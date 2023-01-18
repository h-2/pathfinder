#include <filesystem>
#include <list>
#include <map>
#include <ranges>
#include <string>
#include <vector>

#include <bio/alphabet/fmt.hpp>
#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/ranges/container/concatenated_sequences.hpp>
#include <bio/ranges/to.hpp>
#include <bio/ranges/views/add_reverse_complement.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>
#include <bio/ranges/views/complement.hpp>

#include <bio/io/genomic_region.hpp>
#include <bio/io/txt/misc.hpp>
#include <bio/io/txt/reader.hpp>

/* DETAIL*/

inline constexpr auto ssize = std::ranges::ssize;

struct hash_string
{
    using is_transparent = void;

    std::size_t operator()(std::string const & v) const { return std::hash<std::string>{}(v); }
    std::size_t operator()(char const * v) const { return std::hash<std::string_view>{}(v); }

    std::size_t operator()(std::string_view const & v) const { return std::hash<std::string_view>{}(v); }
};

void assign_append(std::ranges::forward_range auto &&                                                  source,
                   std::ranges::output_range<std::ranges::range_reference_t<decltype(source)>> auto && sink)
{
    size_t old_sink_size = sink.size();
    sink.resize(old_sink_size + source.size());
    std::ranges::copy(source, sink.begin() + old_sink_size);
}

template <>
struct fmt::formatter<bio::io::genomic_region> : fmt::formatter<std::string>
{
    constexpr auto format(bio::io::genomic_region const & reg, auto & ctx) const
    {
        std::string buf;
        fmt::format_to(std::back_insert_iterator{buf}, "{}:{}-{}", reg.chrom, reg.beg, reg.end);
        return fmt::formatter<std::string>::format(buf, ctx);
    }
};

/* PROGRAM */

enum class orientation : bool
{
    plus,
    minus
};

struct node;

struct arc
{
    size_t      target_node_i = -1;
    orientation orient_self   = orientation::plus;
    orientation orient_target = orientation::plus;
};

struct node
{
    std::vector<bio::alphabet::dna5> seq;
    std::string                      name;

    std::string sn;     // this is set to "chr1:1-4,chr2:3-12" for virtual nodes
    int64_t     so = 0; // this is set to -1 for "virtual nodes"
    int64_t     sr = 0;

    std::vector<arc>                     arcs;
    std::vector<bio::io::genomic_region> regions;
};

inline bio::io::genomic_region node_to_region(node const & n)
{
    return {.chrom = n.sn, .beg = n.so, .end = n.so + ssize(n.seq)};
}

inline void append_region(std::vector<bio::io::genomic_region> &                vec,
                          bio::meta::decays_to<bio::io::genomic_region> auto && reg)
{
    assert(reg.beg != reg.end);

    if (!vec.empty() &&                                               // there is previous region
        vec.back().chrom == reg.chrom && vec.back().end == reg.beg && // they touch
        ((vec.back().beg <= vec.back().end) == (reg.beg <= reg.end))) // same orientation
    {
        vec.back().end = reg.end;
    }
    else
    {
        vec.push_back(std::forward<decltype(reg)>(reg));
    }
}

struct visual_graph
{
    std::vector<node>                                                         nodes;
    std::unordered_map<std::string, size_t, hash_string, std::equal_to<void>> name2node_i;
};

using seq_map_t = std::unordered_map<std::string, size_t, hash_string, std::equal_to<void>>;
using arc_map_t = std::unordered_map<std::string, std::vector<arc>, hash_string, std::equal_to<void>>;

inline void parse_segment_optional(std::string_view const field, auto & checks, node & n)
{
    if (field.starts_with("SN:Z:"))
    {
        std::string_view sn = field.substr(5);
        if (sn.empty())
            throw std::runtime_error{"SN field may not be empty."};

        n.sn = static_cast<std::string>(sn);

        checks.have_SN = true;
    }
    else if (field.starts_with("SO:i:"))
    {
        std::string_view so = field.substr(5);
        if (so.empty())
            throw std::runtime_error{"SO field may not be empty."};

        uint64_t offset = 0;
        auto [ptr, ec]  = std::from_chars(so.data(), so.data() + so.size(), offset);
        if (ec != std::errc{})
            throw std::runtime_error{fmt::format("Expected number, got: {}", so)};

        n.so           = offset;
        checks.have_SO = true;
    }
    else if (field.starts_with("SR:i:"))
    {
        std::string_view sr = field.substr(5);
        if (sr.empty())
            throw std::runtime_error{"SR field may not be empty."};

        uint64_t rank  = 0;
        auto [ptr, ec] = std::from_chars(sr.data(), sr.data() + sr.size(), rank);
        if (ec != std::errc{})
            throw std::runtime_error{fmt::format("Expected number, got: {}", sr)};

        n.sr           = rank;
        checks.have_SR = true;
    }
    // other optionals are ignored
}

inline void read_gfa(std::filesystem::path const & path, visual_graph & graph)
{
    auto reader = path == "-" ? bio::io::txt::reader{std::cin, '\t', bio::io::txt::header_kind::starts_with{'#'}}
                              : bio::io::txt::reader{path, '\t', bio::io::txt::header_kind::starts_with{'#'}};

    // seq_map_t segment_names_map;
    // seq_map_t stable_sequence_names_map;
    // arc_map_t arc_map;

    for (bio::io::txt::record & r : reader)
    {
        auto error = [&r]<typename... arg_ts>(fmt::format_string<arg_ts...> fmt_str, arg_ts... args)
        {
            std::string s = "rGFA parse error.\n" + fmt::format(fmt_str, std::forward<arg_ts>(args)...);
            s += "\nSource line:\n";
            s += r.line;
            throw std::runtime_error{std::move(s)};
        };

        if (r.fields.empty())
            error("Don't know how to handle empty line in rGFA.");

        if (r.fields[0].size() != 1)
            error("Don't know how to handle '{}' (size {}) as first field.", r.fields[0], r.fields[0].size());

        switch (r.fields[0][0])
        {
            case 'S': /* SEGMENT */
                {
                    node n;

                    if (r.fields.size() < 6)
                        error("At least 6 fields required for SEGMENT record.");

                    if (auto it = graph.name2node_i.find(r.fields[1]); it != graph.name2node_i.end())
                    {
                        error("Segment names have to be unique, but this one isn't.");
                    }

                    n.name = static_cast<std::string>(r.fields[1]);
                    n.seq =
                      r.fields[2] | bio::views::char_strictly_to<bio::alphabet::dna5> | bio::ranges::to<std::vector>();

                    struct
                    {
                        bool have_SN = false;
                        bool have_SO = false;
                        bool have_SR = false;
                    } checks;

                    for (size_t i = 3; i < r.fields.size(); ++i)
                        parse_segment_optional(r.fields[i], checks, n);

                    if (!checks.have_SN)
                        error("rGFA requires SN optional but none was present.");
                    if (!checks.have_SO)
                        error("rGFA requires SO optional but none was present.");
                    if (!checks.have_SR)
                        error("rGFA requires SR optional but none was present.");

                    // assert(n.seqs.size() == n.names.size());
                    // assert(n.seqs.size() == n.stable_name_i.size());
                    // assert(n.seqs.size() == n.stable_offsets.size());
                    // assert(n.seqs.size() == n.ranks.size());

                    graph.nodes.push_back(std::move(n));
                    graph.name2node_i[graph.nodes.back().name] = graph.nodes.size() - 1;
                    break;
                }
            case 'L': /* LINKE */
                {
                    if (r.fields.size() < 6)
                        error("At least 6 fields required for LINK record.");

                    std::string_view source_seg_id     = r.fields[1];
                    std::string_view source_seg_orient = r.fields[2];
                    std::string_view target_seg_id     = r.fields[3];
                    std::string_view target_seg_orient = r.fields[4];
                    std::string_view cigar             = r.fields[5];

                    arc a;

                    if (auto it = graph.name2node_i.find(target_seg_id); it == graph.name2node_i.end())
                        error("Target segment ID '{}' unknown.", target_seg_id);
                    else
                        a.target_node_i = it->second;

                    if (source_seg_orient == "+")
                        a.orient_self = orientation::plus;
                    else if (source_seg_orient == "-")
                        a.orient_self = orientation::minus;
                    else
                        error("Source orientation '{}' invalid, must be '+' or '-'.", source_seg_orient);

                    if (target_seg_orient == "+")
                        a.orient_target = orientation::plus;
                    else if (target_seg_orient == "-")
                        a.orient_target = orientation::minus;
                    else
                        error("Source orientation '{}' invalid, must be '+' or '-'.", target_seg_orient);

                    if (cigar != "0M")
                        error("Overlaps are not supported.");

                    size_t source_node_i = -1;

                    if (auto it = graph.name2node_i.find(source_seg_id); it == graph.name2node_i.end())
                        error("Source segment ID '{}' unknown.", source_seg_id);
                    else
                        source_node_i = it->second;

                    graph.nodes[source_node_i].arcs.push_back(a);

                    break;
                }
            case '#': /* COMMENT */
                break;
            default:
                error("Don't know how to handle '{}' as first field.", r.fields[0]);
                break;
        }
    }

    for (node & n : graph.nodes)
    {
        n.regions.push_back(node_to_region(n));
    }
}

//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------

inline void print_all_paths(visual_graph const & graph)
{
    auto fn =
      [&graph](auto self, std::vector<std::pair<size_t, orientation>> path, size_t const i_node, orientation const o)
    {
        path.emplace_back(i_node, o);
        node const & n = graph.nodes[i_node];

        if (n.arcs.empty()) // no children → at end
        {
            auto to_name = [&graph](std::pair<size_t, orientation> const & p) -> auto &
            { return graph.nodes[p.first].name; };
            auto to_seq = [&graph](std::pair<size_t, orientation> const & p)
            {
                node const & n2 = graph.nodes[p.first];
                return n2.seq | bio::ranges::detail::reverse_complement_or_not(p.second == orientation::minus);
            };
            // auto to_sn = [&graph] (std::pair<size_t, orientation> const & p) -> std::string
            // {
            //     return graph.nodes[p.first].so >= 0 ?
            //         fmt::format("{}:{}-{}", graph.nodes[p.first].sn, graph.nodes[p.first].so, graph.nodes[p.first].so + graph.nodes[p.first].seq.size()) : graph.nodes[p.first].sn;
            // };

            std::vector<bio::io::genomic_region> regions;
            for (auto && [i_node, o] : path)
            {
                switch (o)
                {
                    case orientation::plus:
                        for (bio::io::genomic_region const & reg : graph.nodes[i_node].regions)
                            append_region(regions, reg);
                        break;
                    case orientation::minus:
                        for (bio::io::genomic_region const & reg : graph.nodes[i_node].regions | std::views::reverse)
                            append_region(regions, bio::io::genomic_region{reg.chrom, reg.end, reg.beg});
                        break;
                }
            }

            fmt::print("Leaf node\n  names: {}\n  seqs: {}\n  regions: {}\n",
                       path | std::views::transform(to_name),
                       path | std::views::transform(to_seq),
                       regions);
            return;
        }

        for (arc const & a : n.arcs | std::views::take(n.arcs.size() - 1))
            self(self, path, a.target_node_i, a.orient_target); // pass path as copy to n-1 children

        // move path to last child
        self(self, std::move(path), n.arcs.back().target_node_i, n.arcs.back().orient_target);
    };

    fn(fn, {}, 0, orientation::plus);
}

inline void graph_to_dot(visual_graph const & graph)
{
    fmt::print("digraph {}\n", "{");

    for (node const & n : graph.nodes)
    {
        fmt::print("{} [label =\"{}:{}\"]\n", n.name, n.name, n.seq);
        for (arc const & a : n.arcs)
            fmt::print("{} -> {}\n", n.name, graph.nodes[a.target_node_i].name);
    }

    fmt::print("{}\n", "}");
}

inline visual_graph discretise(visual_graph const & in_graph, int64_t const w)
{
    // copy in_graph
    visual_graph out_graph = in_graph;

    /* Step 1:
     *
     * Each node that is traversed in reverse order is copied, reverse-complemented and its traversal
     * order set to normal.
     */

    std::unordered_set<size_t> needs_flipping;
    for (size_t i = 0; i < out_graph.nodes.size(); ++i) // size will grow during iteration!
    {
        node & n = out_graph.nodes[i];

        for (arc & a : n.arcs)
            if (a.orient_target == orientation::minus)
                needs_flipping.insert(a.target_node_i);
    }

    // copy and flip nodes
    for (size_t i : needs_flipping)
    {
        node & n = out_graph.nodes[i];

        node new_n{.seq  = n.seq | std::views::reverse | bio::views::complement | bio::ranges::to<std::vector>(),
                   .name = n.name + "_rc",
                   .sn   = n.sn,
                   .so   = n.so,
                   .sr   = n.sr,
                   .arcs{},
                   .regions = n.regions};

        assert(new_n.regions.size() == 1);
        std::swap(new_n.regions.front().beg, new_n.regions.front().end); // reverse coordinates for RC

        /* copy arcs and reorient TODO replace the following with remove_copy_if*/
        std::ranges::copy_if(n.arcs,
                             std::back_inserter(new_n.arcs),
                             [](arc const & a) { return a.orient_self == orientation::minus; });
        for (arc & a : new_n.arcs)
            a.orient_self = orientation::plus;

        /* remove arcs from old */
        std::erase_if(n.arcs, [](arc const & a) { return a.orient_self == orientation::minus; });

        out_graph.name2node_i[new_n.name] = out_graph.nodes.size();
        out_graph.nodes.push_back(std::move(new_n));
    }

    // re-orient existing arcs
    for (node & n : out_graph.nodes)
    {
        for (arc & a : n.arcs)
        {
            if (a.orient_target == orientation::minus)
            {
                std::string new_target_name = out_graph.nodes[a.target_node_i].name + "_rc";
                assert(out_graph.name2node_i.contains(new_target_name));
                a.target_node_i = out_graph.name2node_i[new_target_name];
                a.orient_target = orientation::plus;
            }
        }
    }

    // TODO mark leftover rc'ed nodes to-be-deleted

    /* Step 2:
     *
     * Split all ref-nodes on block-borders
     */

    // auto fac = [] (node const & n)
    // {
    //     assert(n.regions.size() > 0);
    //     return n.regions.back().beg <= n.regions.back().end ? 1 : -1;
    // };

    for (size_t i = 0; i < out_graph.nodes.size(); ++i) // size will grow during iteration!
    {
        node & n = out_graph.nodes[i];
        if (n.sr == 0ul && n.regions.back().beg <= n.regions.back().end) // ref-node
        {
            if (n.regions.back().beg <= n.regions.back().end)            //  plus-orientation
            {
                if ((n.so / w) != (n.so + ssize(n.seq) - 1) / w)         // needs to be split
                {
                    // we only do one split here, if more are necessary, they happen later automatically
                    int64_t new_offset = n.so / w * w + w; // round up to next multiple of w

                    auto new_seq = n.seq | std::views::drop(new_offset - n.so);
                    node new_n{.seq  = new_seq | bio::ranges::to<std::vector>(),
                               .name = n.name + "_n",
                               .sn   = n.sn,
                               .so   = new_offset,
                               .sr   = 0l,
                               .arcs = std::move(n.arcs),
                               .regions{{.chrom = n.sn, .beg = new_offset, .end = new_offset + ssize(new_seq)}}};

                    n.seq.resize(n.seq.size() - new_seq.size());
                    n.arcs.clear(); // this should be implicit by move but is not guaranteed
                    n.arcs.push_back(arc{.target_node_i = out_graph.nodes.size()});
                    assert(n.regions.size() == 1);
                    n.regions.front().end = new_offset;
                    assert(n.regions.front().end - n.regions.front().beg == ssize(n.seq));

                    out_graph.name2node_i[new_n.name] = out_graph.nodes.size();
                    out_graph.nodes.push_back(std::move(new_n));
                }
            }
            else // minus-orientation; this shouldn't happen on linear reference, but who knows
            {
                assert(n.so + 1 >= ssize(n.seq));
                if ((n.so / w) != (n.so + 1 - ssize(n.seq)) / w) // needs to be split
                {
                    // we only do one split here, if more are necessary, they happen later automatically
                    int64_t new_old_offset = n.so / w; // round down to next multiple of w

                    auto new_seq = n.seq | std::views::drop(n.so - new_old_offset);
                    node new_n{.seq  = new_seq | bio::ranges::to<std::vector>(),
                               .name = n.name + "_n",
                               .sn   = n.sn,
                               .so   = n.so,
                               .sr   = 0l,
                               .arcs = std::move(n.arcs),
                               .regions{{.chrom = n.sn, .beg = n.so, .end = n.so - ssize(new_seq)}}};

                    n.seq.resize(n.seq.size() - new_seq.size());
                    n.arcs.clear(); // this should be implicit by move but is not guaranteed
                    n.arcs.push_back(arc{.target_node_i = out_graph.nodes.size()});
                    assert(n.regions.size() == 1);
                    n.so                  = new_old_offset;
                    n.regions.front().beg = new_old_offset;
                    assert(n.regions.front().beg - n.regions.front().end == ssize(n.seq));

                    out_graph.name2node_i[new_n.name] = out_graph.nodes.size();
                    out_graph.nodes.push_back(std::move(new_n));
                }
            }
        }
    }

    /* Step 2.5:
     *
     * Generate regions for all nodes
     */

    /* Step 3a:
     *
     * Generate all non-ref paths, split those paths and add nodes non-ref nodes
     */

    std::vector<std::vector<size_t>> paths;

    auto fn = [&out_graph, &paths](auto self, std::vector<size_t> path, size_t const i_node)
    {
        path.push_back(i_node);
        node const & n = out_graph.nodes[i_node];

        if (n.sr == 0ul || n.arcs.empty()) // ref_node or no children → at end
        {
            paths.push_back(std::move(path));
            return;
        }

        for (arc const & a : n.arcs | std::views::take(n.arcs.size() - 1))
            self(self, path, a.target_node_i); // pass path as copy to n-1 children

        // move path to last child
        self(self, std::move(path), n.arcs.back().target_node_i);
    };

    // generate paths
    for (size_t i = 0; i < out_graph.nodes.size(); ++i)
        if (node const & n = out_graph.nodes[i]; n.sr == 0ul)                                 // ref-node
            for (arc const & a : n.arcs)
                if (node const & target = out_graph.nodes[a.target_node_i]; target.sr != 0ul) // non-ref node
                    fn(fn, {i}, a.target_node_i);

    // DEBUG
    for (std::vector<size_t> & path : paths)
    {
        fmt::print("Non-Ref-Path: {}\n",
                   path | std::views::transform([&](size_t const i) -> std::string_view
                                                { return out_graph.nodes[i].name; }));
    }

    for (std::vector<size_t> & path : paths)
    {
        auto seqs  = path | std::views::transform([&](size_t i) -> auto & { return out_graph.nodes[i].seq; });
        auto sizes = seqs | std::views::transform([](auto && span) { return span.size(); });

        int64_t const total_length = std::reduce(sizes.begin(), sizes.end(), 0l);
        int64_t const n_new_nodes  = (total_length + (w / 2)) / w;
        int64_t const actual_w     = (total_length + n_new_nodes - 1) / n_new_nodes; // this is w+-0.5w (ideally == w)

        assert(n_new_nodes >= 1);
        assert(double(actual_w) >= 0.5 * w);
        assert(double(actual_w) <= 1.5 * w);

        size_t  i_old_node      = 0ul;
        int64_t old_node_offset = 0ul;

        for (int64_t i_new_node = 0l; i_new_node < n_new_nodes; ++i_new_node)
        {
            node new_node{};
            // name, seq and .sn are set below
            new_node.so = -1;
            new_node.sr = -1;
            if (i_new_node < n_new_nodes - 1) // every iteration but the last
                new_node.arcs.push_back(arc{.target_node_i = out_graph.nodes.size() + 1});

            assert(ssize(new_node.seq) < actual_w && i_old_node < path.size());
            while (true)
            {
                node & old_node = out_graph.nodes[path[i_old_node]];
                assert(old_node_offset < ssize(old_node.seq));
                std::span     old_seq{old_node.seq | std::views::drop(old_node_offset)};
                int64_t const desired_size   = actual_w - ssize(new_node.seq);
                int64_t const available_size = ssize(old_seq);
                int64_t const taken_size     = std::min(desired_size, available_size);

                /* setup new node */
                new_node.name += old_node.name; // TODO we need to create unique names here!!
                if (taken_size != ssize(old_node.seq))
                    new_node.name += "_sub";

                assert(old_node.regions.size() > 0);
                if (old_node.regions.back().beg <= old_node.regions.back().end)
                {
                    append_region(new_node.regions,
                                  bio::io::genomic_region{old_node.sn,
                                                          old_node.regions.back().beg + old_node_offset,
                                                          old_node.regions.back().beg + old_node_offset + taken_size});
                }
                else
                {
                    assert(old_node.regions.back().beg - old_node.regions.back().end == ssize(old_node.seq));
                    append_region(new_node.regions,
                                  bio::io::genomic_region{old_node.sn,
                                                          old_node.regions.back().beg - old_node_offset,
                                                          old_node.regions.back().beg - old_node_offset - taken_size});
                }

                std::ranges::copy(old_seq | std::views::take(taken_size), std::back_inserter(new_node.seq));

                if (taken_size == available_size) // need to go to next node
                {
                    /* handle old node */
                    if (old_node.sr != 0l)
                        old_node.sr = -2; // mark for deletion (only non-ref nodes)

                    /* upkeep */
                    ++i_old_node; // go to next old_node
                    old_node_offset = 0l;
                }
                else // (taken_size < available_size) | same old_node will be used for next new_node
                {
                    /* handle old node */
                    old_node_offset += taken_size;

                    /* upkeep */
                    // n_i unchanged; do not go to next old_node
                }

                if (ssize(new_node.seq) < actual_w && i_old_node < path.size())
                {
                    new_node.name += ","; // there will be more
                }
                else
                {
                    break;
                }
            }

#ifndef NDEBUG
            if (i_new_node < n_new_nodes - 1) // every iteration but the last
                assert(ssize(new_node.seq) == actual_w);

#endif
            out_graph.name2node_i[new_node.name] = out_graph.nodes.size();
            out_graph.nodes.push_back(std::move(new_node));
        }

        // add in-arcs to path
        size_t const first_path_node_i = path.front();
        size_t const first_new_node_i  = out_graph.nodes.size() - n_new_nodes;
        for (node & n : out_graph.nodes)
            for (arc const & a : n.arcs)
                if (a.target_node_i == first_path_node_i)
                    n.arcs.push_back(arc{.target_node_i = first_new_node_i}), ({ break; }); // :-)

        // add out-arcs from path
        out_graph.nodes.back().arcs = out_graph.nodes[path.back()].arcs;
    }

    /* Step 3b:
     *
     * Remove all old non-ref nodes
     */

    // remove all arcs from to-be-deleted nodes and all arcs to to-be-deleted nodes
    for (node & n : out_graph.nodes)
    {
        if (n.sr == -2) // marked for deletion
        {
            n.arcs.clear();
        }
        else
        {
            std::erase_if(n.arcs, [&](arc const & a) { return out_graph.nodes[a.target_node_i].sr == -2; });
        }
    }

    // compute for every graph_node the number of graph_nodes before it that will be removed
    std::vector<size_t> del_sums;
    del_sums.resize(out_graph.nodes.size());
    for (size_t i = 1; i < out_graph.nodes.size(); ++i)
        del_sums[i] = del_sums[i - 1] + (out_graph.nodes[i - 1].sr == -2);

    // substract from every arc target index the specific modifier
    for (node & n : out_graph.nodes)
        for (arc & a : n.arcs)
            a.target_node_i -= del_sums[a.target_node_i];

    // actually remove the nodes
    std::erase_if(out_graph.nodes, [](node const & n) { return n.sr == -2; });

    return out_graph;
}
