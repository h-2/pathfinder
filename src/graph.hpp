#include <list>
#include <map>
#include <ranges>
#include <string>
#include <filesystem>
#include <vector>

#include <bio/alphabet/nucleotide/dna5.hpp>
#include <bio/alphabet/fmt.hpp>
#include <bio/ranges/container/concatenated_sequences.hpp>
#include <bio/ranges/views/char_strictly_to.hpp>
#include <bio/ranges/views/complement.hpp>
#include <bio/ranges/to.hpp>

#include <bio/io/genomic_region.hpp>
#include <bio/io/txt/misc.hpp>
#include <bio/io/txt/reader.hpp>

/* DETAIL*/

struct hash_string
{
    using is_transparent = void;

    std::size_t operator()(const std::string& v) const
    {
        return std::hash<std::string>{}(v);
    }
    std::size_t operator()(const char*v) const
    {
        return std::hash<std::string_view>{}(v);
    }

    std::size_t operator()(const std::string_view &v) const
    {
        return std::hash<std::string_view>{}(v);
    }
};

void assign_append(std::ranges::forward_range auto && source,
                   std::ranges::output_range<std::ranges::range_reference_t<decltype(source)>> auto && sink)
{

    size_t old_sink_size = sink.size();
    sink.resize(old_sink_size+ source.size());
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

// struct store_graph
// {
//     bio::ranges::concatenated_sequences<std::string> stable_sequence_names;
//
//     /* segment lines; all have same size */
//     bio::ranges::concatenated_sequences<std::vector<bio::alphabet::dna5>> seqs;
//     bio::ranges::concatenated_sequences<std::string> names;
//     std::vector<uint64_t> stable_name_i; // -> stable_sequence_names
//     std::vector<uint64_t> stable_offsets;
//     std::vector<uint64_t> ranks;
//
//     // bio::ranges::concatenated_sequences<std::vector<arc>> arcs;
//     std::vector<std::vector<arc>> arcs;
//     std::vector<uint64_t> split_counter;
// };


enum class orientation : bool
{
    plus,
    minus
};

struct node;

struct arc
{
    size_t target_node_i      = -1;
    orientation orient_self   = orientation::plus;
    orientation orient_target = orientation::plus;
};

struct node
{
    std::vector<bio::alphabet::dna5> seq;
    std::string name;

    std::string sn;
    uint64_t so = 0;
    uint64_t sr = 0;

    std::vector<arc> arcs;
};


struct visual_graph
{
    std::vector<node> nodes;
    std::unordered_map<std::string, size_t, hash_string, std::equal_to<void>> name2node_i;
};


using seq_map_t = std::unordered_map<std::string, size_t, hash_string, std::equal_to<void>>;
using arc_map_t = std::unordered_map<std::string, std::vector<arc>, hash_string, std::equal_to<void>>;

inline void parse_segment_optional(std::string_view const field,
                                   auto & checks,
                                   node & n)
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
        auto [ptr, ec] = std::from_chars(so.data(), so.data() + so.size(), offset);
        if (ec != std::errc{})
            throw std::runtime_error{fmt::format("Expected number, got: {}", so)};

        n.so = offset;
        checks.have_SO = true;
    }
    else if (field.starts_with("SR:i:"))
    {
        std::string_view sr = field.substr(5);
        if (sr.empty())
            throw std::runtime_error{"SR field may not be empty."};

        uint64_t rank = 0;
        auto [ptr, ec] = std::from_chars(sr.data(), sr.data() + sr.size(), rank);
        if (ec != std::errc{})
            throw std::runtime_error{fmt::format("Expected number, got: {}", sr)};

        n.sr = rank;
        checks.have_SR = true;
    }
    // other optionals are ignored
}

inline void read_gfa(std::filesystem::path const & path, visual_graph & graph)
{

    auto reader = path == "-" ?
                  bio::io::txt::reader{std::cin, '\t', bio::io::txt::header_kind::starts_with{'#'}} :
                  bio::io::txt::reader{path, '\t', bio::io::txt::header_kind::starts_with{'#'}};

    // seq_map_t segment_names_map;
    // seq_map_t stable_sequence_names_map;
    // arc_map_t arc_map;

    for (bio::io::txt::record & r : reader)
    {
        auto error = [&r] <typename ... arg_ts> (fmt::format_string<arg_ts...> fmt_str, arg_ts ... args)
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
                n.seq = r.fields[2] | bio::views::char_strictly_to<bio::alphabet::dna5> | bio::ranges::to<std::vector>();

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
}

enum class reg_placement
{
    ends_before,
    ends_in,
    lies_in,
    contains,
    begins_in,
    begins_after
};

inline reg_placement overlap(std::pair<size_t, size_t> region1, std::pair<size_t, size_t> region2)
{
    //TODO this can probably be optimised with arithmetics
    if (region1.second <= region2.first)
        return reg_placement::ends_before;
    else if (region1.first >= region2.second)
        return reg_placement::begins_after;
    else if (region1.first < region2.first)
        return (region1.second >= region2.second) ? reg_placement::contains : reg_placement::ends_in;
    else if (region1.second > region2.second)
        return reg_placement::begins_in;
    else
        return reg_placement::lies_in;
}

//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------

inline void print_all_paths(visual_graph const & graph)
{
    auto recu = [&] (auto recu,
                     std::vector<bio::alphabet::dna5> seq,
                     std::string path,
                     std::vector<bio::io::genomic_region> stable_path,
                     node const & n,
                     orientation const o) -> void
    {
        path += "->";
        path += n.name;

        if (!stable_path.empty() && stable_path.back().chrom == n.sn && stable_path.back().end == (int64_t)n.so)
            stable_path.back().end += n.seq.size();
        else
            stable_path.push_back({.chrom = n.sn, .beg = (int64_t)n.so, .end = (int64_t)n.so + (int64_t)n.seq.size()});


        if (o == orientation::plus)
            assign_append(n.seq, seq);
        else
            assign_append(n.seq | std::views::reverse | bio::views::complement, seq);

        if (n.arcs.empty()) // leaf node
        {
            fmt::print("Leaf node\n  path: {}\n  seq: {}\n  stable: {}\n", path, seq, stable_path);
        }
        else
        {
            for (arc const a : n.arcs)
                recu(recu, seq, path, stable_path, graph.nodes[a.target_node_i], a.orient_target);
        }
    };

    recu(recu, {}, {}, {}, graph.nodes.front(), orientation::plus);
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


inline visual_graph discretise(visual_graph const & in_graph, size_t const w)
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
        node & n  = out_graph.nodes[i];

        for (arc & a : n.arcs)
            if (a.orient_target == orientation::minus)
                needs_flipping.insert(a.target_node_i);
    }

    // copy and flip nodes
    for (size_t i : needs_flipping)
    {
        node & n  = out_graph.nodes[i];

        node new_n{.seq = n.seq | std::views::reverse | bio::views::complement | bio::ranges::to<std::vector>(),
                   .name = n.name + "_rc",
                   .sn = n.sn,
                   .so = n.so,
                   .sr = n.sr,
                   .arcs{}};

        /* copy arcs and reorient TODO replace the following with remove_copy_if*/
        std::ranges::copy_if(n.arcs,
                             std::back_inserter(new_n.arcs),
                             [] (arc const & a) { return a.orient_self == orientation::minus; });
        for (arc & a : new_n.arcs)
            a.orient_self = orientation::plus;

        /* remove arcs from old */
        std::erase_if(n.arcs, [] (arc const & a) { return a.orient_self == orientation::minus;});

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

    for (size_t i = 0; i < out_graph.nodes.size(); ++i) // size will grow during iteration!
    {
        node & n  = out_graph.nodes[i];
        if (n.sr == 0ul) // ref-node
        {
            if ((n.so / w) != (n.so + n.seq.size()) / w) // needs to be split
            {
                // we only do one split here, if more are necessary, they happen later automatically
                size_t new_offset = n.so / w * w + w; // round up to next multiple of w

                node new_n{.seq = n.seq | std::views::drop(new_offset - n.so) | bio::ranges::to<std::vector>(),
                           .name = n.name + "_n",
                           .sn = n.sn,
                           .so = new_offset,
                           .sr = 0ul,
                           .arcs = std::move(n.arcs),
                };

                n.seq.resize(new_offset - n.so);
                n.arcs.clear(); // this should be implicit by move but is not guaranteed
                n.arcs.push_back(arc{.target_node_i = out_graph.nodes.size()});

                out_graph.name2node_i[new_n.name] = out_graph.nodes.size();
                out_graph.nodes.push_back(std::move(new_n));
            }
        }
    }

    /* Step 3a:
     *
     * Generate all non-ref paths, split those paths and add nodes non-ref nodes
     */

    std::vector<std::vector<size_t>> paths;

    auto fn = [&out_graph, &paths] (auto self, std::vector<size_t> path, size_t const i_node)
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
        if (node const & n = out_graph.nodes[i]; n.sr == 0ul) // ref-node
            for (arc const & a : n.arcs)
                if (node const & target = out_graph.nodes[a.target_node_i]; target.sr != 0ul) // non-ref node
                    fn(fn, {i}, a.target_node_i);

    // DEBUG
    for (std::vector<size_t> & path : paths)
    {
        fmt::print("Non-Ref-Path: {}\n", path | std::views::transform([&] (size_t const i) -> std::string_view { return out_graph.nodes[i].name; }));
    }

    for (std::vector<size_t> & path : paths)
    {
        auto seqs = path | std::views::transform([&](size_t i) -> auto & { return out_graph.nodes[i].seq; });
        auto sizes = seqs | std::views::transform([](auto && span) { return span.size(); });

        size_t const total_length = std::reduce(sizes.begin(), sizes.end(), 0ul);
        size_t const n_new_nodes = (total_length + (w / 2)) / w;
        size_t const actual_w = (total_length + n_new_nodes - 1) / n_new_nodes; // this is w+-0.5w (ideally == w)

        assert(n_new_nodes >= 1);
        assert(double(actual_w) >= 0.5 * w);
        assert(double(actual_w) <= 1.5 * w);

        size_t i_old_node = 0ul;
        size_t old_node_offset = 0ul;

        for (size_t i_new_node = 0ul; i_new_node < n_new_nodes; ++i_new_node)
        {
            node new_node{};
            // name and seq are set below
            new_node.sn = "virtual";
            new_node.so = -1;
            new_node.sr = -1;
            if (i_new_node < n_new_nodes - 1) // every iteration but the last
                new_node.arcs.push_back(arc{.target_node_i = out_graph.nodes.size() + 1});

            assert(new_node.seq.size() < actual_w && i_old_node < path.size());
            while(true)
            {
                node & old_node = out_graph.nodes[path[i_old_node]];
                assert(old_node_offset < old_node.seq.size());
                std::span old_seq{old_node.seq | std::views::drop(old_node_offset)};

                if (size_t required_size = actual_w - new_node.seq.size(); required_size >= old_seq.size()) // consume old_node
                {
                    /* setup new node */
                    std::ranges::copy(old_seq, std::back_inserter(new_node.seq));
                    new_node.name += old_node.name;
                    if (old_seq.size() != old_node.seq.size())
                        new_node.name += "_sub";

                    /* handle old node */
                    if (old_node.sr != 0ul)
                        old_node.sn = "\0"; // mark for deletion (only non-ref nodes)

                    /* upkeep */
                    ++i_old_node; // go to next old_node
                    old_node_offset = 0ul;
                }
                else // (required_size < old_seq.size())
                {
                    /* setup new node */
                    std::ranges::copy(old_seq | std::views::take(required_size),
                                      std::back_inserter(new_node.seq));
                    new_node.name += old_node.name;
                    new_node.name += "_sub";

                    /* handle old node */
                    old_node_offset += required_size;

                    /* upkeep */
                    // n_i unchanged; do not go to next old_node

                }

                if (new_node.seq.size() < actual_w && i_old_node < path.size())
                    new_node.name += ","; // there will be more
                else
                    break;
            }

#ifndef NDEBUG
            if (i_new_node < n_new_nodes - 1) // every iteration but the last
                assert(new_node.seq.size() == actual_w);

#endif
            out_graph.name2node_i[new_node.name] = out_graph.nodes.size();
            out_graph.nodes.push_back(std::move(new_node));
        }


        // add in-arcs to path
        size_t const first_path_node_i = path.front();
        size_t const first_new_node_i = out_graph.nodes.size() - n_new_nodes;
        for (node & n : out_graph.nodes)
            for (arc const & a : n.arcs)
                if (a.target_node_i == first_path_node_i)
                    n.arcs.push_back(arc{.target_node_i = first_new_node_i}), ({break;}); // :-)

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
        if (n.sn == "\0")
        {
            n.arcs.clear();
        }
        else
        {
            std::erase_if(n.arcs, [&] (arc const & a)
            {
                return out_graph.nodes[a.target_node_i].sn == "\0";
            });
        }
    }

    // compute for every graph_node the number of graph_nodes before it that will be removed
    std::vector<size_t> del_sums;
    del_sums.resize(out_graph.nodes.size());
    for (size_t i = 1; i < out_graph.nodes.size(); ++i)
        del_sums[i] = del_sums[i - 1] + (out_graph.nodes[i - 1].sn == "\0");

    // substract from every arc target index the specific modifier
    for (node & n : out_graph.nodes)
        for (arc & a : n.arcs)
            a.target_node_i -= del_sums[a.target_node_i];

    // actually remove the nodes
    std::erase_if(out_graph.nodes, [] (node const & n) { return n.sn == "\0"; });

    return out_graph;
}



// inline void print_chunks(graph const & grph, size_t const w)
// {
//     auto recu = [&grph, &w] (auto recu,
//                          std::vector<bio::alphabet::dna5> seq,
//                          std::string path,
//                          size_t i,
//                          orientation const o,
//                          size_t stable_offset,
//                          size_t local_offset) -> void
//     {
//         path += "->";
//         path += grph.names[i];
//
//         int64_t needed = (int64_t)w - seq.size();
//         int64_t have = (int64_t)seq.size() + grph.seqs[i].size() - local_offset;
//
//         if (have > needed) // we are entering next chunk
//         {
//             std::span in  = grph.seqs[i] | std::views::drop(local_offset) | std::views::take(needed);
//             std::span out = grph.seqs[i] | std::views::drop(local_offset + needed);
//
//             if (o == orientation::plus)
//                 fmt::print("Chunk\n  path: {}\n  seq: {}|{}\n", path, seq, in);
//             else
//                 fmt::print("Chunk\n  path: {}\n  seq: {}|{}\n", path, seq, in | std::views::reverse | bio::views::complement);
//
//             seq.assign(out.begin(), out.end());
//             recu(recu, std::move(seq), std::move(path), i, o, stable_offset + w, local_offset + w);
//             return;
//         }
//         else if (grph.arcs[i].empty()) // leaf node
//         {
//             std::span in  = grph.seqs[i] | std::views::drop(local_offset);
//
//             if (o == orientation::plus)
//                 fmt::print("Leaf Chunk\n  path: {}\n  seq: {}|{}\n", path, seq, in);
//             else
//                 fmt::print("Leaf Chunk\n  path: {}\n  seq: {}|{}\n", path, seq, in | std::views::reverse | bio::views::complement);
//         }
//         else
//         {
//             std::span in  = grph.seqs[i] | std::views::drop(local_offset);
//
//             if (o == orientation::plus)
//                 assign_append(in, seq);
//             else
//                 assign_append(in | std::views::reverse | bio::views::complement, seq);
//
//             for (arc const a : grph.arcs[i])
//             {
//                 recu(recu, seq, path, a.target_i, a.orient_target, stable_offset, 0ull);
//             }
//         }
//     };
//
//     recu(recu, {}, {}, 0, orientation::plus, 0, 0);
// }


// inline void discretise(graph & grph, size_t const w)
// {
    // size_t max_offset = *std::ranges::max_element(grph.stable_offsets);
    //
    // graph new_graph;

    /* 1st iteration
     *
     * All reference segments are split on window borders.
     */

//     for (size_t i = 0; i < grph.names.size(); ++i)
//     {
//         size_t stable_b = grph.stable_offsets[i];
//         size_t stable_e = grph.stable_offsets[i] + grph.seqs[i].size();
//
//         if (stable_e / w > stable_b / w) // end is in different window than beginning
//         {
//             /* create new node */
//             // these are just copied
//             grph.names.push_back(grph.names[i]);
//             grph.stable_name_i.push_back(grph.stable_name_i[i]);
//             grph.ranks.push_back(grph.ranks[i]);
//
//             // these have modified values
//             size_t this_len = w - (stable_b % w);
//             grph.seqs.push_back(grph.seqs[i] | std::views::drop(this_len));
//             grph.stable_offsets.push_back(stable_b + this_len); //TODO can we do this?
//             grph.split_counter.push_back(grph.split_counter[i] + 1);
//
//             // rewire arcs
//             grph.arcs.push_back(std::move(grph.arcs[i])); // move current arcs to new node
//
//
//             /* update old node */
//
//
//
//
//             grph.arcs[i].clear();
//             grph.arcs[i].push_back({.target_i = grph.arcs.size()}); // add single arc from old node to new node
//
//
//         }
//     }



    /* 2nd iteration
     *
     * All non-reference segments need to be assigned to one reference window:
     *
     * 1. all non-reference-paths between reference segments are generated, all non-ref-segments in those path
     * that are already "within" exactly one reference window, are marked as such.
     * 2. all paths that span multiple consecutive
     *
     * paths that span n reference windows where n > 1, are split evenly into n sub-paths.
     * As a result, all non-reference segments "belong" to exactly one reference window.
     * Annotate all non-reference segment
     */


    /* 3nd iteration
     *
     * Identify reference windows
     */


     /**
     *
     * As a result all segements except the left-most one begin on a window boundary

    */

    // for (size_t i = 0; i < max_offset; i += w)
    // {
    //     for (size_t j = 0; j < grph.names.size(); ++j)
    //     {
    //         size_t stable_b = grph.stable_offsets[j];
    //         size_t stable_e = grph.stable_offsets[j] + grph.seqs[j].size();
    //
    //         reg_placement place = overlap({stable_b, stable_e}, {i, i + w});
    //
    //         switch (place)
    //         {
    //             case reg_placement::ends_in:
    //             {
    //                 size_t in_reach = stable_e - (stable_e % w);
    //
    //
    //
    //                 // TODO
    //                 break;
    //             }
    //             case reg_placement::lies_in:
    //                 // TODO
    //                 /* pad the sequence */
    //                 new_graph.seqs.push_back();
    //
    //
    //
    //
    //                 break;
    //             case reg_placement::begins_in:
    //                 // TODO
    //                 break;
    //             case reg_placement::ends_before:
    //             case reg_placement::begins_after:
    //                 // NOTHING NEEDS TO BE DONE
    //                 break;
    //         }
    //
    //
    //
    //
    //     }
    //
    //
    //
    //
    // }
    //
    //

// }

