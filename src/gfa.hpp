
#pragma once

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

#include "misc.hpp"
#include "pf_graph.hpp"

namespace gfa
{

enum class orientation : bool
{
    plus,
    minus
};

struct arc
{
    size_t      target_node_i = -1;
    orientation orient_self   = orientation::plus;
    orientation orient_target = orientation::plus;
};

struct gfa_node
{
    std::vector<bio::alphabet::dna5> seq;
    std::string                      name;

    std::string sn;
    int64_t     so = 0;
    int64_t     sr = 0;

    std::vector<arc> arcs;
};

struct gfa_graph
{
    std::vector<gfa_node>                                                     nodes;
    std::unordered_map<std::string, size_t, hash_string, std::equal_to<void>> name2node_i;
};

inline void parse_segment_optional(std::string_view const field, auto & checks, gfa_node & n)
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

inline void read_graph(std::filesystem::path const & path, gfa_graph & graph)
{
    auto reader = path == "-" ? bio::io::txt::reader{std::cin, '\t', bio::io::txt::header_kind::starts_with{'#'}}
                              : bio::io::txt::reader{path, '\t', bio::io::txt::header_kind::starts_with{'#'}};

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
                    gfa_node n;

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

inline bio::io::genomic_region node_to_region(gfa_node const & n)
{
    return {.chrom = n.sn, .beg = n.so, .end = n.so + ssize(n.seq)};
}

inline void gfa_graph2pf_graph(gfa_graph const & ingraph, pf_graph & out_graph)
{
    //                  orig_i  rc_i
    std::unordered_map<size_t, size_t> orig_to_rc;

    /* copy nodes and mark nodes that need to be reverse complemented */
    for (gfa_node const & in : ingraph.nodes)
    {
        node on{.name = in.name,
                .seq  = in.seq,
                .arcs{},
                .regions{node_to_region(in)},
                .is_ref       = (in.sr == 0),
                .tobe_deleted = false};

        for (arc const & a : in.arcs)
        {
            if (a.orient_self == orientation::minus)
                orig_to_rc[out_graph.nodes.size()] = -1; // new index is still unknown
            if (a.orient_target == orientation::minus)
                orig_to_rc[a.target_node_i] = -1; // new index is still unknown

            /* we copy only plus→plus arcs here */
            if (a.orient_self == orientation::plus && a.orient_target == orientation::plus)
                on.arcs.push_back(a.target_node_i);
        }

        out_graph.nodes.push_back(std::move(on));
    }

    /* add rc nodes */
    for (std::pair<size_t const, size_t> & p : orig_to_rc)
    {
        size_t const i    = p.first;
        size_t &     i_rc = p.second;

        gfa_node const & in = ingraph.nodes[i];
        node const &     on = out_graph.nodes[i];

        node on_rc{.name = on.name + "_rc",
                   .seq  = on.seq | std::views::reverse | bio::views::complement | bio::ranges::to<std::vector>(),
                   .arcs{},
                   .regions{{on.regions.front().chrom, on.regions.front().end, on.regions.front().beg}},
                   .is_ref       = on.is_ref,
                   .tobe_deleted = false};

        /* we copy only minus→plus arcs here */
        for (arc const & a : in.arcs)
            if (a.orient_self == orientation::minus && a.orient_target == orientation::plus)
                on_rc.arcs.push_back(a.target_node_i);

        i_rc = out_graph.nodes.size(); // set index to index of new node
        out_graph.nodes.push_back(std::move(on_rc));
    }

    /* add plus→minus arcs */
    for (size_t i = 0; i < ingraph.nodes.size(); ++i)
    {
        gfa_node const & in = ingraph.nodes[i];
        node &           on = out_graph.nodes[i];

        for (arc const & a : in.arcs)
            if (a.orient_self == orientation::plus && a.orient_target == orientation::minus)
                on.arcs.push_back(orig_to_rc.at(a.target_node_i));

        if (orig_to_rc.contains(i) && on.arcs.empty()) // rc was added and original is disconnected → remove
            on.tobe_deleted = true;
    }

    erase_todo_nodes(out_graph);
}

// TODO do we need something to print gfa_graph as well?
// inline void print_all_paths(pf_graph const & graph)
// {
//     auto fn =
//       [&graph](auto self, std::vector<std::pair<size_t, orientation>> path, size_t const i_node, orientation const o)
//     {
//         path.emplace_back(i_node, o);
//         node const & n = graph.nodes[i_node];
//
//         if (n.arcs.empty()) // no children → at end
//         {
//             auto to_name = [&graph](std::pair<size_t, orientation> const & p) -> auto &
//             { return graph.nodes[p.first].name; };
//             auto to_seq = [&graph](std::pair<size_t, orientation> const & p)
//             {
//                 node const & n2 = graph.nodes[p.first];
//                 return n2.seq | bio::ranges::detail::reverse_complement_or_not(p.second == orientation::minus);
//             };
//             auto to_region = [&graph] (std::pair<size_t, orientation> const & p) -> auto &
//             {
//                 return graph.nodes[p.first].regions;
//             };
//
//             std::vector<bio::io::genomic_region> regions;
//             for (auto && [i_node, o] : path)
//             {
//                 switch (o)
//                 {
//                     case orientation::plus:
//                         for (bio::io::genomic_region const & reg : graph.nodes[i_node].regions)
//                             append_region(regions, reg);
//                         break;
//                     case orientation::minus:
//                         for (bio::io::genomic_region const & reg : graph.nodes[i_node].regions | std::views::reverse)
//                             append_region(regions, bio::io::genomic_region{reg.chrom, reg.end, reg.beg});
//                         break;
//                 }
//             }
//
//             fmt::print("PATH\nnames:\t{}\nseqs:\t{}\nregions:\t{}\nseq_join:\t{}\nregion_join:\t{}\n",
//                        fmt::join(path | std::views::transform(to_name), "\t"),
//                        fmt::join(path | std::views::transform(to_seq), "\t"),
//                        fmt::join(path | std::views::transform(to_region), "\t"),
//                        fmt::join(path | std::views::transform(to_seq), ""),
//                        regions);
//             return;
//         }
//
//         for (arc const & a : n.arcs | std::views::take(n.arcs.size() - 1))
//             self(self, path, a.target_node_i, a.orient_target); // pass path as copy to n-1 children
//
//         // move path to last child
//         self(self, std::move(path), n.arcs.back().target_node_i, n.arcs.back().orient_target);
//     };
//
//     fn(fn, {}, 0, orientation::plus);
// }

} // namespace gfa
