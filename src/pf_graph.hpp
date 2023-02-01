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

#include "index.hpp"
#include "misc.hpp"

struct node
{
    std::string                          name;
    std::vector<bio::alphabet::dna5>     seq;
    std::vector<size_t>                  arcs;
    std::vector<bio::io::genomic_region> regions;
    bool                                 is_ref       = false;
    bool                                 tobe_deleted = false;
};

struct pf_graph
{
    std::vector<node> nodes;

    std::unordered_set<size_t> in_nodes; // indexes of nodes with in-degree 0

    //TODO pull sequences out of the nodes
    // std::unordered_map<std::string, std::vector<bio::alphabet::dna5>>
};

//-----------------------------------------------------------------------------
// path generation
//-----------------------------------------------------------------------------

/* non-recursive graph traversal */
inline void dfs(pf_graph const & graph, size_t const start_node, auto check_fn, auto callback_fn)
{
    std::vector<size_t> path{start_node};
    std::vector<size_t> arc_path{0ul};

    do
    {
        node const & n = graph.nodes[path.back()];

        if (check_fn(n)) // end-condition met (e.g. leaf-node)
        {
            callback_fn(path); // perform callback, e.g. copy current path onto list of valid paths

            // go up in the tree
            path.pop_back();
            arc_path.pop_back();
        }
        else if (n.arcs.size() <= arc_path.back()) // parent has no unprocessed arcs left → go up further
        {
            path.pop_back();
            arc_path.pop_back();
        }
        else // we are adding next node from saved arc index
        {
            /* next node */
            size_t next_i = n.arcs[arc_path.back()];
            path.push_back(next_i);
            arc_path.push_back(0ul);

            /* make arc_path at old position refer to next arc */
            assert(arc_path.size() > 1);
            ++arc_path[arc_path.size() - 2];
        }
    }
    while (!path.empty());
}

/* non-recursive version of the above */
inline std::vector<std::vector<size_t>> generate_all_paths(pf_graph const & graph)
{
    std::vector<std::vector<size_t>> ret;

    auto check_fn    = [](node const & n) { return n.arcs.empty(); };
    auto callback_fn = [&](std::vector<size_t> const & path)
    {
        ret.push_back(path); // copy current path onto list of valid paths
    };

    dfs(graph, 0ul, check_fn, callback_fn);

    return ret;
}

/* generates all maximal paths of the form R1--NR*--R2, i.e.
 * * they start on a reference node
 * * are followed by 1-n non-reference nodes
 * * and end on a reference node (unless the last non-reference node has no outgoing arcs).
 *
 * Also creates paths that begin on non-reference nodes that have no incoming arcs.
 */
inline std::vector<std::vector<size_t>> generate_all_non_ref_paths(pf_graph const & graph)
{
    std::vector<std::vector<size_t>> ret;

    auto check_fn = [](node const & n) { return n.is_ref || n.arcs.empty(); };

    // generate paths
    for (size_t i = 0; i < graph.nodes.size(); ++i)
    {
        if (node const & n = graph.nodes[i]; n.is_ref) // ref-node
        {
            auto callback_fn = [&, i](std::vector<size_t> const & path)
            {
                std::vector<size_t> elem{i}; // initial ref_node
                assign_append(path, elem);   // path starts at first non-ref
                ret.push_back(std::move(elem));
            };

            for (size_t const target_node_i : n.arcs)
                if (node const & target = graph.nodes[target_node_i]; !target.is_ref) // non-ref node
                    dfs(graph, target_node_i, check_fn, callback_fn);
        }
        else if (graph.in_nodes.contains(i)) // in-nodes that are not ref-nodes (e.g. disconnected components)
        {
            auto callback_fn = [&](std::vector<size_t> const & path)
            {
                // here we don't need to prepend the ref-node, because there is none
                ret.push_back(path);
            };
            dfs(graph, i, check_fn, callback_fn);
        }
    }

    return ret;
}

inline void print_outlier_nodes(pf_graph const & graph, index_options const & opts)
{
    std::vector<size_t> in_degrees;
    in_degrees.resize(graph.nodes.size());
    for (node const & n : graph.nodes)
        for (size_t const i_target_node : n.arcs)
            ++in_degrees[i_target_node];

    for (size_t i = 0; i < graph.nodes.size(); ++i)
    {
        node const & n = graph.nodes[i];

        if (n.seq.size() > 2 * opts.window_size)
        {
            fmt::print("Outlier {}\n", i);
            fmt::print("  name:         {}\n", n.name);
            fmt::print("  regions:      {}\n", n.regions);
            fmt::print("  arcs:         {}\n", n.arcs);
            fmt::print("  seq:          {}\n", n.seq);
            fmt::print("  is_ref:       {}\n", n.is_ref);
            fmt::print("  tobe_deleted: {}\n", n.tobe_deleted);
            fmt::print("  in_degree:    {}\n", in_degrees[i]);
        }
    }
}

inline std::vector<std::vector<bio::io::genomic_region>> all_paths_as_collapsed_regions(pf_graph const & graph)
{
    std::vector<std::vector<bio::io::genomic_region>> ret;

    auto check_fn    = [](node const & n) { return n.arcs.empty(); }; // all leaf-nodes
    auto callback_fn = [&](std::vector<size_t> const & path)
    {
        ret.emplace_back();
        for (size_t i_node : path)
            append_regions(graph.nodes[i_node].regions, ret.back());
    };

    dfs(graph, 0ul, check_fn, callback_fn);

    return ret;
}

//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------

inline void print_paths(pf_graph const & graph, std::vector<std::vector<size_t>> const & paths)
{
    fmt::print("PATHS\n======\n");

    auto to_name = [&graph](size_t const i) -> auto &
    {
        return graph.nodes[i].name;
    };
    auto to_seq = [&graph](size_t const i) -> auto &
    {
        return graph.nodes[i].seq;
    };
    auto to_region = [&graph](size_t const i) -> auto &
    {
        return graph.nodes[i].regions;
    };

    std::vector<bio::io::genomic_region> regions;

    for (std::vector<size_t> const & path : paths)
    {
        regions.clear();

        for (size_t const i : path)
            for (bio::io::genomic_region const & reg : graph.nodes[i].regions)
                append_region(reg, regions);

        fmt::print("path\nnames:\t{}\nseqs:\t{}\nregions:\t{}\nseq_join:\t{}\nregion_join:\t{}\n\n",
                   fmt::join(path | std::views::transform(to_name), "\t"),
                   fmt::join(path | std::views::transform(to_seq), "\t"),
                   fmt::join(path | std::views::transform(to_region), "\t"),
                   fmt::join(path | std::views::transform(to_seq), ""),
                   regions);
    }

    fmt::print("\n");
}

inline void print_all_paths(pf_graph const & graph)
{
    std::vector<std::vector<size_t>> paths = generate_all_paths(graph);
    print_paths(graph, paths);
}

inline void graph_to_dot(pf_graph const & graph)
{
    fmt::print("digraph {}\n", "{");

    for (node const & n : graph.nodes)
    {
        fmt::print("\"{}\" [label =\"{}:{}\"]\n", n.name, n.name, n.seq);
        for (size_t const target_node_i : n.arcs)
            fmt::print("\"{}\" -> \"{}\"\n", n.name, graph.nodes[target_node_i].name);
    }

    fmt::print("{}\n\n", "}");
}

inline void print_stats(pf_graph const & graph)
{
    fmt::print("STATS\n======\n");

    fmt::print("Node count:\t{}\n", graph.nodes.size());

    std::vector<size_t> seq_lengths =
      graph.nodes | std::views::transform([](node const & n) { return n.seq.size(); }) | bio::ranges::to<std::vector>();
    std::ranges::sort(seq_lengths);

    fmt::print(
      "Seq lengths\n| Q1\t| Q5\t| Q25\t| Q50\t | Q75\t | Q95\t | Q99\t |\n  {}\t  {}\t  {}\t  {}\t  {}\t  {}\t  {}\n",
      seq_lengths[seq_lengths.size() * 0.01],
      seq_lengths[seq_lengths.size() * 0.05],
      seq_lengths[seq_lengths.size() * 0.25],
      seq_lengths[seq_lengths.size() * 0.50],
      seq_lengths[seq_lengths.size() * 0.75],
      seq_lengths[seq_lengths.size() * 0.95],
      seq_lengths[seq_lengths.size() * 0.99]);
    fmt::print("| min\t| avg  \t| max\t|\n   {}\t  {:.2f}\t  {}\n",
               seq_lengths.front(),
               std::reduce(seq_lengths.begin(), seq_lengths.end(), 0.0) / seq_lengths.size(),
               seq_lengths.back());
    fmt::print("\n");
}

//-----------------------------------------------------------------------------
// graph modification
//-----------------------------------------------------------------------------

inline void erase_todo_nodes(pf_graph & graph)
{
    // remove all arcs from to-be-deleted nodes and all arcs to to-be-deleted nodes
    for (node & n : graph.nodes)
    {
        if (n.tobe_deleted)
            n.arcs.clear();
        else
            std::erase_if(n.arcs, [&](size_t const i) { return graph.nodes[i].tobe_deleted; });
    }

    // compute for every graph_node the number of graph_nodes before it that will be removed
    std::vector<size_t> del_sums;
    del_sums.resize(graph.nodes.size());
    for (size_t i = 1; i < graph.nodes.size(); ++i)
        del_sums[i] = del_sums[i - 1] + (graph.nodes[i - 1].tobe_deleted);

    // substract from every arc index the specific modifier
    for (node & n : graph.nodes)
        for (size_t & i : n.arcs)
            i -= del_sums[i];

    // actually remove the nodes
    std::erase_if(graph.nodes, [](node const & n) { return n.tobe_deleted; });
}

inline void discretise(pf_graph & graph, index_options const & opts)
{
    int64_t const w = opts.window_size;

    fmt::print(stderr, "source_location: {}:{}\n", __FILE__, __LINE__);
    std::cerr << std::flush;

    /* Step 0:
     * determine nodes with in-degree 0
     *
     */

    graph.in_nodes.clear();
    std::vector<size_t> in_degrees;
    in_degrees.resize(graph.nodes.size());
    for (node & n : graph.nodes)
        for (size_t const i_target_node : n.arcs)
            ++in_degrees[i_target_node];
    for (size_t i = 0; i < in_degrees.size(); ++i)
        if (in_degrees[i] == 0)
            graph.in_nodes.insert(i);

    /* Step 1:
     *
     * Split all ref-nodes on block-borders
     */

    for (size_t i = 0; i < graph.nodes.size(); ++i) // size will grow during iteration!
    {
        node & n = graph.nodes[i];
        assert(n.regions.size() == 1);
        if (bio::io::genomic_region & reg = n.regions.back(); n.is_ref) // ref-node
        {
            // we only do one split here; if more are necessary, they happen later automatically
            int64_t split_point   = 0;
            int64_t old_node_size = 0;

            if (reg.beg <= reg.end) // regular orientation
            {
                if ((reg.beg / w) == (reg.end - 1) / w) // doesn't need to be split
                    continue;

                split_point   = reg.beg / w * w + w; // round up to next multiple of w
                old_node_size = split_point - reg.beg;
            }
            else // reverse complemented
            {
                assert(reg.beg + 1 >= ssize(n.seq));
                if ((reg.end / w) == (reg.beg - 1) / w) // doesn't need to be split
                    continue;

                split_point = reg.beg % w == 0     // if already on block border
                                ? reg.beg - w      // take full block; else
                                : reg.beg / w * w; // round down to next multiple of w

                old_node_size = reg.beg - split_point;
            }

            node new_n{.name{},
                       .seq  = n.seq | std::views::drop(old_node_size) | bio::ranges::to<std::vector>(),
                       .arcs = std::move(n.arcs),
                       .regions{{.chrom = reg.chrom, .beg = split_point, .end = reg.end}},
                       .is_ref       = true,
                       .tobe_deleted = false};
            if (opts.descriptive_node_names)
                fmt::format_to(std::back_inserter(new_n.name), "v{}_{}", graph.nodes.size(), n.name);
            else
                fmt::format_to(std::back_inserter(new_n.name), "v{}", graph.nodes.size());

            /* update old node */
            n.seq.resize(old_node_size);
            n.arcs.clear(); // this should be implicit by move but is not guaranteed
            n.arcs.push_back(graph.nodes.size());
            reg.end = split_point;

#ifndef NDEBUG
            bio::io::genomic_region const & new_reg = new_n.regions.back();
            assert(std::max(new_reg.end, new_reg.beg) - std::min(new_reg.end, new_reg.beg) == ssize(new_n.seq));

            assert(std::max(reg.end, reg.beg) - std::min(reg.end, reg.beg) == ssize(n.seq));
#endif

            /* append new node; this will be processed itself later on in the loop */
            graph.nodes.push_back(std::move(new_n));
        }
    }

    fmt::print(stderr, "source_location: {}:{}\n", __FILE__, __LINE__);
    std::cerr << std::flush;

    /* Step 2:
     *
     * Generate all non-ref paths, split those paths and add new nodes and arcs for them
     */

    {
        std::vector<std::vector<size_t>> paths = generate_all_non_ref_paths(graph);

        // DEBUG
        // for (std::vector<size_t> & path : paths)
        // {
        //     fmt::print(stderr, "Non-Ref-Path: {}\n",
        //                path | std::views::transform([&](size_t const i) -> std::string_view
        //                                             { return graph.nodes[i].name; }));
        // }

        for (std::vector<size_t> & path : paths)
        {
            auto seqs  = path | std::views::transform([&](size_t i) -> auto & { return graph.nodes[i].seq; });
            auto sizes = seqs | std::views::transform([](auto && span) { return span.size(); });

            int64_t const total_length = std::reduce(sizes.begin(), sizes.end(), 0l);
            if (total_length <= w)
                continue;
            int64_t const n_new_nodes = (total_length + (w / 2)) / w;
            int64_t const actual_w = (total_length + n_new_nodes - 1) / n_new_nodes; // this is w+-0.5w (ideally == w)

            assert(n_new_nodes >= 1);
            assert(double(actual_w) >= 0.5 * w);
            assert(double(actual_w) <= 1.5 * w);

            size_t  i_old_node      = 0ul;
            int64_t old_node_offset = 0l;

            for (int64_t i_new_node = 0l; i_new_node < n_new_nodes; ++i_new_node)
            {
                node new_node{};
                fmt::format_to(std::back_inserter(new_node.name), "v{}_", graph.nodes.size() + i_new_node);

                // name, seq and .sn are set below
                if (i_new_node < n_new_nodes - 1) // every iteration but the last
                    new_node.arcs.push_back(graph.nodes.size() + 1);

                assert(ssize(new_node.seq) < actual_w && i_old_node < path.size());
                while (true)
                {
                    node & old_node = graph.nodes[path[i_old_node]];
                    assert(old_node_offset < ssize(old_node.seq));
                    std::span     old_seq{old_node.seq | std::views::drop(old_node_offset)};
                    int64_t const desired_size   = actual_w - ssize(new_node.seq);
                    int64_t const available_size = ssize(old_seq);
                    int64_t const taken_size     = std::min(desired_size, available_size);

                    /* setup new node */
                    if (opts.descriptive_node_names)
                    {
                        new_node.name += old_node.name;
                        if (taken_size != ssize(old_node.seq))
                            new_node.name += "_sub";
                    }

                    assert(old_node.regions.size() > 0);

                    if (bio::io::genomic_region & reg = old_node.regions.back(); reg.beg <= reg.end)
                    {
                        append_region(bio::io::genomic_region{reg.chrom,
                                                              reg.beg + old_node_offset,
                                                              reg.beg + old_node_offset + taken_size},
                                      new_node.regions);
                    }
                    else
                    {
                        assert(reg.beg - reg.end == ssize(old_node.seq));
                        append_region(bio::io::genomic_region{reg.chrom,
                                                              reg.beg - old_node_offset,
                                                              reg.beg - old_node_offset - taken_size},
                                      new_node.regions);
                    }

                    std::ranges::copy(old_seq | std::views::take(taken_size), std::back_inserter(new_node.seq));

                    if (taken_size == available_size) // need to go to next node
                    {
                        /* handle old node */
                        if (!old_node.is_ref)
                            old_node.tobe_deleted = true; // mark for deletion (only non-ref nodes)

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
                        if (opts.descriptive_node_names)
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
                graph.nodes.push_back(std::move(new_node));
            }

            // add in-arcs to path
            size_t const first_path_node_i = path.front();
            size_t const first_new_node_i  = graph.nodes.size() - n_new_nodes;
            for (node & n : graph.nodes)
                for (size_t const target_node_i : n.arcs)
                    if (target_node_i == first_path_node_i)
                        n.arcs.push_back(first_new_node_i), ({ break; }); // :-)

            // add out-arcs from path
            graph.nodes.back().arcs = graph.nodes[path.back()].arcs;
        }
    }

    fmt::print(stderr, "source_location: {}:{}\n", __FILE__, __LINE__);
    std::cerr << std::flush;

    /* Step 3:
     *
     * Remove the left-overs
     */
    erase_todo_nodes(graph);

    // print_all_paths(graph);

    /* Step 4:
     *
     * Merge adjacent small ref nodes
     *
     *          ---·NR1
     *         /
     * R1·---·R2---·R3
     *
     * Step 2 turns the previous into:
     *
     * R1·---------·NR1
     *
     * R1·---·R2---·R3
     *
     * Now we remove R2:
     *
     * R1·--------·NR1
     *
     * R1·--------·R3
     */
    in_degrees.clear();
    in_degrees.resize(graph.nodes.size());
    for (node & n : graph.nodes)
        for (size_t const i_target_node : n.arcs)
            ++in_degrees[i_target_node];

    auto can_consume_next = [&graph, &in_degrees, w](node const & n)
    {
        return n.arcs.size() == 1 &&                                                 // n's out-degree == 1
               in_degrees[n.arcs.front()] == 1 &&                                    // next's in-degree == 1
               graph.nodes[n.arcs.front()].is_ref &&                                 // target is also ref
               (n.regions.size() == 1 && n.regions.front().beg % w == 0) &&          // reg begins on block border
               (n.seq.size() + graph.nodes[n.arcs.front()].seq.size() <= 3 * w / 2); // merged size within limits
    };
    auto consume_next = [&graph](node & n)
    {
        node & next_n = graph.nodes[n.arcs.front()];
        n.name += "+";
        n.name += next_n.name;
        assign_append(next_n.seq, n.seq);
        n.arcs = std::move(next_n.arcs);
        append_regions(std::move(next_n.regions), n.regions);
        next_n.tobe_deleted = true;
    };

    // TODO do we need multiple passes of this on larger graphs?
    for (node & n : graph.nodes)
        if (n.is_ref && !n.tobe_deleted && can_consume_next(n))
            consume_next(n);

    erase_todo_nodes(graph);
}
