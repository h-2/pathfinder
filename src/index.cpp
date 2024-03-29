#include "index.hpp"

#include <filesystem>
#include <sharg/all.hpp>
#include <src/pf_graph.hpp>
#include <string>

#include <fmt/core.h>

#include "gfa.hpp"

void parse_index(sharg::parser & parser)
{
    index_options options{};

    parser.add_positional_option(options.gfa,
                                 sharg::config{.description = "The input graph in GFA format.",
                                               .validator   = sharg::input_file_validator{{".gfa", ".gfa.gz"}}});

    parser.add_positional_option(options.index,
                                 sharg::config{.description = "The index file to create.",
                                               .validator   = sharg::output_file_validator{{".pfi", ".pfi.gz"}}});

    parser.add_option(options.window_size,
                      sharg::config{
                        .short_id    = 'w',
                        .long_id     = "window-size",
                        .description = "The size of the genome buckets.",
                        .validator   = sharg::arithmetic_range_validator{2, 2000}
    });
    parser.add_option(options.k_mer_size,
                      sharg::config{
                        .short_id    = 'k',
                        .long_id     = "k-mer-size",
                        .description = "The size of the k-mer",
                        .validator   = sharg::arithmetic_range_validator{10, 30}
    });

    parser.add_option(options.descriptive_node_names,
                      sharg::config{.long_id     = "descriptive-node-names",
                                    .description = "Whether to include original name in new node names."});

    parser.parse();
    index(options);
}

void index(index_options const & options)
{
    //TODO

    pf_graph graph;

    {
        gfa::gfa_graph prelim_graph;
        gfa::read_graph(options.gfa, prelim_graph);
        gfa::gfa_graph2pf_graph(std::move(prelim_graph), graph);
    }

    // fmt::print("stable_seq_names: {}\n", prelim_graph.stable_sequence_names);
    // fmt::print("seqs: {}\n",             prelim_graph.seqs);
    // fmt::print("#arcs: {}\n",            prelim_graph.arcs.concat_size());

    // print_all_paths(graph);
    // graph_to_dot(graph);

    print_stats(graph);

    auto before_paths = all_paths_as_collapsed_regions(graph);

    discretise(graph, options);

    fmt::print("\n=====\t=====\t======\t======\t=====\t====\t\n\n");

    auto after_paths = all_paths_as_collapsed_regions(graph);

    // print_all_paths(graph);
    // graph_to_dot(graph);
    print_stats(graph);

    print_outlier_nodes(graph, options);

    fmt::print("Sorting before_paths...");
    std::ranges::sort(before_paths);
    fmt::print("done.\n");

    fmt::print("Sorting after_paths...");
    std::ranges::sort(after_paths);
    fmt::print("done.\n");

    bool eq = std::ranges::equal(before_paths, after_paths);
    fmt::print("Paths are equal: {}\n", eq);
    if (!eq)
    {
        FILE * file = std::fopen("before.txt", "w");
        for (auto const & pth : before_paths)
            fmt::print(file, "{}\n", pth);
        std::fclose(file);

        file = std::fopen("after.txt", "w");
        for (auto const & pth : after_paths)
            fmt::print(file, "{}\n", pth);
        std::fclose(file);

        fmt::print("Paths before transformation were written to before.txt\n");
        fmt::print("Paths after transformation were written to after.txt\n");
        fmt::print("Please run a diff on them.\n");
    }
}
