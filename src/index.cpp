#include "index.hpp"

#include <filesystem>
#include <sharg/all.hpp>
#include <src/pf_graph.hpp>
#include <string>

#include <fmt/core.h>

#include "gfa.hpp"

struct index_options
{
    std::filesystem::path gfa;
    std::filesystem::path index;

    size_t window_size = 100;
    size_t k_mer_size  = 20;
};

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


    discretise(graph, options.window_size);

    fmt::print("\n=====\t=====\t======\t======\t=====\t====\t\n\n");

    // print_all_paths(graph);
    // graph_to_dot(graph);
    print_stats(graph);
}
