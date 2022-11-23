#include <filesystem>
#include <string>

#include <fmt/core.h>
#include <sharg/all.hpp>

#include "index.hpp"
#include "graph.hpp"

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
                                               .validator = sharg::input_file_validator{{".gfa", ".gfa.gz"}}});

    parser.add_positional_option(options.index,
                                 sharg::config{.description = "The index file to create.",
                                               .validator = sharg::output_file_validator{{".pfi", ".pfi.gz"}}});

    parser.add_option(options.window_size, sharg::config{.short_id = 'w',
                                                         .long_id = "window-size",
                                                         .description = "The size of the genome buckets.",
                                                         .validator = sharg::arithmetic_range_validator{20, 2000}});
    parser.add_option(options.k_mer_size, sharg::config{.short_id = 'k',
                                                        .long_id = "k-mer-size",
                                                        .description = "The size of the k-mer",
                                                        .validator = sharg::arithmetic_range_validator{10, 30}});


    parser.parse();
    index(options);
}

void index(index_options const & options)
{

    //TODO

    graph prelim_graph;
    read_gfa(options.gfa, prelim_graph);

    fmt::print("stable_seq_names: {}\n", prelim_graph.stable_sequence_names);
    fmt::print("seqs: {}\n",             prelim_graph.seqs);
}
