#include <filesystem>
#include <string>

#include <fmt/core.h>
#include <sharg/all.hpp>

#include "index.hpp"

struct index_options
{
    std::filesystem::path gfa;
    std::filesystem::path index;
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

    parser.parse();
    index(options);
}

void index(index_options const & options)
{

    //TODO

}
