#pragma once

#include <sharg/all.hpp>
#include <string>

struct index_options
{
    std::filesystem::path gfa;
    std::filesystem::path index;

    size_t window_size = 100;
    size_t k_mer_size  = 20;

    bool descriptive_node_names = false;
};

void parse_index(sharg::parser & parser);
void index(index_options const & options);
