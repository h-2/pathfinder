#include <string>

#include <sharg/all.hpp>

struct find_options;

void parse_find(sharg::parser & parser);
void find(find_options const & options);
