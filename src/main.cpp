#include <string_view>

#include <fmt/core.h>
#include <sharg/all.hpp>

#include "find.hpp"
#include "index.hpp"

int main(int argc, char const ** argv)
{
    sharg::parser top_level_parser{"pathfinder", argc, argv, sharg::update_notifications::off, {"index", "find"}};

    top_level_parser.info.description.push_back("Choose one of the subcommands.");

    try
    {
        top_level_parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        fmt::print("[Pathfinder ERROR] {}\n", ext.what());
        return -1;
    }

    sharg::parser & sub_parser = top_level_parser.get_sub_parser(); // hold a reference to the sub_parser


    try
    {
        if (sub_parser.info.app_name == std::string_view{"index"})
            parse_index(sub_parser);
//         else if (sub_parser.info.app_name == std::string_view{"find"})
//             parse_find(sub_parser);
        else
            __builtin_unreachable();
    }
    catch (std::exception const & ext)
    {
        fmt::print("[Pathfinder ERROR] {}\n", ext.what());
        return -1;
    }

    return 0;
}
