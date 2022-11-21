set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable (pathfinder main.cpp index.cpp find.cpp)
#set_target_properties(pathfinder PROPERTIES CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
target_link_libraries (pathfinder sharg::sharg fmt::fmt-header-only cereal::cereal biocpp::core biocpp::io)
