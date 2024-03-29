# ===========================================================================
#                  Pathfinder
# ===========================================================================

cmake_minimum_required (VERSION 3.12.0)
string(ASCII 27 Esc)
set(ColourBold  "${Esc}[1m")
set(ColourReset "${Esc}[m")
set(ColourRed   "${Esc}[31m")

message ("${ColourBold}Compiler Detection${ColourReset}")

project (pathfinder CXX)

# ----------------------------------------------------------------------------
# Make "Release" the default cmake build type
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Warn if cmake build type is not "Release"
if (NOT CMAKE_BUILD_TYPE STREQUAL Release)
    message (STATUS "${ColourRed}CMAKE_BUILD_TYPE is not \"Release\", your binaries will be slow.${ColourReset}")
endif ()


# ----------------------------------------------------------------------------
# Begin of dependency detection
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Dependency detection${ColourReset}")

# ZLIB
find_package(ZLIB QUIET REQUIRED)

# CEREAL
set(BUILD_SANDBOX OFF)
set(BUILD_DOC OFF)
set(BUILD_TESTS OFF)
set(SKIP_PERFORMANCE_COMPARISON ON)
set(CEREAL_INSTALL OFF)
add_subdirectory(submodules/cereal)

# FMT
set(FMT_DOC OFF)
set(FMT_INSTALL OFF)
set(FMT_TEST OFF)
add_subdirectory(submodules/fmt)

# BIO
find_package (biocpp COMPONENTS core io QUIET REQUIRED HINTS "${CMAKE_SOURCE_DIR}/submodules/biocpp-core/build_system")

# SHARG
find_package (sharg QUIET REQUIRED HINTS "${CMAKE_SOURCE_DIR}/submodules/sharg-parser/build_system")

message(STATUS "These dependencies were found:")
message(   "     ZLIB             ${ZLIB_FOUND}      ${ZLIB_VERSION_STRING}")
message(   "     CEREAL           ${CEREAL_FOUND}      ${CEREAL_VERSION}")
message(   "     FMT              ${FMT_FOUND}      ${FMT_VERSION}")
message(   "     BIOCPP-CORE      ${BIOCPP_CORE_FOUND}      ${BIOCPP_CORE_VERSION}")
message(   "     BIOCPP-IO        ${BIOCPP_IO_FOUND}      ${BIOCPP_IO_VERSION}")
message(   "     SHARG            ${SHARG_FOUND}      ${SHARG_VERSION}")

# ----------------------------------------------------------------------------
# Add pathfinder targets
# ----------------------------------------------------------------------------

add_subdirectory(src)

# ----------------------------------------------------------------------------
# Add Tests
# ----------------------------------------------------------------------------

#message ("\n${ColourBold}Setting up unit tests${ColourReset}")

#enable_testing ()
#add_subdirectory (test EXCLUDE_FROM_ALL)

#--------------------------------------------------------------------------------------------------
# Clang format target
#--------------------------------------------------------------------------------------------------

find_program(CLANG_FORMAT "clang-format-14")

if (CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
    find_program(CLANG_FORMAT "clang-format")
endif()

if (NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
    add_custom_target (check_format "find" "${CMAKE_CURRENT_SOURCE_DIR}/src" "-name" "'*.[ch]pp'" "-exec"
                       ${CLANG_FORMAT} "-style=file" "-n" "-Werror" "{}" "+"
                       COMMENT "Checking the source-code with clang-format.")


    add_custom_target (reformat "find" "${CMAKE_CURRENT_SOURCE_DIR}/src" "-name" "'*.[ch]pp'" "-exec"
                       ${CLANG_FORMAT} "-style=file" "-i" "{}" "+"
                       COMMENT "Reformatting the source-code with clang-format.")
else ()
    message(STATUS "clang_format not found; not adding the 'check_format' and 'reformat' targets.")
endif()
