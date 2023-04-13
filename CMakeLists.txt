cmake_minimum_required(VERSION 3.25)
project(mandoline2d)

set(CMAKE_CXX_STANDARD 20)

add_executable(mandoline2d main.cpp half_edge.cpp half_edge.h)

target_compile_options(mandoline2d PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W3 /WX>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Werror>
        )
