cmake_minimum_required(VERSION 3.25)
project(mandoline2d)

set(CMAKE_CXX_STANDARD 20)

add_executable(mandoline2d
        main.cpp
        half_edge.cpp
        deps/imgui/imgui.cpp
        deps/imgui/imgui_demo.cpp
        deps/imgui/imgui_draw.cpp
        deps/imgui/imgui_tables.cpp
        deps/imgui/imgui_widgets.cpp
        deps/imgui/backends/imgui_impl_glfw.cpp
        deps/imgui/backends/imgui_impl_opengl3.cpp
        deps/glad/src/gl.c
        )

target_include_directories(mandoline2d PRIVATE
        deps/imgui
        deps/imgui/backends
        deps/glad/include
        )

find_package(OpenGL REQUIRED)
target_link_libraries(mandoline2d OpenGL::GL)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

add_subdirectory(deps/glfw)
target_link_libraries(mandoline2d glfw)

add_subdirectory(deps/fmt)
target_link_libraries(mandoline2d fmt::fmt)

target_compile_options(mandoline2d PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W3 /WX>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Werror>
        )
