cmake_minimum_required(VERSION 3.25)
project(mandoline2d)

set(CMAKE_CXX_STANDARD 20)

add_executable(mandoline2d
        ${PROJECT_SOURCE_DIR}/main.cpp
        ${PROJECT_SOURCE_DIR}/half_edge.cpp
        ${PROJECT_SOURCE_DIR}/cut_mesh.cpp
        ${PROJECT_SOURCE_DIR}/mpm.cpp
        ${PROJECT_SOURCE_DIR}/deps/imgui/imgui.cpp
        ${PROJECT_SOURCE_DIR}/deps/imgui/imgui_demo.cpp
        ${PROJECT_SOURCE_DIR}/deps/imgui/imgui_draw.cpp
        ${PROJECT_SOURCE_DIR}/deps/imgui/imgui_tables.cpp
        ${PROJECT_SOURCE_DIR}/deps/imgui/imgui_widgets.cpp
        ${PROJECT_SOURCE_DIR}/deps/imgui/backends/imgui_impl_glfw.cpp
        ${PROJECT_SOURCE_DIR}/deps/imgui/backends/imgui_impl_opengl3.cpp
        ${PROJECT_SOURCE_DIR}/deps/glad/src/gl.c
        )

target_include_directories(mandoline2d PRIVATE
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/deps/eigen
        ${PROJECT_SOURCE_DIR}/deps/imgui
        ${PROJECT_SOURCE_DIR}/deps/imgui/backends
        ${PROJECT_SOURCE_DIR}/deps/glad/include
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
        $<$<CXX_COMPILER_ID:MSVC>:/W3 /WX /wd4819>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Werror>
        )
