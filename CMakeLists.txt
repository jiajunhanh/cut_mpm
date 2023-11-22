cmake_minimum_required(VERSION 3.25)
project(mandoline2d)

set(CMAKE_CXX_STANDARD 17)

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
        ${PROJECT_SOURCE_DIR}/deps/imgui/backends/imgui_impl_sdl2.cpp
        ${PROJECT_SOURCE_DIR}/deps/imgui/backends/imgui_impl_vulkan.cpp
)

target_include_directories(mandoline2d PRIVATE
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/deps/eigen
        ${PROJECT_SOURCE_DIR}/deps/imgui
        ${PROJECT_SOURCE_DIR}/deps/imgui/backends
)

find_package(Vulkan REQUIRED)
target_link_libraries(mandoline2d PRIVATE Vulkan::Vulkan)

add_subdirectory(deps/sdl)
target_link_libraries(mandoline2d PRIVATE SDL2::SDL2main)
target_link_libraries(mandoline2d PRIVATE SDL2::SDL2)
target_include_directories(mandoline2d PRIVATE ${PROJECT_SOURCE_DIR}/deps/sdl/include)

target_compile_options(mandoline2d PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W3 /WX /wd4819>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Werror>
        #$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
)
