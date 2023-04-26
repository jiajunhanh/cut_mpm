#include "half_edge.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fstream>
#include <iostream>
#include <vector>
#define GL_SILENCE_DEPRECATION
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <GLES2/gl2.h>
#endif
#if defined(_MSC_VER) && (_MSC_VER >= 1900) &&                                 \
    !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif

static void glfw_error_callback(int error, const char *description) {
    // fprintf(stderr, "GLFW Error %d: %s\n", error, description);
    fmt::print(std::cerr, "GLFW Error {}: {}\n", error, description);
}

constexpr int n_grid = 3;
constexpr int n_grid_nodes = (n_grid + 1) * (n_grid + 1);
constexpr float dx = 1.0f / n_grid;

class Vertex {
  public:
    std::array<int, 2> c{};
    std::array<float, 2> q{};
    std::array<bool, 2> b{};

    std::vector<size_t> edges;

    [[nodiscard]] float coord(int d) const {
        return static_cast<float>(c[d]) + q[d];
    }

    [[nodiscard]] int node_id() const { return c[1] * (n_grid + 1) + c[0]; }
};

class Edge {
  public:
    size_t i{}; // point id
    size_t j{}; // point id
};

auto lerp(Vertex &v, const Vertex &vi, const Vertex &vj, float t, int d) {
    auto coord_i = vi.coord(d);
    auto coord_j = vj.coord(d);
    auto coord = coord_i + (coord_j - coord_i) * t;
    v.c[d] = static_cast<int>(std::floor(coord));
    v.q[d] = coord - std::floor(coord);
    v.b[d] = v.q[d] == 0.0f;
}

void add_grid_edges(std::vector<Vertex> &cut_vertices,
                    std::vector<Edge> &cut_edges) {
    std::vector<size_t> grid_cut_vertices[n_grid * n_grid][2];
    for (size_t i = n_grid_nodes; i < cut_vertices.size(); ++i) {
        const auto &v = cut_vertices[i];
        for (int d = 0; d < 2; ++d) {
            if (v.b[d]) {
                grid_cut_vertices[v.c[1] * n_grid + v.c[0]][d ^ 1].emplace_back(
                    i);
            }
        }
    }
    for (int y = 0; y < n_grid; ++y) {
        for (int x = 0; x < n_grid; ++x) {
            auto node_i = y * (n_grid + 1) + x;
            for (int d = 0; d < 2; ++d) {
                auto &v = grid_cut_vertices[y * n_grid + x][d];
                int node_j = node_i;
                if (d == 0) {
                    node_j += 1;
                } else {
                    node_j += n_grid + 1;
                }
                if (v.empty()) {
                    cut_edges.emplace_back(node_i, node_j);
                    continue;
                }
                std::sort(v.begin(), v.end(), [&](size_t i, size_t j) {
                    return cut_vertices[i].q[d] < cut_vertices[j].q[d];
                });
                cut_edges.emplace_back(node_i, v.front());
                for (size_t i = 0; i < v.size() - 1; ++i) {
                    cut_edges.emplace_back(v[i], v[i + 1]);
                }
                cut_edges.emplace_back(v.back(), node_j);
            }
        }
    }
    for (int i = 0; i < n_grid; ++i) {
        cut_edges.emplace_back(n_grid + (n_grid + 1) * i,
                               n_grid + (n_grid + 1) * (i + 1));
        cut_edges.emplace_back(n_grid * (n_grid + 1) + i,
                               n_grid * (n_grid + 1) + i + 1);
    }
    for (auto v : cut_vertices) {
        v.edges.clear();
    }
    for (size_t i = 0; i < cut_edges.size(); ++i) {
        const auto &e = cut_edges[i];
        cut_vertices[e.i].edges.push_back(i);
        cut_vertices[e.j].edges.push_back(i);
    }
}

std::pair<std::vector<Vertex>, std::vector<Edge>>
compute_cut_vertices_and_edges(const std::vector<Vertex> &vertices,
                               const std::vector<Edge> &edges) {
    std::vector<Vertex> cut_vertices;
    std::vector<Edge> cut_edges;
    for (int y = 0; y <= n_grid; ++y) {
        for (int x = 0; x <= n_grid; ++x) {
            Vertex v{};
            v.c[0] = x;
            v.c[1] = y;
            v.b[0] = v.b[1] = true;
            cut_vertices.emplace_back(v);
        }
    }
    cut_vertices.insert(cut_vertices.end(), vertices.cbegin(), vertices.cend());
    for (auto &e : edges) {
        std::vector<std::pair<float, size_t>> intersection_points;
        const auto &vi = vertices[e.i];
        const auto &vj = vertices[e.j];
        for (int d = 0; d < 2; ++d) {
            int z_min;
            int z_max;
            if (vi.c[d] < vj.c[d] ||
                (vi.c[d] == vj.c[d] && vi.q[d] < vj.q[d])) {
                z_min = vi.b[d] ? vi.c[d] : vi.c[d] + 1;
                z_max = vj.c[d];
            } else {
                z_min = vj.b[d] ? vj.c[d] : vj.c[d] + 1;
                z_max = vi.c[d];
            }
            for (int z = z_min; z <= z_max; ++z) {
                auto t = (static_cast<float>(z - vi.c[d]) - vi.q[d]) /
                         (static_cast<float>(vj.c[d] - vi.c[d]) +
                          (vj.q[d] - vi.q[d]));
                Vertex cut_v{};
                lerp(cut_v, vi, vj, t, d ^ 1);
                cut_v.c[d] = z;
                cut_v.q[d] = 0;
                cut_v.b[d] = true;
                if (cut_v.b[d ^ 1]) {
                    intersection_points.emplace_back(t, cut_v.node_id());
                    continue;
                }
                cut_vertices.emplace_back(cut_v);
                intersection_points.emplace_back(t, cut_vertices.size() - 1);
            }
        }
        std::sort(intersection_points.begin(), intersection_points.end(),
                  [](const std::pair<float, size_t> &a,
                     const std::pair<float, size_t> &b) {
                      return a.first < b.first;
                  });
        if (intersection_points.empty()) {
            cut_edges.emplace_back(e.i + n_grid_nodes, e.j + n_grid_nodes);
            continue;
        }
        cut_edges.emplace_back(e.i + n_grid_nodes,
                               intersection_points.front().second);
        for (size_t i = 0; i + 1 < intersection_points.size(); ++i) {
            if (intersection_points[i].second ==
                intersection_points[i + 1].second) {
                continue;
            }
            cut_edges.emplace_back(intersection_points[i].second,
                                   intersection_points[i + 1].second);
        }
        cut_edges.emplace_back(intersection_points.back().second,
                               e.j + n_grid_nodes);
    }
    add_grid_edges(cut_vertices, cut_edges);
    return {cut_vertices, cut_edges};
}

int main(int argc, char *argv[]) {
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;

        // Decide GL+GLSL versions
#if defined(IMGUI_IMPL_OPENGL_ES2)
    // GL ES 2.0 + GLSL 100
    const char *glsl_version = "#version 100";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(__APPLE__)
    // GL 3.2 + GLSL 150
    const char *glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);           // Required on Mac
#else
    // GL 3.0 + GLSL 130
    const char *glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+
    // only glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // 3.0+ only
#endif

    // Create window with graphics context
    GLFWwindow *window = glfwCreateWindow(
        1280, 720, "Dear ImGui GLFW+OpenGL3 example", nullptr, nullptr);
    if (window == nullptr)
        return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    (void)io;
    io.ConfigFlags |=
        ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
    io.ConfigFlags |=
        ImGuiConfigFlags_NavEnableGamepad; // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    // ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can
    // also load multiple fonts and use ImGui::PushFont()/PopFont() to select
    // them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you
    // need to select the font among multiple.
    // - If the file cannot be loaded, the function will return a nullptr.
    // Please handle those errors in your application (e.g. use an assertion, or
    // display an error and quit).
    // - The fonts will be rasterized at a given size (w/ oversampling) and
    // stored into a texture when calling
    // ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame
    // below will call.
    // - Use '#define IMGUI_ENABLE_FREETYPE' in your imconfig file to use
    // Freetype for higher quality font rendering.
    // - Read 'docs/FONTS.md' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string
    // literal you need to write a double backslash \\ !
    // - Our Emscripten build process allows embedding fonts to be accessible at
    // runtime from the "fonts/" folder. See Makefile.emscripten for details.
    // io.Fonts->AddFontDefault();
    // io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf", 18.0f);
    // io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    // io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    // io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    // ImFont* font =
    // io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f,
    // nullptr, io.Fonts->GetGlyphRangesJapanese()); IM_ASSERT(font != nullptr);

    // Our state
    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        // Poll and handle events (inputs, window resize, etc.)
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to
        // tell if dear imgui wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data
        // to your main application, or clear/overwrite your copy of the mouse
        // data.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input
        // data to your main application, or clear/overwrite your copy of the
        // keyboard data. Generally you may always pass all inputs to dear
        // imgui, and hide them from your application based on those two flags.
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // 1. Show the big demo window (Most of the sample code is in
        // ImGui::ShowDemoWindow()! You can browse its code to learn more about
        // Dear ImGui!).
        if (show_demo_window)
            ImGui::ShowDemoWindow(&show_demo_window);

        // 2. Show a simple window that we create ourselves. We use a Begin/End
        // pair to create a named window.
        {
            static float f = 0.0f;
            static int counter = 0;

            ImGui::Begin("Hello, world!"); // Create a window called "Hello,
                                           // world!" and append into it.

            ImGui::Text(
                "This is some useful text."); // Display some text (you can use
                                              // a format strings too)
            ImGui::Checkbox("Demo Window",
                            &show_demo_window); // Edit bools storing our window
                                                // open/close state
            ImGui::Checkbox("Another Window", &show_another_window);

            ImGui::SliderFloat(
                "float", &f, 0.0f,
                1.0f); // Edit 1 float using a slider from 0.0f to 1.0f
            ImGui::ColorEdit3(
                "clear color",
                (float *)&clear_color); // Edit 3 floats representing a color

            if (ImGui::Button(
                    "Button")) // Buttons return true when clicked (most widgets
                               // return true when edited/activated)
                counter++;
            ImGui::SameLine();
            ImGui::Text("counter = %d", counter);

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                        1000.0f / io.Framerate, io.Framerate);
            ImGui::End();
        }

        // 3. Show another simple window.
        if (show_another_window) {
            ImGui::Begin(
                "Another Window",
                &show_another_window); // Pass a pointer to our bool variable
                                       // (the window will have a closing button
                                       // that will clear the bool when clicked)
            ImGui::Text("Hello from another window!");
            if (ImGui::Button("Close Me"))
                show_another_window = false;
            ImGui::End();
        }

        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x * clear_color.w,
                     clear_color.y * clear_color.w,
                     clear_color.z * clear_color.w, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    if (argc < 2) {
        return 0;
    }
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::ifstream ifs(argv[1]);
    size_t n, m;
    ifs >> n >> m;
    while (n--) {
        float x[2]{};
        Vertex v{};
        for (int d = 0; d < 2; ++d) {
            ifs >> x[d];
            assert(x[d] >= 0.0f && x[d] < 1.0f);
            x[d] *= n_grid;
            v.c[d] = static_cast<int>(std::floor(x[d]));
            v.q[d] = x[d] - std::floor(x[d]);
            v.b[d] = (v.q[d] == 0.0f);
            // std::cout << v.c[d] << ' ' << v.q[d] << ' ' << v.b[d] << '\n';
        }
        vertices.emplace_back(v);
    }
    while (m--) {
        size_t i, j;
        ifs >> i >> j;
        edges.emplace_back(i, j);
    }
    auto [cut_vertices, cut_edges] =
        compute_cut_vertices_and_edges(vertices, edges);
    std::cout << cut_vertices.size() << ' ' << cut_edges.size() << '\n';
    for (size_t i = 0; i < cut_vertices.size(); ++i) {
        const auto &v = cut_vertices[i];
        std::cout << i << ":\n";
        for (int d = 0; d < 2; ++d) {
            std::cout << v.c[d] << ' ' << v.q[d] << ' ' << v.b[d] << '\n';
        }
    }
    for (const auto &e : cut_edges) {
        std::cout << e.i << ' ' << e.j << '\n';
    }
    return 0;
}
