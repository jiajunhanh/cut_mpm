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

constexpr int n_grid = 8;
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

static auto lerp(Vertex &v, const Vertex &vi, const Vertex &vj, float t,
                 int d) {
    auto coord_i = vi.coord(d);
    auto coord_j = vj.coord(d);
    auto coord = coord_i + (coord_j - coord_i) * t;
    v.c[d] = static_cast<int>(std::floor(coord));
    v.q[d] = coord - std::floor(coord);
    v.b[d] = v.q[d] == 0.0f;
}

static void add_grid_edges(std::vector<Vertex> &cut_vertices,
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

static std::pair<std::vector<Vertex>, std::vector<Edge>>
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

static void show_cut_cell() {
    constexpr float canvas_width = 600.0f;
    const ImGuiViewport *main_viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(
        ImVec2(main_viewport->WorkPos.x + 20, main_viewport->WorkPos.y + 20),
        ImGuiCond_FirstUseEver);

    if (!ImGui::Begin("Cut-Cell")) {
        ImGui::End();
        return;
    }

    static std::vector<ImVec2> points;
    static bool opt_enable_grid = true;
    static bool opt_draw_cut_vertices = true;
    static bool opt_draw_cut_edges = true;
    static bool adding_line = false;

    ImGui::Checkbox("Enable grid", &opt_enable_grid);
    ImGui::Checkbox("Draw cut-vertices", &opt_draw_cut_vertices);
    ImGui::Checkbox("Draw cut-edges", &opt_draw_cut_edges);
    ImGui::Text(
        "Mouse Left: drag to add lines,\nMouse Right: click for context menu.");

    // Using InvisibleButton() as a convenience 1) it will advance the layout
    // cursor and 2) allows us to use IsItemHovered()/IsItemActive()
    ImVec2 canvas_p0 =
        ImGui::GetCursorScreenPos(); // ImDrawList API uses screen coordinates!
    ImVec2 canvas_sz(canvas_width, canvas_width);
    ImVec2 canvas_p1 =
        ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);

    // Draw border and background color
    ImGuiIO &io = ImGui::GetIO();
    ImDrawList *draw_list = ImGui::GetWindowDrawList();
    draw_list->AddRectFilled(canvas_p0, canvas_p1,
                             IM_COL32(200, 200, 200, 255));

    // This will catch our interactions
    ImGui::InvisibleButton("canvas", canvas_sz,
                           ImGuiButtonFlags_MouseButtonLeft |
                               ImGuiButtonFlags_MouseButtonRight);
    const bool is_hovered = ImGui::IsItemHovered(); // Hovered
    const ImVec2 origin(canvas_p0.x, canvas_p0.y);  // Lock scrolled origin
    const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x,
                                     io.MousePos.y - origin.y);

    if (is_hovered && !adding_line &&
        ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        if (points.empty()) {
            points.emplace_back(mouse_pos_in_canvas);
        }
        points.emplace_back(mouse_pos_in_canvas);
        adding_line = true;
    }
    if (adding_line) {
        points.back() = mouse_pos_in_canvas;
        if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
            adding_line = false;
    }

    // Context menu (under default mouse threshold)
    ImGui::OpenPopupOnItemClick("context", ImGuiPopupFlags_MouseButtonRight);
    if (ImGui::BeginPopup("context")) {
        if (adding_line)
            points.resize(points.size() - 1);
        adding_line = false;
        if (ImGui::MenuItem("Remove one", nullptr, false, !points.empty())) {
            points.resize(points.size() - 1);
            if (points.size() == 1) {
                points.clear();
            }
        }
        if (ImGui::MenuItem("Remove all", nullptr, false, !points.empty())) {
            points.clear();
        }
        ImGui::EndPopup();
    }

    // Construct cut-cell
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    for (size_t i = 0; i < points.size(); ++i) {
        Vertex v{};
        for (int d = 0; d < 2; ++d) {
            points[i][d] = std::clamp(points[i][d], 0.0f, canvas_width - 0.1f);
            auto x = points[i][d] / canvas_width * static_cast<float>(n_grid);
            v.c[d] = static_cast<int>(std::floor(x));
            v.q[d] = x - std::floor(x);
            v.b[d] = (v.q[d] == 0.0f);
        }
        vertices.emplace_back(v);
        auto j = (i + 1) % points.size();
        edges.emplace_back(Edge{i, j});
    }
    auto [cut_vertices, cut_edges] =
        compute_cut_vertices_and_edges(vertices, edges);

    // Draw grid + all lines in the canvas
    draw_list->PushClipRect(canvas_p0, canvas_p1, true);
    if (opt_enable_grid) {
        const float grid_step = canvas_width / static_cast<float>(n_grid);
        for (int i = 1; i < n_grid; ++i) {
            draw_list->AddLine(
                ImVec2(canvas_p0.x + grid_step * static_cast<float>(i),
                       canvas_p0.y),
                ImVec2(canvas_p0.x + grid_step * static_cast<float>(i),
                       canvas_p1.y),
                IM_COL32(0, 0, 0, 40));
            draw_list->AddLine(
                ImVec2(canvas_p0.x,
                       canvas_p0.y + grid_step * static_cast<float>(i)),
                ImVec2(canvas_p1.x,
                       canvas_p0.y + grid_step * static_cast<float>(i)),
                IM_COL32(0, 0, 0, 40));
        }
    }
    for (size_t n = 0; n + 1 < points.size(); ++n) {
        draw_list->AddLine(
            ImVec2(origin.x + points[n].x, origin.y + points[n].y),
            ImVec2(origin.x + points[n + 1].x, origin.y + points[n + 1].y),
            IM_COL32(0, 0, 0, 255), 2.0f);
    }
    if (points.size() >= 3) {
        draw_list->AddLine(
            ImVec2(origin.x + points[0].x, origin.y + points[0].y),
            ImVec2(origin.x + points[points.size() - 1].x,
                   origin.y + points[points.size() - 1].y),
            IM_COL32(0, 0, 0, 255), 2.0f);
    }
    if (opt_draw_cut_edges) {
        for (const auto &e : cut_edges) {
            if (e.i < n_grid_nodes && e.j < n_grid_nodes) {
                continue;
            }
            const auto &v0 = cut_vertices[e.i];
            const auto &v1 = cut_vertices[e.j];
            auto x0 = (static_cast<float>(v0.c[0]) + v0.q[0]) /
                      static_cast<float>(n_grid) * canvas_width;
            auto y0 = (static_cast<float>(v0.c[1]) + v0.q[1]) /
                      static_cast<float>(n_grid) * canvas_width;
            auto x1 = (static_cast<float>(v1.c[0]) + v1.q[0]) /
                      static_cast<float>(n_grid) * canvas_width;
            auto y1 = (static_cast<float>(v1.c[1]) + v1.q[1]) /
                      static_cast<float>(n_grid) * canvas_width;
            ImVec2 p0(origin.x + x0, origin.y + y0);
            ImVec2 p1(origin.x + x1, origin.y + y1);
            draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255), 2);
        }
    }
    if (opt_draw_cut_vertices) {
        for (size_t i = n_grid_nodes; i < cut_vertices.size(); ++i) {
            const auto &v = cut_vertices[i];
            auto x = (static_cast<float>(v.c[0]) + v.q[0]) /
                     static_cast<float>(n_grid) * canvas_width;
            auto y = (static_cast<float>(v.c[1]) + v.q[1]) /
                     static_cast<float>(n_grid) * canvas_width;
            draw_list->AddCircleFilled(ImVec2(origin.x + x, origin.y + y), 4,
                                       IM_COL32(255, 0, 0, 255));
        }
    }
    draw_list->PopClipRect();
    ImGui::End();
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
    // io.IniFilename = nullptr;

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    // ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Our state
    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        show_cut_cell();

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
    return 0;
}
