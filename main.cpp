#include "cut_mesh.h"
#include "half_edge.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include <algorithm>
#include <cmath>
#include <fmt/format.h>
#include <fmt/ostream.h>
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
    fmt::print(std::cerr, "GLFW Error {}: {}\n", error, description);
}

static void show_cut_mesh() {
    constexpr float canvas_width = 600.0f;
    const ImGuiViewport *main_viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(
        ImVec2(main_viewport->WorkPos.x + 20, main_viewport->WorkPos.y + 20),
        ImGuiCond_FirstUseEver);

    if (!ImGui::Begin("Cut-Cell")) {
        ImGui::End();
        return;
    }

    static std::vector<std::array<float, 2>> vertices;
    static std::vector<std::array<size_t, 2>> edges;
    static HalfEdgeMesh cut_mesh;
    static bool opt_enable_grid = true;
    static bool opt_construct_cut_mesh = false;
    static bool opt_draw_cut_vertices = true;
    static bool opt_draw_cut_edges = false;
    static bool opt_draw_original_lines = true;
    static bool adding_line = false;
    static auto selected_half_edge = cut_mesh.half_edges.end();

    ImGui::Checkbox("Enable grid", &opt_enable_grid);
    ImGui::Checkbox("Construct cut-mesh", &opt_construct_cut_mesh);
    ImGui::Checkbox("Draw cut-vertices", &opt_draw_cut_vertices);
    ImGui::Checkbox("Draw cut-edges", &opt_draw_cut_edges);
    ImGui::Checkbox("Draw original lines", &opt_draw_original_lines);
    if (ImGui::Button("Clear cut-mesh")) {
        cut_mesh = HalfEdgeMesh();
        selected_half_edge = cut_mesh.half_edges.end();
    }
    if (!opt_construct_cut_mesh &&
        selected_half_edge != cut_mesh.half_edges.end()) {
        ImGui::SameLine();
        if (ImGui::Button("Next")) {
            selected_half_edge = selected_half_edge->next;
        }
        ImGui::SameLine();
        if (ImGui::Button("Twin")) {
            selected_half_edge = selected_half_edge->twin;
        }
    }
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
        if (vertices.empty()) {
            vertices.emplace_back(
                std::array{mouse_pos_in_canvas.x / canvas_width,
                           mouse_pos_in_canvas.y / canvas_width});
        }
        vertices.emplace_back(std::array{mouse_pos_in_canvas.x / canvas_width,
                                         mouse_pos_in_canvas.y / canvas_width});
        adding_line = true;
    }
    if (adding_line) {
        vertices.back() = std::array{mouse_pos_in_canvas.x / canvas_width,
                                     mouse_pos_in_canvas.y / canvas_width};
        if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
            adding_line = false;
    }

    // Context menu (under default mouse threshold)
    ImGui::OpenPopupOnItemClick("context", ImGuiPopupFlags_MouseButtonRight);
    if (ImGui::BeginPopup("context")) {
        if (adding_line)
            vertices.resize(vertices.size() - 1);
        adding_line = false;
        if (ImGui::MenuItem("Remove one", nullptr, false, !vertices.empty())) {
            vertices.resize(vertices.size() - 1);
            if (vertices.size() == 1) {
                vertices.clear();
            }
        }
        if (ImGui::MenuItem("Remove all", nullptr, false, !vertices.empty())) {
            vertices.clear();
        }
        ImGui::EndPopup();
    }

    // Construct cut-mesh
    for (auto &v : vertices) {
        v[0] = std::clamp(v[0], 0.0f,
                          1.0f - std::numeric_limits<float>::epsilon());
        v[1] = std::clamp(v[1], 0.0f,
                          1.0f - std::numeric_limits<float>::epsilon());
    }
    if (opt_construct_cut_mesh) {
        edges.clear();
        for (size_t i = 0; i < vertices.size(); ++i) {
            auto j = (i + 1) % vertices.size();
            edges.emplace_back(std::array{i, j});
        }
        cut_mesh = construct_cut_mesh(vertices, edges);
        selected_half_edge = cut_mesh.half_edges.begin();
    }

    // Draw grid + all lines in the canvas
    draw_list->PushClipRect(canvas_p0, canvas_p1, true);
    if (opt_enable_grid) {
        const float grid_step = canvas_width / n_grid;
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
    if (opt_draw_original_lines) {
        for (size_t n = 0; n + 1 < vertices.size(); ++n) {
            draw_list->AddLine(
                ImVec2(origin.x + vertices[n][0] * canvas_width,
                       origin.y + vertices[n][1] * canvas_width),
                ImVec2(origin.x + vertices[n + 1][0] * canvas_width,
                       origin.y + vertices[n + 1][1] * canvas_width),
                IM_COL32(0, 0, 0, 255), 2.0f);
        }
        if (vertices.size() >= 3) {
            draw_list->AddLine(
                ImVec2(origin.x + vertices[0][0] * canvas_width,
                       origin.y + vertices[0][1] * canvas_width),
                ImVec2(
                    origin.x + vertices[vertices.size() - 1][0] * canvas_width,
                    origin.y + vertices[vertices.size() - 1][1] * canvas_width),
                IM_COL32(0, 0, 0, 255), 2.0f);
        }
    }
    if (opt_draw_cut_edges) {
        for (const auto &e : cut_mesh.edges) {
            const auto &v0 = e.half_edge->vertex;
            const auto &v1 = e.half_edge->twin->vertex;
            /*if (v0->id < n_grid_nodes && v1->id < n_grid_nodes) {
                continue;
            }*/
            auto pos0 = v0->position * canvas_width;
            auto pos1 = v1->position * canvas_width;
            ImVec2 p0(origin.x + pos0.x(), origin.y + pos0.y());
            ImVec2 p1(origin.x + pos1.x(), origin.y + pos1.y());
            draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255), 2);
        }
    }
    if (!opt_construct_cut_mesh &&
        selected_half_edge != cut_mesh.half_edges.end()) {
        std::vector<ImVec2> points;
        auto h = selected_half_edge;
        do {
            auto v = h->vertex;
            points.emplace_back(origin.x + v->position.x() * canvas_width,
                                origin.y + v->position.y() * canvas_width);
            h = h->next;
        } while (h != selected_half_edge);
        draw_list->AddPolyline(points.data(), static_cast<int>(points.size()),
                               IM_COL32(0, 255, 0, 255), ImDrawFlags_Closed,
                               3.0f);
        const auto &v0 = selected_half_edge->vertex;
        const auto &v1 = selected_half_edge->twin->vertex;
        draw_list->AddLine(ImVec2(origin.x + v0->position.x() * canvas_width,
                                  origin.y + v0->position.y() * canvas_width),
                           ImVec2(origin.x + v1->position.x() * canvas_width,
                                  origin.y + v1->position.y() * canvas_width),
                           IM_COL32(0, 0, 255, 255), 4.0f);
    }
    if (opt_draw_cut_vertices) {
        for (const auto &v : cut_mesh.vertices) {
            draw_list->AddCircleFilled(
                ImVec2(origin.x + v.position.x() * canvas_width,
                       origin.y + v.position.y() * canvas_width),
                4, IM_COL32(255, 0, 0, 255));
        }
    }
    draw_list->PopClipRect();
    ImGui::End();
}

int main() {
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
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        show_cut_mesh();

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
