#include <SDL.h>
#include <SDL_vulkan.h>
#include <vulkan/vulkan.h>

#include <array>
#include <cmath>
#include <cstdio>   // printf, fprintf
#include <cstdlib>  // abort
#include <memory>
#include <string>
#include <vector>

#include "cut_mesh.h"
#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_vulkan.h"
#include "mpm.h"
// #include <vulkan/vulkan_beta.h>

// #define IMGUI_UNLIMITED_FRAME_RATE
#ifdef _DEBUG
#define IMGUI_VULKAN_DEBUG_REPORT
#endif

// Data
static VkAllocationCallbacks* g_Allocator = nullptr;
static VkInstance g_Instance = VK_NULL_HANDLE;
static VkPhysicalDevice g_PhysicalDevice = VK_NULL_HANDLE;
static VkDevice g_Device = VK_NULL_HANDLE;
static uint32_t g_QueueFamily = static_cast<uint32_t>(-1);
static VkQueue g_Queue = VK_NULL_HANDLE;
// static VkDebugReportCallbackEXT g_DebugReport = VK_NULL_HANDLE;
static VkPipelineCache g_PipelineCache = VK_NULL_HANDLE;
static VkDescriptorPool g_DescriptorPool = VK_NULL_HANDLE;

static ImGui_ImplVulkanH_Window g_MainWindowData;
static uint32_t g_MinImageCount = 2;
static bool g_SwapChainRebuild = false;

static void check_vk_result(VkResult err) {
    if (err == 0) return;
    fprintf(stderr, "[vulkan] Error: VkResult = %d\n", err);
    if (err < 0) abort();
}

#ifdef IMGUI_VULKAN_DEBUG_REPORT
static VKAPI_ATTR VkBool32 VKAPI_CALL
debug_report(VkDebugReportFlagsEXT flags, VkDebugReportObjectTypeEXT objectType,
             uint64_t object, size_t location, int32_t messageCode,
             const char* pLayerPrefix, const char* pMessage, void* pUserData) {
    (void)flags;
    (void)object;
    (void)location;
    (void)messageCode;
    (void)pUserData;
    (void)pLayerPrefix;  // Unused arguments
    fprintf(stderr,
            "[vulkan] Debug report from ObjectType: %i\nMessage: %s\n\n",
            objectType, pMessage);
    return VK_FALSE;
}
#endif  // IMGUI_VULKAN_DEBUG_REPORT

static bool IsExtensionAvailable(
    const ImVector<VkExtensionProperties>& properties, const char* extension) {
    for (const VkExtensionProperties& p : properties)
        if (strcmp(p.extensionName, extension) == 0) return true;
    return false;
}

static VkPhysicalDevice SetupVulkan_SelectPhysicalDevice() {
    uint32_t gpu_count;
    VkResult err = vkEnumeratePhysicalDevices(g_Instance, &gpu_count, nullptr);
    check_vk_result(err);
    IM_ASSERT(gpu_count > 0);

    ImVector<VkPhysicalDevice> gpus;
    gpus.resize(static_cast<int>(gpu_count));
    err = vkEnumeratePhysicalDevices(g_Instance, &gpu_count, gpus.Data);
    check_vk_result(err);

    // If a number >1 of GPUs got reported, find discrete GPU if present, or use
    // first one available. This covers most common cases
    // (multi-gpu/integrated+dedicated graphics). Handling more complicated
    // setups (multiple dedicated GPUs) is out of scope of this sample.
    for (VkPhysicalDevice& device : gpus) {
        VkPhysicalDeviceProperties properties;
        vkGetPhysicalDeviceProperties(device, &properties);
        if (properties.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU)
            return device;
    }

    // Use first GPU (Integrated) is a Discrete one is not available.
    if (gpu_count > 0) return gpus[0];
    return VK_NULL_HANDLE;
}

static void SetupVulkan(ImVector<const char*> instance_extensions) {
    VkResult err;

    // Create Vulkan Instance
    {
        VkInstanceCreateInfo create_info = {};
        create_info.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;

        // Enumerate available extensions
        uint32_t properties_count;
        ImVector<VkExtensionProperties> properties;
        vkEnumerateInstanceExtensionProperties(nullptr, &properties_count,
                                               nullptr);
        properties.resize(static_cast<int>(properties_count));
        err = vkEnumerateInstanceExtensionProperties(nullptr, &properties_count,
                                                     properties.Data);
        check_vk_result(err);

        // Enable required extensions
        if (IsExtensionAvailable(
                properties,
                VK_KHR_GET_PHYSICAL_DEVICE_PROPERTIES_2_EXTENSION_NAME))
            instance_extensions.push_back(
                VK_KHR_GET_PHYSICAL_DEVICE_PROPERTIES_2_EXTENSION_NAME);
#ifdef VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME
        if (IsExtensionAvailable(
                properties, VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME)) {
            instance_extensions.push_back(
                VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME);
            create_info.flags |=
                VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR;
        }
#endif

        // Enabling validation layers
#ifdef IMGUI_VULKAN_DEBUG_REPORT
        const char* layers[] = {"VK_LAYER_KHRONOS_validation"};
        create_info.enabledLayerCount = 1;
        create_info.ppEnabledLayerNames = layers;
        instance_extensions.push_back("VK_EXT_debug_report");
#endif

        // Create Vulkan Instance
        create_info.enabledExtensionCount =
            static_cast<uint32_t>(instance_extensions.Size);
        create_info.ppEnabledExtensionNames = instance_extensions.Data;
        err = vkCreateInstance(&create_info, g_Allocator, &g_Instance);
        check_vk_result(err);

        // Setup the debug report callback
#ifdef IMGUI_VULKAN_DEBUG_REPORT
        auto vkCreateDebugReportCallbackEXT =
            (PFN_vkCreateDebugReportCallbackEXT)vkGetInstanceProcAddr(
                g_Instance, "vkCreateDebugReportCallbackEXT");
        IM_ASSERT(vkCreateDebugReportCallbackEXT != nullptr);
        VkDebugReportCallbackCreateInfoEXT debug_report_ci = {};
        debug_report_ci.sType =
            VK_STRUCTURE_TYPE_DEBUG_REPORT_CALLBACK_CREATE_INFO_EXT;
        debug_report_ci.flags = VK_DEBUG_REPORT_ERROR_BIT_EXT |
                                VK_DEBUG_REPORT_WARNING_BIT_EXT |
                                VK_DEBUG_REPORT_PERFORMANCE_WARNING_BIT_EXT;
        debug_report_ci.pfnCallback = debug_report;
        debug_report_ci.pUserData = nullptr;
        err = vkCreateDebugReportCallbackEXT(g_Instance, &debug_report_ci,
                                             g_Allocator, &g_DebugReport);
        check_vk_result(err);
#endif
    }

    // Select Physical Device (GPU)
    g_PhysicalDevice = SetupVulkan_SelectPhysicalDevice();

    // Select graphics queue family
    {
        uint32_t count;
        vkGetPhysicalDeviceQueueFamilyProperties(g_PhysicalDevice, &count,
                                                 nullptr);
        auto queues = static_cast<VkQueueFamilyProperties*>(
            malloc(sizeof(VkQueueFamilyProperties) * count));
        vkGetPhysicalDeviceQueueFamilyProperties(g_PhysicalDevice, &count,
                                                 queues);
        for (uint32_t i = 0; i < count; i++)
            if (queues[i].queueFlags & VK_QUEUE_GRAPHICS_BIT) {
                g_QueueFamily = i;
                break;
            }
        free(queues);
        IM_ASSERT(g_QueueFamily != static_cast<uint32_t>(-1));
    }

    // Create Logical Device (with 1 queue)
    {
        ImVector<const char*> device_extensions;
        device_extensions.push_back("VK_KHR_swapchain");

        // Enumerate physical device extension
        uint32_t properties_count;
        ImVector<VkExtensionProperties> properties;
        vkEnumerateDeviceExtensionProperties(g_PhysicalDevice, nullptr,
                                             &properties_count, nullptr);
        properties.resize(static_cast<int>(properties_count));
        vkEnumerateDeviceExtensionProperties(
            g_PhysicalDevice, nullptr, &properties_count, properties.Data);
#ifdef VK_KHR_PORTABILITY_SUBSET_EXTENSION_NAME
        if (IsExtensionAvailable(properties,
                                 VK_KHR_PORTABILITY_SUBSET_EXTENSION_NAME))
            device_extensions.push_back(
                VK_KHR_PORTABILITY_SUBSET_EXTENSION_NAME);
#endif

        const float queue_priority[] = {1.0f};
        VkDeviceQueueCreateInfo queue_info[1] = {};
        queue_info[0].sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
        queue_info[0].queueFamilyIndex = g_QueueFamily;
        queue_info[0].queueCount = 1;
        queue_info[0].pQueuePriorities = queue_priority;
        VkDeviceCreateInfo create_info = {};
        create_info.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
        create_info.queueCreateInfoCount =
            sizeof(queue_info) / sizeof(queue_info[0]);
        create_info.pQueueCreateInfos = queue_info;
        create_info.enabledExtensionCount =
            static_cast<uint32_t>(device_extensions.Size);
        create_info.ppEnabledExtensionNames = device_extensions.Data;
        err = vkCreateDevice(g_PhysicalDevice, &create_info, g_Allocator,
                             &g_Device);
        check_vk_result(err);
        vkGetDeviceQueue(g_Device, g_QueueFamily, 0, &g_Queue);
    }

    // Create Descriptor Pool
    // The example only requires a single combined image sampler descriptor for
    // the font image and only uses one descriptor set (for that) If you wish to
    // load e.g. additional textures you may need to alter pools sizes.
    {
        VkDescriptorPoolSize pool_sizes[] = {
            {VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1},
        };
        VkDescriptorPoolCreateInfo pool_info = {};
        pool_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        pool_info.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
        pool_info.maxSets = 1;
        pool_info.poolSizeCount =
            static_cast<uint32_t>(IM_ARRAYSIZE(pool_sizes));
        pool_info.pPoolSizes = pool_sizes;
        err = vkCreateDescriptorPool(g_Device, &pool_info, g_Allocator,
                                     &g_DescriptorPool);
        check_vk_result(err);
    }
}

// All the ImGui_ImplVulkanH_XXX structures/functions are optional helpers used
// by the demo. Your real engine/app may not use them.
static void SetupVulkanWindow(ImGui_ImplVulkanH_Window* wd,
                              VkSurfaceKHR surface, int width, int height) {
    wd->Surface = surface;

    // Check for WSI support
    VkBool32 res;
    vkGetPhysicalDeviceSurfaceSupportKHR(g_PhysicalDevice, g_QueueFamily,
                                         wd->Surface, &res);
    if (res != VK_TRUE) {
        fprintf(stderr, "Error no WSI support on physical device 0\n");
        exit(-1);
    }

    // Select Surface Format
    const VkFormat requestSurfaceImageFormat[] = {
        VK_FORMAT_B8G8R8A8_UNORM, VK_FORMAT_R8G8B8A8_UNORM,
        VK_FORMAT_B8G8R8_UNORM, VK_FORMAT_R8G8B8_UNORM};
    const VkColorSpaceKHR requestSurfaceColorSpace =
        VK_COLORSPACE_SRGB_NONLINEAR_KHR;
    wd->SurfaceFormat = ImGui_ImplVulkanH_SelectSurfaceFormat(
        g_PhysicalDevice, wd->Surface, requestSurfaceImageFormat,
        static_cast<size_t>(IM_ARRAYSIZE(requestSurfaceImageFormat)),
        requestSurfaceColorSpace);

    // Select Present Mode
#ifdef IMGUI_UNLIMITED_FRAME_RATE
    VkPresentModeKHR present_modes[] = {VK_PRESENT_MODE_MAILBOX_KHR,
                                        VK_PRESENT_MODE_IMMEDIATE_KHR,
                                        VK_PRESENT_MODE_FIFO_KHR};
#else
    VkPresentModeKHR present_modes[] = {VK_PRESENT_MODE_FIFO_KHR};
#endif
    wd->PresentMode = ImGui_ImplVulkanH_SelectPresentMode(
        g_PhysicalDevice, wd->Surface, &present_modes[0],
        IM_ARRAYSIZE(present_modes));
    // printf("[vulkan] Selected PresentMode = %d\n", wd->PresentMode);

    // Create SwapChain, RenderPass, Framebuffer, etc.
    IM_ASSERT(g_MinImageCount >= 2);
    ImGui_ImplVulkanH_CreateOrResizeWindow(
        g_Instance, g_PhysicalDevice, g_Device, wd, g_QueueFamily, g_Allocator,
        width, height, g_MinImageCount);
}

static void CleanupVulkan() {
    vkDestroyDescriptorPool(g_Device, g_DescriptorPool, g_Allocator);

#ifdef IMGUI_VULKAN_DEBUG_REPORT
    // Remove the debug report callback
    auto vkDestroyDebugReportCallbackEXT =
        (PFN_vkDestroyDebugReportCallbackEXT)vkGetInstanceProcAddr(
            g_Instance, "vkDestroyDebugReportCallbackEXT");
    vkDestroyDebugReportCallbackEXT(g_Instance, g_DebugReport, g_Allocator);
#endif  // IMGUI_VULKAN_DEBUG_REPORT

    vkDestroyDevice(g_Device, g_Allocator);
    vkDestroyInstance(g_Instance, g_Allocator);
}

static void CleanupVulkanWindow() {
    ImGui_ImplVulkanH_DestroyWindow(g_Instance, g_Device, &g_MainWindowData,
                                    g_Allocator);
}

static void FrameRender(ImGui_ImplVulkanH_Window* wd, ImDrawData* draw_data) {
    VkResult err;

    VkSemaphore image_acquired_semaphore =
        wd->FrameSemaphores[wd->SemaphoreIndex].ImageAcquiredSemaphore;
    VkSemaphore render_complete_semaphore =
        wd->FrameSemaphores[wd->SemaphoreIndex].RenderCompleteSemaphore;
    err = vkAcquireNextImageKHR(g_Device, wd->Swapchain, UINT64_MAX,
                                image_acquired_semaphore, VK_NULL_HANDLE,
                                &wd->FrameIndex);
    if (err == VK_ERROR_OUT_OF_DATE_KHR || err == VK_SUBOPTIMAL_KHR) {
        g_SwapChainRebuild = true;
        return;
    }
    check_vk_result(err);

    ImGui_ImplVulkanH_Frame* fd = &wd->Frames[wd->FrameIndex];
    {
        err = vkWaitForFences(
            g_Device, 1, &fd->Fence, VK_TRUE,
            UINT64_MAX);  // wait indefinitely instead of periodically checking
        check_vk_result(err);

        err = vkResetFences(g_Device, 1, &fd->Fence);
        check_vk_result(err);
    }
    {
        err = vkResetCommandPool(g_Device, fd->CommandPool, 0);
        check_vk_result(err);
        VkCommandBufferBeginInfo info = {};
        info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        info.flags |= VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
        err = vkBeginCommandBuffer(fd->CommandBuffer, &info);
        check_vk_result(err);
    }
    {
        VkRenderPassBeginInfo info = {};
        info.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
        info.renderPass = wd->RenderPass;
        info.framebuffer = fd->Framebuffer;
        info.renderArea.extent.width = wd->Width;
        info.renderArea.extent.height = wd->Height;
        info.clearValueCount = 1;
        info.pClearValues = &wd->ClearValue;
        vkCmdBeginRenderPass(fd->CommandBuffer, &info,
                             VK_SUBPASS_CONTENTS_INLINE);
    }

    // Record dear imgui primitives into command buffer
    ImGui_ImplVulkan_RenderDrawData(draw_data, fd->CommandBuffer);

    // Submit command buffer
    vkCmdEndRenderPass(fd->CommandBuffer);
    {
        VkPipelineStageFlags wait_stage =
            VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        VkSubmitInfo info = {};
        info.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
        info.waitSemaphoreCount = 1;
        info.pWaitSemaphores = &image_acquired_semaphore;
        info.pWaitDstStageMask = &wait_stage;
        info.commandBufferCount = 1;
        info.pCommandBuffers = &fd->CommandBuffer;
        info.signalSemaphoreCount = 1;
        info.pSignalSemaphores = &render_complete_semaphore;

        err = vkEndCommandBuffer(fd->CommandBuffer);
        check_vk_result(err);
        err = vkQueueSubmit(g_Queue, 1, &info, fd->Fence);
        check_vk_result(err);
    }
}

static void FramePresent(ImGui_ImplVulkanH_Window* wd) {
    if (g_SwapChainRebuild) return;
    VkSemaphore render_complete_semaphore =
        wd->FrameSemaphores[wd->SemaphoreIndex].RenderCompleteSemaphore;
    VkPresentInfoKHR info = {};
    info.sType = VK_STRUCTURE_TYPE_PRESENT_INFO_KHR;
    info.waitSemaphoreCount = 1;
    info.pWaitSemaphores = &render_complete_semaphore;
    info.swapchainCount = 1;
    info.pSwapchains = &wd->Swapchain;
    info.pImageIndices = &wd->FrameIndex;
    VkResult err = vkQueuePresentKHR(g_Queue, &info);
    if (err == VK_ERROR_OUT_OF_DATE_KHR || err == VK_SUBOPTIMAL_KHR) {
        g_SwapChainRebuild = true;
        return;
    }
    check_vk_result(err);
    wd->SemaphoreIndex =
        (wd->SemaphoreIndex + 1) %
        wd->ImageCount;  // Now we can use the next set of semaphores
}

static void show_cut_mesh() {
    constexpr float canvas_width = 600.0f;
    const ImGuiViewport* main_viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(
        ImVec2(main_viewport->WorkPos.x + 20, main_viewport->WorkPos.y + 20),
        ImGuiCond_FirstUseEver);

    std::string title = "CUT_MPM (" + std::to_string(ImGui::GetIO().Framerate) +
                        " FPS)###CUT_MPM";
    if (!ImGui::Begin(title.c_str())) {
        ImGui::End();
        return;
    }

    constexpr Real gap = 0.0039;
    static std::vector<std::vector<std::array<Real, 2>>> scenarios{
        std::vector<std::array<Real, 2>>{},
        std::vector<std::array<Real, 2>>{{Real{0.3} - gap / 2, Real{0.6}},
                                         {Real{0.001}, Real{0.1}},
                                         {Real{0.001}, Real{0.999}},
                                         {Real{0.999}, Real{0.999}},
                                         {Real{0.999}, Real{0.1}},
                                         {Real{0.3} + gap / 2, Real{0.6}},
                                         {Real{0.95}, Real{0.95}},
                                         {Real{0.05}, Real{0.95}}},
        std::vector<std::array<Real, 2>>{{Real{0.375}, Real{0.59}},
                                         {Real{0.325}, Real{0.95}},
                                         {Real{0.425}, Real{0.95}}},
        std::vector<std::array<Real, 2>>{{Real{0.5}, Real{0.95}},
                                         {Real{0.6}, Real{0.5}},
                                         {Real{0.999}, Real{0.4}},
                                         {Real{0.999}, Real{0.999}},
                                         {Real{0.001}, Real{0.999}},
                                         {Real{0.001}, Real{0.4}},
                                         {Real{0.4}, Real{0.5}}},
        std::vector<std::array<Real, 2>>{{Real{0.075}, Real{0.4}},
                                         {Real{0.375}, Real{0.9}},
                                         {Real{0.675}, Real{0.4}},
                                         {Real{0.675}, Real{0.95}},
                                         {Real{0.075}, Real{0.95}}}};
    static int quality = 4;
    static int boundary = 0;
    static int material = 0;
    static std::vector<std::array<Real, 2>> vertices = scenarios[boundary];
    static auto cut_mesh = std::make_shared<CutMesh>(vertices, quality);
    static MPM mpm(cut_mesh, quality, material);
    static bool opt_draw_grid = true;
    static bool opt_construct_cut_mesh = true;
    // static bool opt_draw_cut_vertices = false;
    // static bool opt_draw_cut_edges = false;
    static bool opt_draw_original_lines = true;
    static bool opt_simulation = false;
    static bool adding_line = false;
    static bool simulating = false;
    static auto selected_half_edge = end(cut_mesh->half_edges());
    static auto selected_face = end(cut_mesh->faces());

    ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x * 0.25f);
    ImGui::SliderInt("Quality", &quality, 1, 8);
    ImGui::SameLine();
    const char* boundaries[] = {"None", "Narrow gap", "Cutting",
                                "Sharp angle (liquid)",
                                "Sharp angle (elastic)"};
    ImGui::Combo("Boundary", &boundary, boundaries, IM_ARRAYSIZE(boundaries));
    ImGui::PopItemWidth();
    ImGui::Checkbox("Draw grid", &opt_draw_grid);
    ImGui::SameLine();
    ImGui::Checkbox("Draw lines", &opt_draw_original_lines);
    ImGui::Checkbox("Construct cut-mesh", &opt_construct_cut_mesh);
    // ImGui::Checkbox("Draw cut-vertices", &opt_draw_cut_vertices);
    // ImGui::Checkbox("Draw cut-edges", &opt_draw_cut_edges);

    if (opt_simulation) opt_construct_cut_mesh = false;
    if (!opt_construct_cut_mesh &&
        selected_half_edge != end(cut_mesh->half_edges())) {
        ImGui::SameLine();
        if (ImGui::Button("Next")) {
            selected_half_edge = selected_half_edge->next;
        }
        ImGui::SameLine();
        if (ImGui::Button("Twin")) {
            selected_half_edge = selected_half_edge->twin;
        }
    }
    ImGui::Checkbox("Simulation", &opt_simulation);
    ImGui::SameLine();
    const char* materials[] = {"Liquid", "Elastic"};
    ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x * 0.25f);
    ImGui::Combo("Material", &material, materials, IM_ARRAYSIZE(materials));
    // if (opt_construct_cut_mesh) opt_simulation = false;
    /*if (ImGui::Button("Clear cut-mesh")) {
        cut_mesh = std::make_shared<CutMesh>();
        selected_half_edge = end(cut_mesh->half_edges());
        selected_face = end(cut_mesh->faces());
    }*/
    ImGui::Text(
        "Mouse Left: drag to add lines,\nMouse Right: click for context menu.");

    // Using InvisibleButton() as a convenience 1) it will advance the layout
    // cursor and 2) allows us to use IsItemHovered()/IsItemActive()
    ImVec2 canvas_p0 =
        ImGui::GetCursorScreenPos();  // ImDrawList API uses screen coordinates!
    ImVec2 canvas_sz(canvas_width, canvas_width);
    ImVec2 canvas_p1 =
        ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y);

    // Draw border and background color
    ImGuiIO& io = ImGui::GetIO();
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    draw_list->AddRectFilled(canvas_p0, canvas_p1,
                             IM_COL32(200, 200, 200, 255));

    // This will catch our interactions
    ImGui::InvisibleButton(
        "canvas", canvas_sz,
        ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
    const bool is_hovered = ImGui::IsItemHovered();  // Hovered
    const ImVec2 origin(canvas_p0.x, canvas_p0.y);   // Lock scrolled origin
    const ImVec2 mouse_pos_in_canvas(io.MousePos.x - origin.x,
                                     io.MousePos.y - origin.y);
    const std::array mouse_pos_in_grid{
        std::clamp(Real{mouse_pos_in_canvas.x} / canvas_width, Real{0.0},
                   Real{1.0} - std::numeric_limits<Real>::epsilon()),
        std::clamp(Real{mouse_pos_in_canvas.y} / canvas_width, Real{0.0},
                   Real{1.0} - std::numeric_limits<Real>::epsilon())};

    if (opt_construct_cut_mesh && is_hovered && !adding_line &&
        ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        if (vertices.empty()) {
            vertices.emplace_back(mouse_pos_in_grid);
        }
        boundary = 0;
        vertices.emplace_back(mouse_pos_in_grid);
        adding_line = true;
    }
    if (adding_line) {
        vertices.back() = mouse_pos_in_grid;
        if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) adding_line = false;
    }

    // Context menu (under default mouse threshold)
    ImGui::OpenPopupOnItemClick("context", ImGuiPopupFlags_MouseButtonRight);
    if (ImGui::BeginPopup("context")) {
        if (adding_line) vertices.resize(vertices.size() - 1);
        adding_line = false;
        if (ImGui::MenuItem("Remove one", nullptr, false, !vertices.empty())) {
            boundary = 0;
            vertices.resize(vertices.size() - 1);
            if (vertices.size() == 1) {
                vertices.clear();
            }
        }
        if (ImGui::MenuItem("Remove all", nullptr, false, !vertices.empty())) {
            boundary = 0;
            vertices.clear();
        }
        ImGui::EndPopup();
    }

    if (opt_construct_cut_mesh) {
        if (boundary) vertices = scenarios[boundary];
        *cut_mesh = CutMesh(vertices, quality);
        selected_half_edge = end(cut_mesh->half_edges());
        selected_face = end(cut_mesh->faces());
    }

    if (opt_simulation) {
        if (!simulating) {
            if (boundary) vertices = scenarios[boundary];
            *cut_mesh = CutMesh(vertices, quality);
            selected_half_edge = end(cut_mesh->half_edges());
            selected_face = end(cut_mesh->faces());
            mpm = MPM(cut_mesh, quality, material);
            mpm.initialize();
            simulating = true;
        }
        for (int i = 0; i < static_cast<int>(std::pow(2, quality - 1)) * 5 / 4;
             ++i) {
            mpm.update();
        }
        for (const auto& p : mpm.particles()) {
            float radius = quality > 4 ? 1.5 : 2;
            draw_list->AddCircleFilled(
                ImVec2(origin.x + static_cast<float>(p.x.x()) * canvas_width,
                       origin.y + static_cast<float>(p.x.y()) * canvas_width),
                radius, IM_COL32(6, 133, 135, 255));
        }
    } else {
        simulating = false;
    }

    // Draw grid + all lines in the canvas
    draw_list->PushClipRect(canvas_p0, canvas_p1, true);
    if (opt_draw_grid) {
        int grid_size = 8 * static_cast<int>(std::pow(2, quality - 1));
        const float grid_step = canvas_width / static_cast<float>(grid_size);
        for (int i = 1; i < grid_size; ++i) {
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
                ImVec2(static_cast<float>(origin.x +
                                          vertices[n][0] * canvas_width),
                       static_cast<float>(origin.y +
                                          vertices[n][1] * canvas_width)),
                ImVec2(static_cast<float>(origin.x +
                                          vertices[n + 1][0] * canvas_width),
                       static_cast<float>(origin.y +
                                          vertices[n + 1][1] * canvas_width)),
                IM_COL32(0, 0, 0, 255), 2.0);
        }
        if (vertices.size() >= 3) {
            draw_list->AddLine(
                ImVec2(static_cast<float>(origin.x +
                                          vertices[0][0] * canvas_width),
                       static_cast<float>(origin.y +
                                          vertices[0][1] * canvas_width)),
                ImVec2(static_cast<float>(origin.x +
                                          vertices[vertices.size() - 1][0] *
                                              canvas_width),
                       static_cast<float>(origin.y +
                                          vertices[vertices.size() - 1][1] *
                                              canvas_width)),
                IM_COL32(0, 0, 0, 255), 2.0);
        }
    }
    /*if (opt_draw_cut_edges) {
        for (const auto& e : cut_mesh->edges()) {
            const auto& v0 = e.half_edge->vertex;
            const auto& v1 = e.half_edge->twin->vertex;
            *//*if (v0->id < n_grid_nodes && v1->id < n_grid_nodes) {
                continue;
            }*//*
            auto pos0 = v0->position * canvas_width;
            auto pos1 = v1->position * canvas_width;
            ImVec2 p0(origin.x + static_cast<float>(pos0.x()),
                      origin.y + static_cast<float>(pos0.y()));
            ImVec2 p1(static_cast<float>(origin.x + pos1.x()),
                      static_cast<float>(origin.y + pos1.y()));
            draw_list->AddLine(p0, p1, IM_COL32(0, 0, 0, 255), 2);
        }
    }*/
    if (opt_simulation) {
        for (const auto& p : mpm.particles()) {
            float radius = quality > 4 ? 1.5 : 2;
            draw_list->AddCircleFilled(
                ImVec2(origin.x + static_cast<float>(p.x.x()) * canvas_width,
                       origin.y + static_cast<float>(p.x.y()) * canvas_width),
                radius, IM_COL32(6, 133, 135, 255));
        }
    }
    if (!opt_construct_cut_mesh) {
        if (selected_half_edge != end(cut_mesh->half_edges())) {
            std::vector<ImVec2> points;
            auto h = selected_half_edge;
            do {
                auto v = h->vertex;
                points.emplace_back(origin.x + v->position.x() * canvas_width,
                                    origin.y + v->position.y() * canvas_width);
                h = h->next;
            } while (h != selected_half_edge);
            draw_list->AddPolyline(
                points.data(), static_cast<int>(points.size()),
                IM_COL32(0, 255, 0, 255), ImDrawFlags_Closed, 2.0);
            const auto& v0 = selected_half_edge->vertex;
            const auto& v1 = selected_half_edge->twin->vertex;
            draw_list->AddLine(
                ImVec2(origin.x +
                           static_cast<float>(v0->position.x()) * canvas_width,
                       origin.y +
                           static_cast<float>(v0->position.y()) * canvas_width),
                ImVec2(origin.x +
                           static_cast<float>(v1->position.x()) * canvas_width,
                       origin.y +
                           static_cast<float>(v1->position.y()) * canvas_width),
                IM_COL32(0, 0, 255, 255), 2.0);
        }
        if (is_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            auto face = cut_mesh->get_enclosing_face(
                {mouse_pos_in_grid[0], mouse_pos_in_grid[1]});
            if (face == selected_half_edge->face) {
                selected_half_edge = end(cut_mesh->half_edges());
            } else {
                selected_half_edge = face->half_edge;
            }
        }
    }
    /*if (opt_draw_cut_vertices) {
        for (const auto& v : cut_mesh->vertices()) {
            draw_list->AddCircleFilled(
                ImVec2(static_cast<float>(origin.x + v.position.x()) *
                           canvas_width,
                       static_cast<float>(origin.y + v.position.y()) *
                           canvas_width),
                4, IM_COL32(255, 0, 0, 255));
        }
    }*/
    draw_list->PopClipRect();
    ImGui::End();
}

// Main code
int main(int, char**) {
    // Setup SDL
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) !=
        0) {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }

    // From 2.0.18: Enable native IME.
#ifdef SDL_HINT_IME_SHOW_UI
    SDL_SetHint(SDL_HINT_IME_SHOW_UI, "1");
#endif

    // Create window with Vulkan graphics context
    auto window_flags = static_cast<SDL_WindowFlags>(
        SDL_WINDOW_VULKAN | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_Window* window =
        SDL_CreateWindow("Dear ImGui SDL2+Vulkan", SDL_WINDOWPOS_CENTERED,
                         SDL_WINDOWPOS_CENTERED, 655, 800, window_flags);

    ImVector<const char*> extensions;
    uint32_t extensions_count = 0;
    SDL_Vulkan_GetInstanceExtensions(window, &extensions_count, nullptr);
    extensions.resize(static_cast<int>(extensions_count));
    SDL_Vulkan_GetInstanceExtensions(window, &extensions_count,
                                     extensions.Data);
    SetupVulkan(extensions);

    // Create Window Surface
    VkSurfaceKHR surface;
    VkResult err;
    if (SDL_Vulkan_CreateSurface(window, g_Instance, &surface) == 0) {
        printf("Failed to create Vulkan surface.\n");
        return 1;
    }

    // Create Framebuffers
    int w, h;
    SDL_GetWindowSize(window, &w, &h);
    ImGui_ImplVulkanH_Window* wd = &g_MainWindowData;
    SetupVulkanWindow(wd, surface, w, h);

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;
    io.ConfigFlags |=
        ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls
    io.ConfigFlags |=
        ImGuiConfigFlags_NavEnableGamepad;  // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    // ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplSDL2_InitForVulkan(window);
    ImGui_ImplVulkan_InitInfo init_info = {};
    init_info.Instance = g_Instance;
    init_info.PhysicalDevice = g_PhysicalDevice;
    init_info.Device = g_Device;
    init_info.QueueFamily = g_QueueFamily;
    init_info.Queue = g_Queue;
    init_info.PipelineCache = g_PipelineCache;
    init_info.DescriptorPool = g_DescriptorPool;
    init_info.Subpass = 0;
    init_info.MinImageCount = g_MinImageCount;
    init_info.ImageCount = wd->ImageCount;
    init_info.MSAASamples = VK_SAMPLE_COUNT_1_BIT;
    init_info.Allocator = g_Allocator;
    init_info.CheckVkResultFn = check_vk_result;
    ImGui_ImplVulkan_Init(&init_info, wd->RenderPass);

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
    // io.Fonts->AddFontDefault();
    // io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf", 18.0f);
    // io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    // io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    // io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    // ImFont* font =
    // io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f,
    // nullptr, io.Fonts->GetGlyphRangesJapanese()); IM_ASSERT(font != nullptr);

    // Upload Fonts
    {
        // Use any command queue
        VkCommandPool command_pool = wd->Frames[wd->FrameIndex].CommandPool;
        VkCommandBuffer command_buffer =
            wd->Frames[wd->FrameIndex].CommandBuffer;

        err = vkResetCommandPool(g_Device, command_pool, 0);
        check_vk_result(err);
        VkCommandBufferBeginInfo begin_info = {};
        begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        begin_info.flags |= VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
        err = vkBeginCommandBuffer(command_buffer, &begin_info);
        check_vk_result(err);

        ImGui_ImplVulkan_CreateFontsTexture(command_buffer);

        VkSubmitInfo end_info = {};
        end_info.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
        end_info.commandBufferCount = 1;
        end_info.pCommandBuffers = &command_buffer;
        err = vkEndCommandBuffer(command_buffer);
        check_vk_result(err);
        err = vkQueueSubmit(g_Queue, 1, &end_info, VK_NULL_HANDLE);
        check_vk_result(err);

        err = vkDeviceWaitIdle(g_Device);
        check_vk_result(err);
        ImGui_ImplVulkan_DestroyFontUploadObjects();
    }

    // Our state
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // Main loop
    bool done = false;
    while (!done) {
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
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            ImGui_ImplSDL2_ProcessEvent(&event);
            if (event.type == SDL_QUIT) done = true;
            if (event.type == SDL_WINDOWEVENT &&
                event.window.event == SDL_WINDOWEVENT_CLOSE &&
                event.window.windowID == SDL_GetWindowID(window))
                done = true;
        }

        // Resize swap chain?
        if (g_SwapChainRebuild) {
            int width, height;
            SDL_GetWindowSize(window, &width, &height);
            if (width > 0 && height > 0) {
                ImGui_ImplVulkan_SetMinImageCount(g_MinImageCount);
                ImGui_ImplVulkanH_CreateOrResizeWindow(
                    g_Instance, g_PhysicalDevice, g_Device, &g_MainWindowData,
                    g_QueueFamily, g_Allocator, width, height, g_MinImageCount);
                g_MainWindowData.FrameIndex = 0;
                g_SwapChainRebuild = false;
            }
        }

        // Start the Dear ImGui frame
        ImGui_ImplVulkan_NewFrame();
        ImGui_ImplSDL2_NewFrame();
        ImGui::NewFrame();

        show_cut_mesh();

        // Rendering
        ImGui::Render();
        ImDrawData* draw_data = ImGui::GetDrawData();
        const bool is_minimized = (draw_data->DisplaySize.x <= 0.0f ||
                                   draw_data->DisplaySize.y <= 0.0f);
        if (!is_minimized) {
            wd->ClearValue.color.float32[0] = clear_color.x * clear_color.w;
            wd->ClearValue.color.float32[1] = clear_color.y * clear_color.w;
            wd->ClearValue.color.float32[2] = clear_color.z * clear_color.w;
            wd->ClearValue.color.float32[3] = clear_color.w;
            FrameRender(wd, draw_data);
            FramePresent(wd);
        }
    }

    // Cleanup
    err = vkDeviceWaitIdle(g_Device);
    check_vk_result(err);
    ImGui_ImplVulkan_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();

    CleanupVulkanWindow();
    CleanupVulkan();

    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
