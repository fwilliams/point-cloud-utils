/*
 *    _____   _       _   _     ____   
 *   /  ___| | |     | | | |   /  _ \  
 *   | |     | |     | | | |   | |_\ \ 
 *   | |  _  | |     | | | |   |  __ / 
 *   | |_| | | |___  | |_| |   | |     
 *   \_____/ |_____| \_____/   |_|     
 *
 *    _     _   _   _____   _          __  _____   _____
 *   | |   / / | | | ____| | |        / / | ____| |  _  \
 *   | |  / /  | | | |__   | |  __   / /  | |__   | |_| |
 *   | | / /   | | |  __|  | | /  | / /   |  __|  |  _  /
 *   | |/ /    | | | |___  | |/   |/ /    | |___  | | \ \
 *   |___/     |_| |_____| |___/|___/     |_____| |_|  \_\
 *
 *  Version 1.0
 *  Bruno Levy, April 2016
 *  INRIA, Project ALICE
 *
 *  Used internally by GLUP Viewer for interfacing with ImGUI
 *
 */

#include <geogram_gfx/basic/GL.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui_private.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/full_screen_effects/ambient_occlusion.h>
#include <geogram_gfx/full_screen_effects/unsharp_masking.h>
#include <geogram_gfx/third_party/ImGui/imgui.h>
#include <geogram_gfx/third_party/ImGui/imgui_impl_glfw.h>
#include <geogram_gfx/third_party/ImGui/imgui_impl_opengl2.h>
#include <geogram_gfx/third_party/ImGui/imgui_impl_opengl3.h>
#include <geogram_gfx/third_party/quicktext/glQuickText.h>
#include <geogram_gfx/third_party/imgui_fonts/roboto_medium.h>
#include <geogram_gfx/third_party/imgui_fonts/cousine_regular.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/command_line.h>

/* 
 *Lots of documentation tags in GLFW that are
 * not understood by CLANG.
 */
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wdocumentation"
#endif

#ifdef GEO_OS_EMSCRIPTEN
#include <GLFW/glfw3.h>
#include <emscripten.h>
#include <geogram/basic/file_system.h>
#else
#include <third_party/glfw/include/GLFW/glfw3.h>
#endif

#include <geogram_gfx/GLUP/GLUP.h>

#include <string.h>
#include <iostream>

#ifdef GEO_GL_LEGACY
static bool vanillaGL = false;
#endif

/***************************************************************************/

static GEO::FullScreenEffectImpl_var effect_;


// Commands may need to update the GUI (when using
// the console or the progressbar). The Command class
// does that outside of the ImGui handler, but if client
// code does that directly, then we may have two nested
// invocations of the ImGui handler, which is not correct.
// This variable avoids to have two nested invocations of
// the ImGui handler.
static bool glup_viewer_gui_locked = false;

extern "C" {
    void glup_viewer_one_frame(void);
}

void glup_viewer_gui_update() {
    // It's a pity, but under Emscripten, only the browser can have the
    // control of the rendering loop, therefore it is not possible (or I
    // don't know how) to update the graphics during computations.
#ifndef GEO_OS_EMSCRIPTEN
    glup_viewer_post_redisplay();    
    if(!glup_viewer_gui_locked) {
        glup_viewer_one_frame();
    }
#endif    
}

/***************************************************************************/

void glup_viewer_gui_init(GLFWwindow* w) {

    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); 
    ImGui_ImplGlfw_InitForOpenGL(w,false);
#ifdef GEO_GL_LEGACY        
    vanillaGL = (strcmp(glupCurrentProfileName(), "VanillaGL") == 0);
    if(vanillaGL) {
        GEO::Logger::out("ImGUI") << "Viewer GUI init (Vanilla)"
                                  << std::endl;        
        ImGui_ImplOpenGL2_Init();        
    } else
#endif
    {
        GEO::Logger::out("ImGUI") << "Viewer GUI init (GL3)"
                                  << std::endl;
#if defined(GEO_OS_APPLE)
        ImGui_ImplOpenGL3_Init("#version 330");	
#elif defined(GEO_OS_ANDROID)
	// TODO: have also version for OpenGL ES 2.0.
        ImGui_ImplOpenGL3_Init("#version 300 es");
#else	
        ImGui_ImplOpenGL3_Init("#version 100");
#endif	
    }

    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 10.0f;
    style.FrameRounding = 10.0f;
    style.GrabRounding = 10.0f;
    io.IniFilename = nullptr;

    io.FontDefault = io.Fonts->AddFontFromMemoryCompressedTTF(
	roboto_medium_compressed_data, roboto_medium_compressed_size, 16.0f
    );

    io.Fonts->AddFontFromMemoryCompressedTTF(
	roboto_medium_compressed_data, roboto_medium_compressed_size, 32.0f
    );

    io.FontDefault = io.Fonts->AddFontFromMemoryCompressedTTF(
	cousine_regular_compressed_data, cousine_regular_compressed_size, 16.0f
    );

    io.Fonts->AddFontFromMemoryCompressedTTF(
	cousine_regular_compressed_data, cousine_regular_compressed_size, 32.0f
    );
    
}

void glup_viewer_gui_cleanup() {
#ifdef GEO_GL_LEGACY        
    if(vanillaGL) {    
        ImGui_ImplOpenGL2_Shutdown();
    } else
#endif
    {
        ImGui_ImplOpenGL3_Shutdown();        
    }
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}


void glup_viewer_gui_begin_frame() {
    glup_viewer_gui_locked = true;
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplOpenGL2_NewFrame();        
    } else
#endif        
    {
        ImGui_ImplOpenGL3_NewFrame();
    }
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void glup_viewer_gui_end_frame() {
    GlupViewerDisplayFunc overlay_func = glup_viewer_get_overlay_func();
    if(overlay_func != nullptr) {
        overlay_func();
        ImGui::Render();
#ifdef GEO_GL_LEGACY            	
	if(vanillaGL) {
	    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
	} else
#endif	    
	{
	    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
	}
    }
    glup_viewer_gui_locked = false;
    // We flush the queued command here, once ImGui::Render() was
    // called, so that if it triggers a frame rendering (e.g. through
    // the Logger or the ProgressLogger), then ImGui calls are correctly
    // nested.
    GEO::Command::flush_queue();
    if(glup_viewer_gui_takes_input()) {
        glup_viewer_post_redisplay();
    }
}

int glup_viewer_gui_takes_input() {
    if(!glup_viewer_is_enabled(GLUP_VIEWER_TWEAKBARS)) {
        return 0;
    }
    return (
        ImGui::GetIO().WantCaptureMouse ||
        ImGui::GetIO().WantCaptureKeyboard
    ) ? 1 : 0;
}

void glup_viewer_gui_mouse_button_callback(
    GLFWwindow* window, int button, int action, int mods
) {
    glup_viewer_post_redisplay();
    ImGui_ImplGlfw_MouseButtonCallback(window, button, action, mods);
}

void glup_viewer_gui_scroll_callback(
    GLFWwindow* window, double xoffset, double yoffset
) {
    ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);
}

void glup_viewer_gui_key_callback(
    GLFWwindow* window, int key, int scancode, int action, int mods
) {
    ImGui_ImplGlfw_KeyCallback(window, key, scancode, action, mods);
}

void glup_viewer_gui_char_callback(GLFWwindow* window, unsigned int c) {
    glup_viewer_post_redisplay();    
    ImGui_ImplGlfw_CharCallback(window, c);
}

void glup_viewer_gui_resize(int width, int height) {
    glup_viewer_post_redisplay();    
    ImGui::GetIO().DisplaySize.x = float(width);
    ImGui::GetIO().DisplaySize.y = float(height);
}

GLboolean glup_viewer_get_arg_bool(const char* param) {
    return GEO::CmdLine::get_arg_bool(param) ? GL_TRUE : GL_FALSE;
}

GLboolean glup_viewer_test_arg_string(const char* param, const char* arg) {
    return (GEO::CmdLine::get_arg(param) == arg) ? GL_TRUE : GL_FALSE;
}

void glup_viewer_set_screen_size_from_args() {
    std::string geometry = GEO::CmdLine::get_arg("gfx:geometry");
    int w,h;
    sscanf(geometry.c_str(),"%dx%d",&w,&h);
    glup_viewer_set_screen_size(w,h);
}

GLboolean glup_viewer_set_effect(GLenum effect) {
    switch(effect) {
	case GLUP_VIEWER_NO_EFFECT:
	    effect_.reset();
	    break;
	case GLUP_VIEWER_AMBIENT_OCCLUSION:
	    effect_ = new GEO::AmbientOcclusionImpl();
	    break;
	case GLUP_VIEWER_UNSHARP_MASKING:
	    effect_ = new GEO::UnsharpMaskingImpl();
	    break;
    }
    GLboolean result = (effect_.is_null() || effect_->OK()) ? GL_TRUE : GL_FALSE;
    if(!effect_.is_null() && !effect_->OK()) {
	effect_.reset();
    }
    return result;
}

void glup_viewer_effect_begin_frame() {
    if(!effect_.is_null()) {
	int w,h;
	glup_viewer_get_screen_size(&w,&h);
	effect_->pre_render(GEO::index_t(w), GEO::index_t(h));
    }
}
    
void glup_viewer_effect_end_frame() {
    if(!effect_.is_null()) {
	effect_->post_render();
    }
}

#ifdef GEO_OS_EMSCRIPTEN

extern "C" {
    
    void drag_drop(GLFWwindow* w, int nb, const char** p);
    void file_system_changed_callback();

    /**
     * \brief This function is called by the HTML shell each
     *  time a file is loaded.
     * \details For now, it uses the fact that the latest loaded
     *  file appears first in the list (to be improved).
     */
    EMSCRIPTEN_KEEPALIVE void file_system_changed_callback() {
        std::vector<std::string> all_files;
        GEO::FileSystem::get_directory_entries("/",all_files);
        if(
            all_files.size() > 0 && 
            GEO::FileSystem::is_file(all_files[0])
        ) {
            const char* pname = all_files[0].c_str();
            drag_drop(nullptr, 1, &pname);
        }
    }
}
#endif
