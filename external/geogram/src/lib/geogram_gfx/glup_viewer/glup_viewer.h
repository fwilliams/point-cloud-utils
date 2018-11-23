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
 */

#ifndef GLUP_VIEWER
#define GLUP_VIEWER

/**
 * \file geogram_gfx/glup_viewer/glup_viewer.h
 * \brief Low-level C API for the GLUP viewer.
 */

/* #define WITH_PNG  */

#ifdef __EMSCRIPTEN__
#define GLFW_INCLUDE_ES2 
#include <GLFW/glfw3.h>
#define GL_GLEXT_PROTOTYPES
#include <GLES2/gl2ext.h>
typedef double GLdouble;
#else
#include <geogram_gfx/third_party/glad/glad.h>
#endif

#include <geogram_gfx/GLUP/GLUP.h>

#ifdef __cplusplus
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#endif

#if defined(_MSC_VER) && defined(GEO_DYNAMIC_LIBS)
#ifdef geogram_gfx_EXPORTS
#define GLUP_VIEWER_API __declspec(dllexport) 
#else
#define GLUP_VIEWER_API __declspec(dllimport) 
#endif
#else
#if defined(__GNUC__) && defined(GEO_DYNAMIC_LIBS)
#define GLUP_VIEWER_API __attribute__ ((visibility("default")))
#else
#define GLUP_VIEWER_API
#endif

#endif


#ifdef __cplusplus
#include <string>
#include <cassert>
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum GlupViewerEvent {
    GLUP_VIEWER_DOWN,
    GLUP_VIEWER_MOVE,
    GLUP_VIEWER_UP
};
    
/**
 * \brief Callback function to be called to display each frame.
 * \see glup_viewer_set_display_func()
 */
typedef void (* GlupViewerDisplayFunc)(void);

/**
 * \brief Callback function to be called when a key is pressed.
 * \details For high-level input, typically key shortcuts for command.
 *  Videogames should use glup_viewer_set_keyboard_func_ext. The function
 *  is supposed to return GL_TRUE if the key was taken into account, GL_FALSE
 *  otherwise.
 * \see glup_viewer_set_keyboard_func()
 */
typedef GLboolean (* GlupViewerKeyboardFunc)(char key);

/**
 * \brief Callback function to be called when a key is pressed.
 * \details For low-level input, typically videogames. The function
 *  is supposed to return GL_TRUE if the key was taken into account, GL_FALSE
 *  otherwise. The key argument is either a single-character string with the
 *  key, or the symbolic name of the key ("left", "right", "alt", "control",
 *  "tab", "esc"). The ev argument is one of GLUP_VIEWER_DOWN, GLUP_VIEWER_UP.
 * \see glup_viewer_set_keyboard_func()
 */
typedef GLboolean (* GlupViewerKeyboardFuncExt)(
    const char* key, enum GlupViewerEvent ev
);
    
/**
 * \brief Callback function associated to an individual key.
 * \see glup_viewer_set_key_func()
 */
typedef void (* GlupViewerKeyFunc)(void);

/**
 * \brief Callback function called when the first frame is displayed.
 * \details Can be used to setup OpenGL objects, create textures etc...
 * \see glup_viewer_set_init_func()
 */
typedef void (* GlupViewerInitFunc)(void);

/**
 * \brief Callback function called when a file is dropped in the window.
 * \see glup_viewer_set_drag_drop_func()
 */
typedef void (* GlupViewerDragDropFunc)(char*);


typedef GLboolean (* GlupViewerMouseFunc)(
    float x, float y, int button, enum GlupViewerEvent event
);

#define GLUP_VIEWER_IDLE_REDRAW 0
#define GLUP_VIEWER_DRAW_SCENE 1
#define GLUP_VIEWER_SHOW_HELP 2
#define GLUP_VIEWER_BACKGROUND 3

#define GLUP_VIEWER_ROTATE_LIGHT 5
#define GLUP_VIEWER_3D 6
#define GLUP_VIEWER_TWEAKBARS 7

#define GLUP_VIEWER_STEREOSCOPIC_DISPLAY 11
#define GLUP_VIEWER_CLIP 12
#define GLUP_VIEWER_SHOW_CLIP 13
#define GLUP_VIEWER_EDIT_CLIP 14
#define GLUP_VIEWER_FIXED_CLIP 15
#define GLUP_VIEWER_FULL_SCREEN 16

#define GLUP_VIEWER_NO_EFFECT 0
#define GLUP_VIEWER_AMBIENT_OCCLUSION 1
#define GLUP_VIEWER_UNSHARP_MASKING 2   
#define GLUP_VIEWER_TEST_EFFECT 3
    
extern GLUP_VIEWER_API void glup_viewer_enable(int cap);
extern GLUP_VIEWER_API void glup_viewer_disable(int cap);
extern GLUP_VIEWER_API GLboolean glup_viewer_is_enabled(int cap);
extern GLUP_VIEWER_API GLboolean* glup_viewer_is_enabled_ptr(int cap);
extern GLUP_VIEWER_API void glup_viewer_toggle(int cap);

#define GLUP_VIEWER_STEREOSCOPIC_EYE_DISTANCE 5
#define GLUP_VIEWER_ZOOM 6
    
extern GLUP_VIEWER_API void glup_viewer_set_float(int param, GLfloat value);
extern GLUP_VIEWER_API GLfloat glup_viewer_get_float(int param);
extern GLUP_VIEWER_API GLfloat* glup_viewer_float_ptr(int param);

extern GLUP_VIEWER_API void glup_viewer_main_loop(int argc, char** argv);
extern GLUP_VIEWER_API void glup_viewer_exit_main_loop(void);
extern GLUP_VIEWER_API void glup_viewer_set_window_title(const char* title);
    
extern GLUP_VIEWER_API void glup_viewer_set_display_func(
    GlupViewerDisplayFunc f
);

extern GLUP_VIEWER_API void glup_viewer_set_overlay_func(
    GlupViewerDisplayFunc f
);

extern GLUP_VIEWER_API void glup_viewer_set_keyboard_func(
    GlupViewerKeyboardFunc f
);

extern GLUP_VIEWER_API void glup_viewer_set_keyboard_func_ext(
    GlupViewerKeyboardFuncExt f
);

extern GLUP_VIEWER_API void glup_viewer_set_mouse_func(GlupViewerMouseFunc f);
extern GLUP_VIEWER_API void glup_viewer_set_init_func(GlupViewerInitFunc f);
extern GLUP_VIEWER_API void glup_viewer_set_drag_drop_func(
    GlupViewerDragDropFunc f
);

extern GLUP_VIEWER_API GlupViewerDisplayFunc glup_viewer_get_display_func(void);
extern GLUP_VIEWER_API GlupViewerDisplayFunc glup_viewer_get_overlay_func(void);
extern GLUP_VIEWER_API GlupViewerKeyboardFunc glup_viewer_get_keyboard_func(
    void
);
    
extern GLUP_VIEWER_API GlupViewerMouseFunc glup_viewer_get_mouse_func(void);
extern GLUP_VIEWER_API GlupViewerInitFunc glup_viewer_get_init_func(void);
extern GLUP_VIEWER_API GlupViewerDragDropFunc glup_viewer_get_drag_drop_func(
    void
);
    
extern GLUP_VIEWER_API void glup_viewer_add_toggle(
    char key, GLboolean* pointer, const char* description
);
extern GLUP_VIEWER_API void glup_viewer_add_key_func(
    char key, GlupViewerKeyFunc f, const char* description
);
extern GLUP_VIEWER_API void glup_viewer_unbind_key(char key);
extern GLUP_VIEWER_API void glup_viewer_set_region_of_interest(
    float xmin, float ymin, float zmin, float xmax, float ymax, float zmax
);
extern GLUP_VIEWER_API void glup_viewer_get_region_of_interest(
    float* xmin, float* ymin, float* zmin, float* xmax, float* ymax, float* zmax
);
extern GLUP_VIEWER_API void glup_viewer_set_screen_size(int w, int h);
extern GLUP_VIEWER_API void glup_viewer_get_screen_size(int* w, int* h);

extern GLUP_VIEWER_API void glup_viewer_clear_text(void);
extern GLUP_VIEWER_API void glup_viewer_printf(const char* format, ...);
extern GLUP_VIEWER_API void glup_viewer_set_background(int texture);
extern GLUP_VIEWER_API int glup_viewer_get_background(void);
extern GLUP_VIEWER_API void glup_viewer_set_background_color(
    GLfloat r, GLfloat g, GLfloat b
);
extern GLUP_VIEWER_API void glup_viewer_set_background_color2(
    GLfloat r, GLfloat g, GLfloat b
);
extern GLUP_VIEWER_API GLfloat* glup_viewer_get_background_color(void);
extern GLUP_VIEWER_API GLfloat* glup_viewer_get_background_color2(void);

extern GLUP_VIEWER_API float* glup_viewer_get_scene_quaternion(void);
extern GLUP_VIEWER_API float* glup_viewer_get_scene_translation(void);
extern GLUP_VIEWER_API float* glup_viewer_get_light_matrix(void);
extern GLUP_VIEWER_API float* glup_viewer_get_light_quaternion(void);
extern GLUP_VIEWER_API float* glup_viewer_get_clip_quaternion(void);

extern GLUP_VIEWER_API void glup_viewer_set_scene_rotation(
    float xyz[3], float angle
);
    
extern GLUP_VIEWER_API void glTexImage2DXPM(char const* const* xpm_data);
extern GLUP_VIEWER_API void glTexImage2Dfile(const char* filename);

extern GLUP_VIEWER_API GLboolean glup_viewer_load_image(
    const char* filename, GLuint* width, GLuint* height, GLuint* bpp,
    GLvoid** pixels
);

extern GLUP_VIEWER_API int glup_viewer_fps(void);

extern GLUP_VIEWER_API void glup_viewer_redraw(void);

extern GLUP_VIEWER_API void glup_viewer_save_transform_for_picking(void);
extern GLUP_VIEWER_API void glup_viewer_get_picked_ray(
    GLdouble* p, GLdouble* v
);
extern GLUP_VIEWER_API void glup_viewer_get_picked_point(
    GLdouble* p, GLboolean* hit_background
);
extern GLUP_VIEWER_API void glup_viewer_project(
    double x_in, double y_in, double z_in,
    double* x_out, double* y_out, double* z_out
);
extern GLUP_VIEWER_API void glup_viewer_unproject(
    double x_in, double y_in, double z_in,
    double* x_out, double* y_out, double* z_out
);
    
extern GLUP_VIEWER_API void glup_viewer_random_color(void);
    
extern GLUP_VIEWER_API void glup_viewer_randomize_colors(int);
    
extern GLUP_VIEWER_API void glup_viewer_random_color_from_index(int index);    

extern GLUP_VIEWER_API void glup_viewer_post_redisplay(void);

extern GLUP_VIEWER_API void glup_viewer_home(void);

extern GLUP_VIEWER_API GLboolean glup_viewer_is_high_dpi(void);

/**
 * \brief Saves a copy of the current window to a file.
 * \details Uses the PPM file format.
 * \param[in] filename name of the file
 */
extern GLUP_VIEWER_API void glup_viewer_snapshot(const char* filename);

/**
 * \brief Sets the effect.
 * \param[in] effect one of GLUP_VIEWER_NO_EFFECT, GLUP_VIEWER_AMBIENT_OCCLUSION, 
 *  GLUP_VIEWER_UNSHARP_MASKING
 */
extern GLUP_VIEWER_API GLboolean glup_viewer_set_effect(GLenum effect);
    
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
inline void glup_viewer_printf(const std::string& s) {
    glup_viewer_printf(s.c_str());
}

inline void glup_viewer_add_toggle(
    char key, bool* toggle, const char* description
) {
    assert(sizeof(bool) == sizeof(GLboolean));
    glup_viewer_add_toggle(key, (GLboolean*)(toggle), description);
}


#endif

#endif

