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

/* 
 * Many documentation tags in GLFW that are
 * not understood by CLANG.
 */
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wdocumentation"
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#endif

#include "glup_viewer.h"

#include <stdlib.h>
#include <stdarg.h>


#ifdef GEO_USE_SYSTEM_GLFW3
#include <GLFW/glfw3.h>
#else
#include <third_party/glfw/include/GLFW/glfw3.h>
#endif

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

#include <geogram_gfx/glup_viewer/glup_viewer_gui_private.h>
#include <geogram_gfx/third_party/quicktext/glQuickText.h>

#include <stdio.h>

#include <math.h>
#include <string.h>
#include <ctype.h>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/types.h>
#include <sys/times.h>
#endif

#ifdef WITH_PNG
#include <png/png.h>
#endif

#ifndef __APPLE__
#include <malloc.h>
#endif

#if defined(__APPLE__) || defined(__linux__) || defined(__EMSCRIPTEN__)
#include <unistd.h>
static void glup_viewer_pause() {
    usleep(1000);
}
#elif defined(WIN32)
static void glup_viewer_pause() {
    Sleep(1);
}
#else
static void glup_viewer_pause() {
}
#endif


/*
 * CLANG complains with a const -> non-const cast
 * for a callback parameter.
 */
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#ifdef __EMSCRIPTEN__
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#endif

static int glup_viewer_needs_redraw = 1;
static double glup_viewer_start_time = 0.0;

void glup_viewer_post_redisplay() {
    /* 
     * We observed that triggering only one frame redraw can cause some
     * flickering under MacOS/X (maybe it uses triple buffering ?)
     */
    glup_viewer_needs_redraw += 2;
    /*
     * Do not queue more than 10 frame redraws
     */
    if(glup_viewer_needs_redraw > 10) {
        glup_viewer_needs_redraw = 10;        
    }
}

#define glup_viewer_assert(x)                                                  \
    if(!(x)) {                                                                 \
        fprintf(stderr, "%s: assertion fail %s -- %d\n",#x,__FILE__, __LINE__);\
        abort();                                                               \
    }

#define glup_viewer_assert_not_reached()                                       \
    {                                                                          \
        fprintf(stderr, "should not go there!! %s -- %d\n",__FILE__,__LINE__); \
        abort();                                                               \
    }

#define glup_viewer_argused(x) (void) x

static GlupViewerDisplayFunc display_func = NULL;
static GlupViewerDisplayFunc overlay_func = NULL;
static GlupViewerKeyboardFunc keyboard_func = NULL;
static GlupViewerKeyboardFuncExt keyboard_func_ext = NULL;
static GlupViewerMouseFunc mouse_func = NULL;
static GlupViewerInitFunc init_func = NULL;
static GlupViewerDragDropFunc drag_drop_func = NULL;

static GLboolean* toggle[256];
static char* key_description[256];
static GlupViewerKeyFunc key_func[256];

static GLboolean caps[17] = {
    GL_FALSE, /* IDLE_REDRAW  */
    GL_FALSE, /* DRAW_SCENE   */
    GL_FALSE, /* SHOW_HELP    */
    GL_FALSE, /* BACKGROUND   */
    GL_FALSE, /* unused       */
    GL_FALSE, /* ROTATE_LIGHT */
    GL_TRUE,  /* 3D           */
    GL_TRUE,  /* TWEAKBARS    */
    GL_FALSE, /* unused */
    GL_FALSE, /* unused */
    GL_FALSE, /* unused */
    GL_FALSE, /* STEREOSCOPIC_DISPLAY         */
    GL_FALSE, /* CLIP         */
    GL_TRUE,  /* SHOW_CLIP    */
    GL_FALSE, /* EDIT_CLIP    */
    GL_FALSE, /* FIXED_CLIP   */
    GL_FALSE  /* FULL_SCREEN  */   
};
static int nb_caps = 17;

static GLfloat params[7] = {
    0.0f,  /* unused */
    0.0f,  /* unused */
    0.0f,  /* unused */
    0.0f,  /* unused */
    0.0f,  /* unused */
    0.09f, /* STEREOSCOPIC_EYE_DISTANCE */
    1.0f   /* ZOOM                      */
};
static int nb_params = 7;

static const char* title = "g33>|< Viewer";

static GLuint background_tex = 0;

static GLFWwindow* glup_viewer_window = NULL;


static GLboolean transparent = GL_FALSE;

/* ========================== Trackball prototypes ========================= */
/* (from SGI, see copyright below) */

/*
 * trackball.h
 * A virtual trackball implementation
 * Written by Gavin Bell for Silicon Graphics, November 1988.
 */

/*
 * Pass the x and y coordinates of the last and current positions of
 * the mouse, scaled so they are from (-1.0 ... 1.0).
 *
 * The resulting rotation is returned as a quaternion rotation in the
 * first paramater.
 */
static void trackball(float q[4], float p1x, float p1y, float p2x, float p2y);

/*
 * Given two quaternions, add them together to get a third quaternion.
 * Adding quaternions to get a compound rotation is analagous to adding
 * translations to get a compound translation.  When incrementally
 * adding rotations, the first argument here should be the new
 * rotation, the second and third the total rotation (which will be
 * over-written with the resulting new total rotation).
 */
static void add_quats(float* q1, float* q2, float* dest);

/*
 * A useful function, builds a rotation matrix in Matrix based on
 * given quaternion.
 */
static void build_rotmatrix(float m[4][4], float q[4]);

/*
 * This function computes a quaternion based on an axis (defined by
 * the given vector) and an angle about which to rotate.  The angle is
 * expressed in radians.  The result is put into the third argument.
 */
static void axis_to_quat(float a[3], float phi, float q[4]);

/* ==================== Timing ========================================= */

static double now() {
#ifdef WIN32
    return (double)(GetTickCount()) / 1000.0;
#else
   struct tms now_tms;
   return (double)(times(&now_tms)) / 100.0;
/*    return (double)(clock()) / (double) (CLOCKS_PER_SEC); */
#endif
}

int glup_viewer_fps() {
    static int init = 1;
    static double ref = 0;
    static int frame = 0;
    static double result = 0;
    if(init) {
        ref = now();
        init = 0;
    }
    frame++;
    if((frame % 20) == 0) {
        double new_t = now();
        result = frame / (new_t - ref);
        ref = new_t;
        frame = 0;
    }
    return (int) result;
}

/* ==================== GlupViewer implementation ========================= */

static int glup_viewer_W = 800;
static int glup_viewer_H = 800;
static int last_x, last_y;
static float cur_rot[4] = {0.0, 0.0, 0.0, 0.0};
static float cur_xlat[3] = {0.0, 0.0, 0.0};
static float cur_rot_light[4] = {0.0, 0.0, 0.0, 0.0};
static float cur_rot_clip[4] = {0.0, 0.0, 0.0, 0.0};
static float cur_xlat_clip[3] = {0.0, 0.0, 0.0};

static enum {
    NONE, ROTATE, PAN, ZOOM
} mode = NONE;
static float xmin = -1, ymin = -1, zmin = -1, xmax = 1, ymax = 1, zmax = 1;
static float roi_radius = 1.4f;
static int window_w = 0, window_h = 0;

float* glup_viewer_get_scene_quaternion() {
    return cur_rot;
}

float* glup_viewer_get_scene_translation() {
    return cur_xlat;
}

float* glup_viewer_get_light_quaternion() {
    return cur_rot_light;
}

float* glup_viewer_get_clip_quaternion() {
    return cur_rot_clip;
}

void glup_viewer_enable(int cap) {
    glup_viewer_assert(cap < nb_caps);
    caps[cap] = GL_TRUE;
}

void glup_viewer_disable(int cap) {
    glup_viewer_assert(cap < nb_caps);
    caps[cap] = GL_FALSE;
}

GLboolean glup_viewer_is_enabled(int cap) {
    glup_viewer_assert(cap < nb_caps);
    return caps[cap];
}

GLboolean* glup_viewer_is_enabled_ptr(int cap) {
    glup_viewer_assert(cap < nb_caps);
    return &(caps[cap]);
}

void glup_viewer_toggle(int cap) {
    glup_viewer_assert(cap < nb_caps);
    caps[cap] = !caps[cap];
}

void glup_viewer_set_float(int param, GLfloat value) {
    glup_viewer_assert(param < nb_params);
    params[param] = value;
}

GLfloat glup_viewer_get_float(int param) {
    glup_viewer_assert(param < nb_params);
    return params[param];
}

GLfloat* glup_viewer_float_ptr(int param) {
    glup_viewer_assert(param < nb_params);
    return &(params[param]);
}

void glup_viewer_get_screen_size(int* w, int* h) {
    *w = glup_viewer_W;
    *h = glup_viewer_H;
}

void glup_viewer_set_screen_size(int w, int h) {
    glup_viewer_W = w;
    glup_viewer_H = h;
}

static void reshape(int w, int h) {

    int viewport_x = 0;
    int viewport_y = 0;
    int viewport_width = w;
    int viewport_height = h;
    
    glViewport(
        viewport_x, viewport_y, viewport_width, viewport_height
    );
        
    glup_viewer_W = w;
    glup_viewer_H = h;
    window_w = w;
    window_h = h;

    if(glup_viewer_is_enabled(GLUP_VIEWER_TWEAKBARS)) {
        glup_viewer_gui_resize(w,h);
    }
    
    glup_viewer_post_redisplay();
}

static GLboolean transform_saved = GL_FALSE;
static GLdouble modelview_save[16];
static GLdouble project_save[16];
static GLint viewport_save[4];

void glup_viewer_save_transform_for_picking() {
    transform_saved = GL_TRUE;
    glupGetMatrixdv(GLUP_MODELVIEW_MATRIX, modelview_save);
    glupGetMatrixdv(GLUP_PROJECTION_MATRIX, project_save);
    glGetIntegerv(GL_VIEWPORT, viewport_save);
}

static GLdouble ray_p_save[3];
static GLdouble ray_v_save[3];
static GLdouble ray_p3d_save[3];
static GLboolean ray_hit_background;

void glup_viewer_get_picked_ray(GLdouble* p, GLdouble* v) {
    unsigned int i;
    for(i = 0; i < 3; i++) {
        p[i] = ray_p_save[i];
        v[i] = ray_v_save[i];
    }
}

void glup_viewer_get_picked_point(GLdouble* p, GLboolean* hit_background) {
    unsigned int i;
    for(i = 0; i < 3; i++) {
        p[i] = ray_p3d_save[i];
    }
    *hit_background = ray_hit_background;
}

void glup_viewer_project(
    double x_in, double y_in, double z_in,
    double* x_out, double* y_out, double* z_out
) {
    glupProject(
        x_in, y_in, z_in,
        modelview_save, project_save, viewport_save,
        x_out, y_out, z_out
    );
}

void glup_viewer_unproject(
    double x_in, double y_in, double z_in,
    double* x_out, double* y_out, double* z_out
) {
    glupUnProject(
        x_in, y_in, z_in,
        modelview_save, project_save, viewport_save,
        x_out, y_out, z_out
    );
}

static void save_picked_ray(int x_in, int y_in) {
    
   if(transform_saved) {
       unsigned int i;
       GLdouble temp[3];
       GLdouble x = (GLdouble) x_in;
       GLdouble y = (GLdouble) (viewport_save[3] - y_in);
       GLfloat z;

        /*
         * Get the correct value of Z
         * (since even 2D mode uses perspective xform, if we do not do
         * that, we cannot get correct x,y).
         */

       glupProject(
           x, y, 0.0, modelview_save, project_save, viewport_save,
           &(temp[0]), &(temp[1]), &(temp[2])
       );
       
       glupUnProject(
           x, y, temp[2], modelview_save, project_save, viewport_save,
           &(ray_p_save[0]), &(ray_p_save[1]), &(ray_p_save[2])
       );
        
       glupUnProject(
           x, y, 1.0, modelview_save, project_save, viewport_save,
           &(ray_v_save[0]), &(ray_v_save[1]), &(ray_v_save[2])
       );
        
       for(i = 0; i < 3; i++) {
           ray_v_save[i] -= ray_p_save[i];
       }

#ifdef __EMSCRIPTEN__
       /* Cannot read depth buffer with OpenGL/ES2.0, zutalors !! */
       z = 0.5f;
#else       
       glReadPixels(
           (GLint) x, (GLint) y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z
       );
#endif
       
       glupUnProject(
           x, y, (double)(z),
           modelview_save, project_save, viewport_save,
           &(ray_p3d_save[0]), &(ray_p3d_save[1]), &(ray_p3d_save[2])
       );
       
       ray_hit_background = (z == 1.0f);
   }
}

static GLboolean call_mouse_func(
    int x, int y, int button, enum GlupViewerEvent event
) {
    double l = (double) (window_w > window_h ? window_w : window_h) / 2.0;
    double fx = x - (double) window_w / 2.0;
    double fy = y - (double) window_h / 2.0;
    int mx = (int) (3000.0 * fx / l);
    int my = -(int) (3000.0 * fy / l);
    GLboolean result = GL_FALSE;
    save_picked_ray(x, y);
    if(mouse_func == NULL) {
        return result;
    }
    result = mouse_func((float) mx, (float) my, button, event);
    glup_viewer_post_redisplay();
    return result;
}

static int mouse_x = 0;
static int mouse_y = 0;
static int mouse_pressed = 0;

static void mouse(GLFWwindow* w, int button, int action, int mods) {
    glup_viewer_post_redisplay();
    glup_viewer_gui_mouse_button_callback(w, button, action, mods);
    if(glup_viewer_gui_takes_input()) {
        return;
    }
    
    mouse_pressed = (action == GLFW_PRESS);
    if(call_mouse_func(
           mouse_x, mouse_y, button,
           (action == GLFW_PRESS) ? GLUP_VIEWER_DOWN : GLUP_VIEWER_UP)
    ) {
        return;
    }
    mode = NONE;
    if(action != GLFW_PRESS) {
        return;
    }
    switch(button) {
        case GLFW_MOUSE_BUTTON_1:
            mode = PAN;
            break;
        case GLFW_MOUSE_BUTTON_3:
            mode = ZOOM;
            break;
        case GLFW_MOUSE_BUTTON_2:
            mode = ROTATE;
            break;
    }
    last_x = mouse_x;
    last_y = mouse_y;
}

static void scroll(GLFWwindow* w, double xoffset, double yoffset) {
    glup_viewer_post_redisplay();
    glup_viewer_gui_scroll_callback(w, xoffset, yoffset);    
    if(glup_viewer_gui_takes_input()) {
        return;
    }
    
/* 
 * Under emscripten and apple, mouse wheel is inversed 
 * as compared to the other platforms. 
 */
#if defined(__EMSCRIPTEN__) || defined(__APPLE__)
    yoffset *= -1.0;
#endif    
    
    if(yoffset > 0.0) {
        *glup_viewer_float_ptr(GLUP_VIEWER_ZOOM) /= 1.1f;
        glup_viewer_post_redisplay();
    } else if(yoffset < 0.0) {
        *glup_viewer_float_ptr(GLUP_VIEWER_ZOOM) *= 1.1f;
        glup_viewer_post_redisplay();
    }
}

static void motion(int x, int y);

static void passive_mouse(GLFWwindow* w, double xf, double yf) {
    glup_viewer_argused(w);    
    glup_viewer_post_redisplay();
    if(glup_viewer_gui_takes_input()) {
        return;
    }
    
    mouse_x = (int)(xf);
    mouse_y = (int)(yf);

    if(mouse_pressed) {
        motion(mouse_x,mouse_y);
    }
}

static void motion(int x, int y) {

    int W = window_w;
    int H = window_h;

    if(call_mouse_func(x, y, -1, GLUP_VIEWER_MOVE)) {
        return;
    }

    switch(mode) {
        case ROTATE: {
            float delta_rot[4];
            trackball(delta_rot,
                (float) (2 * last_x - W) / (float) W,
                (float) (H - 2 * last_y) / (float) H,
                (float) (2 * x - W) / (float) W,
                (float) (H - 2 * y) / (float) H
            );
            if(glup_viewer_is_enabled(GLUP_VIEWER_ROTATE_LIGHT)) {
                add_quats(delta_rot, cur_rot_light, cur_rot_light);
            } else {
                if(
                    glup_viewer_is_enabled(GLUP_VIEWER_EDIT_CLIP) ||
                    !glup_viewer_is_enabled(GLUP_VIEWER_FIXED_CLIP)
                ) {
                    add_quats(delta_rot, cur_rot_clip, cur_rot_clip);
                }
                if(
                    !glup_viewer_is_enabled(GLUP_VIEWER_CLIP) ||
                    !glup_viewer_is_enabled(GLUP_VIEWER_EDIT_CLIP)
                ) {
                    add_quats(delta_rot, cur_rot, cur_rot);
                }
            }
        } break;
        case PAN: {
            if(!glup_viewer_is_enabled(GLUP_VIEWER_ROTATE_LIGHT)) {
                float delta_x = (float) (last_x - x) / (float) W;
                float delta_y = (float) (y - last_y) / (float) H;
                if(!glup_viewer_is_enabled(GLUP_VIEWER_EDIT_CLIP)) {
                    cur_xlat[0] -= 2.0f * delta_x / params[GLUP_VIEWER_ZOOM];
                    cur_xlat[1] -= 2.0f * delta_y / params[GLUP_VIEWER_ZOOM];
                }
                if(
                    glup_viewer_is_enabled(GLUP_VIEWER_EDIT_CLIP) ||
                    !glup_viewer_is_enabled(GLUP_VIEWER_FIXED_CLIP)
                ) {
                    cur_xlat_clip[0] -=
                        2.0f * delta_x / params[GLUP_VIEWER_ZOOM];
                    cur_xlat_clip[1] -=
                        2.0f * delta_y / params[GLUP_VIEWER_ZOOM];
                }
            }
        } break;
        case ZOOM: {
            if(!glup_viewer_is_enabled(GLUP_VIEWER_ROTATE_LIGHT)) {
                params[GLUP_VIEWER_ZOOM] *=
                    (1.0f + (float) (y - last_y) / (float) H);
            }
        } break;
        case NONE:
            break;
    }

    last_x = x;
    last_y = y;
    glup_viewer_post_redisplay();
}

void drag_drop(GLFWwindow* w, int nb, const char** p);

void drag_drop(GLFWwindow* w, int nb, const char** p) {
    int i;
    glup_viewer_post_redisplay();
    glup_viewer_argused(w);    
    if(drag_drop_func != NULL) {
        for(i=0; i<nb; ++i) {
            drag_drop_func((char*)(p[i]));
        }
        glup_viewer_post_redisplay();
    }
}

static float cur_text_x = 0;
static float cur_text_y = 0;

void glup_viewer_clear_text() {
    cur_text_x = -2800;
    cur_text_y = 2800;
    if(glup_viewer_W > glup_viewer_H) {
        cur_text_y *= ((float) glup_viewer_H / (float) glup_viewer_W);
    } else {
        cur_text_x *= ((float) glup_viewer_W / (float) glup_viewer_H);
    }
}

void glup_viewer_printf(const char* format, ...) {
    va_list args;
    char buffer[1024];
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    glupSetColor3f(GLUP_FRONT_AND_BACK_COLOR, 1.0, 1.0, 1.0);
    glQuickTextPrintString(
        (double)(cur_text_x), (double)(cur_text_y), 0.0, 10.0, buffer
    );
    cur_text_y -= 200;
}
 

static void draw_foreground() {
    float aspect;
    GLboolean vertex_colors_save;
    GLboolean texturing_save;
    GLboolean mesh_save;
    GLUPfloat shrink_save;
    GLUPfloat side;
    
    glup_viewer_clear_text();
    glupMatrixMode(GLUP_PROJECTION_MATRIX);
    glupLoadIdentity();
    aspect = (float) (glup_viewer_W) / (float) (glup_viewer_H);
    if(aspect > 1.0f) {
        side = 3000 * aspect;
    } else {
        side = 3000 / aspect;
    }

    glupOrtho2D(
        (double)(-side), (double)(side), (double)(-side), (double)(side)
    );
    
    glupMatrixMode(GLUP_MODELVIEW_MATRIX);
    glupLoadIdentity();
    glDisable(GL_DEPTH_TEST);
    glupDisable(GLUP_LIGHTING);

    if(
        glup_viewer_is_enabled(GLUP_VIEWER_SHOW_HELP) &&
        !glup_viewer_is_enabled(GLUP_VIEWER_TWEAKBARS)         
    ) {
        int i;
        
        glupDisable(GLUP_CLIPPING);
        
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glupSetColor4f(GLUP_FRONT_AND_BACK_COLOR, 0, 0, 0.5, 0.5);

        vertex_colors_save = glupIsEnabled(GLUP_VERTEX_COLORS);
        texturing_save = glupIsEnabled(GLUP_TEXTURING);
        mesh_save = glupIsEnabled(GLUP_DRAW_MESH);
        shrink_save = glupGetCellsShrink();
        
        glupDisable(GLUP_VERTEX_COLORS);
        glupDisable(GLUP_TEXTURING);
        glupDisable(GLUP_DRAW_MESH);
        glupSetCellsShrink(0.0f);
        
        glupBegin(GLUP_TRIANGLES);
        glupVertex2f(-side, -side);
        glupVertex2f(-side, side);
        glupVertex2f(side, side);
        glupVertex2f(-side, -side);
        glupVertex2f(side, side);        
        glupVertex2f(side, -side);
        glupEnd();
        
        glDisable(GL_BLEND);
        glupColor4f(5, 5, 5, 1);

        glupSetMeshWidth(1);

        glupColor3f(5, 5, 5);

        if(glup_viewer_is_enabled(GLUP_VIEWER_IDLE_REDRAW)) {
            glup_viewer_printf(
                " --- %s help [%4d FPS ] --- ", title, glup_viewer_fps()
            );
        } else {
            glup_viewer_printf(" --- %s help --- ", title);
        }
        glup_viewer_printf("");

        for(i = 0; i < 256; i++) {
            if(!isprint(i)) {
                continue;
            }
            if(key_description[i] != NULL) {
                if(toggle[i] != NULL) {
                    glup_viewer_printf(
                        "%c:toggle %s (%s)\n",
                        (char) i, key_description[i],
                        *toggle[i] ? "on" : "off"
                    );
                }
                if(key_func[i] != NULL) {
                    glup_viewer_printf("%c:%s\n", (char) i, key_description[i]);
                }
            }
        }
        glDisable(GL_BLEND);
        if(vertex_colors_save) {
            glupEnable(GLUP_VERTEX_COLORS);
        }
        if(texturing_save) {
            glupEnable(GLUP_TEXTURING);
        }
        if(mesh_save) {
            glupEnable(GLUP_DRAW_MESH);
        }
        glupSetCellsShrink(shrink_save);
    }
}

static void display(void);

static GLboolean in_main_loop_ = GL_FALSE;
void glup_viewer_redraw() {
    if(in_main_loop_) {
        display();
    }
}

#ifndef M_PI
#define M_PI 3.14159
#endif

static float to_radians(float x) {
    return x * 2.0f * (float) M_PI / 360.0f;
}

GLfloat* glup_viewer_get_light_matrix() {
    static GLfloat m[4][4];
    build_rotmatrix(m, cur_rot_light);
    return &m[0][0];
}

static float bkg1[3] = {1.0, 1.0, 1.0};
static float bkg2[3] = {0.0, 0.0, 0.5};

void glup_viewer_set_background_color(GLfloat r, GLfloat g, GLfloat b) {
    bkg1[0] = r;
    bkg1[1] = g;
    bkg1[2] = b;
}

void glup_viewer_set_background_color2(GLfloat r, GLfloat g, GLfloat b) {
    bkg2[0] = r;
    bkg2[1] = g;
    bkg2[2] = b;
}

GLfloat* glup_viewer_get_background_color() {
    return bkg1;
}

GLfloat* glup_viewer_get_background_color2() {
    return bkg2;
}

static void draw_background() {
    float z = 1.0f;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_POLYGON_OFFSET_FILL);    
    glClear((GLbitfield) (GL_DEPTH_BUFFER_BIT));
    glupMatrixMode(GLUP_PROJECTION_MATRIX);
    glupLoadIdentity();
    glupMatrixMode(GLUP_MODELVIEW_MATRIX);
    glupLoadIdentity();
    
    glupDisable(GLUP_LIGHTING);
    glupDisable(GLUP_TEXTURING);
    glupDisable(GLUP_CLIPPING);
    glupDisable(GLUP_DRAW_MESH);
    glupDisable(GLUP_PICKING);
    glupSetCellsShrink(0.0f);

    if(background_tex == 0) {

	if(transparent) {
	    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	    glClear(GL_COLOR_BUFFER_BIT);
	} else {
	    glupEnable(GLUP_VERTEX_COLORS);        
	    glupBegin(GLUP_QUADS);
	    glupColor3fv(bkg1);
	    glupVertex3f(-1, -1, z);
	    glupVertex3f(1, -1, z);
	    glupColor3fv(bkg2);
	    glupVertex3f(1, 1, z);
	    glupVertex3f(-1, 1, z);
	    glupEnd();
	}
        
    } else {
        
        float w = 1.0f;
        float h = 1.0f;
        
        glupEnable(GLUP_TEXTURING);
        glupTextureType(GLUP_TEXTURE_2D);
        glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
        glBindTexture(GL_TEXTURE_2D, background_tex);
        glupColor4f(1, 1, 1, 1);
        glupBegin(GLUP_QUADS);
        glupTexCoord2f(0, 0);
        glupVertex3f(-w, -h, z);
        glupTexCoord2f(1, 0);
        glupVertex3f(w, -h, z);
        glupTexCoord2f(1, 1);
        glupVertex3f(w, h, z);
        glupTexCoord2f(0, 1);
        glupVertex3f(-w, h, z);
        glupEnd();
        glupDisable(GLUP_TEXTURING);
    }
    glEnable(GL_DEPTH_TEST);
    glupEnable(GLUP_LIGHTING);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glupDisable(GLUP_VERTEX_COLORS);
}

static void transform_vector(
    GLfloat* vtrans, const GLfloat* M, const GLfloat* v
) {
    unsigned int i,j;
    for(i=0; i<4; ++i) {
        vtrans[i] = 0.0;
        for(j=0; j<4; ++j) {
            vtrans[i] += M[j*4+i]*v[j];
        }
    }
}

static void actually_render_display(double offset) {
    GLfloat m[4][4];

    /* Infinite light location (normalize (1,1,1) */
    static GLfloat light_position[] = {0.57735f, 0.57735f, 0.57735f, 0.0f};

    /* Transformed light vector */
    GLfloat light_pos_transformed[4];
    
    /* field of view of the larger dimension in degrees */
    float camera_aperture = 9.0f;

    float zoom = params[GLUP_VIEWER_ZOOM];
    float zNear = 1.0f;                /* near clipping plane     */
    float zFar = 10.0f;                /* far clipping plane      */
    float zScreen = 5.0f;              /* screen projection plane */

    /* half of the distance between the eyes, if in stereo mode */
    float eye_offset = (float) offset;

    /* aspect ratio of the window */
    float aspect = (float) glup_viewer_W / (float) glup_viewer_H;

    /* half the width of the screen */
    float vue_max_size = (float) (
        zScreen * (float) tan((double)(to_radians(camera_aperture) / 2.0f))
    );

    /* from the central point of view
       shift of the vue from the current point of view */
    float vue_shift = eye_offset * zNear / zScreen;

    int i;
    
    double clip_eqn[4];

    glup_viewer_effect_begin_frame();
    
    if(glup_viewer_is_enabled(GLUP_VIEWER_BACKGROUND)) {
        draw_background();
    } else {
        glClearColor(bkg1[0], bkg1[1], bkg1[2], 1.0f);
        glClear((GLbitfield) (GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT));
    }

    glEnable(GL_DEPTH_TEST);

    glupMatrixMode(GLUP_PROJECTION_MATRIX);
    glupLoadIdentity();

    if(glup_viewer_is_enabled(GLUP_VIEWER_3D)) {
        float right;
        float top;
        
        if(aspect < 1) {
            top = vue_max_size;
            right = top * aspect;
        } else {
            right = vue_max_size;
            top = right / aspect;
        }
        right /= zoom;
        top /= zoom;
        glupFrustum(
            (double)(-right - vue_shift), (double)(right - vue_shift),
            (double)(-top), (double)(top),
            (double)(zNear), (double)(zFar)
        );
        glupTranslatef(-eye_offset, 0.0, 0.0);
    } else {
        float x = 1.0f / zoom;
        float y = 1.0f / zoom;
        if(aspect > 1.0f) {
            x *= aspect;
        } else {
            y /= aspect;
        }
        glupOrtho(
            (double)(-x), (double)(x),
            (double)(-y), (double)(y),
            (double)(zNear), (double)(zFar)
        );
    }


    build_rotmatrix(m, cur_rot_light);
    transform_vector(light_pos_transformed, &m[0][0], light_position);
    glupLightVector3fv(light_pos_transformed);

    
    glupMatrixMode(GLUP_MODELVIEW_MATRIX);    
    glupLoadIdentity();

    glupDisable(GLUP_CLIPPING);
    if(glup_viewer_is_enabled(GLUP_VIEWER_CLIP)) {
        float* background_color;
        
        glupPushMatrix();

        /* translate the world of the distance between eye and center */
        glupTranslatef(0.0, 0.0, -zScreen);
        glupTranslatef(cur_xlat_clip[0], cur_xlat_clip[1], cur_xlat_clip[2]);
        build_rotmatrix(m, cur_rot_clip);
        glupMultMatrixf(&m[0][0]);

        background_color = glup_viewer_get_background_color();
        glupSetColor3f(
            GLUP_FRONT_AND_BACK_COLOR,
            1.0f - background_color[0],
            1.0f - background_color[1],
            1.0f - background_color[2]
        );


        if(glup_viewer_is_enabled(GLUP_VIEWER_SHOW_CLIP)) {

            float sq_w = 1.25f / params[GLUP_VIEWER_ZOOM];

            GLboolean vertex_colors_save = glupIsEnabled(GLUP_VERTEX_COLORS);
            GLboolean texturing_save = glupIsEnabled(GLUP_TEXTURING);
            glupDisable(GLUP_VERTEX_COLORS);
            glupDisable(GLUP_TEXTURING);

            /* Draw the cross */
            glupBegin(GLUP_LINES);
            glupVertex3f(-sq_w, 0.0f, 0.0f);
            glupVertex3f(sq_w, 0.0f, 0.0f);
            glupVertex3f(0.0f, -sq_w, 0.0f);
            glupVertex3f(0.0f, sq_w, 0.0f);

            /* Draw the square around the cross */
            for(i=0; i<3; ++i) {
                glupVertex3f(sq_w, -sq_w, 0.0f);
                glupVertex3f(sq_w, sq_w, 0.0f);
                glupVertex3f(sq_w, sq_w, 0.0f);            
                glupVertex3f(-sq_w, sq_w, 0.0f);
                glupVertex3f(-sq_w, sq_w, 0.0f);            
                glupVertex3f(-sq_w, -sq_w, 0.0f);
                glupVertex3f(-sq_w, -sq_w, 0.0f);
                glupVertex3f(sq_w, -sq_w, 0.0f);
                sq_w = sq_w * 1.01f;
            }
            
            glupEnd();
            
            if(vertex_colors_save) {
                glupEnable(GLUP_VERTEX_COLORS);
            }

            if(texturing_save) {
                glupEnable(GLUP_TEXTURING);
            }
            
        }

        clip_eqn[0] = 0.0;
        clip_eqn[1] = 0.0;
        clip_eqn[2] = -1.0;
        clip_eqn[3] = 0.0;
        glupEnable(GLUP_CLIPPING); 
        glupClipPlane(clip_eqn);
        glupPopMatrix();
    } else {
        clip_eqn[0] = 0.0;
        clip_eqn[1] = 0.0;
        clip_eqn[2] = 0.0;
        clip_eqn[3] = 0.0;
        glupClipPlane(clip_eqn);
        glupDisable(GLUP_CLIPPING);
    }

    /* translate the world of the distance between eye and center */
    glupTranslatef(0.0, 0.0, -zScreen);

    
    /* apply pan */
    glupTranslatef(cur_xlat[0], cur_xlat[1], cur_xlat[2]);

    /* if enabled, apply rotation */
    if(glup_viewer_is_enabled(GLUP_VIEWER_3D)) {
        build_rotmatrix(m, cur_rot);
        glupMultMatrixf(&m[0][0]);
    }

    /* restrict view to the Region Of Interest */
    glupScalef(1.0f / roi_radius, 1.0f / roi_radius, 1.0f / roi_radius);
    glupTranslatef(
        -0.5f * (xmin + xmax),
        -0.5f * (ymin + ymax),
        -0.5f * (zmin + zmax)
    );

    if(display_func != NULL) {
        glup_viewer_save_transform_for_picking();
        display_func();
    }

    glup_viewer_effect_end_frame();
    
    draw_foreground();
}

static void display() {
#ifndef __EMSCRIPTEN__
    double eye_dist =
        (double)(glup_viewer_get_float(GLUP_VIEWER_STEREOSCOPIC_EYE_DISTANCE));
#endif
    

#ifndef __EMSCRIPTEN__    
    if(glup_viewer_is_enabled(GLUP_VIEWER_STEREOSCOPIC_DISPLAY)) {
        /* left eye  */
        glDrawBuffer(GL_BACK_LEFT);
        glClear(GL_COLOR_BUFFER_BIT);
        actually_render_display(-eye_dist);

        /* right eye */
        glDrawBuffer(GL_BACK_RIGHT);
        glClear(GL_COLOR_BUFFER_BIT);
        actually_render_display(eye_dist);
    } else
#endif
        
    {
#ifndef __EMSCRIPTEN__        
        glDrawBuffer(GL_BACK); 
#endif
        glClear(GL_COLOR_BUFFER_BIT);
        actually_render_display(0.0);
    }
}

static void copy_image_to_clipboard(void);

void glup_viewer_char_callback(GLFWwindow* w, unsigned int c);

void glup_viewer_char_callback(GLFWwindow* w, unsigned int c) {

    glup_viewer_gui_char_callback(w,c);
    
    if(glup_viewer_gui_takes_input()) {
        return;
    }
    
    if(c == 3) {   /* 3 = <Ctrl> C */
        copy_image_to_clipboard();
        return;
    }
    if(keyboard_func != NULL) {
        keyboard_func((char) c);
    }
    if(toggle[c] != NULL) {
        *toggle[c] = !*toggle[c];
        glup_viewer_post_redisplay();
    }
    if(key_func[c] != NULL) {
        key_func[c]();
        glup_viewer_post_redisplay();
    }
    glup_viewer_post_redisplay();
}

static void toggle_clip() {
    glup_viewer_toggle(GLUP_VIEWER_CLIP);
    if(!glup_viewer_is_enabled(GLUP_VIEWER_CLIP)) {
        glup_viewer_disable(GLUP_VIEWER_EDIT_CLIP);
    }
}

static void toggle_edit_clip() {
    if(glup_viewer_is_enabled(GLUP_VIEWER_CLIP)) {
        glup_viewer_toggle(GLUP_VIEWER_EDIT_CLIP);
    } else {
        glup_viewer_disable(GLUP_VIEWER_EDIT_CLIP);
    }
}

static void toggle_fixed_clip() {
    if(glup_viewer_is_enabled(GLUP_VIEWER_CLIP)) {
        glup_viewer_toggle(GLUP_VIEWER_FIXED_CLIP);
        glup_viewer_disable(GLUP_VIEWER_EDIT_CLIP);
    }
}

static void glup_viewer_key_callback(
    GLFWwindow* w, int key, int scancode, int action, int mods
) {
    GLboolean handled = GL_FALSE;
    enum GlupViewerEvent ev;
    const char* keyname = NULL;

#ifdef __EMSCRIPTEN__    
    static char buffer[2];
#endif
    
    glup_viewer_gui_key_callback(w, key, scancode, action, mods);

    if(glup_viewer_gui_takes_input()) {
	/* 
	 * We continue capturing keypresses on function keys even if
	 * ImGUI window is active, else we cannot "run program" with
	 * F5 / "save file" with F2 in geocod !
	 */
	if(key < GLFW_KEY_F1 || key > GLFW_KEY_F12) {
	    return;
	}
    }
    

    if(keyboard_func_ext != NULL && action != GLFW_REPEAT) {
	ev = (action == GLFW_PRESS) ? GLUP_VIEWER_DOWN : GLUP_VIEWER_UP;

#ifdef __EMSCRIPTEN__
	/* Argh, glfwGetKeyName is not implemented in Emscripten */
	if(key < 128) {
	    buffer[1] = '\0';
	    buffer[0] = (char)(key);
	    keyname = buffer;
	}
#else	
	keyname = glfwGetKeyName(key,0);
#endif	
	if(keyname != NULL) {
	    handled = keyboard_func_ext(keyname,ev);	    
	} else {
	    if(key == GLFW_KEY_LEFT) {
		handled = keyboard_func_ext("left",ev);		
	    } else if(key == GLFW_KEY_RIGHT) {
		handled = keyboard_func_ext("right",ev);
	    } else if(key == GLFW_KEY_UP) {
		handled = keyboard_func_ext("up",ev);
	    } else if(key == GLFW_KEY_DOWN) {
		handled = keyboard_func_ext("down",ev);
	    } else if(key == GLFW_KEY_SPACE) {
		handled = keyboard_func_ext(" ",ev);		
	    } else if(key == GLFW_KEY_F1) {
		handled = keyboard_func_ext("F1",ev);				
	    } else if(key == GLFW_KEY_F2) {
		handled = keyboard_func_ext("F2",ev);				
	    } else if(key == GLFW_KEY_F3) {
		handled = keyboard_func_ext("F3",ev);				
	    } else if(key == GLFW_KEY_F4) {
		handled = keyboard_func_ext("F4",ev);				
	    } else if(key == GLFW_KEY_F5) {
		handled = keyboard_func_ext("F5",ev);				
	    } else if(key == GLFW_KEY_F6) {
		handled = keyboard_func_ext("F6",ev);				
	    } else if(key == GLFW_KEY_F7) {
		handled = keyboard_func_ext("F7",ev);				
	    } else if(key == GLFW_KEY_F18) {
		handled = keyboard_func_ext("F8",ev);				
	    } else if(key == GLFW_KEY_F9) {
		handled = keyboard_func_ext("F9",ev);				
	    } else if(key == GLFW_KEY_F10) {
		handled = keyboard_func_ext("F10",ev);				
	    } else if(key == GLFW_KEY_F11) {
		handled = keyboard_func_ext("F11",ev);				
	    } else if(key == GLFW_KEY_F12) {
		handled = keyboard_func_ext("F12",ev);				
	    } else if(key == GLFW_KEY_LEFT_CONTROL) {
		handled = keyboard_func_ext("left_control",ev);
	    } else if(key == GLFW_KEY_RIGHT_CONTROL) {
		handled = keyboard_func_ext("right_control",ev);
	    } else if(key == GLFW_KEY_LEFT_ALT) {
		handled = keyboard_func_ext("left_alt",ev);
	    } else if(key == GLFW_KEY_RIGHT_ALT) {
		handled = keyboard_func_ext("right_alt",ev);
	    } else if(key == GLFW_KEY_LEFT_SHIFT) {
		handled = keyboard_func_ext("left_shift",ev);
	    } else if(key == GLFW_KEY_RIGHT_SHIFT) {
		handled = keyboard_func_ext("right_shift",ev);
	    } else if(key == GLFW_KEY_ESCAPE) {
		handled = keyboard_func_ext("escape",ev);		
	    } else if(key == GLFW_KEY_TAB) {
		handled = keyboard_func_ext("tab",ev);		
	    } else if(key == GLFW_KEY_BACKSPACE) {
		handled = keyboard_func_ext("backspace",ev);		
	    }
	}
	if(handled) {
	    glup_viewer_post_redisplay();
	    return;
	}
    }
    
    if(action != GLFW_PRESS) {
        return;
    }
    
    switch(key) {
        case GLFW_KEY_F1:
            toggle_clip();
            break;
        case GLFW_KEY_F2:
            toggle_edit_clip();
            break;
        case GLFW_KEY_F3:
            toggle_fixed_clip();
            break;
        case GLFW_KEY_F4:
            glup_viewer_toggle(GLUP_VIEWER_SHOW_CLIP);
            break;
    }
    
    glup_viewer_post_redisplay();
}


void glup_viewer_home(void) {
    cur_xlat[0] = cur_xlat[1] = cur_xlat[2] = 0.0;
    params[GLUP_VIEWER_ZOOM] = 1.0;
    trackball(cur_rot, 0.0, 0.0, 0.0, 0.0);
    trackball(cur_rot_light, 0.0, 0.0, 0.0, 0.0);
    trackball(cur_rot_clip, 0.0, 0.0, 0.0, 0.0);
}

static void init() {

#ifndef __EMSCRIPTEN__    
    if(!gladLoadGL()) {
        printf("GLAD: could not load OpenGL\n");
        exit(-1);
    }
#endif

    if(glupCurrentContext() == NULL) {
        glupMakeCurrent(glupCreateContext());
    }

    if(glupCurrentContext() == NULL) {
	exit(-1);
    }
    
    glup_viewer_gui_init(glup_viewer_window);
    
    glup_viewer_disable(GLUP_VIEWER_IDLE_REDRAW);
    glup_viewer_enable(GLUP_VIEWER_DRAW_SCENE);

#ifndef __EMSCRIPTEN__    
    glup_viewer_add_key_func('q', glup_viewer_exit_main_loop, "quit");
    glup_viewer_add_key_func(27, glup_viewer_exit_main_loop, "quit");
#endif
    
    glup_viewer_add_key_func('H', glup_viewer_home, "home");

    glup_viewer_add_toggle('h', &caps[GLUP_VIEWER_SHOW_HELP], "help");
    glup_viewer_add_toggle(
        'l', &caps[GLUP_VIEWER_ROTATE_LIGHT], "light rotation"
    );

    glupEnable(GLUP_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    
    glPolygonOffset(1.0, 2.0);
    glEnable(GL_POLYGON_OFFSET_FILL);
    
    glup_viewer_home();
}

void glup_viewer_set_window_title(const char* s) {
    title = s;
}

/***************************************************************/

void glup_viewer_set_display_func(GlupViewerDisplayFunc f) {
    display_func = f;
}

void glup_viewer_set_overlay_func(GlupViewerDisplayFunc f) {
    overlay_func = f;
}

void glup_viewer_set_keyboard_func(GlupViewerKeyboardFunc f) {
    keyboard_func = f;
}

void glup_viewer_set_keyboard_func_ext(GlupViewerKeyboardFuncExt f) {
    keyboard_func_ext = f;
}

void glup_viewer_set_mouse_func(GlupViewerMouseFunc f) {
    mouse_func = f;
}

void glup_viewer_set_init_func(GlupViewerInitFunc f) {
    init_func = f;
}

void glup_viewer_set_drag_drop_func(GlupViewerDragDropFunc f) {
    drag_drop_func = f;
}

GlupViewerDisplayFunc glup_viewer_get_display_func() {
    return display_func;
}

GlupViewerDisplayFunc glup_viewer_get_overlay_func() {
    return overlay_func;
}

GlupViewerKeyboardFunc glup_viewer_get_keyboard_func() {
    return keyboard_func;
}

GlupViewerMouseFunc glup_viewer_get_mouse_func() {
    return mouse_func;
}

GlupViewerInitFunc glup_viewer_get_init_func() {
    return init_func;
}

GlupViewerDragDropFunc glup_viewer_get_drag_drop_func() {
    return drag_drop_func;
}

/***************************************************************/

static void init_keys_if_needed() {
    static int first = 1;
    if(first) {
        int i;
        first = 0;
        for(i = 0; i < 256; i++) {
            toggle[i] = NULL;
            key_description[i] = NULL;
            key_func[i] = NULL;
        }
    }
}

static void set_key_description(char key, const char* description) {
    if(key_description[(int) key] != NULL) {
        printf(
            "Warning: overriding key %c with %s (was %s)\n",
            key, description ? description : "NULL", key_description[(int) key]
        );
    }
    free(key_description[(int) key]);
    if(description == NULL) {
        key_description[(int) key] = NULL;
    } else {
        key_description[(int) key] = malloc(strlen(description) + 1);
        strcpy(key_description[(int) key], description);
    }
}

void glup_viewer_add_toggle(
    char key, GLboolean* pointer, const char* description
) {
    init_keys_if_needed();
    toggle[(int) key] = pointer;
    key_func[(int) key] = NULL;
    set_key_description(key, description);
}

void glup_viewer_add_key_func(
    char key, GlupViewerKeyFunc f, const char* description
) {
    init_keys_if_needed();
    key_func[(int) key] = f;
    toggle[(int) key] = NULL;
    set_key_description(key, description);
}

void glup_viewer_unbind_key(char key) {
    key_func[(int) key] = NULL;
    toggle[(int) key] = NULL;
    set_key_description(key, NULL);
}

static void cleanup(void) {
    int i;
    glup_viewer_gui_cleanup();
    for(i = 0; i < 256; ++i) {
        if(key_description[i] != NULL) {
            free(key_description[i]);
        }
    }
}

void glup_viewer_exit_main_loop() {
    in_main_loop_ = GL_FALSE;
}


/* exported, used in glup_viewer_gui.cpp */
void glup_viewer_one_frame(void);

void glup_viewer_one_frame() {
    int cur_width;
    int cur_height;
    GlupViewerInitFunc init = init_func;

    /*
       can happen when ImGUI Graphite application 
       triggers update too soon.
    */
    if(glup_viewer_window == NULL) {
	return;
    }
    
    if(init != NULL) {
	init_func = NULL;
        init();
        glup_viewer_start_time = now();
    }

    if(!glfwWindowShouldClose(glup_viewer_window) && in_main_loop_) {
        glfwGetFramebufferSize(glup_viewer_window, &cur_width, &cur_height);
        if(glup_viewer_W != cur_width || glup_viewer_H != cur_height) {
            reshape(cur_width, cur_height);
        }
        glfwPollEvents();
        if(
            glup_viewer_is_enabled(GLUP_VIEWER_IDLE_REDRAW) ||
            glup_viewer_needs_redraw > 0 ||
            (now()-glup_viewer_start_time) < 1.0
            /* overcomes missing update at startup */
        ) {
            if(glup_viewer_is_enabled(GLUP_VIEWER_TWEAKBARS)) {
                glup_viewer_gui_begin_frame();
                display();
                glup_viewer_gui_end_frame();
            } else {
                display();
            }
            if(glup_viewer_needs_redraw > 0) {
                --glup_viewer_needs_redraw;
            }
        } else {
            glup_viewer_pause();
        }
#ifdef __EMSCRIPTEN__
	/* 
	 * Set alpha channel to 1 else the image is composited 
	 * with the background 
	 */
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_TRUE);
	glClear(GL_COLOR_BUFFER_BIT);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
#endif	
        glfwSwapBuffers(glup_viewer_window);
    } 
}

void glup_viewer_main_loop(int argc, char** argv) {
    (void)argc; /* suppresses a warning */
    (void)argv; /* suppresses a warning */
    
    if(!glfwInit()) {
        fprintf(stderr, "Could not initialize GLFW\n");
        exit(-1);
    }
    
    in_main_loop_ = GL_TRUE;

    if(glup_viewer_test_arg_string("gfx:GL_profile", "core")) {
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    }

    
    if(glup_viewer_get_arg_bool("gfx:GL_debug")) {
        glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);            
    }

    if(glup_viewer_is_high_dpi()) {
        glup_viewer_W *= 2;
        glup_viewer_H *= 2;
    }

    glup_viewer_set_screen_size_from_args();

    transparent = glup_viewer_get_arg_bool("gfx:transparent");
    
    if(transparent) {
/*	glfwWindowHint(GLFW_ALPHA_BITS,8);  */
#ifdef GLFW_ALPHA_MASK	    
	glfwWindowHint(GLFW_ALPHA_MASK,GL_TRUE);
#else
	fprintf(stderr,"Transparent not supported by this version of GLFW\n");
#endif	    
    }
   
    if( 
        glup_viewer_get_arg_bool("gfx:fullscreen") ||
        glup_viewer_is_enabled(GLUP_VIEWER_FULL_SCREEN) 
    ) {
        GLFWmonitor* monitor = glfwGetPrimaryMonitor();
        const GLFWvidmode* vidmode = glfwGetVideoMode(monitor);
        glup_viewer_W = vidmode->width;
        glup_viewer_H = vidmode->height;

	if(glup_viewer_get_arg_bool("gfx:no_decoration")) {
	    glfwWindowHint(GLFW_FOCUSED,GL_TRUE);	
	    glfwWindowHint(GLFW_DECORATED,GL_FALSE);
	    glfwWindowHint(GLFW_RESIZABLE,GL_FALSE);
	    glfwWindowHint(GLFW_AUTO_ICONIFY,GL_FALSE);
	    glfwWindowHint(GLFW_FLOATING,GL_FALSE);
	    glfwWindowHint(GLFW_MAXIMIZED,GL_TRUE);
	}				
	
        glup_viewer_window = glfwCreateWindow(
            glup_viewer_W, glup_viewer_H, title,
	    (
		glup_viewer_get_arg_bool("gfx:no_decoration") ?
		glfwGetPrimaryMonitor() : NULL
	    ), 
	    NULL
        );
    } else {
        glup_viewer_window = glfwCreateWindow(
            glup_viewer_W, glup_viewer_H, title, NULL, NULL
        );
    }
   

    if(glup_viewer_window == NULL) {
        fprintf(stderr, "Could not create GLFW window\n");
        exit(-1);
    }

    glfwMakeContextCurrent(glup_viewer_window);
    glfwSwapInterval(1);
    

    glfwSetMouseButtonCallback(glup_viewer_window,mouse);
    glfwSetCursorPosCallback(glup_viewer_window,passive_mouse);
    glfwSetCharCallback(glup_viewer_window,glup_viewer_char_callback);
    glfwSetKeyCallback(glup_viewer_window,glup_viewer_key_callback);

#ifndef __EMSCRIPTEN__    
    glfwSetDropCallback(glup_viewer_window,drag_drop);
#endif
    
    glfwSetScrollCallback(glup_viewer_window,scroll);
    
    init();
    glfwGetFramebufferSize(glup_viewer_window, &glup_viewer_W, &glup_viewer_H);    
    reshape(glup_viewer_W, glup_viewer_H);

    atexit(cleanup);

    glup_viewer_start_time = now();
   
#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(glup_viewer_one_frame, 0, 1);
#else
    while (!glfwWindowShouldClose(glup_viewer_window) && in_main_loop_) {
        glup_viewer_one_frame();
    }
#endif
    
    glfwDestroyWindow(glup_viewer_window);
    glup_viewer_window = NULL;
    glfwTerminate();
    
    in_main_loop_ = GL_FALSE;
}

void glup_viewer_set_region_of_interest(
    float xm, float ym, float zm, float xM, float yM, float zM
) {
    xmin = xm;
    ymin = ym;
    zmin = zm;
    xmax = xM;
    ymax = yM;
    zmax = zM;
    roi_radius = (float) sqrt((double)(
        0.25f * (xmax - xmin) * (xmax - xmin) +
        0.25f * (ymax - ymin) * (ymax - ymin) +
        0.25f * (zmax - zmin) * (zmax - zmin)
    ));
}

void glup_viewer_get_region_of_interest(
    float* xm, float* ym, float* zm, float* xM, float* yM, float* zM
) {
    *xm = xmin;
    *ym = ymin;
    *zm = zmin;
    *xM = xmax;
    *yM = ymax;
    *zM = zmax;
}

void glup_viewer_set_background(int tex) {
    background_tex = (GLuint) tex;
}

int glup_viewer_get_background() {
    return (int) background_tex;
}

/* =============================== Trackball code ========================= */

/*
 * (c) Copyright 1993, 1994, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED
 * Permission to use, copy, modify, and distribute this software for
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.
 *
 * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * US Government Users Restricted Rights
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */
/*
 * Trackball code:
 *
 * Implementation of a virtual trackball.
 * Implemented by Gavin Bell, lots of ideas from Thant Tessman and
 *   the August '88 issue of Siggraph's "Computer Graphics," pp. 121-129.
 *
 * Vector manip code:
 *
 * Original code from:
 * David M. Ciemiewicz, Mark Grossman, Henry Moreton, and Paul Haeberli
 *
 * Much mucking with by:
 * Gavin Bell
 */

/*
 * This size should really be based on the distance from the center of
 * rotation to the point on the object underneath the mouse.  That
 * point would then track the mouse as closely as possible.  This is a
 * simple example, though, so that is left as an Exercise for the
 * Programmer.
 */
#define TRACKBALLSIZE (0.8f)

/*
 * Local function prototypes (not defined in trackball.h)
 */
static float tb_project_to_sphere(float, float, float);
static void normalize_quat(float[4]);

static void vzero(float* v) {
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
}

static void vset(float* v, float x, float y, float z) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

static void vsub(const float* src1, const float* src2, float* dst) {
    dst[0] = src1[0] - src2[0];
    dst[1] = src1[1] - src2[1];
    dst[2] = src1[2] - src2[2];
}

static void vcopy(const float* v1, float* v2) {
    int i;
    for(i = 0; i < 3; i++) {
        v2[i] = v1[i];
    }
}

static void vcross(const float* v1, const float* v2, float* cross) {
    float temp[3];

    temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
    temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
    temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
    vcopy(temp, cross);
}

static float vlength(const float* v) {
    return (float) sqrt((double)(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
}

static void vscale(float* v, float div) {
    v[0] *= div;
    v[1] *= div;
    v[2] *= div;
}

static void vnormal(float* v) {
    vscale(v, 1.0f / vlength(v));
}

static float vdot(const float* v1, const float* v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

static void vadd(const float* src1, const float* src2, float* dst) {
    dst[0] = src1[0] + src2[0];
    dst[1] = src1[1] + src2[1];
    dst[2] = src1[2] + src2[2];
}

/*
 * Ok, simulate a track-ball.  Project the points onto the virtual
 * trackball, then figure out the axis of rotation, which is the cross
 * product of P1 P2 and O P1 (O is the center of the ball, 0,0,0)
 * Note:  This is a deformed trackball-- is a trackball in the center,
 * but is deformed into a hyperbolic sheet of rotation away from the
 * center.  This particular function was chosen after trying out
 * several variations.
 *
 * It is assumed that the arguments to this routine are in the range
 * (-1.0 ... 1.0)
 */
static void trackball(float q[4], float p1x, float p1y, float p2x, float p2y) {
    float a[3]; /* Axis of rotation */
    float phi;  /* how much to rotate about axis */
    float p1[3], p2[3], d[3];
    float t;

    if(p1x == p2x && p1y == p2y) {
        /* Zero rotation */
        vzero(q);
        q[3] = 1.0;
        return;
    }

    /*
     * First, figure out z-coordinates for projection of P1 and P2 to
     * deformed sphere
     */
    vset(p1, p1x, p1y, tb_project_to_sphere(TRACKBALLSIZE, p1x, p1y));
    vset(p2, p2x, p2y, tb_project_to_sphere(TRACKBALLSIZE, p2x, p2y));

    /*
     *  Now, we want the cross product of P1 and P2
     */
    vcross(p2, p1, a);

    /*
     *  Figure out how much to rotate around that axis.
     */
    vsub(p1, p2, d);
    t = vlength(d) / (2.0f * TRACKBALLSIZE);

    /*
     * Avoid problems with out-of-control values...
     */
    if(t > 1.0f) {
        t = 1.0;
    }
    if(t < -1.0f) {
        t = -1.0;
    }
    phi = 2.0f * (float) asin((double)t);

    axis_to_quat(a, phi, q);
}

/*
 *  Given an axis and angle, compute quaternion.
 */
void axis_to_quat(float a[3], float phi, float q[4])
{
    vnormal(a);
    vcopy(a, q);
    vscale(q, (float) sin(((double)phi) / 2.0));
    q[3] = (float) cos(((double)phi) / 2.0);
}

void glup_viewer_set_scene_rotation(
    float xyz[3], float angle
) {
    axis_to_quat(xyz, (float)((double)angle * M_PI / 180.0), cur_rot);
}

/*
 * Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
 * if we are away from the center of the sphere.
 */
static float tb_project_to_sphere(float r, float x, float y) {
    float d, z;

    d = (float) sqrt((double)(x * x + y * y));
    if(d < r * 0.70710678118654752440f) {     /* Inside sphere */
        z = (float) sqrt((double)(r * r - d * d));
    } else {           /* On hyperbola */
        float t = r / 1.41421356237309504880f;
        z = t * t / d;
    }
    return z;
}

/*
 * Given two rotations, e1 and e2, expressed as quaternion rotations,
 * figure out the equivalent single rotation and stuff it into dest.
 *
 * NOTE: This routine is written so that q1 or q2 may be the same
 * as dest (or each other).
 */

void
add_quats(float q1[4], float q2[4], float dest[4])
{
    float t1[4], t2[4], t3[4];
    float tf[4];

    vcopy(q1, t1);
    vscale(t1, q2[3]);

    vcopy(q2, t2);
    vscale(t2, q1[3]);

    vcross(q2, q1, t3);
    vadd(t1, t2, tf);
    vadd(t3, tf, tf);
    tf[3] = q1[3] * q2[3] - vdot(q1, q2);

    dest[0] = tf[0];
    dest[1] = tf[1];
    dest[2] = tf[2];
    dest[3] = tf[3];

    normalize_quat(dest);
}

/*
 * Quaternions always obey:  a^2 + b^2 + c^2 + d^2 = 1.0
 * If they don't add up to 1.0, dividing by their magnitued will
 * renormalize them.
 *
 * Note: See the following for more information on quaternions:
 *
 * - Shoemake, K., Animating rotation with quaternion curves, Computer
 *   Graphics 19, No 3 (Proc. SIGGRAPH'85), 245-254, 1985.
 * - Pletinckx, D., Quaternion calculus as a basic tool in computer
 *   graphics, The Visual Computer 5, 2-13, 1989.
 */
static void
normalize_quat(float q[4])
{
    int i;
    float mag;

    mag = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    for(i = 0; i < 4; i++) {
        q[i] /= mag;
    }
}

/*
 * Build a rotation matrix, given a quaternion rotation.
 *
 */
void
build_rotmatrix(float m[4][4], float q[4])
{
    m[0][0] = 1.0f - 2.0f * (q[1] * q[1] + q[2] * q[2]);
    m[0][1] = 2.0f * (q[0] * q[1] - q[2] * q[3]);
    m[0][2] = 2.0f * (q[2] * q[0] + q[1] * q[3]);
    m[0][3] = 0.0f;

    m[1][0] = 2.0f * (q[0] * q[1] + q[2] * q[3]);
    m[1][1] = 1.0f - 2.0f * (q[2] * q[2] + q[0] * q[0]);
    m[1][2] = 2.0f * (q[1] * q[2] - q[0] * q[3]);
    m[1][3] = 0.0f;

    m[2][0] = 2.0f * (q[2] * q[0] - q[1] * q[3]);
    m[2][1] = 2.0f * (q[1] * q[2] + q[0] * q[3]);
    m[2][2] = 1.0f - 2.0f * (q[1] * q[1] + q[0] * q[0]);
    m[2][3] = 0.0f;

    m[3][0] = 0.0f;
    m[3][1] = 0.0f;
    m[3][2] = 0.0f;
    m[3][3] = 1.0f;
}

/*---------------------------------------------------------------------------*/

#ifdef WIN32

static void copy_image_to_clipboard() {
    int w = glup_viewer_W;
    int h = glup_viewer_H;
    /*    int row_len = w * 3 ;  */
    int image_size, size;
    HANDLE handle;
    char* base_mem;
    char* pData;
    BITMAPINFOHEADER* header;

    /* Step 1: Try to open Window's clipboard
     *   NULL -> bind to current process
     */
    if(!OpenClipboard(NULL)) {
        return;
    }

    /* Step 2: create a shared memory segment, with
     * a DIB (Device Independent Bitmap) in it.
     */

    image_size = 3 * w * h;
    size = sizeof(BITMAPINFOHEADER) + image_size;

    handle = (HANDLE) GlobalAlloc(GHND, size);
    if(handle == NULL) {
        return;
    }

    pData = GlobalLock((HGLOBAL) handle);
    header = (BITMAPINFOHEADER*) pData;

    header->biSize = sizeof(BITMAPINFOHEADER);
    header->biWidth = w;
    header->biHeight = h;
    header->biPlanes = 1;
    header->biBitCount = 24;
    header->biCompression = BI_RGB;
    header->biSizeImage = 0;
    header->biXPelsPerMeter = 1000000;
    header->biYPelsPerMeter = 1000000;
    header->biClrUsed = 0;
    header->biClrImportant = 0;

    base_mem = pData + sizeof(BITMAPINFOHEADER);
    glReadPixels(0, 0, w, h, GL_BGR_EXT, GL_UNSIGNED_BYTE, base_mem);

    /* Step 3: Detach the clipboard from current process. */
    GlobalUnlock((HGLOBAL) handle);
    EmptyClipboard();
    SetClipboardData(CF_DIB, handle);
    CloseClipboard();
}

#else
static void copy_image_to_clipboard() {
}

#endif

/* -------------------------------------------------------------------------- */

static int htoi(char digit) {
    if(digit >= '0' && digit <= '9') {
        return digit - '0';
    }
    if(digit >= 'a' && digit <= 'f') {
        return digit - 'a' + 10;
    }
    if(digit >= 'A' && digit <= 'F') {
        return digit - 'A' + 10;
    }
    fprintf(stderr, "xpm: unknown digit\n");
    return 0;
}

/* The colormap. */
static unsigned char i2r[1024];
static unsigned char i2g[1024];
static unsigned char i2b[1024];
static unsigned char i2a[1024];

/* Converts a two-digit XPM color code into
 *  a color index.
 */
static int char_to_index[256][256];

void glTexImage2DXPM(char const* const* xpm_data) {
    int width, height, nb_colors, chars_per_pixel;
    int line = 0;
    int color = 0;
    int key1 = 0, key2 = 0;
    char* colorcode;
    int x, y;
    unsigned char* rgba;
    unsigned char* pixel;

    sscanf(
        xpm_data[line], "%6d%6d%6d%6d",
        &width, &height, &nb_colors, &chars_per_pixel
    );
    line++;
    if(nb_colors > 1024) {
        fprintf(stderr, "xpm with more than 1024 colors\n");
        return;
    }
    if(chars_per_pixel != 1 && chars_per_pixel != 2) {
        fprintf(stderr, "xpm with more than 2 chars per pixel\n");
        return;
    }
    for(color = 0; color < nb_colors; color++) {
        int r, g, b;
        int none ;
        
        key1 = xpm_data[line][0];
        key2 = (chars_per_pixel == 2) ? xpm_data[line][1] : 0;
        colorcode = strstr(xpm_data[line], "c #");
        none = 0;
        if(colorcode == NULL) {
            colorcode = "c #000000";
            if(strstr(xpm_data[line], "None") != NULL) {
                none = 1;
            } else {
                fprintf(
                    stderr, "unknown xpm color entry (replaced with black)\n"
                );
            }
        }
        colorcode += 3;

	if(strlen(colorcode) == 12) {
	    r = 16 * htoi(colorcode[0]) + htoi(colorcode[1]);
	    g = 16 * htoi(colorcode[4]) + htoi(colorcode[5]);
	    b = 16 * htoi(colorcode[8]) + htoi(colorcode[9]);
	} else {
	    r = 16 * htoi(colorcode[0]) + htoi(colorcode[1]);
	    g = 16 * htoi(colorcode[2]) + htoi(colorcode[3]);
	    b = 16 * htoi(colorcode[4]) + htoi(colorcode[5]);
	}
	    
        i2r[color] = (unsigned char) r;
        i2g[color] = (unsigned char) g;
        i2b[color] = (unsigned char) b;
	if(none) {
	    i2a[color] = 0;
	} else {
	    i2a[color] = 255;
	}
        char_to_index[key1][key2] = color;
        line++;
    }
    rgba = (unsigned char*) malloc((size_t) (width * height * 4));
    pixel = rgba;
    for(y = 0; y < height; y++) {
        for(x = 0; x < width; x++) {
            if(chars_per_pixel == 2) {
                key1 = xpm_data[line][2 * x];
                key2 = xpm_data[line][2 * x + 1];
            } else {
                key1 = xpm_data[line][x];
                key2 = 0;
            }
            color = char_to_index[key1][key2];
            pixel[0] = i2r[color];
            pixel[1] = i2g[color];
            pixel[2] = i2b[color];
            pixel[3] = i2a[color];
            pixel += 4;
        }
        line++;
    }
    
    glTexImage2D(
        GL_TEXTURE_2D, 0,
        GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba
    );
#ifndef __EMSCRIPTEN__    
    glGenerateMipmap(GL_TEXTURE_2D);
#endif    
    free(rgba);
}

/*----------------------------------------------------------------------------*/

#ifdef WITH_PNG
static GLboolean glup_viewer_load_image_PNG(
    const char* filename,
    GLuint* width_in, GLuint* height_in, GLuint* bpp, GLvoid** pixels_in
) {

    FILE* in;
    png_structp png_ptr;
    png_infop info_ptr;
    unsigned long width, height;
    int bit_depth, color_type, interlace_type;
    png_bytep row_pointer;
    unsigned int row;
    char** pixels = (char**) pixels_in;
    png_uint_32 w, h;

    in = fopen(filename, "rb");
    if(in == NULL) {
        fprintf(stderr, "PNG loader: failed opening file %s\n", filename);
        return GL_FALSE;
    }

    png_ptr = png_create_read_struct(
        PNG_LIBPNG_VER_STRING,
        (png_voidp) NULL, (png_error_ptr) NULL, (png_error_ptr) NULL
    );
    if(png_ptr == NULL) {
        fclose(in);
        fprintf(stderr, "PNG loader: error while loading %s\n", filename);
        return GL_FALSE;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if(info_ptr == NULL) {
        png_destroy_read_struct(
            &png_ptr, (png_infopp) NULL, (png_infopp) NULL
        );
        fprintf(stderr, "PNG loader: error while loading %s\n", filename);
        fclose(in);
        return GL_FALSE;
    }

    png_init_io(png_ptr, in);

    /* read header  */
    png_read_info(png_ptr, info_ptr);
    png_get_IHDR(
        png_ptr, info_ptr, &w, &h, &bit_depth, &color_type,
        &interlace_type, NULL, NULL
    );
    width = w;
    height = h;

    if(color_type == PNG_COLOR_TYPE_GRAY) {
        *bpp = 1;
    } else if(color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
        *bpp = 4;
    } else {
        *bpp = 3;
    }

    if(color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_expand(png_ptr);
    }

    *pixels = malloc(width * height * *bpp);
    if(*pixels == NULL) {
        png_destroy_read_struct(
            &png_ptr, (png_infopp) NULL, (png_infopp) NULL
        );
        fprintf(
            stderr,
    "PNG loader: error while loading %s, unable to allocate sufficient memory\n"
            , filename
        );
        fclose(in);
        return GL_FALSE;
    }

    /* Read the image one line at a time */
    row_pointer = (png_bytep) malloc(width * *bpp);

    for(row = 0; row < height; row++) {
        png_read_rows(png_ptr, &row_pointer, NULL, 1);
        memcpy(
            *pixels + (height - 1 - row) * width * *bpp,
            row_pointer,
            width * *bpp
        );
    }
    free(row_pointer);

    png_read_end(png_ptr, info_ptr);
    png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp) NULL);

    *width_in = (GLuint) width;
    *height_in = (GLuint) height;

    fclose(in);
    return GL_TRUE;
}

#endif

/*----------------------------------------------------------------------------*/

#define GLUP_VIEWER_BMP_RLE_8 1
#define GLUP_VIEWER_BMP_RLE_4 2

typedef unsigned int bmp_int32;
typedef unsigned short bmp_int16;

#define readH(x) fread(&(header.x), sizeof(header.x), 1, f);

typedef struct {
    bmp_int16 sType;
    bmp_int32 iSizeOfFile;
    bmp_int16 sReserved1;
    bmp_int16 sReserved2;
    bmp_int32 iOffBits;
    bmp_int32 iSize;
    bmp_int32 iWidth;
    bmp_int32 iHeight;
    bmp_int16 sPlanes;
    bmp_int16 sBitCount;
    bmp_int32 iCompression;
    bmp_int32 iSizeImage;
    bmp_int32 iXpelsPerMeter;
    bmp_int32 iYpelsPerMeter;
    bmp_int32 iClrUsed;
    bmp_int32 iClrImportant;
} GlupViewerBMPHeader;

static void rgb_to_bgr(
    GLuint width, GLuint height, GLuint bpp, GLvoid* pixels
) {
    char* p = (char*) pixels;
    int i;
    if(bpp != 3 && bpp != 4) {
        return;
    }
    for(i = 0; i < (int) (width * height); i++) {
        char tmp = p[0];
        p[0] = p[2];
        p[2] = tmp;
        p += bpp;
    }
}

static GLboolean glup_viewer_load_image_BMP(
    const char* filename,
    GLuint* width, GLuint* height, GLuint* bpp, GLvoid** pixels_in
) {
    GlupViewerBMPHeader header;
    char padding[3];
    size_t rowlen;
    int i;
    char* pixels;
    FILE* f = fopen(filename, "rb");
    int ok = 1;

    if(f == NULL) {
        fprintf(stderr, "cannot open %s\n", filename);
        return GL_FALSE;
    }


   
    /*
     * Argh, I need to read the header field by field, seems
     * that the compiler aligns the fields in a way that makes
     * the memory layout mismatch the file layout.
     */
   
    memset(&header, 0, sizeof(GlupViewerBMPHeader));
    ok = ok && readH(sType);
    ok = ok && readH(iSizeOfFile);
    ok = ok && readH(sReserved1);
    ok = ok && readH(sReserved2);
    ok = ok && readH(iOffBits);
    ok = ok && readH(iSize);
    ok = ok && readH(iWidth);
    ok = ok && readH(iHeight);
    ok = ok && readH(sPlanes);
    ok = ok && readH(sBitCount);
    ok = ok && readH(iCompression);
    ok = ok && readH(iSizeImage);
    ok = ok && readH(iXpelsPerMeter);
    ok = ok && readH(iYpelsPerMeter);
    ok = ok && readH(iClrUsed);
    ok = ok && readH(iClrImportant);

    if(!ok) {
        fprintf(stderr, "Error while reading BMP header of %s\n", filename);
        fclose(f);
        return GL_FALSE;
    }

    if(
        (header.iCompression & GLUP_VIEWER_BMP_RLE_8) ||
        (header.iCompression & GLUP_VIEWER_BMP_RLE_4)
    ) {
        fprintf(
            stderr,
            "BMP loader: cannot load compressed BMP files, sorry\n"
        );
        fclose(f);
        return GL_FALSE;
    }
    *width = header.iWidth;
    *height = header.iHeight;
    *bpp = header.sBitCount / 8;

    pixels = malloc(*width * *height * *bpp);

    rowlen = (size_t) (*width * *bpp);
    for(i = 0; i < (int) *height; i++) {
        ok = ok && fread(pixels + i * (int) rowlen, rowlen, (size_t) 1, f);
        /* BMP rows are padded to 4-bytes multiples */
        ok = ok && fread(padding, rowlen % 4, (size_t) 1, f);
    }

    rgb_to_bgr(*width, *height, *bpp, pixels);
    *pixels_in = pixels;

    fclose(f);

    if(!ok) {
        fprintf(stderr, "Error while reading BMP data from %s\n", filename);
        free(pixels);
        return GL_FALSE;
    }

    return GL_TRUE;
}

/*----------------------------------------------------------------------------*/

static const char* extension(const char* filename) {
    const char* point_loc = strrchr(filename, '.');
    return (point_loc == NULL) ? NULL : (point_loc + 1);
}

GLboolean glup_viewer_load_image(
    const char* filename, GLuint* width, GLuint* height,
    GLuint* bpp, GLvoid** pixels
) {
    const char* ext = extension(filename);
#ifdef WITH_PNG
    if(!strcmp(ext, "png") || !strcmp(ext, "PNG")) {
        return glup_viewer_load_image_PNG(filename, width, height, bpp, pixels);
    }
#endif
    if(!strcmp(ext, "bmp") || !strcmp(ext, "BMP")) {
        return glup_viewer_load_image_BMP(filename, width, height, bpp, pixels);
    }
    return GL_FALSE;
}

void glTexImage2Dfile(const char* filename) {
    GLuint width, height, bpp;
    GLvoid* pixels;
    if(!glup_viewer_load_image(filename, &width, &height, &bpp, &pixels)) {
        return;
    }
    switch(bpp) {
        case 1:
            glTexImage2D(
                GL_TEXTURE_2D, 0, GL_LUMINANCE,
                (GLsizei) width, (GLsizei) height, 0,
                GL_LUMINANCE, GL_UNSIGNED_BYTE, pixels
            );
            break;
        case 3:
            glTexImage2D(
                GL_TEXTURE_2D, 0, GL_RGB,
                (GLsizei) width, (GLsizei) height, 0,
                GL_RGB, GL_UNSIGNED_BYTE, pixels
            );
            break;
        case 4:
            glTexImage2D(
                GL_TEXTURE_2D, 0, GL_RGBA,
                (GLsizei) width, (GLsizei) height, 0,
                GL_RGBA, GL_UNSIGNED_BYTE, pixels
            );
            break;
        default:
            glup_viewer_assert_not_reached();
    }
#ifndef __EMSCRIPTEN__    
    glGenerateMipmap(GL_TEXTURE_2D);
#endif    
    free(pixels);
}


#define c1 0.35 
#define c2 0.5 
#define c3 1.0 

static double color_table[12][3] = {
    {c3, c2, c2},
    {c2, c3, c2},
    {c2, c2, c3},
    {c2, c3, c3},
    {c3, c2, c3},
    {c3, c3, c2},
    
    {c1, c2, c2},
    {c2, c1, c2},
    {c2, c2, c1},
    {c2, c1, c1},
    {c1, c2, c1},
    {c1, c1, c2}
};

static int random_color_index_ = 0 ;


static float white[4] = {
    0.8f, 0.8f, 0.3f, 1.0f
};

void glup_viewer_random_color() {
    glupColor3d(
        color_table[random_color_index_][0], 
        color_table[random_color_index_][1], 
        color_table[random_color_index_][2]
    );
    random_color_index_ = (random_color_index_ + 1) % 12 ;
}
    
void glup_viewer_randomize_colors(int index) {
    random_color_index_ = (index % 12) ;
}

void glup_viewer_random_color_from_index(int index) {
    if(index >= 0) {
        glup_viewer_randomize_colors(index) ;
        glup_viewer_random_color() ;
    } else {
        glupColor4fv(white) ;
    }
}

GLboolean glup_viewer_is_high_dpi(void) {
#ifdef __EMSCRIPTEN__
    return GL_FALSE;
#else    
    GLFWmonitor* monitor = glfwGetPrimaryMonitor();
    const GLFWvidmode* vidmode = glfwGetVideoMode(monitor);
    return (vidmode->width > 1920);
#endif    
}

void glup_viewer_snapshot(const char* filename) {
    FILE* f = fopen(filename, "wb");
    int line_size = glup_viewer_W*3;
    void* data = NULL;
    char* line1;
    char* line2;
    int y;
    int xx;
    char tmp;
    
    if(f == NULL) {
	return;
    }

    /* Read OpenGL frame buffer */
    data = malloc((size_t)(glup_viewer_W*glup_viewer_H*3));
    glReadPixels(
	0, 0, glup_viewer_W, glup_viewer_H,
	GL_RGB, GL_UNSIGNED_BYTE, data
    );

    /* Flip image along vertical axis */
    for(y=0; y<glup_viewer_W/2; ++y) {
	line1 = (char*)data + line_size*y;
	line2 = (char*)data + line_size*(glup_viewer_H-1-y);
	for(xx=0; xx<line_size; ++xx) {
	    tmp = line1[xx];
	    line1[xx] = line2[xx];
	    line2[xx] = tmp;
	}
    }

    /* Save as PPM */
    fprintf(f,"P6\n%d %d\n255\n",glup_viewer_W,glup_viewer_H);
    fwrite(data,1,(size_t)(glup_viewer_W*glup_viewer_H*3),f);
    
    free(data);
    fclose(f);
}

