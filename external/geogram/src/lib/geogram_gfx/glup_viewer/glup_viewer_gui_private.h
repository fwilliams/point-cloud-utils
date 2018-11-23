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


#ifndef GLUP_VIEWER_GUI
#define GLUP_VIEWER_GUI

#include <geogram_gfx/glup_viewer/glup_viewer.h>

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * \file geogram_gfx/glup_viewer/glup_viewer_gui_private.h
     * \short These functions are used internally by
     *  glup_viewer for interfacing with ImGUI.
     */
    
    struct GLFWwindow;

    /**
     * \brief Initializes GLUP viewer GUI structures.
     * \details Called at application startup.
     * \param[in] w a pointer to the GLFWwindow
     */
    void GLUP_VIEWER_API glup_viewer_gui_init(GLFWwindow* w);

    /**
     * \brief Deletes all the variables allocated for the GUI.
     * \details Called at application shutdown.
     */
    void GLUP_VIEWER_API glup_viewer_gui_cleanup(void);

    /**
     * \brief Starts a new frame.
     */
    void GLUP_VIEWER_API glup_viewer_gui_begin_frame(void);

    /**
     * \brief Ends a frame.
     */
    void GLUP_VIEWER_API glup_viewer_gui_end_frame(void);

    /**
     * \brief Tests whether callbacks are directed towards
     *  GUI elements or application.
     * \retval non-zero if callbacks are directed towards GUI
     * \retval zero if callbacks should be taken into account
     *  by the application.
     */
    int GLUP_VIEWER_API glup_viewer_gui_takes_input(void);

    /**
     * \brief Callback for mouse button events.
     * \param[in] window a pointer to the GLFWwindow
     * \param[in] button the mouse button
     * \param[in] action the action
     * \param[in] mods modifiers
     */
    void GLUP_VIEWER_API glup_viewer_gui_mouse_button_callback(
        GLFWwindow* window, int button, int action, int mods
    );

    /**
     * \brief Callback for mouse wheel events.
     * \param[in] window a pointer to the GLFWwindow
     * \param[in] xoffset horizontal displacement
     * \param[in] yoffset vertical displacement
     */
    void GLUP_VIEWER_API glup_viewer_gui_scroll_callback(
        GLFWwindow* window, double xoffset, double yoffset
    );

    /**
     * \brief Callback for low-level key events
     * \param[in] window a pointer to the GLFWwindow
     * \param[in] key the key code
     * \param[in] scancode the scan code
     * \param[in] action the action
     * \param[in] mods modifiers
     */
    void GLUP_VIEWER_API glup_viewer_gui_key_callback(
        GLFWwindow* window, int key, int scancode, int action, int mods
    );

    /**
     * \brief Callback for high-level key events.
     * \param[in] window a pointer to the GLFWwindow
     * \param[in] c the character that corresponds to the
     *  pushed key
     */
    void GLUP_VIEWER_API glup_viewer_gui_char_callback(
	GLFWwindow* window, unsigned int c
    );

    /**
     * \brief Callback for window resize events.
     * \param[in] width the new width
     * \param[in] height the new height
     */
    void GLUP_VIEWER_API glup_viewer_gui_resize(int width, int height);

    /**
     * \brief Redraws the GUI and the scene.
     * \details It can be used from a command to update the graphics during
     *  a computation. It is called whenever a message is displayed in the
     *  console or whenever the progress bar is updated. It is ignored if called
     *  from the redraw or from the overlay callbacks.
     */
    void GLUP_VIEWER_API glup_viewer_gui_update(void);

    /**
     * \brief Tests a boolean command line argument.
     * \param[in] param the name of the argument
     * \return the value of the argument
     */
    GLboolean GLUP_VIEWER_API glup_viewer_get_arg_bool(const char* param);

    /**
     * \brief Tests a string command line argument.
     * \param[in] param the name of the argument
     * \param[in] arg the value to be tested
     * \retval GLtrue if command line argument \p param has value \p arg
     * \retval GLfalse otherwise
     */
    GLboolean GLUP_VIEWER_API glup_viewer_test_arg_string(
	const char* param, const char* arg
    );


    /**
     * \brief Gets the screen size from the arglist and sends it to GLUP.
     */
    void GLUP_VIEWER_API glup_viewer_set_screen_size_from_args(void);


    void GLUP_VIEWER_API glup_viewer_effect_begin_frame(void);
    
    void GLUP_VIEWER_API glup_viewer_effect_end_frame(void);

#ifdef __cplusplus
}
#endif

#endif
