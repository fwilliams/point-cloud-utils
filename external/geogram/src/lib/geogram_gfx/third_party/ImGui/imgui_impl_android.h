/*
 * ImGui Platform Binding for: Android
 * Author: Bruno Levy   Sun Aug 19 08:01:39 CEST 2018
 * Note: not part (yet) of the official ImGui distribution
 */ 

#ifdef __ANDROID__
#include <android_native_app_glue.h>

// param app: if non-null, registers input handler to specified app.
IMGUI_IMPL_API bool     ImGui_ImplAndroid_Init(struct android_app* app = nullptr);

IMGUI_IMPL_API void     ImGui_ImplAndroid_Shutdown();

// Needs to be called at the beginning of each frame,
// before ImGui::NewFrame().
IMGUI_IMPL_API void     ImGui_ImplAndroid_NewFrame();

// Needs to be called at the end of each frame,
// after all other ImGui functions.
IMGUI_IMPL_API void     ImGui_ImplAndroid_EndFrame();

// The event handler, to be used if not registered by ImGui_ImplAndroid_Init().
IMGUI_IMPL_API int32_t  ImGui_ImplAndroid_InputEvent(
    struct android_app* app, AInputEvent* event
);
#endif