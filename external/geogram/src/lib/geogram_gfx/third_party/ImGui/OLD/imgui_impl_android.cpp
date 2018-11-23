// ImGui Platform Binding for: Android
// [Bruno Levy] Sun Aug 19 08:01:39 CEST 2018
// Note: not part (yet) of the official ImGui distribution
// 
// Note: to use, include in the CMakeLists.txt that compiles this file:
// if(ANDROID)
//   target_include_directories(geogram_gfx_third_party PRIVATE
//      ${ANDROID_NDK}/sources/android/native_app_glue
//   )
// endif()

// What works:
//   Rendering with OpenGL ES 2.x
//   Fingers/Stylus/Mouse interaction
//   Virtual and physical keyboard interaction

// TODO (Bugs to be fixed):
// ------------------------
//  - soft keyboard directional keys do not always work (it depends on the used keyboard,
//     for some keyboards, they work a little bit, randomly, for some others they work...)
//
//  - app is restarted when connecting/disconnecting a physical keyboard
//     while application is running (I do not understand, I have:
//        android:configChanges="orientation|keyboardHidden|keyboard"
//     in AndroidManifest.xml) -> this one is more related to android_main.cpp

// TODO (Improvements):
// --------------------
//  - UTF8 text input (probably not very difficult to add).
//
//  - mouse cursors (https://developer.android.com/about/versions/nougat/android-7.0#custom_pointer_api)
//    (need to overload Java function, cannot do that with native_glue I think, unless we can change
//     methods of an existing Java object with JNI)
//
//  - setMousePos



#ifdef __ANDROID__

#include "imgui.h"
#include "imgui_impl_android.h"

#include <EGL/egl.h>
#include <GLES/gl.h>
#include <time.h>
#include <cassert>
#include <stdexcept>
#include <cctype>

namespace {
    double g_Time = 0.0;
    float g_mouseX = 0.0f;
    float g_mouseY = 0.0f;
    bool  g_mousePressed[5] = {false, false, false, false, false};
    bool  g_resetKeys = false;
}

// Some utilities functions that interact with Android.
namespace AndroidUtils {
    
    // Shows or hides the software keyboard.
    void set_soft_keyboard_visibility(struct android_app* app, bool show);

    // Converts a keycode to a unicode.
    // deviceId, keyCode, metaState can be obtained from the InputEvent.
    jint keycode_to_unicode(
	struct android_app* app, int32_t deviceId, int32_t keyCode, int32_t metaState
    );
}


bool ImGui_ImplAndroid_Init(struct android_app* app) {
    g_Time = 0.0;
    g_mouseX = 0.0f;
    g_mouseY = 0.0f;
    for (int i = 0; i < IM_ARRAYSIZE(g_mousePressed); i++) {
	g_mousePressed[i] = false;
    }
    // TODO: mouse cursor
    // TODO: setmousepos ?

    ImGuiIO& io = ImGui::GetIO();
    io.KeyMap[ImGuiKey_Tab] = AKEYCODE_TAB;
    io.KeyMap[ImGuiKey_LeftArrow] = AKEYCODE_DPAD_LEFT;
    io.KeyMap[ImGuiKey_RightArrow] = AKEYCODE_DPAD_RIGHT;
    io.KeyMap[ImGuiKey_UpArrow] = AKEYCODE_DPAD_UP;
    io.KeyMap[ImGuiKey_DownArrow] = AKEYCODE_DPAD_DOWN;
    io.KeyMap[ImGuiKey_PageUp] = AKEYCODE_PAGE_UP;
    io.KeyMap[ImGuiKey_PageDown] = AKEYCODE_PAGE_DOWN;
    io.KeyMap[ImGuiKey_Home] = AKEYCODE_MOVE_HOME;
    io.KeyMap[ImGuiKey_End] = AKEYCODE_MOVE_END;
    io.KeyMap[ImGuiKey_Insert] = AKEYCODE_INSERT;
    io.KeyMap[ImGuiKey_Delete] = AKEYCODE_FORWARD_DEL;
    io.KeyMap[ImGuiKey_Backspace] = AKEYCODE_DEL;
    io.KeyMap[ImGuiKey_Space] = AKEYCODE_SPACE;
    io.KeyMap[ImGuiKey_Enter] = AKEYCODE_ENTER;
    io.KeyMap[ImGuiKey_Escape] = AKEYCODE_ESCAPE;
    io.KeyMap[ImGuiKey_A] = AKEYCODE_A;
    io.KeyMap[ImGuiKey_C] = AKEYCODE_C;
    io.KeyMap[ImGuiKey_V] = AKEYCODE_V;
    io.KeyMap[ImGuiKey_X] = AKEYCODE_X;
    io.KeyMap[ImGuiKey_Y] = AKEYCODE_Y;
    io.KeyMap[ImGuiKey_Z] = AKEYCODE_Z;
    

    // Install callbacks
    if(app != nullptr) {
	app->onInputEvent = ImGui_ImplAndroid_InputEvent;
    }
    
    return true;
}

void ImGui_ImplAndroid_Shutdown() {
    // TODO: destroy mouse cursors.
}

static void ImGui_ImplAndroid_UpdateMousePosAndButtons()
{
    ImGuiIO& io = ImGui::GetIO();

    for (int i = 0; i < IM_ARRAYSIZE(io.MouseDown); i++) {
        io.MouseDown[i] = g_mousePressed[i];
    }
    io.MousePos = ImVec2(g_mouseX, g_mouseY);
}

static void ImGui_ImplAndroid_UpdateMouseCursor() {
    // TODO...
}

void ImGui_ImplAndroid_NewFrame() {
    ImGuiIO& io = ImGui::GetIO();
    // Font atlas needs to be built, call renderer _NewFrame() function
    // e.g. ImGui_ImplOpenGL3_NewFrame()     
    IM_ASSERT(io.Fonts->IsBuilt());     

    // Get current display size
    EGLDisplay display = eglGetDisplay(EGL_DEFAULT_DISPLAY);
    EGLSurface surface = eglGetCurrentSurface(EGL_DRAW);
    int w=0;
    int h=0;
    eglQuerySurface(display, surface, EGL_WIDTH, &w);
    eglQuerySurface(display, surface, EGL_HEIGHT, &h);
    int display_w = w;
    int display_h = h;

    io.DisplaySize = ImVec2((float)w, (float)h);
    io.DisplayFramebufferScale =
	ImVec2(
	    w > 0 ?((float)display_w / w) : 0,
	    h > 0 ? ((float)display_h / h) : 0
    );

    // Setup time step
    timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    double current_time = double(now.tv_sec) + double(now.tv_nsec) * 1e-9;
    
    io.DeltaTime = g_Time > 0.0 ? float(current_time - g_Time) : 1.0f/60.0f;
    g_Time = current_time;

    ImGui_ImplAndroid_UpdateMousePosAndButtons();
    ImGui_ImplAndroid_UpdateMouseCursor();

    // TODO: Gamepad navigation mapping ?
}


void ImGui_ImplAndroid_EndFrame() {
    // g_resetKeys is set when the latest key event came from the soft keyboard,
    // then we need to reset the keys.
    if(g_resetKeys) {
	ImGuiIO& io = ImGui::GetIO();
	for(int key = 0; key < IM_ARRAYSIZE(io.KeysDown); ++key) {
	    io.KeysDown[key] = false;	    
	}
	io.KeyShift = false;
	io.KeyCtrl = false;
	io.KeyAlt = false;
	io.KeySuper = false;
	g_resetKeys = false;
    }
}

// Emulates mouse buttons using multiple fingers:
//   emulated mouse button is determined by number of fingers
//   coordinates are defined by last finger
int32_t  ImGui_ImplAndroid_FingerEvent(
    struct android_app* app, AInputEvent* event
) {
    int32_t action = AMotionEvent_getAction(event);
    bool down_or_move = (action == AMOTION_EVENT_ACTION_DOWN ||
			 action == AMOTION_EVENT_ACTION_MOVE );

    int nb_fingers = int(AMotionEvent_getPointerCount(event));
    int btn = nb_fingers-1;
    for(int i=0; i<IM_ARRAYSIZE(g_mousePressed); ++i) {
	if(i == btn) {
	    g_mousePressed[i] = down_or_move;
	} else {
	    g_mousePressed[i] = false;
	}
    }
    g_mouseX = AMotionEvent_getX(event, nb_fingers-1);
    g_mouseY = AMotionEvent_getY(event, nb_fingers-1);
    
    return 1;
}

// Handles stylus input, like Galaxy SPen. Uses the tiny button on
// the pen to emulate second mouse button.
int32_t ImGui_ImplAndroid_StylusEvent(
    struct android_app* app, AInputEvent* event
) {
    int32_t action = AMotionEvent_getAction(event);
    bool down_or_move = (action == AMOTION_EVENT_ACTION_DOWN ||
			 action == AMOTION_EVENT_ACTION_MOVE );
    int btn = (
	(AMotionEvent_getButtonState(event) &
	    AMOTION_EVENT_BUTTON_STYLUS_PRIMARY) != 0
    ) ? 1 : 0;
    for(int i=0; i<IM_ARRAYSIZE(g_mousePressed); ++i) {
	if(i == btn) {
	    g_mousePressed[i] = down_or_move;
	} else {
	    g_mousePressed[i] = false;
	}
    }
    g_mouseX = AMotionEvent_getX(event, 0);
    g_mouseY = AMotionEvent_getY(event, 0);
    return 1;    
}

// Handles a standard USB or bluetooth mouse connected to the phone.
int32_t  ImGui_ImplAndroid_MouseEvent(
    struct android_app* app, AInputEvent* event
) {
    int32_t buttons = AMotionEvent_getButtonState(event);
    for(int i=0; i<IM_ARRAYSIZE(g_mousePressed); ++i) {
	g_mousePressed[i] = false;
    }
    g_mousePressed[0] = (buttons &  AMOTION_EVENT_BUTTON_PRIMARY) != 0;
    g_mousePressed[1] = (buttons &  AMOTION_EVENT_BUTTON_SECONDARY) != 0;
    g_mousePressed[2] = (buttons &  AMOTION_EVENT_BUTTON_TERTIARY) != 0;
    g_mouseX = AMotionEvent_getX(event, 0);
    g_mouseY = AMotionEvent_getY(event, 0);

    // Mouse wheel
    int32_t action = AMotionEvent_getAction(event);
    if(action == AMOTION_EVENT_ACTION_SCROLL) {
	float hscroll = AMotionEvent_getAxisValue(
	    event, AMOTION_EVENT_AXIS_HSCROLL, 0
        );
	float vscroll = AMotionEvent_getAxisValue(
	    event, AMOTION_EVENT_AXIS_VSCROLL, 0
        );
	ImGuiIO& io = ImGui::GetIO();
	io.MouseWheelH += hscroll;
	io.MouseWheel  += vscroll;
    }
    
    return 1;    
}

int32_t ImGui_ImplAndroid_MotionEvent(
    struct android_app* app, AInputEvent* event
) {
    int32_t result = 0;
    switch(AMotionEvent_getToolType(event,0)) {
	case AMOTION_EVENT_TOOL_TYPE_FINGER:
	    result = ImGui_ImplAndroid_FingerEvent(app, event);
	    break;
	case AMOTION_EVENT_TOOL_TYPE_STYLUS:
	    result = ImGui_ImplAndroid_StylusEvent(app, event);	    
	    break;
	case AMOTION_EVENT_TOOL_TYPE_MOUSE:
	    result = ImGui_ImplAndroid_MouseEvent(app, event);	    	    
	    break;
	default:
	    break;
    }
    return result;
}

int32_t ImGui_ImplAndroid_KeyEvent(
    struct android_app* app, AInputEvent* event
) {
    // Note: important to return 1 on BACK key pressed,
    // else this triggers default behavior that stops the
    // application brutally.

    ImGuiIO& io = ImGui::GetIO();
    
    int32_t action = AKeyEvent_getAction(event);
    int32_t key = AKeyEvent_getKeyCode(event);
    int32_t modifiers = AKeyEvent_getMetaState(event);
    int32_t device = AInputEvent_getDeviceId(event);    

    if(key >= 0 && key < IM_ARRAYSIZE(io.KeysDown)) {
	if((AKeyEvent_getFlags(event) & AKEY_EVENT_FLAG_SOFT_KEYBOARD)) {
	    // The soft keyboard generates Push/Release events when the
	    // key is released. Thus we mark the key as pushed, and
	    // set g_resetKeys so that ImGui_ImplAndroid_EndFrame()
	    // will mark the key as released after ImGui could do what
	    // it has to do with the key.
	    io.KeysDown[key] = true;
	    g_resetKeys = true;
	} else {
	    io.KeysDown[key] = (action == AKEY_EVENT_ACTION_DOWN);
	    g_resetKeys = false;	    
	}
	io.KeyShift = ((modifiers & AMETA_SHIFT_ON) != 0);
	io.KeyCtrl = ((modifiers & AMETA_CTRL_ON) != 0);
	io.KeyAlt = ((modifiers & AMETA_ALT_ON) != 0);
	io.KeySuper = ((modifiers & AMETA_META_ON) != 0);
    }

    if(action == AKEY_EVENT_ACTION_DOWN) {
	if(key == AKEYCODE_BACK) {
	    ImGui::debug_printf("Show virtual keyboard\n");
	    AndroidUtils::set_soft_keyboard_visibility(app, true);
	} else {
	    jint unicode = AndroidUtils::keycode_to_unicode(
		app, device, key, modifiers
	    );
	    // TODO: use AddInputCharactersUTF8()
	    char c = char(unicode);
	    if(isprint(c)) {
		io.AddInputCharacter(c);
	    }
	}
    }
    
    return 1;
}

int32_t  ImGui_ImplAndroid_InputEvent(
    struct android_app* app, AInputEvent* event
) {
    int32_t result = 0;
    switch(AInputEvent_getType(event)) {
	case AINPUT_EVENT_TYPE_MOTION:
	    result = ImGui_ImplAndroid_MotionEvent(app, event);
	    break;
	case AINPUT_EVENT_TYPE_KEY:
	    result = ImGui_ImplAndroid_KeyEvent(app, event);
	    break;
	default:
	    break;
    }
    return result;
}
    
/*****************************************************************/

// Functions that interact with Java.

namespace AndroidUtils {

   // Display soft keyboard programmatically
   //https://groups.google.com/forum/?fromgroups=#!topic/android-ndk/Tk3g00wLKhk
   //   (see alto messages about how to attach/detach thread).
   // There is a function supposed to do that: 
   // ANativeActivity_showSoftInput(
   //   mApplication->activity,ANATIVEACTIVITY_SHOW_SOFT_INPUT_FORCED
   // );
   // (or ANATIVEACTIVITY_SHOW_SOFT_INPUT_IMPLICIT) 
   // but I did not manage to make it work.
    
    void set_soft_keyboard_visibility(struct android_app* app, bool pShow) {

	JavaVM* lJavaVM = app->activity->vm;
	JNIEnv* lJNIEnv = nullptr; 
	bool lThreadAttached = false;

	// Get JNIEnv from lJavaVM using GetEnv to test whether
	// thread is attached or not to the VM. If not, attach it
	// (and note that it will need to be detached at the end
	//  of the function).
	switch (lJavaVM->GetEnv((void**)&lJNIEnv, JNI_VERSION_1_6)) {
	    case JNI_OK:
		break;
	    case JNI_EDETACHED: {
		jint lResult = lJavaVM->AttachCurrentThread(&lJNIEnv, nullptr);
		if(lResult == JNI_ERR) {
		    throw std::runtime_error("Could not attach current thread");
		}
		lThreadAttached = true;
	    } break;
	    case JNI_EVERSION:
		throw std::runtime_error("Invalid java version");
	}
    
	// Retrieves NativeActivity.
	jobject lNativeActivity = app->activity->clazz;
	jclass ClassNativeActivity = lJNIEnv->GetObjectClass(lNativeActivity);

	// Retrieves Context.INPUT_METHOD_SERVICE.
	jclass ClassContext = lJNIEnv->FindClass("android/content/Context");
	jfieldID FieldINPUT_METHOD_SERVICE =
	    lJNIEnv->GetStaticFieldID(
		ClassContext,"INPUT_METHOD_SERVICE", "Ljava/lang/String;"
	);
	jobject INPUT_METHOD_SERVICE =
	    lJNIEnv->GetStaticObjectField(
		ClassContext, FieldINPUT_METHOD_SERVICE
	);

	// Runs getSystemService(Context.INPUT_METHOD_SERVICE).
	jclass ClassInputMethodManager = lJNIEnv->FindClass(
	    "android/view/inputmethod/InputMethodManager"
	);
	jmethodID MethodGetSystemService = lJNIEnv->GetMethodID(
	    ClassNativeActivity, "getSystemService",
	    "(Ljava/lang/String;)Ljava/lang/Object;"
	);
	jobject lInputMethodManager = lJNIEnv->CallObjectMethod(
	    lNativeActivity, MethodGetSystemService,
	    INPUT_METHOD_SERVICE
	);

	// Runs getWindow().getDecorView().
	jmethodID MethodGetWindow = lJNIEnv->GetMethodID(
	    ClassNativeActivity, "getWindow",
	    "()Landroid/view/Window;"
	);
	jobject lWindow = lJNIEnv->CallObjectMethod(
	    lNativeActivity, MethodGetWindow
	);
	jclass ClassWindow = lJNIEnv->FindClass("android/view/Window");
	jmethodID MethodGetDecorView = lJNIEnv->GetMethodID(
	    ClassWindow, "getDecorView", "()Landroid/view/View;"
	);
	jobject lDecorView = lJNIEnv->CallObjectMethod(
	    lWindow, MethodGetDecorView
	);

	if (pShow) {
	    // Runs lInputMethodManager.showSoftInput(...).
	    jmethodID MethodShowSoftInput = lJNIEnv->GetMethodID(
		ClassInputMethodManager, "showSoftInput",
		"(Landroid/view/View;I)Z"
	    );
	    lJNIEnv->CallBooleanMethod(
		lInputMethodManager, MethodShowSoftInput,
		lDecorView, 0
	    );
	} else {
	    // Runs lWindow.getViewToken()
	    jclass ClassView = lJNIEnv->FindClass(
		"android/view/View"
	    );
	    jmethodID MethodGetWindowToken = lJNIEnv->GetMethodID(
		ClassView, "getWindowToken", "()Landroid/os/IBinder;"
	    );
	    jobject lBinder = lJNIEnv->CallObjectMethod(
		lDecorView, MethodGetWindowToken
	    );

	    // lInputMethodManager.hideSoftInput(...).
	    jmethodID MethodHideSoftInput = lJNIEnv->GetMethodID(
		ClassInputMethodManager, "hideSoftInputFromWindow",
		"(Landroid/os/IBinder;I)Z"
	    );
	    lJNIEnv->CallBooleanMethod(
		lInputMethodManager, MethodHideSoftInput,
		lBinder, 0
	    );
	}

	if(lThreadAttached) {
	    lJavaVM->DetachCurrentThread();
	}
    }

    jint keycode_to_unicode(
	struct android_app* app,
	int32_t pDeviceId, int32_t pKeyCode, int32_t pMetaState
    ) {
	jint result = 0;

	// Early exit for special keys
	// (works without it, but well, why calling all that
	//  Java stuff if we now in advance that we do not need
	//  to ?).
	if(
	    pKeyCode == AKEYCODE_TAB ||
	    pKeyCode == AKEYCODE_DPAD_LEFT ||
	    pKeyCode == AKEYCODE_DPAD_RIGHT ||
	    pKeyCode == AKEYCODE_DPAD_UP ||
	    pKeyCode == AKEYCODE_DPAD_DOWN ||
	    pKeyCode == AKEYCODE_PAGE_UP ||
	    pKeyCode == AKEYCODE_PAGE_DOWN ||
	    pKeyCode == AKEYCODE_MOVE_HOME ||
	    pKeyCode == AKEYCODE_MOVE_END ||
	    pKeyCode == AKEYCODE_INSERT ||
	    pKeyCode == AKEYCODE_FORWARD_DEL ||
	    pKeyCode == AKEYCODE_DEL ||
	    pKeyCode == AKEYCODE_ENTER ||
	    pKeyCode == AKEYCODE_ESCAPE
	) {
	    return result;
	}

	
	JavaVM* lJavaVM = app->activity->vm;
	JNIEnv* lJNIEnv = nullptr; 
	bool lThreadAttached = false;

	// Get JNIEnv from lJavaVM using GetEnv to test whether
	// thread is attached or not to the VM. If not, attach it
	// (and note that it will need to be detached at the end
	//  of the function).
	switch (lJavaVM->GetEnv((void**)&lJNIEnv, JNI_VERSION_1_6)) {
	    case JNI_OK:
		break;
	    case JNI_EDETACHED: {
		jint lResult = lJavaVM->AttachCurrentThread(&lJNIEnv, nullptr);
		if(lResult == JNI_ERR) {
		    throw std::runtime_error("Could not attach current thread");
		}
		lThreadAttached = true;
	    } break;
	    case JNI_EVERSION:
		throw std::runtime_error("Invalid java version");
	}

	jclass ClassKeyCharacterMap = lJNIEnv->FindClass(
	    "android/view/KeyCharacterMap"
	);

	jmethodID MethodLoad = lJNIEnv->GetStaticMethodID(
	    ClassKeyCharacterMap, "load",
	    "(I)Landroid/view/KeyCharacterMap;"
	);

	jobject lKeyCharacterMap = lJNIEnv->CallStaticObjectMethod(
	    ClassKeyCharacterMap, MethodLoad, jint(pDeviceId)
	);

	jmethodID MethodGet = lJNIEnv->GetMethodID(
	    ClassKeyCharacterMap, "get",
	    "(II)I"
	);

	result = lJNIEnv->CallIntMethod(
	    lKeyCharacterMap, MethodGet,
	    jint(pKeyCode), jint(pMetaState)
	);

	if(lThreadAttached) {
	    lJavaVM->DetachCurrentThread();
	}

	return result;
    }
}

#endif

/********************************************************************/

/*

 Notes, links etc...
 ===================

 In pre-4.3 Androids, there was a bug on some devices making the 
 app. crash when hiding the soft keyboard. Normally it was fixed in Android 4.3.

 https://stackoverflow.com/questions/15913080/crash-when-closing-soft-keyboard-while-using-native-activity

*/

