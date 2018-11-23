/**
 * Macros and functions for facilitating interfacing
 * ImGui with WebGL when used with emscripten.
 * This file is not part of the initial ImGui distribution.
 * [Bruno Levy 11/26/2017]
 */

#ifndef H_GLUP_COMPAT_H
#define H_GLUP_COMPAT_H

#ifndef GL_POLYGON_MODE
#  define GL_POLYGON_MODE 0x0B40
#endif

#ifndef GL_FILL
#  define GL_FILL 0x1B02
#endif

#ifndef GL_VERTEX_ARRAY_BINDING
#  define GL_VERTEX_ARRAY_BINDING 0x85B5
#endif

#ifndef GL_SAMPLER_BINDING
#  define GL_SAMPLER_BINDING 0x8919
#endif

#ifndef GL_UNPACK_ROW_LENGTH
#  define GL_UNPACK_ROW_LENGTH 0x0CF2
#endif

namespace {

    inline void glup_glGetIntegerv(GLenum name, GLint* value) {
	if(name == GL_VERTEX_ARRAY_BINDING) {
	    *value = glupGetVertexArrayBinding();
	}
#if defined(__EMSCRIPTEN__)  || defined(__APPLE__)
	else if(name == GL_SAMPLER_BINDING) {
	    *value = 0;
	}
#endif

#if defined(__EMSCRIPTEN__) || defined(__ANDROID__)
	else if(name == GL_POLYGON_MODE) {
	    *value = 0;
	}
#endif
        else {
	    glGetIntegerv(name, value);
	}
    }

#if defined(__EMSCRIPTEN__) || defined(__ANDROID__)
#  ifdef glPolygonMode
#    undef glPolygonMode
#  endif    
   inline void glPolygonMode(GLenum, GLenum) {
   }
#endif    

#if defined(__EMSCRIPTEN__)  || defined(__APPLE__)
#  ifdef glBindSampler
#     undef glBindSampler
#  endif    
    inline void glBindSampler(GLenum, GLuint) {
    }
#endif
}


#ifdef glGetIntegerv
#  undef glGetIntegerv
#endif
#define glGetIntegerv glup_glGetIntegerv

#ifdef glBindVertexArray
#  undef glBindVertexArray
#endif
#define glBindVertexArray glupBindVertexArray

#ifdef glGenVertexArrays
#  undef glGenVertexArrays
#endif
#define glGenVertexArrays glupGenVertexArrays

#ifdef __EMSCRIPTEN__
inline void glup_glPixelStorei(GLenum name, GLint value) {
    if(name != GL_UNPACK_ROW_LENGTH) {
	glPixelStorei(name, value);
    }
}
#  define glPixelStorei glup_glPixelStorei
#endif

#endif // Include guard
