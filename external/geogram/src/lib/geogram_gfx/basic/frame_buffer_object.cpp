/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <geogram_gfx/basic/frame_buffer_object.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram/basic/logger.h>
#include <string>

#ifdef GEO_OS_EMSCRIPTEN
#include <emscripten.h>
#endif

namespace GEO {

    FrameBufferObject::FrameBufferObject() : 
        frame_buffer_id(0),
        depth_buffer_id(0),
        offscreen_id(0),
        width(0),
        height(0),
        previous_frame_buffer_id(0)
    {
    }

    FrameBufferObject::~FrameBufferObject() {
        if(offscreen_id != 0) {
            glDeleteTextures(1, &offscreen_id);
            offscreen_id = 0;
        }
        if(depth_buffer_id != 0) {
            glDeleteTextures(1, &depth_buffer_id);
            depth_buffer_id = 0;
        }
        if(frame_buffer_id != 0) {
            glBindFramebuffer(GL_FRAMEBUFFER, previous_frame_buffer_id);
            glDeleteFramebuffers(1, &frame_buffer_id);
        }
    }

    void FrameBufferObject::resize(index_t new_width, index_t new_height) {
        if(width != new_width || height != new_height) {
            width  = new_width; 
            height = new_height;

	    GEO_CHECK_GL();
            glBindTexture(GL_TEXTURE_2D, offscreen_id);	    
            glTexImage2D(
                GL_TEXTURE_2D, 0, internal_storage,
                GLsizei(width), GLsizei(height), 0,
                GL_RGBA, GL_FLOAT, nullptr
            );
	    GEO_CHECK_GL();	    
            if(depth_buffer_id != 0) {
                glBindTexture(GL_TEXTURE_2D, depth_buffer_id);
		GEO_CHECK_GL();		
                glTexImage2D(
                    GL_TEXTURE_2D, 0,
		    GL_DEPTH_COMPONENT,
                    GLsizei(width), GLsizei(height), 0,
                    GL_DEPTH_COMPONENT,
		    GL_UNSIGNED_SHORT,
		    nullptr
                );
		GEO_CHECK_GL();				
            }
            glBindTexture(GL_TEXTURE_2D, 0);
        }
    }

    bool FrameBufferObject::initialize(
        index_t width_in, index_t height_in, 
        bool with_depth_buffer, GLint internal_storage_in,
        bool mipmaps
    ) {

#ifdef GEO_OS_EMSCRIPTEN
	// Initialize WEBGL extension for floating point textures.
	// Note: the object returned by getExtension() needs to be
	// stored (else when it is garbage collected the extension
	// is deinitialized). It is stored in a new field of the
	// Emscripten Module object.
	if(with_depth_buffer) {
	    EM_ASM(
		if(!Module.WEBGL_depth_texture_ext) {
		    Module.WEBGL_depth_texture_ext = Module.ctx.getExtension("WEBGL_depth_texture");
		    if(!Module.WEBGL_depth_texture_ext) {
			alert("No WEBGL_depth_texture support");
		    }
		}
		// In older browser: maybe try as well:
		//   WEBKIT_WEBGL_depth_texture
		//   MOZ_WEBGL_depth_texture"
	    );
	}
#endif	
	
	GEO_CHECK_GL();
	
        //  Get the id of the default frame buffer used by the
        // context. It will be 0 if it is a regular OpenGL context,
        // or it can be a bound frame buffer object since Qt5.4
        // QOpenGLWidget uses a frame buffer to implement light weight
        // OpenGL contexts. Note that it may change when the Qt rendering
	// context is resized (see bind_as_framebuffer()).
	glGetIntegerv(
	    GL_FRAMEBUFFER_BINDING, (GLint*)(&previous_frame_buffer_id)
	);
	
        width = width_in;
        height = height_in;
        internal_storage = internal_storage_in;

        // Generate frame buffer object then bind it.
        glGenFramebuffers(1, &frame_buffer_id);
        glBindFramebuffer(GL_FRAMEBUFFER, frame_buffer_id);

	GEO_CHECK_GL();
	
        if(with_depth_buffer) {
            glGenTextures(1, &depth_buffer_id);
        }

	GEO_CHECK_GL();
	
        // Create the texture we will be using to render to.
        glGenTextures(1, &offscreen_id);
        glBindTexture(GL_TEXTURE_2D, offscreen_id);

	GEO_CHECK_GL();
	
        if(!mipmaps) {
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        }

	GEO_CHECK_GL();
	
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	GEO_CHECK_GL();
	
        glTexImage2D(
            GL_TEXTURE_2D, 0, internal_storage,
            GLsizei(width), GLsizei(height), 0,
            GL_RGBA,	    
	    GL_FLOAT,
	    nullptr
        );

	GEO_CHECK_GL();
	
        // Bind the texture to the frame buffer.
        glFramebufferTexture2D(
            GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
            GL_TEXTURE_2D, offscreen_id, 0
        );

	GEO_CHECK_GL();
	
        // Initialize the depth buffer.
        if(with_depth_buffer) {

            glBindTexture(GL_TEXTURE_2D, depth_buffer_id);
            glTexImage2D(
                GL_TEXTURE_2D, 0,
		GL_DEPTH_COMPONENT,
                GLsizei(width), GLsizei(height), 0,
                GL_DEPTH_COMPONENT,
		GL_UNSIGNED_SHORT,
		nullptr
            );

	    GEO_CHECK_GL();
	
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	    GEO_CHECK_GL();
	
            glFramebufferTexture2D(
                GL_FRAMEBUFFER,
		GL_DEPTH_ATTACHMENT,		
                GL_TEXTURE_2D, depth_buffer_id, 0
            );

	    GEO_CHECK_GL();	    
        }

	GEO_CHECK_GL();

//      Heard it said that DrawBuffer and ReadBuffer
//      should be set to NONE before checking frame buffer
//      status, but it seems to be unnecessary (commented-out
//      for now)
//	
//	glDrawBuffer(GL_NONE);
//	glReadBuffer(GL_NONE);
	
        // Make sure we have not errors.
        if(
            glCheckFramebufferStatus(GL_FRAMEBUFFER) !=
            GL_FRAMEBUFFER_COMPLETE
        ) {
//	    glDrawBuffer(GL_BACK);
//	    glReadBuffer(GL_BACK);
            return false;
        }

	GEO_CHECK_GL();
	
//	glDrawBuffer(GL_BACK);
//	glReadBuffer(GL_BACK);
	
	GEO_CHECK_GL();
	
        // Restore previously bound frame buffer object.
        glBindFramebuffer(GL_FRAMEBUFFER, previous_frame_buffer_id);

	GEO_CHECK_GL();
	
        return true;
    }

    void FrameBufferObject::bind_as_texture() {
        glBindTexture(GL_TEXTURE_2D, offscreen_id);
    }

    void FrameBufferObject::bind_depth_buffer_as_texture() {
        glBindTexture(GL_TEXTURE_2D, depth_buffer_id);
    }

    void FrameBufferObject::bind_as_framebuffer() {
	// Current frame buffer ID may have changed,
	// for instance under Qt if rendering area was resized.
	glGetIntegerv(
	    GL_FRAMEBUFFER_BINDING, (GLint*)(&previous_frame_buffer_id)
	);
        glBindFramebuffer(GL_FRAMEBUFFER, frame_buffer_id);
    }

    void FrameBufferObject::unbind() {
        glBindFramebuffer(GL_FRAMEBUFFER, previous_frame_buffer_id);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    bool FrameBufferObject::is_bound_as_framebuffer() const {
	GLuint current_frame_buffer_id;
	glGetIntegerv(
	    GL_FRAMEBUFFER_BINDING, (GLint*)(&current_frame_buffer_id)
	);
	return (current_frame_buffer_id == frame_buffer_id);
    }
}

