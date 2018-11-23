/*
 *  Copyright (c) 2012-2016, Bruno Levy
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

#ifndef GEOGRAM_GFX_GLUP_GLUP_CONTEXT_VANILLAGL
#define GEOGRAM_GFX_GLUP_GLUP_CONTEXT_VANILLAGL

#include <geogram_gfx/basic/common.h>
#include <geogram_gfx/GLUP/GLUP.h>
#include <geogram_gfx/GLUP/GLUP_context.h>

/**
 * \file geogram_gfx/GLUP/GLUP_context_VanillaGL.h
 * \brief Internal implementation of GLUP using plain old OpenGL.
 */

#ifdef GEO_GL_LEGACY

namespace GLUP {
    using namespace GEO;

    /*********************************************************************/

    /**
     * \brief Implementation of GLUP using Vanilla (old-style) OpenGL.
     * \details This implementation does not use any shader. It is used
     *  as a fallback when the initialization of the other ones fails.
     *  Some primitive may be not implemented, degraded or of very low
     *  performance.
     */
    class Context_VanillaGL : public Context {
    public:
        /**
         * \brief Context_VanillaGL constructor.
         */
        Context_VanillaGL();

        /**
         * \brief Context_VanillaGL destructor.
         */
        virtual ~Context_VanillaGL();
        
        /**
         * \copydoc Context::profile_name()
         */
        virtual const char* profile_name() const;
        
        /**
         * \copydoc Context::primitive_supports_array_mode()
         */
        virtual bool primitive_supports_array_mode(GLUPprimitive prim) const;

        /**
         * \copydoc Context::begin()
         */
        virtual void begin(GLUPprimitive primitive);

        /**
         * \copydoc Context::end()
         */
        virtual void end();
        
    protected:

        /**
         * \copydoc Context::prepare_to_draw()
         */
        virtual void prepare_to_draw(GLUPprimitive primitive);
        
        /**
         * \copydoc Context::do_update_uniform_buffer()
         */
        virtual void do_update_uniform_buffer();
        
        /**
         * \copydoc Context::update_matrices()
         */
        virtual void update_matrices();

        /**
         * \copydoc Context::update_lighting()
         */
        virtual void update_lighting();
        
        /**
         * \brief Configures texturing-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin()
         */
        void configure_OpenGL_texturing();

        /**
         * \brief Configures lighting-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin()
         */
        void configure_OpenGL_lighting();

        /**
         * \brief Configures clipping-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin(). 
         */
        void configure_OpenGL_clipping();
        
        /**
         * \brief Configures lighting-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin(). It needs
         *  to be called after configure_OpenGL_texturing() and
         *  configure_OpenGL_lighting() since it overrides texturing and
         *  lighting settings.
         */
        void configure_OpenGL_picking();        

        
        /**
         * \copydoc Context::setup()
         */
        virtual void setup();
        
        /**
         * \copydoc Context::flush_immediate_buffers()
         */
        virtual void flush_immediate_buffers();


        /**
         * \brief Flushes the immediate buffer with the
         *  current drawing modes. 
         * \details This function is separated from 
         *  flush_immediate_buffers(), since we need to
         *  flush the buffer twice when mesh drawing is
         *  enabled.
         */
        virtual void flush_immediate_buffers_once();
        
        /**
         * \copydoc Context::setup_immediate_buffers()
         */
        virtual void setup_immediate_buffers();


        /**
         * \copydoc Context::setup_primitives()
         */
        virtual void setup_primitives();
        
        /**
         * \copydoc Context::get_state_variable_address()
         */
        Memory::pointer get_state_variable_address(const char* name);

        /**
         * \brief Tests whether cells should be sliced.
         * \retval true if cells should be sliced
         * \retval false otherwise
         */
        bool clip_slice_cells() const {
            return (
                uniform_state_.toggle[GLUP_CLIPPING].get() &&
                uniform_state_.clipping_mode.get() == GLUP_CLIP_SLICE_CELLS &&
                immediate_state_.primitive() >= GLUP_TETRAHEDRA
            ) ;
        }
        
        /**
         * \brief Sends a vertex and its optional attributes to OpenGL.
         * \param[in] v the index of the vertex from the immediate buffer.
         */
        void output_vertex(index_t v) {
            if(immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].is_enabled()) {
                glColor4fv(
                    immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].element_ptr(v)
                );
            }
            if(immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].is_enabled()) {
                glTexCoord4fv(
                    immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].
                    element_ptr(v)
                );
            }
	    if(immediate_state_.buffer[GLUP_NORMAL_ATTRIBUTE].is_enabled()) {
                glNormal3fv(
                    immediate_state_.buffer[GLUP_NORMAL_ATTRIBUTE].
                    element_ptr(v)
                );
	    }
            glVertex4fv(
                immediate_state_.buffer[GLUP_VERTEX_ATTRIBUTE].element_ptr(v)
            );
        }

        /**
         * \brief Sends a triangle normal to OpenGL
         * \param[in] v1 , v2 , v3 the indices of the three vertices from
         *  the immediate buffer.
         */
        void output_normal(index_t v1, index_t v2, index_t v3);

        /**
         * \brief Sends a quad normal to OpenGL
         * \param[in] v1 , v2 , v3 , v4 the indices of the four vertices from
         *  the immediate buffer.
         */
        void output_normal(index_t v1, index_t v2, index_t v3, index_t v4);


        /**
         * \brief Sends a picking id to OpenGL and encodes it as a color.
         * \details The current base picking id is added to the id.
         *  If picking is deactivated or constant by object, 
         *  it does nothing.
         */
        void output_picking_id(index_t id) {
            if(pick_primitives_) {
                glPickingIdAsColor(
                    index_t(uniform_state_.base_picking_id.get()) + id
                );
            }
        }

        /**
         * \brief Encodes a 32 bits integer into the 
         *  current OpenGL color.
         * \details This function is used by the 
         *  picking mechanism.
         * \param[in] id the picking id to be encoded as 
         *  the current OpenGL color
         */
        static void glPickingIdAsColor(index_t id) {
            GLubyte r = GLubyte( id        & 255);
            GLubyte g = GLubyte((id >> 8)  & 255);
            GLubyte b = GLubyte((id >> 16) & 255);
            GLubyte a = GLubyte((id >> 24) & 255);
            glColor4ub(r,g,b,a);
        }
        
        /**
         * \brief Sends a flat-shaded triangle to OpenGL
         * \param[in] v1 , v2 , v3 the indices of the three vertices from the
         *  immediate buffer.
         */
        void draw_triangle(index_t v1, index_t v2, index_t v3) {
            if(
		uniform_state_.toggle[GLUP_LIGHTING].get() &&
		!uniform_state_.toggle[GLUP_VERTEX_NORMALS].get()
	    ) {
                output_normal(v1,v2,v3);
            }
            output_vertex(v1);
            output_vertex(v2);
            output_vertex(v3);
        }

        /**
         * \brief Sends a flat-shaded quad to OpenGL
         * \param[in] v1 , v2 , v3 , v4 the indices of the three 
         *  vertices from the immediate buffer.
         */
        void draw_quad(index_t v1, index_t v2, index_t v3, index_t v4) {
            if(uniform_state_.toggle[GLUP_LIGHTING].get() &&
	       !uniform_state_.toggle[GLUP_VERTEX_NORMALS].get()	       
	    ) {
                output_normal(v1,v2,v3,v4);
            }
            output_vertex(v1);
            output_vertex(v2);
            output_vertex(v4);
            output_vertex(v3);            
        }

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as point primitives.
         */
        void draw_immediate_buffer_GLUP_POINTS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as line primitives.
         */
        void draw_immediate_buffer_GLUP_LINES();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as triangle primitives.
         */
        void draw_immediate_buffer_GLUP_TRIANGLES();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as quad primitives.
         */
        void draw_immediate_buffer_GLUP_QUADS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as tetrahedra primitives.
         */
        void draw_immediate_buffer_GLUP_TETRAHEDRA();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as hexahedra primitives.
         */
        void draw_immediate_buffer_GLUP_HEXAHEDRA();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as prism primitives.
         */
        void draw_immediate_buffer_GLUP_PRISMS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as pyramid primitives.
         */
        void draw_immediate_buffer_GLUP_PYRAMIDS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as connectors primitives.
         */
        void draw_immediate_buffer_GLUP_CONNECTORS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as spheres primitives.
         */
        void draw_immediate_buffer_GLUP_SPHERES();
	
        /**
         * \brief Draws all the primitives from the immediate buffer using
         *  the marching cells algorithm.
         * \details This function is used when clipping is enabled and when
         *  clippping mode is GLUP_CLIP_SLICE_CELLS
         */
        void draw_immediate_buffer_with_marching_cells(
            const MarchingCell& cell
        );

        void begin_indirect_texturing();

        void end_indirect_texturing();
        
    private:
        /**
         * \brief Indicates whether a picking id should be send to 
         *  OpenGL for each primitive.
         */
        bool pick_primitives_;

        /**
         * \brief The program to be used for indirect
         *  texturing.
         */
        GLuint indirect_texturing_program_;

	/**
	 * \brief A display list that draws a tesselated
	 *  unit sphere.
	 */
	GLuint sphere_display_list_;
    };

    /*********************************************************************/

    
}

#endif

#endif
