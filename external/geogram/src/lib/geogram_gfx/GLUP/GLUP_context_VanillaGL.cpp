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

#include <geogram_gfx/GLUP/GLUP_context_VanillaGL.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/basic/GL.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/logger.h>

#ifdef GEO_GL_LEGACY

namespace GLUP {
    using namespace GEO;

    /***********************************************************************/
    
    GLuint create_sphere_display_list();
	
    Context_VanillaGL::Context_VanillaGL() {
        use_core_profile_ = false;
        use_ES_profile_ = false;
        indirect_texturing_program_ = 0;
	sphere_display_list_ = create_sphere_display_list();
    }

    Context_VanillaGL::~Context_VanillaGL() {
        if(indirect_texturing_program_ != 0) {
            glDeleteProgramsARB(1, &indirect_texturing_program_);
            indirect_texturing_program_ = 0;
        }
	if(sphere_display_list_ != 0) {
	    glDeleteLists(sphere_display_list_, 1);
	    sphere_display_list_ = 0;
	}
    }

    const char* Context_VanillaGL::profile_name() const {
        return "VanillaGL";
    }
    
    void Context_VanillaGL::begin(GLUPprimitive primitive) {

	if(!primitive_info_[primitive].implemented) {
	    Logger::warn("GLUP")
		<< "glupBegin(): "
		<< glup_primitive_name(primitive)
		<< " not implemented in this profile" << std::endl;
	}
	
        update_uniform_buffer();

        if(uniform_state_.toggle[GLUP_VERTEX_COLORS].get()) {
            immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].enable();
        } else {
            immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].disable();
        }
        
        if(uniform_state_.toggle[GLUP_TEXTURING].get()) {
            immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].enable();
        } else {
            immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].disable();
        }

        if(
	    uniform_state_.toggle[GLUP_LIGHTING].get() &&
	    uniform_state_.toggle[GLUP_VERTEX_NORMALS].get()
	) {
            immediate_state_.buffer[GLUP_NORMAL_ATTRIBUTE].enable();
        } else {
            immediate_state_.buffer[GLUP_NORMAL_ATTRIBUTE].disable();
        }
        
        immediate_state_.begin(primitive);
        prepare_to_draw(primitive);
        configure_OpenGL_texturing();
        configure_OpenGL_lighting();
        configure_OpenGL_clipping();        
        configure_OpenGL_picking();
    }

    void Context_VanillaGL::prepare_to_draw(GLUPprimitive primitive) {
        if(primitive == GLUP_POINTS) {
            glPointSize(uniform_state_.point_size.get());
        }
    }
    
    void Context_VanillaGL::do_update_uniform_buffer() {
        copy_to_GL_state(GLUP_CLIPPING_ATTRIBUTES_BIT);
        copy_to_GL_state(GLUP_COLORS_ATTRIBUTES_BIT);
        Context::do_update_uniform_buffer();
    }
    
    void Context_VanillaGL::update_matrices() {
        if(matrices_dirty_) {
            Context::update_matrices();
            copy_to_GL_state(GLUP_MATRICES_ATTRIBUTES_BIT);
        }
    }

    void Context_VanillaGL::update_lighting() {
        if(lighting_dirty_) {
            Context::update_matrices();
            copy_to_GL_state(GLUP_LIGHTING_ATTRIBUTES_BIT);
        }
    }
    
    void Context_VanillaGL::configure_OpenGL_texturing() {
        glDisable(GLUP_TEXTURE_1D_TARGET);
        glDisable(GLUP_TEXTURE_2D_TARGET);
        glDisable(GLUP_TEXTURE_3D_TARGET);            
        if( !uniform_state_.toggle[GLUP_PICKING].get() &&
            uniform_state_.toggle[GLUP_TEXTURING].get()
        ) {
            if(uniform_state_.toggle[GLUP_INDIRECT_TEXTURING].get()) {
                begin_indirect_texturing();
            } else {
                switch(uniform_state_.texture_type.get()) {
                case GLUP_TEXTURE_1D:
                    glEnable(GLUP_TEXTURE_1D_TARGET);
                    break;
                case GLUP_TEXTURE_2D: {
                    //  Copy the 2D texture bound to unit 1 to unit 0
                    // else the fixed functionality pipeline cannot
                    // see it.
                    glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
                    GLint tex;
                    glGetIntegerv(GL_TEXTURE_BINDING_2D, &tex);
                    glActiveTexture(GL_TEXTURE0);
                    glBindTexture(GLUP_TEXTURE_2D_TARGET, GLuint(tex));
                    glEnable(GLUP_TEXTURE_2D_TARGET);
                } break;
                case GLUP_TEXTURE_3D: {
#ifdef GEO_GL_TEXTURE_3D                    
                    //  Copy the 3D texture bound to unit 2 to unit 0
                    // else the fixed functionality pipeline cannot
                    // see it.
                    glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_3D_UNIT);
                    GLint tex;
                    glGetIntegerv(GL_TEXTURE_BINDING_3D, &tex);
                    glActiveTexture(GL_TEXTURE0);
                    glBindTexture(GLUP_TEXTURE_3D_TARGET, GLuint(tex));
                    glEnable(GLUP_TEXTURE_3D_TARGET);
#endif                                    
                } break;
                }
            }

            glMatrixMode(GL_TEXTURE);
            glLoadMatrixf(get_matrix(GLUP_TEXTURE_MATRIX));
            glMatrixMode(GL_MODELVIEW);
            
            switch(uniform_state_.texture_mode.get()) {
            case GLUP_TEXTURE_REPLACE:
                // Yes, it's GL_MODULATE also for texture replace,
                // because GL_TEXTURE_REPLACE removes the shading !
                if(uniform_state_.toggle[GLUP_LIGHTING].get()) {
                    glTexEnvi(
                        GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE
                    );
                } else {
                    glTexEnvi(
                        GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE
                    );
                }
                break;
            case GLUP_TEXTURE_MODULATE:
                glTexEnvi(
                    GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE
                );
                break;
            case GLUP_TEXTURE_ADD:
                glTexEnvi(
                    GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_ADD
                );
                break;
            }
        }
    }
   
    void Context_VanillaGL::configure_OpenGL_lighting() {
        GLUPprimitive primitive = immediate_state_.primitive();

        if(primitive == GLUP_POINTS || primitive == GLUP_LINES) {
            glDisable(GL_LIGHTING);
        } else {
            if(uniform_state_.toggle[GLUP_LIGHTING].get()) {
                glEnable(GL_LIGHTING);
                glEnable(GL_NORMALIZE);
            } else {
                glDisable(GL_LIGHTING);                
            }
        }
        
        if(
            uniform_state_.toggle[GLUP_VERTEX_COLORS].get() ||
            primitive == GLUP_POINTS ||
            primitive == GLUP_LINES
        ) {
            glEnable(GL_COLOR_MATERIAL);
        } else {
            glEnable(GL_COLOR_MATERIAL);
            glColor3f(1.0f, 1.0f, 1.0f);
            glDisable(GL_COLOR_MATERIAL);            
        }

        
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
        glColor4fv(uniform_state_.color[GLUP_FRONT_COLOR].get_pointer());
        
        glMaterialfv(
            GL_FRONT, GL_DIFFUSE,
            uniform_state_.color[GLUP_FRONT_COLOR].get_pointer()
        );        
        glMaterialfv(
            GL_BACK, GL_DIFFUSE,
            uniform_state_.color[GLUP_BACK_COLOR].get_pointer()
        );

        static GLfloat ambient[4]  = { 0.2f, 0.2f, 0.2f, 1.0f };
        static GLfloat zero[4]     = { 0.0f, 0.0f, 0.0f, 0.0f };
	
        glMaterialfv(
            GL_FRONT_AND_BACK, GL_AMBIENT, ambient
        );
	if(uniform_state_.specular.get() != 0.0f) {
	    GLfloat specular[4];
	    specular[0] = uniform_state_.specular.get();
	    specular[1] = uniform_state_.specular.get();
	    specular[2] = uniform_state_.specular.get();
	    specular[3] = uniform_state_.specular.get();	    
	    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);        
	    glMaterialf(
		GL_FRONT_AND_BACK, GL_SHININESS, 30.0
	    );
	} else {
	    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, zero);
	    glMaterialf(
		GL_FRONT_AND_BACK, GL_SHININESS, 0.0
	    );
	}
        glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
    }

    void Context_VanillaGL::configure_OpenGL_picking() {
        if(uniform_state_.toggle[GLUP_PICKING].get()) {
            
            glEnable(GL_COLOR_MATERIAL);
            glDisable(GL_LIGHTING);
            
            // Disable colors and texture coordinates
            immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].disable();
            immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].disable();

            //   Disable buffers for interpolated clipping
            // attributes.
            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
            
            uniform_state_.base_picking_id.set(0);
            
            switch(uniform_state_.picking_mode.get()) {
            case GLUP_PICK_PRIMITIVE: {
                pick_primitives_ = true;
            } break;
            case GLUP_PICK_CONSTANT: {
                pick_primitives_ = false;                
                glPickingIdAsColor(index_t(uniform_state_.picking_id.get()));
            } break;
            }
        } else {
            pick_primitives_ = false;
        }
    }

    void Context_VanillaGL::configure_OpenGL_clipping() {
        if(clip_slice_cells()) {
            glNormal3f(
                -world_clip_plane_[0],
                -world_clip_plane_[1],
                -world_clip_plane_[2]            
            );
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
            
            // Bind the intersection buffers as vertex array, colors and texture
            // coordinates (if enabled for the last two ones).
            glEnableClientState(GL_VERTEX_ARRAY);
        
            glVertexPointer(
                4,
                GL_FLOAT,
                0,  // stride
                &(isect_vertex_attribute_[0][0])
            );
        
            if(immediate_state_.buffer[1].is_enabled()) {
                glEnableClientState(GL_COLOR_ARRAY);
                glColorPointer(
                    4,
                    GL_FLOAT,
                    0,  // stride
                    &(isect_vertex_attribute_[1][0])                    
                );
            }
        
            if(immediate_state_.buffer[2].is_enabled()) {
                glEnableClientState(GL_TEXTURE_COORD_ARRAY);
                glTexCoordPointer(
                    4,
                    GL_FLOAT,
                    0,  // stride
                    &(isect_vertex_attribute_[2][0])                    
                );
            }
        }
    }
    
    void Context_VanillaGL::end() {
        flush_immediate_buffers();
        if(clip_slice_cells()) {
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        }
        if(uniform_state_.toggle[GLUP_INDIRECT_TEXTURING].get()) {
            end_indirect_texturing();
        }
    }

    void Context_VanillaGL::begin_indirect_texturing() {

        // TODO: there is one texture matrix per texture unit,
        // make sure that we do what we meant to do (associate
        // the correct matrix with the correct unit)
        // Here, texture matrix is associated with unit 2, whereas
        // it acts on the colormap (needs to be fixed...)

        // TODO: indirect texturing is only implemented for 3D textures
        // for now. Implement it for 2D textures also if needed.

        // Note: indirect texturing uses a fragment shader, not very
        // "Vanilla", but not very "spicy" either (just an ARBfp1.0
        // shader, that most graphic boards should support without
        // any problem, it has been there since OpenGL 2.0).

        //   This program reads the 3D texture (bound to unit 2),
        // transforms the read value with the texture matrix,
        // and then reads the colormap (bound to unit 0) using
        // the transformed value.
        
        static const char* fshader_source =
            "!!ARBfp1.0                                    "
            "TEMP R0;                                      "
            "TEX R0, fragment.texcoord[0], texture[2], 3D; "
            "DP4 R0.x, state.matrix.texture[2].row[0], R0; "
            "TEX result.color, R0, texture[0], 2D;         "
            "END                                           "
            ;
        
        if(indirect_texturing_program_ == 0) {
            glGenProgramsARB(1, &indirect_texturing_program_) ;
            glEnable(GL_FRAGMENT_PROGRAM_ARB);
            glBindProgramARB(
                GL_FRAGMENT_PROGRAM_ARB, indirect_texturing_program_
            );
            glProgramStringARB(
                GL_FRAGMENT_PROGRAM_ARB, GL_PROGRAM_FORMAT_ASCII_ARB,
                GLsizei(strlen(fshader_source)), fshader_source
            );
            GLint errpos ;
            glGetIntegerv(GL_PROGRAM_ERROR_POSITION_ARB, &errpos) ;
            bool ok = ( errpos == -1) ;
            if(!ok) {
                const char* s = (const char*)(
                    glGetString(GL_PROGRAM_ERROR_STRING_ARB)
                );
                Logger::err("ARBfp")
                    << ":" << errpos << ": " << s << std::endl ;
            } 
        }
        glBindProgramARB(
            GL_FRAGMENT_PROGRAM_ARB, indirect_texturing_program_
        );
        glEnable(GL_FRAGMENT_PROGRAM_ARB);
    }

    void Context_VanillaGL::end_indirect_texturing() {
        glDisable(GL_FRAGMENT_PROGRAM_ARB);
    }
    
    void Context_VanillaGL::setup() {
        create_CPU_side_uniform_buffer();
    }



    void Context_VanillaGL::flush_immediate_buffers() {
        shrink_cells_in_immediate_buffers();
        classify_vertices_in_immediate_buffers();        

        glPushAttrib(
            GL_ENABLE_BIT | GL_LIGHTING_BIT |
            GL_POLYGON_BIT | GL_TEXTURE_BIT |
            GL_CURRENT_BIT
        );
        
        if(
            uniform_state_.clipping_mode.get() != GLUP_CLIP_STANDARD &&
            immediate_state_.primitive() != GLUP_POINTS &&
            immediate_state_.primitive() != GLUP_LINES &&
            immediate_state_.primitive() != GLUP_TRIANGLES &&
            immediate_state_.primitive() != GLUP_QUADS                       
        ) {
            glDisable(GL_CLIP_PLANE0);
        }

        flush_immediate_buffers_once();

        if(uniform_state_.toggle[GLUP_PICKING].get()) {
            if(pick_primitives_) {
                uniform_state_.base_picking_id.set(
                    uniform_state_.base_picking_id.get() +
                    int(immediate_state_.nb_primitives())
                );
            }
        } else if(
            uniform_state_.toggle[GLUP_DRAW_MESH].get() &&
            immediate_state_.primitive() != GLUP_POINTS &&
            immediate_state_.primitive() != GLUP_LINES  &&
	    immediate_state_.primitive() != GLUP_SPHERES
        ) {
            // Do it one more time for the mesh
            
            glDisable(GL_LIGHTING);
            glDisable(GLUP_TEXTURE_1D_TARGET);
            glDisable(GLUP_TEXTURE_2D_TARGET);
            glDisable(GLUP_TEXTURE_3D_TARGET);
            glDisable(GL_COLOR_MATERIAL);
            
            // Disable vertex attributes.
            bool va1_enabled = immediate_state_.buffer[1].is_enabled();
            bool va2_enabled = immediate_state_.buffer[2].is_enabled();
            immediate_state_.buffer[1].disable();
            immediate_state_.buffer[2].disable();

            if(clip_slice_cells()) {
                glDisableClientState(GL_COLOR_ARRAY);
                glDisableClientState(GL_TEXTURE_COORD_ARRAY);
            }
            
            glColor3fv(uniform_state_.color[GLUP_MESH_COLOR].get_pointer());
            glLineWidth(uniform_state_.mesh_width.get());
            glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

            flush_immediate_buffers_once();

            // Restore previous state for vertex attributes.
            if(va1_enabled) {
                immediate_state_.buffer[1].enable();
                if(clip_slice_cells()) {
                    glEnableClientState(GL_COLOR_ARRAY);                    
                }
            }
            if(va2_enabled) {
                immediate_state_.buffer[2].enable();
                if(clip_slice_cells()) {
                    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
                }
            }
        }

        glColor3f(1.0f, 1.0f, 1.0f);
        glPopAttrib();        

        immediate_state_.reset();        
    }
    
    void Context_VanillaGL::flush_immediate_buffers_once() {
        switch(immediate_state_.primitive()) {
        case GLUP_POINTS:
            draw_immediate_buffer_GLUP_POINTS();
            break;
        case GLUP_LINES:
            draw_immediate_buffer_GLUP_LINES();
            break;
        case GLUP_TRIANGLES: 
            draw_immediate_buffer_GLUP_TRIANGLES();
            break;
        case GLUP_QUADS:
            draw_immediate_buffer_GLUP_QUADS();
            break;
        case GLUP_TETRAHEDRA:
            draw_immediate_buffer_GLUP_TETRAHEDRA();
            break;
        case GLUP_HEXAHEDRA:
            draw_immediate_buffer_GLUP_HEXAHEDRA();
            break;
        case GLUP_PRISMS:
            draw_immediate_buffer_GLUP_PRISMS();
            break;
        case GLUP_PYRAMIDS:
            draw_immediate_buffer_GLUP_PYRAMIDS();
            break;
        case GLUP_CONNECTORS:
            draw_immediate_buffer_GLUP_CONNECTORS();
            break;
	case GLUP_SPHERES:
            draw_immediate_buffer_GLUP_SPHERES();
            break;
        case GLUP_NB_PRIMITIVES:
            geo_assert_not_reached;
        }
    }

    void Context_VanillaGL::setup_immediate_buffers() {
    }

    void Context_VanillaGL::setup_primitives() {
        primitive_info_.resize(GLUP_NB_PRIMITIVES);
        primitive_info_[GLUP_POINTS].implemented = true;
        primitive_info_[GLUP_LINES].implemented = true;
        primitive_info_[GLUP_TRIANGLES].implemented = true;
        primitive_info_[GLUP_QUADS].implemented = true;
        primitive_info_[GLUP_TETRAHEDRA].implemented = true;
        primitive_info_[GLUP_HEXAHEDRA].implemented = true;
        primitive_info_[GLUP_PYRAMIDS].implemented = true;
        primitive_info_[GLUP_PRISMS].implemented = true;
        primitive_info_[GLUP_CONNECTORS].implemented = true;
        primitive_info_[GLUP_SPHERES].implemented = true;	
    }
    
    Memory::pointer Context_VanillaGL::get_state_variable_address(
        const char* name
    ) {
        geo_assert(variable_to_offset_.find(name) != variable_to_offset_.end());
        return uniform_buffer_data_ + variable_to_offset_[name];
    }
    
    bool Context_VanillaGL::primitive_supports_array_mode(
        GLUPprimitive prim
    ) const {
        geo_argused(prim);
        return false;
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_POINTS() {
        glDisable(GL_LIGHTING);
        glEnable(GL_POINT_SMOOTH);
        glBegin(GL_POINTS);
        for(index_t v=0; v<immediate_state_.nb_vertices(); ++v) {
            output_picking_id(v);
            output_vertex(v);
        }
        glEnd();
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_LINES() {
        glBegin(GL_LINES);
        index_t v = 0;
        while(v < immediate_state_.nb_vertices()) {
            output_picking_id(v/2);
            output_vertex(v);
            output_vertex(v+1);
            v += 2;
        }
        glEnd();
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_TRIANGLES() {
        glBegin(GL_TRIANGLES);
        index_t v = 0;
        while(v < immediate_state_.nb_vertices()) {
            output_picking_id(v/3);
            draw_triangle(v,v+1,v+2);
            v += 3;
        }
        glEnd();
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_QUADS() {
        glBegin(GL_QUADS);
        index_t v = 0;
        while(v < immediate_state_.nb_vertices()) {
            output_picking_id(v/4);
            draw_quad(v,v+1,v+3,v+2);
            v += 4;
        }
        glEnd();
    }


    void Context_VanillaGL::draw_immediate_buffer_with_marching_cells(
        const MarchingCell& cell
    ) {
        
        index_t v0=0;
        while(v0 < immediate_state_.nb_vertices()) {
            index_t config = get_config(v0, cell.nb_vertices());

            //   Compute all the intersection vertices (plus their
            // attributes if enabled).
            for(index_t i=0; i<cell.config_size(config); ++i) {
                index_t e = cell.config_edges(config)[i];
                compute_intersection(
                    v0+cell.edge_vertex(e,0), v0+cell.edge_vertex(e,1), e
                );
            }            


            if(cell.config_size(config) != 0) {
                output_picking_id(v0/cell.nb_vertices());
                
                //   With the bound intersection buffer,
                // we can draw the intersection polygon with
                // a single OpenGL call ! The marching cells table directly
                // refers to the vertices in the intersection buffer.
                glDrawElements(
                    GL_POLYGON,
                    GLsizei(cell.config_size(config)),
                    GL_UNSIGNED_INT,
                    cell.config_edges(config)
                );
            }
            
            v0 += cell.nb_vertices();
        }

    }
    
    void Context_VanillaGL::draw_immediate_buffer_GLUP_TETRAHEDRA() {
        if(clip_slice_cells()) {
            draw_immediate_buffer_with_marching_cells(marching_tet_);
        } else {
            glBegin(GL_TRIANGLES);
            index_t v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {            
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/4);                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    draw_triangle(v0,v1,v2);                           
                    draw_triangle(v1,v0,v3);                             
                    draw_triangle(v0,v2,v3);                             
                    draw_triangle(v2,v1,v3);
                }
                v0 += 4;
            }
            glEnd();
        }
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_HEXAHEDRA() {
        if(clip_slice_cells()) {
            draw_immediate_buffer_with_marching_cells(marching_hex_);
        } else {
            glBegin(GL_QUADS);
            index_t v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/8);                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    index_t v4 = v0+4;
                    index_t v5 = v0+5;
                    index_t v6 = v0+6;
                    index_t v7 = v0+7;
                    draw_quad(v0,v2,v4,v6);
                    draw_quad(v3,v1,v7,v5);
                    draw_quad(v1,v0,v5,v4);
                    draw_quad(v2,v3,v6,v7);
                    draw_quad(v1,v3,v0,v2);
                    draw_quad(v4,v6,v5,v7);
                }
                v0 += 8;
            }
            glEnd();
        }
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_PRISMS() {
        if(clip_slice_cells()) {
            draw_immediate_buffer_with_marching_cells(marching_prism_);
        } else {
            glBegin(GL_QUADS);
            index_t v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/6);                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    index_t v4 = v0+4;
                    index_t v5 = v0+5;
                    draw_quad(v0,v3,v1,v4);
                    draw_quad(v0,v2,v3,v5);
                    draw_quad(v1,v4,v2,v5);
                }
                v0 += 6;
            }
            glEnd();
            glBegin(GL_TRIANGLES);
            v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/6);                                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    index_t v4 = v0+4;
                    index_t v5 = v0+5;
                    draw_triangle(v0,v1,v2);
                    draw_triangle(v5,v4,v3);
                }
                v0 += 6;
            }
            glEnd();
        }
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_PYRAMIDS() {
        if(clip_slice_cells()) {
            draw_immediate_buffer_with_marching_cells(marching_pyramid_);
        } else {
            glBegin(GL_QUADS);
            index_t v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/5);                                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    draw_quad(v0,v1,v3,v2);
                }
                v0 += 5;
            }
            glEnd();
            glBegin(GL_TRIANGLES);
            v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/5);                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    index_t v4 = v0+4;
                    draw_triangle(v0,v4,v1);
                    draw_triangle(v0,v3,v4);
                    draw_triangle(v2,v4,v3);
                    draw_triangle(v2,v1,v4);
                }
                v0 += 5;
            }
            glEnd();
        }
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_CONNECTORS() {
        if(clip_slice_cells()) {
            draw_immediate_buffer_with_marching_cells(marching_connector_);
        } else {
            glBegin(GL_QUADS);
            index_t v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/4);                                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    draw_quad(v0,v1,v3,v2);
                }
                v0 += 4;
            }
            glEnd();
            glBegin(GL_TRIANGLES);
            v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/4);                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    draw_triangle(v2,v1,v0);
                    draw_triangle(v3,v2,v0);
                }
                v0 += 4;
            }
            glEnd();
        }
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_SPHERES() {
	// Note: this is not very optimized. We could probably use
	// a form of instanced drawing supported by early OpenGL,
	// or use a vertex shader in an old VS dialect supported
	// everywhere, but anyway VanillaGL is seldom used, so I
	// will not spend too much on optimizing it (unless I
	// receive a query for that).
	// This would require re-implementing shading in the vertex
	// shader (not too difficult, but annoying ...)
	glMatrixMode(GL_MODELVIEW);
	for(index_t v=0; v<immediate_state_.nb_vertices(); ++v) {
            output_picking_id(v);	    
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
	    const float* p =
		immediate_state_.buffer[GLUP_VERTEX_ATTRIBUTE].element_ptr(v);
	    glPushMatrix();
	    glTranslatef(p[0], p[1], p[2]);
	    glScalef(p[3], p[3], p[3]);	    	    
	    glCallList(sphere_display_list_);
	    glPopMatrix();
	}
    }    
    
    void Context_VanillaGL::output_normal(index_t v1, index_t v2, index_t v3) {
        GLfloat* p1 = immediate_state_.buffer[0].element_ptr(v1);
        GLfloat* p2 = immediate_state_.buffer[0].element_ptr(v2);
        GLfloat* p3 = immediate_state_.buffer[0].element_ptr(v3);

        // scale vector components, else it can generate floating
        // point exceptions when manipulating very small vectors
        // (e.g. at the poles of the sphere generated in Graphite).
        const float s = 100.0f;
        
        GLfloat U[3];
        U[0] = s*(p2[0] - p1[0]);
        U[1] = s*(p2[1] - p1[1]);
        U[2] = s*(p2[2] - p1[2]);
        
        GLfloat V[3];            
        V[0] = s*(p3[0] - p1[0]);
        V[1] = s*(p3[1] - p1[1]);
        V[2] = s*(p3[2] - p1[2]);
        
        glNormal3f(
            U[1]*V[2] - U[2]*V[1],
            U[2]*V[0] - U[0]*V[2],
            U[0]*V[1] - U[1]*V[0]                
        );
    }


    void Context_VanillaGL::output_normal(
        index_t v1, index_t v2, index_t v3, index_t v4
    ) {
        GLfloat* p1 = immediate_state_.buffer[0].element_ptr(v1);
        GLfloat* p2 = immediate_state_.buffer[0].element_ptr(v2);
        GLfloat* p3 = immediate_state_.buffer[0].element_ptr(v3);
        GLfloat* p4 = immediate_state_.buffer[0].element_ptr(v4);        

        // scale vector components, else it can generate floating
        // point exceptions when manipulating very small vectors
        // (e.g. at the poles of the sphere generated in Graphite).
        const float s = 100.0f;
        
        GLfloat U1[3];
        U1[0] = s*(p2[0] - p1[0]);
        U1[1] = s*(p2[1] - p1[1]);
        U1[2] = s*(p2[2] - p1[2]);
        
        GLfloat V1[3];            
        V1[0] = s*(p4[0] - p1[0]);
        V1[1] = s*(p4[1] - p1[1]);
        V1[2] = s*(p4[2] - p1[2]);

        GLfloat U2[3];
        U2[0] = s*(p4[0] - p3[0]);
        U2[1] = s*(p4[1] - p3[1]);
        U2[2] = s*(p4[2] - p3[2]);
        
        GLfloat V2[3];            
        V2[0] = s*(p2[0] - p3[0]);
        V2[1] = s*(p2[1] - p3[1]);
        V2[2] = s*(p2[2] - p3[2]);
        
        glNormal3f(
            (U1[1]*V1[2]-U1[2]*V1[1]) - (U2[1]*V2[2]-U2[2]*V2[1]),
            (U1[2]*V1[0]-U1[0]*V1[2]) - (U2[2]*V2[0]-U2[0]*V2[2]),
            (U1[0]*V1[1]-U1[1]*V1[0]) - (U2[0]*V2[1]-U2[1]*V2[0])               
        );
    }

}

#endif
