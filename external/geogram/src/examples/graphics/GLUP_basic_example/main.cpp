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

#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/GLUP/GLUP.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include "uv.xpm"

namespace {

    bool lighting = true;
    
    /**
     * \brief Zooms in.
     * \details Zooming factor is 1.1x.
     */
    void zoom_in() {
        *glup_viewer_float_ptr(GLUP_VIEWER_ZOOM) *= 1.1f;
    }

    /**
     * \brief Zooms out.
     * \details De-zooming factor is (1/1.1)x.
     */
    void zoom_out() {
        *glup_viewer_float_ptr(GLUP_VIEWER_ZOOM) /= 1.1f;
    }

    GEO::index_t n = 30;
    
    GLint prim = GLUP_POINTS;
    
    bool VBO_mode = false;
    GLuint index_VBO = 0;
    GLuint vertex_VBO = 0;
    GLuint color_VBO = 0;
    GLuint tex_coord_VBO = 0;
//  GLuint VAO = 0; // TODO...
    GLsizei VBO_nb_vertices = 0;
    GLuint texture;

    void reset_VBOs() {
        if(!VBO_mode) {
            return;
        }
        if(index_VBO != 0) {
            glDeleteBuffers(1, &index_VBO);
        }
        index_VBO = 0;
        if(vertex_VBO != 0) {
            glDeleteBuffers(1, &vertex_VBO);
        }
        vertex_VBO = 0;
        if(color_VBO != 0) {
            glDeleteBuffers(1, &color_VBO);
        }
        color_VBO = 0;
        if(tex_coord_VBO != 0) {
            glDeleteBuffers(1, &tex_coord_VBO);
        }
        tex_coord_VBO = 0;
        VBO_mode = false;
        GEO::Logger::out("VBO") << "Deactivated" << std::endl;        
    }

    void setup_VBO_vertices_grid() {
        GEO::index_t nb_vertices = n*n*n;
        size_t size = sizeof(float)*3*nb_vertices;
        float* data = new float[3*nb_vertices];
        float* ptr = data;
        for(GEO::index_t i=0; i<n; ++i) {
            for(GEO::index_t j=0; j<n; ++j) {
                for(GEO::index_t k=0; k<n; ++k) {
                    ptr[0] = float(i)/float(n);
                    ptr[1] = float(j)/float(n);
                    ptr[2] = float(k)/float(n);
                    ptr += 3;
                }
            }
        }
        GEO::update_buffer_object(vertex_VBO, GL_ARRAY_BUFFER, size, data);
        GEO::update_buffer_object(
            color_VBO, GL_ARRAY_BUFFER, size, data
        );
        GEO::update_buffer_object(
            tex_coord_VBO, GL_ARRAY_BUFFER, size, data
        );
        delete[] data;
    }

    void setup_VBO_vertices_sphere() {
        GEO::index_t nb_vertices = n*n;
        size_t size = sizeof(float)*3*nb_vertices;
        float* data = new float[3*nb_vertices];
        float* ptr = data;
        for(GEO::index_t i=0; i<n; ++i) {
            for(GEO::index_t j=0; j<n; ++j) {
                double theta = double(i) * 2.0 * M_PI / double(n-1);
                double phi = -M_PI/2.0 + double(j) * M_PI / double(n-1);
                
                double x = (cos(theta)*cos(phi) + 1.0)/2.0;
                double y = (sin(theta)*cos(phi) + 1.0)/2.0;
                double z = (sin(phi) + 1.0) / 2.0;
                
                ptr[0] = float(x);
                ptr[1] = float(y);
                ptr[2] = float(z);
                ptr += 3;
            }
        }
        GEO::update_buffer_object(vertex_VBO, GL_ARRAY_BUFFER, size, data);
        GEO::update_buffer_object(
            color_VBO, GL_ARRAY_BUFFER, size, data
        );
        GEO::update_buffer_object(
            tex_coord_VBO, GL_ARRAY_BUFFER, size, data
        );
        delete[] data;
    }

    
    void setup_VBO_GLUP_POINTS() {
        VBO_mode = true;
        setup_VBO_vertices_grid();
        VBO_nb_vertices = GLsizei(n*n*n);
    }


    void setup_VBO_GLUP_TRIANGLES() {
        VBO_mode = true;
        setup_VBO_vertices_sphere();
        GEO::index_t nb_triangles = (n-1)*(n-1)*2;
        GLuint* data = new GLuint[nb_triangles*3];
        GLuint* ptr = data;
        for(GEO::index_t i=0; i<n-1; ++i) {
            for(GEO::index_t j=0; j<n-1; ++j) {
                GEO::index_t idx00 = i*n+j;
                GEO::index_t idx01 = i*n+j+1;
                GEO::index_t idx10 = (i+1)*n+j;
                GEO::index_t idx11 = (i+1)*n+j+1;
                ptr[0] = GLuint(idx00);
                ptr[1] = GLuint(idx11);
                ptr[2] = GLuint(idx01);
                ptr[3] = GLuint(idx11);
                ptr[4] = GLuint(idx00);
                ptr[5] = GLuint(idx10);                    
                ptr += 6;
            }
        }
        GEO::update_buffer_object(
            index_VBO, GL_ELEMENT_ARRAY_BUFFER,
            sizeof(GLuint)*3*nb_triangles, data
        );
        delete[] data;
        VBO_nb_vertices = GLsizei(nb_triangles * 3);
    }

    void setup_VBO_GLUP_QUADS() {
        VBO_mode = true;
        setup_VBO_vertices_sphere();        
        GEO::index_t nb_quads = (n-1)*(n-1);
        GLuint* data = new GLuint[nb_quads*4];
        GLuint* ptr = data;
        for(GEO::index_t i=0; i<n-1; ++i) {
            for(GEO::index_t j=0; j<n-1; ++j) {
                GEO::index_t idx00 = i*n+j;
                GEO::index_t idx01 = i*n+j+1;
                GEO::index_t idx10 = (i+1)*n+j;
                GEO::index_t idx11 = (i+1)*n+j+1;
                ptr[0] = GLuint(idx00);
                ptr[1] = GLuint(idx10);
                ptr[2] = GLuint(idx11);
                ptr[3] = GLuint(idx01);
                ptr += 4;
            }
        }
        GEO::update_buffer_object(
            index_VBO, GL_ELEMENT_ARRAY_BUFFER,
            sizeof(GLuint)*4*nb_quads, data
        );
        delete[] data;
        VBO_nb_vertices = GLsizei(nb_quads * 4);
    }

    void setup_VBO_GLUP_TETRAHEDRA() {
        VBO_mode = true;
        setup_VBO_vertices_grid();
        GEO::index_t nb_tets = (n-1)*(n-1)*(n-1);
        GLuint* data = new GLuint[nb_tets*4];
        GLuint* ptr = data;
        for(GEO::index_t i=0; i<n-1; ++i) {
            for(GEO::index_t j=0; j<n-1; ++j) {
                for(GEO::index_t k=0; k<n-1; ++k) {                
                    GEO::index_t idx000 = (i  )*n*n+(j  )*n+(k  );
                    GEO::index_t idx100 = (i+1)*n*n+(j  )*n+(k  );
                    GEO::index_t idx010 = (i  )*n*n+(j+1)*n+(k  );
                    GEO::index_t idx001 = (i  )*n*n+(j  )*n+(k+1);
                    ptr[0] = GLuint(idx000);
                    ptr[1] = GLuint(idx001);
                    ptr[2] = GLuint(idx010);
                    ptr[3] = GLuint(idx100);
                    ptr += 4;
                }
            }
        }
        GEO::update_buffer_object(
            index_VBO, GL_ELEMENT_ARRAY_BUFFER,
            sizeof(GLuint)*4*nb_tets, data
        );
        delete[] data;
        VBO_nb_vertices = GLsizei(nb_tets * 4);
    }

    void setup_VBO_GLUP_HEXAHEDRA() {
        VBO_mode = true;
        setup_VBO_vertices_grid();
        GEO::index_t nb_hex = (n-1)*(n-1)*(n-1);
        GLuint* data = new GLuint[nb_hex*8];
        GLuint* ptr = data;
        for(GEO::index_t i=0; i<n-1; ++i) {
            for(GEO::index_t j=0; j<n-1; ++j) {
                for(GEO::index_t k=0; k<n-1; ++k) {                
                    GEO::index_t idx000 = (i  )*n*n+(j  )*n+(k  );
                    GEO::index_t idx001 = (i  )*n*n+(j  )*n+(k+1);
                    GEO::index_t idx010 = (i  )*n*n+(j+1)*n+(k  );
                    GEO::index_t idx011 = (i  )*n*n+(j+1)*n+(k+1);
                    GEO::index_t idx100 = (i+1)*n*n+(j  )*n+(k  );
                    GEO::index_t idx101 = (i+1)*n*n+(j  )*n+(k+1);
                    GEO::index_t idx110 = (i+1)*n*n+(j+1)*n+(k  );
                    GEO::index_t idx111 = (i+1)*n*n+(j+1)*n+(k+1);
                    ptr[0] = GLuint(idx000);
                    ptr[1] = GLuint(idx010);
                    ptr[2] = GLuint(idx100);
                    ptr[3] = GLuint(idx110);
                    ptr[4] = GLuint(idx001);
                    ptr[5] = GLuint(idx011);
                    ptr[6] = GLuint(idx101);
                    ptr[7] = GLuint(idx111);
                    ptr += 8;
                }
            }
        }
        GEO::update_buffer_object(
            index_VBO, GL_ELEMENT_ARRAY_BUFFER,
            sizeof(GLuint)*8*nb_hex, data
        );
        delete[] data;
        VBO_nb_vertices = GLsizei(nb_hex * 8);
    }
    
    void toggle_VBOs() {
        if(VBO_mode) {
            reset_VBOs();
            return;
        }

        if(!glupPrimitiveSupportsArrayMode(GLUPprimitive(prim))) {
            GEO::Logger::out("VBO") << "Not supported for current primitive"
                                    << std::endl;
            return;
        }
        
        GEO::Logger::out("VBO") << "Activating" << std::endl;
        switch(prim) {
        case GLUP_POINTS: {
            setup_VBO_GLUP_POINTS();
        } break;
        case GLUP_TRIANGLES: {
            setup_VBO_GLUP_TRIANGLES();
        } break;
        case GLUP_QUADS: {
            setup_VBO_GLUP_QUADS();
        } break;
        case GLUP_TETRAHEDRA: {
            setup_VBO_GLUP_TETRAHEDRA();
        } break;
        case GLUP_HEXAHEDRA: {
            setup_VBO_GLUP_HEXAHEDRA();
        }
        };
        GEO::Logger::out("VBO") << "Activated" << std::endl;        
    }
    
    
    /**
     * \brief Toggles mesh display.
     */
    void toggle_mesh() {
        if(glupIsEnabled(GLUP_DRAW_MESH)) {
            glupDisable(GLUP_DRAW_MESH);
        } else {
            glupEnable(GLUP_DRAW_MESH);
        }
    }

    /**
     * \brief Increases cell shrinking factor.
     */
    void inc_shrink() {
        glupSetCellsShrink(glupGetCellsShrink()+0.1f);
    }

    /**
     * \brief Decreases cell shrinking factor.
     */
    void dec_shrink() {
        glupSetCellsShrink(glupGetCellsShrink()-0.1f);
    }
    
    /**
     * \brief Toggles lighting / constant color mode.
     */
    void toggle_lighting() {
        lighting = !lighting;
    }

    void toggle_picking() {
        if(glupIsEnabled(GLUP_PICKING)) {
            glupDisable(GLUP_PICKING);
        } else {
            glupEnable(GLUP_PICKING);
            glupPickingMode(GLUP_PICK_PRIMITIVE);
        }
    }

    void toggle_texturing() {
        if(glupIsEnabled(GLUP_TEXTURING)) {
            glupDisable(GLUP_TEXTURING);
        } else {
            glupEnable(GLUP_TEXTURING);
        }
    }

    void cycle_texturing_mode() {
        static GLUPtextureMode mode = GLUP_TEXTURE_REPLACE;
        mode = GLUPtextureMode(mode + 1);
        if(mode > GLUP_TEXTURE_ADD) {
            mode = GLUP_TEXTURE_REPLACE;
        }
        glupTextureMode(mode);
    }

    void cycle_clipping_mode() {
        static GLUPclipMode mode = GLUP_CLIP_WHOLE_CELLS;
        mode = GLUPclipMode(mode + 1);
        if(mode > GLUP_CLIP_SLICE_CELLS) {
            mode = GLUP_CLIP_STANDARD;
        }
        glupClipMode(mode);
    }
    
    /**
     * \brief Toggles per-vertex colors.
     */
    void toggle_vertex_colors() {
        if(glupIsEnabled(GLUP_VERTEX_COLORS)) {
            glupDisable(GLUP_VERTEX_COLORS);
        } else {
            glupEnable(GLUP_VERTEX_COLORS);
        }
    }

    /**
     * \brief Toggles per-vertex normals.
     */
    void toggle_vertex_normals() {
        if(glupIsEnabled(GLUP_VERTEX_NORMALS)) {
            glupDisable(GLUP_VERTEX_NORMALS);
        } else {
            glupEnable(GLUP_VERTEX_NORMALS);
        }
    }

    /**
     * \brief Toggles fragment discard.
     */
    void toggle_alpha_discard() {
        if(glupIsEnabled(GLUP_ALPHA_DISCARD)) {
            glupDisable(GLUP_ALPHA_DISCARD);
        } else {
            glupEnable(GLUP_ALPHA_DISCARD);
        }
    }
    
    void inc_n() {
        reset_VBOs();
        GEO::index_t new_n = n * 3;
        if(new_n < 100000000) {
            n = new_n;
        }
        GEO::Logger::out("GLUP") << "n  = " << n << std::endl;
        GEO::Logger::out("GLUP") << "n2 = " << n*n << std::endl;
        GEO::Logger::out("GLUP") << "n3 = " << n*n*n << std::endl;        
    }

    void dec_n() {
        reset_VBOs();
        GEO::index_t new_n = n / 3;
        if(new_n >= 2) {
            n = new_n;
        }
        GEO::Logger::out("GLUP") << "n  = " << n << std::endl;
        GEO::Logger::out("GLUP") << "n2 = " << n*n << std::endl;        
        GEO::Logger::out("GLUP") << "n3 = " << n*n*n << std::endl;
    }

    inline void draw_vertex_grid(
        GEO::index_t i, GEO::index_t j, GEO::index_t k, float R=0.0f
    ) {
	glupColor4f(
	    float(i)/float(n),
	    float(j)/float(n),
	    float(k)/float(n),
	    0.0f
	);
        glupTexCoord3f(
            float(i)/float(n),
            float(j)/float(n),
            float(k)/float(n)                        
        );
        glupVertex4f(
            float(i)/float(n),
            float(j)/float(n),
            float(k)/float(n),
	    (R == 0.0f) ? 1.0f : R
        );
    }

    inline void draw_vertex_sphere(GEO::index_t i, GEO::index_t j) {
        double theta = double(i) * 2.0 * M_PI / double(n-1);
        double phi = -M_PI/2.0 + double(j) * M_PI / double(n-1);
        
        double x = (cos(theta)*cos(phi) + 1.0)/2.0;
        double y = (sin(theta)*cos(phi) + 1.0)/2.0;
        double z = (sin(phi) + 1.0) / 2.0;

        glupColor4d(x,y,z,0.0);
        glupTexCoord3d(x,y,z);
	glupNormal3d(x-0.5,y-0.5,z-0.5);
        glupVertex3d(x,y,z);
    }
    
    GLint point_size = 4;

    void inc_point_size() {
        ++point_size;
    }

    void dec_point_size() {
        --point_size;
        if(point_size <= 0) {
            point_size = 0;
        }
    }
    
    void next_primitive() {
        reset_VBOs();        
        ++prim;

	if(prim == GLUP_CONNECTORS) {
	    prim = GLUP_SPHERES;
	}
	
        if(prim > GLUP_SPHERES) {
            prim = GLUP_POINTS;
        }
    }
    
    /**
     * \brief Initializes OpenGL objects.
     * \details Specifed as glup_viewer_set_init_func() callback.
     */
    void init() {
        GEO::Graphics::initialize();
        
        glup_viewer_set_background_color(1.0, 1.0, 1.0);
        glup_viewer_add_toggle(
            'T', glup_viewer_is_enabled_ptr(GLUP_VIEWER_TWEAKBARS),
            "Toggle tweakbars"
        );
        glup_viewer_add_key_func('z', zoom_in, "Zoom in");
        glup_viewer_add_key_func('Z', zoom_out, "Zoom out");
        glup_viewer_disable(GLUP_VIEWER_TWEAKBARS);
        glup_viewer_disable(GLUP_VIEWER_BACKGROUND);
        glup_viewer_add_key_func('m', toggle_mesh, "mesh");

        glupEnable(GLUP_VERTEX_COLORS);
        glupEnable(GLUP_DRAW_MESH);

        glupMatrixMode(GLUP_TEXTURE_MATRIX);
        glupLoadIdentity();

        glGenTextures(1, &texture);
        glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2DXPM(uv);

        glupTextureType(GLUP_TEXTURE_2D);
        glupTextureMode(GLUP_TEXTURE_REPLACE);
        glupClipMode(GLUP_CLIP_WHOLE_CELLS);
    }

    void display_points() {
        glupSetPointSize(GLfloat(point_size));

        if(VBO_mode) {
            glBindBuffer(GL_ARRAY_BUFFER, vertex_VBO);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

            if(glupIsEnabled(GLUP_VERTEX_COLORS)) {
                glBindBuffer(GL_ARRAY_BUFFER, color_VBO);
                glEnableVertexAttribArray(1);
                glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
            }

            if(glupIsEnabled(GLUP_TEXTURING)) {
                glBindBuffer(GL_ARRAY_BUFFER, tex_coord_VBO);
                glEnableVertexAttribArray(2);
                glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
            }
            
            glupDrawArrays(GLUP_POINTS, 0, VBO_nb_vertices);
            
            glDisableVertexAttribArray(0);
            glDisableVertexAttribArray(1);
            glDisableVertexAttribArray(2);            
            glBindBuffer(GL_ARRAY_BUFFER,0);
            
        } else {
            glupBegin(GLUP_POINTS);
            for(GEO::index_t i=0; i<n; ++i) {
                for(GEO::index_t j=0; j<n; ++j) {
                    for(GEO::index_t k=0; k<n; ++k) {
                        draw_vertex_grid(i,j,k);
                    }
                }
            }
            glupEnd();
        }
    }

    void display_spheres() {
	glupBegin(GLUP_SPHERES);
	for(GEO::index_t i=0; i<n; i+=4) {
	    for(GEO::index_t j=0; j<n; j+=4) {
		for(GEO::index_t k=0; k<n; k+=4) {
		    double waved =
			sin(3.0*double(i)/double(n)*2.0*M_PI)*
			sin(3.0*double(j)/double(n)*2.0*M_PI)*
			sin(3.0*double(k)/double(n)*2.0*M_PI);
		    float wave = float(0.5*(1.0 + waved));
		    float Rmin = 0.2f * float(point_size) / float(n);
		    float Rmax = 0.8f * float(point_size) / float(n);
		    float R = wave*Rmin + (1.0f-wave)*Rmax;
		    draw_vertex_grid(i,j,k,R);
		}
	    }
	}
	glupEnd();
    }

    
    void display_lines() {
        glLineWidth(1);
        glupBegin(GLUP_LINES);
        for(GEO::index_t i=0; i<n-1; ++i) {
            for(GEO::index_t j=0; j<n-1; ++j) {
                for(GEO::index_t k=0; k<n-1; ++k) {
                    draw_vertex_grid(i,j,k);
                    draw_vertex_grid(i+1,j,k);

                    draw_vertex_grid(i,j,k);
                    draw_vertex_grid(i,j+1,k);                    

                    draw_vertex_grid(i,j,k);
                    draw_vertex_grid(i,j,k+1);                    
                }
            }
        }
        glupEnd();
    }


    void display_VBO_with_index(GLUPprimitive prim_in) {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_VBO);            
        glBindBuffer(GL_ARRAY_BUFFER, vertex_VBO);
        
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        
        if(glupIsEnabled(GLUP_VERTEX_COLORS)) {
            glBindBuffer(GL_ARRAY_BUFFER, color_VBO);
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        }

        if(glupIsEnabled(GLUP_TEXTURING)) {
            glBindBuffer(GL_ARRAY_BUFFER, tex_coord_VBO);
            glEnableVertexAttribArray(2);
            glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        }
        
        glupDrawElements(
            prim_in, VBO_nb_vertices, GL_UNSIGNED_INT, nullptr
        );
            
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glDisableVertexAttribArray(2);        
        glBindBuffer(GL_ARRAY_BUFFER,0);
    }
    
    void display_triangles() {

        if(VBO_mode) {
            display_VBO_with_index(GLUP_TRIANGLES);            
        } else {
            glupBegin(GLUP_TRIANGLES);            
            for(GEO::index_t i=0; i<n-1; ++i) {
                for(GEO::index_t j=0; j<n-1; ++j) {
                    draw_vertex_sphere(i,j);
                    draw_vertex_sphere(i+1,j+1);                
                    draw_vertex_sphere(i,j+1);
                    
                    draw_vertex_sphere(i+1,j+1);                
                    draw_vertex_sphere(i,j);                
                    draw_vertex_sphere(i+1,j);                
                }
            }
            glupEnd();            
        }
    }


    void display_quads() {

        if(VBO_mode) {
            display_VBO_with_index(GLUP_QUADS);            
        } else {
            glupBegin(GLUP_QUADS);
            
            for(GEO::index_t i=0; i<n-1; ++i) {
                for(GEO::index_t j=0; j<n-1; ++j) {
                    draw_vertex_sphere(i,j);
                    draw_vertex_sphere(i+1,j);
                    draw_vertex_sphere(i+1,j+1);                
                    draw_vertex_sphere(i,j+1);                    
                }
            }

            glupEnd();            
        }
        
    }

    void display_tetrahedra() {
        if(VBO_mode) {
            display_VBO_with_index(GLUP_TETRAHEDRA);            
        } else {
            glupBegin(GLUP_TETRAHEDRA);
            for(GEO::index_t i=0; i<n-1; ++i) {
                for(GEO::index_t j=0; j<n-1; ++j) {
                    for(GEO::index_t k=0; k<n-1; ++k) {
                        draw_vertex_grid(i,j,k);
                        draw_vertex_grid(i,j,k+1);
                        draw_vertex_grid(i,j+1,k);              
                        draw_vertex_grid(i+1,j,k);
                    }
                }
            }
            glupEnd();            
        }
    }

    void display_hexahedra() {
        if(VBO_mode) {
            display_VBO_with_index(GLUP_HEXAHEDRA);            
        } else {
            glupBegin(GLUP_HEXAHEDRA);
            for(GEO::index_t i=0; i<n-1; ++i) {
                for(GEO::index_t j=0; j<n-1; ++j) {
                    for(GEO::index_t k=0; k<n-1; ++k) {
                        draw_vertex_grid(i,  j,    k);
                        draw_vertex_grid(i  ,j+1,  k);
                        draw_vertex_grid(i+1,j,    k);
                        draw_vertex_grid(i+1,j+1,  k);
                        draw_vertex_grid(i,  j,    k+1);
                        draw_vertex_grid(i  ,j+1,  k+1);
                        draw_vertex_grid(i+1,j,    k+1);
                        draw_vertex_grid(i+1,j+1,  k+1);
                    }
                }
            }
            glupEnd();            
        }
    }

    
    void display_prisms() {
        if(VBO_mode) {
            // display_VBO_with_index(GLUP_PRISMS);            
        } else {
            glupBegin(GLUP_PRISMS);
            for(GEO::index_t i=0; i<n-1; ++i) {
                for(GEO::index_t j=0; j<n-1; ++j) {
                    for(GEO::index_t k=0; k<n-1; ++k) {
                        draw_vertex_grid(i,  j,    k);
                        draw_vertex_grid(i  ,j+1,  k);
                        draw_vertex_grid(i+1,j  ,  k);
                        draw_vertex_grid(i,  j,    k+1);
                        draw_vertex_grid(i  ,j+1,  k+1);
                        draw_vertex_grid(i+1,j  ,  k+1);              
                    }
                }
            }
            glupEnd();            
        }
    }

    void display_pyramids() {
        if(VBO_mode) {
            // display_VBO_with_index(GLUP_HEXAHEDRA);            
        } else {
            glupBegin(GLUP_PYRAMIDS);
            for(GEO::index_t i=0; i<n-1; ++i) {
                for(GEO::index_t j=0; j<n-1; ++j) {
                    for(GEO::index_t k=0; k<n-1; ++k) {
                        draw_vertex_grid(i,  j,    k);
                        draw_vertex_grid(i  ,j+1,  k);
                        draw_vertex_grid(i+1,j+1,  k);
                        draw_vertex_grid(i+1,j  ,  k);
                        draw_vertex_grid(i,  j,    k+1);
                    }
                }
            }
            glupEnd();            
        }
    }

    
    /**
     * \brief Draws a test scene.
     * \details Specified as glup_viewer_set_display_func() callback.
     */
    void display() {
        static bool recurse = false;
        
        if(lighting) {
            glupEnable(GLUP_LIGHTING);
        } else {
            glupDisable(GLUP_LIGHTING);            
        }

        if(glupIsEnabled(GLUP_TEXTURING)) {
            glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
            glBindTexture(GL_TEXTURE_2D, texture);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glupTextureType(GLUP_TEXTURE_2D);
        }
        
        glupMatrixMode(GLUP_TEXTURE_MATRIX);
        glupLoadMatrixf(glup_viewer_get_light_matrix());
        
        glupSetColor4f(GLUP_FRONT_COLOR, 1.0f, 1.0f, 0.0f, 0.0f);
        glupSetColor4f(GLUP_BACK_COLOR,  1.0f, 0.0f, 1.0f, 0.0f);        

        if(!recurse)
        switch(prim) {
        case GLUP_POINTS:
            display_points();
            break;
        case GLUP_LINES:
            display_lines();
            break;
        case GLUP_TRIANGLES:
            display_triangles();
            break;
        case GLUP_QUADS:
            display_quads();
            break;
        case GLUP_TETRAHEDRA:
            display_tetrahedra();
            break;
        case GLUP_HEXAHEDRA:
            display_hexahedra();
            break;
        case GLUP_PRISMS:
            display_prisms();
            break;
        case GLUP_PYRAMIDS:
            display_pyramids();
            break;
	case GLUP_SPHERES:
	    display_spheres();
	    break;
        default:
            break;
        }
    }
}

static void overlay() {
    // Displays ImGui demo window.
    ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver);
    ImGui::ShowDemoWindow();
}

int main(int argc, char** argv) {

    GEO::initialize();
    GEO::Logger::instance()->set_quiet(false);

    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("gfx");    

    GEO::CmdLine::set_arg("sys:assert","abort");
    
    std::vector<std::string> filenames;
    if(!GEO::CmdLine::parse(argc, argv, filenames, "")) {
        return 1;
    }

    glup_viewer_set_region_of_interest(
        0.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 1.0f
    );

    glup_viewer_set_window_title(
        "GLUP eyecandy test"
    );
    glup_viewer_set_init_func(init);
    glup_viewer_set_display_func(display);
    glup_viewer_set_overlay_func(overlay);
    glup_viewer_add_key_func('L', toggle_lighting, "toggle lighting");
    glup_viewer_add_key_func('c', toggle_vertex_colors, "toggle vrtx colors");
    glup_viewer_add_key_func('x', dec_shrink, "unshrink cells");
    glup_viewer_add_key_func('w', inc_shrink, "shrink cells");
    glup_viewer_add_key_func('n', inc_n, "increment n");
    glup_viewer_add_key_func('N', dec_n, "decrement n");
    glup_viewer_add_key_func('p', inc_point_size, "increment point size");
    glup_viewer_add_key_func('P', dec_point_size, "decrement point size");
    glup_viewer_add_key_func(' ', next_primitive, "cycle GLUP primitives");
    glup_viewer_add_key_func('o', toggle_picking, "toggle picking");
    glup_viewer_add_key_func('v', toggle_VBOs, "create VBO");
    glup_viewer_add_key_func('t', toggle_texturing, "toggle texturing");
    glup_viewer_add_key_func('y', cycle_texturing_mode, "texturing mode");
    glup_viewer_add_key_func('C', cycle_clipping_mode, "clipping mode");    
    glup_viewer_add_key_func('W', toggle_vertex_normals, "toggle vrtx normals");
    glup_viewer_add_key_func('a', toggle_alpha_discard, "toggle alpha discard");
    
    if(GEO::CmdLine::get_arg_bool("gfx:full_screen")) {
       glup_viewer_enable(GLUP_VIEWER_FULL_SCREEN);
    }

    glup_viewer_main_loop(argc, argv);

    return 0;
}
