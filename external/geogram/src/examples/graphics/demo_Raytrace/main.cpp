/*
 *  Copyright (c) 2012-2014, Bruno Levy All rights reserved.
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

/*
 * GEOGRAM example program:
 * simple mesh raytracing using AABB tree.
 * GUI version of geogram/simple_raytrace
 * disclaimer: not the fastest, lighting model is ridiculous, there is no
 *  antialiasing etc... (just a demo program for the mesh AABB class).
 * mouse interaction is nearly unusable...
 * usage: geogram_demo_Raytrace
 *        geogram_demo_Raytrace meshfile.(obj|ply|mesh|...)
 */ 

#include <geogram_gfx/glup_viewer/glup_viewer.h>
#define RAYTRACE_GUI
#include "../../geogram/simple_raytrace/raytracing.h"

namespace {
    using namespace GEO;

    
    /**
     * \brief An application that demonstrates both
     *  GLUP primitives and glup_viewer application
     *  framework.
     */
    class DemoRaytraceApplication : public Application {
    public:

        /**
         * \brief DemoRaytraceApplication constructor.
         */
        DemoRaytraceApplication(
            int argc, char** argv,
            const std::string& usage
	) :
	    Application(argc, argv, usage),
	    camera_(
		vec3(2.0, 2.0, 1.5), // Position
		vec3(0.5, 0.5, 0.5), // Target
#ifdef GEO_OS_EMSCRIPTEN
		512, 512, // OpenGL ES needs power of two textures
#else
		800, 800,            // Image size		
#endif		
		50.0                 // zoom angle (in degrees)
	    )
	{
            texture_ = 0;
	    scene_.add_object(new HorizontalCheckerboardPlane(0.0))  // The tradition !
   	          ->rename("Checkerboard");
	    
	    scene_.add_object(                                       // A sphere
		new Sphere(vec3(-0.7, -0.7, 1.0),0.7)
		)->set_reflection_coefficient(vec3(0.8, 0.8, 0.8))
		 ->set_diffuse_coefficient(vec3(0.2, 0.2, 0.2))
		 ->rename("Sphere 1");


	    scene_.add_object(                                       // Another sphere
		new Sphere(vec3( 0.0, 0.0, 0.0),0.25)
		)->set_diffuse_coefficient(vec3(0.0, 1.0, 0.5))
		 ->rename("Sphere 2");		

	    
	    // Let there be (two) lights !
	    scene_.add_object(new Light(
				 vec3(0.5, 1.3, 0.5), // Position
				 0.02,                // Radius
				 vec3(0.0, 0.5, 1.0)  // Color
				 )
	    )->rename("Light 1");
	    
	    scene_.add_object(new Light(
				 vec3(1.0, 0.2, 1.0), // Position
				 0.02,                // Radius
				 vec3(1.5, 1.5, 0.5)  // Color
			         )
	    )->rename("Light 2");
	    
	    scene_changed_ = true;
	    zoom_ = 0.0f;
	    left_pane_visible_ = false;
	    background_color_1_ = vec4f(0.0f, 0.0f, 0.2f, 1.0f);	
	    background_color_2_ = vec4f(0.0f, 0.0f, 0.2f, 1.0f);

	    for(index_t i=0; i<4; ++i) {
		quaternion_[i] = 0.0f;
	    }

	    for(index_t i=0; i<3; ++i) {
		translation_[i] = 0.0f;
	    }
	    
        }

        /**
         * \brief DemoRaytraceApplication destructor.
         */
        ~DemoRaytraceApplication() override {
            if(texture_ != 0) {
                glDeleteTextures(1,&texture_);
            }
        }

	void update() {

	    float* new_quaternion = glup_viewer_get_scene_quaternion();
	    for(index_t i=0; i<4; ++i) {
		if(new_quaternion[i] != quaternion_[i]) {
		    scene_changed_ = true;
		    quaternion_[i] = new_quaternion[i];		    
		}
	    }

	    float* new_translation = glup_viewer_get_scene_translation();
	    for(index_t i=0; i<3; ++i) {
		if(new_translation[i] != translation_[i]) {
		    scene_changed_ = true;
		    translation_[i] = new_translation[i];		    
		}
	    }

	    float new_zoom = glup_viewer_get_float(GLUP_VIEWER_ZOOM);
	    if(new_zoom != zoom_) {
		scene_changed_ = true;
		zoom_ = new_zoom;
	    }
	    
	    if(scene_changed_) {
		vec3 pos(2.0, 2.0, 1.5);
		mat4 M;
		set_mat4_from_translation_and_quaternion(
		    M,
		    double(translation_[0]),
		    double(translation_[1]),
		    double(translation_[2]),		    
		    double(quaternion_[2]),
		    double(quaternion_[0]),
		    double(quaternion_[1]),
		    double(quaternion_[3])
		);
		pos = transform_point(pos,M);
		
		vec3 target(
		    0.5 + double(translation_[0]),
		    0.5 + double(translation_[1]),
		    0.5 + double(translation_[2])		    
		);
		
		camera_.update(
		    pos, target, 50.0 / double(zoom_)
		);
#ifdef GEO_OPENMP	
#pragma omp parallel for
#endif
		for(index_t Y=0; Y<camera_.image_height(); ++Y) {
		    for(index_t X=0; X<camera_.image_width(); ++X) {
			Ray R = camera_.launch_ray(X,Y);
			vec3 K = scene_.raytrace(R);
			camera_.set_pixel(X,Y,K);
		    }
		}
	    }
	    if(scene_changed_ && texture_ != 0) {
		glTexImage2D(
		    GL_TEXTURE_2D,
		    0,
		    GL_RGB,
		    GLsizei(camera_.image_width()),
		    GLsizei(camera_.image_height()),
		    0,
		    GL_RGB,
		    GL_UNSIGNED_BYTE,
		    camera_.image_data()
		);
	    }
	    scene_changed_ = false;
	}

	
        /**
         * \brief Displays and handles the GUI for object properties.
         * \details Overloads Application::draw_object_properties().
         */
        void draw_object_properties() override {
	    if(ImGui::Button("Home",ImVec2(-1.0f, 0.0f))) {
		glup_viewer_home();
		zoom_ = 40.0;
		scene_changed_ = true;
	    }
	    if(scene_.draw_gui()) {
		scene_changed_ = true;
	    }
	    update();
        }

        /**
         * \brief Draws the scene according to currently set primitive and
         *  drawing modes.
         */
         void draw_scene() override {
	    update();
	    glupMatrixMode(GLUP_PROJECTION_MATRIX);
	    glupLoadIdentity();	    
	    glupMatrixMode(GLUP_MODELVIEW_MATRIX);
	    glupLoadIdentity();
	    int w,h;
	    glup_viewer_get_screen_size(&w,&h);
	    
	    glupScalef(2.0f*float(h)/float(w), 2.0f, 1.0f);	    	    
	    glupTranslatef(-0.5f, -0.5f, 0.0f);
	    glupEnable(GLUP_TEXTURING);
	    glupDisable(GLUP_LIGHTING);
	    glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
	    glBindTexture(GL_TEXTURE_2D, texture_);
	    glupBegin(GLUP_QUADS);
	    glupTexCoord2d(0.0, 0.0);
	    glupVertex2d(0.0, 1.0);
	    glupTexCoord2d(1.0, 0.0);
	    glupVertex2d(1.0, 1.0);	    
	    glupTexCoord2d(1.0, 1.0);
	    glupVertex2d(1.0, 0.0);
	    glupTexCoord2d(0.0, 1.0);
	    glupVertex2d(0.0, 0.0);	    
	    glupEnd();
	    glupDisable(GLUP_TEXTURING);
        }

        /**
         * \brief Creates the texture.
         * \details This function overloads Application::init_graphics(). It
         *  is called as soon as the OpenGL context is ready for rendering. It
         *  is meant to initialize the graphic objects used by the application.
         */
        void init_graphics() override {
            Application::init_graphics();
	    
            // Create the texture and initialize its texturing modes
            glGenTextures(1, &texture_);
            glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
            glBindTexture(GL_TEXTURE_2D, texture_);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glupTextureType(GLUP_TEXTURE_2D);
            glupTextureMode(GLUP_TEXTURE_REPLACE);
        }

	bool load(const std::string& filename) override {
	    mesh_load(filename, mesh_);
	    normalize_mesh(mesh_);
	    scene_.add_object(new MeshObject(mesh_));
	    scene_changed_ = true;
	    return true;
	}

	
    private:
	Camera camera_;
	Scene scene_;
	Mesh mesh_;
        GLuint texture_;
	bool scene_changed_;
	float translation_[3];
	float quaternion_[4];
	float zoom_;
    };
      
}

int main(int argc, char** argv) {
    DemoRaytraceApplication app(argc, argv, "<filename>");
    app.start();
    return 0;
}
