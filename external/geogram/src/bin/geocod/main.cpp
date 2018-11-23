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

#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/lua/lua_glup.h>
#include <geogram/lua/lua_io.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>
#include <algorithm>

extern "C" {
#include <geogram/third_party/lua/lua.h>
#include <geogram/third_party/lua/lauxlib.h>
#include <geogram/third_party/lua/lualib.h>
}

extern void register_embedded_lua_files(void);

namespace {

    using namespace GEO;

    /**
     * \brief An application that demonstrates both
     *  GLUP primitives and glup_viewer application
     *  framework.
     */
    class GeoCodApplication : public Application {
    public:

        /**
         * \brief GeoCodApplication constructor.
         */
        GeoCodApplication(
            int argc, char** argv,
            const std::string& usage
        ) : Application(argc, argv, usage) {
	    save_dialog_.set_default_filename("hello.lua");	    
	    init_graphics_called_ = false;
	    text_editor_was_visible_ = false;
            // Define the 3d region that we want to display
            // (xmin, ymin, zmin, xmax, ymax, zmax)
            glup_viewer_set_region_of_interest(
                0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f
            );
	    register_embedded_lua_files();
	    exec_command("require(\"preamble\")");
        }

        /**
         * \brief GeoCodApplication destructor.
         */
        virtual ~GeoCodApplication() {
        }
        
        /**
         * \brief Displays and handles the GUI for object properties.
         * \details Overloads Application::draw_object_properties().
         */
        virtual void draw_object_properties() {
	    if(!lua_error_occured_) {
		exec_command("imgui.draw_object_properties()");
	    }
        }

        /**
         * \brief Draws the scene according to currently set primitive and
         *  drawing modes.
         */
        virtual void draw_scene() {
	    if(text_editor_was_visible_ && !text_editor_visible_) {
		run_program();
	    }
	    if(!lua_error_occured_) {
		exec_command("GLUP.draw_scene()");
	    }
	    text_editor_was_visible_ = text_editor_visible_;
        }

        /**
         * \brief Initializes graphics objects.
         * \details This function overloads Application::init_graphics(). It
         *  is called as soon as the OpenGL context is ready for rendering. It
         *  is meant to initialize the graphic objects used by the application.
         */
        virtual void init_graphics() {
            Application::init_graphics();
	    if(!lua_error_occured_) {
		exec_command("GLUP.init_graphics()");
	    }
	    init_graphics_called_ = true;
        }

	void embedded_files_menu(const std::string& prefix) {
	    std::vector<std::string> embedded_files;
	    list_embedded_lua_files(embedded_files);
	    std::sort(embedded_files.begin(), embedded_files.end());
	    for(index_t i=0; i<embedded_files.size(); ++i) {
		if(String::string_starts_with(
		       embedded_files[i],prefix) &&
		   ImGui::MenuItem(embedded_files[i].c_str())
		) {
		    const char* data;
		    get_embedded_lua_file(embedded_files[i].c_str(), &data);
		    text_editor_.load_data(data);
		    run_program();
		    current_file_ = "";		    			
		}
	    }
	}

	
	virtual void draw_fileops_menu() {
	    if(ImGui::BeginMenu("New...")) {
		if(ImGui::MenuItem("empty file")) {
		    text_editor_.clear();
		    current_file_ = "";
		    exec_command("require(\"preamble\")");
		}
		ImGui::Separator();
		ImGui::MenuItem("From template...", nullptr, false, false);
		embedded_files_menu("templates/");
		ImGui::EndMenu();
	    }
	    if(ImGui::BeginMenu("Load example...")) {
		embedded_files_menu("examples/");		
		ImGui::EndMenu();
	    }
	    if(ImGui::BeginMenu("Load game...")) {
		embedded_files_menu("games/");		
		ImGui::EndMenu();
	    }
	    if(ImGui::BeginMenu("Load Shift and Tab\'s adventure...")) {
		embedded_files_menu("book/");		
		ImGui::EndMenu();
	    }
	    if(ImGui::BeginMenu("Load internal lib...")) {
		ImGui::MenuItem("These files are those", nullptr, false, false);
		ImGui::MenuItem(
		    "that are read by require(\"...\")", nullptr, false, false
	        );
		ImGui::MenuItem("You cannot modify them",    nullptr, false, false);
		ImGui::MenuItem("(but you can take a look)", nullptr, false, false);
		ImGui::Separator();
		embedded_files_menu("lib/");		
		ImGui::EndMenu();
	    }
	    if(ImGui::MenuItem("Run program [F5]")) {
		run_program();
	    }
	}

	void run_program() {
	    exec_command("require(\"preamble\")");	    
	    if(exec_command(text_editor_.text().c_str())) {
		Logger::out("LUA") << "Program is OK." << std::endl;
	    } else {
		adjust_lua_glup_state(lua_state_);		
	    }
	    if(!lua_error_occured_ && init_graphics_called_) {
		exec_command("GLUP.init_graphics()");
	    }
	    glup_viewer_enable(GLUP_VIEWER_IDLE_REDRAW);
	}
	
        /**
         * \brief Draws the application menus.
         * \details This function overloads 
         *  Application::draw_application_menus(). It can be used to create
         *  additional menus in the main menu bar.
         */
        virtual void draw_application_menus() {
	    if(!lua_error_occured_) {
		exec_command("imgui.draw_application_menus()");
	    }
	}

	/**
	 * \copydoc Application::load()
	 */
	virtual bool load(const std::string& filename) {
	    text_editor_.load(filename);
	    run_program();
	    current_file_ = filename;
	    return true;
	}

	/**
	 * \copydoc Application::save()
	 */
        virtual bool save(const std::string& filename) {
	    text_editor_.save(filename);
	    current_file_ = filename;
	    return true;
	}

	/**
	 * \copydoc Application::supported_read_file_extensions()
	 */
	virtual std::string supported_read_file_extensions() {
	    return "lua";
	}

	/**
	 * \copydoc Application::supported_write_file_extensions()
	 */
	virtual std::string supported_write_file_extensions() {
	    return "lua";
	}

	/**
	 * \copydoc Application::on_key_pressed()
	 */
	virtual bool on_key_pressed(const char* c) {
	    if(!strcmp(c,"F5")) {
		run_program();
	    } else if(!strcmp(c,"F2") && current_file_ != "") {
		if(save(current_file_)) {
		    Logger::out("I/O")
			<< "Saved " << current_file_ << std::endl;
		}
	    } else {
		exec_command(
		    (
			std::string("imgui.on_key_pressed(\"") + c + "\")"
			).c_str()
		    );
	    }
	    return true;
	}

	/**
	 * \copydoc Application::on_key_released()
	 */
	virtual bool on_key_released(const char* c) {
	    exec_command(
		(
		    std::string("imgui.on_key_released(\"") + c + "\")"
		).c_str()
	    );
	    return true;
	}
	
    protected:
	virtual void draw_viewer_properties() {
	    if(ImGui::Button("run program [F5]", ImVec2(-1,0))) {
		run_program();
	    }
	    Application::draw_viewer_properties();
	}
        
    private:
	bool init_graphics_called_;
	bool text_editor_was_visible_;
    };
      
}

int main(int argc, char** argv) {
    GeoCodApplication app(argc, argv, "<filename>");
    app.start();
    return 0;
}
