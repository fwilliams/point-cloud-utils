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

#ifndef GEOGRAM_GFX_IMGUI_EXT_FILE_DIALOG
#define GEOGRAM_GFX_IMGUI_EXT_FILE_DIALOG

#include <geogram_gfx/basic/common.h>
#include <geogram_gfx/ImGui_ext/imgui_ext.h>

/**
 * \file geogram_gfx/ImGui_ext/file_dialog.h
 * \brief Implementation of a file dialog for ImGui.
 */

namespace GEO {

    // TODO: make it independent of GEO::Application
    // through callbacks.
    
    class Application;
    
    /**
     * \brief Implementation of GLUP viewer's file dialog.
     */
    class GEOGRAM_GFX_API FileDialog {
    public:

        /**
         * \brief FileDialog constructor.
         * \param[in] application a pointer to the Application
         * \param[in] save_mode if true, FileDialog is used to create files
         * \param[in] default_filename the default file name used if save_mode
         *  is set
         */
        FileDialog(
            Application* application,
            bool save_mode=false,
            const std::string& default_filename=""
        );

	/**
	 * \brief FileDialog constructor.
	 */
	FileDialog();
	
	/** 
	 * \brief Sets the default file.
	 * \details Only valid if save_mode is set.
         * \param[in] default_filename the default file name.
	 */	
	void set_default_filename(const std::string& default_filename);
	
        /**
         * \brief Makes this FileDialog visible.
         */
        void show() {
            update_files();
            visible_ = true;
        }

        /**
         * \brief Makes this FileDialog invisibile.
         */
        void hide() {
            visible_ = false;
        }

        /**
         * \brief Tests whether this FileDialog is visible.
         * \retval true if this FileDialog is visible
         * \retval false otherwise
         */
        bool is_visible() const {
            return visible_;
        }

        /**
         * \brief Draws the console and handles the gui.
         */
        void draw();

	/**
	 * \brief Sets whether this file dialog is for 
	 *  saving file.
	 * \details If this file dialog is for saving file,
	 *  then the user can enter the name of a non-existing
	 *  file, else he can only select existing files.
	 * \param[in] x true if this file dialog is for
	 *  saving file.
	 */
	void set_save_mode(bool x) {
	    save_mode_ = x;
	}

	/**
	 * \brief Gets the selected file if any and resets it
	 *  to the empty string.
	 * \return the selected file if there is any or the 
	 *  empty string otherwise.
	 */
	std::string get_and_reset_selected_file() {
	    std::string result;
	    std::swap(result,selected_file_);
	    return result;
	}

	/**
	 * \brief Defines the file extensions managed by this 
	 *  FileDialog.
	 * \param[in] extensions a ';'-separated list of extensions
	 */
	void set_extensions(const std::string& extensions);
	
    protected:

	/**
	 * \brief Tests whether a file can be read.
	 * \param[in] filename the file name to be tested.
	 * \retval true if this file can be read.
	 * \retval false otherwise.
	 */
	bool can_load(const std::string& filename);
	
        /**
         * \brief Updates the list of files and directories
         *  displayed by this FileDialog.
         */
        void update_files();

        /**
         * \brief Changes the current directory.
         * \param[in] directory either the path relative to the
         *  current directory or an absolute path
         */
        void set_directory(const std::string& directory);

        /**
         * \brief The callback for handling the text input.
         * \param[in,out] data a pointer to the callback data
         */
        static int text_input_callback(ImGuiTextEditCallbackData* data);

        /**
         * \brief Called whenever the up or down arrows are pressed.
         * \param[in,out] data a pointer to the callback data
         * \param[in] direction -1 if the up arrow was pressed, 1 if the
         *  down arrow was pressed
         */
        void updown_callback(ImGuiTextEditCallbackData* data, int direction);

        /**
         * \brief Called whenever the tab key is pressed.
         * \param[in,out] data a pointer to the callback data
         */
        void tab_callback(ImGuiTextEditCallbackData* data);

        /**
         * \brief Copies the currently selected file into the 
         *  string currently manipulated by InputText.
         * \param[out] data a pointer to the callback data
         */
        void update_text_edit_callback_data(
            ImGuiTextEditCallbackData* data
        );
        
        /**
         * \brief Called whenever a file is selected.
         * \param[in] force in save_mode, if set, 
         *  overwrites the file even if it already 
         *  exists.
         */
        void file_selected(bool force=false);

	/**
	 * \brief Handles the "are you sure ?" dialog
	 *  when a file is about to be overwritten.
	 */
        void draw_are_you_sure();

	/**
	 * \brief Under Windows, add buttons to change
	 *  disk drive.
	 */
	void draw_disk_drives();
	
    private:
        Application* application_;
        bool visible_;
        std::string directory_;
        index_t current_directory_index_;
        index_t current_file_index_;
        std::vector<std::string> directories_;
        std::vector<std::string> files_;
        std::vector<std::string> extensions_;
        index_t current_write_extension_index_;
        char current_file_[geo_imgui_string_length];
        bool pinned_;
        bool show_hidden_;
        bool scroll_to_file_;
        bool save_mode_;
        bool are_you_sure_;
	std::string selected_file_;
    };
    
}

#endif
