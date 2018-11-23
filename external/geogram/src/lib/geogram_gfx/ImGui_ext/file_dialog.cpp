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

#include <geogram_gfx/ImGui_ext/file_dialog.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/string.h>

#include <algorithm>

namespace {
    
    /**
     * \brief Safer version of strncpy()
     * \param[in] dest a pointer to the destination string
     * \param[in] source a pointer to the source string
     * \param[in] max_dest_size number of characters available in
     *  destination string
     * \return the length of the destination string after copy. If
     *  the source string + null terminator was greater than max_dest_size,
     *  then it is cropped. On exit, dest is always null-terminated (in 
     *  contrast with strncpy()).
     */ 
    size_t safe_strncpy(
        char* dest, const char* source, size_t max_dest_size
    ) {
        strncpy(dest, source, max_dest_size-1);
        dest[max_dest_size-1] = '\0';
        return strlen(dest);
    }

   /**
    * \brief Converts a complete path to a file to a label
    *  displayed in the file browser.
    * \details Strips viewer_path from the input path.
    * \param[in] path the complete path, can be either a directory or
    *  a file
    * \return the label to be displayed in the menu
    */
    std::string path_to_label(
        const std::string& viewer_path, const std::string& path
    ) {
        std::string result = path;
        if(GEO::String::string_starts_with(result, viewer_path)) {
            result = result.substr(
                viewer_path.length(), result.length()-viewer_path.length()
            );
        }
        return result;
    }
}

namespace GEO {

    FileDialog::FileDialog(
        Application* app, bool save_mode, const std::string& default_filename
    ) :
        application_(app),
        visible_(false),
        current_write_extension_index_(0),        
        pinned_(false),
        show_hidden_(false),
        scroll_to_file_(false),
        save_mode_(save_mode),
        are_you_sure_(false)
    {
#if defined(GEO_OS_WINDOWS) || defined(GEO_OS_ANDROID)
	directory_ = FileSystem::documents_directory() + "/";
#else	
        directory_ = FileSystem::get_current_working_directory() + "/";
#endif	
	set_default_filename(default_filename);
	current_file_index_ = 0;
	current_directory_index_ = 0;
	current_write_extension_index_ = 0;
    }

    FileDialog::FileDialog() :
        application_(nullptr),
        visible_(false),
        current_write_extension_index_(0),        
        pinned_(false),
        show_hidden_(false),
        scroll_to_file_(false),
        save_mode_(false),
        are_you_sure_(false)
    {
#if defined(GEO_OS_WINDOWS) || defined(GEO_OS_ANDROID)
	directory_ = FileSystem::documents_directory() + "/";
#else	
        directory_ = FileSystem::get_current_working_directory() + "/";
#endif	
	current_file_index_ = 0;
	current_directory_index_ = 0;
	current_write_extension_index_ = 0;
    }
    
    void FileDialog::set_default_filename(const std::string& default_filename) {
        safe_strncpy(
            current_file_, default_filename.c_str(), sizeof(current_file_)
        );
    }
    
    void FileDialog::update_files() {
        directories_.clear();
        files_.clear();

        directories_.push_back("../");
        
        std::vector<std::string> entries;
        FileSystem::get_directory_entries(directory_, entries);
        std::sort(entries.begin(), entries.end());
        for(index_t i=0; i<entries.size(); ++i) {
            if(can_load(entries[i])) {
                files_.push_back(path_to_label(directory_,entries[i]));
            } else if(FileSystem::is_directory(entries[i])) {
                std::string subdir =
                    path_to_label(directory_,entries[i]) + "/";
                if(show_hidden_ || subdir[0] != '.') {
                    directories_.push_back(subdir);
                }
            }
        }
        if(current_directory_index_ >= directories_.size()) {
            current_directory_index_ = 0;
        }
        if(current_file_index_ >= files_.size()) {
            current_file_index_ = 0;
        }
        if(!save_mode_) {
            if(current_file_index_ >= files_.size()) {
                current_file_[0] = '\0';
            } else {
                safe_strncpy(
                    current_file_,
                    files_[current_file_index_].c_str(),
                    sizeof(current_file_)
                );
            }
        }
    }

    void FileDialog::set_directory(const std::string& directory) {
        current_directory_index_ = 0;
        current_file_index_ = 0;
        if(directory[0] == '/' || directory[1] == ':') {
            directory_ = directory;
        } else {
            directory_ = FileSystem::normalized_path(
                directory_ + "/" +
                directory 
            );
        }
        if(directory_[directory_.length()-1] != '/') {
            directory_ += "/";
        }
        update_files();
    }

    int FileDialog::text_input_callback(ImGuiTextEditCallbackData* data) {
        FileDialog* dlg = static_cast<FileDialog*>(data->UserData);
        if((data->EventFlag & ImGuiInputTextFlags_CallbackCompletion) != 0) {
            dlg->tab_callback(data);
        } else if(
            (data->EventFlag & ImGuiInputTextFlags_CallbackHistory) != 0
        ) {
            if(data->EventKey == ImGuiKey_UpArrow) {
                dlg->updown_callback(data,-1);
            } else if(data->EventKey == ImGuiKey_DownArrow) {
                dlg->updown_callback(data,1);                
            }
        } 
        return 0;
    }


    void FileDialog::updown_callback(
        ImGuiTextEditCallbackData* data, int direction
    ) {
        int next = int(current_file_index_) + direction;
        if(next < 0) {
            if(files_.size() == 0) {
                current_file_index_ = 0;
            } else {
                current_file_index_ = index_t(files_.size()-1);
            }
        } else if(next >= int(files_.size())) {
            current_file_index_ = 0;
        } else {
            current_file_index_ = index_t(next);
        }

        if(files_.size() == 0) {
            current_file_[0] = '\0';
        } else {
            safe_strncpy(
                current_file_,
                files_[current_file_index_].c_str(),
                sizeof(current_file_)
            );
        }
        update_text_edit_callback_data(data);
        scroll_to_file_ = true;        
    }

    void FileDialog::update_text_edit_callback_data(
        ImGuiTextEditCallbackData* data
    ) {
        data->BufTextLen = int(
            safe_strncpy(
                data->Buf, current_file_, (size_t)data->BufSize
            )
        );
        data->CursorPos = data->BufTextLen;
        data->SelectionStart = data->BufTextLen;
        data->SelectionEnd = data->BufTextLen;
        data->BufDirty = true;
    }
    
    void FileDialog::tab_callback(ImGuiTextEditCallbackData* data) {
        std::string file(current_file_);
        bool found = false;
        for(index_t i=0; i<files_.size(); ++i) {
            if(String::string_starts_with(files_[i],file)) {
                current_file_index_ = i;
                found = true;
                break;
            }
        }
        if(found) {
            safe_strncpy(
                current_file_,
                files_[current_file_index_].c_str(),
                sizeof(current_file_)
            );
            update_text_edit_callback_data(data);
            scroll_to_file_ = true;
        }
    }

    void FileDialog::file_selected(bool force) {
        std::string file =
            FileSystem::normalized_path(
                directory_ + "/" +
                current_file_
            );
        
        if(save_mode_) {
            if(!force && FileSystem::is_file(file)) {
                are_you_sure_ = true;
                return;
            } else {
		selected_file_ = file;
                if(
		    application_ != nullptr &&
		    application_->save(file)
		) {
		    Logger::out("I/O") << "Saved "
				       << current_file_ << std::endl;
		}
            }
        } else {
	    selected_file_ = file;
	    if(application_ != nullptr) {
		application_->load(file);
	    }
        }
        
        if(!pinned_) {
            hide();
        }
    }
    
    void FileDialog::draw() {
        
        if(!visible_) {
            return;
        }

        ImGui::SetNextWindowSize(
            ImVec2(
                ImGui::scaling()*400.0f,
                ImGui::scaling()*400.0f
            ),
            ImGuiCond_Once
        );

        ImGui::Begin(
            (std::string(
                save_mode_ ? "Save as...##" : "Load...##"
             )+String::to_string(this)).c_str(),
            &visible_,
            ImGuiWindowFlags_NoCollapse 
        );
        
        if(ImGui::Button("parent")) {
            set_directory("../");
        }
        ImGui::SameLine();
        if(ImGui::Button("home")) {
            set_directory(FileSystem::documents_directory());
            update_files();
        }
        ImGui::SameLine();            
        if(ImGui::Button("refresh")) {
            update_files();
        }

	if(!save_mode_) {
	    ImGui::SameLine();        
	    ImGui::Text("pin");
	    ImGui::SameLine();
	    ImGui::Checkbox("##pin", &pinned_);
	    if(ImGui::IsItemHovered()) {
		ImGui::SetTooltip("Keeps this dialog open.");
	    }
	}

	draw_disk_drives();
        ImGui::Separator();
	
        {
            std::vector<std::string> path;
            String::split_string(directory_, '/', path);
            for(index_t i=0; i<path.size(); ++i) {
                if(i != 0) {
                    ImGui::SameLine();
                }
                // We need to generate a unique id, else there is an id
                // clash with the "home" button right before !!
                if(ImGui::SmallButton(
                       (path[i] + "##path" + String::to_string(i)).c_str())
                ) {
                    std::string new_dir;
                    if(path[0].length() >= 2 && path[0][1] == ':') {
                        new_dir = path[0];
                    } else {
                        new_dir += "/";
                        new_dir += path[0];
                    }
                    for(index_t j=1; j<=i; ++j) {
                        new_dir += "/";
                        new_dir += path[j];
                    }
                    set_directory(new_dir);
                }
            }
        }

        const float footer_size = 35.0f*ImGui::scaling();
        {
            ImGui::BeginChild(
                "##directories",
                ImVec2(ImGui::GetWindowWidth()*0.5f-10.0f*ImGui::scaling(), -footer_size),
                true
            );
            for(index_t i=0; i<directories_.size(); ++i) {
                if(ImGui::Selectable(
                       directories_[i].c_str(),
                       (i == current_directory_index_)
                       )
                ) {
                    current_directory_index_ = i;
                    set_directory(directories_[current_directory_index_]);
                }
            }
            ImGui::EndChild();
        }
        ImGui::SameLine();
        {
            ImGui::BeginChild(
                "##files",
                ImVec2(ImGui::GetWindowWidth()*0.5f-10.0f*ImGui::scaling(), -footer_size),
                true
            );
            for(index_t i=0; i<files_.size(); ++i) {
                if(ImGui::Selectable(
                       files_[i].c_str(),
                       (i == current_file_index_)
                   )
                ) {
                    safe_strncpy(
                        current_file_,files_[i].c_str(),sizeof(current_file_)
                    );
                    if(current_file_index_ == i) {
                        file_selected();
                    } else {
                        current_file_index_ = i;
                    }
                }
                if(scroll_to_file_ && i == current_file_index_) {
                    ImGui::SetScrollHere();
                    scroll_to_file_ = false;
                }
            }
            ImGui::EndChild();

            {
                if(ImGui::Button(save_mode_ ? "Save as" : "Load")) {
                    file_selected();
                }
                ImGui::SameLine();
                ImGui::PushItemWidth(
                    save_mode_ ? -80.0f*ImGui::scaling() : -5.0f*ImGui::scaling()
                );
                if(ImGui::InputText(
                       "##filename",
                       current_file_, geo_imgui_string_length,
                       ImGuiInputTextFlags_EnterReturnsTrue    |
                       ImGuiInputTextFlags_CallbackHistory     |
                       ImGuiInputTextFlags_CallbackCompletion ,
                       text_input_callback,
                       this 
                       )
                ) {
                    scroll_to_file_ = true;
                    std::string file = current_file_;
                    for(index_t i=0; i<files_.size(); ++i) {
                        if(files_[i] == file) {
                            current_file_index_ = i;
                        }
                    }
                    file_selected();
                }
                ImGui::PopItemWidth();
                // Keep auto focus on the input box
                if (ImGui::IsItemHovered()) {
                    // Auto focus previous widget                
                    ImGui::SetKeyboardFocusHere(-1); 
                }

                if(save_mode_) {
                    ImGui::SameLine();
                    ImGui::PushItemWidth(-5.0f*ImGui::scaling());

                    if(
		       application_ != nullptr &&
		       extensions_.size() == 0
		    ) {
                        String::split_string(
                            application_->supported_write_file_extensions(),
                            ';', extensions_
                        );
                        std::string ext =
                            FileSystem::extension(current_file_);
                        for(index_t i=0; i<extensions_.size(); ++i) {
                            if(extensions_[i] == ext) {
                                current_write_extension_index_ = i;
                            }
                        }
                    }
		    
                    std::vector<const char*> write_extensions;
                    for(index_t i=0; i<extensions_.size(); ++i) {
                        write_extensions.push_back(&extensions_[i][0]);
                    }
                    if(ImGui::Combo(
                           "##extension",
                           (int*)(&current_write_extension_index_),
                           &write_extensions[0],
                           int(write_extensions.size())
                       )
                    ) {
                        std::string file = current_file_;
                        file = FileSystem::base_name(file) + "." +
                            extensions_[current_write_extension_index_];
                        safe_strncpy(
                            current_file_, file.c_str(), sizeof(current_file_)
                        );
                    }
                    ImGui::PopItemWidth();
                }
            }
        }
        ImGui::End();
        draw_are_you_sure();
    }

    void FileDialog::draw_are_you_sure() {
        if(are_you_sure_) {
            ImGui::OpenPopup("File exists");
        }
        if(
            ImGui::BeginPopupModal(
                "File exists", nullptr, ImGuiWindowFlags_AlwaysAutoResize
            )
        ) {
            ImGui::Text(
                "%s",
                 (std::string("File ") + current_file_ +
                  " already exists\nDo you want to overwrite it ?"
                 ).c_str()
            );
            ImGui::Separator();
            float w = 120.0f*ImGui::scaling();
            if (ImGui::Button("Overwrite", ImVec2(w,0))) {
                are_you_sure_ = false;
                ImGui::CloseCurrentPopup();
                file_selected(true);
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel",ImVec2(w,0))) {
                are_you_sure_ = false;                
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }
    }

    bool FileDialog::can_load(const std::string& filename) {
	if(!FileSystem::is_file(filename)) {
	    return false;
	}
	if(application_ != nullptr) {
	    return application_->can_load(filename);
	}
	std::string ext = FileSystem::extension(filename);
	for(size_t i=0; i<extensions_.size(); ++i) {
	    if(extensions_[i] == ext || extensions_[i] == "*") {
		return true;
	    }
	}
	return false;
    }


    void FileDialog::set_extensions(const std::string& extensions) {
	extensions_.clear();
	GEO::String::split_string(extensions, ';', extensions_);
    }

    void FileDialog::draw_disk_drives() {
#ifdef GEO_OS_WINDOWS	
	DWORD drives = GetLogicalDrives();
	for(DWORD b=0; b<16; ++b) {
	    if((drives & (1u << b)) != 0) {
		std::string drive;
		drive += char('A' + char(b));
		drive += ":";
		if(ImGui::Button(drive.c_str())) {
		    set_directory(drive);
		}
		ImGui::SameLine();
	    }
	}
#endif	
    }
    
}

