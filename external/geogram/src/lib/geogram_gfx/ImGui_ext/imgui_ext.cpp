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

#include <geogram_gfx/ImGui_ext/imgui_ext.h>
#include <geogram_gfx/ImGui_ext/file_dialog.h>
#include <geogram_gfx/third_party/ImGui/imgui.h>
#include <geogram_gfx/third_party/ImGui/imgui_internal.h>
#include <geogram/basic/logger.h>
#include <map>

namespace {
    bool initialized = false;

    std::map<std::string, GEO::FileDialog*> file_dialogs;
    
    void terminate() {
	for(auto& it : file_dialogs) {
	    delete it.second;
	}
    }
    
    void initialize() {
	if(!initialized) {
	    initialized = true;
	    atexit(terminate);
	}
    }

    /**
     * \brief Manages the GUI of a color editor.
     * \details This creates a custom dialog with the color editor and
     *  a default palette, as in ImGUI example.
     * \param[in] label the label of the widget, passed to ImGUI
     * \param[in,out] color_in a pointer to an array of 3 floats if 
     *  with_alpha is false or 4 floats if with_alpha is true
     * \param[in] with_alpha true if transparency is edited, false otherwise
     * \retval true if the color was changed
     * \retval false otherwise
     */
    bool ColorEdit3or4WithPalette(
	const char* label, float* color_in, bool with_alpha
    ) {
	bool result = false;
	static bool saved_palette_initialized = false;
	static ImVec4 saved_palette[40];
	static ImVec4 backup_color;
	ImGui::PushID(label);
	int flags =
	    ImGuiColorEditFlags_PickerHueWheel |
	    ImGuiColorEditFlags_Float;

	if(!with_alpha) {
	    flags |= ImGuiColorEditFlags_NoAlpha ; 
	}
	
	ImVec4& color = *(ImVec4*)color_in;

	if (!saved_palette_initialized) {

	    for (int n = 0; n < 8; n++) {
		saved_palette[n].x = 0.0f;
		saved_palette[n].y = 0.0f;
		saved_palette[n].z = 0.0f;
	    }

	    saved_palette[0] = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
	    saved_palette[1] = ImVec4(0.5f, 0.5f, 0.5f, 1.0f);
	    saved_palette[2] = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
	    saved_palette[3] = ImVec4(1.0f, 0.0f, 0.0f, 1.0f);
	    saved_palette[4] = ImVec4(1.0f, 1.0f, 0.0f, 1.0f);
	    saved_palette[5] = ImVec4(0.0f, 0.0f, 1.0f, 1.0f);
	    saved_palette[6] = ImVec4(0.0f, 1.0f, 0.0f, 1.0f);
	    saved_palette[7] = ImVec4(0.0f, 1.0f, 1.0f, 1.0f);	    
	    
	    for (int n = 0; n < 32; n++) {
		ImGui::ColorConvertHSVtoRGB(
		    float(n) / 31.0f, 0.8f, 0.8f,
		    saved_palette[n+8].x,
		    saved_palette[n+8].y,
		    saved_palette[n+8].z
		);
	    }
	    saved_palette_initialized = true;
	}
	
	bool open_popup = ImGui::ColorButton(label, color, flags);
	
	if(label[0] != '#') {
	    ImGui::SameLine();	    
	    ImGui::Text("%s",label);
	}
	if (open_popup) {
	    ImGui::OpenPopup("##PickerPopup");
	    backup_color = color;
	}
	if (ImGui::BeginPopup("##PickerPopup")) {
	    if(label[0] != '#') {
		ImGui::Text("%s",label);
		ImGui::Separator();
	    }
	    if(ImGui::ColorPicker4(
		"##picker", (float*)&color,
		flags | ImGuiColorEditFlags_NoSidePreview
		      | ImGuiColorEditFlags_NoSmallPreview
		   )
	    ) {
		result = true;
	    }
	    ImGui::SameLine();
	    ImGui::BeginGroup();
	    ImGui::Text("Current");
	    ImGui::ColorButton(
		"##current", color,
		ImGuiColorEditFlags_NoPicker |
		ImGuiColorEditFlags_AlphaPreviewHalf,
		ImVec2(60,40)
	    );
	    ImGui::Text("Previous");
	    if (ImGui::ColorButton(
		    "##previous", backup_color,
		    ImGuiColorEditFlags_NoPicker |
		    ImGuiColorEditFlags_AlphaPreviewHalf,
		    ImVec2(60,40))
		) {
		color = backup_color;
		result = true;
	    }
	    ImGui::Separator();
	    ImGui::Text("Palette");
	    for (int n = 0; n < 40; n++) {
		ImGui::PushID(n);
		if ( (n % 6) != 0 ) {
		    ImGui::SameLine(0.0f, ImGui::GetStyle().ItemSpacing.y);
		}
		if (ImGui::ColorButton(
			"##palette",
			saved_palette[n],
			ImGuiColorEditFlags_NoPicker |
			ImGuiColorEditFlags_NoTooltip,
			ImVec2(20,20))
		    ) {
		    color = ImVec4(
			saved_palette[n].x,
			saved_palette[n].y,
			saved_palette[n].z,
			color.w
		    ); // Preserve alpha!
		    result = true;
		}
		ImGui::PopID();
	    }
	    ImGui::Separator();
	    if(ImGui::Button(
		   "OK", ImVec2(-1, -1)
	       )
	    ) {
		ImGui::CloseCurrentPopup();
	    }
	    ImGui::EndGroup();
	    ImGui::EndPopup();
	}
	ImGui::PopID();
	return result;
    }
}

namespace ImGui {

    float scaling() {
	ImGuiContext& g = *GImGui;
	float s = 1.0;
	if(g.Font->FontSize > 40.0f) {
	    s = g.Font->FontSize / 30.0f;
	} else {
	    s = g.Font->FontSize / 20.0f;	    
	}
	return s * ImGui::GetIO().FontGlobalScale;
    }

    void set_scaling(float x) {
	ImGui::GetIO().FontGlobalScale = x;
    }
    
    /*******************************************************************/
    
    bool ColorEdit3WithPalette(const char* label, float* color_in) {
	return ColorEdit3or4WithPalette(label, color_in, false);
    }

    bool ColorEdit4WithPalette(const char* label, float* color_in) {
	return ColorEdit3or4WithPalette(label, color_in, true);	
    }

    /*******************************************************************/
    
    void OpenFileDialog(
	const char* label,
	const char* extensions,
	const char* filename,
	ImGuiExtFileDialogFlags flags
    ) {
	initialize();	
	GEO::FileDialog* dlg = nullptr;
	if(file_dialogs.find(label) == file_dialogs.end()) {
	    file_dialogs[label] = new GEO::FileDialog();
	}
	dlg = file_dialogs[label];
	dlg->set_extensions(extensions); 
	if(flags == ImGuiExtFileDialogFlags_Save) {
	    dlg->set_save_mode(true);
	    dlg->set_default_filename(filename);
	} else {
	    dlg->set_save_mode(false);
	}
	dlg->show();
    }

    bool FileDialog(
	const char* label, char* filename, size_t filename_buff_len
    ) {
	if(file_dialogs.find(label) == file_dialogs.end()) {
	    filename[0] = '\0';
	    return false;
	}
	GEO::FileDialog* dlg = file_dialogs[label];
	dlg->draw();
	std::string result = dlg->get_and_reset_selected_file();
	if(result != "") {
	    if(result.length() + 1 >= filename_buff_len) {
		GEO::Logger::err("ImGui_ext") << "filename_buff_len exceeded"
					      << std::endl;
		return false;
	    } else {
		strcpy(filename, result.c_str());
		return true;
	    }
	} else {
	    return false;
	}
    }

    /****************************************************************/
    
}

