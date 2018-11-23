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

#ifndef GEOGRAM_GFX_IMGUI_EXT
#define GEOGRAM_GFX_IMGUI_EXT

#include <geogram_gfx/basic/common.h>

#ifdef GEO_COMPILER_CLANG
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#include <geogram_gfx/third_party/ImGui/imgui.h>

#ifdef GEO_COMPILER_CLANG
#pragma GCC diagnostic pop
#endif


/**
 * \file geogram_gfx/ImGui_ext/imgui_ext.h
 * \brief Extension functions for ImGui.
 */

namespace GEO {
    /**
     * \brief Maximum string length for ImGUI.
     * \details ImGUI uses plain old C strings (I'd prefer it to use
     *  std::string, but it is life...).
     */
    enum { geo_imgui_string_length = 4096 };
}

typedef int ImGuiExtFileDialogFlags;


enum ImGuiExtFileDialogFlags_ {
    ImGuiExtFileDialogFlags_Load  = 1,
    ImGuiExtFileDialogFlags_Save  = 2
};


namespace ImGui {

    /**
     * \brief Gets the global application scaling.
     * \details The global application scaling can be changed
     *   for high DPI displays. Some GUI elements need to
     *   be scaled accordingly.
     * \return The global application scaling.
     */
    float GEOGRAM_GFX_API scaling();

    /**
     * \brief Sets the global application scaling.
     * \details The global application scaling can be changed
     *   for high DPI displays. Some GUI elements need to
     *   be scaled accordingly.
     * \param x The global application scaling.
     */
    void GEOGRAM_GFX_API set_scaling(float x);
    
    /**
     * \brief Manages the GUI of a color editor.
     * \details This creates a custom dialog with the color editor and
     *  a default palette, as in ImGUI example.
     * \param[in] label the label of the widget, passed to ImGUI
     * \param[in,out] color a pointer to an array of 3 floats
     * \retval true if the color was changed
     * \retval false otherwise
     */
    bool GEOGRAM_GFX_API ColorEdit3WithPalette(
	const char* label, float* color
    );

    /**
     * \brief Manages the GUI of a color editor.
     * \details This creates a custom dialog with the color editor and
     *  a default palette, as in ImGUI example.
     * \param[in] label the label of the widget, passed to ImGUI
     * \param[in,out] color a pointer to an array of 4 floats
     * \retval true if the color was changed
     * \retval false otherwise
     */
    bool GEOGRAM_GFX_API ColorEdit4WithPalette(
	const char* label, float* color
    );

    // extensions: ';'-separated list of extensions, whitout '.'
    void GEOGRAM_GFX_API OpenFileDialog(
	const char* label,
	const char* extensions,
	const char* filename,
	ImGuiExtFileDialogFlags flags
    );

    bool GEOGRAM_GFX_API FileDialog(
	const char* label,
	char* filename, size_t filename_buff_len
    );
    
}

#endif

