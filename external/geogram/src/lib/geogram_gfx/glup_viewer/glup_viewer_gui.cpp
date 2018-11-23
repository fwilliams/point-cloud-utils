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

#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui_private.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/glup_viewer/geogram_logo_256.xpm>

#ifdef GEOGRAM_WITH_LUA
#include <geogram_gfx/lua/lua_glup.h>
#include <geogram_gfx/lua/lua_glup_viewer.h>
#include <geogram_gfx/lua/lua_imgui.h>
#include <geogram/lua/lua_io.h>
#endif

#include <geogram/mesh/mesh_io.h>

#include <geogram/basic/string.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/stopwatch.h>

#include <geogram/bibliography/bibliography.h>

#ifdef GEO_OS_EMSCRIPTEN
#include <emscripten.h>
#endif

#include <algorithm>
#include <stdio.h>

#include <geogram_gfx/glup_viewer/colormaps/french.xpm>
#include <geogram_gfx/glup_viewer/colormaps/black_white.xpm>
#include <geogram_gfx/glup_viewer/colormaps/viridis.xpm>
#include <geogram_gfx/glup_viewer/colormaps/rainbow.xpm>
#include <geogram_gfx/glup_viewer/colormaps/cei_60757.xpm>
#include <geogram_gfx/glup_viewer/colormaps/inferno.xpm>
#include <geogram_gfx/glup_viewer/colormaps/magma.xpm>
#include <geogram_gfx/glup_viewer/colormaps/parula.xpm>
#include <geogram_gfx/glup_viewer/colormaps/plasma.xpm>
#include <geogram_gfx/glup_viewer/colormaps/blue_red.xpm>

extern "C" {
    void glup_viewer_char_callback(void* w, unsigned int c);
}


#ifdef GEOGRAM_WITH_LUA
extern "C" {
#include <geogram/third_party/lua/lua.h>    
#include <geogram/third_party/lua/lauxlib.h>
#include <geogram/third_party/lua/lualib.h>
}
#endif

namespace {
    /**
     * \brief Removes the underscores from a string and
     *  replaces them with spaces.
     * \details Utility for the prototype parser used by
     *  Command.
     * \param[in] s a const reference to the input string
     * \return a copy of \p s with underscores replaced with
     *  spaces
     */
    std::string remove_underscores(const std::string& s) {
        std::string result = s;
        for(GEO::index_t i=0; i<result.size(); ++i) {
            if(result[i] == '_') {
                result[i] = ' ';
            }
        }
        return result;
    }



    /**
     * \brief Removes leading and trailing spaces from
     *  a string.
     * \details Utility for the prototype parser used by
     *  Command.
     * \param[in,out] line the line from which spaces should
     *  be trimmed.
     */
    void trim_spaces(std::string& line) {
        size_t i=0;
        for(i=0; i<line.length(); ++i) {
            if(line[i] != ' ') {
                break;
            }
        }
        size_t j=line.length();
        for(j=line.length(); j>0; --j) {
            if(line[j-1] != ' ') {
                break;
            }
        }
        if(j>i) {
            line = line.substr(i, j-i);
        } else {
            line = "";
        }
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

    /*****************************************************************/

    StatusBar::StatusBar() {
        step_ = 0;
        percent_ = 0;
        progress_ = false;
        canceled_ = false;
        nb_active_ = 0;
    }
    
    void StatusBar::begin() {
        progress_ = true;
        canceled_ = false;
        ++nb_active_;
    }

    void StatusBar::progress(GEO::index_t step, GEO::index_t percent) {
        step_ = step;
        percent_ = percent;
	update();
    }

    void StatusBar::end(bool canceled) {
        geo_argused(canceled);
        step_ = 0;
        percent_ = 0;
        progress_ = false;
        --nb_active_;
    }

    void StatusBar::draw() {
        ImGui::Begin(
            "##Status", nullptr,
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse |
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoScrollbar	    
        );
        if(progress_) {
            if(ImGui::Button("cancel")) {
                Progress::cancel();
            }
            ImGui::SameLine();            
            ImGui::Text(
                "%s", Progress::current_task()->task_name().c_str()
            );
            ImGui::SameLine();
            
            std::string overlay =
                String::to_string(step_) + "/" +
                String::to_string(
                    Progress::current_task()->max_steps()
                ) + " (" + 
                String::to_string(percent_) +
                "%)";
            
            ImGui::ProgressBar(
                std::max(0.001f, float(percent_)/float(100.0)),
                ImVec2(-1,0.0),
                overlay.c_str()
            );
        }
        ImGui::End();
    }

    void StatusBar::update() {
	glup_viewer_gui_update();
    }
    
    /*****************************************************************/

    Console::Console(bool* visible_flag) :
        visible_flag_(visible_flag),
        completion_callback_(nullptr),
        history_callback_(nullptr),
        history_index_(0),
        max_history_index_(0),
        fixed_layout_(true) {
	input_buf_[0] = '\0';
	console_font_ = nullptr;
	scroll_to_bottom_ = false;
    }

    void Console::notify_error(const std::string& err) {
	geo_argused(err);
	return;
    }
    
    void Console::div(const std::string& value) {
        this->printf("========== %s", value.c_str());
    }
        
    void Console::out(const std::string& value) {
        this->printf("    %s", value.c_str());
    }
        
    void Console::warn(const std::string& value) {
        this->printf("[W] %s", value.c_str());
        if(visible_flag_ != nullptr) {
            *visible_flag_ = true;
        }
    }
        
    void Console::err(const std::string& value) {
        this->printf("[E] %s", value.c_str());
        if(visible_flag_ != nullptr) {
            *visible_flag_ = true;
        }
	notify_error(value);
    }

    void Console::status(const std::string& value) {
	// Do not display error messages twice.
	if(String::string_starts_with(value, "Error:")) {
	    return;
	}
	this->printf("[status] %s", value.c_str());
    }
    
    void Console::clear() {
        buf_.clear();
        line_offsets_.clear();
    }
    
    void Console::printf(const char* fmt, ...) {
        va_list args;
        va_start(args, fmt);
        int old_size = buf_.size();
        buf_.appendfv(fmt, args); 
        va_end(args);
        for (int new_size = buf_.size(); old_size < new_size; old_size++) {
            if (buf_[old_size] == '\n') {
                line_offsets_.push_back(old_size);
            }
        }
        scroll_to_bottom_ = true;
	update();
    }

    void Console::update() {
        glup_viewer_gui_update();	
    }

    static int TextEditCallbackStub(ImGuiTextEditCallbackData* data) {
        Console* console = (Console*)data->UserData;
	return console->TextEditCallback(data);
    }


    int Console::TextEditCallback(ImGuiTextEditCallbackData* data)  {
        switch (data->EventFlag) {
        case ImGuiInputTextFlags_CallbackCompletion: {
	    if(completion_callback_ != nullptr) {
		static const char*
		    completer_word_break_characters = " .(){},+-*/=";
		// Locate beginning of current word
		const char* word_end = data->Buf + data->CursorPos;
		const char* word_start = word_end;
		while (word_start > data->Buf)
		{
		    const char c = word_start[-1];
		    if(strchr(completer_word_break_characters,c) != nullptr) {
			break;
		    }
		    word_start--;
		}  
		index_t startw = index_t(word_start - data->Buf);
		index_t endw   = index_t(word_end - data->Buf);
		std::string cmpword(word_start, size_t(word_end - word_start));
		std::vector<std::string> matches;
		completion_callback_(
		    this, std::string(data->Buf),
		    startw, endw, cmpword, matches
		);
		if(matches.size() == 0) {
		    this->printf("Completions: no match\n");
		} else if(matches.size() == 1) {
		    // Single match. Delete the beginning of the word and
		    // replace it entirely so we've got nice casing
                    data->DeleteChars(
			(int)(word_start-data->Buf),
			(int)(word_end-word_start)
		    );
                    data->InsertChars(data->CursorPos, matches[0].c_str());
		} else {
		    // Several matches, find longest common prefix
		    std::string longest_prefix;
		    size_t cur_char = 0;
		    bool finished = false;
		    while(!finished) {
			char c = '\0';
			for(size_t i=0; i<matches.size(); ++i) {
			    if(
				cur_char >= matches[i].length() ||
				(i != 0 && matches[i][cur_char] != c)
			    ) {
				finished = true;
				break;
			    }
			    c = matches[i][cur_char];
			}
			if(!finished) {
			    longest_prefix.push_back(c);
			}
			++cur_char;
		    }
		    // Replace edited text with longest prefix
		    if(longest_prefix.length() != 0) {
			data->DeleteChars(
			    (int)(word_start-data->Buf),
			    (int)(word_end-word_start)
			);
			data->InsertChars(
			    data->CursorPos, longest_prefix.c_str()
			);
		    }
		    this->printf("Completions:\n");
		    for(size_t i=0; i<matches.size(); ++i) {
			this->printf(
			    "[%d] ... %s\n", int(i), matches[i].c_str()
			);
		    }
		}
	    }
	} break;
        case ImGuiInputTextFlags_CallbackHistory: {
	    if(history_callback_ != nullptr) {
		std::string history_command;
		//   Call the callback first, to give it the opportunity to
		// declare the history size.
		history_callback_(this, history_index_, history_command);
		if(max_history_index_ > 0) {
		    int h = int(history_index_);
		    if(data->EventKey == ImGuiKey_UpArrow) {
			--h;
			if(h < 0) {
			    h = int(max_history_index_);
			}
		    } else if(data->EventKey == ImGuiKey_DownArrow) {
			++h;
			if(h > int(max_history_index_)) {
			    h = 0;
			}
		    }
		    {
			history_index_ = index_t(h);
			if(history_index_ == max_history_index_) {
			    history_command = "";
			} else {
			    history_callback_(this, history_index_, history_command);
			}
			int newpos = std::min(
			    data->BufSize-1, int(history_command.length())
			);
			strncpy(data->Buf, history_command.c_str(), size_t(newpos));
			data->Buf[newpos] = '\0';
			data->CursorPos = newpos;
			data->SelectionStart = newpos;
			data->SelectionEnd = newpos;
			data->BufTextLen = newpos;
			data->BufDirty = true;   
		    }
		}
	    }
        } break;
	}
        return 0;
    }

    bool Console::exec_command(const char* command) {
	// Note: history_index_ and max_history_index_ are
	// not managed here. They are managed by the callback.
	return Application::instance()->exec_command(command);
    }
    
    void Console::draw(bool* visible, bool with_window) {
	if(!*visible) {
	    return;
	}
	if(with_window) {
	    ImGui::Begin(
		"Console", visible,
		fixed_layout_ ? (
		    ImGuiWindowFlags_NoResize |
		    ImGuiWindowFlags_NoMove |
		    ImGuiWindowFlags_NoCollapse
		) : 0
	    );
	}
        if (ImGui::Button("Clear")) {
            clear();
        }
        ImGui::SameLine();
        bool copy = ImGui::Button("Copy");
        ImGui::SameLine();
        filter_.Draw("Filter", -100.0f);
        ImGui::Separator();

	if(console_font_ != nullptr) {
	    ImGui::PushFont(console_font_);
	} else if(
	    ImGui::GetIO().Fonts->Fonts.size() >= 4 &&
	    Application::instance() != nullptr
	) {
	    if(Application::instance()->scaling() == 2.0f) {
		console_font_ = ImGui::GetIO().Fonts->Fonts[3];
	    } else {
		console_font_ = ImGui::GetIO().Fonts->Fonts[2];
	    }
	    ImGui::PushFont(console_font_);	    
	}

	
	float scaling = ImGui::GetIO().FontDefault->FontSize / 16.0f;
	
        ImGui::BeginChild(
            "scrolling",
	    ImVec2(
		0.0f,
		-20.0f * scaling
	    ),
	    false,
            ImGuiWindowFlags_HorizontalScrollbar
        );
        if (copy) {
            ImGui::LogToClipboard();
        }
            
        {
            const char* buf_begin = buf_.begin();
            const char* line = buf_begin;
            for (int line_no = 0; line != nullptr; line_no++) {
                const char* line_end =
                    (line_no < line_offsets_.Size) ?
                    buf_begin + line_offsets_[line_no] : nullptr;
                
                if (!filter_.IsActive() ||
		    filter_.PassFilter(line, line_end)
		) {
		    bool is_error = ((line_end - line) >= 3 &&
			line[0] == '[' &&
			line[2] == ']' &&
			(line[1] == 'E' || line[1] == 'W')
		    );
		    if(is_error) {
			ImGui::PushStyleColor(
			    ImGuiCol_Text,ImVec4(1.0f,0.0f,0.0f,1.0f)
			);
		    }
                    ImGui::TextUnformatted(line, line_end);
		    if(is_error) {
			ImGui::PopStyleColor();
		    }
                }
                line = line_end && line_end[1] ? line_end + 1 : nullptr;
            }
        } 
        if (scroll_to_bottom_) {
            ImGui::SetScrollHere();
        }
        scroll_to_bottom_ = false;
        ImGui::EndChild();

	ImGui::Text(">>>");
	ImGui::SameLine();
	ImGui::PushItemWidth(-20);
        if(ImGui::InputText(
	       "##CommandInput", input_buf_, geo_imgui_string_length,
	       ImGuiInputTextFlags_EnterReturnsTrue |
	       ImGuiInputTextFlags_CallbackCompletion |
	       ImGuiInputTextFlags_CallbackHistory,
	       &TextEditCallbackStub, (void*)this)
	) {
            char* input_end = input_buf_+strlen(input_buf_);
            while (input_end > input_buf_ && input_end[-1] == ' ') {
		input_end--;
	    }
            *input_end = 0;
            if (input_buf_[0]) {
                exec_command(input_buf_);
	    }
            strcpy(input_buf_, "");
        }
	ImGui::PopItemWidth();
	
        // Keeping auto focus on the input box
        if (
	    ImGui::IsItemHovered() /* || (
		ImGui::IsRootWindowOrAnyChildFocused() &&
		!ImGui::IsAnyItemActive() && !ImGui::IsMouseClicked(0)
	    ) */
        ) {
            ImGui::SetKeyboardFocusHere(-1); // Auto focus previous widget
	}

	if(console_font_ != nullptr) {
	    ImGui::PopFont();
	}

	if(with_window) {
	    ImGui::End();
	}
    }
    
    /*****************************************************************/
    
    CommandInvoker::CommandInvoker() {
    }

    CommandInvoker::~CommandInvoker() {
    }
    
    /*****************************************************************/

    void Command::flush_queue() {
        if(queued_ != nullptr) {
            // Steal the queued command to avoid
            // infinite recursion.
            SmartPointer<Command> queued = queued_;
            queued_ = nullptr;
            queued->apply();
        }
    }
    
    Command::~Command() {
    }

    Command::Command(const std::string& prototype_in) : visible_(false) {

        //   If there is no brace, then protype only has function
        // name, and then the invoker will construct the arguments
        // with names arg1, arg2 ...
        auto_create_args_ =
            (prototype_in.find('(') == std::string::npos);

        if(auto_create_args_) {
            name_ = prototype_in;
            if(name_ == "") {
                name_ = "command";
            }
            help_ = "Hay You Programmer ! No prototype was specified, \n"
                "see Command::make_current() documentation\n"
                "in geogram_gfx/glup_viewer/command.h\n"
                "to specify parameter names (and tooltips)";
            return;
        }
        
        // Parsing the prototype...
        
        std::string prototype = prototype_in;

        // Transform carriage returns into spaces
        {
            for(index_t i=0; i<prototype.size(); ++i) {
                if(prototype[i] == '\n') {
                    prototype[i] = ' ';
                }
            }
        }
        
        // Separate function name from argument list
        size_t p1 = std::string::npos;
        size_t p2 = std::string::npos;
        {
            int level = 0;
            for(size_t i=0; i<prototype.length(); ++i) {
                switch(prototype[i]) {
                case '[':
                    ++level;
                    break;
                case ']':
                    --level;
                    break;
                case '(':
                    if(level == 0) {
                        p1 = i;
                    }
                    break;
                case ')':
                    if(level == 0) {
                        p2 = i;
                    }
                    break;
                }
            }
        }

        
        geo_assert(p1 != std::string::npos && p2 != std::string::npos);
        name_ = prototype.substr(0,p1);
        
        // Trim spaces, and remove return type if it was specified.
        {
            std::vector<std::string> name_parts;
            String::split_string(name_, ' ', name_parts);
            name_ = name_parts[name_parts.size()-1];
        }

        // Find help if any
        {
            size_t bq1 = prototype.find('[',p2);
            size_t bq2 = prototype.find(']',p2);
            if(bq1 != std::string::npos && bq2 != std::string::npos) {
                help_ = prototype.substr(bq1+1, bq2-bq1-1);
            }
        }
        
        std::string args_string = prototype.substr(p1+1,p2-p1-1);
        std::vector<std::string> args;
        String::split_string(args_string, ',', args);
        
        for(index_t i=0; i<args.size(); ++i) {
            std::string arg = args[i];
            std::string default_value;
            std::string help;
            
            // Find help if any
            {
                size_t bq1 = arg.find('[');
                size_t bq2 = arg.find(']');
                if(bq1 != std::string::npos && bq2 != std::string::npos) {
                    help = arg.substr(bq1+1, bq2-bq1-1);
                    arg = arg.substr(0, bq1);
                }
            }
            
            // Find default value if any (to the right of the '=' sign)
            {
                size_t eq = arg.find('=');
                if(eq != std::string::npos) {
                   default_value = arg.substr(eq+1, arg.length()-eq-1);
                   trim_spaces(default_value);
                   arg = arg.substr(0, eq);
                }
            }

            // Analyze argument type and name
            std::vector<std::string> arg_words;
            String::split_string(arg, ' ', arg_words);

            // Argument name is the last word
            const std::string& arg_name = arg_words[arg_words.size()-1];
            int type = -1;
            if(arg_words.size() > 1) {
                bool is_unsigned = false;
                for(index_t w=0; w<arg_words.size()-1; ++w) {
                    if(arg_words[w] == "unsigned") {
                        is_unsigned = true;
                    } else if(arg_words[w] == "bool") {
                        type = Arg::ARG_BOOL;
                    } else if(arg_words[w] == "int") {
                        type = (is_unsigned) ? Arg::ARG_UINT : Arg::ARG_INT;
                    } else if(
                        arg_words[w] == "index_t" ||
                        arg_words[w] == "GEO::index_t"                        
                    ) {
                        type = Arg::ARG_UINT;
                    } else if(
                        arg_words[w] == "float" ||
                        arg_words[w] == "double"                        
                    ) {
                        type = Arg::ARG_FLOAT;
                    } else if(
                        arg_words[w] == "string" ||
                        arg_words[w] == "std::string" ||
                        arg_words[w] == "string&" ||
                        arg_words[w] == "std::string&"
                    ) {
                        type = Arg::ARG_STRING;
                    }
                }
            }

            switch(type) {
            case Arg::ARG_BOOL: {
                bool val = false;
                if(default_value != "") {
                    String::from_string(default_value, val);
                }
                add_arg(arg_name, val, help);
            } break;
            case Arg::ARG_INT: {
                int val = 0;
                if(default_value != "") {
                    String::from_string(default_value, val);                    
                }
                add_arg(arg_name, val, help);
            } break;
            case Arg::ARG_UINT: {
                unsigned int val = 0;
                if(default_value != "") {
                    String::from_string(default_value, val);
                }
                add_arg(arg_name, val, help);
            } break;
            case Arg::ARG_FLOAT: {
                float val = 0.0f;
                if(default_value != "") {
                    String::from_string(default_value, val);
                }
                add_arg(arg_name, val, help);
            } break;
            case Arg::ARG_STRING: {
                if(default_value != "") {
                    // Remove quotes
                    default_value = default_value.substr(
                        1, default_value.length()-2
                    );
                    add_arg(arg_name, default_value, help);
                }
            } break;
            default: {
                geo_assert_not_reached;
            }
            }
        }
        name_ = remove_underscores(name_);
    }
    
    void Command::apply() {
        if(invoker_ != nullptr) {
            invoker_->invoke();
        }
    }
    
    int Command::int_arg_by_index(index_t i) const {
        const Arg& arg = find_arg_by_index(i);
        geo_assert(
            arg.type == Arg::ARG_INT ||
            arg.type == Arg::ARG_UINT
        );
        int result = arg.val.int_val;
        if(
            arg.type == Arg::ARG_UINT &&
            result < 0
        ) {
            Logger::warn("Cmd")
                << "Argument " << arg.name
                << "Of type UNIT had a negative value"
                << std::endl;
            result = 0;
        }
        return result;
    }

    unsigned int Command::uint_arg_by_index(index_t i) const {
        const Arg& arg = find_arg_by_index(i);
        geo_assert(
            arg.type == Arg::ARG_INT ||
            arg.type == Arg::ARG_UINT
        );
        int result = arg.val.int_val;
        if(result < 0) {
            Logger::warn("Cmd")
                << "Argument " << arg.name
                << "queried as uint had a negative value"
                << std::endl;
            result = 0;
        }
        return (unsigned int)(result);
    }
    
    void Command::draw() {
        if(ImGui::Button("apply")) {
            queued_ = this;
        }
        if(ImGui::IsItemHovered()) {
            if(help_ == "") {
                ImGui::SetTooltip("apply command");
            } else {
                ImGui::SetTooltip("%s",help_.c_str());                
            }
        }
        ImGui::SameLine();
        if(ImGui::Button("default")) {
            reset_factory_settings();
        }
        if(ImGui::IsItemHovered()) {
            ImGui::SetTooltip("reset factory settings");
        }
        ImGui::Separator();            
        for(index_t i=0; i<args_.size(); ++i) {
            args_[i].draw();
        }
        ImGui::Separator();            
    }

    void Command::reset_factory_settings() {
        for(index_t i=0; i<args_.size(); ++i) {
            args_[i].val = args_[i].default_val;
        }
    }

    void Command::ArgVal::clear() {
        bool_val = false;
        int_val = 0;
        float_val = 0.0f;
        string_val[0] = '\0';
    }

    Command::ArgVal::ArgVal(const Command::ArgVal& rhs) {
        bool_val = rhs.bool_val;
        int_val = rhs.int_val;
        float_val = rhs.float_val;
        Memory::copy(string_val, rhs.string_val, 64);
    }

    Command::ArgVal& Command::ArgVal::operator=(const ArgVal& rhs) {
        if(&rhs != this) {
            bool_val = rhs.bool_val;
            int_val = rhs.int_val;
            float_val = rhs.float_val;
            Memory::copy(string_val, rhs.string_val, 64);
        }
        return *this;
    }

    Command::Arg::Arg() {
        name = "unnamed";
        type = ARG_BOOL;
        val.clear();
        default_val.clear();
        val.bool_val = false;
        default_val.bool_val = false;
    }
    
    Command::Arg::Arg(
        const std::string& name_in, bool x,
        const std::string& help_in
    ) {
        name = name_in;
        help = help_in;
        type = ARG_BOOL;
        val.clear();
        default_val.clear();
        val.bool_val = x;
        default_val.bool_val = x;
    }

    Command::Arg::Arg(
        const std::string& name_in, int x,
        const std::string& help_in
    ) {
        name = name_in;
        help = help_in;                
        type = ARG_INT;
        val.clear();
        default_val.clear();
        val.int_val = x;
        default_val.int_val = x;
    }

    Command::Arg::Arg(
        const std::string& name_in, unsigned int x,
        const std::string& help_in
    ) {
        name = name_in;
        help = help_in;                
        type = ARG_UINT;
        val.clear();
        default_val.clear();
        val.int_val = int(x);
        default_val.int_val = int(x);
    }

    Command::Arg::Arg(
        const std::string& name_in, float x,
        const std::string& help_in                
    ) {
        name = name_in;
        help = help_in;
        type = ARG_FLOAT;
        val.clear();
        default_val.clear();
        val.float_val = x;
        default_val.float_val = x;
    }

    Command::Arg::Arg(
        const std::string& name_in, double x,
        const std::string& help_in                                
    ) {
        name = name_in;
        help = help_in;                
        type = ARG_FLOAT;
        val.clear();
        default_val.clear();
        val.float_val = float(x);
        default_val.float_val = float(x);
    }

    Command::Arg::Arg(
        const std::string& name_in, const std::string& x,
        const std::string& help_in
    ) {
        name = name_in;
        help = help_in;                
        type = ARG_STRING;
        val.clear();
        default_val.clear();
        geo_assert(x.length() < 63);
        Memory::copy(default_val.string_val,x.c_str(), x.length());
        default_val.string_val[x.length()+1] = '\0';        
        Memory::copy(val.string_val,x.c_str(), x.length());        
        val.string_val[x.length()+1] = '\0';        
    }

    void Command::Arg::draw() {
        // Some widgets with labels are too wide,
        // therefore their label is displayed with
        // a separate ImGui::Text().
        //   Each ImGUI widget requires a unique Id,
        // that is normally generated from the label.
        //   If label starts with "##", then ImGui
        // makes it invisible (and uses what's after
        // "##" to generate the Id).
        switch(type) {
        case ARG_BOOL:
            ImGui::Checkbox(remove_underscores(name).c_str(), &val.bool_val);
            if(help != "" && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("%s",help.c_str());
            }
            break;
        case ARG_INT:
            ImGui::Text("%s",remove_underscores(name).c_str());
            if(help != "" && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("%s",help.c_str());
            }
            ImGui::InputInt(("##" + name).c_str(), &val.int_val);
            break;
        case ARG_UINT:
            ImGui::Text("%s",remove_underscores(name).c_str());            
            if(help != "" && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("%s",help.c_str());
            }
            ImGui::InputInt(("##" + name).c_str(), &val.int_val);
            break;
        case ARG_FLOAT:
            ImGui::Text("%s",remove_underscores(name).c_str());
            if(help != "" && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("%s",help.c_str());
            }
            ImGui::InputFloat(("##" + name).c_str(), &val.float_val);
            break;
        case ARG_STRING:
            ImGui::Text("%s",remove_underscores(name).c_str());
            if(help != "" && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("%s",help.c_str());
            }
            ImGui::InputText(("##" + name).c_str(), val.string_val, 64);
            break;
        }
    }

    SmartPointer<Command> Command::current_;
    SmartPointer<Command> Command::queued_;

    /**********************************************************************/

    namespace {
    
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

    }

    /**********************************************************************/

    TextEditor::TextEditor(bool* visible) : visible_(visible) {
	impl_.SetText("\n");
	impl_.SetCursorPosition(
	    ::TextEditor::Coordinates(0,0)
	);
	impl_.SetLanguageDefinition(
	    ::TextEditor::LanguageDefinition::Lua()
	);
	impl_.SetPalette(::TextEditor::GetDarkPalette());	    	    
    }

    std::string TextEditor::text() const {
	return impl_.GetText();
    }
    
    void TextEditor::draw() {
	ImGui::Begin(
	    "Text Editor", visible_,
	    ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse 
	);

	if(Application::instance()->scaling() == 2.0f) {
	    ImGui::PushFont(ImGui::GetIO().Fonts->Fonts[3]);
	} else {
	    ImGui::PushFont(ImGui::GetIO().Fonts->Fonts[2]);
	}

	impl_.Render("##source");
	
	ImGui::PopFont();
	
	ImGui::End();
    }

    void TextEditor::load(const std::string& filename) {
	std::ifstream in(filename.c_str());
	std::string text;
	std::string line;
	while(std::getline(in,line)) {
	    text += line;
	    text += "\n";
	}
	impl_.SetText(text);
	impl_.SetCursorPosition(
	    ::TextEditor::Coordinates(0,0)
	);
    }

    void TextEditor::save(const std::string& filename) {
	std::ofstream out(filename.c_str());
	out << impl_.GetText();
    }

    void TextEditor::clear() {
	impl_.SetText("\n");
	impl_.SetCursorPosition(
	    ::TextEditor::Coordinates(0,0)
	);
    }

    void TextEditor::load_data(const char* data) {
	impl_.SetText(data);
	impl_.SetCursorPosition(
	    ::TextEditor::Coordinates(0,0)
	);
    }
    
    /**********************************************************************/
    
    Application* Application::instance_ = nullptr;
    
    Application::Application(
        int argc, char** argv, const std::string& usage, lua_State* lua_state
    ) :
        argc_(argc),
        argv_(argv),
        usage_(usage),
        load_dialog_(this,false),
        save_dialog_(this,true),
	text_editor_visible_(false),
	text_editor_(&text_editor_visible_)
    {

        name_ = (argc == 0) ? "" : FileSystem::base_name(argv[0]);
        geo_assert(instance_ == nullptr);
        instance_ = this;

        GEO::initialize();
        Logger::instance()->set_quiet(false);

        // Import the arg groups needed by graphics.
        CmdLine::import_arg_group("standard");
        CmdLine::import_arg_group("algo");
        CmdLine::import_arg_group("gfx");

        left_pane_visible_ = true;
        right_pane_visible_ = true;
        console_visible_ = false;
	text_editor_visible_ = false;

        console_ = new Console(&console_visible_);
        status_bar_ = new StatusBar;

        lighting_ = true;
	effect_ = GLUP_VIEWER_NO_EFFECT;

        clip_mode_ = GLUP_CLIP_WHOLE_CELLS;

        geogram_logo_texture_ = 0;

        scaling_ = 1.0f;
	fixed_layout_ = true;
        retina_mode_ = false;

	GEO::CmdLine::declare_arg(
	    "gfx:gui", true, "show/hide GUI"
	);

	GEO::CmdLine::declare_arg(
	    "gfx:background", true, "fancy/plain background"
	);
	
	GEO::CmdLine::declare_arg(
	    "gfx:rotate", "",
	    "initial rotation (four ','-separated doubles: x,y,z,angle)"
	);

	GEO::CmdLine::declare_arg(
	    "gfx:snapshot", "", "snapshot image file name (.ppm)"
	);

	GEO::CmdLine::declare_arg(
	    "gfx:keypress", "", "initial simulated sequence of pressed keys"
	);

#ifdef GEOGRAM_WITH_LUA	
	lua_error_occured_ = false;
	if(lua_state != nullptr) {
	    lua_state_ = lua_state;
	    owns_lua_state_ = false;
	} else {
	    lua_state_ = luaL_newstate();
	    owns_lua_state_ = true;
	    luaL_openlibs(lua_state_);
	    init_lua_io(lua_state_);
	}
	init_lua_glup(lua_state_);
	init_lua_glup_viewer(lua_state_);		
	init_lua_imgui(lua_state_);
#endif
	geo_cite_with_info(
	    "WEB:ImGUI",
	    "Used to create the GUI of GEOGRAM utilities "
	    "(vorpaview, geobox, geocod)."
	);

	background_color_1_ = vec4f(1.0f, 1.0f, 1.0f, 1.0f);	
	background_color_2_ = vec4f(0.0f, 0.0f, 0.7f, 1.0f);
    }

    Application::~Application() {
        if(geogram_logo_texture_ != 0) {
            glDeleteTextures(1, &geogram_logo_texture_);
        }
        for(index_t i=0; i<colormaps_.size(); ++i) {
            if(colormaps_[i].texture != 0) {
                glDeleteTextures(1, &colormaps_[i].texture);                
            }
        }
#ifdef GEOGRAM_WITH_LUA
	if(owns_lua_state_) {
	    lua_close(lua_state_);
	    lua_state_ = nullptr;
	}
#endif	
        geo_assert(instance_ == this);        
        instance_ = nullptr;
    }

    static GLboolean keyboard_func_ext(const char* q, GlupViewerEvent ev) {
	GLboolean result = GL_FALSE;
	switch(ev) {
	    case GLUP_VIEWER_DOWN:
		result = GLboolean(Application::instance()->on_key_pressed(q));
		break;
	    case GLUP_VIEWER_MOVE:
		break;
	    case GLUP_VIEWER_UP:
		result = GLboolean(Application::instance()->on_key_released(q));
		break;
	}
	return result;
    }
    
    void Application::start() {
        std::vector<std::string> filenames;
        if(!GEO::CmdLine::parse(argc_, argv_, filenames, usage_)) {
            return;
        }

        if(filenames.size() == 1 && FileSystem::is_directory(filenames[0])) {
            path_ = filenames[0];
        } else if(filenames.size() > 0) {
            for(index_t i=0; i<filenames.size(); ++i) {
                load(filenames[i]);
            }
            if(filenames.size() > 0) {
                path_ = FileSystem::dir_name(filenames[filenames.size()-1]);
            }
        } else {
            path_ = FileSystem::documents_directory();
        }

        glup_viewer_set_window_title(name_.c_str());
        glup_viewer_set_init_func(init_graphics_callback);
        glup_viewer_set_display_func(draw_scene_callback);
        glup_viewer_set_overlay_func(draw_gui_callback);
        glup_viewer_set_drag_drop_func(dropped_file_callback);
	glup_viewer_set_keyboard_func_ext(keyboard_func_ext);
	
        if(GEO::CmdLine::get_arg_bool("gfx:full_screen")) {
            glup_viewer_enable(GLUP_VIEWER_FULL_SCREEN);
        }

	if(GEO::CmdLine::get_arg_bool("gfx:gui")) {
	    glup_viewer_enable(GLUP_VIEWER_TWEAKBARS);
	} else {
	    glup_viewer_disable(GLUP_VIEWER_TWEAKBARS);	    
	}

	if(GEO::CmdLine::get_arg_bool("gfx:background")) {
	    glup_viewer_enable(GLUP_VIEWER_BACKGROUND);	    
	} else {
	    glup_viewer_disable(GLUP_VIEWER_BACKGROUND);	    	    
	}

        glup_viewer_main_loop(argc_, argv_);
    }

    void Application::quit() {
	glup_viewer_exit_main_loop();
    }
    
    bool Application::save(const std::string& filename) {
        Logger::warn("GLUP") << "Could not save " << filename << std::endl;
        Logger::warn("GLUP") << "Application::save() needs to be overloaded"
                             << std::endl;
        return false;
    }
    
    bool Application::load(const std::string& filename) {
        Logger::warn("GLUP") << "Could not load " << filename << std::endl;
        Logger::warn("GLUP") << "Application::load() needs to be overloaded"
                             << std::endl;
        return false;
    }

    bool Application::can_load(const std::string& filename) {
        std::string extensions_str = supported_read_file_extensions();
        if(extensions_str == "") {
            return false;
        }
        if(extensions_str == "*") {
            return true;
        }
        std::string extension = FileSystem::extension(filename);
        std::vector<std::string> extensions;
        String::split_string(extensions_str, ';', extensions);
        for(index_t i=0; i<extensions.size(); ++i) {
            if(extensions[i] == extension) {
                return true;
            }
        }
        return false;
    }

    std::string Application::supported_read_file_extensions() {
        return "";
    }

    std::string Application::supported_write_file_extensions() {
        return "";
    }

    ImTextureID Application::convert_to_ImTextureID(GLuint gl_texture_id_in) {
        // It is not correct to directly cast a GLuint into a void*
        // (generates warnings), therefore I'm using a union.
        union {
            GLuint gl_texture_id;
            ImTextureID imgui_texture_id;
        };
        imgui_texture_id = nullptr;
        gl_texture_id = gl_texture_id_in;
        return imgui_texture_id;
    }
    
    void Application::browse(const std::string& path, bool subdirs) {
        std::vector<std::string> files;
        GEO::FileSystem::get_directory_entries(path,files);
        
        for(GEO::index_t i=0; i<files.size(); ++i) {
            if(GEO::FileSystem::is_directory(files[i]) && subdirs) {
                if(ImGui::BeginMenu(path_to_label(path_,files[i]).c_str())) {
                    browse(files[i]);
                    ImGui::EndMenu();
                }
            } else {
                if(can_load(files[i])) {
                    if(ImGui::MenuItem(path_to_label(path_,files[i]).c_str())) {
                        load(files[i]);
                    }
                }
            }
        }
    }
    
    void Application::init_graphics() {
        Logger::instance()->register_client(console_);
        Progress::set_client(status_bar_);
        
        GEO::Graphics::initialize();
        
        glup_viewer_add_toggle(
            'a', glup_viewer_is_enabled_ptr(GLUP_VIEWER_IDLE_REDRAW), "animate"
        );
        glup_viewer_add_toggle(
            'T', glup_viewer_is_enabled_ptr(GLUP_VIEWER_TWEAKBARS), "tweakbars"
        );

        glup_viewer_add_key_func('z', zoom_in, "Zoom in");
        glup_viewer_add_key_func('Z', zoom_out, "Zoom out");
        glup_viewer_add_toggle('L', &lighting_, "lighting");        
        
#ifdef GEO_OS_EMSCRIPTEN
        {
            std::vector<std::string> all_files;
            GEO::FileSystem::get_directory_entries("/",all_files);
            if(all_files.size() > 0 && can_load(all_files[0])) {
                load(all_files[0]);
            }
        }
#endif

        glGenTextures(1, &geogram_logo_texture_);
        glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
        glBindTexture(GL_TEXTURE_2D, geogram_logo_texture_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2DXPM(geogram_logo_256_xpm);
	
        if(glup_viewer_is_high_dpi()) {
            retina_mode_ = true;                
            scaling_ = 2.0f;
        }

	std::string rotate = GEO::CmdLine::get_arg("gfx:rotate");
	if(rotate != "") {
	    float axis[3];
	    float angle;
	    sscanf(
		rotate.c_str(),
		"%f,%f,%f,%f",
		&axis[0], &axis[1], &axis[2], &angle
	    );
	    glup_viewer_set_scene_rotation(axis,angle);
	}

	std::string keys = GEO::CmdLine::get_arg("gfx:keypress");
	for(size_t i=0; i<keys.length(); ++i) {
	    glup_viewer_char_callback(nullptr,(unsigned int)(keys[i]));
	}
    }

    void Application::init_graphics_callback() {
        if(instance() != nullptr) {
            instance()->init_graphics();
        }
    }
    
    void Application::draw_scene() {
    }

    void Application::draw_scene_callback() {
        if(instance() != nullptr) {
	    
	    glup_viewer_set_background_color(
		instance()->background_color_1_.x,
		instance()->background_color_1_.y,
		instance()->background_color_1_.z		
	    );
	
	    glup_viewer_set_background_color2(
		instance()->background_color_2_.x,
		instance()->background_color_2_.y,
		instance()->background_color_2_.z		
	    );
	    
            if(instance()->lighting_) {
                glupEnable(GLUP_LIGHTING);
            } else {
                glupDisable(GLUP_LIGHTING);                
            }
            glupClipMode(instance()->clip_mode_);
            instance()->draw_scene();
	    std::string snapshot = CmdLine::get_arg("gfx:snapshot");
	    if(snapshot != "") {
		glup_viewer_snapshot(snapshot.c_str());
		exit(0);
	    }
        }
    }

    void Application::draw_gui() {
	// TODO: better management of font size (for now we use two
	// hardwired sizes)
	// Fonts are created in glup_viewer_gui_init(), in glup_viewer_gui_private.cpp
	
	//ImGui::GetIO().FontGlobalScale = scaling_;
	if(scaling_ == 2.0f) {
	    ImGui::GetIO().FontDefault = ImGui::GetIO().Fonts->Fonts[1];
	} else {
	    ImGui::GetIO().FontDefault = ImGui::GetIO().Fonts->Fonts[0];
	}
	
        draw_menu_bar();
	if(text_editor_visible_) {
	    int w,h;
	    glup_viewer_get_screen_size(&w, &h);
	    if(status_bar_->active()) {
		h -= (STATUS_HEIGHT()+1);
	    }
	    if(console_visible_) {
		h -= (CONSOLE_HEIGHT()+1);
	    }
	    h -= MENU_HEIGHT();
	    ImGui::SetNextWindowPos(
		ImVec2(0.0f, float(MENU_HEIGHT())),
		fixed_layout_ ? ImGuiCond_Always : ImGuiCond_Once
	    );
	    ImGui::SetNextWindowSize(
		ImVec2(float(w), float(h)),
		fixed_layout_ ? ImGuiCond_Always : ImGuiCond_Once		
	    );
	    text_editor_.draw();
	    if(console_visible_) {
		draw_console();
	    }
	    if(status_bar_->active()) {
		draw_status_bar();
	    }
	} else {
	    if(left_pane_visible_) {
		draw_left_pane();
	    }
	    if(right_pane_visible_) {
		draw_right_pane();
	    }
	    if(console_visible_) {
		draw_console();
	    }
	    if(status_bar_->active()) {
		draw_status_bar();
	    }
	}
        load_dialog_.draw();
        save_dialog_.draw();
    }

    void Application::draw_left_pane() {
        int w,h;
        glup_viewer_get_screen_size(&w, &h);
        if(status_bar_->active()) {
            h -= (STATUS_HEIGHT()+1);
        }
        if(console_visible_) {
            h -= (CONSOLE_HEIGHT()+1);
        }
        h -= MENU_HEIGHT();

        if(Command::current() != nullptr) {
            h /= 2;
        }
        
	if(fixed_layout_) {
	    ImGui::SetNextWindowPos(
		ImVec2(0.0f, float(MENU_HEIGHT())),
		ImGuiCond_Always 
	    );
	    ImGui::SetNextWindowSize(
		ImVec2(float(PANE_WIDTH()), float(h)),
		ImGuiCond_Always 
	    );
	} else {
	    ImGui::SetNextWindowPos(
		ImVec2(float(MENU_HEIGHT()), 2.0f*float(MENU_HEIGHT())),
		ImGuiCond_Once			    
	    );
	    ImGui::SetNextWindowSize(
		ImVec2(float(PANE_WIDTH()), float(h)/2.0f),
		ImGuiCond_Once			    	    
	    );
	}

        draw_viewer_properties_window();

        if(Command::current() != nullptr) {
            ImGui::SetNextWindowPos(
                ImVec2(0.0f, float(MENU_HEIGHT()+h+1)),
                ImGuiCond_Always
            );
            ImGui::SetNextWindowSize(
                ImVec2(float(PANE_WIDTH()), float(h-1)),
                ImGuiCond_Always
            );
            draw_command();
        }
    }

    void Application::draw_right_pane() {
        int w,h;
        glup_viewer_get_screen_size(&w, &h);
        if(status_bar_->active()) {
            h -= (STATUS_HEIGHT()+1);
        }
        if(console_visible_) {
            h -= (CONSOLE_HEIGHT()+1);
        }
        h -= MENU_HEIGHT();

        ImGui::SetNextWindowPos(
            ImVec2(float(w-PANE_WIDTH()), float(MENU_HEIGHT())),
	    fixed_layout_ ? ImGuiCond_Always : ImGuiCond_Once	    
        );
        
        ImGui::SetNextWindowSize(
            ImVec2(float(PANE_WIDTH()), float(h)),
	    fixed_layout_ ? ImGuiCond_Always : ImGuiCond_Once	    
        );
        
        draw_object_properties_window();
    }

    void Application::draw_viewer_properties_window() {
        ImGui::Begin(
            "Viewer",
	    fixed_layout_ ? &left_pane_visible_ : nullptr,
            fixed_layout_ ? (
		ImGuiWindowFlags_NoResize |
		ImGuiWindowFlags_NoMove |
		ImGuiWindowFlags_NoCollapse
	    ) : 0
        );
        draw_viewer_properties();
        ImGui::End();
    }

    void Application::draw_viewer_properties() {
        if(ImGui::Button("home [H]", ImVec2(-1,0))) {
            glup_viewer_home();
        }
        ImGui::Separator();
        ImGui::Checkbox(
            "Lighting [L]", &lighting_
        );
        if(lighting_) {
            ImGui::Checkbox(
                "edit light [l]",
                (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_ROTATE_LIGHT)
            );
        }
    
        ImGui::Separator();
        ImGui::Checkbox(
            "Clipping [F1]", (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_CLIP)
        );
        if(glup_viewer_is_enabled(GLUP_VIEWER_CLIP)) {
            ImGui::Combo(
                "mode",
                (int*)&clip_mode_,                
                "std. GL\0cells\0straddle\0slice\0\0"
            );
            ImGui::Checkbox(
                "edit clip [F2]",
                (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_EDIT_CLIP)
            );
            ImGui::Checkbox(
                "fixed clip [F3]",
                (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_FIXED_CLIP)
            );
        }
    
        ImGui::Separator();
        ImGui::Text("Style");
	if(ImGui::ColorEdit3WithPalette(
	       "Background 1", background_color_1_.data())
	) {
	    glup_viewer_enable(GLUP_VIEWER_BACKGROUND);
	}
	if(ImGui::ColorEdit3WithPalette(
	       "Background 2", background_color_2_.data())
	) {
	    glup_viewer_enable(GLUP_VIEWER_BACKGROUND);	    
	}
	if(
	  ImGui::Combo(
	    "sfx",
	    (int*)&effect_,
	    "none\0SSAO\0cartoon\0\0"
	  )
	) {
	    if(!glup_viewer_set_effect(effect_)) {
		Logger::err("GLUP") << "Could not activate effect"
				    << std::endl;
	    }
	}
    }
    
    void Application::draw_object_properties_window() {
        ImGui::Begin(
            "Object", &right_pane_visible_,
	    fixed_layout_ ? (
		ImGuiWindowFlags_NoResize |
		ImGuiWindowFlags_NoMove |
		ImGuiWindowFlags_NoCollapse
	    ) : 0
        );
        draw_object_properties();
        ImGui::End();
    }

    void Application::draw_object_properties() {
        ImGui::Separator();
        ImGui::Text("Object properties...");
    }
    
    void Application::draw_command() {
        geo_assert(Command::current() != nullptr);
        if(!Command::current()->is_visible()) {
            Command::reset_current();
            return;
        }
        ImGui::Begin(
            Command::current()->name().c_str(),
            Command::current()->is_visible_ptr(),
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse 
        );
        Command::current()->draw();
        ImGui::End();
    }
    
    void Application::draw_console() {
        int w,h;
        glup_viewer_get_screen_size(&w, &h);
        h = h - CONSOLE_HEIGHT();
        if(!fixed_layout_ || status_bar_->active()) {
            h -= (STATUS_HEIGHT() + 1);
        }
	if(fixed_layout_) {
	    ImGui::SetNextWindowPos(
		ImVec2(0.0f, float(h)),
		ImGuiCond_Always 
	    );
	    ImGui::SetNextWindowSize(
		ImVec2(float(w),float(CONSOLE_HEIGHT())),
		ImGuiCond_Always 
	    );
	} else {
	    ImGui::SetNextWindowPos(
		ImVec2(float(MENU_HEIGHT()), float(h) - float(MENU_HEIGHT())),
		ImGuiCond_Once
	    );
	    ImGui::SetNextWindowSize(
		ImVec2(
		    float(w) - 2.0f * float(MENU_HEIGHT()),
		    float(CONSOLE_HEIGHT())
		),
		ImGuiCond_Once	    
	    );
	}
	console_->set_fixed_layout(fixed_layout_);
        console_->draw(&console_visible_);
    }

    void Application::draw_status_bar() {
        int w,h;
        glup_viewer_get_screen_size(&w, &h);
        ImGui::SetNextWindowPos(
            ImVec2(0.0f, float(h-STATUS_HEIGHT())),
            ImGuiCond_Always
        );
        ImGui::SetNextWindowSize(
            ImVec2(float(w),float(STATUS_HEIGHT()-1)),
            ImGuiCond_Always
        );
        status_bar_->draw();
    }
    
    void Application::draw_menu_bar() {
        if(ImGui::BeginMainMenuBar()) {
            if(ImGui::BeginMenu("File")) {
                if(supported_read_file_extensions() != "") {
                    draw_load_menu();
                }
#ifndef GEO_OS_EMSCRIPTEN		
		if(current_file_ != "") {
		    if(ImGui::MenuItem("Save")) {
			if(save(current_file_)) {
			    Logger::out("I/O") << "Saved "
					       << current_file_ << std::endl;
			} else {
			    Logger::out("I/O") << "Could not save "
					       << current_file_ << std::endl;
			}
		    }
		}
#endif		
                if(supported_write_file_extensions() != "") {
                    draw_save_menu();
                }
		draw_fileops_menu();
                draw_about();
#ifndef GEO_OS_EMSCRIPTEN                        
                ImGui::Separator();
                if(ImGui::MenuItem("quit [q]")) {
                    glup_viewer_exit_main_loop();
                }
#endif
                ImGui::EndMenu();
            }
            draw_windows_menu();
            draw_application_menus();
            
            ImGui::EndMainMenuBar();            
        }
    }

    void Application::draw_application_menus() {
        // Meant to be overloaded in derived classes.
    }
    
    void Application::draw_load_menu() {
#ifdef GEO_OS_EMSCRIPTEN
            ImGui::Text("To load a file,");
            ImGui::Text("use the \"Browse\"");
            ImGui::Text("button on the top");
            ImGui::Text("(or \"recent files\"");
            ImGui::Text("below)");
            ImGui::Separator();
            if(ImGui::BeginMenu("Recent files...")) {
                browse(path_);
                ImGui::EndMenu(); 
            }
#else
        if(ImGui::MenuItem("Load...")) {
            load_dialog_.show();
        }
#endif        
    }

    void Application::draw_save_menu() {
#ifdef GEO_OS_EMSCRIPTEN
        if(ImGui::BeginMenu("Save as...")) {
	    ImGui::MenuItem("Supported extensions:", nullptr, false, false);
            std::vector<std::string> extensions;
            String::split_string(
                supported_write_file_extensions(), ';', extensions
            );
            for(index_t i=0; i<extensions.size(); ++i) {
		ImGui::MenuItem(
		    (" ." + extensions[i]).c_str(), nullptr, false, false
		);	    		
	    }
	    ImGui::Separator();
	    static char buff[geo_imgui_string_length];
	    if(current_file_ != "") {
		strcpy(buff, current_file_.c_str());
	    } else if (extensions.size() != 0) {
		strcpy(buff, ("out." + extensions[0]).c_str());		
	    }

	    if(ImGui::InputText(
		   "##MenuFileName",buff,geo_imgui_string_length,
		   ImGuiInputTextFlags_EnterReturnsTrue)
	    ) {
		current_file_ = buff;
		if(String::string_starts_with(current_file_, "/")) {
		    current_file_ = current_file_.substr(
			1,current_file_.length()-1
		    );
		}
		if(save(current_file_)) {
		    std::string command =
			"saveFileFromMemoryFSToDisk(\'" +
			current_file_ +
			"\');" ;
		    emscripten_run_script(command.c_str());
		}
	    }
            ImGui::EndMenu();
        }
#else        
        if(ImGui::MenuItem("Save as...")) {
            save_dialog_.show();
        }
#endif        
    }

    void Application::draw_fileops_menu() {
    }
    
    void Application::draw_about() {
        ImGui::Separator();
        if(ImGui::BeginMenu("About...")) {
            ImGui::Text("%s : a GEOGRAM application", name_.c_str());
            ImGui::Image(
                convert_to_ImTextureID(geogram_logo_texture_),
                ImVec2(256.0f * scaling(), 256.0f * scaling())
            );
            ImGui::Text("\n");            
            ImGui::Separator();
            ImGui::Text("\n");
            ImGui::Text("GEOGRAM website: ");
            ImGui::Text("http://alice.loria.fr/software/geogram");

            ImGui::Text("\n");
            ImGui::Separator();
            ImGui::Text(
                "%s",
                (
                    "GEOGRAM version:" +
                    Environment::instance()->get_value("version")
                ).c_str()
            );
            ImGui::EndMenu();
        }
    }

    void Application::draw_windows_menu() {
        if(ImGui::BeginMenu("Windows")) {
            ImGui::MenuItem("object properties", nullptr, &right_pane_visible_);
            ImGui::MenuItem("viewer properties", nullptr, &left_pane_visible_);
            ImGui::MenuItem("console", nullptr, &console_visible_);
            ImGui::MenuItem("text editor", nullptr, &text_editor_visible_);
	    
            if(ImGui::MenuItem(
                "show/hide GUI [T]", nullptr,
                (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_TWEAKBARS)
            )) {
                glup_viewer_post_redisplay();
            }
            ImGui::Separator();
            if(ImGui::MenuItem("Big text", nullptr, &retina_mode_)) {
                scaling_ = retina_mode_ ? 2.0f : 1.0f;
                glup_viewer_post_redisplay();
            }
            ImGui::EndMenu();
        }
    }
    
    void Application::draw_gui_callback() {
        if(instance() != nullptr) {
            instance()->draw_gui();
        }
    }

    void Application::dropped_file_callback(char* filename) {
        if(instance() != nullptr) {
            instance()->load(std::string(filename));
        }
    }

    void Application::init_colormap(
        const std::string& name, const char** xpm_data
    ) {
        colormaps_.push_back(ColormapInfo());
        colormaps_.rbegin()->name = name;
        glGenTextures(1, &colormaps_.rbegin()->texture);
        glBindTexture(GL_TEXTURE_2D, colormaps_.rbegin()->texture);
        glTexImage2DXPM(xpm_data);
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(
            GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR
        );
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void Application::init_colormaps() {
        init_colormap("french", french_xpm);
        init_colormap("black_white", black_white_xpm);
        init_colormap("viridis", viridis_xpm);
        init_colormap("rainbow", rainbow_xpm);
        init_colormap("cei_60757", cei_60757_xpm);
        init_colormap("inferno", inferno_xpm);
        init_colormap("magma", magma_xpm);
        init_colormap("parula", parula_xpm);
        init_colormap("plasma", plasma_xpm);
        init_colormap("blue_red", blue_red_xpm);
    }

    /**********************************************************************/

    bool Application::exec_command(const char* command) {
#ifdef GEOGRAM_WITH_LUA	
	if(luaL_dostring(lua_state_,command)) {
	    adjust_lua_glup_state(lua_state_);
	    const char* msg = lua_tostring(lua_state_,-1);
	    const char* msg2 = strchr(msg,']');
	    if(msg2 != nullptr) {
		msg = msg2+2;
	    }
	    Logger::err("LUA") << "line " << msg << std::endl;
	    lua_error_occured_ = true;
	} else {
	    lua_error_occured_ = false;
	}
	return !lua_error_occured_;
#else
	geo_argused(command);
	Logger::err("LUA") << "Compiled without LUA support"
			   << std::endl;
	return false;
#endif	
    }

    bool Application::on_key_pressed(const char* q) {
	geo_argused(q);
	return false;
    }

    bool Application::on_key_released(const char* q) {
	geo_argused(q);
	return false;
    }

    /**********************************************************************/

    int Application::MENU_HEIGHT() const {
	return int(20 * scaling());
    }

    int Application::PANE_WIDTH() const {
	return int(140 * scaling());
    }

    int Application::CONSOLE_HEIGHT() const {
	return int(200 * scaling());
    }

    int Application::STATUS_HEIGHT() const {
	return retina_mode() ? 48 : 35;
    }
    
    /**********************************************************************/

    SimpleMeshApplication::SimpleMeshApplication(
        int argc, char** argv, const std::string& usage
    ) : Application(argc, argv, usage) {
        std::vector<std::string> extensions;
        GEO::MeshIOHandlerFactory::list_creators(extensions);
        file_extensions_ = String::join_strings(extensions, ';');

        save_dialog_.set_default_filename("out.meshb");
	
        anim_speed_ = 1.0f;
        anim_time_ = 0.0f;

        show_vertices_ = false;
        show_vertices_selection_ = true;
        vertices_size_ = 1.0f;
	vertices_color_ = vec4f(0.0f, 1.0f, 0.0f, 1.0f);
	
        show_surface_ = true;
        show_surface_sides_ = false;
        show_mesh_ = true;
	mesh_color_ = vec4f(0.0f, 0.0f, 0.0f, 1.0f);
	mesh_width_ = 0.1f;
	
        show_surface_borders_ = false;
	surface_color_ =   vec4f(0.0f, 0.5f, 1.0f, 1.0f);
	surface_color_2_ = vec4f(1.0f, 0.5f, 0.0f, 1.0f); 
	
        show_volume_ = false;
	volume_color_ = vec4f(0.9f, 0.9f, 0.9f, 1.0f);	
        cells_shrink_ = 0.0f;
        show_colored_cells_ = false;
        show_hexes_ = true;
	show_connectors_ = true;

        GEO::CmdLine::declare_arg(
            "attributes", true, "load mesh attributes"
        );

        GEO::CmdLine::declare_arg(
            "single_precision", true, "use single precision vertices (FP32)"
        );

        show_attributes_ = false;
        current_colormap_texture_ = 0;
        attribute_min_ = 0.0f;
        attribute_max_ = 0.0f;
        attribute_ = "vertices.point_fp32[0]";
        attribute_name_ = "point_fp32[0]";
        attribute_subelements_ = MESH_VERTICES;
    }

    std::string SimpleMeshApplication::supported_read_file_extensions() {
        return file_extensions_;
    }
    
    std::string SimpleMeshApplication::supported_write_file_extensions() {
        return file_extensions_;
    }

    void SimpleMeshApplication::autorange() {
        if(attribute_subelements_ != MESH_NONE) {
            attribute_min_ = 0.0;
            attribute_max_ = 0.0;
            const MeshSubElementsStore& subelements =
                mesh_.get_subelements_by_type(attribute_subelements_);
            ReadOnlyScalarAttributeAdapter attribute(
                subelements.attributes(), attribute_name_
            );
            if(attribute.is_bound()) {
                attribute_min_ = Numeric::max_float32();
                attribute_max_ = Numeric::min_float32();
                for(index_t i=0; i<subelements.nb(); ++i) {
                    attribute_min_ =
                        std::min(attribute_min_, float(attribute[i]));
                    attribute_max_ =
                        std::max(attribute_max_, float(attribute[i]));
                }
            } 
        }
    }

    std::string SimpleMeshApplication::attribute_names() {
        return mesh_.get_scalar_attributes();
    }

    void SimpleMeshApplication::set_attribute(const std::string& attribute) {
        attribute_ = attribute;
        std::string subelements_name;
        String::split_string(
            attribute_, '.',
            subelements_name,
            attribute_name_
        );
        attribute_subelements_ =
            mesh_.name_to_subelements_type(subelements_name);
        if(attribute_min_ == 0.0f && attribute_max_ == 0.0f) {
            autorange();
        } 
    }
    
    void SimpleMeshApplication::draw_object_properties() {
        ImGui::Checkbox("attributes", &show_attributes_);
        if(show_attributes_) {
            if(attribute_min_ == 0.0f && attribute_max_ == 0.0f) {
                autorange();
            } 
            if(ImGui::Button(
                   (attribute_ + "##Attribute").c_str(), ImVec2(-1,0))
            ) {
                ImGui::OpenPopup("##Attributes");                
            }
            if(ImGui::BeginPopup("##Attributes")) {
                std::vector<std::string> attributes;
                String::split_string(attribute_names(), ';', attributes);
                for(index_t i=0; i<attributes.size(); ++i) {
                    if(ImGui::Button(attributes[i].c_str())) {
                        set_attribute(attributes[i]);
                        ImGui::CloseCurrentPopup();
                    }
                }
                ImGui::EndPopup();                
            }
            ImGui::InputFloat("min",&attribute_min_);
            ImGui::InputFloat("max",&attribute_max_);
            if(ImGui::Button("autorange", ImVec2(-1,0))) {
                autorange();
            }
            if(ImGui::ImageButton(
                   convert_to_ImTextureID(current_colormap_texture_),
                   ImVec2(115.0f*scaling(),8.0f*scaling()))
            ) {
                ImGui::OpenPopup("##Colormap");
            }
            if(ImGui::BeginPopup("##Colormap")) {
                for(index_t i=0; i<colormaps_.size(); ++i) {
                    if(ImGui::ImageButton(
                           convert_to_ImTextureID(colormaps_[i].texture),
                           ImVec2(100.0f*scaling(),8.0f*scaling()))
                    ) {
                        current_colormap_texture_ = colormaps_[i].texture;
                        ImGui::CloseCurrentPopup();
                    }
                }
                ImGui::EndPopup();
            }
        }
        
        if(mesh_.vertices.dimension() >= 6) {
            ImGui::Separator();
            ImGui::Checkbox(
                "Animate [a]",
                (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_IDLE_REDRAW)
            );
            ImGui::SliderFloat("spd.", &anim_speed_, 1.0f, 10.0f, "%.1f");
            ImGui::SliderFloat("t.", &anim_time_, 0.0f, 1.0f, "%.2f");
        }
    
        ImGui::Separator();    
        ImGui::Checkbox("##VertOnOff", &show_vertices_);
	ImGui::SameLine();
	ImGui::ColorEdit3WithPalette("Vert. [p]", vertices_color_.data());

        if(show_vertices_) {
            ImGui::Checkbox("selection", &show_vertices_selection_);            
            ImGui::SliderFloat("sz.", &vertices_size_, 0.1f, 5.0f, "%.1f");
        }

        if(mesh_.facets.nb() != 0) {
            ImGui::Separator();
            ImGui::Checkbox("##SurfOnOff", &show_surface_);
	    ImGui::SameLine();
	    ImGui::ColorEdit3WithPalette(
		"Surf. [S]", surface_color_.data()
	    );
            if(show_surface_) {
		ImGui::Checkbox("##SidesOnOff", &show_surface_sides_);
		ImGui::SameLine();
		ImGui::ColorEdit3WithPalette(
		    "2sided [c]", surface_color_2_.data()
		);
		
                ImGui::Checkbox("##MeshOnOff", &show_mesh_);
		ImGui::SameLine();
		ImGui::ColorEdit3WithPalette("mesh [m]", mesh_color_.data());

		if(show_mesh_) {
		    ImGui::SliderFloat(
			"wid.", &mesh_width_, 0.1f, 2.0f, "%.1f"
		    );
		}
		
                ImGui::Checkbox("borders [B]", &show_surface_borders_);
            }
        }

        if(mesh_.cells.nb() != 0) {
            ImGui::Separator();
            ImGui::Checkbox("##VolumeOnOff", &show_volume_);
	    ImGui::SameLine();
	    ImGui::ColorEdit3WithPalette("Volume [V]", volume_color_.data());
            if(show_volume_) {
                ImGui::SliderFloat(
                    "shrk.", &cells_shrink_, 0.0f, 1.0f, "%.2f"
                );        
                if(!mesh_.cells.are_simplices()) {
                    ImGui::Checkbox("colored cells [C]", &show_colored_cells_);
                    ImGui::Checkbox("hexes [j]", &show_hexes_);
                }
            }
        }
    }

    void SimpleMeshApplication::increment_anim_time_callback() {
        instance()->anim_time_ = std::min(
            instance()->anim_time_ + 0.05f, 1.0f
        );
    }
        
    void SimpleMeshApplication::decrement_anim_time_callback() {
        instance()->anim_time_ = std::max(
            instance()->anim_time_ - 0.05f, 0.0f
        );
    }

    void SimpleMeshApplication::increment_cells_shrink_callback() {
        instance()->cells_shrink_ = std::min(
            instance()->cells_shrink_ + 0.05f, 1.0f
        );
    }
        
    void SimpleMeshApplication::decrement_cells_shrink_callback() {
        instance()->cells_shrink_ = std::max(
            instance()->cells_shrink_ - 0.05f, 0.0f
        );
    }

    void SimpleMeshApplication::init_graphics() {
        glup_viewer_add_toggle('p', &show_vertices_, "vertices");
        glup_viewer_add_toggle('S', &show_surface_, "surface");
        glup_viewer_add_toggle('c', &show_surface_sides_, "2sided");
        glup_viewer_add_toggle('B', &show_surface_borders_,"borders");
        glup_viewer_add_toggle('m', &show_mesh_, "mesh");
        glup_viewer_add_toggle('V', &show_volume_, "volume");
        glup_viewer_add_toggle('j', &show_hexes_, "hexes");
        glup_viewer_add_toggle('k', &show_connectors_, "connectors");	
        glup_viewer_add_toggle('C', &show_colored_cells_, "colored cells");

        glup_viewer_add_key_func(
            'r', decrement_anim_time_callback, "Decrement time"
        );
        
        glup_viewer_add_key_func(
            't', increment_anim_time_callback, "Increment time"
        );
        
        glup_viewer_add_key_func(
            'x', decrement_cells_shrink_callback, "Decrement shrink"
        );
        
        glup_viewer_add_key_func(
            'w', increment_cells_shrink_callback, "Increment shrink"
        );

        init_colormaps();
        current_colormap_texture_ = colormaps_[3].texture;
        Application::init_graphics();
    }
    
    void SimpleMeshApplication::draw_scene() {

        if(mesh_gfx_.mesh() == nullptr) {
            return;
        }
        
        if(glup_viewer_is_enabled(GLUP_VIEWER_IDLE_REDRAW)) {
            anim_time_ = float(
                sin(double(anim_speed_) * GEO::SystemStopwatch::now())
            );
            anim_time_ = 0.5f * (anim_time_ + 1.0f);
        }
        
        mesh_gfx_.set_lighting(lighting_);
        mesh_gfx_.set_time(double(anim_time_));

        if(show_attributes_) {
            mesh_gfx_.set_scalar_attribute(
                attribute_subelements_, attribute_name_,
                double(attribute_min_), double(attribute_max_),
                current_colormap_texture_, 1
            );
        } else {
            mesh_gfx_.unset_scalar_attribute();
        }
        
        if(show_vertices_) {
            mesh_gfx_.set_points_color(
		vertices_color_.x, vertices_color_.y, vertices_color_.z
	    );
            mesh_gfx_.set_points_size(vertices_size_);
            mesh_gfx_.draw_vertices();
        }

        if(show_vertices_selection_) {
            mesh_gfx_.set_points_color(1.0, 0.0, 0.0);
            mesh_gfx_.set_points_size(2.0f * vertices_size_);
            mesh_gfx_.set_vertices_selection("selection");
            mesh_gfx_.draw_vertices();
            mesh_gfx_.set_vertices_selection("");            
        }

	mesh_gfx_.set_mesh_color(0.0, 0.0, 0.0);

	mesh_gfx_.set_surface_color(
	    surface_color_.x, surface_color_.y, surface_color_.z
	);
        if(show_surface_sides_) {
	    mesh_gfx_.set_backface_surface_color(
		surface_color_2_.x, surface_color_2_.y, surface_color_2_.z
	    );
        }
	
        mesh_gfx_.set_show_mesh(show_mesh_);
	mesh_gfx_.set_mesh_color(mesh_color_.x, mesh_color_.y, mesh_color_.z);
	mesh_gfx_.set_mesh_width(index_t(mesh_width_*10.0f));
	
        if(show_surface_) {
	    float specular_backup = glupGetSpecular();
	    glupSetSpecular(0.4f);
            mesh_gfx_.draw_surface();
	    glupSetSpecular(specular_backup);	    
        }
        
        if(show_surface_borders_) {
            mesh_gfx_.draw_surface_borders();
        }

        if(show_mesh_) {
            mesh_gfx_.draw_edges();
        }

        if(show_volume_) {

            if(
                glupIsEnabled(GLUP_CLIPPING) &&
                glupGetClipMode() == GLUP_CLIP_SLICE_CELLS
            ) {
                mesh_gfx_.set_lighting(false);
            }
            
            mesh_gfx_.set_shrink(double(cells_shrink_));
            mesh_gfx_.set_draw_cells(GEO::MESH_HEX, show_hexes_);
            mesh_gfx_.set_draw_cells(GEO::MESH_CONNECTOR, show_connectors_);
	    
            if(show_colored_cells_) {
                mesh_gfx_.set_cells_colors_by_type();
            } else {
                mesh_gfx_.set_cells_color(
		    volume_color_.x, volume_color_.y, volume_color_.z
		);
            }
            mesh_gfx_.draw_volume();

            mesh_gfx_.set_lighting(lighting_);
        }
    }

    bool SimpleMeshApplication::load(const std::string& filename) {
        if(!FileSystem::is_file(filename)) {
            Logger::out("I/O") << "is not a file" << std::endl;
        }
        mesh_gfx_.set_mesh(nullptr);

        mesh_.clear(false,false);
        
        if(GEO::CmdLine::get_arg_bool("single_precision")) {
            mesh_.vertices.set_single_precision();
        }
        
        MeshIOFlags flags;
        if(CmdLine::get_arg_bool("attributes")) {
            flags.set_attribute(MESH_FACET_REGION);
            flags.set_attribute(MESH_CELL_REGION);            
        } 
        if(!mesh_load(filename, mesh_, flags)) {
            return false;
        }

        if(
            FileSystem::extension(filename) == "obj6" ||
            FileSystem::extension(filename) == "tet6"
        ) {
            Logger::out("Vorpaview")
                << "Displaying mesh animation." << std::endl;

            glup_viewer_enable(GLUP_VIEWER_IDLE_REDRAW);
            
            mesh_gfx_.set_animate(true);
            double xyzmin[3];
            double xyzmax[3];
            get_bbox(mesh_, xyzmin, xyzmax, true);
            glup_viewer_set_region_of_interest(
                float(xyzmin[0]), float(xyzmin[1]), float(xyzmin[2]),
                float(xyzmax[0]), float(xyzmax[1]), float(xyzmax[2])
            );
        } else {
            mesh_gfx_.set_animate(false);            
            mesh_.vertices.set_dimension(3);
            double xyzmin[3];
            double xyzmax[3];
            get_bbox(mesh_, xyzmin, xyzmax, false);
            glup_viewer_set_region_of_interest(
                float(xyzmin[0]), float(xyzmin[1]), float(xyzmin[2]),
                float(xyzmax[0]), float(xyzmax[1]), float(xyzmax[2])
            );
        }

        show_vertices_ = (mesh_.facets.nb() == 0);
        mesh_gfx_.set_mesh(&mesh_);

	current_file_ = filename;
        return true;
    }

    bool SimpleMeshApplication::save(const std::string& filename) {
        MeshIOFlags flags;
        if(CmdLine::get_arg_bool("attributes")) {
            flags.set_attribute(MESH_FACET_REGION);
            flags.set_attribute(MESH_CELL_REGION);            
        }
	if(mesh_save(mesh_, filename, flags)) {
	    current_file_ = filename;
	    return true;
	}
	return false;
    }
    
    void SimpleMeshApplication::get_bbox(
        const Mesh& M_in, double* xyzmin, double* xyzmax, bool animate
    ) {
        geo_assert(M_in.vertices.dimension() >= index_t(animate ? 6 : 3));
        for(index_t c = 0; c < 3; c++) {
            xyzmin[c] = Numeric::max_float64();
            xyzmax[c] = Numeric::min_float64();
        }

        for(index_t v = 0; v < M_in.vertices.nb(); ++v) {
            if(M_in.vertices.single_precision()) {
                const float* p = M_in.vertices.single_precision_point_ptr(v);
                for(coord_index_t c = 0; c < 3; ++c) {
                    xyzmin[c] = std::min(xyzmin[c], double(p[c]));
                    xyzmax[c] = std::max(xyzmax[c], double(p[c]));
                    if(animate) {
                        xyzmin[c] = std::min(xyzmin[c], double(p[c+3]));
                        xyzmax[c] = std::max(xyzmax[c], double(p[c+3]));
                    }
                }
            } else {
                const double* p = M_in.vertices.point_ptr(v);
                for(coord_index_t c = 0; c < 3; ++c) {
                    xyzmin[c] = std::min(xyzmin[c], p[c]);
                    xyzmax[c] = std::max(xyzmax[c], p[c]);
                    if(animate) {
                        xyzmin[c] = std::min(xyzmin[c], p[c+3]);
                        xyzmax[c] = std::max(xyzmax[c], p[c+3]);
                    }
                }
            }
        }
    }
}

