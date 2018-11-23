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

#ifndef GEOGRAM_GFX_GLUP_VIEWER_GLUP_VIEWER_GUI
#define GEOGRAM_GFX_GLUP_VIEWER_GLUP_VIEWER_GUI

#include <geogram_gfx/basic/common.h>
#include <geogram_gfx/mesh/mesh_gfx.h>

#include <geogram/mesh/mesh.h>

#include <geogram/basic/counted.h>
#include <geogram/basic/smart_pointer.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/string.h>

#include <geogram_gfx/ImGui_ext/imgui_ext.h>
#include <geogram_gfx/ImGui_ext/file_dialog.h>

#include <geogram_gfx/third_party/ImGuiColorTextEdit/TextEditor.h>

struct lua_State;

/**
 * \file geogram_gfx/glup_viewer/glup_viewer_gui.h
 * \brief Utilities and C++ classes for easily creating 
 *  applications with glup_viewer and ImGui.
 */

namespace GEO {

    /*****************************************************************/

    /**
     * \brief Implementation of the GLUP viewer's status bar, that
     *  displays the progress bar.
     */
    class GEOGRAM_GFX_API StatusBar : public GEO::ProgressClient {
    public:

        /**
         * \brief StatusBar constructor.
         */
        StatusBar();
        
        /**
         * \copydoc GEO::ProgressClient::begin()
         */
        virtual void begin();
        
        /**
         * \copydoc GEO::ProgressClient::progress()
         */
        virtual void progress(GEO::index_t step, GEO::index_t percent);
        
        /**
         * \copydoc GEO::ProgressClient::end()
         */
        virtual void end(bool canceled);

        /**
         * \brief Draws the status bar and handles the GUI.
         */
        void draw();

        /**
         * \brief Tests whether this status bar should be displayed.
         * \details The status bar is displayed whenever there is an
         *  active progress bar.
         */
        bool active() const {
            return (nb_active_ > 0);
        }

	/**
	 * \brief Redraws the GUI.
	 */
	virtual void update();
	
      private:
        bool progress_;
        index_t step_;
        index_t percent_;
        bool canceled_;
        index_t nb_active_;
    };
    
    /*****************************************************************/

    /**
     * \brief Implementation of the GLUP viewer's console, where
     *  Logger message are displayed.
     * \details Inspired from ImGui AppLog example.
     */
    class GEOGRAM_GFX_API Console : public GEO::LoggerClient {
    
    public:

        /**
         * \brief Console constructor.
         * \param[in] visible_flag an optional pointer to application's
         *  variable that controls the visibility of this Console.
         */
        Console(bool* visible_flag = nullptr);
        
        /**
         * \copydoc GEO::LoggerClient::div()
         */
        virtual void div(const std::string& value);
        
        /**
         * \copydoc GEO::LoggerClient::out()
         */
        virtual void out(const std::string& value);
        
        /**
         * \copydoc GEO::LoggerClient::warn()
         */
        virtual void warn(const std::string& value);
        
        /**
         * \copydoc GEO::LoggerClient::err()
         */
        virtual void err(const std::string& value);
        
        /**
         * \copydoc GEO::LoggerClient::status()
         */
        virtual void status(const std::string& value);

        /**
         * \brief Clears the contents of the console.
         */
        void clear();

        /**
         * \brief Displays a formatted string to the console.
         */
        virtual void printf(const char* fmt, ...) /* IM_FMTARGS(1) */;

        /**
         * \brief Draws the console and handles the gui.
         * \param[in] visible an optional pointer to a visibility
         *  flag, controlled by a close button if different from nullptr.
	 * \param[in] with_window if true, then creates a new window
	 *  using imgui::Begin() / imgui::End(), else caller is responsible
	 *  for doing that.
         */
        virtual void draw(bool* visible=nullptr, bool with_window=true);

	int TextEditCallback(ImGuiTextEditCallbackData* data);

	void show() {
	    *visible_flag_ = true;
	}

	void hide() {
	    *visible_flag_ = false;
	}

	typedef void (*CompletionCallback)(
	    Console* console,
	    const std::string& line, index_t startw, index_t endw,
	    const std::string& cmpword, std::vector<std::string>& matches
	);

	void set_completion_callback(CompletionCallback CB) {
	    completion_callback_ = CB;
	}

	typedef void (*HistoryCallback)(
	    Console* console,
	    index_t index,
	    std::string& command
	);

	void set_history_callback(HistoryCallback CB) {
	    history_callback_ = CB;
	}

	void set_history_size(index_t n) {
	    if(n != max_history_index_) {
		history_index_  = n;
	    }
	    max_history_index_ = n;
	}

	void set_fixed_layout(bool x) {
	    fixed_layout_ = x;
	}

	void set_console_font(ImFont* font) {
	    console_font_ = font;
	}
	
      protected:
	/**
	 * \brief This function is called whenever an error is
	 *  displayed using err()
	 * \details Base implementation does nothing. This function
	 *  is meant to be overloaded in derived classes.
	 * \param[in] err the error message sent to err()
	 */
	virtual void notify_error(const std::string& err);

	virtual bool exec_command(const char* command);
	
	/**
	 * \brief Redraws the GUI.
	 */
	virtual void update();
	
        ImGuiTextBuffer buf_;
        ImGuiTextFilter filter_;
        /** \brief Index to lines offset */
        ImVector<int>      line_offsets_;   
        bool               scroll_to_bottom_;
        bool*              visible_flag_;
	char               input_buf_[geo_imgui_string_length];
	CompletionCallback completion_callback_;
	HistoryCallback    history_callback_;
	index_t            history_index_;
	index_t            max_history_index_;
	bool               fixed_layout_;
	ImFont*            console_font_;
    };
    
    /*****************************************************************/

    class Application;
    
    /**
     * \brief Abstract class for calling functions or 
     *  calling member functions.
     * \details Used internally by Command.
     */
    class GEOGRAM_GFX_API CommandInvoker : public Counted {
    public:
        /**
         * \brief CommandInvoker constructor.
         */
        CommandInvoker();

        /**
         * \brief CommandInvoker destructor.
         */
        virtual ~CommandInvoker();

        /**
         * \brief Invokes the target function.
         */
        virtual void invoke() = 0;

        /**
         * \brief Creates the arguments in the target
         *  command.
         * \details This function is used when client code
         *  did not provide a function prototype to the
         *  constructor of Command. 
         */
        virtual void auto_create_args() = 0;
    };

    /**
     * \brief Automatic reference-counted pointer to a CommandInvoker.
     * \details Used internally by Command.
     */
    typedef SmartPointer<CommandInvoker> CommandInvoker_var;
    
    /*****************************************************************/

    /**
     * \brief Manages the GUI of a command with ImGUI.
     * \details Client code will only need to use set_current()
     */
    class GEOGRAM_GFX_API Command : public Counted {
    public:

        /**
         * \brief Binds the current command to a function.
         * \details This makes the command dialog box display 
         *  parameters that correspond to the arguments of the function,
         *  and whenever the 'apply' button is pushed, the function is
         *  invoked with the arguments.
         *  Example:
         *  \code
         *     void my_command_impl(float x, float y, bool normalize) {
         *     ... do something
         *     }
         *     ...
         *     ... 
         *     if(ImGui::MenuItem("my command")) {
         *         GEO::Command::set_current(
         *            "void my_command_impl(float x, float y, bool normalize)",
         *            &my_command_impl
         *         )
         *     }
         *  \endcode
         *  The first argument (the string with the function prototype) is 
         *  necessary to retrieve the names of the parameters. In addition,
         *  default values and tooltips may be specified, as follows:
         *  \code
         *         GEO::Command::set_current(
         *            "void my_command_impl(                        "
         *            "   float x=0 [this is the x coordinate],     "
         *            "   float y=0 [this is the y coordinate],     "
         *            "   bool normalize=true [normalize the vector]"
         *            ") [does something with a vector]             ",
         *            &my_command_impl
         *         )
         *  \endcode
         *  Each text in the square brackets corresponds to a tooltip attached
         *  to an argument. There can be also one for the function. The default
         *  values are those that are obtained at initialization, or those that
         *  are set when the 'default' button is pushed.
         * \tparam FPTR function pointer type
         * \param[in] prototype a string with the prototype of the function,
         *  as written in C++. In addition, the function and each parameter 
         *  can be documented in square brackets.
         * \param[in] tfun the function pointer.
         */
        template<class FPTR> static void set_current(
            const std::string& prototype, FPTR tfun
        );
        
        /**
         * \brief Binds the current command to a member function
         *  of an object.
         * \details This makes the command dialog box display 
         *  parameters that correspond to the arguments of the function,
         *  and whenever the 'apply' button is pushed, the function is
         *  invoked with the arguments.
         *  Example:
         *  \code
         *     class MyCommands {
         *        void my_command_impl(float x, float y, bool normalize) {
         *        ... do something
         *        }
         *     };
         *     MyCommands my_commands;
         *     ...
         *     ... 
         *     if(ImGui::MenuItem("my command")) {
         *         GEO::Command::set_current(
         *            "void my_command_impl(float x, float y, bool normalize)",
         *            &my_commands, &MyCommands::my_command_impl
         *         )
         *     }
         *  \endcode
         *  The first argument (the string with the function prototype) is 
         *  necessary to retrieve the names of the parameters. In addition,
         *  default values and tooltips may be specified, as follows:
         *  \code
         *         GEO::Command::set_current(
         *            "void my_command_impl(                        "
         *            "   float x=0 [this is the x coordinate],     "
         *            "   float y=0 [this is the y coordinate],     "
         *            "   bool normalize=true [normalize the vector]"
         *            ") [does something with a vector]             ",
         *            &my_commands, &MyCommands::my_command_impl
         *         )
         *  \endcode
         *  Each text in the square brackets corresponds to a tooltip attached
         *  to an argument. There can be also one for the function. The default
         *  values are those that are obtained at initialization, or those that
         *  are set when the 'default' button is pushed.
         * \tparam T object class
         * \tparam TFPTR function pointer type, should be a member function
         *  of class T
         * \param[in] prototype a string with the prototype of the function,
         *  as written in C++. In addition, the function and each parameter 
         *  can be documented in square brackets.
         * \param[in] target a pointer to the object, of class T
         * \param[in] tfun the pointer to the member function
         */
        template<class T, class TFPTR> static void set_current(
            const std::string& prototype, T* target, TFPTR tfun
        );

        /**
         * \brief Flushes the potentially queued command invokation.
         * \details When the user pushes the 'apply' button, the command
         *  is not invoked immediatly, because we are still in the ImGUI
         *  handling function. This function is called by the framework
         *  at the end of the frame, when the ImGUI handler is already
         *  finished. It can potentially re-trigger frame rendering operations,
         *  through the Logger and ProgressLogger. Without this mechanism,
         *  it would nest two ImGUI handlers, which is not allowed.
         */
        static void flush_queue();
        
        /**
         * \brief Command constructor.
         * \param[in] prototype a const reference to a string with
         *  the prototype of the function that implements the callback, 
         *  as declared in the C++ sources.
         * \note Regular client code should not need to use this function.
         */
        Command(const std::string& prototype);


        /**
         * \brief Gets the name of this command.
         * \return a const reference to the name of this command.
         */
        const std::string& name() const {
            return name_;
        }
        
        /**
         * \brief Sets the invoker.
         * \details The invoker is used internally to transmit the stored
         *  arguments to the parameters of a function. 
         * \param[in] invoker a pointer to the CommandInvoker. Ownership
         *  is transferred to this Command.
         * \note Regular client code should not need to use this function.
         */
        void set_invoker(CommandInvoker* invoker) {
            invoker_ = invoker;
            if(auto_create_args_) {
                invoker_->auto_create_args();
                auto_create_args_ = false;
            }
        }
        
        /**
         * \brief Command destructor.
         */
        virtual ~Command();

        /**
         * \brief Tests whether this Command is visible.
         * \retval true if this Command is visible
         * \retval false otherwise
         */
        bool is_visible() const {
            return visible_;
        }

        /**
         * \brief Gets a pointer to the visibility flag of
         *   this command.
         * \return a pointer to the visibility flag
         */
        bool* is_visible_ptr() {
            return &visible_;
        }
        
        /**
         * \brief Displays and manages the GUI of this 
         *  Command.
         * \note Regular client code should not need to use this function.
         */
        virtual void draw();

        /**
         * \brief Restores default parameter values for
         *  all parameters.
         * \details This is the function that is called when the user pushes
         *  the 'default' button.
         * \note Regular client code should not need to use this function.
         */
        virtual void reset_factory_settings();

        /**
         * \brief Gets the value of the parameters and
         *  does the task. 
         * \details This is the function that is called when the user pushes
         *  the 'apply' button. It does not invoke the command immediatly,
         *  the command invocation is queued, and executed later by
         *  Command::flush_queue(), once we are no longer in the ImGui 
         *  handler (else we would have two nested ImGui handlers, which is
         *  not correct).
         * \note Regular client code should not need to use this function.
         */
        virtual void apply() ;

        /**
         * \brief Gets the current command.
         * \return a pointer to the current command
         * \note Regular client code should not need to use this function.
         */
        static Command* current() {
            return current_;
        }

        /**
         * \brief Resets the current command.
         */
        static void reset_current() {
            current_.reset();
        }
        
        /**
         * \brief Sets the current command.
         * \param[in] command a pointer to the command
         *  to be set as current
         * \note Regular client code should not need to use this function.
         */
        static void set_current(Command* command) {
            current_ = command;
            command->visible_ = true;
        }

        /**
         * \brief Gets the value of a boolean argument by index.
         * \details If the index is out of range, or if the
         *  argument is not of the correct type, then an assertion
         *  failure is triggered.
         * \param[in] i the index of the argument
         * \return the value of the argument
         * \note Regular client code should not need to use this function.
         */
        bool bool_arg_by_index(index_t i) const {
            const Arg& arg = find_arg_by_index(i);
            geo_assert(arg.type == Arg::ARG_BOOL);
            return arg.val.bool_val;
        }

        /**
         * \brief Gets the value of an integer argument by index.
         * \details If the index is out of range, or if the
         *  argument is not of the correct type, then an assertion
         *  failure is triggered.
         * \param[in] i the index of the argument
         * \return the value of the argument
         * \note Regular client code should not need to use this function.
         */
        int int_arg_by_index(index_t i) const;

        /**
         * \brief Gets the value of an unsigned integer argument by index.
         * \details If the index is out of range, or if the
         *  argument is not of the correct type, then an assertion
         *  failure is triggered. If the stored value is negative, then
         *  it is clamped to 0, and a warning message is displayed on
         *  the console.
         * \param[in] i the index of the argument
         * \return the value of the argument
         * \note Regular client code should not need to use this function.
         */
        unsigned int uint_arg_by_index(index_t i) const; 

        /**
         * \brief Gets the value of a floating-point argument by index.
         * \details If the index is out of range, or if the
         *  argument is not of the correct type, then an assertion
         *  failure is triggered.
         * \param[in] i the index of the argument
         * \return the value of the argument
         * \note Regular client code should not need to use this function.
         */
        float float_arg_by_index(index_t i) const {
            const Arg& arg = find_arg_by_index(i);
            geo_assert(arg.type == Arg::ARG_FLOAT);
            return arg.val.float_val;
        }

        /**
         * \brief Gets the value of a floating-point argument by index
         *  and converts it to a double.
         * \details If the index is out of range, or if the
         *  argument is not of the correct type, then an assertion
         *  failure is triggered.
         * \param[in] i the index of the argument
         * \return the value of the argument
         * \note Regular client code should not need to use this function.
         */
        double double_arg_by_index(index_t i) const {
            return double(float_arg_by_index(i));
        }

        /**
         * \brief Gets the value of a string argument by index.
         * \details If the index is out of range, or if the
         *  argument is not of the correct type, then an assertion
         *  failure is triggered.
         * \param[in] i the index of the argument
         * \return the value of the argument
         * \note Regular client code should not need to use this function.
         */
        std::string string_arg_by_index(index_t i) const {
            const Arg& arg = find_arg_by_index(i);
            geo_assert(arg.type == Arg::ARG_STRING);
            return std::string(arg.val.string_val);
        }

        /**
         * \brief Gets the value of an argument by index.
         * \details This function is generic, and has several
         *  specializations for bool, int, unsigned int, float, double and
         *  std::string. For all other types, an assertion failure is 
         *  triggered.
         * \tparam T type of the argument
         * \param[in] i the index of the argument
         * \param[out] val a reference to the argument
         * \note Regular client code should not need to use this function.
         */
        template<class T> void get_arg_by_index(index_t i, T& val) {
            geo_argused(val);
            Logger::err("Cmd") << "Attempted to read argument #"
                               << i
                               << " to variable of unknown type"
                               << std::endl;
            geo_assert_not_reached;
        }

        /***************************************************************/

        /**
         * \brief Invokes a member function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam T class of the target object
         * \param[in] target the target object
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        template <class T> void invoke(
            T* target, void (T::*fptr)(void)
        ) {
            this->assert_nb_args_matches(0);
            if(target != nullptr && fptr != nullptr) {
                (*target.*fptr)();
            }
        }

        /**
         * \brief Invokes a member function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam T class of the target object
         * \tparam ARG0 type of the argument
         * \param[in] target the target object
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        template <
            class T,
            class ARG0
        > void invoke(
            T* target,
            void (T::*fptr)(ARG0)
        ) {
            this->assert_nb_args_matches(1);            
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            if(target != nullptr && fptr != nullptr) {            
                (*target.*fptr)(a0);
            }
        }

        /**
         * \brief Invokes a member function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam T class of the target object
         * \tparam ARG0 , ARG1 types of the arguments
         * \param[in] target the target object
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        template <
            class T,
            class ARG0, class ARG1
        > void invoke(
            T* target,
            void (T::*fptr)(ARG0,ARG1)
        ) {
            this->assert_nb_args_matches(2);                        
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            if(target != nullptr && fptr != nullptr) {
                (*target.*fptr)(a0,a1);
            }
        }

        /**
         * \brief Invokes a member function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam T class of the target object
         * \tparam ARG0 , ARG1 , ARG2 types of the arguments
         * \param[in] target the target object
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class T,
            class ARG0, class ARG1, class ARG2
        > void invoke(
            T* target,
            void (T::*fptr)(ARG0,ARG1,ARG2)
        ) {
            this->assert_nb_args_matches(3);                                    
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            if(target != nullptr && fptr != nullptr) {            
                (*target.*fptr)(a0,a1,a2);
            }
        }

        /**
         * \brief Invokes a member function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam T class of the target object
         * \tparam ARG0 , ARG1 , ARG2 , ARG3 types of the arguments
         * \param[in] target the target object
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class T,
            class ARG0, class ARG1, class ARG2, class ARG3
        > void invoke(
            T* target,
            void (T::*fptr)(ARG0,ARG1,ARG2,ARG3)
        ) {
            this->assert_nb_args_matches(4);            
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            ARG3 a3;
            this->get_arg_by_index(3,a3);
            if(target != nullptr && fptr != nullptr) {
                (*target.*fptr)(a0,a1,a2,a3);
            }
        }

        /**
         * \brief Invokes a member function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam T class of the target object
         * \tparam ARG0 , ARG4 types of the arguments
         * \param[in] target the target object
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class T,
            class ARG0, class ARG1, class ARG2, class ARG3,
            class ARG4
        > void invoke(
            T* target,
            void (T::*fptr)(ARG0,ARG1,ARG2,ARG3,ARG4)
        ) {
            this->assert_nb_args_matches(5);                        
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            ARG3 a3;
            this->get_arg_by_index(3,a3);
            ARG4 a4;
            this->get_arg_by_index(4,a4);
            if(target != nullptr && fptr != nullptr) {
                (*target.*fptr)(a0,a1,a2,a3,a4);
            }
        }

        /**
         * \brief Invokes a member function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam T class of the target object
         * \tparam ARG0 , ARG5 types of the arguments
         * \param[in] target the target object
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class T,
            class ARG0, class ARG1, class ARG2, class ARG3,
            class ARG4, class ARG5
        > void invoke(
            T* target,
            void (T::*fptr)(ARG0,ARG1,ARG2,ARG3,ARG4,ARG5)
        ) {
            this->assert_nb_args_matches(6);                                    
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            ARG3 a3;
            this->get_arg_by_index(3,a3);
            ARG4 a4;
            this->get_arg_by_index(4,a4);
            ARG5 a5;
            this->get_arg_by_index(5,a5);
            if(target != nullptr && fptr != nullptr) {
                (*target.*fptr)(a0,a1,a2,a3,a4,a5);
            }
        }

        /**
         * \brief Invokes a member function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam T class of the target object
         * \tparam ARG0 , ARG6 types of the arguments
         * \param[in] target the target object
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class T,
            class ARG0, class ARG1, class ARG2, class ARG3,
            class ARG4, class ARG5, class ARG6
        > void invoke(
            T* target,
            void (T::*fptr)(ARG0,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6)
        ) {
            this->assert_nb_args_matches(7);            
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            ARG3 a3;
            this->get_arg_by_index(3,a3);
            ARG4 a4;
            this->get_arg_by_index(4,a4);
            ARG5 a5;
            this->get_arg_by_index(5,a5);
            ARG6 a6;
            this->get_arg_by_index(6,a6);
            if(target != nullptr && fptr != nullptr) {
                (*target.*fptr)(a0,a1,a2,a3,a4,a5,a6);
            }
        }

        /**
         * \brief Invokes a member function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam T class of the target object
         * \tparam ARG0 , ARG7 types of the arguments
         * \param[in] target the target object
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class T,
            class ARG0, class ARG1, class ARG2, class ARG3,
            class ARG4, class ARG5, class ARG6, class ARG7
        > void invoke(
            T* target,
            void (T::*fptr)(ARG0,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7)
        ) {
            this->assert_nb_args_matches(8);                        
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            ARG3 a3;
            this->get_arg_by_index(3,a3);
            ARG4 a4;
            this->get_arg_by_index(4,a4);
            ARG5 a5;
            this->get_arg_by_index(5,a5);
            ARG6 a6;
            this->get_arg_by_index(6,a6);
            ARG7 a7;
            this->get_arg_by_index(7,a7);
            if(target != nullptr && fptr != nullptr) {
                (*target.*fptr)(a0,a1,a2,a3,a4,a5,a6,a7);
            }
        }

        /**************************************************************/

        /**
         * \brief Invokes a function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        void invoke(
            void (*fptr)(void)
        ) {
            this->assert_nb_args_matches(0);                                    
            if(fptr != nullptr) {
                (*fptr)();
            }
        }

        /**
         * \brief Invokes a function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam ARG0 type of the argument
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class ARG0
        > void invoke(
            void (*fptr)(ARG0)
        ) {
            this->assert_nb_args_matches(1);            
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            if(fptr != nullptr) {
                (*fptr)(a0);
            }
        }

        /**
         * \brief Invokes a function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam ARG0 , ARG1 types of the arguments
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class ARG0, class ARG1
        > void invoke(
            void (*fptr)(ARG0,ARG1)
        ) {
            this->assert_nb_args_matches(2);                        
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            if(fptr != nullptr) {
                (*fptr)(a0,a1);
            }
        }

        /**
         * \brief Invokes a function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam ARG0 , ARG1 , ARG2 types of the arguments
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class ARG0, class ARG1, class ARG2
        > void invoke(
            void (*fptr)(ARG0,ARG1,ARG2)
        ) {
            this->assert_nb_args_matches(3);                                    
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            if(fptr != nullptr) {
                (*fptr)(a0,a1,a2);
            }
        }

        /**
         * \brief Invokes a function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam ARG0 , ARG1 , ARG2 , ARG3 types of the arguments
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class ARG0, class ARG1, class ARG2, class ARG3
        > void invoke(
            void (*fptr)(ARG0,ARG1,ARG2,ARG3)
        ) {
            this->assert_nb_args_matches(4);             
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            ARG3 a3;
            this->get_arg_by_index(3,a3);
            if(fptr != nullptr) {
                (*fptr)(a0,a1,a2,a3);
            }
        }

        /**
         * \brief Invokes a function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam ARG0 , ARG1 , ARG2 , ARG3 , ARG4 types of the arguments
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class ARG0, class ARG1, class ARG2, class ARG3,
            class ARG4
        > void invoke(
            void (*fptr)(ARG0,ARG1,ARG2,ARG3,ARG4)
        ) {
            this->assert_nb_args_matches(5);                         
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            ARG3 a3;
            this->get_arg_by_index(3,a3);
            ARG4 a4;
            this->get_arg_by_index(4,a4);
            if(fptr != nullptr) {
                (*fptr)(a0,a1,a2,a3,a4);
            }
        }

        /**
         * \brief Invokes a function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam ARG0 , ARG5 types of the arguments
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class ARG0, class ARG1, class ARG2, class ARG3,
            class ARG4, class ARG5
        > void invoke(
            void (*fptr)(ARG0,ARG1,ARG2,ARG3,ARG4,ARG5)
        ) {
            this->assert_nb_args_matches(6);             
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            ARG3 a3;
            this->get_arg_by_index(3,a3);
            ARG4 a4;
            this->get_arg_by_index(4,a4);
            ARG5 a5;
            this->get_arg_by_index(5,a5);
            if(fptr != nullptr) {
                (*fptr)(a0,a1,a2,a3,a4,a5);
            }
        }

        /**
         * \brief Invokes a function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam ARG0 , ARG6 types of the arguments
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class ARG0, class ARG1, class ARG2, class ARG3,
            class ARG4, class ARG5, class ARG6
        > void invoke(
            void (*fptr)(ARG0,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6)
        ) {
            this->assert_nb_args_matches(7);                         
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            ARG3 a3;
            this->get_arg_by_index(3,a3);
            ARG4 a4;
            this->get_arg_by_index(4,a4);
            ARG5 a5;
            this->get_arg_by_index(5,a5);
            ARG6 a6;
            this->get_arg_by_index(6,a6);
            if(fptr != nullptr) {
                (*fptr)(a0,a1,a2,a3,a4,a5,a6);
            }
        }

        /**
         * \brief Invokes a function with the stored arguments.
         * \details If the type or number of stored arguments do not match
         *  the function pointer, then an assertion failure is triggered.
         * \tparam ARG0 , ARG7 types of the arguments
         * \param[in] fptr the pointer to the member function to be called.
         * \note Used by internal CommandInvoker mechanism.
         */
        
        template <
            class ARG0, class ARG1, class ARG2, class ARG3,
            class ARG4, class ARG5, class ARG6, class ARG7
        > void invoke(
            void (*fptr)(ARG0,ARG1,ARG2,ARG3,ARG4,ARG5,ARG6,ARG7)
        ) {
            this->assert_nb_args_matches(8);             
            ARG0 a0;
            this->get_arg_by_index(0,a0);
            ARG1 a1;
            this->get_arg_by_index(1,a1);
            ARG2 a2;
            this->get_arg_by_index(2,a2);
            ARG3 a3;
            this->get_arg_by_index(3,a3);
            ARG4 a4;
            this->get_arg_by_index(4,a4);
            ARG5 a5;
            this->get_arg_by_index(5,a5);
            ARG6 a6;
            this->get_arg_by_index(6,a6);
            ARG7 a7;
            this->get_arg_by_index(7,a7);
            if(fptr != nullptr) {
                (*fptr)(a0,a1,a2,a3,a4,a5,a6,a7);
            }
        }

        /**************************************************************/
        
    protected:

        /**
         * \brief Tests whether the number of declared arguments
         *  matches a specified number.
         * \details If the number of arguments differs from the expected
         *  number, then an assertion failure is triggered. 
         *  When auto_create_args_ is set, number of stored arguments should
         *  be 0.
         * \param[in] nb expected number of arguments.
         */
        void assert_nb_args_matches(index_t nb) {
            geo_assert(
                (auto_create_args_ && args_.size() == 0) ||
                (args_.size() == nb)
            );
        }

        
        /**
         * \brief Adds a parameter to this command
         * \tparam T type of the parameter, deduced from
         *  \p default_val
         * \param[in] name name of the parameter
         * \param[in] default_val default value of the parameter
         * \param[in] help optionnal text, displayed in
         *  a tooltip
         */
        template<class T> void add_arg(
            const std::string& name, const T& default_val,
            const std::string& help = ""
        ) {
            args_.push_back(Arg(name, default_val, help));
        }

        /**
         * \brief Creates an argument at a given index
         * \details Used by the auto_create_args mechanism, used
         *  when client code did not provide any function prototype
         * \tparam T type of the parameter, deduced from
         *  \p default_val
         * \param[in] default_val default value of the parameter
         */
        template<class T> void create_arg(
            index_t i, const T& default_val
        ) {
            if(i >= args_.size()) {
                args_.resize(i+1);
            }
            args_[i] = Arg("arg " + String::to_string(i), default_val);
        }

	static void set_queued(Command* command) {
	    queued_ = command;
	}
	
    private:
        
        /**
         * \brief Internal representation of an argument
         *   value.
         */
        struct GEOGRAM_GFX_API ArgVal {
            /**
             * \brief Resets all values stored to zero.
             */
            void clear();

            /**
             * \brief Argval constructor.
             */
            ArgVal() {
                clear();
            }

            /**
             * \brief Argval copy constructor.
             * \param[in] rhs a const reference to the ArgVal
             *  to be copied
             */
            ArgVal(const ArgVal& rhs);

            /**
             * \brief Argval assignment operator.
             * \param[in] rhs a const reference to the ArgVal
             *  to be copied
             * \return a reference to this ArgVal after assignment
             */
            ArgVal& operator=(const ArgVal& rhs);
            
            bool bool_val;
            int int_val;
            float float_val;
            char string_val[64];
        };

        /**
         * \brief Internal representation of an argument.
         * \details An argument has a type, a name, a
         *  default value, and an optionnal tooltip.
         */
        struct GEOGRAM_GFX_API Arg {

            /**
             * \brief Arg default constructor.
             */
            Arg();
            
            /**
             * \brief Arg constructor from bool.
             * \param[in] name_in the name of the argument
             * \param[in] x the default value
             * \param[in] help_in an optional text, displayed
             *  as a tooltip
             */
            Arg(
                const std::string& name_in, bool x,
                const std::string& help_in=""
            );

            /**
             * \brief Arg constructor from int.
             * \param[in] name_in the name of the argument
             * \param[in] x the default value
             * \param[in] help_in an optional text, displayed
             *  as a tooltip
             */
            Arg(
                const std::string& name_in, int x,
                const std::string& help_in=""
            );

            /**
             * \brief Arg constructor from unsigned int.
             * \details Stored internally as (signed) int
             * \param[in] name_in the name of the argument
             * \param[in] x the default value
             * \param[in] help_in an optional text, displayed
             *  as a tooltip
             */
            Arg(
                const std::string& name_in, unsigned int x,
                const std::string& help_in=""
            );

            /**
             * \brief Arg constructor from float.
             * \param[in] name_in the name of the argument
             * \param[in] x the default value
             * \param[in] help_in an optional text, displayed
             *  as a tooltip
             */
            Arg(
                const std::string& name_in, float x,
                const std::string& help_in=""                
            );

            /**
             * \brief Arg constructor from double.
             * \details Stored internally as a float
             * \param[in] name_in the name of the argument
             * \param[in] x the default value
             * \param[in] help_in an optional text, displayed
             *  as a tooltip
             */
            Arg(
                const std::string& name_in, double x,
                const std::string& help_in=""                                
            );

            /**
             * \brief Arg constructor from string.
             * \details String length is limited to 64 characters.
             * \param[in] name_in the name of the argument
             * \param[in] x the default value
             * \param[in] help_in an optional text, displayed
             *  as a tooltip
             */
            Arg(
                const std::string& name_in, const std::string& x,
                const std::string& help_in=""
            );

            /**
             * \brief Displays and manages the GUI for this Arg.
             */
            void draw();
            
            enum { ARG_BOOL, ARG_INT, ARG_UINT, ARG_FLOAT, ARG_STRING } type;
            std::string name;
            std::string help;
            ArgVal val;
            ArgVal default_val;
        };

        /**
         * \brief Finds an argument value by name.
         * \details Triggers an assertion failure if no such
         *  argument exists.
         * \param[in] name name of the argument
         * \return a const reference to the Arg
         */
        const Arg& find_arg(const std::string& name) const {
            for(index_t i=0; i<args_.size(); ++i) {
                if(args_[i].name == name) {
                    return args_[i];
                }
            }
            geo_assert_not_reached;
        }

        /**
         * \brief Gets an Arg by its index.
         * \param[in] i the index of the Arg
         * \return a const reference to the Arg
         */
        const Arg& find_arg_by_index(index_t i) const {
            geo_assert(i < args_.size());
            return args_[i];
        }
        
    private:
        std::string name_;
        std::string help_;
        vector<Arg> args_;
        CommandInvoker_var invoker_;
        bool visible_;
        /**
         * \brief If no prototype was specified, then
         *  arguments are automatically created.
         */
        bool auto_create_args_; 
        static SmartPointer<Command> current_;
        static SmartPointer<Command> queued_;
    };

/***********************************************************************/

    /**
     * \copydoc Command::get_arg_by_index()
     */
    template<> inline void Command::get_arg_by_index(
        index_t i, bool& val
    ) {
        if(auto_create_args_) {
            val = false;
            this->create_arg(i, val);
        } else {
            val = this->bool_arg_by_index(i);
        }
    }

    /**
     * \copydoc Command::get_arg_by_index()
     */
    template<> inline void Command::get_arg_by_index(
        index_t i, int& val
    ) {
        if(auto_create_args_) {
            val = 0;
            this->create_arg(i, val);
        } else {        
            val = this->int_arg_by_index(i);
        }
    }

    /**
     * \copydoc Command::get_arg_by_index()
     */
    template<> inline void Command::get_arg_by_index(
        index_t i, unsigned int& val
    ) {
        if(auto_create_args_) {
            val = 0;
            this->create_arg(i, val);
        } else {
            val = this->uint_arg_by_index(i);
        }
    }

    /**
     * \copydoc Command::get_arg_by_index()
     */
    template<> inline void Command::get_arg_by_index(
        index_t i, float& val
    ) {
        if(auto_create_args_) {
            val = 0.0f;
            this->create_arg(i, val);
        } else {
            val = this->float_arg_by_index(i);
        }
    }

    /**
     * \copydoc Command::get_arg_by_index()
     */
    template<> inline void Command::get_arg_by_index(
        index_t i, double& val
    ) {
        if(auto_create_args_) {
            val = 0.0;
            this->create_arg(i, val);
        } else {        
            val = this->double_arg_by_index(i);
        }
    }

    /**
     * \copydoc Command::get_arg_by_index()
     */
    template<> inline void Command::get_arg_by_index(
        index_t i, std::string& val
    ) {
        if(auto_create_args_) {
            val = "";
            this->create_arg(i, val);
        } else {                
            val = this->string_arg_by_index(i);
        }
    }

    /*****************************************************************/

    /**
     * \brief An implementation of CommandInvoker that calls 
     *  a function.
     * \tparam FPTR function pointer type for the function to be called
     */
    template <class FPTR>
    class FunctionCommandInvoker : public CommandInvoker {
    public:

        /**
         * \brief FunctionCommandInvoker constructor.
         * \param[in] command a pointer to the Command object
         * \param[in] fun the function pointer
         */
        FunctionCommandInvoker(
            Command* command,
            FPTR fun
        ) :
            command_(command),
            fun_(fun) {
        }

        /**
         * \copydoc CommandInvoker::invoke()
         */
        virtual void invoke() {
            command_->invoke(fun_);
        }

        /**
         * \copydoc CommandInvoker::auto_create_args()
         */
        virtual void auto_create_args() {
            command_->invoke(FPTR(nullptr));            
        }
        
    private:
        Command* command_;
        FPTR fun_; 
    };

    /*****************************************************************/

    /**
     * \brief An implementation of CommandInvoker that calls 
     *  a member function of an object.
     * \tparam T class of the object
     * \tparam TFPTR function pointer type for the function to be called
     */
    
    template <class T, class TFPTR>
    class MemberFunctionCommandInvoker : public CommandInvoker {
    public:

        /**
         * \brief MemberFunctionCommandInvoker constructor.
         * \param[in] command a pointer to the Command object
         * \param[in] target a pointer to the object
         * \param[in] target_fun the member function pointer
         */
        
        MemberFunctionCommandInvoker(
            Command* command,
            T* target,
            TFPTR target_fun
        ) :
            command_(command),
            target_(target),
            target_fun_(target_fun) {
        }

        /**
         * \copydoc CommandInvoker::invoke()
         */
        
        virtual void invoke() {
            command_->invoke(target_, target_fun_);
        }

        /**
         * \copydoc CommandInvoker::auto_create_args()
         */
        
        virtual void auto_create_args() {
            command_->invoke((T*)(nullptr), (TFPTR)(nullptr));            
        }

        
    private:
        Command* command_;
        T* target_;
        TFPTR target_fun_; 
    };

    /*****************************************************************/
    
    template<class FPTR> inline void Command::set_current(
        const std::string& prototype,
        FPTR fun
    ) {
        set_current(new Command(prototype));
        current()->set_invoker(
            new FunctionCommandInvoker<FPTR>(current(), fun)
        );
    }

    /*****************************************************************/    
    
    template<class T, class TFPTR> inline void Command::set_current(
        const std::string& prototype,
        T* target,
        TFPTR tfun
    ) {
        set_current(new Command(prototype));
        current()->set_invoker(
            new MemberFunctionCommandInvoker<T, TFPTR>(current(), target, tfun)
        );
    }
    
    /*****************************************************************/    

    /**
     * \brief A simple text editor
     */
    class GEOGRAM_GFX_API TextEditor {
    public:
	TextEditor(bool* visible);
	void draw();
	std::string text() const;
	void load(const std::string& filename);
	void save(const std::string& filename);
	void clear();
	void load_data(const char* data);
	
    private:
	::TextEditor impl_;
	bool* visible_;
    };
    
    /*****************************************************************/    

    
    /**
     * \brief Base class for glup_viewer applications with a gui.
     */
    class GEOGRAM_GFX_API Application {
    public:
        /**
         * \brief Application constructor.
         * \param[in] argc , argv command line arguments copied from main()
         * \param[in] usage the usage string
	 * \param[in] lua_state an optional pointer to a LUA state or nullptr. If
	 *  nullptr, then a LUA state is created.
         * \see CmdLine::parse()
         */
        Application(
	    int argc, char** argv, const std::string& usage,
	    lua_State* lua_state = nullptr
	);

        /**
         * \brief Application destructor.
         */
        virtual ~Application();

        /**
         * \brief Starts the application.
         */
        void start();


	/**
	 * \brief Stops the application.
	 */
	void quit();
	
        /**
         * \brief Gets the instance.
         * \return A pointer to the instance.
         */
        static Application* instance() {
            return instance_;
        }

        /**
         * \brief Saves the current content to a file.
         * \details Baseclass implementation does nothing. Derived classes
         *  may overload this function.
         * \retval true if the file could be sucessfully saved.
         * \retval false otherwise
         */
        virtual bool save(const std::string& filename);
        
        /**
         * \brief Loads a file.
         * \details Baseclass implementation does nothing. Derived classes
         *  may overload this function.
         * \retval true if the file could be sucessfully loaded
         * \retval false otherwise
         */
        virtual bool load(const std::string& filename);

        /**
         * \brief Tests whether a file can be loaded.
         * \details This function can be used to filter the files displayed
         *  in the "Load..." menu. Baseclass implementation always return true.
         *  Derived classes may overload it and return false for files with
         *  unknown extensions.
         */  
        virtual bool can_load(const std::string& filename);

        /**
         * \brief Gets the list of supported file extensions for reading.
         * \details This function may be olverloaded by derived class. Base
         *  class implementation returns "". If this function returns "", then
         *  no "Load..." option is displayed in the "File" menu. 
         * \return The semi-colon separated list of supported file extensions,
         *  or "*" if all file extensions are supported.
         */
        virtual std::string supported_read_file_extensions(); 

        /**
         * \brief Gets the list of supported file extensions for writing.
         * \details This function may be olverloaded by derived class. Base
         *  class implementation returns "". If this function returns "", then
         *  no "Save..." option is displayed in the "File" menu.  
         *   If it returns a colon-separated list of extensions, then the
         *  "Save..." option displays a list of possible file names for each
         *  supported extension.
         * \return The semi-colon separated list of supported file extensions,
         *  or "*" if all file extensions are supported.
         */
        virtual std::string supported_write_file_extensions(); 

        /**
         * \brief Gets the scaling applied to all dimensions.
         * \details This function is used by retina displays to ensure that
         *  GUI elements remain visible.
         * \return The scaling used for fonts and all sizes.
         */
        float scaling() const {
            return scaling_;
        }

        /**
         * \brief Sets the scaling applied to all dimensions.
         * \details This function is used by retina displays to ensure that
         *  GUI elements remain visible.
         * \param x The scaling used for fonts and all sizes.
         */
	void set_scaling(float x) {
	    scaling_ = x;
	}


	bool console_visible() const {
	    return console_visible_;
	}

	void set_console_visible(bool x) {
	    console_visible_ = x;
	}

	virtual bool exec_command(const char* command);

	Console* console() {
	    return console_;
	}
	
	bool retina_mode() const {
	    return retina_mode_;
	}

	void set_retina_mode(bool x) {
	    retina_mode_ = x;
	    scaling_ = retina_mode_ ? 2.0f : 1.0f;
	}
	
	virtual bool on_key_pressed(const char* key);
	virtual bool on_key_released(const char* key);

	void set_lighting(bool x) {
	    lighting_ = x;
	}


	void set_background_color_1(float r, float g, float b) {
	    background_color_1_.x = r;
	    background_color_1_.y = g;
	    background_color_1_.z = b;
	    background_color_1_.w = 1.0;
	}

	void set_background_color_2(float r, float g, float b) {
	    background_color_2_.x = r;
	    background_color_2_.y = g;
	    background_color_2_.z = b;
	    background_color_1_.w = 1.0;	    
	}

	void set_background_color_1(const vec4f& value) {
	    background_color_1_ = value;
	}

	void set_background_color_2(const vec4f& value) {
	    background_color_2_ = value;
	}
	
	const vec4f& get_background_color_1() const {
	    return background_color_1_;
	}

	const vec4f& get_background_color_2() const {
	    return background_color_2_;	    
	}
	
    protected:

        /**
         * \brief Converts an OpenGL texture ID into an ImGUI texture ID.
         * \param[in] gl_texture_id_in the OpenGL texture ID
         * \return the corresponding ImGUI texture ID
         */
        ImTextureID convert_to_ImTextureID(GLuint gl_texture_id_in);
        
        /**
         * \brief Recursively browses a directory and generates
         *  menu items.
         * \param[in] path the path to be browsed
         * \param[in] subdirs if true, browse subdirectories as well
         */
        void browse(const std::string& path, bool subdirs=false);
        
        /**
         * \brief Initializes graphic objects.
         * \details This function is called once the OpenGL context is 
         *  ready for rendering. Derived classes may overload it to 
         *  initialize additional graphic objects.
         */
        virtual void init_graphics();

        /**
         * \brief Calls Application's instance init_graphics().
         * \details Used as a callback for glup_viewer / GLFW.
         */
        static void init_graphics_callback();

        /**
         * \brief Draws the scene.
         * \details Derived classes need to overload this function.
         */
        virtual void draw_scene();

        /**
         * \brief Calls Application's instance draw_scene().
         * \details Used as a callback for glup_viewer / GLFW.
         */
        static void draw_scene_callback();

        /**
         * \brief Draws the graphic user interface using ImGUI.
         * \details Derived classes may overload this function, and/or
         *  the functions that draw individual GUI elements.
         */
        virtual void draw_gui();

        /**
         * \brief Draws the menu bar.
         */
        virtual void draw_menu_bar();

        /**
         * \brief Draws the load menu and browser.
         */
        virtual void draw_load_menu();

        /**
         * \brief Draws the save menu.
         */
        virtual void draw_save_menu();

	/**
	 * \brief Draws other file operation menu.
	 * \details Default implementation does nothing.
	 *  It can be overloaded to add other menu
	 *  items in the file menu.
	 */
	virtual void draw_fileops_menu();
	
        /**
         * \brief Draws the about box in the file menu.
         */
        virtual void draw_about();

        /**
         * \brief Draws the windows menu.
         */
        virtual void draw_windows_menu();

        /**
         * \brief Draws the application menus.
         * \details Meant to be overloaded by derived classes.
         */
        virtual void draw_application_menus();
        
        /**
         * \brief Draws the left pane.
         * \details The left pane has the viewer properties and the command.
         */
        virtual void draw_left_pane();

        /**
         * \brief Draws the viewer properties window frame and contents.
         */
        virtual void draw_viewer_properties_window();
        
        /**
         * \brief Draws the contents of viewer properties window.
         */
        virtual void draw_viewer_properties();

        /**
         * \brief Draw the object properties window frame and contents.
         */
        virtual void draw_object_properties_window();
        
        /**
         * \brief Draws the contents of the object properties window.
         */
        virtual void draw_object_properties();
        
        /**
         * \brief Draws the active command if any.
         */
        virtual void draw_command();
        
        /**
         * \brief Draws the right pane.
         * \details The right pane has the object properties. This function
         *  is meant to be overloaded by client code.
         */
        virtual void draw_right_pane();

        /**
         * \brief Draws the console.
         */
        virtual void draw_console();

        /**
         * \brief Draws the status bar.
         */
        virtual void draw_status_bar();
        
        /**
         * \brief Calls Application's instance draw_gui()
         * \details Used as a callback for glup_viewer / GLFW.
         */
        static void draw_gui_callback();

        /**
         * \brief Calls Application's instance load()
         * \details Used as a callback for glup_viewer / GLFW.
         */
        static void dropped_file_callback(char* filename);

        /**
         * \brief Makes the console visible.
         */
        void show_console() {
            console_visible_ = true;
        }

        /**
         * \brief Makes the console invisible.
         */
        void hide_console() {
            console_visible_ = false;
        }

        /**
         * \brief Initializes a new colormap from name and xpm data.
         * \details This function can be called only once the OpenGL
         *  context is ready, for instance in the init_graphics() function.
         * \param[in] name the name of the colormap
         * \param[in] xpm_data the image data of the colormap.
         */
        void init_colormap(const std::string& name, const char** xpm_data);

        /**
         * \brief Initializes all the default colormaps.
         * \details This function can be called only once the OpenGL
         *  context is ready, for instance in the init_graphics() function.
         */
        void init_colormaps();

    protected:
        static Application* instance_;

        int argc_;
        char** argv_;
        
        std::string usage_;
        std::string name_;
        std::string path_;
        
        bool left_pane_visible_;
        bool right_pane_visible_;
        bool console_visible_;

        virtual int MENU_HEIGHT() const;
        virtual int PANE_WIDTH() const;
        virtual int CONSOLE_HEIGHT() const;
        virtual int STATUS_HEIGHT() const;
        
        SmartPointer<Console> console_;
        SmartPointer<StatusBar> status_bar_;

        bool lighting_;
	GLenum effect_;

        GLUPclipMode clip_mode_;
        GLuint geogram_logo_texture_;

        struct ColormapInfo {
            ColormapInfo() : texture(0) {
            }
            GLuint texture;
            std::string name;
        };

        vector<ColormapInfo> colormaps_;

        FileDialog load_dialog_;
        FileDialog save_dialog_;
	std::string current_file_;
	
	bool text_editor_visible_;
	TextEditor text_editor_;

	bool fixed_layout_;
        float scaling_;
        bool retina_mode_;

	vec4f background_color_1_;
	vec4f background_color_2_;

#ifdef GEOGRAM_WITH_LUA	
	lua_State* lua_state_;
#else
	void* lua_state_;
#endif	
	bool lua_error_occured_;
	bool owns_lua_state_;
    };

    /*****************************************************************/

    /**
     * \brief An Application that manipulates a single Mesh.
     */
    class GEOGRAM_GFX_API SimpleMeshApplication : public Application {
    public:
        /**
         * \brief Application constructor.
         * \param[in] argc , argv command line arguments copied from main()
         * \param[in] usage the usage string
         * \see CmdLine::parse()
         */
        SimpleMeshApplication(int argc, char** argv, const std::string& usage);

        /**
         * \copydoc Application::supported_read_file_extensions()
         */
        virtual std::string supported_read_file_extensions();

        /**
         * \copydoc Application::supported_write_file_extensions()
         */
        virtual std::string supported_write_file_extensions();

        /**
         * \copydoc Application::draw_object_properties()
         */
        virtual void draw_object_properties();

        /**
         * \copydoc Application::draw_scene()
         */
        virtual void draw_scene();

        /**
         * \copydoc Application::init_graphics()
         */
        virtual void init_graphics();

        /**
         * \copydoc Application::load()
         */
        virtual bool load(const std::string& filename);

        /**
         * \copydoc Application::save()
         */
        virtual bool save(const std::string& filename);
        
        /**
         * \brief Gets the instance.
         * \return a pointer to the current SimpleMeshApplication.
         */
        static SimpleMeshApplication* instance() {
            SimpleMeshApplication* result =
                dynamic_cast<SimpleMeshApplication*>(Application::instance());
            geo_assert(result != nullptr);
            return result;
        }

    protected:

        /**
         * \brief Gets the bounding box of a mesh animation.
         * \details In animated mode, the mesh animation is stored as 
         *  a mesh with 6d coordinates, that correspond to the geometric 
         *  location at the vertices at time 0 and at time 1.
         * \param[in] M_in the mesh
         * \param[out] xyzmin a pointer to the three minimum coordinates
         * \param[out] xyzmax a pointer to the three maximum coordinates
         * \param[in] animate true if displaying a mesh animation
         */
        void get_bbox(
            const Mesh& M_in, double* xyzmin, double* xyzmax, bool animate
        );
        
        /**
         * \brief increments the animation time in the current instance.
         * \details Callback bound to the 't' key
         */
        static void increment_anim_time_callback();

        /**
         * \brief derements the animation time in the current instance.
         * \details Callback bound to the 'r' key
         */
        static void decrement_anim_time_callback();

        /**
         * \brief increments the cells shrinkage.
         * \details Callback bound to the 'w' key
         */
        static void increment_cells_shrink_callback();

        /**
         * \brief decrements the cells shrinkage.
         * \details Callback bound to the 'x' key
         */
        static void decrement_cells_shrink_callback();

    protected:
        /**
         * \brief Gets the mesh.
         * \return a pointer to the Mesh.
         */
        Mesh* mesh() {
            return &mesh_;
        }

        /**
         * \brief Gets the mesh graphics.
         * \return a pointer to the MeshGfx.
         */
        MeshGfx* mesh_gfx() {
            return &mesh_gfx_;
        }

        /**
         * \brief Makes the vertices visible.
         */
        void show_vertices() {
            show_vertices_ = true;
        }

        /**
         * \brief Makes the vertices invisible.
         */
        void hide_vertices() {
            show_vertices_ = false;
        }

        /**
         * \brief Makes the surface facets visible.
         */
        void show_surface() {
            show_surface_ = true;
        }

        /**
         * \brief Makes the surface facets invisible.
         */
        void hide_surface() {
            show_surface_ = false;            
        }

        /**
         * \brief Makes the volume cells visible.
         */
        void show_volume() {
            show_volume_ = true;
        }

        /**
         * \brief Makes the volume cells invisible.
         */
        void hide_volume() {
            show_volume_ = false;            
        }

        /**
         * \brief Makes the attributes visible.
         */
        void show_attributes() {
            show_attributes_ = true;
        }
        
        /**
         * \brief Makes the attributes invisible.
         */
        void hide_attributes() {
            show_attributes_ = false;
        }
        
        /**
         * \brief Adjusts the current minimum and maximum attribute value
         *  to the currently bound attribute if any.
         */
        void autorange();

        /**
         * \brief Gets the list of attribute names.
         * \return the ';'-separated list of attribute names.
         */
        std::string attribute_names();

        /**
         * \brief Sets the currently displayed attribute.
         * \param[in] attribute the name of the attribute
         *  to be displayed, prefixed by element type.
         */
        void set_attribute(const std::string& attribute);
        
    protected:
        Mesh mesh_;
        MeshGfx mesh_gfx_;
        std::string file_extensions_;

        float anim_speed_;
        float anim_time_;

        bool show_vertices_;
        bool show_vertices_selection_;
        float vertices_size_;
	vec4f vertices_color_;

        bool show_surface_;
        bool show_surface_sides_;        
        bool show_mesh_;
	float mesh_width_;
	vec4f mesh_color_;
	
        bool show_surface_borders_;
	vec4f surface_color_;
	vec4f surface_color_2_;

        bool show_volume_;
        float cells_shrink_;
	vec4f volume_color_;
        bool show_colored_cells_;
        bool show_hexes_;
	bool show_connectors_;

        bool show_attributes_;
        GLuint current_colormap_texture_;
        std::string       attribute_;
        MeshElementsFlags attribute_subelements_;
        std::string       attribute_name_;
        float             attribute_min_;
        float             attribute_max_;
    };

    /*****************************************************************/    
}

#endif
