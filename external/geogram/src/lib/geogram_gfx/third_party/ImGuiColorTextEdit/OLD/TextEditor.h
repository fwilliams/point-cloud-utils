#pragma once

// [Bruno Levy] For the older compilers that does not have
// regex support, uncomment to deactivate syntax highlighting
// #define GEO_OLD_COMPILER

#include <string>
#include <vector>
#include <memory>
#include <map>

#ifndef GEO_OLD_COMPILER
#include <regex>
#include <unordered_set>
#include <unordered_map>
#include <array>
#else
#include <set>
namespace std {

    /**
     * \brief Minimal drop-in replacement for class std::array.
     */
    template <class T, unsigned int SIZE> class array {
      public:
	typedef array<T,SIZE> thisclass;

	array() {
	}

	array(const thisclass& rhs) {
	    for(unsigned int i=0; i<SIZE; ++i) {
		data_[i] = rhs.data_[i];
	    }
	}

	thisclass& operator=(const thisclass& rhs) {
	    for(unsigned int i=0; i<SIZE; ++i) {
		data_[i] = rhs.data_[i];
	    }
	    return *this;
	}
	
	T& operator[](unsigned int idx) {
	    return data_[idx];
	}
	const T& operator[](unsigned int idx) const {
	    return data_[idx];		
	}

      private:
	T data_[SIZE];
    };

    /**
     * \brief Drop-in replacement for std::move(), that does
     *  nothing here.
     */
    template <class T> inline const T& move(const T& x) {
	return x;
    }
    
}
#endif

// [Bruno Levy] include path redirected to geogram.
#include <geogram_gfx/api/defs.h>
#include <geogram_gfx/third_party/ImGui/imgui.h>

// [Bruno Levy] additional callback type.
enum TextEditorAction {
    TEXT_EDITOR_RUN, TEXT_EDITOR_SAVE,
    TEXT_EDITOR_FIND, TEXT_EDITOR_STOP,
    TEXT_EDITOR_TOOLTIP,
    TEXT_EDITOR_COMPLETION,
    TEXT_EDITOR_TEXT_CHANGED
};

typedef void (*TextEditorCallback)(
    TextEditorAction action, void* client_data
);

class GEOGRAM_GFX_API TextEditor
{
public:
	enum class PaletteIndex : uint8_t
	{
		Default,
		Keyword,
		Number,
		String,
		CharLiteral,
		Punctuation,
		Preprocessor,
		Identifier,
		KnownIdentifier,
		PreprocIdentifier,
		Comment,
		MultiLineComment,
		Background,
		Cursor,
		Selection,
		ErrorMarker,
		Breakpoint,
		LineNumber,
		CurrentLineFill,
		CurrentLineFillInactive,
		CurrentLineEdge,
		Max
	};

	struct Breakpoint
	{
		int mLine;
		bool mEnabled;
		std::string mCondition;

		Breakpoint()
			: mLine(-1)
			, mEnabled(false)
		{}
	};
	
	struct Coordinates
	{
		int mLine, mColumn;
		Coordinates() : mLine(0), mColumn(0) {}
		Coordinates(int aLine, int aColumn) : mLine(aLine), mColumn(aColumn)
		{
			assert(aLine >= 0);
			assert(aColumn >= 0);
		}
		static Coordinates Invalid() { static Coordinates invalid(-1, -1); return invalid; }

		bool operator ==(const Coordinates& o) const
		{
			return
				mLine == o.mLine &&
				mColumn == o.mColumn;
		}

		bool operator !=(const Coordinates& o) const
		{
			return
				mLine != o.mLine ||
				mColumn != o.mColumn;
		}

		bool operator <(const Coordinates& o) const
		{
			if (mLine != o.mLine)
				return mLine < o.mLine;
			return mColumn < o.mColumn;
		}

		bool operator >(const Coordinates& o) const
		{
			if (mLine != o.mLine)
				return mLine > o.mLine;
			return mColumn > o.mColumn;
		}

		bool operator <=(const Coordinates& o) const
		{
			if (mLine != o.mLine)
				return mLine < o.mLine;
			return mColumn <= o.mColumn;
		}

		bool operator >=(const Coordinates& o) const
		{
			if (mLine != o.mLine)
				return mLine > o.mLine;
			return mColumn >= o.mColumn;
		}
	};

	struct Identifier
	{
		Coordinates mLocation;
		std::string mDeclaration;
	};

	typedef std::string String;

	// [Bruno Levy] for old compilers that do not have std::unordered_set
	// and map.
#ifdef GEO_OLD_COMPILER
	typedef std::map<std::string, Identifier> Identifiers;
	typedef std::set<std::string> Keywords;
	typedef std::set<int> Breakpoints;
#else	
	typedef std::unordered_map<std::string, Identifier> Identifiers;
	typedef std::unordered_set<std::string> Keywords;
	typedef std::unordered_set<int> Breakpoints;
#endif
	
	typedef std::map<int, std::string> ErrorMarkers;
	typedef std::array<ImU32, (unsigned)PaletteIndex::Max> Palette;
	typedef char Char;
	
	struct Glyph
	{
		Char mChar;
	        PaletteIndex mColorIndex : 8; // [Bruno Levy]: 7 is too small to hold all values
	        bool mMultiLineComment : 1;

		Glyph(Char aChar, PaletteIndex aColorIndex) : mChar(aChar), mColorIndex(aColorIndex), mMultiLineComment(false) {}
	};

	typedef std::vector<Glyph> Line;
	typedef std::vector<Line> Lines;

	struct GEOGRAM_GFX_API LanguageDefinition
	{
		typedef std::pair<std::string, PaletteIndex> TokenRegexString;
		typedef std::vector<TokenRegexString> TokenRegexStrings;

		std::string mName;
		Keywords mKeywords;
		Identifiers mIdentifiers;
		Identifiers mPreprocIdentifiers;
		std::string mCommentStart, mCommentEnd;

		TokenRegexStrings mTokenRegexStrings;

		bool mCaseSensitive;

		static LanguageDefinition CPlusPlus();
		static LanguageDefinition HLSL();
		static LanguageDefinition GLSL();
		static LanguageDefinition C();
		static LanguageDefinition SQL();
		static LanguageDefinition AngelScript();
		static LanguageDefinition Lua();
	};

	TextEditor();
	~TextEditor();

	void SetLanguageDefinition(const LanguageDefinition& aLanguageDef);
	const LanguageDefinition& GetLanguageDefinition() const { return mLanguageDefinition; }

	const Palette& GetPalette() const { return mPalette; }
	void SetPalette(const Palette& aValue);

	void SetErrorMarkers(const ErrorMarkers& aMarkers) { mErrorMarkers = aMarkers; }
	const ErrorMarkers& GetErrorMarkers() const { return mErrorMarkers; } //[Bruno]
	
	void SetBreakpoints(const Breakpoints& aMarkers) { mBreakpoints = aMarkers; }
	const Breakpoints& GetBreakpoints() const { return mBreakpoints; } //[Bruno]
	
	void Render(const char* aTitle, const ImVec2& aSize = ImVec2(), bool aBorder = false);
	void SetText(const std::string& aText);
	std::string GetText() const;
	//[Bruno] made public
	std::string GetText(const Coordinates& aStart, const Coordinates& aEnd) const;	
	std::string GetSelectedText() const;
	Coordinates GetSelectionStart() const {
	    return mState.mSelectionStart;
	}
	Coordinates GetSelectionEnd() const {
	    return mState.mSelectionEnd;
	}
	int nbLines() const {
	    return int(mLines.size());
	}
	std::string GetLine(int line) const; 
	
	int GetTotalLines() const { return (int)mLines.size(); }
	bool IsOverwrite() const { return mOverwrite; }

	void SetReadOnly(bool aValue);
	bool IsReadOnly() const { return mReadOnly; }

	Coordinates GetCursorPosition() const { return GetActualCursorCoordinates(); }

	// [Bruno Levy] added this function.
	Coordinates GetMousePosition() const { return ScreenPosToCoordinates(ImGui::GetMousePos()); }
	
	void SetCursorPosition(const Coordinates& aPosition);

	void InsertText(const std::string& aValue);
	void InsertText(const char* aValue);

	void MoveUp(int aAmount = 1, bool aSelect = false);
	void MoveDown(int aAmount = 1, bool aSelect = false);
	void MoveLeft(int aAmount = 1, bool aSelect = false, bool aWordMode = false);
	void MoveRight(int aAmount = 1, bool aSelect = false, bool aWordMode = false);
	void MoveTop(bool aSelect = false);
	void MoveBottom(bool aSelect = false);
	void MoveHome(bool aSelect = false);
	void MoveEnd(bool aSelect = false);

	void SetSelectionStart(const Coordinates& aPosition);
	void SetSelectionEnd(const Coordinates& aPosition);
	void SetSelection(const Coordinates& aStart, const Coordinates& aEnd, bool awordmode = false);
	void SelectWordUnderCursor();
	bool HasSelection() const;

	void Copy();
	void Cut();
	void Paste();
	void Delete();

	bool CanUndo() const;
	bool CanRedo() const;
	void Undo(int aSteps = 1);
	void Redo(int aSteps = 1);

	static const Palette& GetDarkPalette();
	static const Palette& GetLightPalette();

	// [Bruno Levy] gets the text source element under the pointer
	// when a tooltip should be displayed.
	const std::string& GetWordContext() const {
	    return word_context_;
	}
	
private:

// [Bruno Levy] Made syntax hightlighting optional	
#ifndef GEO_OLD_COMPILER   
	typedef std::vector<std::pair<std::regex, PaletteIndex>> RegexList;
#endif
	
	struct EditorState
	{
		Coordinates mSelectionStart;
		Coordinates mSelectionEnd;
		Coordinates mCursorPosition;
	};

	class UndoRecord
	{
	public:
		UndoRecord() {}
		~UndoRecord() {}

		UndoRecord(
			const std::string& aAdded,
			const TextEditor::Coordinates aAddedStart, 
			const TextEditor::Coordinates aAddedEnd, 
			
			const std::string& aRemoved, 
			const TextEditor::Coordinates aRemovedStart,
			const TextEditor::Coordinates aRemovedEnd,
			
			TextEditor::EditorState& aBefore, 
			TextEditor::EditorState& aAfter);

		void Undo(TextEditor* aEditor);
		void Redo(TextEditor* aEditor);

		std::string mAdded;
		Coordinates mAddedStart;
		Coordinates mAddedEnd;

		std::string mRemoved;
		Coordinates mRemovedStart;
		Coordinates mRemovedEnd;

		EditorState mBefore;
		EditorState mAfter;
	};

	typedef std::vector<UndoRecord> UndoBuffer;

	void ProcessInputs();
	void Colorize(int aFromLine = 0, int aCount = -1);
	void ColorizeRange(int aFromLine = 0, int aToLine = 0);
	void ColorizeInternal();
	int TextDistanceToLineStart(const Coordinates& aFrom) const;
	void EnsureCursorVisible();
	int GetPageSize() const;
	int AppendBuffer(std::string& aBuffer, char chr, int aIndex);

	Coordinates GetActualCursorCoordinates() const;
	Coordinates SanitizeCoordinates(const Coordinates& aValue) const;
	void Advance(Coordinates& aCoordinates) const;
	void DeleteRange(const Coordinates& aStart, const Coordinates& aEnd);
	int InsertTextAt(Coordinates& aWhere, const char* aValue);
	void AddUndo(UndoRecord& aValue);
	Coordinates ScreenPosToCoordinates(const ImVec2& aPosition) const;
	Coordinates FindWordStart(const Coordinates& aFrom) const;
	Coordinates FindWordEnd(const Coordinates& aFrom) const;

	// [Bruno Levy] additional function for finding the context for
	//  tooltips
	bool IsWordContextBoundary(char c) const;
	Coordinates FindWordContextStart(const Coordinates& aFrom) const;
	Coordinates FindWordContextEnd(const Coordinates& aFrom) const;
	std::string GetWordContextAt(const Coordinates & aCoords) const;
	
	bool IsOnWordBoundary(const Coordinates& aAt) const;
	void RemoveLine(int aStart, int aEnd);
	void RemoveLine(int aIndex);
	Line& InsertLine(int aIndex);
	void EnterCharacter(Char aChar);
	void BackSpace();
	void DeleteSelection();
	std::string GetWordUnderCursor() const;
	std::string GetWordAt(const Coordinates& aCoords) const;

	float mLineSpacing;
	Lines mLines;
	EditorState mState;
	UndoBuffer mUndoBuffer;
	int mUndoIndex;
	
	int mTabSize;
	bool mOverwrite;
	bool mReadOnly;
	bool mWithinRender;
	bool mScrollToCursor;
	bool mWordSelectionMode;
	int mColorRangeMin, mColorRangeMax;

	Palette mPalette;
	LanguageDefinition mLanguageDefinition;

// [Bruno Levy] Made syntax hightlighting optional		
#ifndef GEO_OLD_COMPILER
	RegexList mRegexList;
#endif
	bool mCheckMultilineComments;
	Breakpoints mBreakpoints;
	ErrorMarkers mErrorMarkers;
	ImVec2 mCharAdvance;
	Coordinates mInteractiveStart, mInteractiveEnd;

// [Bruno Levy] Additional callback.
  public:
	void set_callback(
	    TextEditorCallback cb, void* cb_cli_data = nullptr
	) {
	    callback_ = cb;
	    callback_client_data_ = cb_cli_data;
	}
  private:
	TextEditorCallback callback_;
	void* callback_client_data_;
	std::string word_context_;
};

