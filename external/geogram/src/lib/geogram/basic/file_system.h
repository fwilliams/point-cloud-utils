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

#ifndef GEOGRAM_BASIC_FILE_SYSTEM
#define GEOGRAM_BASIC_FILE_SYSTEM

#include <geogram/basic/common.h>
#include <geogram/basic/numeric.h>
#include <string>
#include <vector>

/**
 * \file geogram/basic/file_system.h
 * \brief Functions and times for filesystem manipulation
 */

namespace GEO {

    /**
     * \brief Abstraction layer for file-system management.
     */
    namespace FileSystem {

        /**
         * \brief Checks if a path is a regular file.
         * \param[in] path system path to verify.
         * \retval true if \p path is a regular file.
         * \retval false otherwise.
         */
        bool GEOGRAM_API is_file(const std::string& path);

        /**
         * \brief Checks if a path is a directory.
         * \param[in] path system path to verify.
         * \retval true if \p path is a directory.
         * \retval false otherwise.
         */
        bool GEOGRAM_API is_directory(const std::string& path);

        /**
         * \brief Creates a directory
         * \details This recursively creates a new directory given by its \b
         * absolute path \p path, creating any missing intermediate
         * directories on the fly.
         * \param[in] path absolute path to the directory to be created.
         * \retval true if the directory was successfully created.
         * \retval false otherwise.
         */
        bool GEOGRAM_API create_directory(const std::string& path);

        /**
         * \brief Deletes a directory
         * \details This deletes the directory specified by path \p path. The
         * path must specify an empty directory.
         * \param[in] path the path of the directory to be removed.
         * \retval true if the directory was successfully deleted.
         * \retval false otherwise.
         */
        bool GEOGRAM_API delete_directory(const std::string& path);

        /**
         * \brief Deletes a file
         * \param[in] path the path of the file to be deleted.
         * \retval true if the file path was successfully deleted
         * \retval false otherwise
         */
        bool GEOGRAM_API delete_file(const std::string& path);

        /**
         * \brief Lists directory contents
         * \details Lists all the files and sub-directories in the directory
         * specified by \p path, and stores the list in \p result. Special
         * entries "." and ".." are not stored in \p result.
         * \param[in] path the path to the directory to list.
         * \param[in] result output vector of files and sub-directories.
         * \retval true if \p path specifies a readable directory.
         * \retval false otherwise.
         */
        bool GEOGRAM_API get_directory_entries(
            const std::string& path, std::vector<std::string>& result
        );

        /**
         * \brief Gets the current working directory.
         * \return The absolute path to the current directory.
         */
        std::string GEOGRAM_API get_current_working_directory();

        /**
         * \brief Sets the working directory.
         * \param[in] path path to the new working directory.
         * \retval true if the current directory could be changed to \p path.
         * \retval false otherwise.
         */
        bool GEOGRAM_API set_current_working_directory(
            const std::string& path
        );

        /**
         * \brief Renames or moves a file.
         * \details This renames the existing file or directory specified by
         * path \p old_name to the new path \p new_name. The new name must not
         * be the name of an existing file or directory. If \p old_name and \p
         * new_name are not in the same directory, \p old_name is moved to the
         * \p new_name.
         * \param[in] old_name path of the file or directory to be renamed.
         * \param[in] new_name new path of the file or directory.
         * \retval true if the file was renamed successfully.
         * \retval false otherwise.
         */
        bool GEOGRAM_API rename_file(
            const std::string& old_name, const std::string& new_name
        );

        /**
         * \brief Gets a file last modification time.
         * \param[in] path the path to an existing file or directory.
         * \return the last modification time in seconds
         */
        Numeric::uint64 GEOGRAM_API get_time_stamp(
            const std::string& path
        );

        /**
         * \brief Gets a path extension
         * \details Extracts the extension from the path \p path, that is any
         * character that appear after the last dot (.) and after any
         * directory separator character. If \p path has no extension, the
         * empty string is returned.
         *
         * Examples
         * - extension("/dir/file.cpp") -> "cpp"
         * - extension("file") -> ""
         * - extension("/dir.ext/file") -> ""
         *
         * \param[in] path the path to a file or directory
         * \return the path's extension (without the dot), or the empty
         * string if none.
         */
        std::string GEOGRAM_API extension(const std::string& path);

        /**
         * \brief Gets a path base name
         * \details Extracts the base name from the path \p path, that is any
         * character that appear after the last directory separator. If
         * parameter \p remove_extension is \c true (the default), the
         * extension is removed from the base name, otherwise is it kept. If
         * the path does not contain any directory separator, the whole path
         * is returned.
         *
         * Examples
         * - base_name("/dir/file.cpp") -> "file"
         * - base_name("/dir/file.cpp", false) -> "file.cpp"
         * - base_name("file") -> "file"
         *
         * \param[in] path the path to a file or directory
         * \param[in] remove_extension whether to remove the extension from
         * the base name or not.
         */
        std::string GEOGRAM_API base_name(
            const std::string& path, bool remove_extension = true
        );

        /**
         * \brief Gets a path directory
         * \details Extracts the directory from the path \p path, that is any
         * character that appear before the last directory separator. 
         *  If the path does not contain any directory separator, 
         * string "." is returned.
         *
         * Examples
         * - dir_name("/dir/file.cpp") -> "dir"
         * - dir_name("file") -> "."
         * - dir_name("/") -> "/"
         *
         * \param[in] path the path to a file or directory
         * \return the path directory or "." if none
         */
        std::string GEOGRAM_API dir_name(const std::string& path);

        /**
         * \brief Lists directory contents
         * \details Lists all the files and sub-directories in the directory
         * specified by \p path, and stores the list in \p result. Special
         * entries "." and ".." are not stored in \p result. If parameter
         * recursive is set to \c true, \p result will include the entries of
         * all sub-directories in \p path recursively.
         * \param[in] path the path to an existing directory
         * \param[in] result output vector of entries in \p path
         * \param[in] recursive recursively traverses all sub-directories in
         * \p path
         */
        void GEOGRAM_API get_directory_entries(
            const std::string& path,
            std::vector<std::string>& result, bool recursive
        );

        /**
         * \brief Lists files in a directory
         * \details Lists all the files in the directory specified by \p path,
         * and stores the list in \p result. Special entries "." and ".." are
         * not stored in \p result. If parameter recursive is set to \c true,
         * \p result will include the entries of all sub-directories in \p
         * path recursively.
         * \param[in] path the path to an existing directory
         * \param[in] result output vector of files in \p path
         * \param[in] recursive recursively traverses all sub-directories in
         * \p path
         * \see get_directory_entries()
         */
        void GEOGRAM_API get_files(
            const std::string& path,
            std::vector<std::string>& result, bool recursive = false
        );

        /**
         * \brief Lists sub-directories in a directory
         * \details Lists all the sub-directories in the directory specified
         * by \p path, and stores the list in \p result. Special entries "."
         * and ".." are not stored in \p result. If parameter recursive is set
         * to \c true, \p result will include the entries of all
         * sub-directories in \p path recursively.
         * \param[in] path the path to an existing directory
         * \param[in] result output vector of sub-directories in \p path
         * \param[in] recursive recursively traverses all sub-directories in
         * \p path
         * \see get_directory_entries()
         */
        void GEOGRAM_API get_subdirectories(
            const std::string& path,
            std::vector<std::string>& result, bool recursive = false
        );

        /**
         * \brief Converts a path to Unix format
         * \details It changes all Windows "\" directory separators into Unix
         * "/" directory separators.
         * \param[in,out] path the path to be converted
         */
        void GEOGRAM_API flip_slashes(std::string& path);

        /**
         * \brief Copies a file
         * \param[in] from name of the file to be copied
         * \param[out] to name of the copy
         * \retval true if the copy was successful
         * \retval false otherwise
         */
        bool GEOGRAM_API copy_file(
            const std::string& from, const std::string& to
        );

        /**
         * \brief Marks a filename as executable.
         * \details On unix, it chmods the file, on Windows, does nothing.
         * \param[in] filename name of the file to be made executable
         */
        void GEOGRAM_API set_executable_flag(const std::string& filename);


        /**
         * \brief Modifies the last modification time of a file.
         * \param[in] filename name of the file.
         */
        void GEOGRAM_API touch(const std::string& filename);

        /**
         * \brief Normalizes a path.
         * \details A path is normalized if it is absolute and it does not 
         *  contain any "../" component.
         * \param[in] path the path to be normalized. The path can have 
         *  components that do not exist.
         * \return the normalized path
         */
        std::string GEOGRAM_API normalized_path(const std::string& path);


        /**
         * \brief Gets the current user's home directory.
         * \return The path to the current user's home directory as a string.
         */
        std::string GEOGRAM_API home_directory();

        /**
         * \brief Gets the current user's home directory.
         * \details Under unix, it returns the content of the HOME environment
         *  variable. Under Windows, it returns the "My Documents" directory.
         * \return The path to the current user's home directory as a string.
         */
        std::string GEOGRAM_API documents_directory();
    }
}

#endif

