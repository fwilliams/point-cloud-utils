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

#include <geogram/basic/file_system.h>
#include <geogram/basic/line_stream.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/string.h>

#include <iostream>
#include <fstream>
#include <assert.h>

#ifdef GEO_OS_WINDOWS
#include <windows.h>
#include <io.h>
#include <shlobj.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dirent.h>
#include <unistd.h>
#include <stdio.h>
#endif

namespace GEO {

    namespace FileSystem {

        // OS-dependent functions
#ifdef GEO_OS_WINDOWS

        bool is_file(const std::string& path) {
            WIN32_FIND_DATA file;
            HANDLE file_handle = FindFirstFile(path.c_str(), &file);
            if(file_handle == INVALID_HANDLE_VALUE) {
                return false;
            }
            FindClose(file_handle);
            return (file.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) == 0;
        }

        bool is_directory(const std::string& path) {
            WIN32_FIND_DATA file;
            HANDLE file_handle = FindFirstFile(path.c_str(), &file);
            if(file_handle == INVALID_HANDLE_VALUE) {
                return false;
            }
            FindClose(file_handle);
            return (file.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
        }

        bool create_directory(const std::string& path_in) {
            std::vector<std::string> path;
            String::split_string(path_in, '/', path);
            std::string current;
            int start_at = 0;
            if(path_in.at(1) == ':') {
                current += path_in.at(0);
                current += path_in.at(1);
                start_at = 1;
            }
            else if(path_in.at(0) != '/' && path_in.at(0) != '\\') {
                current += get_current_working_directory();
            }
            for(size_t i = start_at; i < path.size(); i++) {
                current += "/";
                current += path[i];
                if(path[i].at(0) == '.' &&
                    path[i].at(1) == '.' &&
                    path[i].length() == 2
                ) {
                    continue;
                }
                if(!is_directory(current)) {
                    if(!::CreateDirectory(current.c_str(), nullptr)) {
                        Logger::err("OS")
                            << "Could not create directory "
                            << current << std::endl;
                        return false;
                    }
                }
            }
            return true;
        }

        bool delete_directory(const std::string& path) {
            return ::RemoveDirectory(path.c_str()) != FALSE;
        }

        bool delete_file(const std::string& path) {
            return ::DeleteFile(path.c_str()) != FALSE;
        }

        bool get_directory_entries(
            const std::string& path, std::vector<std::string>& result
        ) {
            std::string dirname = path;
            if(dirname.at(dirname.size() - 1) != '/' &&
                dirname.at(dirname.size() - 1) != '\\'
            ) {
                dirname += '/';
            }

            std::string current_directory = get_current_working_directory();

            bool dir_found = set_current_working_directory(dirname);
            if(!dir_found) {
                return false;
            }

            WIN32_FIND_DATA file;
            HANDLE file_handle = FindFirstFile("*.*", &file);
            if(file_handle != INVALID_HANDLE_VALUE) {
                do {
                    std::string file_name = file.cFileName;
                    if(file_name != "." && file_name != "..") {
                        file_name = dirname + file_name;
                        flip_slashes(file_name);
                        result.push_back(file_name);
                    }
                } while(FindNextFile(file_handle, &file));
                FindClose(file_handle);
            }
            set_current_working_directory(current_directory);
            return true;
        }

        std::string get_current_working_directory() {
            char buf[2048];
            std::string result = "";
            if(GetCurrentDirectory(sizeof(buf), buf)) {
                result = buf;
                flip_slashes(result);
            }
            return result;
        }

        bool set_current_working_directory(const std::string& path_in) {
            std::string path = path_in;
            if(
		path.at(path.size() - 1) != '/' &&
		path.at(path.size() - 1) != '\\') {
                path += "/";
            }
            return SetCurrentDirectory(path.c_str()) != -1;
        }

        bool rename_file(
            const std::string& old_name, const std::string& new_name
        ) {
            return ::rename(old_name.c_str(), new_name.c_str()) != -1;
        }

        Numeric::uint64 get_time_stamp(
            const std::string& path
        ) {
            WIN32_FILE_ATTRIBUTE_DATA infos;
            if(!GetFileAttributesEx(
		   path.c_str(), GetFileExInfoStandard, &infos)
	    ) {
                return 0;
            }
            return infos.ftLastWriteTime.dwLowDateTime;
        }

#else

        bool is_file(const std::string& path) {
            //   We use 'stat' and not 'lstat' since
            // we want to be able to follow symbolic
            // links (required for instance when testing
            // input path in GOMGEN, there can be some
            // symlinks in system includes)
            struct stat buff;
            if(stat(path.c_str(), &buff)) {
                return false;
            }
            return S_ISREG(buff.st_mode);
        }

        bool is_directory(const std::string& path) {
            //   We use 'stat' and not 'lstat' since
            // we want to be able to follow symbolic
            // links (required for instance when testing
            // input path in GOMGEN, there can be some
            // symlinks in system includes)
            struct stat buff;
            if(stat(path.c_str(), &buff)) {
                return false;
            }
            return S_ISDIR(buff.st_mode);
        }

        bool create_directory(const std::string& path_in) {
            std::vector<std::string> path;
            String::split_string(path_in, '/', path);
            std::string current;
            for(size_t i = 0; i < path.size(); i++) {
                current += "/";
                current += path[i];
                if(!is_directory(current)) {
                    if(mkdir(current.c_str(), 0755) != 0) {
                        Logger::err("OS")
                            << "Could not create directory "
                            << current << std::endl;
                        return false;
                    }
                }
            }
            return true;
        }

        bool delete_directory(const std::string& path) {
            return rmdir(path.c_str()) == 0;
        }

        bool delete_file(const std::string& path) {
            return unlink(path.c_str()) == 0;
        }

        bool get_directory_entries(
            const std::string& path, std::vector<std::string>& result
        ) {
            std::string dirname = path;
            if(dirname[dirname.length() - 1] != '/') {
                dirname += "/";
            }
            DIR* dir = opendir(dirname.c_str());
            if(dir == nullptr) {
                Logger::err("OS")
                    << "Could not open directory " << dirname
                    << std::endl;
                return false;
            }
            struct dirent* entry = readdir(dir);
            while(entry != nullptr) {
                std::string current = std::string(entry->d_name);
                // Ignore . and ..
                if(current != "." && current != "..") {
                    if(dirname != "./") {
                        current = dirname + current;
                    }
                    // Ignore symbolic links and other special Unix stuff
                    if(is_file(current) || is_directory(current)) {
                        result.push_back(current);
                    }
                }
                entry = readdir(dir);
            }
            closedir(dir);
            return true;
        }

        std::string get_current_working_directory() {
            char buff[4096];
            return std::string(getcwd(buff, 4096));
        }

        bool set_current_working_directory(const std::string& path) {
            return chdir(path.c_str()) == 0;
        }

        bool rename_file(
            const std::string& old_name, const std::string& new_name
        ) {
            if(is_file(new_name)) {
                return false;
            }
            return ::rename(old_name.c_str(), new_name.c_str()) == 0;
        }

        Numeric::uint64 get_time_stamp(
            const std::string& path
        ) {
            struct stat buffer;
            if(!stat(path.c_str(), &buffer)) {
                return Numeric::uint64(buffer.st_mtime);
            }
            return 0;
        }

#endif

        // OS-independent functions

        std::string extension(const std::string& path) {
            size_t len = path.length();
            if(len != 0) {
                for(size_t i = len - 1; i != 0; i--) {
                    if(path[i] == '/' || path[i] == '\\') {
                        break;
                    }
                    if(path[i] == '.') {
                        return String::to_lowercase(path.substr(i + 1));
                    }
                }
            }
            return std::string();
        }

        std::string base_name(const std::string& path, bool remove_extension) {
            long int len = (long int)(path.length());
            if(len == 0) {
                return std::string();
            }
            long int dot_pos = len;
            long int i;
            for(i = len - 1; i >= 0; i--) {
                if(path[size_t(i)] == '/' || path[size_t(i)] == '\\') {
                    break;
                }
                if(remove_extension && path[size_t(i)] == '.') {
                    dot_pos = i;
                }
            }
            return path.substr(size_t(i + 1), size_t(dot_pos - i - 1));
        }

        std::string dir_name(const std::string& path) {
            size_t len = path.length();
            if(len != 0) {
                for(size_t i = len - 1; i != 0; i--) {
                    if(path[i] == '/' || path[i] == '\\') {
                        return path.substr(0, i);
                    }
                }
            }
            return ".";
        }

        void get_directory_entries(
            const std::string& path,
            std::vector<std::string>& result, bool recursive
        ) {
	    // TODO: seems to be bugged, enters infinite recursion...
            get_directory_entries(path, result);
            if(recursive) {
                for(size_t i = 0; i < result.size(); i++) {
                    if(is_directory(result[i])) {
                        get_directory_entries(result[i], result, true);
                    }
                }
            }
        }

        void get_files(
            const std::string& path,
            std::vector<std::string>& result, bool recursive
        ) {
            std::vector<std::string> entries;
            get_directory_entries(path, entries, recursive);
            for(size_t i = 0; i < entries.size(); i++) {
                if(is_file(entries[i])) {
                    result.push_back(entries[i]);
                }
            }
        }

        void get_subdirectories(
            const std::string& path,
            std::vector<std::string>& result, bool recursive
        ) {
            std::vector<std::string> entries;
            get_directory_entries(path, entries, recursive);
            for(size_t i = 0; i < entries.size(); i++) {
                if(is_directory(entries[i])) {
                    result.push_back(entries[i]);
                }
            }
        }

        void flip_slashes(std::string& s) {
            for(size_t i = 0; i < s.length(); i++) {
                if(s[i] == '\\') {
                    s[i] = '/';
                }
            }
        }

        bool copy_file(const std::string& from, const std::string& to) {
            FILE* fromf = fopen(from.c_str(), "rb");
            if(fromf == nullptr) {
                Logger::err("FileSyst")
		    << "Could not open source file:" << from << std::endl;
                return false;
            }
            FILE* tof = fopen(to.c_str(),"wb");
            if(tof == nullptr) {
                Logger::err("FileSyst")
		    << "Could not create file:" << to << std::endl;
                fclose(fromf);
                return false;
            }

            bool result = true;
            const size_t buff_size = 4096;
            char buff[buff_size];
            size_t rdsize = 0;
            do {
                rdsize = fread(buff, 1, buff_size, fromf);
                if(fwrite(buff, 1, rdsize, tof) != rdsize) {
                    Logger::err("FileSyst") << "I/O error when writing to file:"
                                            << to << std::endl;
                    result = false;
                    break;
                }
            } while(rdsize == 4096);
            
            fclose(fromf);
            fclose(tof);
            return result;
        }

        void set_executable_flag(const std::string& filename) {
            geo_argused(filename);
#ifdef GEO_OS_UNIX
            if(::chmod(filename.c_str(), 0755) != 0) {
                Logger::err("FileSyst")
                    << "Could not change file permissions for:"
                    << filename << std::endl;
            }
#endif            
        }

        void GEOGRAM_API touch(const std::string& filename) {
#ifdef GEO_OS_APPLE
           {
                struct stat buff;
                int rc = stat(filename.c_str(), &buff); 
                if(rc != 0) {  // FABIEN NOT SURE WE GET THE TOUCH
                    Logger::err("FileSystem")
                        << "Could not touch file:"
                        << filename
                        << std::endl;
                }
            }
#elif defined(GEO_OS_UNIX)
            {
                int rc = utimensat(
                    AT_FDCWD,
                    filename.c_str(),
                    nullptr,
                    0
                );
                if(rc != 0) {
                    Logger::err("FileSystem")
                        << "Could not touch file:"
                        << filename
                        << std::endl;
                }
            }
#elif defined GEO_OS_WINDOWS
            {
                HANDLE hfile = CreateFile(
                    filename.c_str(),
                    GENERIC_READ | GENERIC_WRITE,
                    FILE_SHARE_READ | FILE_SHARE_WRITE,
                    nullptr,
                    OPEN_EXISTING,
                    FILE_ATTRIBUTE_NORMAL,
                    nullptr
                );
                if(hfile == INVALID_HANDLE_VALUE) {
                    Logger::err("FileSystem")
                        << "Could not touch file:"
                        << filename
                        << std::endl;
                }
                SYSTEMTIME now_system;
                FILETIME now_file;
                GetSystemTime(&now_system);
                SystemTimeToFileTime(&now_system, &now_file);
                SetFileTime(hfile,nullptr,&now_file,&now_file);
                CloseHandle(hfile);
            }
#endif            
        }
        
        std::string normalized_path(const std::string& path_in) {

            if(path_in == "") {
                return "";
            }
            
            std::string path = path_in;
            std::string result;

#ifdef GEO_OS_UNIX
            // If this is a relative path, prepend "./"
            if(path[0] != '/') {
                path = "./" + path;
            }

            char buffer[PATH_MAX];
            char* p = realpath(path.c_str(), buffer);
            if(p != nullptr) {
                result = std::string(p);
            } else {
                // realpath() only works for existing paths and existing file,
                // therefore we attempt calling it on the input path by adding
                // one component at a time.
                size_t pos = 1;
                while(pos != std::string::npos) {
                    pos = path.find('/',pos);
                    if(pos != std::string::npos) {
                        std::string path_part = path.substr(0,pos);
                        p = realpath(path_part.c_str(), buffer);
                        if(p == nullptr) {
                            break;
                        } else {
                            result = std::string(p) +
                                path.substr(pos, path.length()-pos);
                        }
                        ++pos;
                        if(pos == path.length()) {
                            break;
                        }
                    } 
                }
            }
#endif
            
#ifdef GEO_OS_WINDOWS
            TCHAR buffer[MAX_PATH];
            GetFullPathName(path.c_str(), MAX_PATH, buffer, nullptr);
            result = std::string(buffer);
#endif
            
            flip_slashes(result);
            return result;
        }


        std::string home_directory() {
            std::string home;
#if defined GEO_OS_WINDOWS
            wchar_t folder[MAX_PATH+1];
            HRESULT hr = SHGetFolderPathW(0, CSIDL_PROFILE, 0, 0, folder);
            if (SUCCEEDED(hr)) {
                char result[MAX_PATH+1];
                wcstombs(result, folder, MAX_PATH);
                home=std::string(result);
                flip_slashes(home);
            }
#elif defined GEO_OS_EMSCRIPTEN
            home="/";
#else            
            char* result = getenv("HOME");
            if(result != nullptr) {
                home=result;
            }
#endif
            return home;
        }

        std::string documents_directory() {
            std::string home;
#if defined GEO_OS_WINDOWS
            wchar_t folder[MAX_PATH+1];
            HRESULT hr = SHGetFolderPathW(0, CSIDL_MYDOCUMENTS, 0, 0, folder);
            if (SUCCEEDED(hr)) {
                char result[MAX_PATH+1];
                wcstombs(result, folder, MAX_PATH);
                home=std::string(result);
                flip_slashes(home);
            }
#elif defined GEO_OS_EMSCRIPTEN
            home="/";
#elif defined GEO_OS_ANDROID
            char* result = getenv("EXTERNAL_STORAGE");
            if(result != nullptr) {
                home=result;
            }
#else            
            char* result = getenv("HOME");
            if(result != nullptr) {
                home=result;
            }
#endif
            return home;
        }
        
    }

    
}

