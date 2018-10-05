/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef __MEMORY_INFO_H
#define __MEMORY_INFO_H

#include <stdexcept>

namespace vcg 
{
	//WARNING: All the classes derived from MemoryInfo has been intended to be instantiated as a singleton in the host application 
	//(i.e. in every application using it just an instance of a class derived from MemoryInfo should be declared).

	class MemoryInfo
	{
	public:
		class MemoryInfoException : public std::exception
		{
		public:
			MemoryInfoException(const char* text)
				:std::exception(),_exctext(text){}

			~MemoryInfoException() throw() {}
			inline const char* what() const throw() {return _exctext;}
		private:
			const char* _exctext;
		};

		MemoryInfo(std::ptrdiff_t originalmem)
			:_originaltotalmemory(originalmem),_currentfreememory(_originaltotalmemory)
		{       
		}

		virtual ~MemoryInfo() {}
		virtual void acquiredMemory(std::ptrdiff_t mem) = 0;
		virtual std::ptrdiff_t usedMemory() const = 0;
		virtual std::ptrdiff_t currentFreeMemory() const = 0;
		virtual void releasedMemory(std::ptrdiff_t mem = 0) = 0;
		virtual bool isAdditionalMemoryAvailable(std::ptrdiff_t mem) = 0;

	protected:
		const std::ptrdiff_t _originaltotalmemory;
		std::ptrdiff_t _currentfreememory;
	};

	//WARNING: this is not a thread safe class. The object derived from MemoryInfo are intended to be used inside GLMeshAttributeFeeder as static variable in order to manage the available GPUMemory.
	//We strongly recommend you to define in your code a thread safe version of the class, defining mutexed access member functions. 
	//This class should be consider just as a basic example for the implementations of the required functionalities. 
	//It is safe to use it just when the user has only one mesh to pass to the GPU.

	class NotThreadSafeMemoryInfo : public MemoryInfo
	{
	public:
		NotThreadSafeMemoryInfo(std::ptrdiff_t originalmem)
			:MemoryInfo(originalmem)
		{
		}

		~NotThreadSafeMemoryInfo() {}

		void acquiredMemory(std::ptrdiff_t mem)
		{
			/*if (mem > _originaltotalmemory)
				throw MemoryInfo::MemoryInfoException("It has been requested more memory than the total one.\\n");
			else 
				if (mem > _currentfreememory)
					throw MemoryInfo::MemoryInfoException("It has been requested more memory than the free available one.\\n");
				else*/
					_currentfreememory -= mem;
		}

		std::ptrdiff_t usedMemory() const
		{
			return _originaltotalmemory - _currentfreememory;
		}

		std::ptrdiff_t currentFreeMemory() const
		{
			return _currentfreememory;
		}

		void releasedMemory(std::ptrdiff_t mem = 0)
		{
			/*if (mem > _originaltotalmemory)
				throw MemoryInfo::MemoryInfoException("It has been released more memory than the total one. Something strange happened!\\n");
			else*/
				_currentfreememory += mem;
		}

		bool isAdditionalMemoryAvailable(std::ptrdiff_t mem)
		{
			return (_currentfreememory >= mem);
		}
	};
};

#endif