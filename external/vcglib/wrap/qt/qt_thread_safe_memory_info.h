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

#ifndef __QT_THREAD_SAFE_MEMORY_INFO_H
#define __QT_THREAD_SAFE_MEMORY_INFO_H

#include <QReadWriteLock>

#include <wrap/system/memory_info.h>

namespace vcg
{
    class QtThreadSafeMemoryInfo : public vcg::NotThreadSafeMemoryInfo
    {
    public:
	    QtThreadSafeMemoryInfo(std::ptrdiff_t originalmem)
            :vcg::NotThreadSafeMemoryInfo(originalmem),lock(QReadWriteLock::Recursive)
        {
        }

	    ~QtThreadSafeMemoryInfo()
        {
        }

	    void acquiredMemory(std::ptrdiff_t mem)
        {
            QWriteLocker locker(&lock);
            vcg::NotThreadSafeMemoryInfo::acquiredMemory(mem);
        }

	    std::ptrdiff_t usedMemory() const
        {
            QReadLocker locker(&lock);
            return vcg::NotThreadSafeMemoryInfo::usedMemory();
        }

	    std::ptrdiff_t currentFreeMemory() const
        {
            QReadLocker locker(&lock);
            return vcg::NotThreadSafeMemoryInfo::currentFreeMemory();
        }

	    void releasedMemory(std::ptrdiff_t mem = 0)
        {
            QWriteLocker locker(&lock);
            vcg::NotThreadSafeMemoryInfo::releasedMemory(mem);
        }
	
	    bool isAdditionalMemoryAvailable(std::ptrdiff_t mem)
        {
            QReadLocker locker(&lock);
            return vcg::NotThreadSafeMemoryInfo::isAdditionalMemoryAvailable(mem);
        }
    private:
	    //mutable objects can be modified from the declared const functions
	    //in this way we have not to modified the basic vcg::MemoryInfo interface for the logically const functions
	    //whose need to lock the mutex for a simple reading operation
	    mutable QReadWriteLock lock;
    };
}

#endif
