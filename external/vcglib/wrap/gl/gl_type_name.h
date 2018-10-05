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

#ifndef __VCG_GL_TYPE_NAME
#define __VCG_GL_TYPE_NAME

namespace vcg 
{
	template <typename T>
	class GL_TYPE_NM
	{public:
	static GLenum SCALAR() { assert(0); return 0;}
	};
	template <> class GL_TYPE_NM<float>
	{ public:
	typedef GLfloat ScalarType;
	static GLenum SCALAR() { return GL_FLOAT; }
	};
	template <> class GL_TYPE_NM<double>
	{public:
	typedef GLdouble ScalarType;
	static GLenum SCALAR() { return GL_DOUBLE; }
	};
}

#endif