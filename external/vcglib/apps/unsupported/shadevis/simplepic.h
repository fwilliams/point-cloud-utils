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
/****************************************************************************
  History

$Log: not supported by cvs2svn $

****************************************************************************/

#ifndef __VCG_SIMPLE_PIC
#define __VCG_SIMPLE_PIC
#include <vcg/math/matrix44.h>

namespace vcg {
  template <class PixType> 
  class SimplePic
  {public:
    std::vector<PixType> img;
   int sx,sy;
   void Create(int tx,int ty)
   {
     sx=tx;sy=ty;
     img.resize(sx*sy);
   }
   PixType &Pix(int x, int y) {return img[sx*y+x];}

   void OpenGLSnap(GLenum format=0)
	{
		int vp[4];
		glGetIntegerv( GL_VIEWPORT,vp );		// Lettura viewport
		glPixelStorei( GL_PACK_ROW_LENGTH, 0);
		glPixelStorei( GL_PACK_ALIGNMENT, 1);
		int tx = vp[2];
		int ty = vp[3];

		Create(tx,ty);

		GLenum mtype  = 0;

		if(format==0) {
				format = GL_RGBA;
        mtype = GL_UNSIGNED_BYTE;
		}
		if(format==GL_DEPTH_COMPONENT) {
				format = GL_DEPTH_COMPONENT;
        mtype = GL_FLOAT;
		}
 		glReadPixels(vp[0],vp[1],vp[2],vp[3],format,mtype,(GLvoid *)&img[0]);
	}
	bool SavePPM( const char * filename )
	{
		FILE * fp = fopen(filename,"wb");
		if(fp==0) return false;


			fprintf(fp,"P6\n%d %d\n255\n",sx,sy);

			for(int i=0;i<sx*sy;++i)
			{
			 fwrite(&(img[i]),3,1,fp);
			}
	
		fclose(fp);
		return true;
	}

};

}
#endif // __VCG_MESH_VISIBILITY