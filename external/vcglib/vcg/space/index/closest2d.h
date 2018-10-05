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


#ifndef __VCG_TRIMESH_CLOSEST_2D
#define __VCG_TRIMESH_CLOSEST_2D

#include <math.h>
#include <vcg/space/point2.h>
#include <vcg/space/segment2.h>
#include <vcg/space/box2.h>
#include <vcg/simplex/face/distance.h>

namespace vcg {
    namespace tri {

        //**MARKER CLASSES**//
        template <class MESH_TYPE,class OBJ_TYPE>
        class Tmark
        {
            MESH_TYPE *m;
        public:
            Tmark(){}
            Tmark(	MESH_TYPE *m) {SetMesh(m);}
            void UnMarkAll(){ vcg::tri::UnMarkAll(*m);}
            bool IsMarked(OBJ_TYPE* obj){return (vcg::tri::IsMarked(*m,obj));}
            void Mark(OBJ_TYPE* obj){ vcg::tri::Mark(*m,obj);}
            void SetMesh(MESH_TYPE *_m)
            {
              m=_m;
            }
        };

        template <class MESH, class GRID>
            typename MESH::FaceType * GetClosestEdgeBase( MESH & mesh,GRID & gr,const typename GRID::CoordType & _p,
                                                          const typename GRID::ScalarType _maxDist,typename GRID::ScalarType & _minDist,
                                                          typename GRID::CoordType &_closestPt)
        {
            typedef typename GRID::ScalarType ScalarType;
            typedef Point3<ScalarType> Point3x;
            typedef FaceTmark<MESH> MarkerFace;
            MarkerFace mf;
            mf.SetMesh(&mesh);
            vcg::PointSegment2DEPFunctor<ScalarType> PDistFunct;
            _minDist=_maxDist;
            return (gr.GetClosest(PDistFunct,mf,_p,_maxDist,_minDist,_closestPt));
        }

    }	 // end namespace tri
}	 // end namespace vcg

#endif
