/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2014                                                \/)\/    *
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

/* Optimizes given UV-mapping with
 * [Least Square Conformal Maps]
 * (minimizes angle distrotions).
 *
 * Needs:
 * (-) per-vertex texture coords
 * (-) some fixed boundary:
 *     Fixed vertices are the flagged ones.
 *     By default: (use fixedMask parameter to customize)
 *     BORDER or SELECTED verts are fixed.
 *
 * Example of usage:
 *    MyMesh m;
 *    vcg::tri::UpdateFlags< MyMesh >::VertexBorderFromNone( m );
 *    vcg::tri::OptimizeUV_LSCM( m );
 *
 */

#ifndef __VCG_IGL_LEAST_SQUARES_CONFORMAL_MAPS
#define __VCG_IGL_LEAST_SQUARES_CONFORMAL_MAPS

#include <igl/lscm.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>

namespace vcg{
namespace tri{

template<class MeshType >
void OptimizeUV_LSCM( MeshType& m ,
                      unsigned int fixedMask =
                        MeshType::VertexType::BORDER |
                        MeshType::VertexType::SELECTED
                     )
{

    // check requirements
    vcg::tri::VertexVectorHasPerVertexTexCoord( m.vert );
    vcg::tri::VertexVectorHasPerVertexFlags( m.vert );

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXi b;
    Eigen::MatrixXd bc;
    Eigen::MatrixXd V_uv;

    vcg::tri::MeshToMatrix< MeshType >::GetTriMeshData( m, F, V );

    // build fixed points data
    int nFixed = 0;
    for (int i=0; i<(int)m.vert.size(); i++) {
        if (m.vert[i].Flags()&fixedMask) nFixed++;
    }

    // all fixed, nothing to do? get out to avoid crashes
    if (nFixed == m.vert.size()) return;

    b.resize(nFixed);
    bc.resize(nFixed,2);
    for (int i=0,k=0; i<(int)m.vert.size(); i++) {
        if (m.vert[i].Flags()&fixedMask)  {
            b(k) = i;
            bc(k,0) = m.vert[i].T().P()[0];
            bc(k,1) = m.vert[i].T().P()[1];
            k++;
        }
    }

    // apply Least Square Conformal Maps
    ::igl::lscm(V,F,b,bc,V_uv);

    // copy results back to mesh
    for (int i=0; i<(int)m.vert.size(); i++) {
        m.vert[i].T().P()[0] = V_uv(i,0);
        m.vert[i].T().P()[1] = V_uv(i,1);
    }
}

}} // namespaces
#endif
