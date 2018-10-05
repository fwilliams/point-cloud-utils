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
Revision 1.1  2006/12/10 19:55:09  ganovelli
first draft. Working but  ugly interface. right mouse of the button to place a pos, then prss buttons.

*/
#pragma once
/** the definition of vertex */
#include<vcg/simplex/vertex/base.h>						
/** the definition of face */
#include<vcg/simplex/face/base.h>
/** definition of triangle mesh */
#include<vcg/complex/complex.h>

/** allocation vertices and faces of triangle mesh */
#include<wrap/io_trimesh/import_PLY.h>


class DummyEdge;
class StraightFace;

/**************************************************************************************************************************/
/*    DEFINITION OF A VERY STRAIGHT MESH. No optional atributes, just normals in the vertices and flags in vertices and faces*/

/** definition of a very simple vertex type. Just coordinates and normal as attributes*/
class StraightVertex: public vcg::VertexSimp2< StraightVertex, DummyEdge, StraightFace, vcg::vert::Coord3f,vcg::vert::VFAdj,vcg::vert::Normal3f,vcg::vert::BitFlags>{};

/** definition of a very simple face type. Just color and reference to vertices as attribute*/
class StraightFace: public vcg::FaceSimp2< StraightVertex, DummyEdge, StraightFace,  vcg::	face::VertexRef,  vcg::	face::FFAdj,  vcg::	face::VFAdj,vcg::	face::Normal3f,vcg::face::BitFlags > {};

/** definition of a very simple mesh*/
class MyStraightMesh: public vcg::tri::TriMesh< std::vector<StraightVertex>,std::vector<StraightFace> >{};

/****************************************************************************************************************************/
