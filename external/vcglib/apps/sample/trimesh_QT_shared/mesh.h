#ifndef MESH_H
#define MESH_H

/// vcg imports
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/create/platonic.h>

#include <wrap/io_trimesh/import.h>

using namespace vcg;
class CFaceO;
class CVertexO;
class CEdgeO;

struct MyUsedTypes : public UsedTypes<Use<CVertexO>		::AsVertexType, vcg::Use<CEdgeO   >::AsEdgeType,
	Use<CFaceO>			::AsFaceType>{};

/// compositing wanted proprieties
class CVertexO : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags>{};
class CFaceO   : public vcg::Face<  MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::BitFlags > {};
class CEdgeO : public vcg::Edge<MyUsedTypes,vcg::edge::BitFlags,vcg::edge::EVAdj,vcg::edge::EEAdj>{};
class CMeshO   : public vcg::tri::TriMesh< std::vector<CVertexO>, std::vector<CFaceO>,std::vector<CEdgeO> > 
{
};

#endif