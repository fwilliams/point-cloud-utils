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
#include <stdio.h>

#include<vcg/complex/complex.h>
#include<vcg/simplex/face/distance.h>
#include<vcg/simplex/face/component_ep.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/component_ep.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/space/intersection3.h>
#include <vcg/space/index/aabb_binary_tree/aabb_binary_tree.h>

typedef float AScalarType;

using namespace vcg;

class AVertex;
class AFace;

struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<AVertex>		::AsVertexType,
											vcg::Use<AFace>			::AsFaceType>{};

class AVertex     : public Vertex< MyUsedTypes, vertex::Normal3f, vertex::Coord3f,vertex::BitFlags >{};
class AFace       : public Face<   MyUsedTypes, face::VertexRef, face::Normal3f, face::EdgePlane, face::BitFlags> {};
class AMesh     : public vcg::tri::TriMesh< std::vector<AVertex>, std::vector<AFace> > { };

typedef vcg::AABBBinaryTreeIndex<AFace, AScalarType, vcg::EmptyClass> AIndex;

static AMesh gMesh;
static AIndex gIndex;

static void CreateMesh(void) {
	vcg::tri::Dodecahedron<AMesh>(gMesh);

	vcg::tri::UpdateFlags<AMesh>::Clear(gMesh);
	vcg::tri::UpdateNormal<AMesh>::PerVertexNormalized(gMesh);
	vcg::tri::UpdateComponentEP<AMesh>::Set(gMesh);
}

static void SetIndex(void) {
	gIndex.Set(gMesh.face.begin(), gMesh.face.end());
}

static void TestClosest(void) {
	vcg::face::PointDistanceEPFunctor<AIndex::ScalarType> getPtDist;
	const AIndex::CoordType queryPoint((AIndex::ScalarType)0, (AIndex::ScalarType)0, (AIndex::ScalarType)0);
	const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();

	AIndex::ObjPtr closestFace;
	AIndex::ScalarType closestDist;
	AIndex::CoordType closestPoint;

	vcg::EmptyClass a;
	closestFace = gIndex.GetClosest(getPtDist, a, queryPoint, maxDist, closestDist, closestPoint);

	printf("GetClosest Test:\n");

	if (closestFace != 0) {
		printf("\tface     : 0x%p\n", closestFace);
		printf("\tdistance : %f\n", closestDist);
		printf("\tpoint    : [%f, %f, %f]\n", closestPoint[0], closestPoint[1], closestPoint[2]);
	}
	else {
		printf("\tno object found (index is probably empty).\n");
	}
}

static void TestKClosest(void) {
	vcg::face::PointDistanceEPFunctor<AIndex::ScalarType> getPtDist;
	const unsigned int k = 10;
	const AIndex::CoordType queryPoint((AIndex::ScalarType)0, (AIndex::ScalarType)0, (AIndex::ScalarType)0);
	const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();

	std::vector<AIndex::ObjPtr> closestObjects;
	std::vector<AIndex::ScalarType> closestDistances;
	std::vector<AIndex::CoordType> closestPoints;

	vcg::EmptyClass a;
	unsigned int rk = gIndex.GetKClosest(getPtDist, a, k, queryPoint, maxDist, closestObjects, closestDistances, closestPoints);

	printf("GetKClosest Test:\n");
	printf("\tfound %d objects\n", rk);
}

static void TestRay(void) {
	const bool TEST_BACK_FACES = true;

	vcg::RayTriangleIntersectionFunctor<TEST_BACK_FACES> rayIntersector;
	const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();
	const AIndex::CoordType rayOrigin((AIndex::ScalarType)0, (AIndex::ScalarType)0, (AIndex::ScalarType)0);
	const AIndex::CoordType rayDirection((AIndex::ScalarType)1, (AIndex::ScalarType)0, (AIndex::ScalarType)0);
	const vcg::Ray3<AIndex::ScalarType, false> ray(rayOrigin, rayDirection);

	AIndex::ObjPtr isectFace;
	AIndex::ScalarType rayT;
	AIndex::CoordType isectPt;

	vcg::EmptyClass a;
	isectFace = gIndex.DoRay(rayIntersector, a , ray, maxDist, rayT);

	printf("DoRay Test:\n");
	if (isectFace != 0) {
		printf("\tface  : 0x%p\n", isectFace);
		printf("\tray t : %f\n", rayT);
	}
	else {
		printf("\tno object found (index is probably empty).\n");
	}
}

int main (int /*argc*/, char ** /*argv*/) {
	CreateMesh();

	SetIndex();

	printf("Spatial Index Tests\n");
	printf("---\n");
	TestClosest();
	printf("---\n");
	TestKClosest();
	printf("---\n");
	TestRay();
	printf("---\n");

	return (0);
}
