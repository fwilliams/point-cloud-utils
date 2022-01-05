#include <npe.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <string>

#include "vcg_utils.h"
//#include "mls_utils/apss.h"
#include "common.h"


namespace {

using namespace vcg;
class VCGMeshEdge;
class VCGMeshFace;
class VCGMeshVertex;
struct VCGMeshUsedTypes : public UsedTypes<	Use<VCGMeshVertex>   ::AsVertexType,
        Use<VCGMeshEdge>     ::AsEdgeType,
        Use<VCGMeshFace>     ::AsFaceType>{};
class VCGMeshVertex  : public Vertex<VCGMeshUsedTypes,
                                     vertex::VFAdj, vertex::Coord3d, vertex::Normal3d, vertex::Radiusd,
                                     vertex::BitFlags, vertex::CurvatureDird, vertex::Curvatured> {};
class VCGMeshFace    : public Face<VCGMeshUsedTypes,
                                   face::VFAdj, face::FFAdj,
                                   face::Normal3d, face::VertexRef, face::BitFlags> {};
class VCGMeshEdge    : public Edge<VCGMeshUsedTypes>{};
class VCGMesh : public tri::TriMesh<std::vector<VCGMeshVertex>, std::vector<VCGMeshFace>, std::vector<VCGMeshEdge>> {};

}


const char* mesh_principal_curvatures_doc = R"Qu8mg5v7(
Estimate principal curvature directions and magnitudes for a mesh

Parameters
----------
v : #v by 3 Matrix of mesh vertex 3D positions
f : #f by 3 Matrix of face (triangle) indices
r : optional floating point radius of neighborhood to consider when estimating curvature

Returns
-------
A tuple (k1, k2, d1, d2) where:
  k1 is an array of shape (#v,) of maximum curvature magnitudes
  k2 is an array of shape (#v,) of minimum curvature magnitudes
  d1 is an array of shape (#v, 3) of maximum curvature directions
  d2 is an array of shape (#v, 3) of minimum curvature directions
)Qu8mg5v7";
npe_function(mesh_principal_curvatures)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_long, dense_longlong)
npe_default_arg(r, double, -1.0)
npe_doc(mesh_principal_curvatures_doc)
npe_begin_code()
{
    VCGMesh m;
    vcg_mesh_from_vf(v, f, m);
	tri::UpdateTopology<VCGMesh>::FaceFace(m);
	tri::UpdateTopology<VCGMesh>::VertexFace(m);
	if (r > 0.0) {
	    tri::UpdateCurvature<VCGMesh>::PrincipalDirectionsPCA(m, r);
	} else {
	    tri::UpdateCurvature<VCGMesh>::PrincipalDirections(m);
	}

    npe_Matrix_v k1(m.vn, 1), k2(m.vn, 1), d1(m.vn, 3), d2(m.vn, 3);

    int vcount = 0;
    for (VCGMesh::VertexIterator vit = m.vert.begin(); vit != m.vert.end(); vit++) {
        k1(vcount, 0) = vit->cK1();
        k2(vcount, 0) = vit->cK2();
        for (int i = 0; i < 3; i++) {
            d1(vcount, i) = vit->cPD1()[i];
            d2(vcount, i) = vit->cPD2()[i];
        }
        vcount += 1;
    }

    return std::make_tuple(npe::move(k1), npe::move(k2), npe::move(d1), npe::move(d2));
}
npe_end_code()


const char* mesh_mean_and_gaussian_curvatures_doc = R"Qu8mg5v7(
Estimate mean and Gaussian curvatures for a mesh

Parameters
----------
v : #v by 3 Matrix of mesh vertex 3D positions
f : #f by 3 Matrix of face (triangle) indices

Returns
-------
A tuple (kg, kh, d1, d2) where:
  kg is an array of shape (#v,) of per-vertex Gaussian curvatures
  kh is an array of shape (#v,) of per-vertex mean curvatures
)Qu8mg5v7";
npe_function(mesh_mean_and_gaussian_curvatures)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_long, dense_longlong)
npe_doc(mesh_mean_and_gaussian_curvatures_doc)
npe_begin_code()
{
    VCGMesh m;
    vcg_mesh_from_vf(v, f, m);
	tri::UpdateTopology<VCGMesh>::FaceFace(m);
	tri::UpdateTopology<VCGMesh>::VertexFace(m);
    tri::UpdateCurvature<VCGMesh>::MeanAndGaussian(m);

    npe_Matrix_v kg(m.vn, 1), kh(m.vn, 1);

    int vcount = 0;
    for (VCGMesh::VertexIterator vit = m.vert.begin(); vit != m.vert.end(); vit++) {
        kg(vcount, 0) = vit->cKg();
        kh(vcount, 0) = vit->cKh();
        vcount += 1;
    }

    return std::make_tuple(npe::move(kg), npe::move(kh));
}
npe_end_code()


//const char* pointcloud_apss_curvature_doc = R"Qu8mg5v7(
//)Qu8mg5v7";
//npe_function(pointcloud_apss_curvature)
//npe_arg(v, dense_float, dense_double)
//npe_arg(f, dense_int, dense_long, dense_longlong)
//npe_doc(pointcloud_apss_curvature_doc)
//npe_begin_code()
//{
//    VCGMesh m;
//    vcg_mesh_from_vf(v, f, m);
//	tri::UpdateTopology<VCGMesh>::FaceFace(m);
//	tri::UpdateTopology<VCGMesh>::VertexFace(m);
////    vcg_mesh_from_v(v, m);
//	GaelMls::APSS<VCGMesh> apss(m);
//
//    npe_Matrix_v r(m.vn, 1);
//
//    int vcount = 0;
//    for (VCGMesh::VertexIterator vit = m.vert.begin(); vit != m.vert.end(); vit++) {
//        auto p = apss.project(vit->cP());
//        r(vcount, 0) = apss.approxMeanCurvature(p);
//        std::cout << p[0] << ", " << p[1] << ", " << p[2] << " -- " << r(vcount, 0) << std::endl;
//        vcount += 1;
//    }
//
//    return npe::move(r);
//}
//npe_end_code()