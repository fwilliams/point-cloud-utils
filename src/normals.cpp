#include <npe.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/pointcloud_normal.h>
#include <igl/per_vertex_normals.h>
#include <string>

#include "vcg_utils.h"
#include "common.h"


namespace {

using namespace vcg;
class VCGMeshEdge;
class VCGMeshFace;
class VCGMeshVertex;
struct VCGMeshUsedTypes : public UsedTypes<	Use<VCGMeshVertex>   ::AsVertexType,
        Use<VCGMeshEdge>     ::AsEdgeType,
        Use<VCGMeshFace>     ::AsFaceType>{};
class VCGMeshVertex  : public Vertex<VCGMeshUsedTypes, vertex::Coord3d, vertex::Normal3d, vertex::BitFlags> {};
class VCGMeshFace    : public Face<VCGMeshUsedTypes, face::FFAdj,  face::Normal3d, face::VertexRef, face::BitFlags> {};
class VCGMeshEdge    : public Edge<VCGMeshUsedTypes>{};
class VCGMesh : public tri::TriMesh<std::vector<VCGMeshVertex>, std::vector<VCGMeshFace>, std::vector<VCGMeshEdge>> {};

}


const char* estimate_point_cloud_normals_doc = R"Qu8mg5v7(
Estimate normals for a point cloud by locally fitting a plane to a small neighborhood of points

Parameters
----------
v : #v by 3 array of vertex positions (each row is a vertex)
k : Number of nearest neighbors to use in the estimate for the normal of a point. Default: 10.
smoothing_iterations : Number of smoothing iterations to apply to the estimated normals. Default: 0.

Returns
-------
A #v x 3 array of normals where each row i is the normal at vertex v[i]

)Qu8mg5v7";
npe_function(estimate_point_cloud_normals)
npe_arg(v, dense_float, dense_double)
npe_default_arg(k, int, 10)
npe_default_arg(smoothing_iterations, int, 0)
npe_doc(estimate_point_cloud_normals_doc)
npe_begin_code()
{
    VCGMesh m;
    vcg_mesh_from_v(v, m);

    tri::PointCloudNormal<VCGMesh>::Param p;
    p.fittingAdjNum = k;
    p.smoothingIterNum = smoothing_iterations;
    p.viewPoint = vcg::Point3d(0.0f, 0.0f, 0.0f);
    p.useViewPoint = false;

    tri::PointCloudNormal<VCGMesh>::Compute(m, p, (vcg::CallBackPos*) nullptr);

    npe_Matrix_v ret(m.vn, 3);

    int vcount = 0;
    for (VCGMesh::VertexIterator vit = m.vert.begin(); vit != m.vert.end(); vit++) {
        ret(vcount, 0) = vit->N()[0];
        ret(vcount, 1) = vit->N()[1];
        ret(vcount, 2) = vit->N()[2];
        vcount += 1;
    }

    return npe::move(ret);

}
npe_end_code()


const char* estimate_mesh_normals_doc = R"igl_Qu8mg5v7(
Compute vertex normals of a mesh from its vertices and faces using face area weighting

Parameters
----------
v : #v by 3 Matrix of mesh vertex 3D positions
f : #f by 3 Matrix of face (triangle) indices
weighting_type : Weighting type must be one of 'uniform', 'angle', or 'area' (default is 'uniform')

Returns
-------
n : list of vertex normals of shape #v by 3

See also
--------

Notes
-----

Examples
--------
)igl_Qu8mg5v7";
npe_function(estimate_mesh_normals)
npe_doc(estimate_mesh_normals_doc)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_default_arg(weighting_type, std::string, std::string("uniform"))
npe_begin_code()
{
    validate_mesh(v, f);
    int weighting = -1;
    if (weighting_type == "uniform") {
        weighting  = 0;
    } else if (weighting_type == "area") {
        weighting = 1;
    } else if (weighting_type == "angle") {
        weighting = 2;
    } else {
        throw pybind11::value_error("Invalid weighting_type must be one of 'uniform', 'angle', or 'area', but got '" +
                                    weighting_type + "'");
    }

    npe_Matrix_v n;
    igl::PerVertexNormalsWeightingType wtype = igl::PerVertexNormalsWeightingType(weighting);
    igl::per_vertex_normals(v, f, wtype, n);

    return npe::move(n);
}
npe_end_code()
