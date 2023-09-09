#include <npe.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/smooth.h>
#include <igl/per_vertex_normals.h>
#include <string>

#include "common/vcg_utils.h"
#include "common/common.h"


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


const char* laplacian_smooth_mesh_doc = R"igl_Qu8mg5v7(
Smooth a mesh using Laplacian smoothing

Args:
    v : \#v by 3 Matrix of mesh vertex 3D positions
    f : \#f by 3 Matrix of face (triangle) indices
    num_iters : Number of smoothing iterations
    use_cotan_weights : Whether to use cotagent weighting (False by default)

Returns:
    n : list of vertex normals of shape #v by 3

)igl_Qu8mg5v7";
npe_function(laplacian_smooth_mesh)
npe_doc(laplacian_smooth_mesh_doc)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64, dense_uint32, dense_uint64)
npe_arg(num_iters, int)
npe_default_arg(use_cotan_weights, bool, false)
npe_begin_code()
{
    validate_mesh(v, f);
    if (num_iters < 0) {
        throw pybind11::value_error("Invalid value for argument num_iters in smooth_mesh_laplacian. Must be a positive integer or 0.");
    }
    if (num_iters == 0) {  // No-op for 0
        return npe::move(v);
    }

    VCGMesh m;
    vcg_mesh_from_vf(v, f, m);

    vcg::tri::Smooth<VCGMesh>::VertexCoordLaplacian(m, num_iters, false /* SmoothSelected */, use_cotan_weights);

    npe_Matrix_v vsmooth(m.vn, 3);
    int vcount = 0;
    for (VCGMesh::VertexIterator vit = m.vert.begin(); vit != m.vert.end(); vit++) {
        vsmooth(vcount, 0) = vit->P()[0];
        vsmooth(vcount, 1) = vit->P()[1];
        vsmooth(vcount, 2) = vit->P()[2];
        vcount += 1;
    }

    return npe::move(vsmooth);
}
npe_end_code()
