#pragma once

#include <Eigen/Core>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>

/*
 * Copy a mesh stored as a #V x 3 matrix of vertices, V, and a #F x 3 matrix of face indices into a VCG mesh
 */
template <typename DerivedV, typename DerivedF, typename DerivedN, typename VCGMesh>
static void vcg_mesh_from_vfn(const Eigen::MatrixBase<DerivedV>& V,
                              const Eigen::MatrixBase<DerivedF>& F,
                              const Eigen::MatrixBase<DerivedN>& N, VCGMesh& m) {
    using namespace vcg;
    typename VCGMesh::VertexIterator vit = Allocator<VCGMesh>::AddVertices(m, V.rows());
    std::vector<typename VCGMesh::VertexPointer> ivp(V.rows());
    for (int i = 0; i < V.rows(); i++) {
        ivp[i] = &*vit;
        vit->P() = typename VCGMesh::CoordType(V(i, 0), V(i, 1), V(i, 2));
        if (N.rows() > 0) {
            vit->N() = typename VCGMesh::CoordType(N(i, 0), N(i, 1), N(i, 2));
        }
        vit++;
    }

    if (F.rows() > 0) {
        typename VCGMesh::FaceIterator fit = Allocator<VCGMesh>::AddFaces(m, F.rows());
        for (int i = 0; i < F.rows(); i++) {
            fit->V(0) = ivp[F(i, 0)];
            fit->V(1) = ivp[F(i, 1)];
            fit->V(2) = ivp[F(i, 2)];
            fit++;
        }
    }

    tri::UpdateBounding<VCGMesh>::Box(m);
}

/*
 * Copy a mesh stored as a #V x 3 matrix of vertices, V, and a #F x 3 matrix of face indices into a VCG mesh
 */
template <typename DerivedV, typename DerivedF, typename VCGMesh>
static void vcg_mesh_from_vf(const Eigen::MatrixBase<DerivedV>& V,
                             const Eigen::MatrixBase<DerivedF>& F,
                             VCGMesh& m) {
    Eigen::MatrixXd N(0, 3);
    vcg_mesh_from_vfn(V, F, N, m);
}


/*
 * Copy a mesh stored as a #V x 3 matrix of vertices, V, and a #F x 3 matrix of face indices into a VCG mesh
 */
template <typename DerivedV, typename VCGMesh>
static void vcg_mesh_from_v(const Eigen::MatrixBase<DerivedV>& V,
                            VCGMesh& m) {
    Eigen::MatrixXi F(0, 3);
    Eigen::MatrixXd N(0, 3);
    vcg_mesh_from_vfn(V, F, N, m);
}


/*
 * Use this to sample vertex indices when we are sampling from a point cloud
 */
template<class MeshType>
class EigenVertexIndexSampler {
public:
    typedef typename MeshType::VertexType VertexType;

    // The mesh we are sampling from
    MeshType &sampled_mesh;

    // Indices into the mesh vertex array, this is an eigen matrix of some type
    typedef Eigen::Matrix<std::ptrdiff_t, Eigen::Dynamic, 1> IndexArray;
    IndexArray &indices;

    // Number of vertices
    int vcount = 0;

    EigenVertexIndexSampler(MeshType &in_mesh, IndexArray &out_inds) :
            sampled_mesh(in_mesh), indices(out_inds) {
    }

    void trim() {
        indices.conservativeResize(vcount, 1);
    }

    void reset() {
        vcount = 0;
    }

    void maybe_resize() {
        // If we are about to overflow indexes, double its size
        if (indices.size() <= vcount) {
            const int n_rows = indices.size() == 0 ? 1024 : indices.size();
            indices.conservativeResize(2 * n_rows, 1);
        }
    }

    void AddVert(const VertexType &p) {
        maybe_resize();
        std::ptrdiff_t p_offset = &p - &*sampled_mesh.vert.begin();
        indices(vcount, 0) = p_offset;
        vcount += 1;
    }
};