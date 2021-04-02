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