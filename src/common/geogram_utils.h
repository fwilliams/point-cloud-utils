#include <Eigen/Core>
#include <geogram/basic/geometry.h>
#include <geogram/basic/common.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/RVD.h>

#include <mutex>

#ifndef GEOGRAM_UTILS_H
#define GEOGRAM_UTILS_H

extern bool geogram_is_initialized;
extern std::mutex geogram_init_mutex;

/*
 * Need to call this at the beginning of any function using Geogram
 */
void init_geogram_only_once();

/*
 * Convert a (V, F) pair to a geogram mesh
 */
template <typename DerivedV, typename DerivedF>
void vf_to_geogram_mesh(const Eigen::MatrixBase<DerivedV> &V, const Eigen::MatrixBase<DerivedF> &F, GEO::Mesh &M) {
    M.clear();
    // Setup vertices
    M.vertices.create_vertices((int) V.rows());
    for (int i = 0; i < (int) M.vertices.nb(); ++i) {
        GEO::vec3 &p = M.vertices.point(i);
        p[0] = V(i, 0);
        p[1] = V(i, 1);
        p[2] = (V.cols() == 2 ? 0 : V(i, 2));
    }
    // Setup faces
    if (F.cols() == 3) {
        M.facets.create_triangles((int) F.rows());
    } else if (F.cols() == 4) {
        M.facets.create_quads((int) F.rows());
    } else {
        throw std::runtime_error("Mesh face type not supported");
    }
    for (int c = 0; c < (int) M.facets.nb(); ++c) {
        for (int lv = 0; lv < F.cols(); ++lv) {
            M.facets.set_vertex(c, lv, F(c, lv));
        }
    }
    M.facets.connect();
}

/*
 * Convert a (V, T) pair to a geogram mesh
 */
template <typename DerivedV, typename DerivedF, typename DerivedT>
void vft_to_geogram_tet_mesh(const Eigen::MatrixBase<DerivedV> &V, const Eigen::MatrixBase<DerivedF> &F, const Eigen::MatrixBase<DerivedT> &T, GEO::Mesh &M) {
    M.clear();

    // Setup vertices
    M.vertices.create_vertices((int) V.rows());
    for (int i = 0; i < (int) M.vertices.nb(); ++i) {
        GEO::vec3 &p = M.vertices.point(i);
        p[0] = V(i, 0);
        p[1] = V(i, 1);
        p[2] = (V.cols() == 2 ? 0 : V(i, 2));
    }

    // Setup faces
    if (F.cols() == 3) {
        M.facets.create_triangles((int) F.rows());
    } else {
        throw std::runtime_error("Mesh face type not supported");
    }
    for (int c = 0; c < (int) M.facets.nb(); ++c) {
        for (int lv = 0; lv < F.cols(); ++lv) {
            M.facets.set_vertex(c, lv, F(c, lv));
        }
    }

    // Setup tets
    if (T.cols() == 4) {
        M.cells.create_tets((int) T.rows());
    } else {
        throw std::runtime_error("Mesh cell type not supported");
    }
    for (int c = 0; c < (int) M.cells.nb(); ++c) {
        for (int lv = 0; lv < T.cols(); ++lv) {
            M.cells.set_vertex(c, lv, T(c, lv));
        }
    }

    M.cells.connect();
    M.facets.connect();
}

#endif // GEOGRAM_UTILS_H
