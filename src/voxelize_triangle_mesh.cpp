#include <npe.h>

#include "common/tribox.h"
#include "common/common.h"

namespace {

template <typename VertexMat, typename FaceMat>
void voxelize_triangle_mesh(const VertexMat& v, const FaceMat& f,
                            Eigen::Vector3d voxel_size,
                            Eigen::Vector3d voxel_origin,
                            Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& out_vox) {

    const Eigen::Vector3d box_half_size(voxel_size[0] / 2, voxel_size[1] / 2, voxel_size[2] / 2);

    std::vector<int> voxel_idxs;
    voxel_idxs.reserve(3*4*f.rows());
    for(int fi = 0; fi < f.rows(); fi += 1) {
        Eigen::Vector3d v0(v(f(fi, 0), 0), v(f(fi, 0), 1), v(f(fi, 0), 2));
        Eigen::Vector3d v1(v(f(fi, 1), 0), v(f(fi, 1), 1), v(f(fi, 1), 2));
        Eigen::Vector3d v2(v(f(fi, 2), 0), v(f(fi, 2), 1), v(f(fi, 2), 2));

        double minx, miny, minz, maxx, maxy, maxz;
        int inix, iniy, iniz, inox, inoy, inoz;
        minx = std::min(v0[0], std::min(v1[0], v2[0]));
        miny = std::min(v0[1], std::min(v1[1], v2[1]));
        minz = std::min(v0[2], std::min(v1[2], v2[2]));
        maxx = std::max(v0[0], std::max(v1[0], v2[0]));
        maxy = std::max(v0[1], std::max(v1[1], v2[1]));
        maxz = std::max(v0[2], std::max(v1[2], v2[2]));

        inix = static_cast<int>(std::floor((minx - voxel_origin[0]) / voxel_size[0]));
        iniy = static_cast<int>(std::floor((miny - voxel_origin[1]) / voxel_size[1]));
        iniz = static_cast<int>(std::floor((minz - voxel_origin[2]) / voxel_size[2]));

        inox = static_cast<int>(std::ceil((maxx - voxel_origin[0]) / voxel_size[0]));
        inoy = static_cast<int>(std::ceil((maxy - voxel_origin[1]) / voxel_size[1]));
        inoz = static_cast<int>(std::ceil((maxz - voxel_origin[2]) / voxel_size[2]));

        for (int widx = inix; widx <= inix; widx++) {
            for (int hidx = iniy; hidx <= inoy; hidx++) {
                for (int didx = iniz; didx <= inoz; didx++) {
                    const Eigen::Vector3d vox_pos = Eigen::Vector3d(widx, hidx, didx).array() * voxel_size.array();
                    const Eigen::Vector3d box_center = (voxel_origin + vox_pos).eval();

                    double* tri_verts[3] = {const_cast<double*>(v0.data()),
                                            const_cast<double*>(v1.data()),
                                            const_cast<double*>(v2.data())};

                    if (triBoxOverlap((double*) box_center.data(), (double*) box_half_size.data(), tri_verts)) {
                        voxel_idxs.push_back(widx);
                        voxel_idxs.push_back(hidx);
                        voxel_idxs.push_back(didx);
                        // break;
                    }
                }
            }
        }
    }

    out_vox.resize(voxel_idxs.size()/3, 3);
    memcpy(out_vox.data(), voxel_idxs.data(), voxel_idxs.size() * sizeof(int));
}

}


npe_function(_voxelize_triangle_mesh_internal)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64, dense_uint32, dense_uint64)
npe_arg(vox_origin, dense_double)
npe_arg(vox_size, npe_matches(vox_origin))
npe_begin_code()
    using MatrixI = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    validate_mesh(v, f);

    if (vox_origin.size() != 3) {
        throw pybind11::value_error("Invalid shape");
    }

    if (vox_size.size() != 3) {
        throw pybind11::value_error("Invalid shape");
    }

    Eigen::Vector3d vorgn = vox_origin;
    Eigen::Vector3d vsize = vox_size;

    if (vsize[0] <= 0.0 || vsize[1] <= 0.0 || vsize[2] <= 0.0) {
        throw pybind11::value_error("Invalid voxel size");
    }

    MatrixI out_vox;

    voxelize_triangle_mesh(v, f, vsize, vorgn, out_vox);

    return npe::move(out_vox);

npe_end_code()