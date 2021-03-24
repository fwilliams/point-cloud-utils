//
// FIXME: This function is not publically exposed.
//
#include <npe.h>
#include <npe_typedefs.h>
#include "common.h"
#include <vector>
#include <tuple>


const char* ds_voxel_mesh = R"igl_Qu8mg5v7(
Given a 3D grid of voxels representing an indicator function
(<= 0 meaning outside and > 0 meaning inside), construct a
triangular mesh of the boundary cells

Parameters
----------
data: Flat list of voxel data
nx: Number of voxels in the x axis
ny: Number of voxels in the y axis
nz: Number of voxels in the z axis

Returns
-------
v : array of shape [n, 3] of n vertices for the mesh
f : array of shape [m, 3] of m faces for the mesh

See also
--------

Notes
-----

Examples
--------
)igl_Qu8mg5v7";

npe_function(voxel_mesh)
npe_doc(ds_voxel_mesh)
npe_arg(data, dense_float, dense_double)
npe_arg(nx, std::uint64_t)
npe_arg(ny, std::uint64_t)
npe_arg(nz, std::uint64_t)
npe_begin_code()
{
    if (data.rows() != 1 && data.cols() != 1) {
        throw pybind11::value_error("Invalid shape for data. Expected a flat array.");
    }

    if (nx*ny*nz != data.rows()*data.cols()) {
        throw pybind11::value_error("Invalid shape for data. Expected len(data) = nx * ny * nz.");
    }

    std::vector<Eigen::RowVector3d> vs;
    std::vector<Eigen::RowVector3i> fs;

    for (std::uint64_t x = 0; x < nx; x++) {
        for (std::uint64_t y = 0; y < ny; y++) {
            for (std::uint64_t z = 0; z < nz; z++) {
                std::uint64_t idx = x*nx*ny + y*ny + z;
                std::uint64_t idx_in = x*nx*ny + y*ny + (z+1);
                std::uint64_t idx_out = x*nx*ny + y*ny + (z-1);
                std::uint64_t idx_left = (x-1)*nx*ny + y*ny + z;
                std::uint64_t idx_right = (x+1)*nx*ny + y*ny + z;
                std::uint64_t idx_up = x*nx*ny + (y+1)*ny + z;
                std::uint64_t idx_down = x*nx*ny + (y-1)*ny + z;

                bool inside = data.data()[idx] > 0.0;

                typedef typename npe_Matrix_data::Scalar Scalar;
                Scalar val_in;
                Scalar val_out;
                Scalar val_left;
                Scalar val_right;
                Scalar val_up;
                Scalar val_down;

                if (inside) {
                    if (z == 0) {
                        val_out = 0.0;
                        val_in = data.data()[idx_in];
                    } else if (z == nz-1) {
                        val_in = 0.0;
                        val_out = data.data()[idx_out];
                    } else {
                        val_in = data.data()[idx_in];
                        val_out = data.data()[idx_out];
                    }

                    if (y == 0) {
                        val_down = 0.0;
                        val_up = data.data()[idx_up];
                    } else if (y == ny - 1) {
                        val_up = data.data()[idx_up];
                        val_down = 0.0;
                    } else  {
                        val_up = data.data()[idx_up];
                        val_down = data.data()[idx_down];
                    }

                    if (x == 0) {
                        val_left = 0.0;
                        val_right = data.data()[idx_right];
                    } else if (x == nx - 1) {
                        val_right = 0.0;
                        val_left = data.data()[idx_left];
                    } else {
                        val_left = data.data()[idx_left];
                        val_right = data.data()[idx_right];
                    }

                    Eigen::RowVector3d ctr((double)x, (double)y, (double)z);
                    Eigen::RowVector3d in(0.0, 0.0, 0.5);
                    Eigen::RowVector3d right(0.5, 0.0, 0.0);
                    Eigen::RowVector3d up(0.0, 0.5, 0.0);

                    Eigen::RowVector3d xp_yp_zp = ctr + right + up + in;
                    Eigen::RowVector3d xn_yp_zp = ctr - right + up + in;
                    Eigen::RowVector3d xp_yn_zp = ctr + right - up + in;
                    Eigen::RowVector3d xp_yp_zn = ctr + right + up - in;
                    Eigen::RowVector3d xn_yn_zp = ctr - right - up + in;
                    Eigen::RowVector3d xn_yp_zn = ctr - right + up - in;
                    Eigen::RowVector3d xp_yn_zn = ctr + right - up - in;
                    Eigen::RowVector3d xn_yn_zn = ctr - right - up - in;


                    if (val_right <= 0.0) {
                        int vsize = vs.size();
                        vs.push_back(xp_yn_zn);
                        vs.push_back(xp_yn_zp);
                        vs.push_back(xp_yp_zp);
                        vs.push_back(xp_yp_zn);
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+2, vsize+1));
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+3, vsize+2));
                    }
                    if (val_left <= 0.0) {
                        int vsize = vs.size();
                        vs.push_back(xn_yn_zn);
                        vs.push_back(xn_yn_zp);
                        vs.push_back(xn_yp_zp);
                        vs.push_back(xn_yp_zn);
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+1, vsize+2));
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+2, vsize+3));
                    }

                    if (val_in <= 0.0) {
                        int vsize = vs.size();
                        vs.push_back(xn_yn_zp);
                        vs.push_back(xn_yp_zp);
                        vs.push_back(xp_yp_zp);
                        vs.push_back(xp_yn_zp);
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+2, vsize+1));
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+3, vsize+2));
                    }
                    if (val_out <= 0.0) {
                        int vsize = vs.size();
                        vs.push_back(xn_yn_zn);
                        vs.push_back(xn_yp_zn);
                        vs.push_back(xp_yp_zn);
                        vs.push_back(xp_yn_zn);
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+1, vsize+2));
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+2, vsize+3));
                    }

                    if (val_up <= 0.0) {
                        int vsize = vs.size();
                        vs.push_back(xn_yp_zn);
                        vs.push_back(xn_yp_zp);
                        vs.push_back(xp_yp_zp);
                        vs.push_back(xp_yp_zn);
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+1, vsize+2));
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+2, vsize+3));
                    }

                    if (val_down <= 0.0) {
                        int vsize = vs.size();
                        vs.push_back(xn_yn_zn);
                        vs.push_back(xn_yn_zp);
                        vs.push_back(xp_yn_zp);
                        vs.push_back(xp_yn_zn);
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+2, vsize+1));
                        fs.push_back(Eigen::RowVector3i(vsize+0, vsize+3, vsize+2));
                    }

                }
            }
        }
    }

    EigenDenseF64 retv(vs.size(), 3);
    for (int i = 0; i < vs.size(); i++) {
        retv(i, 0) = vs[i][0];
        retv(i, 1) = vs[i][1];
        retv(i, 2) = vs[i][2];
    }

    EigenDenseF64 retf(fs.size(), 3);
    for (int i = 0; i < fs.size(); i++) {
        retf(i, 0) = fs[i][0];
        retf(i, 1) = fs[i][1];
        retf(i, 2) = fs[i][2];
    }

    return std::make_tuple(npe::move(retv), npe::move(retf));
}
npe_end_code()
