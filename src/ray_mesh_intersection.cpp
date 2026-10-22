#include <npe.h>
#include <common/embree_intersector.h>
#include <tuple>
#include <numeric>

#include "common/common.h"

namespace py = pybind11;

using Intersector = EmbreeIntersector;
void hack_extra_ray_mesh_bindings(pybind11::module& m) {
    py::class_<Intersector, std::shared_ptr<Intersector>>(m, "_RayMeshIntersectorInternal")
    .def(py::init([]() {
        return std::shared_ptr<Intersector>(new EmbreeIntersector());
    }));
}


npe_function(_populate_ray_intersector_internal)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64, dense_uint32, dense_uint64)
npe_arg(isector, std::shared_ptr<EmbreeIntersector>)
npe_begin_code()
    EmbreeIntersector::PointMatrixType v_copy = v.template cast<float>();
    EmbreeIntersector::FaceMatrixType f_copy = f.template cast<int>();
    isector->init(v_copy, f_copy, true /*is_static*/);
npe_end_code()


npe_function(_intersect_ray_intersector_internal)
npe_arg(ray_o, dense_float, dense_double)
npe_arg(ray_d, npe_matches(ray_o))
npe_arg(isector, std::shared_ptr<Intersector>)
npe_default_arg(ray_near, double, 0.0)
npe_default_arg(ray_far, double, std::numeric_limits<double>::infinity())
npe_begin_code()
    bool use_single_ray_origin = false;
    if (ray_o.size() == 3) {
        use_single_ray_origin = true;
    } else {
        if (ray_o.rows() != ray_d.rows()) {
            throw pybind11::value_error("ray_o and ray_d must have the same number of rows (one ray origin per ray direction). "
                                        "(Note: ray_o can have one row to use the same origin for all directions)");
        }
    }

    if (ray_o.cols() != 3 && !use_single_ray_origin) {
        throw pybind11::value_error("Invalid shape for ray_o must have shape (N, 3) but got (" +
                                    std::to_string(ray_o.rows()) + ", " + std::to_string(ray_o.cols()) + ").");
    }
    if (ray_d.cols() != 3) {
        throw pybind11::value_error("Invalid shape for ray_d must have shape (N, 3) but got (" +
                                    std::to_string(ray_d.rows()) + ", " + std::to_string(ray_d.cols()) + ").");
    }

    EigenDense<int> ret_fid(ray_d.rows(), 1);
    npe_Matrix_ray_o ret_bc(ray_d.rows(), 3);
    npe_Matrix_ray_o ret_t(ray_d.rows(), 1);
    for (int i = 0; i < ray_d.rows(); i++) {
        Eigen::RowVector3f o_i;
        if (use_single_ray_origin) {
            o_i = Eigen::RowVector3f((float)ray_o(0, 0), (float)ray_o(1, 0), (float)ray_o(2, 0));
        } else {
            o_i = Eigen::RowVector3f((float)ray_o(i, 0), (float)ray_o(i, 1), (float)ray_o(i, 2));
        }
        Eigen::RowVector3f d_i((float)ray_d(i, 0), (float)ray_d(i, 1), (float)ray_d(i, 2));
        Hit hit;

        bool is_hit = isector->intersectRay(o_i, d_i, hit, ray_near, ray_far);
        if (is_hit) {
            ret_fid(i, 0) = (int) hit.id;
            ret_bc(i, 0) = 1.0 - hit.u - hit.v;
            ret_bc(i, 1) = hit.u;
            ret_bc(i, 2) = hit.v;
            ret_t(i, 0) = hit.t;
        } else {
            ret_fid(i, 0) = -1;
            ret_t(i, 0) = std::numeric_limits<npe_Scalar_ray_o>::infinity();
            ret_bc(i, 0) = 0.0;
            ret_bc(i, 1) = 0.0;
            ret_bc(i, 2) = 0.0;
        }

    }

    return std::make_tuple(npe::move(ret_fid), npe::move(ret_bc), npe::move(ret_t));
npe_end_code()



const char* ray_mesh_intersection_doc = R"Qu8mg5v7(
Compute intersection between a set of rays and a triangle mesh

Args:
    v : \#v by 3 array of vertex positions (each row is a vertex)
    f : \#f by 3 Matrix of face (triangle) indices
    ray_o : array of shape (#rays, 3) of ray origins (one per row) or a single array of shape (3,) to use
    ray_d : array of shape (#rays, 3) of ray directions (one per row)
    ray_near : an optional floating point value indicating the distance along each ray to start searching (default 0.0)
    ray_far : an optional floating point value indicating the maximum distance along each ray to search (default inf)

Returns:
    f_id : an array of shape (#rays,) representing the face id hit by each ray
    bc : an array of shape (#rays, 3) where each row is the barycentric coordinates within each face of the ray intersection
    t : the distance along each ray to the intersection
)Qu8mg5v7";
npe_function(ray_mesh_intersection)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64, dense_uint32, dense_uint64)
npe_arg(ray_o, npe_matches(v))
npe_arg(ray_d, npe_matches(v))
npe_default_arg(ray_near, double, 0.0)
npe_default_arg(ray_far, double, std::numeric_limits<double>::infinity())
npe_doc(ray_mesh_intersection_doc)
npe_begin_code()
{
    bool use_single_ray_origin = false;
    if (ray_o.size() == 3) {
        use_single_ray_origin = true;
    } else {
        if (ray_o.rows() != ray_d.rows()) {
            throw pybind11::value_error("ray_o and ray_d must have the same number of rows (one ray origin per ray direction). "
                                        "(Note: ray_o can have one row to use the same origin for all directions)");
        }
    }

    if (ray_o.cols() != 3 && !use_single_ray_origin) {
        throw pybind11::value_error("Invalid shape for ray_o must have shape (N, 3) but got (" +
                                    std::to_string(ray_o.rows()) + ", " + std::to_string(ray_o.cols()) + ").");
    }
    if (ray_d.cols() != 3) {
        throw pybind11::value_error("Invalid shape for ray_d must have shape (N, 3) but got (" +
                                    std::to_string(ray_d.rows()) + ", " + std::to_string(ray_d.cols()) + ").");
    }
    validate_mesh(v, f);


    EmbreeIntersector isector;

    EmbreeIntersector::PointMatrixType v_copy = v.template cast<float>();
    EmbreeIntersector::FaceMatrixType f_copy = f.template cast<int>();
    isector.init(v_copy, f_copy, true /*is_static*/);

    npe_Matrix_f ret_fid(ray_d.rows(), 1);
    npe_Matrix_v ret_bc(ray_d.rows(), 3);
    npe_Matrix_v ret_t(ray_d.rows(), 1);
    for (int i = 0; i < ray_d.rows(); i++) {
        Eigen::RowVector3f o_i;
        if (use_single_ray_origin) {
            o_i = Eigen::RowVector3f((float)ray_o(0, 0), (float)ray_o(1, 0), (float)ray_o(2, 0));
        } else {
            o_i = Eigen::RowVector3f((float)ray_o(i, 0), (float)ray_o(i, 1), (float)ray_o(i, 2));
        }
        Eigen::RowVector3f d_i((float)ray_d(i, 0), (float)ray_d(i, 1), (float)ray_d(i, 2));
        Hit hit;

        bool is_hit = isector.intersectRay(o_i, d_i, hit, ray_near, ray_far);
        if (is_hit) {
            ret_fid(i, 0) = (npe_Scalar_f) hit.id;
            ret_bc(i, 0) = 1.0 - hit.u - hit.v;
            ret_bc(i, 1) = hit.u;
            ret_bc(i, 2) = hit.v;
            ret_t(i, 0) = hit.t;
        } else {
            ret_fid(i, 0) = -1;
            ret_t(i, 0) = std::numeric_limits<npe_Scalar_v>::infinity();
            ret_bc(i, 0) = 0.0;
            ret_bc(i, 1) = 0.0;
            ret_bc(i, 2) = 0.0;
        }

    }

    return std::make_tuple(npe::move(ret_fid), npe::move(ret_bc), npe::move(ret_t));

}
npe_end_code()
