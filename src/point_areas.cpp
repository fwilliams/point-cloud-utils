#include <npe.h>
#include <igl/copyleft/cgal/point_areas.h>
#include "common/common.h"

const char* estimate_point_cloud_areas_doc = R"igl_Qu8mg5v7(
Compute a consistent inside/outside field given a triangle soup and evaluate that field
at a set of query points

Parameters
----------
v: (#v, 3)-shaped array of mesh vertex positions (one vertex position per row)
f: (#f, 3)-shaped array of mesh face indexes into v (a row (fi, fj, fk) indicate the 3 vertices of a face)
q: (#q, 3)-shaped array of query positions at which to evaluat the winding number field

Returns
-------
A (#q,)-shaped array with a sign value for each query point (positive for outside and 
negative for inside)
)igl_Qu8mg5v7";
npe_function(estimate_point_cloud_areas)
npe_arg(p, dense_float, dense_double)
npe_arg(i, dense_int, dense_long, dense_longlong)
npe_arg(n, npe_matches(v))
npe_doc(estimate_point_cloud_areas_doc)
npe_begin_code()
{
    validate_point_cloud_normals(p, n)
    if(i.rows() != p.rows()) {
        throw pybind11::value_error("Invalid shape for i (knn indices), must be (#p, k) where #p is " + 
                                    "the number of input points (" + std::to_string(p.rows()) + 
                                    "). Got (" + std::to_string(i.rows()) + ", " + 
                                    std::to_string(i.cols()) + ")");
    }
    EigenDenseLike<npe_Matrix_p> a;
    igl::copyleft::cgal::point_areas(p, i, n, a);

    return npe::move(a);
}
npe_end_code()