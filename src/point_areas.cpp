#include <npe.h>
#include <igl/copyleft/cgal/point_areas.h>
#include "common/common.h"

const char* estimate_point_cloud_areas_doc = R"igl_Qu8mg5v7(
Compute a consistent inside/outside field given a triangle soup and evaluate that field
at a set of query points

Args:
    p : (\#p, 3)-shaped array of 3D points
    i : (\#p, k)-shaped array of indices where i[a, b] is the index of b^th nearest neighbor to point p[a]
    n : (\#p, 3)-shaped array of point normals
Returns:
    areas : (\#o,)-shaped array with estimated areas for each point
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