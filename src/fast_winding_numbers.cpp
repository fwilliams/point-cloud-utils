#include <npe.h>

#include <igl/fast_winding_number.h>

#include "common/common.h"

const char* triangle_soup_fast_winding_number_doc = R"igl_Qu8mg5v7(
Compute a consistent inside/outside field given a triangle soup and evaluate that field
at a set of query points

Args:
    v: (#v, 3)-shaped array of mesh vertex positions (one vertex position per row)
    f: (#f, 3)-shaped array of mesh face indexes into v (a row (fi, fj, fk) indicate the 3 vertices of a face)
    q: (#q, 3)-shaped array of query positions at which to evaluat the winding number field

Returns:
    A (#q,)-shaped array with a sign value for each query point (positive for outside and
    negative for inside)
)igl_Qu8mg5v7";
npe_function(triangle_soup_fast_winding_number)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64, dense_uint32, dense_uint64)
npe_arg(p, npe_matches(v))
npe_doc(triangle_soup_fast_winding_number_doc)
npe_begin_code()
{
    validate_mesh(v, f);
    validate_point_cloud(p, false /* allow_0 */);
    EigenDenseLike<npe_Matrix_v> w;
    igl::fast_winding_number(v, f, p, w);

    return npe::move(w);
}
npe_end_code()



const char* point_cloud_fast_winding_number_doc = R"igl_Qu8mg5v7(
Compute a consistent inside/outside field given a point cloud and evaluate that field
at a set of query points

Args:
    p : (#p, 3)-shaped array of points
    n : (#p, 3)-shaped array of point normals at p
    a : (#p, 1)-shaped array of point areas (e.g. computed with pcu.estimate_point_cloud_areas)
    q : (#q, 3)-shaped array of query positions at which to evaluat the winding number field

Returns:
    signs : (\#q,)-shaped array with a sign value for each query point (positive for outside and negative for inside)
)igl_Qu8mg5v7";
npe_function(point_cloud_fast_winding_number)
npe_arg(p, dense_float, dense_double)
npe_arg(n, npe_matches(p))
npe_arg(a, npe_matches(p))
npe_arg(q, npe_matches(p))
npe_doc(point_cloud_fast_winding_number_doc)
npe_begin_code()
{
    validate_point_cloud_normals(p, n);
    validate_point_cloud(q, false /* allow_0 */);
    EigenDenseLike<npe_Matrix_p> w;
    igl::fast_winding_number(p, n, a, q, w);

    return npe::move(w);
}
npe_end_code()

