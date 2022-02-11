#include <npe.h>
#include <igl/signed_distance.h>
#include <numeric>

#include "common/common.h"

const char* signed_distance_doc = R"igl_Qu8mg5v7(
Computes signed distances of a point cloud with respect to a Mesh using Fast Winding Numbers

Parameters
----------
p : (#p, 3)-shaped array point cloud (one 3D point per row)
v: (#v, 3)-shaped array of mesh vertex positions (one vertex position per row)
f: (#f, 3)-shaped array of mesh face indexes into v (a row (fi, fj, fk) indicate the 3 vertices of a face)
lower_bound: The minimum distance value possible (use this to clamp SDF values). negative infinite by default
upper_bound: The maximum distance value possible (use this to clamp SDF values). negative infinite by default

Returns
-------
A tuple (s, fi, bc) where:
  - s is a (#p,) shaped array of signed distance values for each query point in p
  - fi is a (#p,) shaped array of indices to the closest face for each query point in p
  - bc is a (#p, 3) shaped array of barycentric coordinates for the closest point on
    the mesh to each query point in p
)igl_Qu8mg5v7";
npe_function(signed_distance_to_mesh)
npe_arg(p, dense_float, dense_double)
npe_arg(v, npe_matches(p))
npe_arg(f, dense_int, dense_long, dense_longlong)
npe_default_arg(lower_bound, float, -std::numeric_limits<float>::infinity())
npe_default_arg(upper_bound, float, std::numeric_limits<float>::infinity())
npe_doc(signed_distance_doc)
npe_begin_code()
{
    validate_mesh(v, f);
    validate_point_cloud(p, false /* allow_0 */);

    EigenDense<npe_Scalar_p> s;
    EigenDense<int> fi;
    EigenDense<npe_Scalar_p> c;
    EigenDense<npe_Scalar_p> n;

    EigenDense<npe_Scalar_p> p_cp = p.template cast<npe_Scalar_p>();
    EigenDense<npe_Scalar_p> v_cp = v.template cast<npe_Scalar_p>();
    EigenDense<int> f_cp = f.template cast<int>(); // FIXME: LibIGL only works with int
    igl::signed_distance(p_cp, v_cp, f_cp, igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
                         (npe_Scalar_p) lower_bound, (npe_Scalar_p) upper_bound, s, fi, c, n);

    npe_Matrix_p v_a(p.rows(), 3), v_b(p.rows(), 3), v_c(p.rows(), 3);
    npe_Matrix_p barycentric_coords;
    for (int i = 0; i < fi.rows(); i++) {
        v_a.row(i) = v.row(f(fi(i), 0));
        v_b.row(i) = v.row(f(fi(i), 1));
        v_c.row(i) = v.row(f(fi(i), 2));
    }
    igl::barycentric_coordinates(c, v_a, v_b, v_c, barycentric_coords);

    return std::make_tuple(npe::move(s), npe::move(fi), npe::move(barycentric_coords));
}
npe_end_code()
