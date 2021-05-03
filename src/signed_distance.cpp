#include <npe.h>
#include <igl/signed_distance.h>
#include <numeric>

#include "common.h"

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
A tuple (s, i, c) where
    s is a (#p,) shaped array of signed distance values for each point in p
    i is a (#p,) shaped array of indices to the closest face for each point in p
    c is a (#p, 3) shaped array of the closest point on the mesh for each point in p
)igl_Qu8mg5v7";
npe_function(signed_distance)
npe_arg(p, dense_float, dense_double)
npe_arg(v, npe_matches(p))
npe_arg(f, dense_int, dense_long, dense_longlong)
npe_default_arg(lower_bound, float, -std::numeric_limits<float>::infinity())
npe_default_arg(upper_bound, float, std::numeric_limits<float>::infinity())
npe_begin_code()
{
    validate_mesh(v, f);
    validate_point_cloud(p, false /* allow_0 */);

    EigenDense<npe_Scalar_p> s;
    EigenDense<int> i;
    EigenDense<npe_Scalar_p> c;
    EigenDense<npe_Scalar_p> n;

    EigenDense<npe_Scalar_p> p_cp = p.template cast<npe_Scalar_p>();
    EigenDense<npe_Scalar_p> v_cp = v.template cast<npe_Scalar_p>();
    EigenDense<int> f_cp = f.template cast<int>(); // FIXME: LibIGL only works with int
    igl::signed_distance(p_cp, v_cp, f_cp, igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
                         (npe_Scalar_p) lower_bound, (npe_Scalar_p) upper_bound, s, i, c, n);

    return std::make_tuple(npe::move(s), npe::move(i), npe::move(c));
}
npe_end_code()