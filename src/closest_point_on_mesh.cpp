#include <npe.h>

#include <igl/point_mesh_squared_distance.h>
#include <igl/barycentric_coordinates.h>

#include "common.h"


const char* closest_points_on_mesh_doc = R"Qu8mg5v7(
Compute distances from a set of points p to a triangle mesh (v, f)

Parameters
----------
p : (#p, 3)-shaped array of query point positions
v : (#v, 3)-shaped array of mesh vertex positions
f : (#f, 3)-shaped array of triangle face indices

Returns
-------
A tuple (d, f_idx, bc) where:
 - d is a (#p,)-shaped array of shortest distances for each query point p
 - f_idx is a (#p,)-shaped array of indices into f of the face containing the closest
   point to eaach query point
 - bc is a (#p, 3)-shaped array of barycentric coordinates for each query point

Notes
-----
Known bugs: This only computes distances to given primitives. So
unreferenced vertices are ignored. However, degenerate primitives are
handled correctly: triangle [1 2 2] is treated as a segment [1 2], and
triangle [1 1 1] is treated as a point. So one _could_ add extra
combinatorially degenerate rows to Ele for all unreferenced vertices to
also get distances to points.
)Qu8mg5v7";
npe_function(closest_points_on_mesh)
npe_arg(p, dense_float, dense_double)
npe_arg(v, npe_matches(p))
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_doc(closest_points_on_mesh_doc)
npe_begin_code()
{
    validate_point_cloud(p);
    validate_mesh(v, f);

    npe_Matrix_p sqr_dists, closest_points, barycentric_coords;
    npe_Matrix_f face_idxs;
    igl::point_mesh_squared_distance(p, v, f, sqr_dists, face_idxs, closest_points);

    npe_Matrix_p v_a(p.rows(), 3), v_b(p.rows(), 3), v_c(p.rows(), 3);
    for (int i = 0; i < face_idxs.rows(); i++) {
        v_a.row(i) = v.row(f(face_idxs(i), 0));
        v_b.row(i) = v.row(f(face_idxs(i), 1));
        v_c.row(i) = v.row(f(face_idxs(i), 2));
        sqr_dists(i) = sqrt(sqr_dists(i));
    }
    igl::barycentric_coordinates(closest_points, v_a, v_b, v_c, barycentric_coords);

    return std::make_tuple(npe::move(sqr_dists), npe::move(face_idxs), npe::move(barycentric_coords));
}
npe_end_code()
