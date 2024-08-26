#include <npe.h>

#include <igl/round.h>
#include <igl/unique_rows.h>
#include <igl/colon.h>

#include "common/common.h"

namespace {

template <
    typename DerivedV, 
    typename DerivedSV, 
    typename DerivedSVI, 
    typename DerivedSVJ>
inline void remove_duplicate_vertices(
    const Eigen::MatrixBase<DerivedV>& V,
    const double epsilon,
    Eigen::PlainObjectBase<DerivedSV>& SV,
    Eigen::PlainObjectBase<DerivedSVI>& SVI,
    Eigen::PlainObjectBase<DerivedSVJ>& SVJ) {

    static_assert(
        (DerivedSVI::RowsAtCompileTime == 1 || DerivedSVI::ColsAtCompileTime == 1) &&
        (DerivedSVJ::RowsAtCompileTime == 1 || DerivedSVJ::ColsAtCompileTime == 1),
        "SVI and SVJ need to have RowsAtCompileTime == 1 or ColsAtCompileTime == 1");
    if(epsilon > 0) {
        DerivedV rV,rSV;
        igl::round((V/(epsilon)).eval(),rV);
        igl::unique_rows(rV,rSV,SVI,SVJ);
        SV = V(SVI.derived(),Eigen::all);
    } else {
        igl::unique_rows(V,SV,SVI,SVJ);
    }
}

template <
    typename DerivedV, 
    typename DerivedF,
    typename DerivedSV, 
    typename DerivedSVI, 
    typename DerivedSVJ,
    typename DerivedSF>
inline void remove_duplicate_vertices(
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedF>& F,
    const double epsilon,
    Eigen::PlainObjectBase<DerivedSV>& SV,
    Eigen::PlainObjectBase<DerivedSVI>& SVI,
    Eigen::PlainObjectBase<DerivedSVJ>& SVJ,
    Eigen::PlainObjectBase<DerivedSF>& SF) {
    // SVI and SVJ need to have RowsAtCompileTime == 1 or ColsAtCompileTime == 1
    static_assert(
        (DerivedSVI::RowsAtCompileTime == 1 || DerivedSVI::ColsAtCompileTime == 1) &&
        (DerivedSVJ::RowsAtCompileTime == 1 || DerivedSVJ::ColsAtCompileTime == 1),
        "SVI and SVJ need to have RowsAtCompileTime == 1 or ColsAtCompileTime == 1");
    using namespace Eigen;
    using namespace std;
    remove_duplicate_vertices(V,epsilon,SV,SVI,SVJ);
    SF.resizeLike(F);
    int64_t fcount = 0;
    for(int f = 0; f < F.rows(); f++) {
        bool is_degen = false;
        for(int c = 0; c < F.cols(); c++) {
            for (int c2 = c + 1; c2 < F.cols(); c2++) {
                if (SVJ(F(f,c)) == SVJ(F(f,c2))) {
                    is_degen = true;
                    break;
                }
            }
            if (is_degen) {
                break;
            }
            SF(fcount,c) = SVJ(F(f,c));
        }
        if (!is_degen) {
            fcount++;
        }
    }
    SF.conservativeResize(fcount, F.cols());
}


}





const char* remove_duplicate_points_doc = R"igl_Qu8mg5v7(
Removes duplicated points from a point cloud where two points are considered the same if their distance is below
some threshold

Args:
    x : \#x by 3 Matrix of 3D positions
    epsilon: threshold below which two points are considered equal
    return_index: If true, return indices to map between input and output

Returns:
    x_new : \#x_new x 3 Point cloud with duplicates removed
    if return indices is set, this function also returns:
        svi : \#x x 1 indices so that x_new = x[svi]
        svj : \#x_new x 1 indices so that x = x_new[svj]

See also:
    deduplicate_mesh_vertices
)igl_Qu8mg5v7";
npe_function(deduplicate_point_cloud)
npe_doc(remove_duplicate_points_doc)
npe_arg(points, dense_float, dense_double)
npe_arg(epsilon, double)
npe_default_arg(return_index, bool, true)
npe_begin_code()
{
    validate_point_cloud(points);
    Eigen::Matrix<npe_Scalar_points, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> x_copy = points;
    Eigen::Matrix<npe_Scalar_points, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> x_out;
    Eigen::Matrix<int32_t, Eigen::Dynamic, 1> svj;
    Eigen::Matrix<int32_t, Eigen::Dynamic, 1> svi;

    remove_duplicate_vertices(x_copy, epsilon, x_out, svi, svj);

    if (return_index) {
        return pybind11::cast(std::make_tuple(npe::move(x_out), npe::move(svi), npe::move(svj)));
    } else {
        return npe::move(x_out);
    }
}
npe_end_code()



const char* remove_duplicate_mesh_vertices_doc = R"igl_Qu8mg5v7(
Removes duplicated vertices from a triangle mesh two vertices are considered the same if their distance is below
some threshold

Args:
    v : \#v by 3 Matrix of mesh vertex 3D positions
    f : \#f by 3 Matrix of face (triangle) indices
    epsilon: threshold below which two points are considered equal
    return_index: If true, return indices to map between input and output

Returns:
    v_out : \#v x 3 array of mesh vertices with duplicates removed
    f_out : \#f x 3 array of mesh faces corresponding to the deduplicated mesh
    svi : \#x x 1 indices so that v_out = v[svi] (only returned if return_index is True)
    svj : \#x_new x 1 indices so that v = v_out[svj] (only returned if return_index is True)

See also:
    deduplicate_point_cloud
)igl_Qu8mg5v7";
npe_function(deduplicate_mesh_vertices)
npe_doc(remove_duplicate_mesh_vertices_doc)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64)
npe_arg(epsilon, double)
npe_default_arg(return_index, bool, true)
npe_begin_code()
{
    validate_mesh(v, f);
    Eigen::Matrix<npe_Scalar_v, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> v_copy = v;
    Eigen::Matrix<npe_Scalar_f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> f_copy = f;
    Eigen::Matrix<npe_Scalar_v, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> v_out;
    Eigen::Matrix<npe_Scalar_f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> f_out;
    Eigen::Matrix<int32_t, Eigen::Dynamic, 1> svi;
    Eigen::Matrix<int32_t, Eigen::Dynamic, 1> svj;

    remove_duplicate_vertices(v_copy, f_copy, epsilon, v_out, svi, svj, f_out);

    if (return_index) {
        return pybind11::cast(std::make_tuple(npe::move(v_out), npe::move(f_out), npe::move(svi), npe::move(svj)));
    } else {
        return pybind11::cast(std::make_tuple(npe::move(v_out), npe::move(f_out)));
    }
}
npe_end_code()
