#include <npe.h>
#include <igl/per_vertex_normals.h>
#include <string>

#include "common/common.h"
#include <igl/cotmatrix.h>
#include <igl/barycenter.h>



const char* laplacian_smooth_mesh_doc = R"igl_Qu8mg5v7(
Smooth a mesh using Laplacian smoothing

Args:
    v : \#v by 3 Matrix of mesh vertex 3D positions
    f : \#f by 3 Matrix of face (triangle) indices
    num_iters : Number of smoothing iterations to perform
    step_size : Step size per iteration (default is 1e-4)

Returns:
    n : list of vertex normals of shape #v by 3

)igl_Qu8mg5v7";
npe_function(laplacian_smooth_mesh)
npe_doc(laplacian_smooth_mesh_doc)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64)
npe_arg(num_iters, int)
npe_default_arg(step_size, double, 1e-1)
npe_default_arg(use_cotan_weights, bool, false)
npe_begin_code()
{
    validate_mesh(v, f);
    if (num_iters < 0) {
        throw pybind11::value_error("Invalid value for argument num_iters in smooth_mesh_laplacian. Must be a positive integer or 0.");
    }
    if (num_iters == 0) {  // No-op for 0
        return npe::move(v);
    }

    // Heuristic but kind of works maybe???
    if (use_cotan_weights) {
        step_size *= 1e-3;
    }
    Eigen::SparseMatrix<npe_Scalar_v> L;
    Eigen::SparseMatrix<npe_Scalar_v> M(v.rows(), v.rows());
    M.setIdentity();
    EigenDenseLike<npe_Matrix_v> u = v;
    EigenDenseLike<npe_Matrix_v> barycenter;
    Eigen::Matrix<npe_Scalar_v, Eigen::Dynamic, 1> dblA;
    igl::cotmatrix(v, f, L);
    
    for (int iter = 0; iter < num_iters; iter += 1) {
        if (use_cotan_weights) {
            igl::massmatrix(u, f, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
        }
        const auto & S = (M - step_size*L);

        Eigen::SimplicialLLT<Eigen::SparseMatrix<npe_Scalar_v>> solver(S);
        u = solver.solve(M*u).eval();

        if (use_cotan_weights) {
            igl::doublearea(u, f, dblA);
            double area = 0.5*dblA.sum();
            igl::barycenter(u, f, barycenter);
            
            Eigen::Matrix<npe_Scalar_v, 1, 3> centroid(0,0,0);
            for(int i = 0; i < barycenter.rows();i++) {
                centroid += 0.5*dblA(i)/area*barycenter.row(i);
            }
            u.rowwise() -= centroid;
            u.array() /= sqrt(area);
            u.rowwise() += centroid;
        }
    }

    return npe::move(u);
}
npe_end_code()
