#include <npe.h>

#include "common/common.h"



const char* triangle_mesh_face_areas_doc = R"Qu8mg5v7(
Compute the areas of each face of a triangle mesh

Parameters
----------
v : #v by 3 array of vertex positions (each row is a vertex)
f : #f by 3 Matrix of face (triangle) indices

Returns
-------
A numpy array, areas, of shape (#faces,) where areas[i] is the area of the face f[i] 
)Qu8mg5v7";
npe_function(mesh_face_areas)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_long, dense_longlong)
npe_default_arg(num_threads, int, -1)
npe_doc(triangle_mesh_face_areas_doc)
npe_begin_code()
    using Real = npe_Scalar_v;
    using Vec3 = Eigen::Matrix<Real, 1, 3>;

    validate_mesh(v, f);
    auto set_parallel = OmpSetParallelism(num_threads);
    EigenDenseLike<npe_Matrix_v> areas(f.rows(), 1);

    bool threw_exception = false;

    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
        #if defined(_OPENMP)
        #pragma omp for nowait
        #endif
        for (int i = 0; i < f.rows(); i += 1) {
            if (PyErr_CheckSignals() != 0) {
                #if defined(_OPENMP)
                    if (threw_exception) {
                        continue;
                    }
                    #pragma omp critical
                    {
                        threw_exception = true;
                    }
                    // This doesn't work on Windows :(
                    //#pragma omp cancel for
                #else
                    threw_exception = true;
                    break;
                #endif
            }

            const Vec3 v1 = v.row(f(i, 0));
            const Vec3 v2 = v.row(f(i, 1));
            const Vec3 v3 = v.row(f(i, 2));

            const Real a = (v2 - v1).norm();
            const Real b = (v3 - v2).norm();
            const Real c = (v1 - v3).norm();

            const Real p = 0.5 * (a + b + c);

            areas(i, 0) = sqrt(p * std::max(p - a, (Real) 0.0) * 
                                std::max(p - b, (Real) 0.0) * 
                                std::max(p - c, (Real) 0.0));
        }
    }


    if (threw_exception) {
        throw pybind11::error_already_set();
    }

    return npe::move(areas);
npe_end_code()