#include <npe.h>

#include <fstream>
#include <iostream>
#include <functional>
#include <cstdlib>


#include <igl/random_points_on_mesh.h>
#include <igl/blue_noise.h>
#include <igl/doublearea.h>

#include "common/common.h"


const char* sample_mesh_poisson_disk_doc = R"Qu8mg5v7(
Downsample a point set (possibly on a mesh) so that samples are approximately evenly spaced. This function uses the method in ("Parallel Poisson Disk Sampling with Spectrum Analysis on Surface")[http://graphics.cs.umass.edu/pubs/sa_2010.pdf]

Args:
    v : \#v by 3 array of mesh vertex positions
    f : \#f by 3 array of mesh face indices
    num_samples: desired number of Poisson Disk samples. Note that the actual number of returned samples will not be exactly this value (see sample_num_tolerance) to control the range of possible returned samples. Note: If this value <= 0, then the parameter radius is used to decide the number of samples
    radius : desired separation between points, if num_samples <= 0, then this value is used to determine the sampling (-1.0, by default).
    use_geodesic_distance : Use geodesic distance on the mesh downsampling. (True by default).
    best_choice_sampling : When downsampling, always keep the sample that will remove the fewest number of samples. (True by default).
    random_seed : A random seed used to generate the samples. Passing in 0 will use the current time. (0 by default).
    sample_num_tolerance: If you requested a target number of samples, by passsing num_samples > 0, then this function will return between (1 - sample_num_tolerance) * num_samples and (1 + sample_num_tolerance) * num_samples. Setting a very small value for this parameter will increase convergence time. (0.04 by default).
    oversampling_factor: To generate Poisson disk samples, we first generate a very dense (uniform) random sampling of the mesh, then prune these down to have the Poisson disk property. This parameter controls how many dense samples are generated. i.e. we generate oversampling_factor * num_samples samples (if you passed in radius, we estimate num_samples from the input points and radius). This parameter must be >= 1.0. (Default 40.0).
Returns:
    f_idx : a (m,)-shaped array of face indices into f where m is the number of Poisson-disk samples
    bc : a (m, 3)-shaped array of barycentric coordinates where m is the number of Poisson-disk samples

)Qu8mg5v7";
npe_function(sample_mesh_poisson_disk)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64)
npe_arg(num_samples, int)
npe_default_arg(radius, double, 0.0)
npe_default_arg(use_geodesic_distance, bool, true)
npe_default_arg(best_choice_sampling, bool, true)
npe_default_arg(random_seed, unsigned int, 0)
npe_default_arg(sample_num_tolerance, float, 0.04)
npe_default_arg(oversampling_factor, float, 40.0)
npe_doc(sample_mesh_poisson_disk_doc)
npe_begin_code()
{
    validate_mesh(v, f);

    if (num_samples <= 0 && radius <= 0.0) {
        throw pybind11::value_error("Cannot have both num_samples <= 0 and radius <= 0");
    }
    if (sample_num_tolerance > 1.0 || sample_num_tolerance <= 0.0) {
        throw pybind11::value_error("sample_num_tolerance must be in (0, 1]");
    }

    if (random_seed != 0) {
        srand(random_seed);
    }

    if(radius <= 0 && num_samples > 0) {
        const double total_area = [&](){Eigen::VectorXd A; igl::doublearea(v,f,A);return A.array().sum()/2;}();
        if (total_area <= 0) {
            throw pybind11::value_error("Mesh has zero area");
        }
        radius = sqrt(total_area / (0.7 * M_PI * num_samples)); // 0.7 is a density factor
    }

    EigenDenseLike<npe_Matrix_v> ret_bc, ret_p;
    Eigen::Matrix<npe_Scalar_f, Eigen::Dynamic, 1> ret_fi;
    igl::blue_noise(v, f, (npe_Scalar_v) radius, ret_bc, ret_fi, ret_p);


    return std::make_tuple(npe::move(ret_fi), npe::move(ret_bc));
}
npe_end_code()


const char* sample_mesh_random_doc = R"Qu8mg5v7(
Generate uniformly distributed random point samples on a mesh

Args:
    v : (\#v, 3)-shaped array of mesh vertex positions
    f : (\#f, 3)-shaped array of mesh face indices
    num_samples : The number of samples to generate
    random_seed : A random seed used to generate the samples. Passing in 0 will use the current time. (0 by default).
Returns:
    f_idx : (num_samples,) shaped array of face indices into f where
    bc : (num_samples, 3) shaped array of barycentric coordinates

)Qu8mg5v7";
npe_function(sample_mesh_random)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int32, dense_int64, dense_uint32, dense_uint64)
npe_arg(num_samples, int)
npe_default_arg(random_seed, unsigned int, 0)
npe_doc(sample_mesh_random_doc)
npe_begin_code()
{
    validate_mesh(v, f);

    if (num_samples <= 0) {
        throw pybind11::value_error("num_samples must be positive");
    }

    if (random_seed != 0) {
        srand(random_seed);
    }

    EigenDenseLike<npe_Matrix_v> ret_bc, ret_p;
    Eigen::Matrix<npe_Scalar_f, Eigen::Dynamic, 1> ret_fi;

    igl::random_points_on_mesh(num_samples, v, f, ret_bc, ret_fi, ret_p);
    return std::make_tuple(npe::move(ret_fi), npe::move(ret_bc));
}
npe_end_code()


