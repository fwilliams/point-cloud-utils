#include <igl/readOBJ.h>
#include <igl/readPLY.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>

#include <npe.h>
#include <npe_typedefs.h>

const int IglDefaultOptions = Eigen::RowMajor;
typedef Eigen::Matrix<std::float_t, Eigen::Dynamic, Eigen::Dynamic, IglDefaultOptions, Eigen::Dynamic, Eigen::Dynamic> EigenDenseF32;
typedef Eigen::Matrix<std::double_t, Eigen::Dynamic, Eigen::Dynamic, IglDefaultOptions, Eigen::Dynamic, Eigen::Dynamic> EigenDenseF64;
typedef Eigen::Matrix<std::int32_t, Eigen::Dynamic, Eigen::Dynamic, IglDefaultOptions, Eigen::Dynamic, Eigen::Dynamic> EigenDenseI32;
typedef Eigen::Matrix<std::int64_t, Eigen::Dynamic, Eigen::Dynamic, IglDefaultOptions, Eigen::Dynamic, Eigen::Dynamic> EigenDenseI64;

const char* ds_read_obj = R"igl_Qu8mg5v7(
Read a mesh from an ascii obj file, filling in vertex positions, normals
and texture coordinates. Mesh may have faces of any number of degree.

Parameters
----------
filename : string, path to .obj file
dtype : data-type of the returned faces, texture coordinates and normals, optional. Default is `float64`.
        (returned faces always have type int32.)

Returns
-------
v : array of vertex positions #v by 3
tc : array of texture coordinats #tc by 2
n : array of corner normals #n by 3
f : #f array of face indices into vertex positions
ftc : #f array of face indices into vertex texture coordinates
fn : #f array of face indices into vertex normals

See also
--------

Notes
-----

Examples
--------
>>> v, f, n = read_obj("my_model.obj")

)igl_Qu8mg5v7";

npe_function(read_obj)
npe_doc(ds_read_obj)
npe_arg(filename, std::string)
npe_default_arg(dtype, npe::dtype, "float64")
npe_begin_code()

  if (dtype.type() == npe::type_f32) {
    EigenDenseF32 v, tc, n;
    EigenDenseI32 f, ftc, fn;
    bool ret = igl::readOBJ(filename, v, tc, n, f, ftc, fn);
    if (!ret) {
      throw std::invalid_argument("File '" + filename + "' not found.");
    }
    return std::make_tuple(npe::move(v), npe::move(f), npe::move(n));
  } else if (dtype.type() == npe::type_f64) {
    EigenDenseF64 v, tc, n;
    EigenDenseI32 f, ftc, fn;
    bool ret = igl::readOBJ(filename, v, tc, n, f, ftc, fn);
    if (!ret) {
      throw std::invalid_argument("File '" + filename + "' not found.");
    }
    return std::make_tuple(npe::move(v), npe::move(f), npe::move(n));
  } else {
    throw pybind11::type_error("Only float32 and float64 dtypes are supported.");
  }

npe_end_code()




const char* ds_write_obj = R"igl_Qu8mg5v7(
Write a mesh in an ascii obj file.

Parameters
----------
filename : path to outputfile
v : #v by 3 array of vertex positions
f : #f x 3 array of face indices into vertex positions
n : #v x 3 array of vertex normals

Returns
-------
ret : bool if output was successful

See also
--------
read_obj

Notes
-----
None

Examples
--------
# Mesh in (v, f, n)
>>> success = write_obj(v, f, n)
)igl_Qu8mg5v7";

npe_function(write_obj)
npe_doc(ds_write_obj)
npe_arg(filename, std::string)
npe_arg(v, dense_f64, dense_f32)
npe_arg(f, dense_i32, dense_i64)
npe_arg(n, dense_f64, dense_f32)
npe_begin_code()

  npe_Matrix_n fn;
  npe_Matrix_n tc;
  npe_Matrix_n ftc;
  return igl::writeOBJ(filename, v, f, n, fn, tc, ftc);

npe_end_code()




const char* ds_read_off = R"igl_Qu8mg5v7(
Read a mesh from an ascii OFF file, filling in vertex positions, normals
and texture coordinates. Mesh may have faces of any number of degree.

Parameters
----------
filename : string, path to .off file
read_normals : bool, determines whether normals are read. If false, returns []
dtype : data-type of the returned vertices, faces, and normals, optional. Default is `float64`.
        (returned faces always have type int32.)

Returns
-------
v : array of vertex positions #v by 3
f : #f list of face indices into vertex positions
n : list of vertex normals #v by 3

See also
--------
read_triangle_mesh, read_obj

Notes
-----
None

Examples
--------
>>> v, f, n = read_off("my_model.off")
)igl_Qu8mg5v7";

npe_function(read_off)
npe_doc(ds_read_off)
npe_arg(filename, std::string)
npe_default_arg(read_normals, bool, true)
npe_default_arg(dtype, npe::dtype, "float64")
npe_begin_code()

  if (dtype.type() == npe::type_f32) {
    EigenDenseF32 v, n;
    EigenDenseI32 f;
    bool ret;
    if (read_normals) {
      ret = igl::readOFF(filename, v, f, n);
    }
    else {
      ret = igl::readOFF(filename, v, f);
    }

    if (!ret) {
      throw std::invalid_argument("File '" + filename + "' not found.");
    }

    return std::make_tuple(npe::move(v), npe::move(f), npe::move(n));
  } else if (dtype.type() == npe::type_f64) {
    EigenDenseF64 v, n;
    EigenDenseI32 f;
    bool ret;
    if (read_normals) {
      ret = igl::readOFF(filename, v, f, n);
    }
    else {
      ret = igl::readOFF(filename, v, f);
    }

    if (!ret) {
      throw std::invalid_argument("File '" + filename + "' not found.");
    }

    return std::make_tuple(npe::move(v), npe::move(f), npe::move(n));
  } else {
    throw pybind11::type_error("Only float32 and float64 dtypes are supported.");
  }

npe_end_code()




const char* ds_read_ply = R"igl_Qu8mg5v7(
Read a mesh from a PLY file, filling in vertex positions, normals
and texture coordinates. Mesh may have faces of any number of degree.

Parameters
----------
filename : string, path to .off file
dtype : data-type of the returned vertices, faces, and normals, optional. Default is `float64`.
        (returned faces always have type int32.)

Returns
-------
v : array of vertex positions #v by 3
f : #f list of face indices into vertex positions
n : #v by 3 or empt list of vertex normals
uv : #v by 2 (or empty ) list of uv coordinates

See also
--------
read_off, read_obj

Notes
-----
None

Examples
--------
>>> v, f, n, uv = read_ply("my_model.ply")
)igl_Qu8mg5v7";

npe_function(read_ply)
npe_doc(ds_read_ply)
npe_arg(filename, std::string)
npe_default_arg(dtype, npe::dtype, "float64")
npe_begin_code()

  if (dtype.type() == npe::type_f32) {
    EigenDenseF32 v, n, uv;
    EigenDenseI32 f;
    bool ret = igl::readPLY(filename, v, f, n, uv);

    if (!ret) {
      throw std::invalid_argument("Failed to read file '" + filename + "'");
    }

    return std::make_tuple(npe::move(v), npe::move(f), npe::move(n), npe::move(uv));
  } else if (dtype.type() == npe::type_f64) {
    EigenDenseF64 v, n, uv;
    EigenDenseI32 f;
    bool ret = igl::readPLY(filename, v, f, n, uv);

    if (!ret) {
      throw std::invalid_argument("Failed to read file '" + filename + "'");
    }

    return std::make_tuple(npe::move(v), npe::move(f), npe::move(n), npe::move(uv));
  } else {
    throw pybind11::type_error("Only float32 and float64 dtypes are supported.");
  }

npe_end_code()
