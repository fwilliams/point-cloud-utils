#include <npe.h>

#include <string>
#include <fstream>
#include <unordered_map>

#include <igl/readSTL.h>
#include <igl/writeSTL.h>


std::unordered_map<std::string, pybind11::object> load_mesh_stl(const std::string& filename) {
    Eigen::MatrixXd v, n;
    Eigen::MatrixXi f;

    std::ifstream file;
    file.open(filename);
    if (file.fail()) {
        throw pybind11::value_error("Failed to open file " + filename);
    }

    igl::readSTL(file, v, f, n);
    
    std::unordered_map<std::string, pybind11::object> vertex_ret;
    vertex_ret.insert(std::make_pair("positions", npe::move(v)));

    std::unordered_map<std::string, pybind11::object> face_ret;
    face_ret.insert(std::make_pair("vertex_ids", npe::move(f)));
    if (n.rows() > 0) {
        face_ret.insert(std::make_pair("normals", npe::move(n)));
    }

    std::unordered_map<std::string, pybind11::object> ret;
    ret["vertex_data"] = pybind11::cast(vertex_ret);
    ret["face_data"] = pybind11::cast(face_ret);   
    ret["textures"] = pybind11::cast(std::vector<std::string>());
    ret["normal_maps"] = pybind11::cast(std::vector<std::string>());

    return ret;
}


void save_mesh_stl(std::string filename,
                   pybind11::array& v_positions, pybind11::array& v_normals,
                   pybind11::array& v_texcoords, pybind11::array& v_colors, pybind11::array& v_quality,
                   pybind11::array& v_radius, pybind11::array& v_texids, pybind11::array& v_flags,
                   pybind11::array& f_vertex_ids, pybind11::array& f_normals, pybind11::array& f_colors,
                   pybind11::array& f_quality, pybind11::array& f_flags,
                   pybind11::array& w_colors, pybind11::array& w_normals, pybind11::array& w_texcoords,
                   pybind11::array& w_texids,
                   pybind11::dict& custom_v_attribs,
                   pybind11::dict& custom_f_attribs,
                   std::vector<std::string>& textures, std::vector<std::string>& normal_maps,
                   pybind11::dtype dtype_f, pybind11::dtype dtype_i) {
    ssize_t num_vertices = v_positions.shape()[0];
    ssize_t num_faces = f_vertex_ids.shape()[0];

    bool has_v_positions = assert_shape_and_dtype(v_positions, "v_positions", dtype_f, {-num_vertices, 3});
    bool has_v_normals = assert_shape_and_dtype(v_normals, "v_normals", dtype_f, {-num_vertices, 3});
    bool has_v_texcoords = assert_shape_and_dtype(v_texcoords, "v_texcoords", dtype_f, {-num_vertices, 2});
    bool has_v_colors = assert_shape_and_dtype(v_colors, "v_colors", dtype_f, {-num_vertices, 4});
    bool has_v_quality = assert_shape_and_dtype(v_quality, "v_quality", dtype_f, {-num_vertices});
    bool has_v_radius = assert_shape_and_dtype(v_radius, "v_radius", dtype_f, {-num_vertices});
    bool has_v_texids = assert_shape_and_dtype(v_texids, "v_texids", dtype_i, {-num_vertices});
    bool has_v_flags = assert_shape_and_dtype(v_flags, "v_flags", dtype_i, {-num_vertices});

    bool has_f_vertex_ids = assert_shape_and_dtype(f_vertex_ids, "f_vertex_ids", dtype_i, {-num_faces, 3});
    bool has_f_normals = assert_shape_and_dtype(f_normals, "f_normals", dtype_f, {-num_faces, 3});
    bool has_f_colors = assert_shape_and_dtype(f_colors, "f_colors", dtype_f, {-num_faces, 4});
    bool has_f_quality = assert_shape_and_dtype(f_quality, "f_quality", dtype_f, {-num_faces});
    bool has_f_flags = assert_shape_and_dtype(f_flags, "f_flags", dtype_i, {-num_faces});

    bool has_w_texcoords = assert_shape_and_dtype(w_texcoords, "w_texcoords", dtype_f, {-num_faces, 3, 2});
    bool has_w_texids = assert_shape_and_dtype(w_texids, "w_texids", dtype_i, {-num_faces, 3});

    if (!dtype_f.equal(pybind11::dtype("f")) && !dtype_f.equal(pybind11::dtype("d"))) {
        throw pybind11::value_error("Vertex data must be float32 or float64");
    }

    if (!dtype_i.equal(pybind11::dtype("i")) && !dtype_i.equal(pybind11::dtype("l"))) {
        throw pybind11::value_error("Vertex data must be int32 or int64");
    }

    if (!has_v_positions) {
        throw pybind11::value_error("Must write at least vertex positions");
    }

    if (has_v_normals) {
        std::string warning_str = "Vertex normals are ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_f_colors) {
        std::string warning_str = "Face colors are ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_f_quality) {
        std::string warning_str = "Face quality is ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_f_flags) {
        std::string warning_str = "Face flags are ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_w_texcoords) {
        std::string warning_str = "Wedge texcoords are ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_w_texids) {
        std::string warning_str = "Wedge texids are ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_colors) {
        std::string warning_str = "Vertex colors are ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_texcoords) {
        std::string warning_str = "Vertex texcoords are ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_quality) {
        std::string warning_str = "Vertex quality is ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_radius) {
        std::string warning_str = "Vertex radius is ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_texids) {
        std::string warning_str = "Vertex texture ids is ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_flags) {
        std::string warning_str = "Vertex flags are ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> v, fn;
    Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> f;

    v.resize(num_vertices, 3);
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 3; j++) {
            if (dtype_f.equal(pybind11::dtype("f"))) {
                v(i, j) = *((float*) v_positions.data(i, j));
            } else if (dtype_f.equal(pybind11::dtype("d"))) {
                v(i, j) = *((double*) v_positions.data(i, j));
            } else {
                throw std::runtime_error("Bad float dtype -- this should never happen");
            }
        }
    }

    if (has_f_normals && !has_f_vertex_ids) {
        throw pybind11::value_error("Can't have face normals and no faces");
    }

    if (has_f_normals) {
        fn.resize(num_faces, 3);
    }
    if (has_f_vertex_ids) {
        f.resize(num_faces, 3);
        for (int i = 0; i < num_faces; i += 1) {
            for (int j = 0; j < 3; j += 1) {
                if (dtype_i.equal(pybind11::dtype("i"))) {
                    f(i, j) = *((int32_t*) f_vertex_ids.data(i, j));
                } else if (dtype_i.equal(pybind11::dtype("l"))) {
                    f(i, j) = *((int64_t*) f_vertex_ids.data(i, j));
                } else {
                    throw std::runtime_error("Bad integer dtype -- this should never happen");
                }

                if (has_f_normals) {
                    if (dtype_f.equal(pybind11::dtype("f"))) {
                        fn(i, j) = *((float*) f_normals.data(i, j));
                    } else if (dtype_f.equal(pybind11::dtype("d"))) {
                        fn(i, j) = *((double*) f_normals.data(i, j));
                    } else {
                        throw std::runtime_error("Bad float dtype -- this should never happen");
                    }
                }
            }
        }
    }

    if (has_v_normals) {
        if (!igl::writeSTL(filename, v, f, fn)) {
            throw std::runtime_error("Failed to write stl file " + filename);
        }
    } else {
        if (!igl::writeSTL(filename, v, f)) {
            throw std::runtime_error("Failed to write stl file " + filename);
        }
    }
}