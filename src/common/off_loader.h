#include <npe.h>

#include <string>
#include <unordered_map>

#include <igl/readOFF.h>
#include <igl/writeOFF.h>


std::unordered_map<std::string, pybind11::object> load_mesh_off(const std::string& filename) {
    
    std::vector<std::vector<double> > vV, vN, vC;
    std::vector<std::vector<int64_t> > fF;
    igl::readOFF(filename, vV, fF, vN, vC);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> v, vn, vc;
    Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> f;
    
    v.resize(vV.size(), 3);
    std::cerr << "LOADED THIS SHIZ WITH " << v.rows() << " vertices" << std::endl;
    if (vN.size() > 0) {
        vn.resize(vN.size(), 3);
    }
    if (vC.size() > 0) {
        vc.resize(vC.size(), 4);
    }
    for (int i = 0; i < vV.size(); i += 1) {
        if (vV[i].size() != 3) {
            throw std::runtime_error("Invalid vertices have " + std::to_string(vV[i].size()) + " dimensions. Expected 3.");
        }
        v(i, 0) = vV[i][0];
        v(i, 1) = vV[i][1];
        v(i, 2) = vV[i][2];

        if (vN.size() > 0) {
            if (vV[i].size() != 3) {
                throw std::runtime_error("Invalid vertex normals have " + std::to_string(vN[i].size()) + " dimensions. Expected 3.");
            }
            vn(i, 0) = vN[i][0];
            vn(i, 1) = vN[i][1];
            vn(i, 2) = vN[i][2];
        }

        if (vC.size() > 0) {
            if (vV[i].size() != 3) {
                throw std::runtime_error("Invalid vertex colors have " + std::to_string(vC[i].size()) + " dimensions. Expected 3.");
            }
            vc(i, 0) = vC[i][0];
            vc(i, 1) = vC[i][1];
            vc(i, 2) = vC[i][2];
            vc(i, 3) = 1.;
        }
    }

    if (fF.size() > 0) {
        f.resize(fF.size(), 3);
    }
    for (int i = 0; i < fF.size(); i += 1) {
        if (fF[i].size() != 3) {
            throw std::runtime_error("Invalid faces have " + std::to_string(fF[i].size()) + " vertices. Expected 3.");
        }
        f(i, 0) = fF[i][0];
        f(i, 1) = fF[i][1];
        f(i, 2) = fF[i][2];
    }

    std::unordered_map<std::string, pybind11::object> vertex_ret;
    vertex_ret.insert(std::make_pair("positions", npe::move(v)));
    if (vn.rows() > 0) {
        vertex_ret.insert(std::make_pair("normals", npe::move(vn)));
    }
    if (vc.rows() > 0) {
        vertex_ret.insert(std::make_pair("colors", npe::move(vc)));
    }
    
    std::unordered_map<std::string, pybind11::object> face_ret;
    face_ret.insert(std::make_pair("vertex_ids", npe::move(f)));
    
    std::unordered_map<std::string, pybind11::object> ret;
    ret["vertex_data"] = pybind11::cast(vertex_ret);
    ret["face_data"] = pybind11::cast(face_ret);   
    ret["textures"] = pybind11::cast(std::vector<std::string>());
    ret["normal_maps"] = pybind11::cast(std::vector<std::string>());

    return ret;
}




void save_mesh_off(std::string filename,
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

    if (has_f_normals) {
        std::string warning_str = "Face normals are ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_f_colors) {
        std::string warning_str = "Face colors are ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_f_quality) {
        std::string warning_str = "Face quality is ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_f_flags) {
        std::string warning_str = "Face flags are ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_w_texcoords) {
        std::string warning_str = "Wedge texcoords are ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_w_texids) {
        std::string warning_str = "Wedge texids are ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_texcoords) {
        std::string warning_str = "Vertex texcoords are ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_quality) {
        std::string warning_str = "Vertex quality is ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_radius) {
        std::string warning_str = "Vertex radius is ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_texids) {
        std::string warning_str = "Vertex texture ids is ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_flags) {
        std::string warning_str = "Vertex flags are ignored by off files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> v, vn, vc;
    Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> f;

    v.resize(num_vertices, 3);
    if (has_v_normals) {
        vn.resize(num_vertices, 3);
    }
    if (has_v_colors) {
        vc.resize(num_vertices, 3);
    }
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 3; j++) {
            if (dtype_f.equal(pybind11::dtype("f"))) {
                v(i, j) = *((float*) v_positions.data(i, j));
                if (has_v_normals) {
                    vn(i, j) = *((float*) v_normals.data(i, j));
                }
                if (has_v_colors) {
                    vc(i, j) = *((float*) v_colors.data(i, j));
                }
            } else if (dtype_f.equal(pybind11::dtype("d"))) {
                v(i, j) = *((double*) v_positions.data(i, j));
                if (has_v_normals) {
                    vn(i, j) = *((double*) v_normals.data(i, j));
                }
                if (has_v_colors) {
                    vc(i, j) = *((double*) v_colors.data(i, j));
                }
            } else {
                throw std::runtime_error("Bad float dtype -- this should never happen");
            }
        }
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
            }
        }
    }

    if (!igl::writeOFF(filename, v, f, vc)) {
        throw std::runtime_error("Failed to write stl file " + filename);
    }
}