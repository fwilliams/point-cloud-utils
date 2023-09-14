#include <npe.h>

#include <string>
#include <filesystem>
#include <unordered_map>

#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader.h>

#include "common/numpy_utils.h"


std::unordered_map<std::string, pybind11::object> load_mesh_obj(const std::string& filename) {

    std::filesystem::path abspath = std::filesystem::absolute(filename);
    std::string dir = abspath.parent_path().string();

    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = dir; // Path to material files

    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(filename, reader_config)) {
        if (!reader.Error().empty()) {
            throw std::runtime_error("TinyObjReader: " + reader.Error());
        }
    }
    if (!reader.Warning().empty()) {
        std::string warning_str = "TinyObjReader: " + reader.Warning();
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> v, n, tc, c, vw;
    Eigen::Matrix<int64_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> f, mat_ids;

    auto& attrib = reader.GetAttrib();
    auto& shapes = reader.GetShapes();
    auto& materials = reader.GetMaterials();

    // Copy vertex data into Eigen arrays
    {
        const size_t total_vertices = attrib.vertices.size() / 3;
        v.resize(total_vertices, 3);
        if (attrib.normals.size() > 0) {
            n.resize(total_vertices, 3);
        }
        if (attrib.texcoords.size() > 0) {
            tc.resize(total_vertices, 3);
        }
        if (attrib.colors.size() > 0) {
            c.resize(total_vertices, 4);
        }
        if (attrib.vertex_weights.size() > 0) {
            vw.resize(total_vertices, 4);
        }
        for (int i = 0; i < total_vertices; i += 1) {
            v(i, 0) = attrib.vertices[3*i+0];
            v(i, 1) = attrib.vertices[3*i+1];
            v(i, 2) = attrib.vertices[3*i+2];

            if (attrib.normals.size() > 0) {
                n(i, 0) = attrib.normals[3*i+0];
                n(i, 1) = attrib.normals[3*i+1];
                n(i, 2) = attrib.normals[3*i+2];
            }

            if (attrib.texcoords.size() > 0) {
                tc(i, 0) = attrib.texcoords[2*i+0];
                tc(i, 1) = attrib.texcoords[2*i+1];
                if (attrib.texcoord_ws.size() > 0) {
                    tc(i, 2) = attrib.texcoord_ws[i];
                } else {
                    tc(i, 2) = 0;
                }
            }

            if (attrib.colors.size() > 0) {
                c(i, 0) = attrib.colors[3*i+0];
                c(i, 1) = attrib.colors[3*i+1];
                c(i, 2) = attrib.colors[3*i+2];
                if (attrib.colors.size() > 3*total_vertices) {
                    c(i, 3) = attrib.colors[4*i+3];
                } else {
                    c(i, 3) = 1;
                }
            }

            if (attrib.vertex_weights.size() > 0) {
                vw(i, 0) = attrib.vertex_weights[i+0];
            }
        }
    }

    // Copy face data into eigen arrays
    {
        size_t total_faces = 0;
        size_t total_mat_ids = 0;
        for (size_t s = 0; s < shapes.size(); s++) {
            total_faces += shapes[s].mesh.indices.size() / 3;
            total_mat_ids += shapes[s].mesh.material_ids.size();
        }
        f.resize(total_faces, 3);
        mat_ids.resize(total_mat_ids, 1);

        size_t f_offset = 0;
        for (size_t sidx = 0; sidx < shapes.size(); sidx += 1) {  // Loop over shapes
            size_t max_vertex_index = 0;

            size_t index_offset = 0; // Loop over faces(polygon)
            for (size_t fidx = 0; fidx < shapes[sidx].mesh.num_face_vertices.size(); fidx += 1) {
                const size_t num_vertices_in_face = size_t(shapes[sidx].mesh.num_face_vertices[fidx]);
                if (num_vertices_in_face != 3) {
                    throw std::runtime_error("Only triangular faces are supported. Got face with " + std::to_string(num_vertices_in_face) + " vertices.");
                }
                
                for (size_t vidx = 0; vidx < num_vertices_in_face; vidx += 1) {
                    tinyobj::index_t idx = shapes[sidx].mesh.indices[index_offset + vidx];
                    f(f_offset, vidx) = size_t(idx.vertex_index);
                }

                if (shapes[sidx].mesh.material_ids.size() > 0) {
                    mat_ids(f_offset, 0) = shapes[sidx].mesh.material_ids[fidx];
                } else if (mat_ids.size() > 0) {
                    mat_ids(f_offset, 0) = -1;
                }

                f_offset += 1;
                index_offset += num_vertices_in_face;
            }
        }
    }

    std::unordered_map<std::string, pybind11::object> vertex_ret;
    vertex_ret.insert(std::make_pair("positions", npe::move(v)));
    if (n.rows() > 0) {
        vertex_ret.insert(std::make_pair("normals", npe::move(n)));
    }
    if (tc.rows() > 0) {
        vertex_ret.insert(std::make_pair("texcoords", npe::move(tc)));
    }
    if (c.rows() > 0) {
        vertex_ret.insert(std::make_pair("colors", npe::move(c)));
    }

    std::unordered_map<std::string, pybind11::object> vertex_attribs;
    if (vw.rows() > 0) {
        vertex_attribs.insert(std::make_pair("weights", npe::move(vw)));
    }
    vertex_ret.insert(std::make_pair("custom_attributes", pybind11::cast(vertex_attribs)));

    std::unordered_map<std::string, pybind11::object> face_ret;
    face_ret.insert(std::make_pair("vertex_ids", npe::move(f)));

    std::unordered_map<std::string, pybind11::object> face_attribs;
    if (mat_ids.rows() > 0) {
        face_attribs.insert(std::make_pair("material_ids", npe::move(mat_ids)));
    }
    face_ret.insert(std::make_pair("custom_attributes", pybind11::cast(face_attribs)));

    // TODO: Support materials
    std::unordered_map<std::string, pybind11::object> ret;
    ret["vertex_data"] = pybind11::cast(vertex_ret);
    ret["face_data"] = pybind11::cast(face_ret);
    ret["textures"] = pybind11::cast(std::vector<std::string>());
    ret["normal_maps"] = pybind11::cast(std::vector<std::string>());

    return ret;
}



void save_mesh_obj(std::string filename,
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
        std::string warning_str = "Face normals are ignored by obj files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_f_colors) {
        std::string warning_str = "Face colors are ignored by obj files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_f_quality) {
        std::string warning_str = "Face quality is ignored by obj files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_f_flags) {
        std::string warning_str = "Face flags are ignored by obj files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_w_texcoords) {
        std::string warning_str = "Wedge texcoords are ignored by obj files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_w_texids) {
        std::string warning_str = "Wedge texids are ignored by obj files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_quality) {
        std::string warning_str = "Vertex quality is ignored by obj files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_radius) {
        std::string warning_str = "Vertex radius is ignored by obj files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_texids) {
        std::string warning_str = "Vertex texture ids is ignored by stl files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }
    if (has_v_flags) {
        std::string warning_str = "Vertex flags are ignored by obj files";
        PyErr_Warn(PyExc_RuntimeWarning, warning_str.c_str());
    }

    std::ofstream outstream(filename);
    for (int i = 0; i < num_vertices; i += 1) {
        outstream << "v ";
        for (int j = 0; j < 3; j++) {
            if (dtype_f.equal(pybind11::dtype("f"))) {
                float value = *((float*)v_positions.data(i, j));
                outstream << value << " ";
            } else if (dtype_f.equal(pybind11::dtype("d"))) {
                double value = *((double*)v_positions.data(i, j));
                outstream << value << " ";
            }
        }
        if (has_v_colors) {
            if (has_f_colors) {
                throw pybind11::value_error("Cannot specify both vertex colors and face colors");
            }
            for (int j = 0; j < 3; j++) {
                if (dtype_f.equal(pybind11::dtype("f"))) {
                    float value = *((float*)v_colors.data(i, j));
                    outstream << value << " ";
                } else if (dtype_f.equal(pybind11::dtype("d"))) {
                    double value = *((double*)v_colors.data(i, j));
                    outstream << value << " ";
                }
            }
        }
        outstream << std::endl;

        if (has_v_normals) {
            outstream << "vn ";
            for (int j = 0; j < 3; j++) {
                if (dtype_f.equal(pybind11::dtype("f"))) {
                    float value = *((float*)v_normals.data(i, j));
                    outstream << value << " ";
                } else if (dtype_f.equal(pybind11::dtype("d"))) {
                    double value = *((double*)v_normals.data(i, j));
                    outstream << value << " ";
                }
            }
            outstream << std::endl;
        }

        if (has_v_texcoords) {
            outstream << "vt ";
            int num_texcoords = v_texcoords.shape()[1];
            for (int j = 0; j < num_texcoords; j++) {
                if (dtype_f.equal(pybind11::dtype("f"))) {
                    float value = *((float*)v_texcoords.data(i, j));
                    outstream << value << " ";
                } else if (dtype_f.equal(pybind11::dtype("d"))) {
                    double value = *((double*)v_texcoords.data(i, j));
                    outstream << value << " ";
                }
            }
            outstream << std::endl;
        }
    }

    if (has_f_vertex_ids) {
        for (int i = 0; i < num_faces; i += 1) {
            outstream << "f ";
            for (int j = 0; j < 3; j++) {
                if (dtype_i.equal(pybind11::dtype("i"))) {
                    int32_t value = *((int32_t*)f_vertex_ids.data(i, j)) + 1;
                    outstream << value << "/";
                    if (has_v_texcoords) {
                        outstream << value;
                    }
                    outstream << "/";
                    if (has_v_normals) {
                        outstream << value;
                    }
                } else if (dtype_i.equal(pybind11::dtype("l"))) {
                    int64_t value = *((int64_t*)f_vertex_ids.data(i, j)) + 1;
                    outstream << value << "/";
                    if (has_v_texcoords) {
                        outstream << value;
                    }
                    outstream << "/";
                    if (has_v_normals) {
                        outstream << value;
                    }
                }
                outstream << " ";
            }
            outstream << std::endl;
        }
    }
   
}