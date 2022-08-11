#ifndef PLY_LOADER_H
#define PLY_LOADER_H

#include <unordered_map>
#include <fstream>

#define TINYPLY_IMPLEMENTATION
#include <tinyply/source/tinyply.h>
#include <npe.h>
#include <npe_typedefs.h>
#include <pybind11/stl.h>

#include "common/strutil.h"
#include "common/numpy_utils.h"


pybind11::dtype ply_type_to_dtype(const std::shared_ptr<tinyply::PlyData> data) {
    switch(data->t) {
        case tinyply::Type::FLOAT32: {
            return pybind11::dtype::of<std::float_t>();
        }
        case tinyply::Type::FLOAT64: {
            return pybind11::dtype::of<std::double_t>();
        }
        case tinyply::Type::INT8: {
            return pybind11::dtype::of<std::int8_t>();
        }
        case tinyply::Type::INT16: {
            return pybind11::dtype::of<std::int16_t>();
        }
        case tinyply::Type::INT32: {
            return pybind11::dtype::of<std::int32_t>();
        }
        case tinyply::Type::UINT8: {
            return pybind11::dtype::of<std::uint8_t>();
        }
        case tinyply::Type::UINT16: {
            return pybind11::dtype::of<std::uint16_t>();
        }
        case tinyply::Type::UINT32: {
            return pybind11::dtype::of<std::uint32_t>();
        }
        default: {
            throw std::runtime_error("Internal PLY loading error. Cannot determine type.");
        }
    }
}

tinyply::Type dtype_to_ply_type(pybind11::dtype dtype) {
    if (dtype.equal(pybind11::dtype::of<std::float_t>())) {
        return tinyply::Type::FLOAT32;
    } else if (dtype.equal(pybind11::dtype::of<std::double_t>())) {
        return tinyply::Type::FLOAT64;
    } else if (dtype.equal(pybind11::dtype::of<std::int8_t>())) {
        return tinyply::Type::INT8;
    } else if (dtype.equal(pybind11::dtype::of<std::int16_t>())) {
        return tinyply::Type::INT16;
    } else if (dtype.equal(pybind11::dtype::of<std::int32_t>())) {
        return tinyply::Type::INT32;
    } else if (dtype.equal(pybind11::dtype::of<std::uint8_t>())) {
        return tinyply::Type::UINT8;
    } else if (dtype.equal(pybind11::dtype::of<std::uint16_t>())) {
        return tinyply::Type::UINT16;
    } else if (dtype.equal(pybind11::dtype::of<std::uint32_t>())) {
        return tinyply::Type::UINT32;
    } else {
        throw std::runtime_error("Internal PLY loading error. Cannot determine type.");
    }
}


pybind11::array ply_data_to_array(std::shared_ptr<tinyply::PlyData> attrib) {
    pybind11::dtype attrib_dtype = ply_type_to_dtype(attrib);
    if (attrib->count == 0) {
        return pybind11::array(attrib_dtype, std::vector<size_t>({0, 0}));
    }
    size_t bytes_per_scalar = attrib_dtype.elsize();
    if (bytes_per_scalar <= 0) {
        throw std::runtime_error("Internal PLY loading error. Type has no defined byte size.");
    }
    size_t numel = attrib->buffer.size_bytes() / bytes_per_scalar;
    size_t num_rows = attrib->count;
    size_t num_cols = numel / num_rows;
    if (attrib->buffer.size_bytes() != num_rows * num_cols * bytes_per_scalar) {
        throw std::runtime_error("PLY loading internal error");
    }
    if (num_cols == 1) {
        pybind11::array attrib_array(attrib_dtype, std::vector<size_t>({num_rows}));
        std::memcpy(attrib_array.mutable_data(), attrib->buffer.get(), attrib->buffer.size_bytes());
        return attrib_array;
    } else {
        pybind11::array attrib_array(attrib_dtype, std::vector<size_t>({num_rows, num_cols}));
        std::memcpy(attrib_array.mutable_data(), attrib->buffer.get(), attrib->buffer.size_bytes());
        return attrib_array;
    }

}

std::shared_ptr<tinyply::PlyData> request_properties_from_element(
        tinyply::PlyFile& ply_file,
        std::unordered_map<std::string, std::set<std::string>>& loaded_properties,
        std::string elementKey, const std::vector<std::string> property_keys,
        const uint32_t list_size_hint = 0, std::string error_message = "") {
    std::shared_ptr<tinyply::PlyData> ret;
    try {
        ret = ply_file.request_properties_from_element(elementKey, property_keys, list_size_hint);
        for (int i = 0; i < property_keys.size(); i++) {
            loaded_properties[elementKey].insert(property_keys[i]);
        }
    } catch (const std::exception & e) {
        ret = nullptr;
        if (error_message != "") {
            throw std::runtime_error(error_message);
        }
    }
    return ret;
}


void load_mesh_ply(const std::string& filename, std::unordered_map<std::string, pybind11::object>& ret) {
    tinyply::PlyFile plyf;
    std::ifstream filestream (filename, std::ifstream::binary);

    plyf.parse_header(filestream);

    // Track baked in properties which we've already loaded so we don't load them twice
    std::unordered_map<std::string, std::set<std::string>> loaded_properties;
    loaded_properties["vertex"] = std::set<std::string>();
    loaded_properties["face"] = std::set<std::string>();

    // Set up loaders for baked in properties
    auto vertex_pos_ptr = request_properties_from_element(
            plyf, loaded_properties,
            "vertex", { "x", "y", "z" }, 0,
            "PLY file has no vertex positions.");
    auto vertex_normals_ptr = request_properties_from_element(
            plyf, loaded_properties,"vertex", { "nx", "ny", "nz" });
    auto vertex_colors_ptr = request_properties_from_element(
            plyf, loaded_properties,"vertex", { "red", "green", "blue" });
    auto vertex_alpha_ptr = request_properties_from_element(
            plyf, loaded_properties,"vertex", { "alpha" });
    auto vertex_quality_ptr = request_properties_from_element(
            plyf, loaded_properties,"vertex", { "quality" });
    auto vertex_radius_ptr = request_properties_from_element(
            plyf, loaded_properties,"vertex", { "radius" });
    auto vertex_flags_ptr = request_properties_from_element(
            plyf, loaded_properties,"face", { "flags" });
    std::shared_ptr<tinyply::PlyData> vertex_texcoord_ptr;
    try {
        vertex_texcoord_ptr = request_properties_from_element(
                plyf, loaded_properties, "vertex", { "s", "t" }, 0, "no ST");
    } catch (const std::exception & e) {
        try {
            vertex_texcoord_ptr = request_properties_from_element(
                    plyf, loaded_properties, "vertex", { "u", "v" }, 0, "no UV");
        } catch (const std::exception& e) {
            vertex_texcoord_ptr = request_properties_from_element(
                        plyf, loaded_properties, "vertex", { "texcoord" }, 0);
        }
    }

    auto face_indices_ptr = request_properties_from_element(
            plyf, loaded_properties, "face", { "vertex_indices" });
    auto face_flags_ptr = request_properties_from_element(
            plyf, loaded_properties, "face", { "flags" });
    auto face_colors_ptr = request_properties_from_element(
            plyf, loaded_properties, "face", { "red", "green", "blue" });
    auto face_alpha_ptr = request_properties_from_element(
            plyf, loaded_properties, "face", { "alpha" });
    auto face_quality_ptr = request_properties_from_element(
            plyf, loaded_properties, "face", { "quality" });
    auto face_normals_ptr = request_properties_from_element(
            plyf, loaded_properties, "face", { "nx", "ny", "nz" });
    auto wedge_texindex_ptr = request_properties_from_element(
            plyf, loaded_properties, "face", { "texnumber" });
    std::shared_ptr<tinyply::PlyData> wedge_texcoords_ptr;
    try {
        wedge_texcoords_ptr = request_properties_from_element(
                plyf, loaded_properties, "face", { "s", "t" }, 0, "no ST");
    } catch (const std::exception & e) {
        try {
            wedge_texcoords_ptr = request_properties_from_element(
                    plyf, loaded_properties, "face", { "u", "v" }, 0, "no UV");
        } catch (const std::exception & e) {
            wedge_texcoords_ptr = request_properties_from_element(
                    plyf, loaded_properties, "face", { "texcoord" }, 0);
        }

    }

    // Load paths to texture and normal maps files
    std::vector<std::string> texture_files;
    std::vector<std::string> normal_maps;
    for (const auto & c : plyf.get_comments()) {
        std::string trimmed = strutil::trim_copy(c);
        if (strutil::starts_with(strutil::to_lower(trimmed), "texturefile")) {
            std::string tex_filepath = strutil::trim_copy(trimmed.substr(11));
            texture_files.push_back(tex_filepath);
        }
        // std::cout << "\t[ply_header] Comment: " << c << std::endl;
    }

    // Set up loaders for custom vertex and face attributes
    std::vector<std::pair<std::string, std::shared_ptr<tinyply::PlyData>>> custom_vertex_attribs;
    std::vector<std::pair<std::string, std::shared_ptr<tinyply::PlyData>>> custom_face_attribs;
    for (const auto & e : plyf.get_elements()) {
        // std::cout << "\t[ply_header] element: " << e.name << " (" << e.size << ")" << std::endl;
        if (e.name != "vertex" && e.name != "face") {
            std::cerr << "Warning: Element [" << e.name << "] is not supported and will be ignored" << std::endl;
            continue;
        }

        for (const auto & p : e.properties) {
            // std::cout << "\t[ply_header] \tproperty: " << p.name << " (type=" << tinyply::PropertyTable[p.propertyType].str << ")";
            if (loaded_properties[e.name].find(p.name) == loaded_properties[e.name].end()) {
                if (e.name == "vertex") {
                    std::shared_ptr<tinyply::PlyData> dataptr;
                    if (p.isList) {
                        dataptr = plyf.request_properties_from_element(e.name, { p.name }, p.listCount);
                    } else {
                        dataptr = plyf.request_properties_from_element(e.name, { p.name });
                    }
                    custom_vertex_attribs.push_back(std::make_pair(p.name, dataptr));
                } else {
                    std::shared_ptr<tinyply::PlyData> dataptr;
                    if (p.isList) {
                        dataptr = plyf.request_properties_from_element(e.name, { p.name }, p.listCount);
                    } else {
                        dataptr = plyf.request_properties_from_element(e.name, { p.name });
                    }
                    custom_face_attribs.push_back(std::make_pair(p.name, dataptr));
                }
            }
            // if (p.isList) { std::cout << " (list_type=" << tinyply::PropertyTable[p.listType].str << ")"; }
            // std::cout << std::endl;
        }
    }

    plyf.read(filestream);

    std::unordered_map<std::string, pybind11::object> vertex_attribs;
    for (int i = 0; i < custom_vertex_attribs.size(); i += 1) {
        std::string attrib_name = custom_vertex_attribs[i].first;
        std::shared_ptr<tinyply::PlyData> attrib = custom_vertex_attribs[i].second;
        vertex_attribs.insert(std::make_pair(attrib_name, (pybind11::object) ply_data_to_array(attrib)));
    }

    std::unordered_map<std::string, pybind11::object> face_attribs;
    for (int i = 0; i < custom_face_attribs.size(); i += 1) {
        std::string attrib_name = custom_face_attribs[i].first;
        std::shared_ptr<tinyply::PlyData> attrib = custom_face_attribs[i].second;
        face_attribs.insert(std::make_pair(attrib_name, (pybind11::object) ply_data_to_array(attrib)));
    }

    std::unordered_map<std::string, pybind11::object> vertex_ret;
    if (vertex_pos_ptr) {
        vertex_ret.insert(std::make_pair("positions", (pybind11::object) ply_data_to_array(vertex_pos_ptr)));
    }
    if (vertex_normals_ptr) {
        vertex_ret.insert(std::make_pair("normals", (pybind11::object) ply_data_to_array(vertex_normals_ptr)));
    }
    if (vertex_colors_ptr) {
        vertex_ret.insert(std::make_pair("colors", (pybind11::object) ply_data_to_array(vertex_colors_ptr)));
    }
    if (vertex_texcoord_ptr) {
        vertex_ret.insert(std::make_pair("texcoords", (pybind11::object) ply_data_to_array(vertex_texcoord_ptr)));
    }
    if (vertex_alpha_ptr) {
        vertex_ret.insert(std::make_pair("alpha", (pybind11::object) ply_data_to_array(vertex_alpha_ptr)));
    }
    if (vertex_flags_ptr) {
        vertex_ret.insert(std::make_pair("flags", (pybind11::object) ply_data_to_array(vertex_flags_ptr)));
    }
    if (vertex_quality_ptr) {
        vertex_ret.insert(std::make_pair("quality", (pybind11::object) ply_data_to_array(vertex_quality_ptr)));
    }
    if (vertex_radius_ptr) {
        vertex_ret.insert(std::make_pair("radius", (pybind11::object) ply_data_to_array(vertex_radius_ptr)));
    }

    std::unordered_map<std::string, pybind11::object> face_ret;
    if (face_indices_ptr) {
        face_ret.insert(std::make_pair("vertex_ids", (pybind11::object) ply_data_to_array(face_indices_ptr)));
    }
    if (face_flags_ptr) {
        face_ret.insert(std::make_pair("flags", (pybind11::object) ply_data_to_array(face_flags_ptr)));
    }
    if (face_colors_ptr) {
        face_ret.insert(std::make_pair("colors", (pybind11::object) ply_data_to_array(face_colors_ptr)));
    }
    if (face_alpha_ptr) {
        face_ret.insert(std::make_pair("alpha", (pybind11::object) ply_data_to_array(face_alpha_ptr)));
    }
    if (face_quality_ptr) {
        face_ret.insert(std::make_pair("quality", (pybind11::object) ply_data_to_array(face_quality_ptr)));
    }
    if (face_normals_ptr) {
        face_ret.insert(std::make_pair("normals", (pybind11::object) ply_data_to_array(face_normals_ptr)));
    }
    if (wedge_texindex_ptr) {
        face_ret.insert(std::make_pair("wedge_tex_ids", (pybind11::object) ply_data_to_array(wedge_texindex_ptr)));
    }
    if (wedge_texcoords_ptr) {
        pybind11::array wedge_texcoords_array = ply_data_to_array(wedge_texcoords_ptr);
        wedge_texcoords_array.resize(std::vector<size_t>({(size_t) wedge_texcoords_array.shape(0), 3, 2}));
        face_ret.insert(std::make_pair("wedge_texcoords", (pybind11::object) wedge_texcoords_array));
    }

    vertex_ret["custom_attributes"] = pybind11::cast(vertex_attribs);
    pybind11::dict ret_vertex_data = pybind11::cast(vertex_ret);

    face_ret["custom_attributes"] = pybind11::cast(face_attribs);
    pybind11::dict ret_face_data = pybind11::cast(face_ret);
    pybind11::list ret_normalmaps = pybind11::list();
    pybind11::list ret_textures = pybind11::cast(texture_files);

    ret["vertex_data"] = ret_vertex_data;
    ret["face_data"] = ret_face_data;
    ret["textures"] = ret_textures;
    ret["normal_maps"] = ret_normalmaps;
}


void save_mesh_ply(std::string filename,
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

//    bool has_w_normals = assert_shape_and_dtype(w_normals, "w_normals", dtype_f, {-num_faces, 3, 3});
//    bool has_w_colors = assert_shape_and_dtype(w_colors, "w_colors", dtype_f, {-num_faces, 3, 4});
    bool has_w_texcoords = assert_shape_and_dtype(w_texcoords, "w_texcoords", dtype_f, {-num_faces, 3, 2});
    bool has_w_texids = assert_shape_and_dtype(w_texids, "w_texids", dtype_i, {-num_faces, 3});

    tinyply::PlyFile plyf;
    tinyply::Type ply_type_f;
    tinyply::Type ply_type_i = tinyply::Type::INT32;

    if (dtype_f.equal(pybind11::dtype::of<std::float_t>())) {
        ply_type_f = tinyply::Type::FLOAT32;
    } else if (dtype_f.equal(pybind11::dtype::of<std::double_t>())) {
        ply_type_f = tinyply::Type::FLOAT64;
    }

    if (has_w_texcoords) {
        plyf.add_properties_to_element(
                "face", { "texcoords" }, ply_type_f, num_faces,
                reinterpret_cast<std::uint8_t*>(w_texcoords.mutable_data()), tinyply::Type::UINT8, 6);
    }
    if (has_w_texids) {
        plyf.add_properties_to_element(
                "face", { "texnumber" }, ply_type_i, num_faces,
                reinterpret_cast<std::uint8_t*>(w_texcoords.mutable_data()), tinyply::Type::INVALID, 0);
    }


    if (has_f_vertex_ids) {
        plyf.add_properties_to_element(
                "face", { "vertex_indices" }, ply_type_i, num_faces,
                reinterpret_cast<std::uint8_t*>(f_vertex_ids.mutable_data()), tinyply::Type::UINT8, 3);
    }
    if (has_f_normals) {
        plyf.add_properties_to_element(
                "face", { "nx", "ny", "z" }, ply_type_f, num_faces,
                reinterpret_cast<std::uint8_t*>(f_normals.mutable_data()), tinyply::Type::INVALID, 0);
    }
    if (has_f_colors) {
        plyf.add_properties_to_element(
                "face", { "red", "green", "blue", "alpha" }, ply_type_f, num_faces,
                reinterpret_cast<std::uint8_t*>(f_colors.mutable_data()), tinyply::Type::INVALID, 0);
    }
    if (has_f_quality) {
        plyf.add_properties_to_element(
                "face", { "quality" }, ply_type_f, num_faces,
                reinterpret_cast<std::uint8_t*>(f_quality.mutable_data()), tinyply::Type::INVALID, 0);
    }
    if (has_f_flags) {
        plyf.add_properties_to_element(
                "face", { "flags" }, ply_type_i, num_faces,
                reinterpret_cast<std::uint8_t*>(f_flags.mutable_data()), tinyply::Type::INVALID, 0);
    }
    for (auto kv = custom_f_attribs.begin(); kv != custom_f_attribs.end(); kv++) {
        pybind11::str key = (pybind11::str) kv->first;
        pybind11::array value = pybind11::array::ensure(kv->second);
        if (value.shape(0) != num_faces) {
            throw pybind11::value_error("Invalid face attribute " + std::string(key) + ". Must have same number of rows as faces.");
        }
        size_t num_cols = 1;
        for (int i = 1; i < value.ndim(); i += 1) {
            num_cols *= value.shape(i);
        }
        if (num_cols == 0) {
            throw pybind11::value_error("Invalid face attribute " + std::string(key) + " has zero elements.");
        }
        try {
            tinyply::Type ply_dtype = dtype_to_ply_type(value.dtype());
            plyf.add_properties_to_element(
                    "face", { key }, ply_dtype, num_vertices,
                    reinterpret_cast<std::uint8_t*>(value.mutable_data()), tinyply::Type::UINT8, num_cols);
        } catch (const std::runtime_error& e) {
            throw pybind11::value_error("Invalid dtype for custom face attribute "+ std::string(key) + ".");
        }
    }

    if (has_v_positions) {
        plyf.add_properties_to_element(
                "vertex", { "x", "y", "z" }, ply_type_f, num_vertices,
                reinterpret_cast<std::uint8_t*>(v_positions.mutable_data()), tinyply::Type::INVALID, 0);
    }
    if (has_v_normals) {
        plyf.add_properties_to_element(
                "vertex", { "nx", "ny", "nz" }, ply_type_f, num_vertices,
                reinterpret_cast<std::uint8_t*>(v_normals.mutable_data()), tinyply::Type::INVALID, 0);
    }
    if (has_v_texcoords) {
        plyf.add_properties_to_element(
                "vertex", { "s", "t" }, ply_type_f, num_vertices,
                reinterpret_cast<std::uint8_t*>(v_texcoords.mutable_data()), tinyply::Type::INVALID, 0);
    }
    if (has_v_colors) {
        plyf.add_properties_to_element(
                "vertex", { "red", "green", "blue", "alpha" }, ply_type_f, num_vertices,
                reinterpret_cast<std::uint8_t*>(v_colors.mutable_data()), tinyply::Type::INVALID, 0);
    }
    if (has_v_quality) {
        plyf.add_properties_to_element(
                "vertex", { "quality" }, ply_type_f, num_vertices,
                reinterpret_cast<std::uint8_t*>(v_quality.mutable_data()), tinyply::Type::INVALID, 0);
    }
    if (has_v_radius) {
        plyf.add_properties_to_element(
                "vertex", { "radius" }, ply_type_f, num_vertices,
                reinterpret_cast<std::uint8_t*>(v_radius.mutable_data()), tinyply::Type::INVALID, 0);
    }
    if (has_v_texids) {
        plyf.add_properties_to_element(
                "vertex", { "texnumber" }, ply_type_i, num_vertices,
                reinterpret_cast<std::uint8_t*>(v_texids.mutable_data()), tinyply::Type::INVALID, 0);
    }
    if (has_v_flags) {
        plyf.add_properties_to_element(
                "vertex", { "flags" }, ply_type_i, num_vertices,
                reinterpret_cast<std::uint8_t*>(v_flags.mutable_data()), tinyply::Type::INVALID, 0);
    }
    for (auto kv = custom_v_attribs.begin(); kv != custom_v_attribs.end(); kv++) {
        pybind11::str key = (pybind11::str) kv->first;
        pybind11::array value = pybind11::array::ensure(kv->second);
        if (value.shape(0) != num_vertices) {
            throw pybind11::value_error("Invalid vertex attribute " + std::string(key) + 
                                        ". Must have same number of rows as vertices.");
        }
        size_t num_cols = 1;
        for (int i = 1; i < value.ndim(); i += 1) {
            num_cols *= value.shape(i);
        }
        if (num_cols == 0) {
            throw pybind11::value_error("Invalid vertex attribute " + std::string(key) + " has zero elements.");
        }
        try {
            if (num_cols == 1) {
                tinyply::Type ply_dtype = dtype_to_ply_type(value.dtype());
                plyf.add_properties_to_element(
                    "vertex", { key }, ply_dtype, num_vertices,
                    reinterpret_cast<std::uint8_t*>(value.mutable_data()), tinyply::Type::INVALID, 0);
            } else {
                tinyply::Type ply_dtype = dtype_to_ply_type(value.dtype());
                plyf.add_properties_to_element(
                    "vertex", { key }, ply_dtype, num_vertices,
                    reinterpret_cast<std::uint8_t*>(value.mutable_data()), tinyply::Type::UINT8, num_cols);
            }
        }
        catch (const std::runtime_error& e) {
            throw pybind11::value_error("Invalid dtype for custom vertex attribute "+ std::string(key) + ".");
        }
    }


    std::filebuf fb_binary;
    fb_binary.open(filename, std::ios::out | std::ios::binary);
    std::ostream outstream_binary(&fb_binary);
    if (outstream_binary.fail()) throw std::runtime_error("failed to open " + filename);
    plyf.write(outstream_binary, true);
}

#endif // PLY_LOADER_H