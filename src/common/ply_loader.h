#pragma once

#include <unordered_map>
#include <fstream>

#include <tinyply/source/tinyply.h>
#include <npe.h>
#include <npe_typedefs.h>
#include <pybind11/stl.h>

#include "common/strutil.h"


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

    };
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
    pybind11::array attrib_array(attrib_dtype, std::vector<size_t>({num_rows, num_cols}));
    std::memcpy(attrib_array.mutable_data(), attrib->buffer.get(), attrib->buffer.size_bytes());
    return attrib_array;
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

