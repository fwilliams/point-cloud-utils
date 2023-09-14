#include <unordered_map>

#include <npe.h>
#include <npe_typedefs.h>
#include <pybind11/stl.h>

#include "common/common.h"
#include "common/strutil.h"
#include "common/ply_loader.h"
#include "common/off_loader.h"
#include "common/stl_loader.h"
#include "common/obj_loader.h"
#include "common/numpy_utils.h"


npe_function(load_mesh_internal)
npe_arg(filename, std::string)
npe_arg(dtype, npe::dtype)
npe_begin_code()
{
    if (strutil::ends_with(strutil::to_lower(strutil::trim_copy(filename)), "ply")) {
        return load_mesh_ply(filename);
    } else if (strutil::ends_with(strutil::to_lower(strutil::trim_copy(filename)), "off")) {
        return load_mesh_off(filename);
    } else if (strutil::ends_with(strutil::to_lower(strutil::trim_copy(filename)), "stl")) {
        return load_mesh_stl(filename);
    } else if (strutil::ends_with(strutil::to_lower(strutil::trim_copy(filename)), "obj")) {
        return load_mesh_obj(filename);
    } else {
        throw pybind11::value_error("File extension type not supported for file " + filename + " (must be .ply, .off, .stl, or .obj)");
    }
}
npe_end_code()



npe_function(save_mesh_internal)
npe_arg(filename, std::string)

npe_arg(v_positions, pybind11::array)
npe_arg(v_normals, pybind11::array)
npe_arg(v_texcoords, pybind11::array)
npe_arg(v_colors, pybind11::array)
npe_arg(v_quality, pybind11::array)
npe_arg(v_radius, pybind11::array)
npe_arg(v_texids, pybind11::array)
npe_arg(v_flags, pybind11::array)

npe_arg(f_vertex_ids, pybind11::array)
npe_arg(f_normals, pybind11::array)
npe_arg(f_colors, pybind11::array)
npe_arg(f_quality, pybind11::array)
npe_arg(f_flags, pybind11::array)

npe_arg(w_colors, pybind11::array)
npe_arg(w_normals, pybind11::array)
npe_arg(w_texcoords, pybind11::array)
npe_arg(w_texids, pybind11::array)

npe_arg(custom_v_attribs, pybind11::dict)
npe_arg(custom_f_attribs, pybind11::dict)

npe_arg(textures, std::vector<std::string>)
npe_arg(normal_maps, std::vector<std::string>)

npe_arg(dtype_f, npe::dtype)
npe_arg(dtype_i, npe::dtype)

npe_begin_code()
{
    if (strutil::ends_with(strutil::to_lower(strutil::trim_copy(filename)), "ply")) {
        save_mesh_ply(filename,
                      v_positions,
                      v_normals,
                      v_texcoords,
                      v_colors,
                      v_quality,
                      v_radius,
                      v_texids,
                      v_flags,
                      f_vertex_ids,
                      f_normals,
                      f_colors,
                      f_quality,
                      f_flags,
                      w_colors,
                      w_normals,
                      w_texcoords,
                      w_texids,
                      custom_v_attribs,
                      custom_f_attribs,
                      textures,
                      normal_maps,
                      dtype_f,
                      dtype_i);
        return;
    } else if (strutil::ends_with(strutil::to_lower(strutil::trim_copy(filename)), "obj")) {
        save_mesh_obj(filename,
                      v_positions,
                      v_normals,
                      v_texcoords,
                      v_colors,
                      v_quality,
                      v_radius,
                      v_texids,
                      v_flags,
                      f_vertex_ids,
                      f_normals,
                      f_colors,
                      f_quality,
                      f_flags,
                      w_colors,
                      w_normals,
                      w_texcoords,
                      w_texids,
                      custom_v_attribs,
                      custom_f_attribs,
                      textures,
                      normal_maps,
                      dtype_f,
                      dtype_i);
        return;
    } else if (strutil::ends_with(strutil::to_lower(strutil::trim_copy(filename)), "stl")) {
        save_mesh_stl(filename,
                      v_positions,
                      v_normals,
                      v_texcoords,
                      v_colors,
                      v_quality,
                      v_radius,
                      v_texids,
                      v_flags,
                      f_vertex_ids,
                      f_normals,
                      f_colors,
                      f_quality,
                      f_flags,
                      w_colors,
                      w_normals,
                      w_texcoords,
                      w_texids,
                      custom_v_attribs,
                      custom_f_attribs,
                      textures,
                      normal_maps,
                      dtype_f,
                      dtype_i);
        return;
    } else if (strutil::ends_with(strutil::to_lower(strutil::trim_copy(filename)), "off")) {
        save_mesh_off(filename,
                      v_positions,
                      v_normals,
                      v_texcoords,
                      v_colors,
                      v_quality,
                      v_radius,
                      v_texids,
                      v_flags,
                      f_vertex_ids,
                      f_normals,
                      f_colors,
                      f_quality,
                      f_flags,
                      w_colors,
                      w_normals,
                      w_texcoords,
                      w_texids,
                      custom_v_attribs,
                      custom_f_attribs,
                      textures,
                      normal_maps,
                      dtype_f,
                      dtype_i);
        return;
    } else {
        throw pybind11::value_error("File extension type not supported for file " + filename + " (must be .ply, .off, .stl, or .obj)");
    }
}
npe_end_code()