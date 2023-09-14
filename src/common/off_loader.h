#include <npe.h>

#include <string>
#include <unordered_map>

#include <igl/readOFF.h>


std::unordered_map<std::string, pybind11::object> load_mesh_off(const std::string& filename) {
    
    Eigen::MatrixXd v;
    Eigen::MatrixXi f;
    igl::readOFF(filename, v, f);

    std::unordered_map<std::string, pybind11::object> vertex_ret;
    vertex_ret.insert(std::make_pair("positions", npe::move(v)));
    
    std::unordered_map<std::string, pybind11::object> face_ret;
    face_ret.insert(std::make_pair("vertex_ids", npe::move(f)));
    
    std::unordered_map<std::string, pybind11::object> ret;
    ret["vertex_data"] = pybind11::cast(vertex_ret);
    ret["face_data"] = pybind11::cast(face_ret);   
    ret["textures"] = pybind11::cast(std::vector<std::string>());
    ret["normal_maps"] = pybind11::cast(std::vector<std::string>());

    return ret;
}