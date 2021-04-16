#include <npe.h>
#include <octree/octree.cpp>

#include <iostream>


namespace py = pybind11;

void hack_extra_bindings(pybind11::module& m) {
    py::class_<Octree, std::shared_ptr<Octree>>(m, "Octree")
    .def(py::init([](size_t max_depth, double ox, double oy, double oz, double size){
        return std::shared_ptr<Octree>(new Octree(max_depth, Eigen::Vector3d(ox, oy, oz), size));
    }))
    .def("clear", &Octree::Clear)
    .def("is_empty", &Octree::IsEmpty)
    .def("get_min_bound", [](const Octree& octree) {
        Eigen::Vector3d ret = octree.GetMinBound();
        return std::make_tuple(ret[0], ret[1], ret[2]);
    })
    .def("get_max_bound", [](const Octree& octree) {
        Eigen::Vector3d ret = octree.GetMaxBound();
        return std::make_tuple(ret[0], ret[1], ret[2]);
    })
    .def("get_center", [](const Octree& octree) {
        Eigen::Vector3d ret = octree.GetCenter();
        return std::make_tuple(ret[0], ret[1], ret[2]);
    })
    .def("get_aabb", [](const Octree& octree) {
        std::tuple<Eigen::Vector3d, Eigen::Vector3d> ret = octree.GetAxisAlignedBoundingBox();
        return std::make_tuple(
                std::make_tuple(std::get<0>(ret)[0], std::get<0>(ret)[1], std::get<0>(ret)[2]),
                std::make_tuple(std::get<1>(ret)[0], std::get<1>(ret)[1], std::get<0>(ret)[2]));
    });
}




const char* build_octree_from_pointcloud_internal_doc = R"Qu8mg5v7(

)Qu8mg5v7";
npe_function(build_octree_from_pointcloud_internal)
npe_arg(octree, std::shared_ptr<Octree>)
npe_arg(points, dense_float, dense_double)
npe_arg(pad_amount, double)
npe_doc(build_octree_from_pointcloud_internal_doc)
npe_begin_code()
{
    octree->ConvertFromPointCloud(points, pad_amount);
}
npe_end_code()


const char* insert_points_into_octree_internal_doc = R"Qu8mg5v7(

)Qu8mg5v7";
npe_function(insert_points_into_octree_internal)
npe_arg(octree, std::shared_ptr<Octree>)
npe_arg(points, dense_float, dense_double)
npe_arg(base_index, int)
npe_doc(insert_points_into_octree_internal_doc)
npe_begin_code()
{
    for (int i = 0; i < points.rows(); i += 1) {
        octree->InsertPoint(points.row(i).template cast<double>(),
                            OctreePointColorLeafNode::GetInitFunction(),
                            OctreePointColorLeafNode::GetUpdateFunction(base_index + i, Eigen::Vector3d(0.0, 0.0, 0.0)),
                            OctreeInternalPointNode::GetInitFunction(),
                            OctreeInternalPointNode::GetUpdateFunction(base_index + i));
    }
}
npe_end_code()


const char* get_octree_point_leaves_internal_doc = R"Qu8mg5v7(

)Qu8mg5v7";
npe_function(get_octree_point_leaves_internal)
npe_arg(octree, std::shared_ptr<Octree>)
npe_arg(points, dense_float, dense_double)
npe_doc(get_octree_point_leaves_internal_doc)
npe_begin_code()
{
    std::vector<pybind11::tuple> ret_info(points.rows());
    std::vector<pybind11::list> ret_idxs(points.rows());

    #pragma omp parallel_for
    for (int i = 0; i < points.rows(); i += 1) {
        auto ret_i = octree->LocateLeafNode(points.row(i).template cast<double>());
        if (ret_i.second == nullptr) {
            ret_info[i] = pybind11::none();
            ret_idxs[i] = pybind11::list();
        } else {
            ret_info[i] = pybind11::cast(
                    std::make_tuple(npe::move(ret_i.second->origin_), ret_i.second->size_,
                                    ret_i.second->depth_, ret_i.second->child_index_));
            ret_info[i] = pybind11::cast(
                    std::dynamic_pointer_cast<OctreePointColorLeafNode>(ret_i.first)->indices_);
        }
    }
    return std::make_tuple(ret_info, ret_idxs);
}
npe_end_code()



const char* get_octree_point_depths_internal_doc = R"Qu8mg5v7(

)Qu8mg5v7";
npe_function(get_octree_point_depths_internal)
npe_arg(octree, std::shared_ptr<Octree>)
npe_arg(points, dense_float, dense_double)
npe_doc(get_octree_point_depths_internal_doc)
npe_begin_code()
{
    Eigen::VectorXi ret_depth(points.rows());

    #pragma omp parallel_for
    for (int i = 0; i < points.rows(); i += 1) {
        auto ret_i = octree->LocateLeafNode(points.row(i).template cast<double>());
        if (ret_i.second == nullptr) {
            ret_depth[i] = -1;
        } else {
            ret_depth[i] = ret_i.second->depth_;
        }
    }
    return npe::move(ret_depth);
}
npe_end_code()



