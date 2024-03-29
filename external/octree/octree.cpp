// ----------------------------------------------------------------------------
// -                        Open3D: www.open3d.org                            -
// ----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2018 www.open3d.org
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ----------------------------------------------------------------------------

#include "octree.h"

#include <Eigen/Dense>
#include <algorithm>
#include <unordered_map>

//#include "open3d/geometry/PointCloud.h"
//#include "open3d/geometry/VoxelGrid.h"


//
// OctreeInternalNode
//
std::shared_ptr<OctreeNodeInfo> OctreeInternalNode::GetInsertionNodeInfo(
        const std::shared_ptr<OctreeNodeInfo>& node_info,
        const Eigen::Vector3d& point) {
    if (!Octree::IsPointInBound(point, node_info->origin_, node_info->size_)) {
        throw std::runtime_error(
                "Internal error: cannot insert to child since point not in "
                "parent node bound.");
    }

    double child_size = node_info->size_ / 2.0;
    size_t x_index = point(0) < node_info->origin_(0) + child_size ? 0 : 1;
    size_t y_index = point(1) < node_info->origin_(1) + child_size ? 0 : 1;
    size_t z_index = point(2) < node_info->origin_(2) + child_size ? 0 : 1;
    size_t child_index = x_index + y_index * 2 + z_index * 4;
    Eigen::Vector3d child_origin =
            node_info->origin_ + Eigen::Vector3d(x_index * child_size,
                                                 y_index * child_size,
                                                 z_index * child_size);
    auto child_node_info = std::make_shared<OctreeNodeInfo>(
            child_origin, child_size, node_info->depth_ + 1, child_index);
    return child_node_info;
}

std::function<std::shared_ptr<OctreeInternalNode>()>
OctreeInternalNode::GetInitFunction() {
    return []() -> std::shared_ptr<OctreeInternalNode> {
        return std::make_shared<OctreeInternalNode>();
    };
}

std::function<void(std::shared_ptr<OctreeInternalNode>)>
OctreeInternalNode::GetUpdateFunction() {
    return [](std::shared_ptr<OctreeInternalNode> node) -> void {};
}



//
// OctreeInternalPointNode
//
std::function<std::shared_ptr<OctreeInternalNode>()>
OctreeInternalPointNode::GetInitFunction() {
    return []() -> std::shared_ptr<OctreeInternalNode> {
        return std::make_shared<OctreeInternalPointNode>();
    };
}

std::function<void(std::shared_ptr<OctreeInternalNode>)>
OctreeInternalPointNode::GetUpdateFunction(size_t idx) {
    // Here the captured "idx" cannot be a reference, must be a copy,
    // otherwise pybind does not have the correct value
    return [idx](std::shared_ptr<OctreeInternalNode> node) -> void {
        if (auto internal_point_node = std::dynamic_pointer_cast<
                OctreeInternalPointNode>(node)) {
            internal_point_node->indices_.push_back(idx);
        } else {
            throw std::runtime_error(
                    "Internal error: internal node must be "
                    "OctreeInternalPointNode");
        }
    };
}



//
// OctreeColorLeafNode
//
std::function<std::shared_ptr<OctreeLeafNode>()>
OctreeColorLeafNode::GetInitFunction() {
    return []() -> std::shared_ptr<OctreeLeafNode> {
        return std::make_shared<OctreeColorLeafNode>();
    };
}

std::function<void(std::shared_ptr<OctreeLeafNode>)>
OctreeColorLeafNode::GetUpdateFunction(const Eigen::Vector3d& color) {
    // Here the captured "color" cannot be a reference, must be a copy,
    // otherwise pybind does not have the correct value
    return [color](std::shared_ptr<OctreeLeafNode> node) -> void {
        if (auto color_leaf_node = std::dynamic_pointer_cast<OctreeColorLeafNode>(node)) {
            color_leaf_node->color_ = color;
        } else {
            throw std::runtime_error("Internal error: leaf node must be OctreeColorLeafNode");
        }
    };
}

std::shared_ptr<OctreeLeafNode> OctreeColorLeafNode::Clone() const {
    auto cloned_node = std::make_shared<OctreeColorLeafNode>();
    cloned_node->color_ = color_;
    return cloned_node;
}

bool OctreeColorLeafNode::operator==(const OctreeLeafNode& that) const {
    try {
        const OctreeColorLeafNode& that_color_node =
                dynamic_cast<const OctreeColorLeafNode&>(that);
        return this->color_.isApprox(that_color_node.color_);
    } catch (const std::exception&) {
        return false;
    }
}





//
// OctreePointColorLeafNode
//
std::function<std::shared_ptr<OctreeLeafNode>()>
OctreePointColorLeafNode::GetInitFunction() {
    return []() -> std::shared_ptr<OctreeLeafNode> {
        return std::make_shared<OctreePointColorLeafNode>();
    };
}

std::function<void(std::shared_ptr<OctreeLeafNode>)>
OctreePointColorLeafNode::GetUpdateFunction(size_t idx,
                                            const Eigen::Vector3d& color) {
    // Here the captured "idx" cannot be a reference, must be a copy,
    // otherwise pybind does not have the correct value
    return [idx,
            color](std::shared_ptr<OctreeLeafNode> node) -> void {
        if (auto point_color_leaf_node = std::dynamic_pointer_cast<OctreePointColorLeafNode>(node)) {
            OctreeColorLeafNode::GetUpdateFunction(color)(point_color_leaf_node);
            point_color_leaf_node->indices_.push_back(idx);
        } else {
            throw std::runtime_error(
                    "Internal error: leaf node must be "
                    "OctreePointColorLeafNode");
        }
    };
}

std::shared_ptr<OctreeLeafNode> OctreePointColorLeafNode::Clone() const {
    auto cloned_node = std::make_shared<OctreePointColorLeafNode>();
    cloned_node->color_ = color_;
    cloned_node->indices_ = indices_;
    return cloned_node;
}

bool OctreePointColorLeafNode::operator==(const OctreeLeafNode& that) const {
    try {
        const OctreePointColorLeafNode& that_point_color_node =
                dynamic_cast<const OctreePointColorLeafNode&>(that);

        return this->color_.isApprox(that_point_color_node.color_) &&
               this->indices_.size() == that_point_color_node.indices_.size() &&
               this->indices_ == that_point_color_node.indices_;
    } catch (const std::exception&) {
        return false;
    }
}



Octree::Octree(const Octree& src_octree)
        : origin_(src_octree.origin_),
          size_(src_octree.size_),
          max_depth_(src_octree.max_depth_) {
    // First traversal: clone nodes without edges
    std::unordered_map<std::shared_ptr<OctreeNode>, std::shared_ptr<OctreeNode>>
            map_src_to_dst_node;
    auto f_build_map =
            [&map_src_to_dst_node](
                    const std::shared_ptr<OctreeNode>& src_node,
                    const std::shared_ptr<OctreeNodeInfo>& src_node_info)
                    -> bool {
                if (auto src_internal_node =
                        std::dynamic_pointer_cast<OctreeInternalNode>(src_node)) {
                    auto dst_internal_node = std::make_shared<OctreeInternalNode>();
                    map_src_to_dst_node[src_internal_node] = dst_internal_node;
                } else if (auto src_leaf_node =
                        std::dynamic_pointer_cast<OctreeLeafNode>(
                                src_node)) {
                    map_src_to_dst_node[src_leaf_node] = src_leaf_node->Clone();
                } else {
                    throw std::runtime_error("Internal error: unknown node type");
                }
                return false;
            };
    src_octree.Traverse(f_build_map);

    // Second traversal: add edges
    auto f_clone_edges =
            [&map_src_to_dst_node](
                    const std::shared_ptr<OctreeNode>& src_node,
                    const std::shared_ptr<OctreeNodeInfo>& src_node_info)
                    -> bool {
                if (auto src_internal_node =
                        std::dynamic_pointer_cast<OctreeInternalNode>(src_node)) {
                    auto dst_internal_node =
                            std::dynamic_pointer_cast<OctreeInternalNode>(
                                    map_src_to_dst_node.at(src_internal_node));
                    for (size_t child_index = 0; child_index < 8; child_index++) {
                        auto src_child_node = src_internal_node->children_[child_index];
                        if (src_child_node != nullptr) {
                            auto dst_child_node =
                                    map_src_to_dst_node.at(src_child_node);
                            dst_internal_node->children_[child_index] = dst_child_node;
                        }
                    }
                }
                return false;
            };
    src_octree.Traverse(f_clone_edges);

    // Save root node to dst_octree (this octree)
    root_node_ = map_src_to_dst_node.at(src_octree.root_node_);
}

bool Octree::operator==(const Octree& that) const {
    // Check basic properties
    bool rc = true;
    rc = rc && origin_.isApprox(that.origin_);
    rc = rc && size_ == that.size_;
    rc = rc && max_depth_ == that.max_depth_;
    if (!rc) {
        return rc;
    }

    // Assign and check node ids
    std::unordered_map<std::shared_ptr<OctreeNode>, size_t> map_node_to_id;
    std::unordered_map<size_t, std::shared_ptr<OctreeNode>> map_id_to_node;
    size_t next_id = 0;
    auto f_assign_node_id =
            [&map_node_to_id, &map_id_to_node, &next_id](
                    const std::shared_ptr<OctreeNode>& node,
                    const std::shared_ptr<OctreeNodeInfo>& node_info) -> bool {
                map_node_to_id[node] = next_id;
                map_id_to_node[next_id] = node;
                next_id++;
                return false;
            };

    map_node_to_id.clear();
    map_id_to_node.clear();
    next_id = 0;
    Traverse(f_assign_node_id);
    std::unordered_map<std::shared_ptr<OctreeNode>, size_t>
            this_map_node_to_id = map_node_to_id;
    std::unordered_map<size_t, std::shared_ptr<OctreeNode>>
            this_map_id_to_node = map_id_to_node;
    size_t num_nodes = next_id;

    map_node_to_id.clear();
    map_id_to_node.clear();
    next_id = 0;
    that.Traverse(f_assign_node_id);
    std::unordered_map<std::shared_ptr<OctreeNode>, size_t>
            that_map_node_to_id = map_node_to_id;
    std::unordered_map<size_t, std::shared_ptr<OctreeNode>>
            that_map_id_to_node = map_id_to_node;

    rc = rc && this_map_node_to_id.size() == num_nodes &&
         that_map_node_to_id.size() == num_nodes &&
         this_map_id_to_node.size() == num_nodes &&
         that_map_id_to_node.size() == num_nodes;
    if (!rc) {
        return rc;
    }

    // Check nodes
    for (size_t id = 0; id < num_nodes; ++id) {
        std::shared_ptr<OctreeNode> this_node = this_map_id_to_node.at(id);
        std::shared_ptr<OctreeNode> that_node = that_map_id_to_node.at(id);
        // Check if internal node
        auto is_same_node_type = false;
        auto this_internal_node =
                std::dynamic_pointer_cast<OctreeInternalNode>(this_node);
        auto that_internal_node =
                std::dynamic_pointer_cast<OctreeInternalNode>(that_node);
        if (this_internal_node != nullptr && that_internal_node != nullptr) {
            is_same_node_type = true;
            for (size_t child_index = 0; child_index < 8; child_index++) {
                const std::shared_ptr<OctreeNode>& this_child =
                        this_internal_node->children_[child_index];
                int this_child_id = -1;
                if (this_child != nullptr) {
                    this_child_id = int(this_map_node_to_id.at(this_child));
                }
                const std::shared_ptr<OctreeNode>& that_child =
                        that_internal_node->children_[child_index];
                int that_child_id = -1;
                if (that_child != nullptr) {
                    that_child_id = int(that_map_node_to_id.at(that_child));
                }
                rc = rc && this_child_id == that_child_id;
            }
        }
        // Check if leaf node
        auto this_leaf_node =
                std::dynamic_pointer_cast<OctreeLeafNode>(this_node);
        auto that_leaf_node =
                std::dynamic_pointer_cast<OctreeLeafNode>(that_node);
        if (this_leaf_node != nullptr && that_leaf_node != nullptr) {
            is_same_node_type = true;
            rc = rc && (*this_leaf_node) == (*that_leaf_node);
        }
        // Handle case where node types are different
        rc = rc && is_same_node_type;
    }

    return rc;
}

Octree& Octree::Clear() {
    // Inherited Clear function
    root_node_ = nullptr;
    origin_.setZero();
    size_ = 0;
    return *this;
}

bool Octree::IsEmpty() const { return root_node_ == nullptr; }

Eigen::Vector3d Octree::GetMinBound() const {
    if (IsEmpty()) {
        return Eigen::Vector3d::Zero();
    } else {
        return origin_;
    }
}

Eigen::Vector3d Octree::GetMaxBound() const {
    if (IsEmpty()) {
        return Eigen::Vector3d::Zero();
    } else {
        return origin_ + Eigen::Vector3d(size_, size_, size_);
    }
}

Eigen::Vector3d Octree::GetCenter() const {
    return origin_ + Eigen::Vector3d(size_, size_, size_) / 2;
}

std::tuple<Eigen::Vector3d, Eigen::Vector3d> Octree::GetAxisAlignedBoundingBox() const {

    Eigen::Vector3d min_b = GetMinBound();
    Eigen::Vector3d max_b = GetMaxBound();
    return std::make_tuple(min_b, max_b);
}

// FRANCIS TODO: Return pair of vector3ds here
//OrientedBoundingBox Octree::GetOrientedBoundingBox() const {
//    return OrientedBoundingBox::CreateFromAxisAlignedBoundingBox(
//            GetAxisAlignedBoundingBox());
//}

// None of this is implemented. Will probably kill it
Octree& Octree::Transform(const Eigen::Matrix4d& transformation) {
    throw std::runtime_error("Not implemented");
    return *this;
}

Octree& Octree::Translate(const Eigen::Vector3d& translation, bool relative) {
    throw std::runtime_error("Not implemented");
    return *this;
}

Octree& Octree::Scale(const double scale, const Eigen::Vector3d& center) {
    throw std::runtime_error("Not implemented");
    return *this;
}

Octree& Octree::Rotate(const Eigen::Matrix3d& R,
                       const Eigen::Vector3d& center) {
    throw std::runtime_error("Not implemented");
    return *this;
}


template <typename DerivedP>
void Octree::ConvertFromPointCloud(const Eigen::MatrixBase<DerivedP>& point_cloud,
                                   double size_expand) {
    if (size_expand > 1 || size_expand < 0) {
        pybind11::value_error("pad_amount should be between 0 and 1");
    }

    // Set bounds
    Clear();
    Eigen::Array3d min_bound = point_cloud.rowwise().minCoeff().eval().template cast<double>().array();
    Eigen::Array3d max_bound = point_cloud.rowwise().maxCoeff().eval().template cast<double>().array();
    Eigen::Array3d center = (min_bound + max_bound) / 2;
    Eigen::Array3d half_sizes = center - min_bound;
    double max_half_size = half_sizes.maxCoeff();
    origin_ = min_bound.min(center - max_half_size);
    if (max_half_size == 0) {
        size_ = size_expand;
    } else {
        size_ = max_half_size * 2 * (1 + size_expand);
    }

    // Insert points
    for (size_t idx = 0; idx < point_cloud.rows(); idx++) {
        const Eigen::Vector3d& color = Eigen::Vector3d::Zero();
        InsertPoint(point_cloud.row(idx).template cast<double>(),
                    OctreePointColorLeafNode::GetInitFunction(),
                    OctreePointColorLeafNode::GetUpdateFunction(idx, color),
                    OctreeInternalPointNode::GetInitFunction(),
                    OctreeInternalPointNode::GetUpdateFunction(idx));
    }
}

void Octree::InsertPoint(
        const Eigen::Vector3d& point,
        const std::function<std::shared_ptr<OctreeLeafNode>()>& fl_init,
        const std::function<void(std::shared_ptr<OctreeLeafNode>)>& fl_update,
        const std::function<std::shared_ptr<OctreeInternalNode>()>& fi_init,
        const std::function<void(std::shared_ptr<OctreeInternalNode>)>& fi_update) {
    // if missing, create basic internal node init and update functions
    auto _fi_init = fi_init;
    if (_fi_init == nullptr) {
        _fi_init = OctreeInternalNode::GetInitFunction();
    }
    auto _fi_update = fi_update;
    if (_fi_update == nullptr) {
        _fi_update = OctreeInternalNode::GetUpdateFunction();
    }

    if (root_node_ == nullptr) {
        if (max_depth_ == 0) {
            root_node_ = fl_init();
        } else {
            root_node_ = _fi_init();
        }
    }
    auto root_node_info =
            std::make_shared<OctreeNodeInfo>(origin_, size_, 0, 0);

    InsertPointRecurse(root_node_, root_node_info, point, fl_init, fl_update,
                       _fi_init, _fi_update);
}

void Octree::InsertPointRecurse(
        const std::shared_ptr<OctreeNode>& node,
        const std::shared_ptr<OctreeNodeInfo>& node_info,
        const Eigen::Vector3d& point,
        const std::function<std::shared_ptr<OctreeLeafNode>()>& fl_init,
        const std::function<void(std::shared_ptr<OctreeLeafNode>)>& fl_update,
        const std::function<std::shared_ptr<OctreeInternalNode>()>& fi_init,
        const std::function<void(std::shared_ptr<OctreeInternalNode>)>&
        fi_update) {
    if (!IsPointInBound(point, node_info->origin_, node_info->size_)) {
        return;
    }
    if (node_info->depth_ > max_depth_) {
        return;
    } else if (node_info->depth_ == max_depth_) {
        if (auto leaf_node = std::dynamic_pointer_cast<OctreeLeafNode>(node)) {
            fl_update(leaf_node);
        } else {
            throw std::runtime_error("Internal error: leaf node must be OctreeLeafNode");
        }
    } else {
        if (auto internal_node = std::dynamic_pointer_cast<OctreeInternalNode>(node)) {
            // Update internal node with information about the current point
            fi_update(internal_node);

            // Get child node info
            std::shared_ptr<OctreeNodeInfo> child_node_info =
                    internal_node->GetInsertionNodeInfo(node_info, point);

            // Create child node with factory function
            size_t child_index = child_node_info->child_index_;
            if (internal_node->children_[child_index] == nullptr) {
                if (node_info->depth_ == max_depth_ - 1) {
                    internal_node->children_[child_index] = fl_init();
                } else {
                    internal_node->children_[child_index] = fi_init();
                }
            }

            // Insert to the child
            InsertPointRecurse(internal_node->children_[child_index],
                               child_node_info, point, fl_init, fl_update,
                               fi_init, fi_update);
        } else {
            throw std::runtime_error("Internal error: internal node must be OctreeInternalNode");
        }
    }
}

bool Octree::IsPointInBound(const Eigen::Vector3d& point,
                            const Eigen::Vector3d& origin,
                            const double& size) {
    bool rc = (Eigen::Array3d(origin) <= Eigen::Array3d(point)).all() &&
              (Eigen::Array3d(point) < Eigen::Array3d(origin) + size).all();
    return rc;
}

void Octree::Traverse(
        const std::function<bool(const std::shared_ptr<OctreeNode>&,
                                 const std::shared_ptr<OctreeNodeInfo>&)>& f) {
    // root_node_'s child index is 0, though it isn't a child node
    TraverseRecurse(root_node_,
                    std::make_shared<OctreeNodeInfo>(origin_, size_, 0, 0), f);
}

void Octree::Traverse(
        const std::function<bool(const std::shared_ptr<OctreeNode>&,
                                 const std::shared_ptr<OctreeNodeInfo>&)>& f)
const {
    // root_node_'s child index is 0, though it isn't a child node
    TraverseRecurse(root_node_,
                    std::make_shared<OctreeNodeInfo>(origin_, size_, 0, 0), f);
}

void Octree::TraverseRecurse(
        const std::shared_ptr<OctreeNode>& node,
        const std::shared_ptr<OctreeNodeInfo>& node_info,
        const std::function<bool(const std::shared_ptr<OctreeNode>&,
                                 const std::shared_ptr<OctreeNodeInfo>&)>& f) {
    if (node == nullptr) {
        return;
    } else if (auto internal_node =
            std::dynamic_pointer_cast<OctreeInternalNode>(node)) {
        // Allow caller to avoid traversing further down this tree path
        if (f(internal_node, node_info)) return;

        double child_size = node_info->size_ / 2.0;

        for (size_t child_index = 0; child_index < 8; ++child_index) {
            size_t x_index = child_index % 2;
            size_t y_index = (child_index / 2) % 2;
            size_t z_index = (child_index / 4) % 2;

            auto child_node = internal_node->children_[child_index];
            Eigen::Vector3d child_node_origin =
                    node_info->origin_ + Eigen::Vector3d(double(x_index),
                                                         double(y_index),
                                                         double(z_index)) *
                                         child_size;
            auto child_node_info = std::make_shared<OctreeNodeInfo>(
                    child_node_origin, child_size, node_info->depth_ + 1,
                    child_index);
            TraverseRecurse(child_node, child_node_info, f);
        }
    } else if (auto leaf_node =
            std::dynamic_pointer_cast<OctreeLeafNode>(node)) {
        f(leaf_node, node_info);
    } else {
        throw std::runtime_error("Internal error: unknown node type");
    }
}

std::pair<std::shared_ptr<OctreeLeafNode>, std::shared_ptr<OctreeNodeInfo>>
Octree::LocateLeafNode(const Eigen::Vector3d& point) const {
    std::shared_ptr<OctreeLeafNode> target_leaf_node = nullptr;
    std::shared_ptr<OctreeNodeInfo> target_leaf_node_info = nullptr;
    auto f_locate_leaf_node =
            [&target_leaf_node, &target_leaf_node_info, &point](
                    const std::shared_ptr<OctreeNode>& node,
                    const std::shared_ptr<OctreeNodeInfo>& node_info) -> bool {
                bool skip_children = false;
                if (IsPointInBound(point, node_info->origin_, node_info->size_)) {
                    if (auto leaf_node =
                            std::dynamic_pointer_cast<OctreeLeafNode>(node)) {
                        target_leaf_node = leaf_node;
                        target_leaf_node_info = node_info;
                    }
                } else {
                    skip_children = true;
                }
                return skip_children;
            };
    Traverse(f_locate_leaf_node);
    return std::make_pair(target_leaf_node, target_leaf_node_info);
}

