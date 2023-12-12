#include <npe.h>
#include <pybind11/stl.h>

#include <igl/adjacency_list.h>

#include "common/common.h"



const char* adjacency_list_doc = R"igl_Qu8mg5v7(
Compute the adjacency list given a set of mesh faces.

Args:
    f : \#f by 3 Matrix of face (triangle) indices

Returns:
    adj_list : a list of lists such that adj_list[i] contains the indexes of vertices adjacent to vertex i

)igl_Qu8mg5v7";
npe_function(adjacency_list)
npe_doc(adjacency_list_doc)
npe_arg(f, dense_int32, dense_int64)
npe_begin_code()
{
    validate_mesh_faces(f);

    std::vector<std::vector<int>> A;

    igl::adjacency_list(f, A);

    return A;
}
npe_end_code()