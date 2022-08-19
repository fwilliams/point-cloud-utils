#include <npe.h>

#include <igl/bfs_orient.h>

#include <tuple>

#include "common/common.h"


const char* orient_mesh_faces_doc = R"igl_Qu8mg5v7(
Consistently orient faces of a mesh within each connected component

Parameters
----------
f : (#f, 3)-shaped NumPy array of face (triangle) indices

Returns
-------
oriented_faces  : (#f, 3)-shaped NumPy array of faces which are consistently oriented
face_components : (#f,)-shaped NumPy array of connected component ids 
                  (i.e. face_components[i] is the component id of facef[i]) 

See also
--------

)igl_Qu8mg5v7";
npe_function(orient_mesh_faces)
npe_doc(orient_mesh_faces_doc)
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_default_arg(weighting_type, std::string, std::string("uniform"))
npe_begin_code()
{
    validate_mesh_faces(f);

    npe_Matrix_f oriented_faces;
    npe_Matrix_f face_components;
    igl::bfs_orient(f, oriented_faces, face_components);

    return std::make_tuple(npe::move(oriented_faces), npe::move(face_components));
}
npe_end_code()