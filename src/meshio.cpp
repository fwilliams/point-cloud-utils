#include <igl/readOBJ.h>
#include <igl/readPLY.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <igl/writePLY.h>
#include <igl/writeOFF.h>

#include <unordered_map>

#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <npe.h>
#include <npe_typedefs.h>
#include <pybind11/stl.h>

#include "common/common.h"

namespace {

void runtime_warning(const std::string& message, Py_ssize_t level=1) {
    PyErr_WarnEx(PyExc_RuntimeWarning, message.c_str(), level);
}

using namespace vcg;

class CVertex;
class CFace;

struct MyTypes: public UsedTypes< Use<CVertex>::AsVertexType, Use<CFace>::AsFaceType >{};

class CVertex  : public Vertex< MyTypes,
        vertex::Coord3d, vertex::BitFlags, vertex::Normal3d, vertex::Color4b, vertex::TexCoord2d, vertex::Qualityd, vertex::Radiusf > {};
class CFace    : public Face< MyTypes,
        face::VertexRef, face::BitFlags, face::Color4b, face::Qualityd, face::Normal3d, face::WedgeColor4b, face::WedgeTexCoord2d, face::WedgeNormal > {};
class CMesh    : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace> > {};

typedef CMesh::VertexPointer VertexPointer;
typedef CMesh::VertexIterator VertexIterator;
typedef Point3<CMesh::ScalarType> Point3x;
typedef std::vector<Point3x> Hole;

typedef CMesh::VertexPointer VertexPointer;
typedef CMesh::VertexIterator VertexIterator;
typedef CMesh::FaceContainer FaceContainer;
typedef CMesh::ScalarType ScalarType;


bool assert_shape_and_dtype(const pybind11::array& arr, std::string name, pybind11::dtype dtype,
                            const std::vector<ssize_t>& shape) {
    if (!arr.dtype().is(dtype)) {
        throw pybind11::value_error("Invalid dtype for argument '" + name + "'. Expected '" +
                                    dtype.kind() + "' but got '" + arr.dtype().kind() + "'.");
    }
    if (shape.size() != arr.ndim()) {
        throw pybind11::value_error("Invalid number of dimensions for argument '" + name + "'. Expected " +
                                    std::to_string(shape.size()) + " but got " + std::to_string(arr.ndim()) + ".");
    }
    bool nonempty = true;
    for (int i = 0; i < shape.size(); i++) {
        if (arr.shape()[i] <= 0) {
            nonempty = false;
        }
        if (shape[i] < 0) {
            if (arr.shape()[i] == 0) {
                continue;
            } else if (arr.shape()[i] == -shape[i]) {
                continue;
            }
        } else if (shape[i] == arr.shape()[i]) {
            continue;
        }

        throw pybind11::value_error("Invalid  shape for argument '" + name + "' at dimension " +
                                    std::to_string(i) + ". Expected " + std::to_string(shape[i]) +
                                    " but got " + std::to_string(arr.shape()[i]) + ".");
    }

    return nonempty;
}

template <typename Scalar>
void load_mesh_vcg(CMesh& m, int mask, std::unordered_map<std::string, pybind11::object>& ret) {
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 3*2, Eigen::RowMajor> FMatrix32;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 3*3, Eigen::RowMajor> FMatrix33;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 3*4, Eigen::RowMajor> FMatrix34;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 4, Eigen::RowMajor> FMatrix4;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor> FMatrix3;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 2, Eigen::RowMajor> FMatrix2;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> FMatrix1;
    typedef Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> IMatrix3;
    typedef Eigen::Matrix<int, Eigen::Dynamic, 1> IMatrix1;

    const int num_vertices = m.vn;
    const int num_faces = m.fn;

    // Per-vertex attributes
    FMatrix3 vertex_pos, vertex_normals;
    FMatrix2 vertex_texcoord;
    IMatrix1 vertex_texindex;
    FMatrix1 vertex_quality, vertex_radius;
    FMatrix4 vertex_colors;
    IMatrix1 vertex_flags;

    IMatrix3 face_indices;
    IMatrix1 face_flags;
    FMatrix4 face_colors;
    FMatrix1 face_quality;
    FMatrix3 face_normals;

    FMatrix34 wedge_colors;
    FMatrix32 wedge_texcoords;
    IMatrix3 wedge_texindex;
    FMatrix33 wedge_normals;

    const bool has_v_coord = true;  // We always assume there are spatial vertex coordinates
    const bool has_v_flags = mask & tri::io::Mask::IOM_VERTFLAGS;
    const bool has_v_color = mask & tri::io::Mask::IOM_VERTCOLOR;
    const bool has_v_quality = mask & tri::io::Mask::IOM_VERTQUALITY;
    const bool has_v_normal = mask & tri::io::Mask::IOM_VERTNORMAL;
    const bool has_v_texcoord = mask & tri::io::Mask::IOM_VERTTEXCOORD;
    const bool has_v_radius = mask & tri::io::Mask::IOM_VERTRADIUS;

    const bool has_f_index = true; //mask & tri::io::Mask::IOM_FACEINDEX;
    const bool has_f_flags = mask & tri::io::Mask::IOM_FACEFLAGS;
    const bool has_f_color = mask & tri::io::Mask::IOM_FACECOLOR;
    const bool has_f_quality = mask & tri::io::Mask::IOM_FACEQUALITY;
    const bool has_f_normal = mask & tri::io::Mask::IOM_FACENORMAL;

    const bool has_w_color = mask & vcg::tri::io::Mask::IOM_WEDGCOLOR;
    const bool has_w_texcoord = mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD;
    const bool has_w_texindex = mask & vcg::tri::io::Mask::IOM_WEDGTEXMULTI;
    const bool has_w_normal = mask & vcg::tri::io::Mask::IOM_WEDGNORMAL;

    {
        if (has_v_coord) {
//            std::cout << "has vertex position" << std::endl;
            vertex_pos.resize(num_vertices, 3);
        }
        if (has_v_flags) {
//            std::cout << "has vertex flags" << std::endl;
            vertex_flags.resize(num_vertices, 1);
        }
        if (has_v_color) {
//            std::cout << "has vertex colors" << std::endl;
            vertex_colors.resize(num_vertices, 4);
        }
        if (has_v_quality) {
//            std::cout << "has vertex quality" << std::endl;
            vertex_quality.resize(num_vertices, 1);
        }
        if (has_v_normal) {
//            std::cout << "has vertex normals" << std::endl;
            vertex_normals.resize(num_vertices, 3);
        }
        if (has_v_texcoord) {
//            std::cout << "has vertex texcoords" << std::endl;
            vertex_texcoord.resize(num_vertices, 2);
            if (m.textures.size() > 0 || m.normalmaps.size() > 0) {
                vertex_texindex.resize(num_vertices, 1);
            }
        }
        if (has_v_radius) {
//            std::cout << "has vertex radius" << std::endl;
            vertex_radius.resize(num_vertices, 1);
        }
    }


    {
        if (has_f_index) {
//            std::cout << "has face index" << std::endl;
            face_indices.resize(num_faces, 3);
        }
        if (has_f_flags) {
//            std::cout << "has face flags" << std::endl;
            face_flags.resize(num_faces, 3);
        }
        if (has_f_color) {
//            std::cout << "has face colors" << std::endl;
            face_colors.resize(num_faces, 4);
        }
        if (has_f_quality) {
//            std::cout << "has facequality" << std::endl;
            face_quality.resize(num_faces, 1);
        }
        if (has_f_normal) {
//            std::cout << "has face normals" << std::endl;
            face_normals.resize(num_faces, 3);
        }
    }


    {
        if (has_w_color) {
//            std::cout << "has wedge color" << std::endl;
            wedge_colors.resize(num_faces, 3 * 4);
        }
        if (has_w_texcoord) {
//            std::cout << "has wedge texcoord" << std::endl;
            wedge_texcoords.resize(num_faces, 3 * 2);
            if (m.textures.size() > 0 || m.normalmaps.size() > 0) {
                wedge_texindex.resize(num_faces, 3);
            }
        }
        if (has_w_normal) {
//            std::cout << "has wedge normal" << std::endl;
            wedge_normals.resize(num_faces, 3 * 3);
        }
    }

    int vcount = 0;
    for (CMesh::VertexIterator it = m.vert.begin(); it != m.vert.end(); it++) {
        if (has_v_coord) {
            for (int i = 0; i < 3; i++) { vertex_pos(vcount, i) = it->cP()[i]; }
        }
        if (has_v_flags) {
            vertex_flags(vcount, 0) = it->Flags();
        }
        if (has_v_color) {
            for (int i = 0; i < 4; i++) { vertex_colors(vcount, i) = Scalar(it->cC()[i]) / 255.0; }
        }
        if (has_v_quality) {
            vertex_quality(vcount, 0) = it->cQ();
        }
        if (has_v_normal) {
            for (int i = 0; i < 3; i++) { vertex_normals(vcount, i) = it->cN()[i]; }
        }
        if (has_v_texcoord) {
            vertex_texcoord(vcount, 0) = it->cT().U();
            vertex_texcoord(vcount, 1) = it->cT().V();
            if (m.textures.size() > 0) {
                int texindex = int(it->cT().N());
                if (texindex < 0 || texindex >= m.textures.size()) {
                    texindex = -1;
                }
                vertex_texindex(vcount, 0) = texindex;
            }
        }
        if (has_v_radius) {
            vertex_radius(vcount, 0) = it->cR();
        }

        vcount += 1;
    }

    int fcount = 0;
    for (CMesh::FaceIterator it = m.face.begin(); it != m.face.end(); it++) {
        if (has_f_index) {
            for (int i = 0; i < 3; i++) { face_indices(fcount, i) = int(it->cV(i) - &(*m.vert.begin())); }
        }
        if (has_f_flags) {
            face_flags(fcount, 0) = it->cFlags();
        }
        if (has_f_color) {
            for (int i = 0; i < 4; i++) { face_colors(fcount, i) = Scalar(it->cC()[i]) / 255.0; }
        }
        if (has_f_quality) {
            face_quality(fcount, 0) = it->Q();
        }
        if (has_f_normal) {
            for (int i = 0; i < 3; i++) { face_normals(fcount, i) = it->cN()[i]; }
        }

        if (has_w_color) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 4; j++) {
                    wedge_colors(fcount, i * 3 + j) = Scalar(it->cWC(i)[j]) / 255.0;
                }
            }
        }
        if (has_w_texcoord) {
            for (int i = 0; i < 3; i++) {
                wedge_texcoords(fcount, i + 0) = it->cWT(i).U();
                wedge_texcoords(fcount, i + 1) = it->cWT(i).V();
                if (m.textures.size() > 0) {
                    int texindex = int(it->cWT(i).N());
                    if (texindex < 0 || texindex >= m.textures.size()) {
                        texindex = -1;
                    }
                    wedge_texindex(fcount, i) = texindex;
                }
            }
        }
        if (has_w_normal) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    wedge_normals(fcount, i *3 + j) = it->cWN(i)[j];
                }
            }
        }

        fcount += 1;
    }

    // Pack a big dict with everything in it. On the python side, we convert this to a class
    std::unordered_map<std::string, pybind11::object> vertex_ret;
    vertex_ret["positions"] = npe::move(vertex_pos);
    vertex_ret["flags"] = npe::move(vertex_flags);
    vertex_ret["colors"] = npe::move(vertex_colors);
    vertex_ret["quality"] = npe::move(vertex_quality);
    vertex_ret["normals"] = npe::move(vertex_normals);
    vertex_ret["texcoords"] = npe::move(vertex_texcoord);
    vertex_ret["tex_ids"] = npe::move(vertex_texindex);
    vertex_ret["radius"] = npe::move(vertex_radius);

    std::unordered_map<std::string, pybind11::object> face_ret;
    face_ret["vertex_ids"] = npe::move(face_indices);
    face_ret["flags"] = npe::move(face_flags);
    face_ret["colors"] = npe::move(face_colors);
    face_ret["quality"] = npe::move(face_quality);
    face_ret["normals"] = npe::move(face_normals);

    const int num_wedge_colors = wedge_colors.rows();
    pybind11::array ret_wedge_colors = npe::move(wedge_colors);
    ret_wedge_colors.resize({num_wedge_colors, 3, 4});
    face_ret["wedge_colors"] = ret_wedge_colors;

    const int num_wedge_texcoords = wedge_texcoords.rows();
    pybind11::array ret_wedge_texcoords = npe::move(wedge_texcoords);
    ret_wedge_texcoords.resize({num_wedge_texcoords, 3, 2});
    face_ret["wedge_texcoords"] = ret_wedge_texcoords;
    face_ret["wedge_tex_ids"] = npe::move(wedge_texindex);

    const int num_wedge_normals = wedge_normals.rows();
    pybind11::array ret_wedge_normals = npe::move(wedge_normals);
    ret_wedge_normals.resize({num_wedge_normals, 3, 3});
    face_ret["wedge_normals"] = ret_wedge_normals;

    pybind11::dict ret_vertex_data = pybind11::cast(vertex_ret);
    pybind11::dict ret_face_data = pybind11::cast(face_ret);
    pybind11::list ret_normalmaps = pybind11::cast(m.normalmaps);
    pybind11::list ret_textures = pybind11::cast(m.textures);

    ret["vertex_data"] = ret_vertex_data;
    ret["face_data"] = ret_face_data;
    ret["textures"] = ret_textures;
    ret["normal_maps"] = ret_normalmaps;

    // TODO: Handle custom vertex and face attributes
    //    std::cout << "Loaded mesh with " + std::to_string(m.vn) + " vertices" << std::endl;
    //    std::cout << "Mesh vertex attributes (" << m.vert_attr.size() << ")" << std::endl;
    //    std::cout << "Mesh face attributes (" << m.face_attr.size() << ")" << std::endl;
    //    std::cout << "Mesh edge attributes (" << m.edge_attr.size() << ")" << std::endl;
    //    for (auto attrptr : m.vert_attr) {
    //        std::cout << "  " << attrptr._name << std::endl;
    //    }
}


template <typename ScalarF, typename ScalarI>
void write_mesh_vcg(std::string filename, pybind11::array& v_positions, pybind11::array& v_normals,
                    pybind11::array& v_texcoords, pybind11::array& v_colors, pybind11::array& v_quality,
                    pybind11::array& v_radius, pybind11::array& v_texids, pybind11::array& v_flags,
                    pybind11::array& f_vertex_ids, pybind11::array& f_normals, pybind11::array& f_colors,
                    pybind11::array& f_quality, pybind11::array& f_flags,
                    pybind11::array& w_colors, pybind11::array& w_normals, pybind11::array& w_texcoords,
                    pybind11::array& w_texids,
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

    bool has_w_normals = assert_shape_and_dtype(w_normals, "w_normals", dtype_f, {-num_faces, 3, 3});
    bool has_w_colors = assert_shape_and_dtype(w_colors, "w_colors", dtype_f, {-num_faces, 3, 4});
    bool has_w_texcoords = assert_shape_and_dtype(w_texcoords, "w_texcoords", dtype_f, {-num_faces, 3, 2});
    bool has_w_texids = assert_shape_and_dtype(w_texids, "w_texids", dtype_i, {-num_faces, 3});

    int mask = tri::io::Mask::IOM_NONE;
    if (has_v_positions) {
        mask |= tri::io::Mask::IOM_VERTCOORD;
    }
    if (has_v_normals) {
        mask |= tri::io::Mask::IOM_VERTNORMAL;
    }
    if (has_v_texcoords) {
        mask |= tri::io::Mask::IOM_VERTTEXCOORD;
    }
    if (has_v_colors) {
        mask |= tri::io::Mask::IOM_VERTCOLOR;
    }
    if (has_v_quality) {
        mask |= tri::io::Mask::IOM_VERTQUALITY;
    }
    if (has_v_radius) {
        mask |= tri::io::Mask::IOM_VERTRADIUS;
    }
    if (has_v_flags) {
        mask |= tri::io::Mask::IOM_VERTFLAGS;
    }

    if (has_f_vertex_ids) {
        mask |= tri::io::Mask::IOM_FACEINDEX;
    }
    if (has_f_normals) {
        mask |= tri::io::Mask::IOM_FACENORMAL;
    }
    if (has_f_colors) {
        mask |= tri::io::Mask::IOM_FACECOLOR;
    }
    if (has_f_quality) {
        mask |= tri::io::Mask::IOM_FACEQUALITY;
    }
    if (has_f_flags) {
        mask |= tri::io::Mask::IOM_FACEFLAGS;
    }

    if (has_w_normals) {
        mask |= tri::io::Mask::IOM_WEDGNORMAL;
    }
    if (has_w_colors) {
        mask |= tri::io::Mask::IOM_WEDGCOLOR;
    }
    if (has_w_texcoords) {
        mask |= tri::io::Mask::IOM_WEDGTEXCOORD;
    }
    if (has_w_texids) {
        mask |= tri::io::Mask::IOM_WEDGTEXMULTI;
    }


    CMesh m;
    CMesh::VertexIterator vi = vcg::tri::Allocator<CMesh>::AddVertices(m, num_vertices);
    for (int i = 0; i < num_vertices; i++) {
        if (has_v_positions) {
            vi->P() = CMesh::CoordType(*((ScalarF*)v_positions.data(i, 0)),
                                       *((ScalarF*)v_positions.data(i, 1)),
                                       *((ScalarF*)v_positions.data(i, 2)));
        }
        if (has_v_normals) {
            vi->N() = CMesh::VertexType::NormalType(*((ScalarF*)v_normals.data(i, 0)),
                                                    *((ScalarF*)v_normals.data(i, 1)),
                                                    *((ScalarF*)v_normals.data(i, 2)));
        }
        if (has_v_texcoords) {
            vi->T() = CMesh::VertexType::TexCoordType(*((ScalarF*)v_texcoords.data(i, 0)),
                                                      *((ScalarF*)v_texcoords.data(i, 1)));
            if (has_v_texids) {
                vi->T().N() = *((ScalarI*)v_texids.data(i));
            }
        }
        if (has_v_colors) {
            vi->C() = CMesh::VertexType::ColorType((unsigned char) (*((ScalarF*)v_colors.data(i, 0)) * 255.0),
                                                   (unsigned char) (*((ScalarF*)v_colors.data(i, 1)) * 255.0),
                                                   (unsigned char) (*((ScalarF*)v_colors.data(i, 2)) * 255.0),
                                                       (unsigned char) (*((ScalarF*)v_colors.data(i, 3)) * 255.0));
        }
        if (has_v_quality) {
            vi->Q() = CMesh::VertexType::QualityType(*((ScalarF*) v_quality.data(i)));
        }
        if (has_v_radius) {
            vi->R() = CMesh::VertexType::QualityType(*((ScalarF*) v_radius.data(i)));
        }
        if (has_v_flags) {
            vi->Flags() = int(*((ScalarI*) v_flags.data(i)));
        }
        ++vi;
    }

    CMesh::FaceIterator fi = vcg::tri::Allocator<CMesh>::AddFaces(m, num_faces);
    for (int i = 0; i < num_faces; i++) {
        if (has_f_vertex_ids) {
            ScalarI fv1 = *((ScalarI*)f_vertex_ids.data(i, 0));
            ScalarI fv2 = *((ScalarI*)f_vertex_ids.data(i, 1));
            ScalarI fv3 = *((ScalarI*)f_vertex_ids.data(i, 2));
            if (fv1 >= m.vert.size() || fv2 >= m.vert.size() || fv3 >= m.vert.size()) {
                throw pybind11::value_error("Invalid face (" + std::to_string(fv1) + ", " + std::to_string(fv2) +
                                            ", " + std::to_string(fv3) + ") at index " + std::to_string(i) +
                                            " exceeds the number of vertices (" + std::to_string(num_vertices) + ").");
            }

            fi->V(0) = &(m.vert.at(fv1));
            fi->V(1) = &(m.vert.at(fv2));
            fi->V(2) = &(m.vert.at(fv3));
        }
        if (has_f_normals) {
            fi->N() = CMesh::FaceType::NormalType(*((ScalarF*)f_normals.data(i, 0)),
                                                  *((ScalarF*)f_normals.data(i, 1)),
                                                  *((ScalarF*)f_normals.data(i, 2)));
        }
        if (has_f_colors) {
            fi->C() = CMesh::FaceType::ColorType((unsigned char) (*((ScalarF*)f_colors.data(i, 0)) * 255.0),
                                                 (unsigned char) (*((ScalarF*)f_colors.data(i, 1)) * 255.0),
                                                 (unsigned char) (*((ScalarF*)f_colors.data(i, 2)) * 255.0),
                                                 (unsigned char) (*((ScalarF*)f_colors.data(i, 3)) * 255.0));
        }
        if (has_f_quality) {
            fi->Q() = CMesh::FaceType::QualityType(*((ScalarF*) f_quality.data(i)));
        }
        if (has_f_flags) {
            fi->Flags() = int(*((ScalarI*) f_flags.data(i)));
        }

        if (has_w_colors) {
            for (int j = 0; j < 3; j++) {
                fi->WC(j) = CMesh::FaceType::ColorType((unsigned char) (*((ScalarF*)w_colors.data(i, j, 0)) * 255.0),
                                                       (unsigned char) (*((ScalarF*)w_colors.data(i, j, 1)) * 255.0),
                                                       (unsigned char) (*((ScalarF*)w_colors.data(i, j, 2)) * 255.0),
                                                       (unsigned char) (*((ScalarF*)w_colors.data(i, j, 3)) * 255.0));
            }
        }
        if (has_w_texcoords) {
            for (int j = 0; j < 3; j++) {
                fi->WT(j) = CMesh::FaceType::TexCoordType(*((ScalarF*)w_texcoords.data(i, j, 0)),
                                                          *((ScalarF*)w_texcoords.data(i, j, 1)));
                if (has_w_texids) {
                    ScalarI texindex = *((ScalarI*)w_texids.data(i, j));
                    fi->WT(j).N() = texindex;
                }
            }
        }
        if (has_w_normals) {
            for (int j = 0; j < 3; j++) {
                fi->WN(j) = CMesh::FaceType::NormalType(*((ScalarF*)w_normals.data(i, j, 0)),
                                                        *((ScalarF*)w_normals.data(i, j, 1)),
                                                        *((ScalarF*)w_normals.data(i, j, 2)));
            }
        }

        ++fi;
    }

    int err = tri::io::Exporter<CMesh>::Save(m, filename.c_str(), mask);
    if (err) {
        throw pybind11::value_error("Error during loading " + filename + ": '" +
                                    tri::io::Exporter<CMesh>::ErrorMsg(err) + "'");
    }
}




}

const char* ds_load_mesh = R"igl_Qu8mg5v7(
Load a mesh
)igl_Qu8mg5v7";
npe_function(load_mesh_internal)
npe_doc(ds_load_mesh)
npe_arg(filename, std::string)
npe_default_arg(dtype, npe::dtype, "float64")
npe_begin_code()
{
    CMesh m;
    int mask = 0;
    tri::io::Importer<CMesh>::LoadMask(filename.c_str(), mask);

    int err = tri::io::Importer<CMesh>::Open(m, filename.c_str());
    if (err) {
        if(tri::io::Importer<CMesh>::ErrorCritical(err)) {
            throw pybind11::value_error("Error during loading " + filename + ": '" + tri::io::Importer<CMesh>::ErrorMsg(err) + "'");
        } else{
            std::string warning_str = "Noncritical error (" + std::to_string(err) + ") during loading " +
                    filename + ": '" + tri::io::Importer<CMesh>::ErrorMsg(err) + "'";
            runtime_warning(warning_str);
        }
    }

    std::unordered_map<std::string, pybind11::object> ret;
    if (dtype.equal(npe::dtype("float64"))) {
        load_mesh_vcg<double>(m, mask, ret);
    } else if (dtype.equal(npe::dtype("float32"))) {
        load_mesh_vcg<float>(m, mask, ret);
    } else{
        throw pybind11::value_error("Invalid dtype. Must be one of float32 or float64");
    }

    return ret;


}
npe_end_code()



const char* ds_save_mesh = R"igl_Qu8mg5v7(
Save a mesh
)igl_Qu8mg5v7";
npe_function(save_mesh_internal)
npe_doc(ds_save_mesh)
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

npe_arg(textures, std::vector<std::string>)
npe_arg(normal_maps, std::vector<std::string>)

npe_arg(dtype_f, npe::dtype)
npe_arg(dtype_i, npe::dtype)

npe_begin_code()
{
    write_mesh_vcg<double, int>(filename,
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
                   textures,
                   normal_maps,
                   dtype_f,
                   dtype_i);
}
npe_end_code()