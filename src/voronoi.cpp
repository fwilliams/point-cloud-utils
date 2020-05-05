#include <Eigen/Core>

#include <vector>
#include <mutex>
#include <sstream>
#include <chrono>

#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_fill_holes.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_degree3_vertices.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/voronoi/RVD_callback.h>
#include <geogram/voronoi/RVD_mesh_builder.h>
#include <geogram/numerics/predicates.h>

#include <npe.h>
#include <pybind11/stl.h>

#include "common.h"
#include "geogram_utils.h"


///////////////////////////////////////
//
// When you come back to this, look here for example code on how to extrac the RVD:
// external/geogram/src/examples/geogram/compute_RVD/main.cpp
//
///////////////////////////////////////

namespace {
/**
     * \brief The callback called for each RVD polyhedron. Constructs a
     *  mesh with the boundary of all cells.
     * \details Its member functions are called for each RVD polyhedron,
     *  i.e. the intersections between the volumetric mesh tetrahedra and
     *  the Voronoi cells. Based on set_simplify_xxx(), a smaller number of
     *  polyhedra can be generated.
     */
class SaveRVDCells : public GEO::RVDPolyhedronCallback {
public:

    /**
     * \brief SaveRVDCells constructor.
     * \param[out] output_mesh a reference to the generated mesh
     */
    SaveRVDCells(GEO::Mesh& output_mesh, bool simplify_tets, bool simplify_voronoi, bool simplify_boundary, double shrink_cells=0.0)
        : output_mesh_(output_mesh), shrink_(shrink_cells) {
        my_vertex_map_ = nullptr;

        // If set, then only one polyhedron per (connected component of) restricted Voronoi
        // cell is generated.
        set_simplify_internal_tet_facets(simplify_tets);

        // If set, then only one polygon per Voronoi facet is generated.
        set_simplify_voronoi_facets(simplify_voronoi);

        // If set, then the intersection between a Voronoi cell and the boundary surface is
        // replaced with a single polygon whenever possible (i.e. when its topology is a
        // disk and when it has at least 3 corners).
        set_simplify_boundary_facets(simplify_boundary);

        // If set, then the intersections are available as Mesh objects through the function
        // process_polyhedron_mesh(). Note that this is implied by simplify_voronoi_facets
        // or simplify_boundary.
        if(shrink_cells != 0.0) {
            set_use_mesh(true);
        }
    }

    ~SaveRVDCells() {
        delete my_vertex_map_;
        my_vertex_map_ = nullptr;
    }

    /**
     * \brief Called at the beginning of RVD traversal.
     */
    virtual void begin() {
        GEO::RVDPolyhedronCallback::begin();
        output_mesh_.clear();
        output_mesh_.vertices.set_dimension(3);
    }

    /**
     * \brief Called at the end of RVD traversal.
     */
    virtual void end() {
        GEO::RVDPolyhedronCallback::end();
        output_mesh_.facets.connect();
    }

    /**
     * \brief Called at the beginning of each RVD polyhedron.
     * \param[in] seed , tetrahedron the (seed,tetrahedron) pair that
     *  defines the RVD polyhedron, as the intersection between the Voronoi
     *  cell of the seed and the tetrahedron.
     */
    virtual void begin_polyhedron(GEO::index_t seed, GEO::index_t tetrahedron) {
        GEO::geo_argused(tetrahedron);
        GEO::geo_argused(seed);

        //   The RVDVertexMap is used to map the symbolic representation of vertices
        // to indices. Here we reset indexing for each new cell, so that vertices shared
        // by the faces of two different cells will be duplicated. We do that because
        // we construct the boundary of the cells in a surfacic mesh (for visualization
        // purposes). Client code that has a data structure for polyhedral volumetric mesh
        // will not want to reset indexing (and will comment-out the following three lines).
        // It will also construct the RVDVertexMap in the constructor.

        delete my_vertex_map_;
        my_vertex_map_ = new GEO::RVDVertexMap;
        my_vertex_map_->set_first_vertex_index(output_mesh_.vertices.nb());
    }

    /**
     * \brief Called at the beginning of each RVD polyhedron.
     * \param[in] facet_seed if the facet is on a Voronoi bisector,
     *  the index of the Voronoi seed on the other side of the bisector,
     *  else index_t(-1)
     * \param[in] facet_tet if the facet is on a tethedral facet, then
     *  the index of the tetrahedron on the other side, else index_t(-1)
     */
    virtual void begin_facet(GEO::index_t facet_seed, GEO::index_t facet_tet) {
        GEO::geo_argused(facet_seed);
        GEO::geo_argused(facet_tet);
        current_facet_.resize(0);
    }

    virtual void vertex(const double* geometry, const GEOGen::SymbolicVertex& symb) {
        // Find the index of the vertex associated with its symbolic representation.
        GEO::index_t vid = my_vertex_map_->find_or_create_vertex(seed(), symb);

        // If the vertex does not exist in the mesh, create it.
        if(vid >= output_mesh_.vertices.nb()) {
            output_mesh_.vertices.create_vertex(geometry);
        }

        // Memorize the current facet.
        current_facet_.push_back(vid);
    }

    virtual void end_facet() {
        // Create the facet from the memorized indices.
        GEO::index_t f = output_mesh_.facets.nb();
        output_mesh_.facets.create_polygon(current_facet_.size());
        for(GEO::index_t i=0; i<current_facet_.size(); ++i) {
            output_mesh_.facets.set_vertex(f,i,current_facet_[i]);
        }
    }

    virtual void end_polyhedron() {
        // Nothing to do.
    }

    virtual void process_polyhedron_mesh() {
        // This function is called for each cell if set_use_mesh(true) was called.
        // It is the case if simplify_voronoi_facets(true) or
        // simplify_boundary_facets(true) was called.
        //   Note1: most users will not need to overload this function (advanded use
        //   only).
        //   Note2: mesh_ is managed internally by RVDPolyhedronCallback class, as an
        // intermediary representation to store the cell before calling the callbacks.
        // It is distinct from the output_mesh_ constructed by the callbacks.

        //   The current cell represented by a Mesh can be
        // filtered/modified/post-processed (member variable mesh_)
        // here, before calling base class's implementation.
        //   As an example, we shrink the cells. More drastic modifications/
        // transformations of the mesh can be done (see base class's implementation
        // in geogram/voronoi/RVD_polyhedron_callback.cpp).

        double shrink = shrink_;
        if(shrink != 0.0 && mesh_.vertices.nb() != 0) {
            GEO::vec3 center(0.0, 0.0, 0.0);
            for(GEO::index_t v=0; v<mesh_.vertices.nb(); ++v) {
                center += GEO::vec3(mesh_.vertices.point_ptr(v));
            }
            center = (1.0 / double(mesh_.vertices.nb())) * center;
            for(GEO::index_t v=0; v<mesh_.vertices.nb(); ++v) {
                GEO::vec3 p(mesh_.vertices.point_ptr(v));
                p = shrink * center + (1.0 - shrink) * p;
                mesh_.vertices.point_ptr(v)[0] = p.x;
                mesh_.vertices.point_ptr(v)[1] = p.y;
                mesh_.vertices.point_ptr(v)[2] = p.z;
            }
        }

        //  The default implementation simplifies Voronoi facets
        // and boundary mesh facets based on the boolean flags
        // defined by set_simplify_xxx(). Then it calls the callbacks
        // for each mesh facet.
        RVDPolyhedronCallback::process_polyhedron_mesh();

    }

private:
    GEO::vector<GEO::index_t> current_facet_;
    GEO::Mesh& output_mesh_;
    double shrink_;
    GEO::RVDVertexMap* my_vertex_map_;
};


/*
 * Compute the restricted Voronoi Diagram
 */
template <typename DerivedV, typename DerivedF, typename DerivedT, typename DerivedP>
void voronoi_diagram(const Eigen::MatrixBase<DerivedV> &V,
                     const Eigen::MatrixBase<DerivedF> &F,
                     const Eigen::MatrixBase<DerivedT> &T,
                     const Eigen::MatrixBase<DerivedP> &centers,
                     Eigen::PlainObjectBase<DerivedP> &retV,
                     std::vector<std::vector<int>> &retF,
                     std::vector<std::vector<int>> &retC) {

    // Why are we copying here?
    std::vector<double> centers_vec;
    for (GEO::index_t i = 0; i < centers.rows(); i++) {
        for (GEO::index_t j = 0; j < centers.cols(); j++) {
            centers_vec.push_back(centers(i, j));
        }
    }
    GEO::Mesh M;
    vft_to_geogram_tet_mesh(V, F, T, M);

    GEO::Delaunay_var delaunay_ = GEO::Delaunay::create(3);
    GEO::RestrictedVoronoiDiagram_var RVD_ = GEO::RestrictedVoronoiDiagram::create(delaunay_, &M);

    GEO::index_t nb_points = centers.rows();
    delaunay_->set_vertices(nb_points, centers_vec.data());

    RVD_->set_volumetric(true);

    GEO::Mesh retM;
    SaveRVDCells callback(retM, false, false, false, 0.0);
    RVD_->for_each_polyhedron(callback);

    retV.resize(retM.vertices.nb(), retM.vertices.dimension());
    for (int i = 0; i < retM.vertices.nb(); i++) {
        for (int j = 0; j < retM.vertices.dimension(); j++) {
            retV(i, j) = retM.vertices.point_ptr(i)[j];
        }
    }

//    retF.resize(retM.facets.nb());
//    for (int i = 0; i < retM.facets.nb(); i++) {
//        for (int j = 0; j < retM.facets.nb_vertices(i); j++) {
//            retF[i].push_back(retM.facets.vertex(i, j));
//        }
//    }

//    retC.resize(retM.cells.nb());
//    for (int i = 0; i < retM.cells.nb(); i++) {
//        for (int j = 0; j < retM.cells.nb_vertices(i); j++) {
//            retC[i].push_back(retM.cells.vertex(i, j));
//        }
//        pybind11::print("Cell " + std::to_string(i));
//        switch(retM.cells.type(0)) {
//        case GEO::MeshCellType::MESH_TET:
//            pybind11::print("  Cell Type: MESH_TET");
//            break;
//        case GEO::MeshCellType::MESH_HEX:
//            pybind11::print("  Cell Type: MESH_HEX");
//            break;
//        case GEO::MeshCellType::MESH_PRISM:
//            pybind11::print("  Cell Type: MESH_PRISM");
//            break;
//        case GEO::MeshCellType::MESH_PYRAMID:
//            pybind11::print("  Cell Type: MESH_PYRAMID");
//            break;
//        case GEO::MeshCellType::MESH_CONNECTOR:
//            pybind11::print("  Cell Type: MESH_CONNECTOR");
//            break;
//        default:
//            pybind11::print("  Cell Type: WTF...");
//        }
//        pybind11::print("  NB Vertices: " + std::to_string(retM.cells.nb_vertices(i)));
//        pybind11::print("  NB Corners: " + std::to_string(retM.cells.nb_corners(i)));
//    }

//    std::stringstream ss_v;
//    ss_v << "Vertex attributes: ";
//    GEO::vector<std::string> v_attrib_names;
//    M.vertices.attributes().list_attribute_names(v_attrib_names);
//    for (int i = 0; i < v_attrib_names.size(); i++) {
//        ss_v << v_attrib_names[i] << " ";
//    }
//    pybind11::print(ss_v.str());

//    std::stringstream ss_e;
//    ss_e << "Edge attributes: ";
//    GEO::vector<std::string> e_attrib_names;
//    M.edges.attributes().list_attribute_names(e_attrib_names);
//    for (int i = 0; i < e_attrib_names.size(); i++) {
//        ss_e << e_attrib_names[i] << " ";
//    }
//    pybind11::print(ss_e.str());

//    std::stringstream ss_f;
//    ss_f << "Facet attributes: ";
//    GEO::vector<std::string> f_attrib_names;
//    M.facets.attributes().list_attribute_names(f_attrib_names);
//    for (int i = 0; i < f_attrib_names.size(); i++) {
//        ss_f << f_attrib_names[i] << " ";
//    }
//    pybind11::print(ss_f.str());

//    std::stringstream ss_c;
//    ss_c << "Cell attributes: ";
//    GEO::vector<std::string> c_attrib_names;
//    M.cells.attributes().list_attribute_names(c_attrib_names);
//    for (int i = 0; i < c_attrib_names.size(); i++) {
//        ss_c << c_attrib_names[i] << " ";
//    }
//    pybind11::print(ss_c.str());

//    pybind11::print("There are " + std::to_string(retM.vertices.nb()) + " vertices, " +
//                    std::to_string(retM.edges.nb()) + " edges, " +
//                    std::to_string(retM.facets.nb()) + " facets, and " +
//                    std::to_string(retM.cells.nb()) + " cells.");

//    switch(retM.cells.type(0)) {
//    case GEO::MeshCellType::MESH_TET:
//        pybind11::print("Cell Type: MESH_TET");
//        break;
//    case GEO::MeshCellType::MESH_HEX:
//        pybind11::print("Cell Type: MESH_HEX");
//        break;
//    case GEO::MeshCellType::MESH_PRISM:
//        pybind11::print("Cell Type: MESH_PRISM");
//        break;
//    case GEO::MeshCellType::MESH_PYRAMID:
//        pybind11::print("Cell Type: MESH_PYRAMID");
//        break;
//    case GEO::MeshCellType::MESH_CONNECTOR:
//        pybind11::print("Cell Type: MESH_CONNECTOR");
//        break;
//    default:
//        pybind11::print("Cell Type: WTF...");
//    }
}

}

const char* voronoi_diagram_3d_doc = R"Qu8mg5v7(
Compute the 3D voronoi diagram

Parameters
----------
centers: A [n, 3] array where each row is a point in a voronoi diagram restricted to
         the unit cube [0, 1]^3

Returns
-------
v, f, c where:
        v is a nv-by-3 array of voronoi vertices
        f is a list of nf lists of voronoi faces (indexing into v)
        c is a list of nc lists of voronoi cells (indexing into f)

)Qu8mg5v7";
npe_function(voronoi_diagram_3d)
npe_arg(centers, dense_float, dense_double)
npe_doc(voronoi_diagram_3d_doc)
npe_begin_code()

  if (centers.maxCoeff() > 1.0 || centers.minCoeff() < 0.0) {
    throw pybind11::value_error("Centers of voronoi diagram must lie in the unit cube [0, 1]^3");
  }

  init_geogram_only_once();

  Eigen::MatrixXd V(9, 3);
  V << 0. , 1. , 0. ,
       0. , 0. , 1. ,
       1. , 0. , 1. ,
       1. , 0. , 0. ,
       0. , 0. , 0. ,
       1. , 1. , 1. ,
       1. , 1. , 0. ,
       0. , 1. , 1. ,
       0.5, 0.5, 0.5;

  Eigen::MatrixXi T(12, 4);
  T << 8, 0, 7, 5,
       6, 2, 5, 8,
       3, 1, 2, 8,
       0, 8, 7, 4,
       7, 8, 5, 1,
       8, 0, 5, 6,
       3, 1, 8, 4,
       2, 1, 5, 8,
       8, 0, 6, 4,
       7, 8, 1, 4,
       6, 4, 3, 8,
       2, 6, 3, 8;
  Eigen::MatrixXi F(12, 3);
  F << 0, 5, 7,
       6, 2, 5,
       3, 1, 2,
       0, 7, 4,
       7, 5, 1,
       0, 6, 5,
       3, 4, 1,
       2, 1, 5,
       0, 4, 6,
       7, 1, 4,
       6, 4, 3,
       2, 6, 3;

//  auto start = std::chrono::high_resolution_clock::now();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ctrs_copy = centers.template cast<double>();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> voronoi_v;
  std::vector<std::vector<int>> voronoi_f;
  std::vector<std::vector<int>> voronoi_c;

  GEO::Mesh input_mesh;
  vft_to_geogram_tet_mesh(V, F, T, input_mesh);

  GEO::Delaunay_var delaunay_ = GEO::Delaunay::create(3);
  GEO::RestrictedVoronoiDiagram_var RVD_ = GEO::RestrictedVoronoiDiagram::create(delaunay_, &input_mesh);

  GEO::index_t nb_points = centers.rows();
  delaunay_->set_vertices(nb_points, ctrs_copy.data());

  RVD_->set_volumetric(true);

//  auto init_end = std::chrono::high_resolution_clock::now();
  GEO::Mesh output_mesh;
  SaveRVDCells callback(output_mesh, false, false, false, 0.0);
  RVD_->for_each_polyhedron(callback);
//  auto rvd_end = std::chrono::high_resolution_clock::now();

  voronoi_v.resize(output_mesh.vertices.nb(), output_mesh.vertices.dimension());
  for (int i = 0; i < output_mesh.vertices.nb(); i++) {
      for (int j = 0; j < output_mesh.vertices.dimension(); j++) {
          voronoi_v(i, j) = output_mesh.vertices.point_ptr(i)[j];
      }
  }

  voronoi_f.resize(output_mesh.facets.nb());
  for (int i = 0; i < output_mesh.facets.nb(); i++) {
      for (int j = 0; j < output_mesh.facets.nb_vertices(i); j++) {
          voronoi_f[i].push_back(output_mesh.facets.vertex(i, j));
      }
  }

  voronoi_c.resize(output_mesh.cells.nb());
  for (int i = 0; i < output_mesh.cells.nb(); i++) {
      voronoi_c[i].resize(output_mesh.cells.nb_vertices(i));
      for (int j = 0; j < output_mesh.cells.nb_vertices(i); j++) {
          voronoi_c[i].push_back(output_mesh.cells.vertex(i, j));
      }
  }

//  auto write_end = std::chrono::high_resolution_clock::now();

//  std::chrono::duration<double> elapsed_init = init_end - start;
//  std::chrono::duration<double> elapsed_rvd = rvd_end - init_end;
//  std::chrono::duration<double> elapsed_write = write_end - rvd_end;
//  std::chrono::duration<double> elapsed_tot = write_end - start;

//  pybind11::print("Init time: " +std::to_string(elapsed_init.count()) + "s");
//  pybind11::print("RVD time: " +std::to_string(elapsed_rvd.count()) + "s");
//  pybind11::print("Write time: " +std::to_string(elapsed_write.count()) + "s");
//  pybind11::print("Total time: " +std::to_string(elapsed_tot.count()) + "s");

  return std::make_tuple(npe::move(voronoi_v));

npe_end_code()
