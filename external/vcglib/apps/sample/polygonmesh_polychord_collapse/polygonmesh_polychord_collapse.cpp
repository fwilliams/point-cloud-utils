#include <iostream>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/complex/algorithms/polygon_support.h>
#include <vcg/complex/algorithms/polygon_polychord_collapse.h>

using namespace vcg;

class PolyVertex;
class PolyFace;

struct PolyUsedTypes : public vcg::UsedTypes<
        vcg::Use<PolyVertex>::AsVertexType,
        vcg::Use<PolyFace>::AsFaceType
        > {};

class PolyVertex : public vcg::Vertex<
        PolyUsedTypes,
        vcg::vertex::Coord3f,
        vcg::vertex::Normal3f,
        vcg::vertex::BitFlags
        > {};

class PolyFace : public vcg::Face<
        PolyUsedTypes,
        vcg::face::PolyInfo,
        vcg::face::Normal3f,
        vcg::face::BitFlags,
        vcg::face::PFVAdj,
        vcg::face::PFFAdj
        > {};

class PolyMesh : public vcg::tri::TriMesh<
        std::vector<PolyVertex>,
        std::vector<PolyFace>
        > {

public:
    /**
     * @brief open Loads a polygonal mesh from file.
     * @param mesh The mesh object into which to extract the mesh.
     * @param filename The filename where the mesh is stored.
     * @return An error code.
     */
    static int openMesh (PolyMesh &mesh, const char *filename)
    {
        // try to load the mesh from the file
        int err = vcg::tri::io::Importer<PolyMesh>::Open(mesh, filename);

        // check if successfully loaded
        if (err == 0)
        {
            // update bounding box
            vcg::tri::UpdateBounding<PolyMesh>::Box(mesh);
            // update topology
            vcg::tri::UpdateTopology<PolyMesh>::FaceFace(mesh);
            // update normals TODO: compute average normal in the polygon
            vcg::tri::UpdateNormal<PolyMesh>::PerFaceNormalized(mesh);
            // update flags
            vcg::tri::UpdateFlags<PolyMesh>::Clear(mesh);
        }

        return err;
    }

    /**
     * @brief saveMesh Writes a polygonal mesh into a file.
     * @param mesh The mesh to write.
     * @param filename The filename where to store the mesh.
     * @return
     */
    static int saveMesh (PolyMesh &mesh, const char *filename) {
        // try to write the mesh into the file
        return vcg::tri::io::Exporter<PolyMesh>::Save(mesh, filename);
    }
};

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cout << "Error. Usage: " << argv[0] << " meshfilename" << std::endl;
    return -1;
  }

  PolyMesh mesh;

  // open mesh
  int err = PolyMesh::openMesh(mesh, argv[1]);
  if (err != 0)
    return err;

  // collapse **********************************************************************************************
  vcg::tri::PolychordCollapse<PolyMesh>::CollapseAllPolychords(mesh, true);

  // these don't work with polygonal meshes:
  vcg::tri::Allocator<PolyMesh>::CompactFaceVector(mesh);
  vcg::tri::Allocator<PolyMesh>::CompactVertexVector(mesh);

  // save mesh
  PolyMesh::saveMesh(mesh, "output.obj");

  return 0;
}

