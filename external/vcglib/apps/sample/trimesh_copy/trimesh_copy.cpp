
// stuff to define the mesh
#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>
// io
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

#include <cstdlib>

#include <sys/timeb.h>
#include <iostream>
#include <string>


class MyVertex;
class MyEdge;
class MyFace;

struct MyUsedTypes: public vcg::UsedTypes<vcg::Use<MyVertex>::AsVertexType,vcg::Use<MyEdge>::AsEdgeType,vcg::Use<MyFace>::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes,vcg::vertex::VFAdj,vcg::vertex::Coord3f,vcg::vertex::Normal3f,vcg::vertex::Mark,vcg::vertex::BitFlags  >
{
};

class MyEdge : public vcg::Edge< MyUsedTypes> {};

class MyFace    : public vcg::Face< MyUsedTypes,
    vcg::face::VFAdj,
    vcg::face::VertexRef,
    vcg::face::BitFlags > {};

// the main mesh class
class MyMesh    : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> > {};

class OcfVertex;
class OcfEdge;
class OcfFace;

// Declaration of the semantic of the used types
class OcfUsedTypes: public vcg::UsedTypes < vcg::Use<OcfVertex>::AsVertexType,
    vcg::Use<OcfEdge   >::AsEdgeType,
    vcg::Use<OcfFace  >::AsFaceType >{};


// The Main Vertex Class
// Most of the attributes are optional and must be enabled before use.
// Each vertex needs 40 byte, on 32bit arch. and 44 byte on 64bit arch.

class OcfVertex  : public vcg::Vertex< OcfUsedTypes,vcg::vertex::InfoOcf,vcg::vertex::Coord3f,vcg::vertex::BitFlags,vcg::vertex::Normal3fOcf,vcg::vertex::VFAdjOcf,vcg::vertex::MarkOcf>
{
};


// The Main Edge Class
// Currently it does not contains anything.
class OcfEdge : public vcg::Edge<OcfUsedTypes>
{
};

// Each face needs 32 byte, on 32bit arch. and 48 byte on 64bit arch.
class OcfFace    : public vcg::Face<  OcfUsedTypes,vcg::face::InfoOcf,vcg::face::VertexRef,vcg::face::BitFlags,vcg::face::VFAdjOcf> {};

class OcfMesh    : public vcg::tri::TriMesh< vcg::vertex::vector_ocf<OcfVertex>, vcg::face::vector_ocf<OcfFace> >
{
};

void Usage()
{
    printf(
        "\nUsage:  "\
        "trimeshcopy fileIn -(n|o) [fileOut]\n"\
        "trimeshcopy test vcg::MeshCopy efficiency.\nIt imports a fileIn file into a user defined mesh and test how long vcg::MeshCopy needs to copy the imported mesh in a second one.The copy time is expressed in milliseconds.\nIf the -n flag is used a non-optional attributes mesh will be tested, defining -o, instead, the target mesh will be an ocf one.\nA fileOut file can be passed to the tool in order to check if the mesh was successfully copied.\nThe file will be exported in PLY file format.\n"
        );
    exit(-1);
}

template <class MeshType>
bool UnitTest_Append(const char *filename1, const char *filename2)
{
  MeshType mr;
  MeshType ml;

  int startOpen=clock();
  int err=vcg::tri::io::Importer<MeshType>::Open(mr,filename1);
  if(err)
  {
      std::cerr << "Unable to open mesh " << filename1 << " : " << vcg::tri::io::Importer<MyMesh>::ErrorMsg(err) << std::endl;
      exit(-1);
  }
  int endOpen = clock();
  std::cout << "mesh loaded in " << float(endOpen-startOpen)/CLOCKS_PER_SEC << " msecs. Verts: " << mr.VN() << " Faces: " << mr.FN() << "\n";

  int startCopy = clock();
  vcg::tri::Append<MeshType,MeshType>::Mesh(ml,mr,false,true);
  int endCopy = clock();
  std::cout << "mesh copied in " << float(endCopy-startCopy)/CLOCKS_PER_SEC << " msecs." << std::endl;

  assert(ml.VN()==mr.VN());
  assert(ml.en==mr.en);
  assert(ml.FN()==mr.FN());

  int startSave = clock();
  vcg::tri::io::ExporterPLY<MeshType>::Save(ml,filename2);
  int endSave = clock();
  std::cout << "mesh saved in " << float(endSave-startSave)/CLOCKS_PER_SEC << " msecs." << std::endl;
  return true;
}

int main(int /*argc*/ ,char**argv)
{
    UnitTest_Append<MyMesh>(argv[1],"out.ply");
    UnitTest_Append<OcfMesh>(argv[1],"out.ply");
    return 0;
}
