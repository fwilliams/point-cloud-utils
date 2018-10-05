#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/attribute_seam.h>

#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

/*
  this sample shows how to transfer per wedge attributes from wedges to vertices.
  during the process new vertices could be created.
*/
using namespace vcg;

#define TEST_IN_PLACE_SPLIT

#ifdef TEST_IN_PLACE_SPLIT

class SrcVertex;
class SrcFace;

struct ScrUsedTypes : public UsedTypes<	Use<SrcVertex>::AsVertexType,
																				Use<SrcFace>::AsFaceType>{};

class SrcVertex : public vcg::Vertex
<	ScrUsedTypes,
	vcg::vertex::InfoOcf,
	vcg::vertex::Coord3f,
	vcg::vertex::TexCoordfOcf,
	vcg::vertex::BitFlags
> { };

class SrcFace : public vcg::Face
<	ScrUsedTypes,
	vcg::face::InfoOcf,
	vcg::face::VertexRef,
	vcg::face::WedgeTexCoordfOcf
> { };

class SrcMesh : public vcg::tri::TriMesh  <vcg::vertex::vector_ocf<SrcVertex>, vcg::face::vector_ocf<SrcFace> > { };

typedef SrcVertex DstVertex;
typedef SrcFace   DstFace;
typedef SrcMesh   DstMesh;

#else

// source mesh type: per-wedge texture coordinates
class SrcVertex;
class SrcFace;


struct SrcUsedTypes : public UsedTypes<	Use<SrcVertex>::AsVertexType,
																				Use<SrcFace>::AsFaceType>{};

class SrcVertex : public vcg::Vertex   <SrcUsedTypes, vcg::vertex::Coord3f, vcg::vertex::TexCoord2f, vcg::vertex::BitFlags> { };
class SrcFace   : public vcg::Face     <SrcUsedTypes, vcg::face::VertexRef, vcg::face::WedgeTexCoord2f> { };
class SrcMesh   : public vcg::tri::TriMesh  <std::vector<SrcVertex>, std::vector<SrcFace> > { };


// destination mesh type: per-vertex texture coordinates
class DstVertex; 
class DstFace;

struct DstUsedTypes : public UsedTypes<	Use<SrcVertex>::AsVertexType,
																				Use<SrcFace>::AsFaceType>{};

class DstVertex : public vcg::Vertex   <DstUsedTypes, vcg::vertex::Coord3f, vcg::vertex::TexCoord2f, vcg::vertex::BitFlags> { };
class DstFace   : public vcg::Face     <DstUsedTypes, vcg::face::VertexRef> { };
class DstMesh   : public vcg::tri::TriMesh  <std::vector<DstVertex>, std::vector<DstFace> > { };

#endif

// extract wedge attributes functor.
// given a source face and a wedge index, this functor extracts all the relevant attributes from the wedge
// and transfer them to the destination vertex.
// source and destination meshes are provided to allow for attribute presence checking (.Is*Enabled()).
inline void ExtractVertex(const SrcMesh & srcMesh, const SrcFace & f, int whichWedge, const DstMesh & dstMesh, DstVertex & v)
{
	(void)srcMesh;
	(void)dstMesh;

	v.P() = f.cP(whichWedge);
	v.T() = f.cWT(whichWedge);
}

// sample compare functor.
// given two destination vertices, this functor tells if they are identical in all relevan attributes.
// source and destination meshes are provided to allow for attribute presence checking (.Is*Enabled()).
inline bool CompareVertex(const DstMesh & m, const DstVertex & vA, const DstVertex & vB)
{
	(void)m;

	return (vA.cT() == vB.cT());
}

// sample copy functor.
// given two destination vertices, this functor is asked to copy all relevan attributes.
// source and destination meshes are provided to allow for attribute presence checking (.Is*Enabled()).
inline void CopyVertex(const DstMesh & m, const DstVertex & vSrc, DstVertex & vDst)
{
	(void)m;

	vDst.P() = vSrc.cP();
	vDst.T() = vSrc.cT();
}

void usage(void)
{
	printf("usage : trimesh_split_vertex <src_ply_file_name> <dst_ply_file_name>\n");
	printf("where : <src_ply_file_name> : source PLY trimesh file name with texture coordinates per wedge\n");
	printf("        <dst_ply_file_name> : destination PLY trimesh file name with texture coordinates per vertex\n");
	printf("exit.\n");
}

int main(int argc, char ** argv)
{
	if (argc != 3)
	{
		usage();
		return -1;
	}

	SrcMesh srcMesh;
#ifdef TEST_IN_PLACE_SPLIT
	srcMesh.face.EnableWedgeTexCoord();
#endif
	vcg::tri::io::ImporterPLY<SrcMesh>::Open(srcMesh, argv[1]);
	if ((srcMesh.VN() <= 0) || (srcMesh.FN() <= 0))
	{
		printf("invalid source mesh file.\n");
		return -1;
	}
	const int srcVN = srcMesh.VN();
	const int srcFN = srcMesh.FN();
	printf("source mesh succesfully loaded.\n");

#ifdef TEST_IN_PLACE_SPLIT
	DstMesh & dstMesh = srcMesh;
	dstMesh.vert.EnableTexCoord();
	vcg::tri::AttributeSeam::SplitVertex(dstMesh, ExtractVertex, CompareVertex);
#else
	DstMesh dstMesh;
	vcg::tri::AttributeSeam::SplitVertex(srcMesh, dstMesh, ExtractVertex, CompareVertex, CopyVertex);
	dstMesh.textures = srcMesh.textures;
#endif
	if (vcg::tri::io::ExporterPLY<DstMesh>::Save(dstMesh, argv[2], vcg::tri::io::Mask::IOM_VERTCOORD | vcg::tri::io::Mask::IOM_VERTTEXCOORD) != 0)
	{
		printf("cannot save destination mesh file.\n");
		return -1;
	}
	printf("destination mesh succesfully saved.\n");
	const int dstVN = dstMesh.VN();
	const int dstFN = dstMesh.FN();

	printf("\n");
	printf("statistics:\n");
	printf("  input mesh vertices count    : %d\n", srcVN);
	printf("  input mesh faces count       : %d\n", srcFN);
	printf("  splitted mesh vertices count : %d\n", dstVN);
	printf("  splitted mesh faces count    : %d\n", dstFN);
	printf("\n");

	return 0;
}
