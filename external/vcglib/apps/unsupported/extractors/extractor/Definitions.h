#ifndef ___DEFINITIONS
#define ___DEFINITIONS

#include <vcg/simplex/vertex/with/afvmvn.h>
#include <vcg/simplex/face/with/afavfn.h>
#include <vcg/simplex/edge/with/ae.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/allocate.h>

typedef float ScalarType;

class Edge;
class Face;
class Vertex : public vcg::VertexAFVMVN< ScalarType, Edge, Face > {};
class Face		: public vcg::FaceAFAVFN< Vertex, Edge, Face> {};
class Mesh		: public vcg::tri::TriMesh< std::vector< Vertex>, std::vector< Face > > {};

typedef vcg::tri::Allocator< Mesh > Allocator;
typedef vcg::Box3< int >						BoundingBox;
typedef Vertex* VertexPointer;

#endif //___DEFINITIONS