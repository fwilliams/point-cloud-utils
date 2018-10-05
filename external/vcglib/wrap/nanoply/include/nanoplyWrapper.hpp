/****************************************************************************
* NanoPLY                                                                   *
* NanoPLY is a C++11 header-only library to read and write PLY file         *
*                                                                           *
* Copyright(C) 2014-2015                                                    *
* Visual Computing Lab                                                      *
* ISTI - Italian National Research Council                                  *
*                                                                           *
* This Source Code Form is subject to the terms of the Mozilla Public       *
* License, v. 2.0. If a copy of the MPL was not distributed with this       *
* file, You can obtain one at http://mozilla.org/MPL/2.0/.                  *
*                                                                           *
****************************************************************************/

#ifndef NANOPLY_WRAPPER_VCG_H
#define NANOPLY_WRAPPER_VCG_H

#include <wrap/nanoply/include/nanoply.hpp>
#include <vcg/space/point.h>
#include <map>

namespace nanoply
{

	template <class MeshType>
	class NanoPlyWrapper{

				
	private:

    typedef typename MeshType::PointerToAttribute PointerToAttribute;
    typedef typename MeshType::ScalarType ScalarType;

    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::VertexType::ScalarType VertexCoordScalar;
    typedef typename MeshType::VertexType::NormalType::ScalarType VertexNormScalar;
    typedef typename MeshType::VertexType::ColorType::ScalarType VertexColorScalar;
    typedef typename MeshType::VertexType::QualityType VertexQuality;
    typedef typename MeshType::VertexType::RadiusType VertexRadius;
    typedef typename MeshType::VertexType::FlagType VertexFlag;
    typedef typename MeshType::VertexType::TexCoordType::ScalarType VertexTexScalar;
    typedef typename MeshType::VertexType::CurvatureType::ScalarType VertexCurScalar;
    typedef typename MeshType::VertexType::CurScalarType VertexDirCurScalar;
    typedef typename MeshType::VertexType::CurVecType::ScalarType VertexDirCurVecScalar;

    typedef typename MeshType::EdgeType EdgeType;
    typedef typename MeshType::EdgeType::ColorType::ScalarType EdgeColorScalar;
    typedef typename MeshType::EdgeType::QualityType EdgeQuality;
    typedef typename MeshType::EdgeType::FlagType EdgeFlag;

    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::FaceType::NormalType::ScalarType FaceNormScalar;
    typedef typename MeshType::FaceType::ColorType::ScalarType FaceColorScalar;
    typedef typename MeshType::FaceType::QualityType FaceQuality;
    typedef typename MeshType::FaceType::FlagType FaceFlag;
    typedef typename MeshType::FaceType::TexCoordType::ScalarType FaceTexScalar;
    typedef typename MeshType::FaceType::CurScalarType FaceDirCurScalar;
    typedef typename MeshType::FaceType::CurVecType::ScalarType FaceDirCurVecScalar;
    typedef typename MeshType::FaceType::WedgeColorType::ScalarType WedgeColorScalar;
    typedef typename MeshType::FaceType::WedgeNormalType::ScalarType WedgeNormalScalar;

       

		template<class T> static PlyType getEntity() { return NNP_UNKNOWN_TYPE };
		template<> static PlyType getEntity<unsigned char>(){ return NNP_UINT8; };
		template<> static PlyType getEntity<char>(){ return NNP_INT8; };
		template<> static PlyType getEntity<unsigned short>(){ return NNP_UINT16; };
		template<> static PlyType getEntity<short>(){ return NNP_INT16; };
		template<> static PlyType getEntity<unsigned int>(){ return NNP_UINT32; };
		template<> static PlyType getEntity<int>(){ return NNP_INT32; };
		template<> static PlyType getEntity<float>(){ return NNP_FLOAT32; };
		template<> static PlyType getEntity<double>(){ return NNP_FLOAT64; };

    template<class T> static PlyType getEntityList() { return NNP_UNKNOWN_TYPE; };
		template<> static PlyType getEntityList<unsigned char>(){ return NNP_LIST_UINT8_UINT8; };
		template<> static PlyType getEntityList<char>(){ return NNP_LIST_UINT8_INT8; };
		template<> static PlyType getEntityList<unsigned short>(){ return NNP_LIST_UINT8_UINT16; };
		template<> static PlyType getEntityList<short>(){ return NNP_LIST_UINT8_INT16; };
		template<> static PlyType getEntityList<unsigned int>(){ return NNP_LIST_UINT8_UINT32; };
		template<> static PlyType getEntityList<int>(){ return NNP_LIST_UINT8_INT32; };
		template<> static PlyType getEntityList<float>(){ return NNP_LIST_UINT8_FLOAT32; };
		template<> static PlyType getEntityList<double>(){ return NNP_LIST_UINT8_FLOAT64; };


		template<class Container, class Type, int n>
		inline static void PushDescriport(std::vector<PlyProperty>& prop, ElementDescriptor& elem, PlyEntity entity, void* ptr)
		{
			prop.push_back(PlyProperty(getEntity<Type>(), entity));
			DescriptorInterface* di = new DataDescriptor<Container, n, Type>(entity, ptr);
			elem.dataDescriptor.push_back(di);
		}


		template<class Container, class Type, int n>
		inline static void PushDescriportList(std::vector<PlyProperty>& prop, ElementDescriptor& elem, PlyEntity entity, void* ptr)
		{
			prop.push_back(PlyProperty(getEntityList<Type>(), entity));
			DescriptorInterface* di = new DataDescriptor<Container, n, Type>(entity, ptr);
			elem.dataDescriptor.push_back(di);
		}


		template<class Container, class Type, int n>
		inline static void PushDescriport(std::vector<PlyProperty>& prop, ElementDescriptor& elem, std::string& name, void* ptr)
		{
			prop.push_back(PlyProperty(getEntity<Type>(), name));
			DescriptorInterface* di = new DataDescriptor<Container, n, Type>(name, ptr);
			elem.dataDescriptor.push_back(di);
		}


		template<class Container, class Type, int n>
		inline static void PushDescriportList(std::vector<PlyProperty>& prop, ElementDescriptor& elem, std::string& name, void* ptr)
		{
			prop.push_back(PlyProperty(getEntityList<Type>(), name));
			DescriptorInterface* di = new DataDescriptor<Container, n, Type>(name, ptr);
			elem.dataDescriptor.push_back(di);
		}
	
	public:


		typedef enum {
			IO_NONE = 0x00000000,

			IO_VERTCOORD = 0x00000001,
			IO_VERTFLAGS = 0x00000002,
			IO_VERTCOLOR = 0x00000004,
			IO_VERTQUALITY = 0x00000008,
			IO_VERTNORMAL = 0x00000010,
			IO_VERTTEXCOORD = 0x00000020,
			IO_VERTRADIUS = 0x00000040,
			IO_VERTCURV = 0x00000080,
			IO_VERTCURVDIR = 0x00000100,
			IO_VERTATTRIB = 0x00000200,

			IO_FACEINDEX = 0x00000400,
			IO_FACEFLAGS = 0x00000800,
			IO_FACECOLOR = 0x00001000,
			IO_FACEQUALITY = 0x00002000,
			IO_FACENORMAL = 0x00004000,
			IO_FACECURVDIR = 0x00008000,
			IO_FACEATTRIB = 0x00010000,

			IO_EDGEINDEX = 0x00020000,
			IO_EDGEQUALITY = 0x00040000,
			IO_EDGECOLOR = 0x00080000,
			IO_EDGEFLAGS = 0x00100000,
			IO_EDGEATTRIB = 0x00200000,

			IO_WEDGCOLOR = 0x00400000,
			IO_WEDGTEXCOORD = 0x00800000,
			IO_WEDGTEXMULTI = 0x01000000, // when textrue index is explicit
			IO_WEDGNORMAL = 0x02000000,

			IO_BITPOLYGONAL = 0x04000000, // loads explicit polygonal mesh

			IO_CAMERA = 0x08000000,
			IO_MESHATTRIB = 0x10000000,

			IO_FLAGS = IO_VERTFLAGS | IO_FACEFLAGS,

			IO_ALL = 0xFFFFFFFF
		}BitMask;



		class CustomAttributeDescriptor
		{
		public:

			typedef std::map<std::string, ElementDescriptor::PropertyDescriptor> MapMeshAttrib;
			typedef std::map<std::string, ElementDescriptor::PropertyDescriptor>::iterator MapMeshAttribIter;
			typedef std::map<std::string, std::vector<PlyProperty>> MapMeshAttribProp;
			typedef std::map<std::string, std::vector<PlyProperty>>::iterator MapMeshAttribPropIter;
			
			ElementDescriptor::PropertyDescriptor vertexAttrib;
			ElementDescriptor::PropertyDescriptor faceAttrib;
			ElementDescriptor::PropertyDescriptor edgeAttrib;
			std::vector<PlyProperty> vertexAttribProp;
			std::vector<PlyProperty> faceAttribProp;
			std::vector<PlyProperty> edgeAttribProp;

			MapMeshAttrib meshAttrib;
			MapMeshAttribProp meshAttribProp;
			std::map<std::string, int> meshAttribCnt;


			CustomAttributeDescriptor::~CustomAttributeDescriptor()
			{
				for (int i = 0; i < vertexAttrib.size(); i++)
					delete vertexAttrib[i];
				for (int i = 0; i < edgeAttrib.size(); i++)
					delete edgeAttrib[i];
				for (int i = 0; i < faceAttrib.size(); i++)
					delete faceAttrib[i];
				CustomAttributeDescriptor::MapMeshAttribIter iter = meshAttrib.begin();
				for (; iter != meshAttrib.end(); iter++)
					for (int i = 0; i < (*iter).second.size(); i++)
						delete (*iter).second[i];
			}


			template<class Container, class Type, int n>
			void AddVertexAttribDescriptor(std::string& name, PlyType type, void* ptr)
			{
				vertexAttrib.push_back(new DataDescriptor<Container, n, Type>(name, ptr));
				vertexAttribProp.push_back(PlyProperty(type, name));
			}

			template<class Container, class Type, int n>
			void AddEdgeAttribDescriptor(std::string& name, PlyType type, void* ptr)
			{
				edgeAttrib.push_back(new DataDescriptor<Container, n, Type>(name, ptr));
				edgeAttribProp.push_back(PlyProperty(type, name));
			}

			template<class Container, class Type, int n>
			void AddFaceAttribDescriptor(std::string& name, PlyType type, void* ptr)
			{
				faceAttrib.push_back(new DataDescriptor<Container, n, Type>(name, ptr));
				faceAttribProp.push_back(PlyProperty(type, name));
			}

			template<class Container, class Type, int n>
			void AddMeshAttribDescriptor(std::string& nameAttrib, std::string& nameProp, PlyType type, void* ptr)
			{
				meshAttrib[nameAttrib].push_back(new DataDescriptor<Container, n, Type>(nameProp, ptr));
				meshAttribProp[nameAttrib].push_back(PlyProperty(type, nameProp));
			}

			void AddMeshAttrib(std::string& name, int cnt)
			{
				meshAttribCnt[name] = cnt;
			}

			void GetMeshAttrib(std::string filename)
			{
				nanoply::Info info(filename);
				if (info.errInfo == nanoply::NNP_OK)
				{
					for (int i = 0; i < info.elemVec.size(); i++)
					{
						if (info.elemVec[i].plyElem == NNP_UNKNOWN_ELEM && info.elemVec[i].name != "camera")
							meshAttribCnt[info.elemVec[i].name] = info.elemVec[i].cnt;
					}
				}
			}

		};
		


		static int LoadModel(const char* filename, MeshType& mesh, unsigned int bitMask, CustomAttributeDescriptor& custom)
		{
			nanoply::Info info(filename);
			if (info.errInfo != nanoply::NNP_OK)
				return info.errInfo;

			//Camera
			ElementDescriptor cameraDescr(std::string("camera"));
			vcg::Point3<ScalarType> tra;
			vcg::Matrix44<ScalarType> rot;
			size_t count = info.GetElementCount(std::string("camera"));
			if (count > 0 && (bitMask & BitMask::IO_CAMERA))
			{
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("view_px"), &tra[0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("view_py"), &tra[1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("view_pz"), &tra[2]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("x_axisx"), &rot[0][0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("x_axisy"), &rot[0][1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("x_axisz"), &rot[0][2]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("y_axisx"), &rot[1][0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("y_axisy"), &rot[1][1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("y_axisz"), &rot[1][2]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("z_axisx"), &rot[2][0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("z_axisy"), &rot[2][1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("z_axisz"), &rot[2][2]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("focal"), &mesh.shot.Intrinsics.FocalMm));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("scalex"), &mesh.shot.Intrinsics.PixelSizeMm[0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("scaley"), &mesh.shot.Intrinsics.PixelSizeMm[1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("centerx"), &mesh.shot.Intrinsics.CenterPx[0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("centery"), &mesh.shot.Intrinsics.CenterPx[1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<int, 1, int>(std::string("viewportx"), &mesh.shot.Intrinsics.ViewportPx[0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<int, 1, int>(std::string("viewporty"), &mesh.shot.Intrinsics.ViewportPx[1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("k1"), &mesh.shot.Intrinsics.k[0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("k2"), &mesh.shot.Intrinsics.k[1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("k3"), &mesh.shot.Intrinsics.k[2]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("k4"), &mesh.shot.Intrinsics.k[3]));
			}
			
			//Vertex
			std::vector<std::string> nameList;
			VertexType::Name(nameList);
			ElementDescriptor vertexDescr(NNP_VERTEX_ELEM);
			count = info.GetVertexCount();
			if (nameList.size() > 0 && count > 0)
			{
				vcg::tri::Allocator<MeshType>::AddVertices(mesh, count);
				if ((bitMask & BitMask::IO_VERTCOORD) && VertexType::HasCoord())
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 3, VertexCoordScalar>(NNP_PXYZ, (*mesh.vert.begin()).P().V()));
				if ((bitMask & BitMask::IO_VERTNORMAL) && vcg::tri::HasPerVertexNormal(mesh))
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 3, VertexNormScalar>(NNP_NXYZ, (*mesh.vert.begin()).N().V()));
				if ((bitMask & BitMask::IO_VERTCOLOR) && vcg::tri::HasPerVertexColor(mesh))
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 4, VertexColorScalar>(NNP_CRGBA, (*mesh.vert.begin()).C().V()));
				if ((bitMask & BitMask::IO_VERTQUALITY) && vcg::tri::HasPerVertexQuality(mesh))
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexQuality>(NNP_QUALITY, &(*mesh.vert.begin()).Q()));
				if ((bitMask & BitMask::IO_VERTFLAGS) && vcg::tri::HasPerVertexFlags(mesh))
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexFlag>(NNP_BITFLAG, &(*mesh.vert.begin()).Flags()));
				if ((bitMask & BitMask::IO_VERTRADIUS) && vcg::tri::HasPerVertexRadius(mesh))
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexRadius>(NNP_DENSITY, &(*mesh.vert.begin()).R()));
				if ((bitMask & BitMask::IO_VERTTEXCOORD) && vcg::tri::HasPerVertexTexCoord(mesh))
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 2, VertexTexScalar>(NNP_TEXTURE2D, (*mesh.vert.begin()).T().P().V()));
				if ((bitMask & BitMask::IO_VERTCURV) && vcg::tri::HasPerVertexCurvature(mesh))
				{
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexCurScalar>(NNP_KG, &(*mesh.vert.begin()).Kg()));
          vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexCurScalar>(NNP_KH, &(*mesh.vert.begin()).Kh()));
				}
				if ((bitMask & BitMask::IO_VERTCURVDIR) && vcg::tri::HasPerVertexCurvatureDir(mesh))
				{
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexDirCurScalar>(NNP_K1, &(*mesh.vert.begin()).K1()));
          vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexDirCurScalar>(NNP_K2, &(*mesh.vert.begin()).K2()));
          vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 3, VertexDirCurVecScalar>(NNP_K1DIR, (*mesh.vert.begin()).PD1().V()));
          vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 3, VertexDirCurVecScalar>(NNP_K2DIR, (*mesh.vert.begin()).PD2().V()));
				}
				if ((bitMask & BitMask::IO_VERTATTRIB) && custom.vertexAttrib.size() > 0)
				{
					for (int i = 0; i < custom.vertexAttrib.size(); i++)
					{
            std::set<PointerToAttribute>::iterator ai;
						for (ai = mesh.vert_attr.begin(); ai != mesh.vert_attr.end(); ++ai)
						{
							if ((*custom.vertexAttrib[i]).name == (*ai)._name)
							{
								custom.vertexAttrib[i]->base = ai->_handle->DataBegin();
								break;
							}
						}
						vertexDescr.dataDescriptor.push_back(custom.vertexAttrib[i]);
					}
				}
			}

			//Edge
			nameList.clear();
			EdgeType::Name(nameList);
			ElementDescriptor edgeDescr(NNP_EDGE_ELEM);
			count = info.GetEdgeCount();
			std::vector<vcg::Point2i> edgeIndex;
			if (nameList.size() > 0 && count > 0)
			{
				vcg::tri::Allocator<MeshType>::AddEdges(mesh, count);
				if ((bitMask & BitMask::IO_EDGEINDEX) && MeshType::EdgeType::HasVertexRef())
				{
					edgeIndex.resize(count);
					edgeDescr.dataDescriptor.push_back(new DataDescriptor<vcg::Point2i, 1, int>(NNP_EDGE_V1, &(*edgeIndex.begin()).V()[0]));
					edgeDescr.dataDescriptor.push_back(new DataDescriptor<vcg::Point2i, 1, int>(NNP_EDGE_V2, &(*edgeIndex.begin()).V()[1]));
				}
				if ((bitMask & BitMask::IO_EDGEQUALITY) && vcg::tri::HasPerEdgeQuality(mesh))
					edgeDescr.dataDescriptor.push_back(new DataDescriptor<EdgeType, 1, EdgeQuality>(NNP_QUALITY, &(*mesh.edge.begin()).Q()));
				if ((bitMask & BitMask::IO_EDGECOLOR) && vcg::tri::HasPerEdgeColor(mesh))
					edgeDescr.dataDescriptor.push_back(new DataDescriptor<EdgeType, 4, EdgeColorScalar>(NNP_CRGBA, (*mesh.edge.begin()).C().V()));
				if ((bitMask & BitMask::IO_EDGEFLAGS) && vcg::tri::HasPerEdgeFlags(mesh))
					edgeDescr.dataDescriptor.push_back(new DataDescriptor<EdgeType, 1, EdgeFlag>(NNP_BITFLAG, &(*mesh.edge.begin()).Flags()));
				if ((bitMask & BitMask::IO_EDGEATTRIB) && custom.edgeAttrib.size() > 0)
				{
          for (int i = 0; i < custom.edgeAttrib.size(); i++)
          {
            std::set<PointerToAttribute>::iterator ai;
            for (ai = mesh.edge_attr.begin(); ai != mesh.edge_attr.end(); ++ai)
            {
              if ((*custom.edgeAttrib[i]).name == (*ai)._name)
              {
                custom.edgeAttrib[i]->base = ai->_handle->DataBegin();
                break;
              }
            }
            edgeDescr.dataDescriptor.push_back(custom.edgeAttrib[i]);
          }
				}
			}

			//Face
			nameList.clear();
			FaceType::Name(nameList);
			ElementDescriptor faceDescr(NNP_FACE_ELEM);
			count = info.GetFaceCount();
			std::vector<vcg::Point3i> faceIndex;
			std::vector<vcg::ndim::Point<6, FaceTexScalar>> wedgeTexCoord;
			if (nameList.size() > 0 && count > 0)
			{
				vcg::tri::Allocator<MeshType>::AddFaces(mesh, count);
				if ((bitMask & BitMask::IO_FACEINDEX) && FaceType::HasVertexRef())
				{
					faceIndex.resize(count);
					faceDescr.dataDescriptor.push_back(new DataDescriptor<vcg::Point3i, 3, int>(NNP_FACE_VERTEX_LIST, (*faceIndex.begin()).V()));
				}
				if ((bitMask & BitMask::IO_FACEFLAGS) && vcg::tri::HasPerFaceFlags(mesh))
					faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 1, FaceFlag>(NNP_BITFLAG, &(*mesh.face.begin()).Flags()));
				if ((bitMask & BitMask::IO_FACECOLOR) && vcg::tri::HasPerFaceColor(mesh))
					faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 4, FaceColorScalar>(NNP_CRGBA, (*mesh.face.begin()).C().V()));
				if ((bitMask & BitMask::IO_FACEQUALITY) && vcg::tri::HasPerFaceQuality(mesh))
					faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 1, FaceQuality>(NNP_QUALITY, &(*mesh.face.begin()).Q()));
				if ((bitMask & BitMask::IO_FACENORMAL) && vcg::tri::HasPerFaceNormal(mesh))
					faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 3, FaceNormScalar>(NNP_NXYZ, (*mesh.face.begin()).N().V()));
				if ((bitMask & BitMask::IO_VERTCURVDIR) && vcg::tri::HasPerFaceCurvatureDir(mesh))
				{
					faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 1, FaceDirCurScalar>(NNP_K1, &(*mesh.face.begin()).K1()));
          faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 1, FaceDirCurScalar>(NNP_K2, &(*mesh.face.begin()).K2()));
					faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 3, FaceDirCurVecScalar>(NNP_K1DIR, (*mesh.face.begin()).PD1().V()));
          faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 3, FaceDirCurVecScalar>(NNP_K2DIR, (*mesh.face.begin()).PD2().V()));
				}
				if (((bitMask & BitMask::IO_WEDGTEXCOORD) || (bitMask & BitMask::IO_WEDGTEXMULTI)) && vcg::tri::HasPerWedgeTexCoord(mesh))
				{
					wedgeTexCoord.resize(count);
          faceDescr.dataDescriptor.push_back(new DataDescriptor<vcg::ndim::Point<6, FaceTexScalar>, 6, FaceTexScalar>(NNP_FACE_WEDGE_TEX, (*wedgeTexCoord.begin()).V()));
				}
				if ((bitMask & BitMask::IO_WEDGTEXMULTI) && vcg::tri::HasPerWedgeTexCoord(mesh))
					faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 1, short>(NNP_TEXTUREINDEX, &(*mesh.face.begin()).WT(0).N()));
				if ((bitMask & BitMask::IO_WEDGCOLOR) && vcg::tri::HasPerWedgeColor(mesh))
					faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 12, WedgeColorScalar>(NNP_FACE_WEDGE_COLOR, (*mesh.face.begin()).WC(0).V()));
				if ((bitMask & BitMask::IO_WEDGNORMAL) && vcg::tri::HasPerWedgeNormal(mesh))
					faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 9, WedgeNormalScalar>(NNP_FACE_WEDGE_NORMAL, (*mesh.face.begin()).WN(0).V()));
				if ((bitMask & BitMask::IO_FACEATTRIB) && custom.faceAttrib.size() > 0)
				{
					for (int i = 0; i < custom.faceAttrib.size(); i++)
					{
            std::set<PointerToAttribute>::iterator ai;
						for (ai = mesh.face_attr.begin(); ai != mesh.face_attr.end(); ++ai)
						{
							if ((*custom.faceAttrib[i]).name == (*ai)._name)
							{
								custom.faceAttrib[i]->base = ai->_handle->DataBegin();
								break;
							}
						}
						faceDescr.dataDescriptor.push_back(custom.faceAttrib[i]);
							
					}
				}

			}

			std::vector<ElementDescriptor*> meshDescr;
			meshDescr.push_back(&cameraDescr);
			meshDescr.push_back(&vertexDescr);
			meshDescr.push_back(&edgeDescr);
			meshDescr.push_back(&faceDescr);

			//Mesh attribute
			if ((bitMask & BitMask::IO_MESHATTRIB))
			{
				CustomAttributeDescriptor::MapMeshAttribIter iter = custom.meshAttrib.begin();
				for (; iter != custom.meshAttrib.end(); iter++)
				{
					std::string name((*iter).first);
					meshDescr.push_back(new ElementDescriptor(name));
					count = info.GetElementCount(name);
					if (count > 1)
					{
						meshDescr.back()->dataDescriptor = (*iter).second;
					}
				}

			}
			if (!OpenModel(info, meshDescr))
				return info.errInfo;

			mesh.shot.SetViewPoint(tra);
			mesh.shot.Extrinsics.SetRot(rot);
      for (int i = 0; i < faceIndex.size(); i++)
				for (int j = 0; j < 3; j++)
					mesh.face[i].V(j) = &mesh.vert[faceIndex[i][j]];
      for (int i = 0; i < wedgeTexCoord.size(); i++)
      {
        for (int j = 0; j < 3; j++)
        {
          mesh.face[i].WT(j).U() = wedgeTexCoord[i][j * 2];
          mesh.face[i].WT(j).V() = wedgeTexCoord[i][j * 2 + 1];
        }
      }
			for (int i = 0; i < edgeIndex.size(); i++)
			{
				mesh.edge[i].V(0) = &mesh.vert[edgeIndex[i].X()];
				mesh.edge[i].V(1) = &mesh.vert[edgeIndex[i].Y()];
			}

			for (int i = 0; i < cameraDescr.dataDescriptor.size(); i++)
				delete cameraDescr.dataDescriptor[i];
			for (int i = 0; i < vertexDescr.dataDescriptor.size(); i++)
				if (vertexDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete vertexDescr.dataDescriptor[i];
			for (int i = 0; i < edgeDescr.dataDescriptor.size(); i++)
				if (edgeDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete edgeDescr.dataDescriptor[i];
			for (int i = 0; i < faceDescr.dataDescriptor.size(); i++)
				if (faceDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete faceDescr.dataDescriptor[i];
			mesh.textures = info.textureFile;
			return info.errInfo;
		}


		static int LoadModel(const char* filename, MeshType& mesh, unsigned int bitMask)
		{
			CustomAttributeDescriptor custom;
			return LoadModel(filename, mesh, bitMask, custom);
		}
	

		static bool SaveModel(const char* filename, MeshType& mesh, unsigned int bitMask, CustomAttributeDescriptor& custom, bool binary)
		{
			//Camera
			std::vector<PlyProperty> cameraProp;
			ElementDescriptor cameraDescr(std::string("camera"));
			vcg::Point3<ScalarType> tra = mesh.shot.Extrinsics.Tra();
			vcg::Matrix44<ScalarType> rot = mesh.shot.Extrinsics.Rot();
			if (bitMask & BitMask::IO_CAMERA)
			{
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("view_px"), &tra[0]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("view_py"), &tra[1]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("view_pz"), &tra[2]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("x_axisx"), &rot[0][0]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("x_axisy"), &rot[0][1]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("x_axisz"), &rot[0][2]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("y_axisx"), &rot[1][0]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("y_axisy"), &rot[1][1]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("y_axisz"), &rot[1][2]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("z_axisx"), &rot[2][0]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("z_axisy"), &rot[2][1]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("z_axisz"), &rot[2][2]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("focal"), &mesh.shot.Intrinsics.FocalMm);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("scalex"), &mesh.shot.Intrinsics.PixelSizeMm[0]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("scaley"), &mesh.shot.Intrinsics.PixelSizeMm[1]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("centerx"), &mesh.shot.Intrinsics.CenterPx[0]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("centery"), &mesh.shot.Intrinsics.CenterPx[1]);
				PushDescriport<int, int, 1>(cameraProp, cameraDescr, std::string("viewportx"), &mesh.shot.Intrinsics.ViewportPx[0]);
				PushDescriport<int, int, 1>(cameraProp, cameraDescr, std::string("viewporty"), &mesh.shot.Intrinsics.ViewportPx[1]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("k1"), &mesh.shot.Intrinsics.k[0]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("k2"), &mesh.shot.Intrinsics.k[1]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("k3"), &mesh.shot.Intrinsics.k[2]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("k4"), &mesh.shot.Intrinsics.k[3]);
			}

			//Vertex
			std::vector<std::string> nameList;
			VertexType::Name(nameList);
			std::vector<PlyProperty> vertexProp;
			ElementDescriptor vertexDescr(NNP_VERTEX_ELEM);
			if (nameList.size() > 0 && mesh.vert.size() > 0)
			{
				if ((bitMask & BitMask::IO_VERTCOORD) && VertexType::HasCoord())
					PushDescriport<VertexType, VertexCoordScalar, 3>(vertexProp, vertexDescr, NNP_PXYZ, (*mesh.vert.begin()).P().V());
				if ((bitMask & BitMask::IO_VERTNORMAL) && vcg::tri::HasPerVertexNormal(mesh))
					PushDescriport<VertexType, VertexNormScalar, 3>(vertexProp, vertexDescr, NNP_NXYZ, (*mesh.vert.begin()).N().V());
				if ((bitMask & BitMask::IO_VERTCOLOR) && vcg::tri::HasPerVertexColor(mesh))
					PushDescriport<VertexType, VertexColorScalar, 4>(vertexProp, vertexDescr, NNP_CRGBA, (*mesh.vert.begin()).C().V());
				if ((bitMask & BitMask::IO_VERTQUALITY) && vcg::tri::HasPerVertexQuality(mesh))
					PushDescriport<VertexType, VertexQuality, 1>(vertexProp, vertexDescr, NNP_QUALITY, &(*mesh.vert.begin()).Q());
				if ((bitMask & BitMask::IO_VERTFLAGS) && vcg::tri::HasPerVertexFlags(mesh))
					PushDescriport<VertexType, VertexFlag, 1>(vertexProp, vertexDescr, NNP_BITFLAG, &(*mesh.vert.begin()).Flags());
				if ((bitMask & BitMask::IO_VERTRADIUS) && vcg::tri::HasPerVertexRadius(mesh))
					PushDescriport<VertexType, VertexRadius, 1>(vertexProp, vertexDescr, NNP_DENSITY, &(*mesh.vert.begin()).R());
				if ((bitMask & BitMask::IO_VERTTEXCOORD) && vcg::tri::HasPerVertexTexCoord(mesh))
					PushDescriport<VertexType, VertexTexScalar, 2>(vertexProp, vertexDescr, NNP_TEXTURE2D, (*mesh.vert.begin()).T().P().V());
				if ((bitMask & BitMask::IO_VERTCURV) && vcg::tri::HasPerVertexCurvature(mesh))
				{
					PushDescriport<VertexType, VertexCurScalar, 1>(vertexProp, vertexDescr, NNP_KG, &(*mesh.vert.begin()).Kg());
          PushDescriport<VertexType, VertexCurScalar, 1>(vertexProp, vertexDescr, NNP_KH, &(*mesh.vert.begin()).Kh());
				}
				if ((bitMask & BitMask::IO_VERTCURVDIR) && vcg::tri::HasPerVertexCurvatureDir(mesh))
				{
					PushDescriport<VertexType, VertexDirCurScalar, 1>(vertexProp, vertexDescr, NNP_K1, &(*mesh.vert.begin()).K1());
          PushDescriport<VertexType, VertexDirCurScalar, 1>(vertexProp, vertexDescr, NNP_K2, &(*mesh.vert.begin()).K2());
					PushDescriportList<VertexType, VertexDirCurVecScalar, 3>(vertexProp, vertexDescr, NNP_K1DIR, (*mesh.vert.begin()).PD1().V());
          PushDescriportList<VertexType, VertexDirCurVecScalar, 3>(vertexProp, vertexDescr, NNP_K2DIR, (*mesh.vert.begin()).PD2().V());
				}
				if ((bitMask & BitMask::IO_VERTATTRIB) && custom.vertexAttrib.size() > 0)
				{
					for (int i = 0; i < custom.vertexAttrib.size(); i++)
					{
						vertexProp.push_back(custom.vertexAttribProp[i]);
						vertexDescr.dataDescriptor.push_back(custom.vertexAttrib[i]);
					}
				}
			}

			//Edge
			nameList.clear();
			EdgeType::Name(nameList);
			std::vector<PlyProperty> edgeProp;
			ElementDescriptor edgeDescr(NNP_VERTEX_ELEM);
			std::vector<vcg::Point2i> edgeIndex;
			for (int i = 0; i < mesh.edge.size(); i++)
				edgeIndex.push_back(vcg::Point2i(vcg::tri::Index(mesh, mesh.edge[i].V(0)), vcg::tri::Index(mesh, mesh.edge[i].V(1))));
			if (nameList.size() > 0 && mesh.edge.size() > 0)
			{
				if ((bitMask & BitMask::IO_EDGEINDEX) && EdgeType::HasVertexRef())
				{
					PushDescriport<vcg::Point2i, int, 1>(edgeProp, edgeDescr, NNP_EDGE_V1, &(*edgeIndex.begin()).V()[0]);
					PushDescriport<vcg::Point2i, int, 1>(edgeProp, edgeDescr, NNP_EDGE_V2, &(*edgeIndex.begin()).V()[1]);
				}
				if ((bitMask & BitMask::IO_EDGEQUALITY) && vcg::tri::HasPerEdgeQuality(mesh))
					PushDescriport<EdgeType, EdgeQuality, 1>(edgeProp, edgeDescr, NNP_QUALITY, &(*mesh.edge.begin()).Q());
				if ((bitMask & BitMask::IO_EDGECOLOR) && vcg::tri::HasPerEdgeColor(mesh))
					PushDescriport<EdgeType, EdgeColorScalar, 4>(edgeProp, edgeDescr, NNP_CRGBA, (*mesh.edge.begin()).C().V());
				if ((bitMask & BitMask::IO_EDGEFLAGS) && vcg::tri::HasPerEdgeFlags(mesh))
					PushDescriport<EdgeType, EdgeFlag, 1>(edgeProp, edgeDescr, NNP_BITFLAG, &(*mesh.edge.begin()).Flags());
				if ((bitMask & BitMask::IO_EDGEATTRIB) && custom.edgeAttrib.size() > 0)
				{
					for (int i = 0; i < custom.edgeAttrib.size(); i++)
					{
						edgeProp.push_back(custom.edgeAttribProp[i]);
						edgeDescr.dataDescriptor.push_back(custom.edgeAttrib[i]);
					}
				}
			}
			
			//Face
			nameList.clear();
			FaceType::Name(nameList);
			std::vector<PlyProperty> faceProp;
			ElementDescriptor faceDescr(NNP_FACE_ELEM);
			std::vector<vcg::Point3i> faceIndex;
			std::vector<vcg::ndim::Point<6, FaceTexScalar>> wedgeTexCoord;
			for (int i = 0; i < mesh.face.size(); i++)
				faceIndex.push_back(vcg::Point3i(vcg::tri::Index(mesh, mesh.face[i].V(0)), vcg::tri::Index(mesh, mesh.face[i].V(1)), vcg::tri::Index(mesh, mesh.face[i].V(2))));
			
			if (((bitMask & BitMask::IO_WEDGTEXCOORD) || (bitMask & BitMask::IO_WEDGTEXMULTI)) && vcg::tri::HasPerWedgeTexCoord(mesh))
			{
				for (int i = 0; i < mesh.face.size(); i++)
				{
          wedgeTexCoord.push_back(vcg::ndim::Point<6, FaceTexScalar>());
					for (int j = 0; j < 3; j++)
					{
						wedgeTexCoord.back()[j * 2] = mesh.face[i].WT(j).U();
						wedgeTexCoord.back()[j * 2 + 1] = mesh.face[i].WT(j).V();
					}
				}
			}
			if (nameList.size() > 0 && mesh.face.size() > 0)
			{
				if ((bitMask & BitMask::IO_FACEINDEX) && FaceType::HasVertexRef())
					PushDescriportList<vcg::Point3i, int, 3>(faceProp, faceDescr, NNP_FACE_VERTEX_LIST, (*faceIndex.begin()).V());
				if ((bitMask & BitMask::IO_FACEFLAGS) && vcg::tri::HasPerFaceFlags(mesh))
					PushDescriport<FaceType, FaceFlag, 1>(faceProp, faceDescr, NNP_BITFLAG, &(*mesh.face.begin()).Flags());
				if ((bitMask & BitMask::IO_FACECOLOR) && vcg::tri::HasPerFaceColor(mesh))
					PushDescriport<FaceType, FaceColorScalar, 4>(faceProp, faceDescr, NNP_CRGBA, (*mesh.face.begin()).C().V());
				if ((bitMask & BitMask::IO_FACEQUALITY) && vcg::tri::HasPerFaceQuality(mesh))
					PushDescriport<FaceType, FaceQuality, 1>(faceProp, faceDescr, NNP_QUALITY, &(*mesh.face.begin()).Q());
				if ((bitMask & BitMask::IO_FACENORMAL) && vcg::tri::HasPerFaceNormal(mesh))
					PushDescriport<FaceType, FaceNormScalar, 3>(faceProp, faceDescr, NNP_NXYZ, (*mesh.face.begin()).N().V());
				if ((bitMask & BitMask::IO_VERTCURVDIR) && vcg::tri::HasPerFaceCurvatureDir(mesh))
				{
					PushDescriport<FaceType, FaceDirCurScalar, 1>(faceProp, faceDescr, NNP_K1, &(*mesh.face.begin()).K1());
          PushDescriport<FaceType, FaceDirCurScalar, 1>(faceProp, faceDescr, NNP_K2, &(*mesh.face.begin()).K2());
					PushDescriportList<FaceType, FaceDirCurVecScalar, 3>(faceProp, faceDescr, NNP_K1DIR, (*mesh.face.begin()).PD1().V());
          PushDescriportList<FaceType, FaceDirCurVecScalar, 3>(faceProp, faceDescr, NNP_K2DIR, (*mesh.face.begin()).PD2().V());
				}
				if (((bitMask & BitMask::IO_WEDGTEXCOORD) || (bitMask & BitMask::IO_WEDGTEXMULTI)) && vcg::tri::HasPerWedgeTexCoord(mesh))
          PushDescriportList<vcg::ndim::Point<6, FaceTexScalar>, FaceTexScalar, 6>(faceProp, faceDescr, NNP_FACE_WEDGE_TEX, (*wedgeTexCoord.begin()).V());
				if ((bitMask & BitMask::IO_WEDGTEXMULTI) && vcg::tri::HasPerWedgeTexCoord(mesh))
					PushDescriport<FaceType, short, 1>(faceProp, faceDescr, NNP_TEXTUREINDEX, &(*mesh.face.begin()).WT(0).N());
				if ((bitMask & BitMask::IO_WEDGCOLOR) && vcg::tri::HasPerWedgeColor(mesh))
					PushDescriportList<FaceType, WedgeColorScalar, 12>(faceProp, faceDescr, NNP_FACE_WEDGE_COLOR, (*mesh.face.begin()).WC(0).V());
				if ((bitMask & BitMask::IO_WEDGNORMAL) && vcg::tri::HasPerWedgeNormal(mesh))
					PushDescriportList<FaceType, WedgeNormalScalar, 9>(faceProp, faceDescr, NNP_FACE_WEDGE_NORMAL, (*mesh.face.begin()).WN(0).V());
				if ((bitMask & BitMask::IO_FACEATTRIB) && custom.faceAttrib.size() > 0)
				{
					for (int i = 0; i < custom.faceAttrib.size(); i++)
					{
						faceProp.push_back(custom.faceAttribProp[i]);
						faceDescr.dataDescriptor.push_back(custom.faceAttrib[i]);
					}
				}
				
			}
			
			Info infoSave;
			infoSave.filename = filename;
			infoSave.binary = binary;
			PlyElement cameraElem(std::string("camera"), cameraProp, 1);
			PlyElement vertexElem(NNP_VERTEX_ELEM, vertexProp, mesh.vert.size());
			PlyElement edgeElem(NNP_EDGE_ELEM, edgeProp, mesh.edge.size());
			PlyElement faceElem(NNP_FACE_ELEM, faceProp, mesh.face.size());
			infoSave.AddPlyElement(cameraElem);
			infoSave.AddPlyElement(vertexElem);
			infoSave.AddPlyElement(edgeElem);
			infoSave.AddPlyElement(faceElem);
			infoSave.textureFile = mesh.textures;
			std::vector<ElementDescriptor*> meshDescr;
			meshDescr.push_back(&cameraDescr);
			meshDescr.push_back(&vertexDescr);
			meshDescr.push_back(&edgeDescr);
			meshDescr.push_back(&faceDescr);

			//Mesh attribute
			if ((bitMask & BitMask::IO_MESHATTRIB))
			{
				CustomAttributeDescriptor::MapMeshAttribIter iter = custom.meshAttrib.begin();
				CustomAttributeDescriptor::MapMeshAttribPropIter iterProp = custom.meshAttribProp.begin();
				for (; iter != custom.meshAttrib.end(); iter++, iterProp++)
				{
					std::string name((*iter).first);
					PlyElement customElem(name, (*iterProp).second, custom.meshAttribCnt[(*iter).first]);
					infoSave.AddPlyElement(customElem);
					meshDescr.push_back(new ElementDescriptor(name));
					meshDescr.back()->dataDescriptor = (*iter).second;
				}
			}
			
			bool flag = nanoply::SaveModel(infoSave.filename, meshDescr, infoSave);

			for (int i = 0; i < cameraDescr.dataDescriptor.size(); i++)
				delete cameraDescr.dataDescriptor[i];
			for (int i = 0; i < vertexDescr.dataDescriptor.size(); i++)
				if (vertexDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete vertexDescr.dataDescriptor[i];
			for (int i = 0; i < edgeDescr.dataDescriptor.size(); i++)
				if (edgeDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete edgeDescr.dataDescriptor[i];
			for (int i = 0; i < faceDescr.dataDescriptor.size(); i++)
				if (faceDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete faceDescr.dataDescriptor[i];
			
			return flag;
		}
		

		static bool SaveModel(const char* filename, MeshType& mesh, unsigned int bitMask, bool binary)
		{
			CustomAttributeDescriptor custom;
			return SaveModel(filename, mesh, bitMask, custom, binary);
		}


	};

}

#endif