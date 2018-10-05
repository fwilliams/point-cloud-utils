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

#include <iostream>
#include <wrap/nanoply/include/nanoplyWrapper.hpp>

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/normal.h>

struct Material
{
  vcg::Point3f kd;
  vcg::Point3f ks;
  float rho;
};


class MyVertex;
class MyFace;

class MyUsedTypes : public vcg::UsedTypes < vcg::Use< MyVertex >::AsVertexType, vcg::Use< MyFace >::AsFaceType>{};
class MyVertex : public vcg::Vertex <  MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::Radiusf, vcg::vertex::BitFlags>{};
class MyFace : public vcg::Face <  MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::BitFlags>{};

class MyMesh : public vcg::tri::TriMesh < std::vector< MyVertex >, std::vector< MyFace > >
{
public:
  MyMesh::PerVertexAttributeHandle<int> vertexMaterial;
  MyMesh::PerFaceAttributeHandle<MyVertex::CoordType> faceBarycenter;
  MyMesh::PerMeshAttributeHandle<std::vector<Material>> material;
    
  MyMesh()
  {
    vertexMaterial = vcg::tri::Allocator<MyMesh>::AddPerVertexAttribute<int>(*this, std::string("materialId"));
    faceBarycenter = vcg::tri::Allocator<MyMesh>::AddPerFaceAttribute<MyVertex::CoordType>(*this, std::string("barycenter"));
    material = vcg::tri::Allocator<MyMesh>::AddPerMeshAttribute<std::vector<Material>>(*this, std::string("material"));
  }

  void FillMesh()
  {
    vcg::tri::Icosahedron(*this);
    vcg::tri::UpdateNormal<MyMesh>::PerFaceNormalized(*this);
    vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalized(*this);

    int tempC = 255 / vert.size();
    for (int i = 0; i < vert.size(); i++)
    {
      if (i < 2)
        vertexMaterial[i] = -1;
      else if (i < vert.size() / 2)
        vertexMaterial[i] = 0;
      else
        vertexMaterial[i] = 1;
      vert[i].R() = i*i / 2.0f;
      vert[i].C() = vcg::Color4b(tempC*i, tempC*i, tempC*i, 255);
    }

    for (int i = 0; i < face.size(); i++)
     faceBarycenter[i] = vcg::Barycenter(face[i]);

    material().resize(2);
    material()[0] = { vcg::Point3f(0.1f, 0.2f, 0.3f), vcg::Point3f(0.3f, 0.3f, 0.3f), 5.0f };
    material()[1] = { vcg::Point3f(0.1f, 0.1f, 0.1f), vcg::Point3f(0.5f, 0.3f, 0.4f), 50.0f };
  }

};


bool Load(const char* filename, MyMesh& mesh)
{
  //Create the data descriport for the custom attributes
  nanoply::NanoPlyWrapper<MyMesh>::CustomAttributeDescriptor customAttrib;
  customAttrib.GetMeshAttrib(filename);
  int count = customAttrib.meshAttribCnt["material"];
  mesh.material().resize(count);
  customAttrib.AddVertexAttribDescriptor<int, int, 1>(std::string("materialId"), nanoply::NNP_INT32, NULL);
  customAttrib.AddFaceAttribDescriptor<vcg::Point3f, float, 3>(std::string("barycenter"), nanoply::NNP_LIST_UINT8_FLOAT32, NULL);
	if (count > 0)
	{
		customAttrib.AddMeshAttribDescriptor<Material, float, 3>(std::string("material"), std::string("kd"), nanoply::NNP_FLOAT32, mesh.material()[0].kd.V());
		customAttrib.AddMeshAttribDescriptor<Material, float, 3>(std::string("material"), std::string("ks"), nanoply::NNP_FLOAT32, mesh.material()[0].ks.V());
		customAttrib.AddMeshAttribDescriptor<Material, float, 1>(std::string("material"), std::string("rho"), nanoply::NNP_FLOAT32, &mesh.material()[0].rho);
	}

  //Load the ply file
  unsigned int mask = 0;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_VERTCOORD;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_VERTNORMAL;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_VERTCOLOR;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_VERTRADIUS;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_VERTATTRIB;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_FACEINDEX;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_FACENORMAL;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_FACEATTRIB;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_MESHATTRIB;
	return (nanoply::NanoPlyWrapper<MyMesh>::LoadModel(filename, mesh, mask, customAttrib) != 0);
}



bool Save(const char* filename, MyMesh& mesh, bool binary)
{
  //Create the data descriport for the custom attributes
  nanoply::NanoPlyWrapper<MyMesh>::CustomAttributeDescriptor customAttrib;
  customAttrib.AddVertexAttribDescriptor<int, int, 1>(std::string("materialId"), nanoply::NNP_INT32, &mesh.vertexMaterial[0]);
  customAttrib.AddFaceAttribDescriptor<vcg::Point3f, float, 3>(std::string("barycenter"), nanoply::NNP_LIST_UINT8_FLOAT32, mesh.faceBarycenter[0].V());
	if (mesh.material().size() > 0)
	{
		customAttrib.AddMeshAttrib(std::string("material"), mesh.material().size());
		customAttrib.AddMeshAttribDescriptor<Material, float, 3>(std::string("material"), std::string("kd"), nanoply::NNP_LIST_UINT8_FLOAT32, mesh.material()[0].kd.V());
		customAttrib.AddMeshAttribDescriptor<Material, float, 3>(std::string("material"), std::string("ks"), nanoply::NNP_LIST_UINT8_FLOAT32, mesh.material()[0].ks.V());
		customAttrib.AddMeshAttribDescriptor<Material, float, 1>(std::string("material"), std::string("rho"), nanoply::NNP_FLOAT32, &mesh.material()[0].rho);
	}

  //Save the ply file
  unsigned int mask = 0;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_VERTCOORD;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_VERTNORMAL;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_VERTCOLOR;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_VERTRADIUS;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_VERTATTRIB;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_FACEINDEX;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_FACENORMAL;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_FACEATTRIB;
  mask |= nanoply::NanoPlyWrapper<MyMesh>::IO_MESHATTRIB;
  return nanoply::NanoPlyWrapper<MyMesh>::SaveModel(filename, mesh, mask, customAttrib, binary);
}



int main()
{
  MyMesh mesh1;
  mesh1.FillMesh();
  Save("example_ascii.ply", mesh1, false);
  Save("example_binary.ply", mesh1, true);
  MyMesh mesh2, mesh3;
  Load("example_ascii.ply", mesh2);
  Load("example_binary.ply", mesh3);
  return true;
}
