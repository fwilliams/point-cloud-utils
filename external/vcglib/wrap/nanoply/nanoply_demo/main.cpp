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
#include <nanoply.hpp>

template<typename T, int N>
struct Container
{

public:
  T data[N];

  Container(){}

  Container(T* temp, int n)
  {
    for (int i = 0; i < std::min(n, N); i++)
      data[i] = temp[i];
  }

  T* V()
  {
    return data;
  }

  bool operator == (Container<T, N> const & m) const
  {
    bool flag = true;
    for (int i = 0; i < N; i++)
      flag = flag && (data[i] == m.data[i]);
    return flag;
  }
};


typedef Container<float, 3> Point3f;
typedef Container<unsigned char, 4> Color4f;
typedef Container<int, 3> VertexIndex;

struct MyVertexInfo
{
  Color4f c;
  float density;
  int materialId;

  bool operator == (MyVertexInfo const & m) const
  {
    return (c == m.c && m.density == density && m.materialId == materialId);
  }

};

struct MyMaterialInfo
{
  Point3f kd;
  Point3f ks;
  float rho;

  bool operator == (MyMaterialInfo const & m) const
  {
    return (kd == m.kd && ks == m.ks && rho == m.rho);
  }
};


class MyMesh
{
public:
  std::vector<Point3f> coordVec;
  std::vector<Point3f> normalVec;
  std::vector<MyVertexInfo> infoVec;
  std::vector<VertexIndex> faceIndex;
  std::vector<MyMaterialInfo> material;

  void FillMesh()
  {
    float pos[] = { 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0 };
    int index[] = { 0, 1, 2, 0, 2, 3, 0, 3, 1, 3, 2, 1 };
    float norm[] = { 0.57735, 0.57735, 0.57735, -0.57735, 0.57735, -0.57735, -0.57735, -0.57735, 0.57735, 0.57735, -0.57735, -0.57735 };
    unsigned char color[] = { 68, 68, 68, 255, 177, 68, 177, 255, 177, 177, 68, 255, 68, 177, 177 };
    float density[] = { 3.5, 2.0, 4.0, 3.0 };
    float materialId[] = { 1, 0, -1, 1 };
    float materialValue[] = { 0.2, 0.3, 0.2, 0.5, 0.5, 0.6, 20.0, 0.1, 0.1, 0.1, 0.7, 0.5, 0.4, 1.0 };
    coordVec.push_back(Point3f(pos, 3)); coordVec.push_back(Point3f(&pos[3], 3)); coordVec.push_back(Point3f(&pos[6], 3)); coordVec.push_back(Point3f(&pos[9], 3));
    normalVec.push_back(Point3f(norm, 3)); normalVec.push_back(Point3f(&norm[3], 3)); normalVec.push_back(Point3f(&norm[6], 3)); normalVec.push_back(Point3f(&norm[9], 3));
    MyVertexInfo info1 = { Color4f(color, 4), density[0], materialId[0] }; infoVec.push_back(info1);
    MyVertexInfo info2 = { Color4f(&color[4], 4), density[1], materialId[1] }; infoVec.push_back(info2);
    MyVertexInfo info3 = { Color4f(&color[8], 4), density[2], materialId[2] }; infoVec.push_back(info3);
    MyVertexInfo info4 = { Color4f(&color[12], 4), density[3], materialId[3] }; infoVec.push_back(info4);
    faceIndex.push_back(VertexIndex(index, 3));	faceIndex.push_back(VertexIndex(&index[3], 3)); faceIndex.push_back(VertexIndex(&index[6], 3)); faceIndex.push_back(VertexIndex(&index[9], 3));
    MyMaterialInfo mat1 = { Point3f(materialValue, 3), Point3f(&materialValue[3], 3), materialValue[6] }; material.push_back(mat1);
    MyMaterialInfo mat2 = { Point3f(&materialValue[7], 3), Point3f(&materialValue[10], 3), materialValue[13] }; material.push_back(mat2);
  }

  bool operator == (MyMesh& m)
  {
    bool flag = (coordVec == m.coordVec);
    flag = flag && (normalVec == m.normalVec);
    flag = flag && (infoVec == m.infoVec);
    flag = flag && (faceIndex == m.faceIndex);
    flag = flag && (material == m.material);
    return flag;
  }
};


bool Load(const char* filename, MyMesh& mesh)
{
  //Get file info
  nanoply::Info info(filename);
  if (info.errInfo != nanoply::NNP_OK)
  {
    std::cout << "Invalid file format" << std::endl;
    return false;
  }

  //Resize the element containers
  int vertCnt = info.GetVertexCount();
  if (vertCnt <= 0)
  {
    std::cout << "The file does't contain any vertex." << std::endl;
    return false;
  }
  mesh.coordVec.resize(vertCnt);
  mesh.normalVec.resize(vertCnt);
  mesh.infoVec.resize(vertCnt);
  int faceCnt = info.GetFaceCount();
  mesh.faceIndex.resize(faceCnt);
  int materialCnt = info.GetElementCount(std::string("material"));
  mesh.material.resize(2);

  //Create the vertex properties descriptor (what ply property and where to save its data)
  nanoply::ElementDescriptor vertex(nanoply::NNP_VERTEX_ELEM);
  if (vertCnt > 0)
  {
    vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<Point3f, 3, float>(nanoply::NNP_PXYZ, (*mesh.coordVec.begin()).V()));
    vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<Point3f, 3, float>(nanoply::NNP_NXYZ, (*mesh.normalVec.begin()).V()));
    vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<MyVertexInfo, 4, unsigned char>(nanoply::NNP_CRGBA, (*mesh.infoVec.begin()).c.V()));
    vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<MyVertexInfo, 1, float>(nanoply::NNP_DENSITY, &(*mesh.infoVec.begin()).density));
    vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<MyVertexInfo, 1, int>(std::string("materialId"), &(*mesh.infoVec.begin()).materialId));
  }

  //Create the face properties descriptor (what ply property and where the data is stored)
  nanoply::ElementDescriptor face(nanoply::NNP_FACE_ELEM);
  if (mesh.faceIndex.size() > 0)
    face.dataDescriptor.push_back(new nanoply::DataDescriptor<VertexIndex, 3, int>(nanoply::NNP_FACE_VERTEX_LIST, (*mesh.faceIndex.begin()).V()));

  //Create the material properties descriptor (what ply property and where the data is stored)
  nanoply::ElementDescriptor material(std::string("material"));
  if (mesh.material.size() > 0)
  {
    material.dataDescriptor.push_back(new nanoply::DataDescriptor<MyMaterialInfo, 3, float>(std::string("kd"), (*mesh.material.begin()).kd.V()));
    material.dataDescriptor.push_back(new nanoply::DataDescriptor<MyMaterialInfo, 3, float>(std::string("ks"), (*mesh.material.begin()).ks.V()));
    material.dataDescriptor.push_back(new nanoply::DataDescriptor<MyMaterialInfo, 1, float>(std::string("rho"), &(*mesh.material.begin()).rho));
  }

  //Create the mesh descriptor
  std::vector<nanoply::ElementDescriptor*> meshDescr;
  meshDescr.push_back(&vertex);
  meshDescr.push_back(&face);
  meshDescr.push_back(&material);

  //Open the file and save the element data according the relative element descriptor
  OpenModel(info, meshDescr);
  for (int i = 0; i < vertex.dataDescriptor.size(); i++)
    delete vertex.dataDescriptor[i];
  for (int i = 0; i < face.dataDescriptor.size(); i++)
    delete face.dataDescriptor[i];
  for (int i = 0; i < material.dataDescriptor.size(); i++)
    delete material.dataDescriptor[i];
  return (info.errInfo == nanoply::NNP_OK);
}



bool Save(const char* filename, MyMesh& mesh, bool binary)
{
  //Create the vector of vertex properties to save in the file
  std::vector<nanoply::PlyProperty> vertexProp;
  vertexProp.push_back(nanoply::PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_PXYZ));
  vertexProp.push_back(nanoply::PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_NXYZ));
  vertexProp.push_back(nanoply::PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_DENSITY));
  vertexProp.push_back(nanoply::PlyProperty(nanoply::NNP_FLOAT32, nanoply::NNP_CRGBA));
  vertexProp.push_back(nanoply::PlyProperty(nanoply::NNP_INT32, "materialId"));

  //Create the vector of face properties to save in the file
  std::vector<nanoply::PlyProperty> faceProp;
  faceProp.push_back(nanoply::PlyProperty(nanoply::NNP_LIST_UINT8_UINT32, nanoply::NNP_FACE_VERTEX_LIST));

  //Create the vector of material properties to save in the file
  std::vector<nanoply::PlyProperty> materialProp;
  materialProp.push_back(nanoply::PlyProperty(nanoply::NNP_LIST_UINT8_FLOAT32, "kd"));
  materialProp.push_back(nanoply::PlyProperty(nanoply::NNP_LIST_UINT8_FLOAT32, "ks"));
  materialProp.push_back(nanoply::PlyProperty(nanoply::NNP_FLOAT32, "rho"));

  //Create the PlyElement
  nanoply::PlyElement vertexElem(nanoply::NNP_VERTEX_ELEM, vertexProp, mesh.coordVec.size());
  nanoply::PlyElement faceElem(nanoply::NNP_FACE_ELEM, faceProp, mesh.faceIndex.size());
  nanoply::PlyElement materialElem(std::string("material"), materialProp, mesh.material.size());

  //Create the Info object with the data to save in the header
  nanoply::Info infoSave;
  infoSave.filename = filename;
  infoSave.binary = binary;
  infoSave.AddPlyElement(vertexElem);
  infoSave.AddPlyElement(faceElem);
  infoSave.AddPlyElement(materialElem);

  //Create the vertex properties descriptor (what ply property and where the data is stored)
  nanoply::ElementDescriptor vertex(nanoply::NNP_VERTEX_ELEM);
  if (mesh.coordVec.size() > 0)
  {
    vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<Point3f, 3, float>(nanoply::NNP_PXYZ, (*mesh.coordVec.begin()).V()));
    vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<Point3f, 3, float>(nanoply::NNP_NXYZ, (*mesh.normalVec.begin()).V()));
    vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<MyVertexInfo, 4, unsigned char>(nanoply::NNP_CRGBA, (*mesh.infoVec.begin()).c.V()));
    vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<MyVertexInfo, 1, float>(nanoply::NNP_DENSITY, &(*mesh.infoVec.begin()).density));
    vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<MyVertexInfo, 1, int>(std::string("materialId"), &(*mesh.infoVec.begin()).materialId));
  }

  //Create the face properties descriptor (what ply property and where the data is stored)
  nanoply::ElementDescriptor face(nanoply::NNP_FACE_ELEM);
  if (mesh.faceIndex.size() > 0)
    face.dataDescriptor.push_back(new nanoply::DataDescriptor<VertexIndex, 3, int>(nanoply::NNP_FACE_VERTEX_LIST, (*mesh.faceIndex.begin()).V()));

  //Create the material properties descriptor (what ply property and where the data is stored)
  nanoply::ElementDescriptor material(std::string("material"));
  if (mesh.material.size() > 0)
  {
    material.dataDescriptor.push_back(new nanoply::DataDescriptor<MyMaterialInfo, 3, float>(std::string("kd"), (*mesh.material.begin()).kd.V()));
    material.dataDescriptor.push_back(new nanoply::DataDescriptor<MyMaterialInfo, 3, float>(std::string("ks"), (*mesh.material.begin()).ks.V()));
    material.dataDescriptor.push_back(new nanoply::DataDescriptor<MyMaterialInfo, 1, float>(std::string("rho"), &(*mesh.material.begin()).rho));
  }

  //Create the mesh descriptor
  std::vector<nanoply::ElementDescriptor*> meshDescr;
  meshDescr.push_back(&vertex);
  meshDescr.push_back(&face);
  meshDescr.push_back(&material);

  //Save the file
  bool result = nanoply::SaveModel(infoSave.filename, meshDescr, infoSave);

  for (int i = 0; i < vertex.dataDescriptor.size(); i++)
    delete vertex.dataDescriptor[i];
  for (int i = 0; i < face.dataDescriptor.size(); i++)
    delete face.dataDescriptor[i];
  for (int i = 0; i < material.dataDescriptor.size(); i++)
    delete material.dataDescriptor[i];
  return result;
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
  if (mesh2 == mesh1)
    std::cout << "Write and read ASCII ply file:  SUCCESS\n";
  else
    std::cout << "Write and read ASCII ply file:  FAIL\n";
  if (mesh3 == mesh1)
    std::cout << "Write and read binary ply file:  SUCCESS\n";
  else
    std::cout << "Write and read binary ply file:  FAIL\n";
  return true;
}
