/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#ifndef MESH_TO_MATRIX
#define MESH_TO_MATRIX

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/harmonic.h>

using namespace std;

namespace vcg {
namespace tri {
template < typename MeshType >
class MeshToMatrix
{

    // define types

    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::FaceIterator FaceIterator;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::VertexIterator VertexIterator;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> MatrixXm;

    static void GetTriEdgeAdjacency(const MatrixXm& V,
                                    const Eigen::MatrixXi& F,
                                    Eigen::MatrixXi& EV,
                                    Eigen::MatrixXi& FE,
                                    Eigen::MatrixXi& EF)
    {
        (void)V;
        //assert(igl::is_manifold(V,F));
        std::vector<std::vector<int> > ETT;
        for(int f=0;f<F.rows();++f)
            for (int i=0;i<3;++i)
            {
                // v1 v2 f vi
                int v1 = F(f,i);
                int v2 = F(f,(i+1)%3);
                if (v1 > v2) std::swap(v1,v2);
                std::vector<int> r(4);
                r[0] = v1; r[1] = v2;
                r[2] = f;  r[3] = i;
                ETT.push_back(r);
            }
        std::sort(ETT.begin(),ETT.end());

        // count the number of edges (assume manifoldness)
        int En = 1; // the last is always counted
        for(unsigned i=0;i<ETT.size()-1;++i)
            if (!((ETT[i][0] == ETT[i+1][0]) && (ETT[i][1] == ETT[i+1][1])))
                ++En;

        EV = Eigen::MatrixXi::Constant((int)(En),2,-1);
        FE = Eigen::MatrixXi::Constant((int)(F.rows()),3,-1);
        EF = Eigen::MatrixXi::Constant((int)(En),2,-1);
        En = 0;

        for(unsigned i=0;i<ETT.size();++i)
        {
            if (i == ETT.size()-1 ||
                    !((ETT[i][0] == ETT[i+1][0]) && (ETT[i][1] == ETT[i+1][1]))
                    )
            {
                // Border edge
                std::vector<int>& r1 = ETT[i];
                EV(En,0)     = r1[0];
                EV(En,1)     = r1[1];
                EF(En,0)    = r1[2];
                FE(r1[2],r1[3]) = En;
            }
            else
            {
                std::vector<int>& r1 = ETT[i];
                std::vector<int>& r2 = ETT[i+1];
                EV(En,0)     = r1[0];
                EV(En,1)     = r1[1];
                EF(En,0)    = r1[2];
                EF(En,1)    = r2[2];
                FE(r1[2],r1[3]) = En;
                FE(r2[2],r2[3]) = En;
                ++i; // skip the next one
            }
            ++En;
        }

        // Sort the relation EF, accordingly to EV
        // the first one is the face on the left of the edge

        for(unsigned i=0; i<EF.rows(); ++i)
        {
            int fid = EF(i,0);
            bool flip = true;
            // search for edge EV.row(i)
            for (unsigned j=0; j<3; ++j)
            {
                if ((F(fid,j) == EV(i,0)) && (F(fid,(j+1)%3) == EV(i,1)))
                    flip = false;
            }

            if (flip)
            {
                int tmp = EF(i,0);
                EF(i,0) = EF(i,1);
                EF(i,1) = tmp;
            }
        }
    }

public:

    // return mesh as vector of vertices and faces
    static void GetTriMeshData(const MeshType &mesh,
                               Eigen::MatrixXi &faces,
                               MatrixXm &vert)
    {
        tri::RequireCompactness(mesh);
        // create eigen matrix of vertices
        vert=MatrixXm(mesh.VN(), 3);

        // copy vertices
        for (int i = 0; i < mesh.VN(); i++)
            for (int j = 0; j < 3; j++)
                vert(i,j) = mesh.vert[i].cP()[j];

        // create eigen matrix of faces
        faces=Eigen::MatrixXi(mesh.FN(), 3);

        // copy faces
        for (int i = 0; i < mesh.FN(); i++)
            for (int j = 0; j < 3; j++)
                faces(i,j) = (int)tri::Index(mesh,mesh.face[i].cV(j));
    }

    // return normals of the mesh
    static void GetNormalData(const MeshType &mesh,
                              MatrixXm &Nvert,
                              MatrixXm &Nface)
    {
        // create eigen matrix of vertices
        Nvert=MatrixXm(mesh.VN(), 3);
        Nface=MatrixXm(mesh.FN(), 3);

        // per vertices normals
        for (int i = 0; i < mesh.VN(); i++)
            for (int j = 0; j < 3; j++)
                Nvert(i,j) = mesh.vert[i].cN()[j];

        // per vertices normals
        for (int i = 0; i < mesh.FN(); i++)
            for (int j = 0; j < 3; j++)
                Nface(i,j) = mesh.face[i].cN()[j];
    }

    // get face to face adjacency
    static void GetTriFFAdjacency(MeshType &mesh,
                                  Eigen::MatrixXi &FFp,
                                  Eigen::MatrixXi &FFi)
    {
        tri::UpdateTopology<MeshType>::FaceFace(mesh);
        FFp = Eigen::MatrixXi(mesh.FN(),3);
        FFi = Eigen::MatrixXi(mesh.FN(),3);

        for (int i = 0; i < mesh.FN(); i++)
            for (int j = 0; j < 3; j++)
            {
                FaceType *AdjF=mesh.face[i].FFp(j);
                if (AdjF==&mesh.face[i])
                {
                    FFp(i,j)=-1;
                    FFi(i,j)=-1;
                }
                else
                {
                    FFp(i,j)=tri::Index(mesh,AdjF);
                    FFi(i,j)=mesh.face[i].FFi(j);
                }
            }
    }

    // get edge to face and edge to vertex adjacency
    static void GetTriEdgeAdjacency(const MeshType &mesh,
                                    Eigen::MatrixXi& EV,
                                    Eigen::MatrixXi& FE,
                                    Eigen::MatrixXi& EF)
    {
        Eigen::MatrixXi faces;
        MatrixXm vert;
        GetTriMeshData(mesh,faces,vert);
        GetTriEdgeAdjacency(vert,faces,EV,FE,EF);
    }

    static Eigen::Vector3d VectorFromCoord(CoordType v)
    {
        Eigen::Vector3d ret(v[0],v[1],v[2]);
        return ret;
    }

    template< class VecType >
    static void PerVertexArea(MeshType &m, VecType &h)
    {
        tri::RequireCompactness(m);
        h.resize(m.vn);
        for(int i=0;i<m.vn;++i) h[i]=0;
        for(FaceIterator fi=m.face.begin(); fi!=m.face.end();++fi)
        {
            ScalarType a = DoubleArea(*fi)/6.0;
            for(int j=0;j<fi->VN();++j)
                h[tri::Index(m,fi->V(j))] += a;
        }
    }

    template< class VecType >
    static void PerFaceArea(MeshType &m, VecType &h)
    {
        tri::RequireCompactness(m);
        h.resize(m.fn);
        for(int i=0;i<m.fn;++i)
            h[i] =DoubleArea(m.face[i])/2.0;
    }


    static void MassMatrixEntry(MeshType &m,
                                std::vector<std::pair<int,int> > &index,
                                std::vector<ScalarType> &entry,
                                bool vertexCoord=true)
    {
        tri::RequireCompactness(m);

        typename MeshType::template PerVertexAttributeHandle<ScalarType> h =
                tri::Allocator<MeshType>:: template GetPerVertexAttribute<ScalarType>(m, "area");
        for(int i=0;i<m.vn;++i) h[i]=0;

        for(FaceIterator fi=m.face.begin(); fi!=m.face.end();++fi)
        {
            ScalarType a = DoubleArea(*fi);
            for(int j=0;j<fi->VN();++j)
                h[tri::Index(m,fi->V(j))] += a;
        }
        ScalarType maxA=0;
        for(int i=0;i<m.vn;++i)
            maxA = max(maxA,h[i]);

        //store the index and the scalar for the sparse matrix
        for (size_t i=0;i<m.vert.size();i++)
        {
            if (vertexCoord)
            {
                for (size_t j=0;j<3;j++)
                {
                    int currI=(i*3)+j;
                    index.push_back(std::pair<int,int>(currI,currI));
                    entry.push_back(h[i]/maxA);
                }
            }
            else
            {
                int currI=i;
                index.push_back(std::pair<int,int>(currI,currI));
                entry.push_back(h[i]/maxA);
            }
        }
        tri::Allocator<MeshType>::template DeletePerVertexAttribute<ScalarType>(m,h);
    }


    static void GetLaplacianEntry(MeshType &mesh,
                                  FaceType &f,
                                  std::vector<std::pair<int,int> > &index,
                                  std::vector<ScalarType> &entry,
                                  bool cotangent,
                                  ScalarType weight = 1,
                                  bool vertexCoord=true)
    {
        if (cotangent) vcg::tri::MeshAssert<MeshType>::OnlyTriFace(mesh);

        for (int i=0;i<f.VN();i++)
        {

            if (cotangent)
            {
                weight=Harmonic<MeshType>::template CotangentWeight<ScalarType>(f,i);
            }

            //get the index of the vertices
            int indexV0=Index(mesh,f.V0(i));
            int indexV1=Index(mesh,f.V1(i));

            if (vertexCoord)
            {
                //then assemble the matrix
                for (int j=0;j<3;j++)
                {
                    //multiply by 3 and add the component
                    int currI0=(indexV0*3)+j;
                    int currI1=(indexV1*3)+j;

                    index.push_back(std::pair<int,int>(currI0,currI0));
                    entry.push_back(weight);
                    index.push_back(std::pair<int,int>(currI0,currI1));
                    entry.push_back(-weight);

                    index.push_back(std::pair<int,int>(currI1,currI1));
                    entry.push_back(weight);
                    index.push_back(std::pair<int,int>(currI1,currI0));
                    entry.push_back(-weight);
                }
            }
            else
            {
                int currI0=(indexV0);
                int currI1=(indexV1);

                index.push_back(std::pair<int,int>(currI0,currI0));
                entry.push_back(weight);
                index.push_back(std::pair<int,int>(currI0,currI1));
                entry.push_back(-weight);

                index.push_back(std::pair<int,int>(currI1,currI1));
                entry.push_back(weight);
                index.push_back(std::pair<int,int>(currI1,currI0));
                entry.push_back(-weight);
            }
        }
    }


    static void GetLaplacianMatrix(MeshType &mesh,
                                   std::vector<std::pair<int,int> > &index,
                                   std::vector<ScalarType> &entry,
                                   bool cotangent,
                                   ScalarType weight = 1,
                                   bool vertexCoord=true )
    {
        //store the index and the scalar for the sparse matrix
        for (size_t i=0;i<mesh.face.size();i++)
            GetLaplacianEntry(mesh,mesh.face[i],index,entry,cotangent,weight,vertexCoord);
    }



};

} // end namespace tri
} // end namespace vcg
#endif // MESH_TO_MATRIX_CONVERTER
