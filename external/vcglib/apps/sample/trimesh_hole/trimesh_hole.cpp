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

#include<vcg/space/triangle3.h>
#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/hole.h>
#include<vcg/complex/algorithms/local_optimization.h>
#include<vcg/complex/algorithms/local_optimization/tri_edge_flip.h>
#include<vcg/complex/algorithms/smooth.h>
#include<vcg/complex/algorithms/refine.h>

#include<vcg/complex/algorithms/update/selection.h>

// topology computation
#include<vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/normal.h>

// half edge iterators
#include<vcg/simplex/face/pos.h>

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

using namespace vcg;
using namespace std;

class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>		::AsVertexType,
    Use<MyFace>			::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags, vertex::Normal3f, vertex::Mark, vertex::Color4b >{};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef,face::FFAdj, face::Mark, face::BitFlags, face::Normal3f> {};

class MyMesh : public tri::TriMesh< vector<MyVertex>, vector<MyFace > >{};

//Delaunay
class MyDelaunayFlip: public vcg::tri::TriEdgeFlip< MyMesh, MyDelaunayFlip > {
public:
    typedef  vcg::tri::TriEdgeFlip< MyMesh,  MyDelaunayFlip > TEF;
    inline MyDelaunayFlip(  const TEF::PosType &p, int i,BaseParameterClass *pp) :TEF(p,i,pp){}
};

bool callback(int percent, const char *str) {
  cout << "str: " << str << " " << percent << "%\r";
  return true;
}

template <class MESH>
bool NormalTest(typename face::Pos<typename MESH::FaceType> pos)
{
    //giro intorno al vertice e controllo le normali
    typename MESH::ScalarType thr = 0.0f;
        typename MESH::CoordType NdP = vcg::TriangleNormal<typename MESH::FaceType>(*pos.f);
    typename MESH::CoordType tmp, oop, soglia = typename MESH::CoordType(thr,thr,thr);
    face::Pos<typename MESH::FaceType> aux=pos;
    do{
        aux.FlipF();
        aux.FlipE();
                oop = Abs(tmp - ::vcg::TriangleNormal<typename MESH::FaceType>(*pos.f));
        if(oop < soglia )return false;
    }while(aux != pos && !aux.IsBorder());

    return true;
}

int main(int argc,char ** argv){

    if(argc<5)
    {
        printf(
            "\n     HoleFilling (" __DATE__ ")\n"
            "Visual Computing Group I.S.T.I. C.N.R.\n"
      "Usage: trimesh_hole #algorithm #size filein.ply fileout.ply \n"
            "#algorithm: \n"
            " 1) Trivial Ear \n"
            " 2) Minimum weight Ear \n"
            " 3) Selfintersection Ear \n"
            " 4) Minimum weight \n"
            );
        exit(0);
    }

    int algorithm = atoi(argv[1]);
    int holeSize  = atoi(argv[2]);
    if(algorithm < 0 && algorithm > 4)
    {
    printf("Error in algorithm's selection %i\n",algorithm);
        exit(0);
    }

    MyMesh m;

    if(tri::io::ImporterPLY<MyMesh>::Open(m,argv[3])!=0)
    {
        printf("Error reading file  %s\n",argv[2]);
        exit(0);
    }

    //update the face-face topology
    tri::UpdateTopology<MyMesh>::FaceFace(m);
    tri::UpdateNormal<MyMesh>::PerVertexPerFace(m);
    tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
  assert(tri::Clean<MyMesh>::IsFFAdjacencyConsistent(m));

    //compute the average of face area
    float AVG,sumA=0.0f;
    int numA=0,indice;
    indice = m.face.size();
    MyMesh::FaceIterator fi;
    for(fi=m.face.begin();fi!=m.face.end();++fi)
    {
            sumA += DoubleArea(*fi)/2;
            numA++;
            for(int ind =0;ind<3;++ind)
                fi->V(ind)->InitIMark();
    }
    AVG=sumA/numA;

  //tri::Hole<MyMesh> holeFiller;
    switch(algorithm)
    {
  case 1:			tri::Hole<MyMesh>::EarCuttingFill<tri::TrivialEar<MyMesh> >(m,holeSize,false);                	        break;
  case 2:   	tri::Hole<MyMesh>::EarCuttingFill<tri::MinimumWeightEar< MyMesh> >(m,holeSize,false,callback);          break;
  case 3: 		tri::Hole<MyMesh>::EarCuttingIntersectionFill<tri::SelfIntersectionEar< MyMesh> >(m,holeSize,false);		break;
  case 4: 		tri::Hole<MyMesh>::MinimumWeightFill(m,holeSize, false); tri::UpdateTopology<MyMesh>::FaceFace(m);      break;
    }

    tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);

  assert(tri::Clean<MyMesh>::IsFFAdjacencyConsistent(m));

  printf("\nStart refinig...\n");

/*start refining */
    MyMesh::VertexIterator vi;
    MyMesh::FaceIterator f;
    std::vector<MyMesh::FacePointer> vf;
    f =  m.face.begin();
    f += indice;
    for(; f != m.face.end();++f)
    {
        if(!f->IsD())
        {
            f->SetS();
        }
    }

    std::vector<MyMesh::FacePointer *> FPP;
    std::vector<MyMesh::FacePointer> added;
    std::vector<MyMesh::FacePointer>::iterator vfit;
    int i=1;
    printf("\n");

    for(f =  m.face.begin();f!=m.face.end();++f) if(!(*f).IsD())
    {
        if( f->IsS() )
        {
            f->V(0)->IsW();
            f->V(1)->IsW();
            f->V(2)->IsW();
        }
        else
        {
            f->V(0)->ClearW();
            f->V(1)->ClearW();
            f->V(2)->ClearW();
        }
    }
    BaseParameterClass pp;
                vcg::LocalOptimization<MyMesh> Fs(m,&pp);
                Fs.SetTargetMetric(0.0f);
                Fs.Init<MyDelaunayFlip >();
                Fs.DoOptimization();


    do
    {
        vf.clear();
        f =  m.face.begin();
        f += indice;
        for(; f != m.face.end();++f)
        {
            if(f->IsS())
            {
                bool test= true;
                for(int ind =0;ind<3;++ind)
                    f->V(ind)->InitIMark();
                test = (DoubleArea<MyMesh::FaceType>(*f)/2) > AVG;
                if(test)
                {
                    vf.push_back(&(*f));
                }
            }
        }

        //info print
    printf("\r Refining [%d] - > %d",i,int(vf.size()));
        i++;

        FPP.clear();
        added.clear();

        for(vfit = vf.begin(); vfit!=vf.end();++vfit)
        {
            FPP.push_back(&(*vfit));
        }
        int toadd= vf.size();
        MyMesh::FaceIterator f1,f2;
        f2 = tri::Allocator<MyMesh>::AddFaces(m,(toadd*2),FPP);
        MyMesh::VertexIterator vertp = tri::Allocator<MyMesh>::AddVertices(m,toadd);
        std::vector<MyMesh::FacePointer> added;
        added.reserve(toadd);
        vfit=vf.begin();

        for(int i = 0; i<toadd;++i,f2++,vertp++)
        {
            f1=f2;
            f2++;
            TriSplit<MyMesh,CenterPointBarycenter<MyMesh> >::Apply(vf[i],&(*f1),&(*f2),&(*vertp),CenterPointBarycenter<MyMesh>() );
            f1->SetS();
            f2->SetS();
            for(int itr=0;itr<3;itr++)
            {
                f1->V(itr)->SetW();
                f2->V(itr)->SetW();
            }
            added.push_back( &(*f1) );
            added.push_back( &(*f2) );
        }

        BaseParameterClass pp;
        vcg::LocalOptimization<MyMesh> FlippingSession(m,&pp);
        FlippingSession.SetTargetMetric(0.0f);
        FlippingSession.Init<MyDelaunayFlip >();
        FlippingSession.DoOptimization();

    }while(!vf.empty());

    vcg::LocalOptimization<MyMesh> Fiss(m,&pp);
    Fiss.SetTargetMetric(0.0f);
    Fiss.Init<MyDelaunayFlip >();
    Fiss.DoOptimization();

/*end refining */

    tri::io::ExporterPLY<MyMesh>::Save(m,"PreSmooth.ply",false);

    int UBIT = MyMesh::VertexType::NewBitFlag();
    f =  m.face.begin();
    f += indice;
    for(; f != m.face.end();++f)
    {
        if(f->IsS())
        {
            for(int ind =0;ind<3;++ind){
                if(NormalTest<MyMesh>(face::Pos<MyMesh::FaceType>(&(*f),ind )))
                {
                    f->V(ind)->SetUserBit(UBIT);
                }
            }
            f->ClearS();
        }
    }

    for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD())
    {
        if( vi->IsUserBit(UBIT) )
        {
            (*vi).SetS();
            vi->ClearUserBit(UBIT);
        }
    }

    tri::Smooth<MyMesh>::VertexCoordLaplacian(m,1,true);

    printf("\nCompleted. Saving....\n");

  tri::io::ExporterPLY<MyMesh>::Save(m,argv[4],false);
    return 0;
}

