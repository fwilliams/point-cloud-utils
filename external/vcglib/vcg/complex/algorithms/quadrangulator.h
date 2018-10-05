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
#ifndef MIQ_QUADRANGULATOR_H
#define MIQ_QUADRANGULATOR_H

#include <vcg/complex/complex.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/jumping_pos.h>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/complex/algorithms/update/texture.h>

#define precisionQ 0.0000000001

namespace vcg {
namespace tri {
template <class TriMesh,class PolyMesh>
class Quadrangulator
{

public:
    typedef typename TriMesh::FaceType TriFaceType;
    typedef typename TriMesh::VertexType TriVertexType;
    typedef typename TriMesh::CoordType CoordType;
    typedef typename TriMesh::ScalarType ScalarType;

    typedef typename PolyMesh::FaceType PolyFaceType;
    typedef typename PolyMesh::VertexType PolyVertexType;
    typedef typename PolyMesh::CoordType PolyCoordType;
    typedef typename PolyMesh::ScalarType PolyScalarType;


    struct InterpolationInfo
    {
        CoordType Pos3D;
        vcg::Point2<ScalarType> PosUV;
        ScalarType alpha;
        bool to_split;

        InterpolationInfo()
        {
            Pos3D=CoordType(0,0,0);
            PosUV=vcg::Point2<ScalarType>(0,0);
            to_split=false;
            alpha=-1;
        }
    };

    //the interpolation map that is saved once to be univoque per edge
    typedef std::pair<CoordType,CoordType > KeyEdgeType;

    std::map<KeyEdgeType,InterpolationInfo> InterpMap;

    //ScalarType UVtolerance;

private:

    bool ToSplit(const vcg::Point2<ScalarType> &uv0,
                 const vcg::Point2<ScalarType> &uv1,
                 int Dir,
                 ScalarType &alpha)
    {
        ScalarType val0=uv0.V(Dir);
        ScalarType val1=uv1.V(Dir);
        int IntegerLine0=floor(val0);
        int IntegerLine1=floor(val1);
        if (IntegerLine0==IntegerLine1)
            return false;//no integer line pass throught the edge

        bool swapped=false;
        if (IntegerLine0>IntegerLine1)
        {
            std::swap(IntegerLine0,IntegerLine1);
            std::swap(val0,val1);
            assert(val1>=val0);
            swapped=true;
        }

        //then get the first if extist that overcome the threshold
        int IntegerSplit=IntegerLine0+1;
        bool found=false;
        ScalarType dist1,dist0;
        for (int i=IntegerSplit;i<=IntegerLine1;i++)
        {
            dist1=fabs(val1-IntegerSplit);
            dist0=fabs(val0-IntegerSplit);

//            if ((dist0>=UVtolerance)&&
//                (dist1>=UVtolerance))
            if ((val0!=IntegerSplit)&&
                (val1!=IntegerSplit))
            {
                found=true;
                break;
            }
            IntegerSplit++;
        }
        if (!found)return false;

        //have to check distance also in opposite direction
        ScalarType lenght=val1-val0;
        assert(lenght>=0);
        //alpha=1.0-(dist/lenght);
        alpha=(dist1/lenght);
        if (swapped)alpha=1-alpha;
        assert((alpha>0)&&(alpha<1));
        return true;
    }

    void RoundInitial(TriMesh &to_split)
    {
        ScalarType minTolerance=precisionQ;
        //first add all eddge
        for (int i=0;i<to_split.face.size();i++)
        {
            TriFaceType *f=&to_split.face[i];
            for (int j =0;j<3;j++)
            {
                vcg::Point2<ScalarType> UV=f->WT(j).P();

                int int0=floor(UV.X()+0.5);
                int int1=floor(UV.Y()+0.5);

                ScalarType diff0=(fabs(UV.X()-(ScalarType)int0));
                ScalarType diff1=(fabs(UV.Y()-(ScalarType)int1));

                if (diff0<minTolerance)
                    UV.X()=(ScalarType)int0;
                if (diff1<minTolerance)
                    UV.Y()=(ScalarType)int1;

                f->WT(j).P()=UV;
            }
        }
    }

    void RoundSplits(TriMesh &to_split,int dir)
    {
        ScalarType minTolerance=precisionQ;
        //first add all eddge
        for (size_t i=0;i<to_split.face.size();i++)
        {
            TriFaceType *f=&to_split.face[i];
            for (int j =0;j<3;j++)
            {
                CoordType p0=f->P0(j);
                CoordType p1=f->P1(j);
                KeyEdgeType k(p0,p1);
                assert(InterpMap.count(k)==1);
                if (!InterpMap[k].to_split)continue;
                //then get the intepolated value
                vcg::Point2<ScalarType> UV=InterpMap[k].PosUV;

                int int0=floor(UV.X()+0.5);
                int int1=floor(UV.Y()+0.5);

                ScalarType diff0=(fabs(UV.X()-(ScalarType)int0));
                ScalarType diff1=(fabs(UV.Y()-(ScalarType)int1));

                if (diff0<minTolerance)
                    UV.X()=(ScalarType)int0;
                if (diff1<minTolerance)
                    UV.Y()=(ScalarType)int1;

                InterpMap[k].PosUV=UV;
            }
        }
    }

    void InitSplitMap(TriMesh &to_split,
                      int dir)
    {
       assert((dir==0)||(dir==1));
       InterpMap.clear();
       //printf("direction %d\n",dir );
       //first add all eddge
       for (int i=0;i<to_split.face.size();i++)
       {
           TriFaceType *f=&to_split.face[i];
           for (int j =0;j<3;j++)
           {
               CoordType p0=f->P0(j);
               CoordType p1=f->P1(j);
               vcg::Point2<ScalarType> Uv0=f->V0(j)->T().P();
               vcg::Point2<ScalarType> Uv1=f->V1(j)->T().P();
               KeyEdgeType k(p0,p1);
//               printf("p0 (%5.5f,%5.5f,%5.5f) p1(%5.5f,%5.5f,%5.5f) \n",p0.X(),p0.Y(),p0.Z(),p1.X(),p1.Y(),p1.Z());
//               printf("uv0 (%5.5f,%5.5f) uv1(%5.5f,%5.5f) \n",Uv0.X(),Uv0.Y(),Uv1.X(),Uv1.Y());
//               fflush(stdout);
               assert(InterpMap.count(k)==0);
               InterpMap[k]=InterpolationInfo();
           }
       }

       //then set the ones to be splitted
       for (size_t i=0;i<to_split.face.size();i++)
       {
           TriFaceType *f=&to_split.face[i];
           for (int j =0;j<3;j++)
           {
               CoordType p0=f->P0(j);
               CoordType p1=f->P1(j);
               vcg::Point2<ScalarType> uv0=f->V0(j)->T().P();
               vcg::Point2<ScalarType> uv1=f->V1(j)->T().P();

               ScalarType alpha;
               if (!ToSplit(uv0,uv1,dir,alpha))continue;

               KeyEdgeType k(p0,p1);
               assert(InterpMap.count(k)==1);
               InterpMap[k].Pos3D=p0*alpha+p1*(1-alpha);
               InterpMap[k].PosUV=uv0*alpha+uv1*(1-alpha);
               InterpMap[k].to_split=true;
               InterpMap[k].alpha=alpha;
           }
       }

       //then make them coherent
       for (size_t i=0;i<to_split.face.size();i++)
       {
           TriFaceType *f=&to_split.face[i];
           for (int j =0;j<3;j++)
           {
               CoordType p0=f->P0(j);
               CoordType p1=f->P1(j);
               vcg::Point2<ScalarType> uv0=f->V0(j)->T().P();
               vcg::Point2<ScalarType> uv1=f->V1(j)->T().P();
//               if (p0>p1)continue; //only one verse of coherence

               KeyEdgeType k0(p0,p1);
               assert(InterpMap.count(k0)==1);//there should be already in the
                                              //table and it should be coherent

               KeyEdgeType k1(p1,p0);
               if(InterpMap.count(k1)==0)continue;//REAL border, no need for update

               bool to_split0=InterpMap[k0].to_split;
               bool to_split1=InterpMap[k1].to_split;

               //the find all possible cases
               if ((!to_split0)&&(!to_split1))continue;

               if ((to_split0)&&(to_split1))
               {
                   CoordType Pos3D=InterpMap[k1].Pos3D;
                   InterpMap[k0].Pos3D=Pos3D;

                   //check if need to make coherent also the UV Position
                   //skip the fake border and do the rest
                   bool IsBorderFF=(f->FFp(j)==f);

                   if (!IsBorderFF) //in this case they should have same UVs
                       InterpMap[k0].PosUV=InterpMap[k1].PosUV;
                   else
                   {
                       ScalarType alpha=InterpMap[k1].alpha;
                       assert((alpha>=0)&&(alpha<=1));
                       alpha=1-alpha;
                       InterpMap[k0].PosUV=alpha*uv0+(1-alpha)*uv1;
                       InterpMap[k0].alpha=alpha;
                   }

               }
               else
               if ((!to_split0)&&(to_split1))
               {
                   CoordType Pos3D=InterpMap[k1].Pos3D;
                   InterpMap[k0].Pos3D=Pos3D;

                   //check if need to make coherent also the UV Position
                   //skip the fake border and do the rest
                   bool IsBorderFF=(f->FFp(j)==f);

                   InterpMap[k0].to_split=true;

                   if (!IsBorderFF) //in this case they should have same UVs
                       InterpMap[k0].PosUV=InterpMap[k1].PosUV;
                   else //recalculate , it pass across a seam
                   {
                       ScalarType alpha=InterpMap[k1].alpha;
                       assert((alpha>=0)&&(alpha<=1));
                       alpha=1-alpha;
                       InterpMap[k0].PosUV=alpha*uv0+(1-alpha)*uv1;
                       InterpMap[k0].alpha=alpha;
                   }
              }
           }
       }

       RoundSplits(to_split,dir);
    }

    // Basic subdivision class
    // This class must provide methods for finding the position of the newly created vertices
    // In this implemenation we simply put the new vertex in the MidPoint position.
    // Color and TexCoords are interpolated accordingly.
    template<class MESH_TYPE>
    struct SplitMidPoint : public   std::unary_function<vcg::face::Pos<typename MESH_TYPE::FaceType> ,  typename MESH_TYPE::CoordType >
    {
        typedef typename MESH_TYPE::VertexType VertexType;
        typedef typename MESH_TYPE::FaceType FaceType;
        typedef typename MESH_TYPE::CoordType CoordType;

        std::map<KeyEdgeType,InterpolationInfo> *MapEdge;

        void operator()(typename MESH_TYPE::VertexType &nv,
                        vcg::face::Pos<typename MESH_TYPE::FaceType>  ep)
        {
            VertexType* v0=ep.f->V0(ep.z);
            VertexType* v1=ep.f->V1(ep.z);
            assert(v0!=v1);

            CoordType p0=v0->P();
            CoordType p1=v1->P();
            assert(p0!=p1);

            KeyEdgeType k(p0,p1);
            bool found=(MapEdge->count(k)==1);
            assert(found);
            bool to_split=(*MapEdge)[k].to_split;
            assert(to_split);

            //get the value on which the edge must be splitted
            nv.P()= (*MapEdge)[k].Pos3D;
            //nv.N()= v0->N()*alpha+v1->N()*(1.0-alpha);
            nv.T().P()=(*MapEdge)[k].PosUV;
        }

        vcg::TexCoord2<ScalarType> WedgeInterp(vcg::TexCoord2<ScalarType> &t0,
                                               vcg::TexCoord2<ScalarType> &t1)
        {
            return (vcg::TexCoord2<ScalarType>(0,0));
        }

        SplitMidPoint(std::map<KeyEdgeType,InterpolationInfo> *_MapEdge){MapEdge=_MapEdge;}
    };

    template <class MESH_TYPE>
    class EdgePredicate
    {
        typedef typename MESH_TYPE::VertexType VertexType;
        typedef typename MESH_TYPE::FaceType FaceType;
        typedef typename MESH_TYPE::ScalarType ScalarType;

        std::map<KeyEdgeType,InterpolationInfo> *MapEdge;

    public:

        bool operator()(vcg::face::Pos<typename MESH_TYPE::FaceType> ep) const
        {
            VertexType* v0=ep.f->V0(ep.z);
            VertexType* v1=ep.f->V1(ep.z);
            assert(v0!=v1);

            CoordType p0=v0->P();
            CoordType p1=v1->P();
            assert(p0!=p1);

            KeyEdgeType k(p0,p1);
            bool found=(MapEdge->count(k)==1);
            assert(found);
            bool to_split=(*MapEdge)[k].to_split;
            return(to_split);
        }

        EdgePredicate(std::map<KeyEdgeType,InterpolationInfo> *_MapEdge){MapEdge=_MapEdge;}
    };

    void SplitTrisDir(TriMesh &to_split,
                       int dir)
    {
        bool done=true;
        //int step=0;
        while (done)
        {
            printf("Number of Vertices %d \n",to_split.vn);
            fflush(stdout);

            InitSplitMap(to_split,dir);

            SplitMidPoint<TriMesh> splMd(&InterpMap);
            EdgePredicate<TriMesh> eP(&InterpMap);

            done=vcg::tri::RefineE<TriMesh,SplitMidPoint<TriMesh>,EdgePredicate<TriMesh> >(to_split,splMd,eP);

        }
        printf("Number of Vertices %d \n",to_split.vn);
        fflush(stdout);
         fflush(stdout);
    }


    bool IsOnIntegerLine(vcg::Point2<ScalarType> uv0,
                         vcg::Point2<ScalarType> uv1)
    {
        for (int dir=0;dir<2;dir++)
        {
            ScalarType val0=uv0.V(dir);
            ScalarType val1=uv1.V(dir);
            int integer0=floor(uv0.V(dir)+0.5);
            int integer1=floor(uv1.V(dir)+0.5);
            if (integer0!=integer1)continue;
//            if ((fabs(val0-(ScalarType)integer0))>=UVtolerance)continue;
//            if ((fabs(val1-(ScalarType)integer1))>=UVtolerance)continue;
            if (val0!=(ScalarType)floor(val0))continue;
            if (val1!=(ScalarType)floor(val1))continue;
            return true;
        }
        return false;
    }

    bool IsOnIntegerVertex(vcg::Point2<ScalarType> uv,
                           bool IsB)
    {
        int onIntegerL=0;
        for (int dir=0;dir<2;dir++)
        {
            ScalarType val0=uv.V(dir);
            int integer0=floor(val0+0.5);
            //if ((fabs(val0-(ScalarType)integer0))<UVtolerance)onIntegerL++;
            if (val0==(ScalarType)floor(val0))onIntegerL++;
        }
        if ((IsB)&&(onIntegerL>0))return true;
        return (onIntegerL==2);
    }


    void InitIntegerEdgesVert(TriMesh &Tmesh)
    {
        //IntegerEdges.clear();
        vcg::tri::UpdateFlags<TriMesh>::FaceSetF(Tmesh);
        vcg::tri::UpdateFlags<TriMesh>::FaceClearS(Tmesh);
        vcg::tri::UpdateFlags<TriMesh>::VertexClearS(Tmesh);

        for (unsigned int i=0;i<Tmesh.face.size();i++)
        {
            TriFaceType *f=&Tmesh.face[i];
            if (f->IsD())continue;
            for (int j=0;j<3;j++)
            {
                bool IsBorder=f->IsB(j);
                if (IsBorder)
                    f->ClearF(j);
                else
                {
                    vcg::Point2<ScalarType> uv0=f->WT(j).P();
                    vcg::Point2<ScalarType> uv1=f->WT((j+1)%3).P();

                    if (IsOnIntegerLine(uv0,uv1))
                    {
                        f->ClearF(j);
                        TriFaceType *f1=f->FFp(j);
                        int z=f->FFi(j);
                        assert(f1!=f);
                        f1->ClearF(z);
                    }
                }

                bool BorderV=f->V(j)->IsB();

                if (IsOnIntegerVertex(f->WT(j).P(),BorderV))
                    f->V(j)->SetS();
            }
        }
    }

    short int AlignmentEdge(TriFaceType *f,
                      int edge_index)
    {
        vcg::Point2<ScalarType> uv0=f->WT(edge_index).P();
        vcg::Point2<ScalarType> uv1=f->WT((edge_index+1)%3).P();
        if (uv0.X()==uv1.X())return 0;
        if (uv0.Y()==uv1.Y())return 1;
        return -1;
    }

    void FindPolygon(vcg::face::Pos<TriFaceType> &currPos,
                     std::vector<TriVertexType *> &poly,
                     std::vector<short int> &UVpoly)
    {
        currPos.F()->SetV();
        currPos.F()->C()=vcg::Color4b(255,0,0,255);
        poly.clear();
        assert(currPos.V()->IsS());
        TriVertexType *v_init=currPos.V();
        poly.push_back(currPos.V());

        //retrieve UV
        int indexV0=currPos.E();

        short int Align=AlignmentEdge(currPos.F(),currPos.E());

        std::vector<short int> TempUVpoly;
        TempUVpoly.push_back(Align);

        do
        {
            currPos.NextNotFaux();
            currPos.F()->SetV();
            currPos.F()->C()=vcg::Color4b(255,0,0,255);

            if ((currPos.V()->IsS())&&(currPos.V()!=v_init))
            {
                poly.push_back(currPos.V());

                short int Align=AlignmentEdge(currPos.F(),currPos.E());

                TempUVpoly.push_back(Align);
            }

        }while (currPos.V()!=v_init);

        //then shift the order of UV by one
        //to be consistent with edge ordering
        int size=TempUVpoly.size();
        for (int i=0;i<size;i++)
            UVpoly.push_back(TempUVpoly[(i+1)%size]);
    }

    void FindPolygons(TriMesh &Tmesh,
                      std::vector<std::vector<TriVertexType *> > &polygons,
                      std::vector<std::vector<short int> > &UV)
    {
        vcg::tri::UpdateFlags<TriMesh>::FaceClearV(Tmesh);
        for (unsigned int i=0;i<Tmesh.face.size();i++)
        {
            TriFaceType * f=&Tmesh.face[i];
            if (f->IsV())continue;

            for (int j=0;j<3;j++)
            {
                TriVertexType* v0=f->V0(j);
                if (!v0->IsS())continue;
                if (f->IsF(j))continue;

                vcg::face::Pos<TriFaceType> startPos(f,j);

                std::vector<TriVertexType *> poly;
                std::vector< short int> UVpoly;

                FindPolygon(startPos,poly,UVpoly);

                if (poly.size()>2)
                {
                   assert(poly.size()==UVpoly.size());
                   std::reverse(poly.begin(),poly.end());
//                    std::reverse(UVpoly.begin(),UVpoly.end());

                    polygons.push_back(poly);
                    UV.push_back(UVpoly);

                }
                //only one polygon per initial face
                break;
            }
        }

    }

    //FUNCTIONS NEEDED BY "UV WEDGE TO VERTEX" FILTER
    static void ExtractVertex(const TriMesh & srcMesh,
                              const TriFaceType & f,
                              int whichWedge,
                              const TriMesh & dstMesh,
                              TriVertexType & v)
    {
        (void)srcMesh;
        (void)dstMesh;
        // This is done to preserve every single perVertex property
        // perVextex Texture Coordinate is instead obtained from perWedge one.
        v.ImportData(*f.cV(whichWedge));
        v.T() = f.cWT(whichWedge);
    }

    static bool CompareVertex(const TriMesh & m,
                              TriVertexType & vA,
                              TriVertexType & vB)
    {
        (void)m;
        return (vA.cT() == vB.cT());
    }

    void ConvertWTtoVT(TriMesh &Tmesh)
    {
        int vn = Tmesh.vn;
        vcg::tri::AttributeSeam::SplitVertex(Tmesh, ExtractVertex, CompareVertex);
        vcg::tri::UpdateTopology<TriMesh>::FaceFace(Tmesh);
       // vcg::tri::UpdateFlags<TriMesh>::FaceBorderFromFF(Tmesh);
    }

    void ConvertVTtoWT(TriMesh &Tmesh)
    {
        vcg::tri::UpdateTexture<TriMesh>::WedgeTexFromVertexTex(Tmesh);
        vcg::tri::Clean<TriMesh>::RemoveDuplicateVertex(Tmesh);
    }

    void ReupdateMesh(TriMesh &Tmesh)
    {
        vcg::tri::UpdateNormal<TriMesh>::PerFaceNormalized(Tmesh);	 // update Normals
        vcg::tri::UpdateNormal<TriMesh>::PerVertexNormalized(Tmesh);// update Normals
        //compact the mesh
        vcg::tri::Allocator<TriMesh>::CompactVertexVector(Tmesh);
        vcg::tri::Allocator<TriMesh>::CompactFaceVector(Tmesh);
        vcg::tri::UpdateTopology<TriMesh>::FaceFace(Tmesh);			 // update Topology
        vcg::tri::UpdateTopology<TriMesh>::TestFaceFace(Tmesh);		 //and test it
        //set flags
        vcg::tri::UpdateFlags<TriMesh>::VertexClearV(Tmesh);
        vcg::tri::UpdateFlags<TriMesh>::FaceBorderFromFF(Tmesh);
        vcg::tri::UpdateFlags<TriMesh>::VertexBorderFromFaceBorder(Tmesh);
    }

public:


    void TestIsProper(TriMesh &Tmesh)
    {


        //test manifoldness
        int test=vcg::tri::Clean<TriMesh>::CountNonManifoldVertexFF(Tmesh);
        //assert(test==0);
        if (test != 0)
            cerr << "Assertion failed: TestIsProper NonManifoldVertices!" << endl;

        test=vcg::tri::Clean<TriMesh>::CountNonManifoldEdgeFF(Tmesh);
        //assert(test==0);
        if (test != 0)
            cerr << "Assertion failed: TestIsProper NonManifoldEdges" << endl;

        for (unsigned int i=0;i<Tmesh.face.size();i++)
        {
            TriFaceType *f=&Tmesh.face[i];
            assert (!f->IsD());
            for (int z=0;z<3;z++)
            {
                //int indexOpp=f->FFi(z);
                TriFaceType *Fopp=f->FFp(z);
                if (Fopp==f) continue;
                //assert( f->FFp(z)->FFp(f->FFi(z))==f );
                if (f->FFp(z)->FFp(f->FFi(z))!=f)
                    cerr << "Assertion failed: TestIsProper f->FFp(z)->FFp(f->FFi(z))!=f " << endl;
            }
        }
    }



    void Quadrangulate(TriMesh &Tmesh,
                       PolyMesh &Pmesh,
                       std::vector< std::vector< short int> > &UV)
    {
        UV.clear();
        Pmesh.Clear();

        vcg::tri::UpdateTopology<TriMesh>::FaceFace(Tmesh);

        TestIsProper(Tmesh);

        RoundInitial(Tmesh);

        //UVtolerance=tolerance;

        //split to per vert
        ConvertWTtoVT(Tmesh);


        vcg::tri::Allocator<TriMesh>::CompactVertexVector(Tmesh);
        vcg::tri::Allocator<TriMesh>::CompactFaceVector(Tmesh);
        vcg::tri::UpdateTopology<TriMesh>::FaceFace(Tmesh);
        (void)Pmesh;
        //TestIsProper(Tmesh);

        //then split the tris along X
        SplitTrisDir(Tmesh,0);
        SplitTrisDir(Tmesh,1);

        //merge back the mesh and WT coords
        ConvertVTtoWT(Tmesh);

        //CleanMesh(Pmesh);

        //update properties of the mesh
        ReupdateMesh(Tmesh);

        //test manifoldness
        TestIsProper(Tmesh);

        InitIntegerEdgesVert(Tmesh);

        for (int i=0;i<Tmesh.face.size();i++)
            Tmesh.face[i].C()=vcg::Color4b(255,255,255,255);

        std::vector<std::vector<TriVertexType *> > polygons;
        FindPolygons(Tmesh,polygons,UV);

        //then add to the polygonal mesh
        Pmesh.Clear();

        int numV=vcg::tri::UpdateSelection<TriMesh>::VertexCount(Tmesh);

        //first create vertices
        vcg::tri::Allocator<PolyMesh>::AddVertices(Pmesh,numV);

        std::map<CoordType,int> VertMap;
        int index=0;
        for(unsigned int i=0;i<Tmesh.vert.size();i++)
        {
            if (!Tmesh.vert[i].IsS())continue;

            CoordType pos=Tmesh.vert[i].P();
            CoordType norm=Tmesh.vert[i].N();
            vcg::Point2<ScalarType> UV=Tmesh.vert[i].T().P();
            Pmesh.vert[index].P()=typename PolyMesh::CoordType(pos.X(),pos.Y(),pos.Z());
            Pmesh.vert[index].N()=typename PolyMesh::CoordType(norm.X(),norm.Y(),norm.Z());
            Pmesh.vert[index].T().P()=UV;
            VertMap[pos]=index;
            index++;
        }

        //then add polygonal mesh
        vcg::tri::Allocator<PolyMesh>::AddFaces(Pmesh,polygons.size());
        for (unsigned int i=0;i<polygons.size();i++)
        {
            int size=polygons[i].size();
            Pmesh.face[i].Alloc(size);
            for (int j=0;j<size;j++)
            {
                CoordType pos=(polygons[i][j])->P();
                assert(VertMap.count(pos)==1);
                int index=VertMap[pos];
                Pmesh.face[i].V(j)=&(Pmesh.vert[index]);
            }
        }


    }
};
}
}
#endif
