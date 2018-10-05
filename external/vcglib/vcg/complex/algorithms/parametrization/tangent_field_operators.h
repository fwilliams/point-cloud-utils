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
#include <vcg/math/histogram.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <algorithm>

#ifndef VCG_TANGENT_FIELD_OPERATORS
#define VCG_TANGENT_FIELD_OPERATORS

namespace vcg {
namespace tri{




template < typename ScalarType >
vcg::Point2<ScalarType> InterpolateNRosy2D(const std::vector<vcg::Point2<ScalarType> > &V,
                                           const std::vector<ScalarType>  &W,
                                           const int N)
{
    // check parameter
    assert(V.size() == W.size());
    assert( N > 0);

    // create a vector of angles
    std::vector<ScalarType> angles(V.size(), 0);

    // make angle mod 2pi/N by multiplying times N
    for (size_t i = 0; i < V.size(); i++)
        angles[i] = std::atan2(V[i].Y(), V[i].X() ) * N;

    // create vector of directions
    std::vector<vcg::Point2<ScalarType> > VV(V.size(), vcg::Point2<ScalarType>(0,0));

    // compute directions
    for (size_t i = 0; i < V.size(); i++) {
        VV[i].X() = std::cos(angles[i]);
        VV[i].Y() = std::sin(angles[i]);
    }

    // average vector
    vcg::Point2<ScalarType> Res(0,0);

    // compute average of the unit vectors
    ScalarType Sum=0;
    for (size_t i = 0; i < VV.size(); i++)
    {
        Res += VV[i] * W[i];
        Sum+=W[i];
    }
    assert(Sum>0);
    Res /=Sum;

    //R /= VV.rows();

    // scale them back
    ScalarType a = std::atan2(Res.Y(), Res.X()) / N;
    Res.X() = std::cos(a);
    Res.Y() = std::sin(a);

    return Res;
}

template < typename ScalarType >
vcg::Point3<ScalarType> InterpolateNRosy3D(const std::vector<vcg::Point3<ScalarType> > &V,
                                           const std::vector<vcg::Point3<ScalarType> > &Norm,
                                           const std::vector<ScalarType>  &W,
                                           const int N,
                                           const vcg::Point3<ScalarType> &TargetN)
{
    typedef typename vcg::Point3<ScalarType> CoordType;
    ///create a reference frame along TargetN
    CoordType TargetZ=TargetN;
    TargetZ.Normalize();
    CoordType U=CoordType(1,0,0);
    if (fabs(TargetZ*U)>0.999)
        U=CoordType(0,1,0);

    CoordType TargetX=TargetZ^U;
    CoordType TargetY=TargetX^TargetZ;
    TargetX.Normalize();
    TargetY.Normalize();
    vcg::Matrix33<ScalarType> RotFrame=vcg::TransformationMatrix(TargetX,TargetY,TargetZ);
    vcg::Matrix33<ScalarType> RotFrameInv=vcg::Inverse(RotFrame);
    std::vector<vcg::Point2<ScalarType> > Cross2D;
    ///rotate each vector to transform to 2D
    for (size_t i=0;i<V.size();i++)
    {
        CoordType NF=Norm[i];
        NF.Normalize();
        CoordType Vect=V[i];
        Vect.Normalize();
        //ScalarType Dot=fabs(Vect*NF);
        //std::cout << "V[i] " << V[i].X() << " " << V[i].Y() << std::endl << std::flush;

        ///rotate the vector to become tangent to the reference plane
        vcg::Matrix33<ScalarType> RotNorm=vcg::RotationMatrix(Norm[i],TargetN);
        //std::cout << "Norm[i] " << Norm[i].X() << " " << Norm[i].Y() << " " << Norm[i].Z()<< std::endl;
        //std::cout << "TargetN " << TargetN.X() << " " << TargetN.Y() << " " << TargetN.Z()<< std::endl<< std::flush;

        CoordType rotV=RotNorm*V[i];
        //assert(fabs(rotV*TargetN)<0.000001);
        rotV.Normalize();
        //std::cout << "rotV " << rotV.X() << " " << rotV.Y() << " " << rotV.Z()<< std::endl<< std::flush;

        ///trassform to the reference frame
        rotV=RotFrame*rotV;
        assert(!isnan(rotV.X()));
        assert(!isnan(rotV.Y()));

        //it's 2D from now on
        Cross2D.push_back(vcg::Point2<ScalarType>(rotV.X(),rotV.Y()));

    }

    vcg::Point2<ScalarType> AvDir2D=InterpolateNRosy2D(Cross2D,W,N);
    assert(!isnan(AvDir2D.X()));
    assert(!isnan(AvDir2D.Y()));
    CoordType AvDir3D=CoordType(AvDir2D.X(),AvDir2D.Y(),0);
    //transform back to 3D
    AvDir3D=RotFrameInv*AvDir3D;
    return AvDir3D;
}

template <class MeshType>
class CrossField
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename vcg::face::Pos<FaceType> PosType;
    typedef typename vcg::Triangle3<ScalarType> TriangleType;

private:

    static void SubDivideDir(const CoordType &Edge0,
                             const CoordType &Edge1,
                             std::vector<CoordType> &SubDEdges,
                             int Nsub)
    {
        CoordType Dir0=Edge0;
        CoordType Dir1=Edge1;
        ScalarType a=Edge0.Norm();
        Dir0.Normalize();
        Dir1.Normalize();
        /*
         *
         *
         *                     A
         *                   /  |
         *    Dir1         /    |
         *     ^         /      |
         *     |    c  /        |
         *    /      /          |b
         *   |     /            |
         *  /    /              |
         *     /                |
         *   /                  |
         * B____________________C       ->Dir0
         *          a
         *
         * */
        ScalarType BTotal=vcg::Angle(Dir0,Dir1);
        //get angle step
        ScalarType StepAngle=BTotal/(ScalarType)Nsub;
        //get the other edge
        CoordType Edge2=Edge1-Edge0;

        //and its direction
        CoordType Dir2=Edge2;
        Dir2.Normalize();
        //find angle C
        ScalarType C=vcg::Angle(Dir2,-Dir0);
        //safety checks
        assert(BTotal<=(M_PI));
        assert(BTotal>=0);
        assert(C<=(M_PI));
        assert(C>=0);

        //push the first one
        SubDEdges.push_back(Edge0);
        for (size_t i=1;i<Nsub;i++)
        {
            //find angle interval
            ScalarType B=StepAngle*(ScalarType)i;
            ScalarType A=(M_PI-C-B);
            assert(A<=(M_PI));
            assert(A>=0);
            //find current lenght
            ScalarType b=a*sin(B)/sin(A);
            //then move along the direction of edge b
            CoordType intervB=Dir2*b;
            //finally find absolute position summing edge 0
            intervB+=Edge0;
            SubDEdges.push_back(intervB);
        }
        //push the last one
        SubDEdges.push_back(Edge1);
    }

    static void FindSubDir(vcg::Triangle3<ScalarType> T3,
                             size_t Nvert,
                             std::vector<CoordType> &SubDEdges,
                             int Nsub)
    {
        CoordType P0=T3.P0(Nvert);
        CoordType P1=T3.P1(Nvert);
        CoordType P2=T3.P2(Nvert);
        P1-=P0;
        P2-=P0;
        SubDivideDir(P1,P2,SubDEdges,Nsub);
        for (size_t j=0;j<SubDEdges.size();j++)
            SubDEdges[j]+=P0;
    }

    static void SubdivideTris(vcg::Triangle3<ScalarType> T3,
                             size_t Nvert,
                             std::vector<vcg::Triangle3<ScalarType> > &SubTris,
                             int Nsub)
    {
        std::vector<CoordType> SplittedPos;
        FindSubDir(T3,Nvert,SplittedPos,Nsub);
        CoordType P0=T3.P(Nvert);
        //then create the triangles
        for (size_t j=0;j<SplittedPos.size()-1;j++)
        {
           TriangleType T(P0,SplittedPos[j+1],SplittedPos[j]);
           SubTris.push_back(T);
        }
    }

    static ScalarType Sign(ScalarType a){return (ScalarType)((a>0)?+1:-1);}

    static void FindSubTriangles(const typename vcg::face::Pos<FaceType> &vPos,
                                 std::vector<TriangleType> &SubFaces,
                                 std::vector<FaceType*> &OriginalFace,
                                 int numSub=3)
    {
        SubFaces.clear();
        OriginalFace.clear();

        //not ct on border
        assert(!vPos.IsBorder());
        //the vertex should be the one on the pos
        assert(vPos.F()->V(vPos.E())==vPos.V());

        // collect all faces around the star of the vertex
        std::vector<PosType> StarPos;
        vcg::face::VFOrderedStarFF(vPos, StarPos);

        //all the infos for the strip of triangles

        VertexType* v0=vPos.V();
        CoordType P0=v0->P();
        //create the strip of triangles
        for (size_t i = 0; i < StarPos.size(); ++i)
        {
            PosType currPos=StarPos[i];
            FaceType *CurrF=currPos.F();
            OriginalFace.push_back(CurrF);

            VertexType *v1=CurrF->V2(currPos.E());
            VertexType *v2=CurrF->V1(currPos.E());
            CoordType P1=v1->P();
            CoordType P2=v2->P();
            assert(v1!=v2);
            assert(v0!=v1);
            assert(v0!=v2);

            SubdivideTris(vcg::Triangle3<ScalarType>(P0,P1,P2),0,SubFaces,numSub);
        }
    }

    static void InterpolateCrossSubTriangles(const std::vector<TriangleType> &SubFaces,
                                             const std::vector<FaceType*> &OriginalFace,
                                             std::vector<CoordType> &Dir,
                                             std::vector<CoordType> &NormSubF)
    {
        Dir.clear();
        //check
        assert(SubFaces.size()>OriginalFace.size());

        int SubDivisionSize=SubFaces.size()/OriginalFace.size();
        assert(SubDivisionSize>=3);
        assert(SubDivisionSize%2==1);//should be odd
        int MiddlePos=SubDivisionSize/2+1;//the one in the middle that should preserve face's cross field

        //calculate the interpolation weights and intervals
        std::vector<std::pair<FaceType*,FaceType*> > InterpFaces;
        std::vector<ScalarType> InterpWeights;

        //the one in the middle should spond to one of the
        //two faces the rest is interpolated
        for (size_t i=0;i<SubFaces.size();i++)
        {
            //find the interval and get the index in the sub_interval
            int interval=i/SubDivisionSize;
            int sub_int=i%SubDivisionSize;

            int IndexF0=-1,IndexF1=-1;
            ScalarType alpha;

            //get the index of the two faces upon which to interpolate
            if (sub_int<MiddlePos)
            {
                IndexF0=(interval+(OriginalFace.size()-1)) % OriginalFace.size();
                IndexF1=interval;
                alpha=1-(ScalarType)(sub_int+MiddlePos)/(ScalarType)SubDivisionSize;
            }
            else
            if (sub_int>MiddlePos)
            {
                IndexF0=interval;
                IndexF1=(interval+1) % OriginalFace.size();
                alpha=1-(sub_int-MiddlePos)/(ScalarType)SubDivisionSize;
            }
            else //sub_int==MiddlePos
            {
                IndexF0=interval;
                IndexF1=(interval+1) % OriginalFace.size();
                alpha=1;
            }
            assert((IndexF0>=0)&&(IndexF0<OriginalFace.size()));
            assert((IndexF1>=0)&&(IndexF1<OriginalFace.size()));

            FaceType* F0=OriginalFace[IndexF0];
            FaceType* F1=OriginalFace[IndexF1];

            InterpFaces.push_back(std::pair<FaceType*,FaceType*>(F0,F1));
            InterpWeights.push_back(alpha);
        }
        assert(InterpFaces.size()==InterpWeights.size());
        //then calculate the interpolated cross field
        for (size_t i=0;i<InterpFaces.size();i++)
        {
            std::vector<vcg::Point3<ScalarType> > TangVect;
            std::vector<vcg::Point3<ScalarType> > Norms;
            std::vector<ScalarType>  W;
            Norms.push_back(InterpFaces[i].first->N());
            Norms.push_back(InterpFaces[i].second->N());
            CoordType Dir0=InterpFaces[i].first->PD1();
            CoordType Dir1=InterpFaces[i].second->PD1();
            TangVect.push_back(Dir0);
            TangVect.push_back(Dir1);
            ScalarType CurrW=InterpWeights[i];
            W.push_back(CurrW);
            W.push_back(1-CurrW);

            CoordType TargetN=InterpFaces[i].first->N();
            if (CurrW<0.5)
                TargetN=InterpFaces[i].second->N();

            NormSubF.push_back(TargetN);
            //CoordType TargetN=vcg::Normal(SubFaces[i].P(0),SubFaces[i].P(1),SubFaces[i].P(2));
            TargetN.Normalize();
            CoordType InterpD0=InterpolateNRosy3D(TangVect,Norms,W,4,TargetN);
            //CoordType InterpD0=Dir1;
            //if (CurrW>0.5)InterpD0=Dir0;

            Dir.push_back(InterpD0);
        }
    }

    static ScalarType turn (const CoordType &norm, const CoordType &d0, const CoordType &d1)
    {
        //first check if they are coplanar up to a certain delta
        return ((d0 ^ d1).normalized() * norm);
    }

    static void InterpolateDir(const CoordType &Dir0,
                               const CoordType &Dir1,
                               const CoordType &Sep0,
                               const CoordType &Sep1,
                               const TriangleType &t0,
                               const TriangleType &t1,
                               CoordType &Interpolated,
                               size_t &Face)
    {
       //find smallest edge
       ScalarType smallestE=std::numeric_limits<ScalarType>::max();
       for (int j=0;j<3;j++)
       {
          ScalarType L0=(t0.P0(j)-t0.P1(j)).Norm();
          ScalarType L1=(t1.P0(j)-t1.P1(j)).Norm();
          if (L0<smallestE) smallestE=L0;
          if (L1<smallestE) smallestE=L1;
       }

       //safety check
       assert(t0.P(0)==t1.P(0));

       CoordType Origin=t0.P(0);
       TriangleType T0Rot(CoordType(0,0,0),t0.P(1)-Origin,t0.P(2)-Origin);
       TriangleType T1Rot(CoordType(0,0,0),t1.P(1)-Origin,t1.P(2)-Origin);

       //then rotate normal of T0 to match with normal of T1
       CoordType N0=vcg::Normal(T0Rot.cP(0),T0Rot.cP(1),T0Rot.P(2));
       CoordType N1=vcg::Normal(T1Rot.cP(0),T1Rot.cP(1),T1Rot.cP(2));
       N0.Normalize();
       N1.Normalize();
       vcg::Matrix33<ScalarType> rotation=vcg::RotationMatrix(N0,N1);

       //transform the first triangle
       T0Rot.P(1)=rotation*T0Rot.P(1);
       T0Rot.P(2)=rotation*T0Rot.P(2);

       //also rotate the directions
       CoordType Dir0Rot=rotation*Dir0;
       CoordType Dir1Rot=Dir1;
       CoordType Sep0Rot=rotation*Sep0;
       CoordType Sep1Rot=Sep1;

       //find the interpolation angles
       ScalarType Angle0=vcg::Angle(Dir0Rot,Sep0Rot);
       ScalarType Angle1=vcg::Angle(Dir1Rot,Sep1Rot);
       assert(Angle0>=0);
       assert(Angle1>=0);
       assert(Angle0<=(M_PI/2));
       assert(Angle1<=(M_PI/2));
       ScalarType alpha=0.5;//Angle0/(Angle0+Angle1);

       //then interpolate the direction
       //CoordType Interp=Dir0Rot*(1-alpha)+Dir1Rot*alpha;
       Interpolated=Sep0Rot*(1-alpha)+Sep1Rot*alpha;
       Interpolated.Normalize();
       Interpolated*=smallestE;

       //then find the triangle which falls into
       CoordType bary0,bary1;
       bool Inside0=vcg::InterpolationParameters(T0Rot,Interpolated,bary0);
       bool Inside1=vcg::InterpolationParameters(T1Rot,Interpolated,bary1);
       assert(Inside0 || Inside1);
//       if (!(Inside0 || Inside1))
//       {
//           std::cout << "Not Inside " << Interpolated.X() << "," << Interpolated.Y() << "," << Interpolated.Z() << std::endl;
//           std::cout << "bary0 " << bary0.X() << "," << bary0.Y() << "," << bary0.Z() << std::endl;
//           std::cout << "bary1 " << bary1.X() << "," << bary1.Y() << "," << bary1.Z() << std::endl;
//           std::cout << "Diff0 " << fabs(bary0.Norm() - 1) << std::endl;
//           std::cout << "Diff1 " << fabs(bary1.Norm() - 1) << std::endl;
//       }

       if (Inside0)
       {
           Interpolated=t0.P(0)*bary0.X()+t0.P(1)*bary0.Y()+t0.P(2)*bary0.Z();
           Interpolated-=Origin;
           Face=0;
       }
       else
           Face=1;

       //otherwise is already in the tangent space of t0
       Interpolated.Normalize();
    }

    static void ReduceOneDirectionField(std::vector<CoordType> &directions,
                                        std::vector<FaceType*> &faces)
    {
        //compute distribution and find closest pait
        std::vector<ScalarType> DirProd;
        ScalarType MaxProd=-1;
        int MaxInd0=-1,MaxInd1=-1;
        for (size_t i=0;i<directions.size();i++)
        {
            ScalarType prod=0;
            for (size_t j=1;j<directions.size();j++)
            {
                int Index0=i;
                int Index1=(i+j)%directions.size();
                CoordType Dir0=directions[Index0];
                CoordType Dir1=directions[Index1];
                ScalarType currP=(Dir0*Dir1);
                if (currP>MaxProd)
                {
                    MaxProd=currP;
                    MaxInd0=Index0;
                    MaxInd1=Index1;
                }
                prod+=currP;
            }
            DirProd.push_back(prod);
        }
        assert(MaxInd0!=MaxInd1);

        //then find the one to be deleted
        int IndexDel=MaxInd0;
        if (DirProd[MaxInd1]>DirProd[MaxInd0])IndexDel=MaxInd1;

        std::vector<CoordType> SwapV(directions.begin(),directions.end());
        std::vector<FaceType*> SwapF(faces.begin(),faces.end());

        directions.clear();
        faces.clear();

        for (size_t i=0;i<SwapV.size();i++)
        {
            if (i==IndexDel)continue;
            directions.push_back(SwapV[i]);
            faces.push_back(SwapF[i]);
        }
    }

public:

    static size_t FindSeparatrices(const typename vcg::face::Pos<FaceType> &vPos,
                                std::vector<CoordType> &directions,
                                std::vector<FaceType*> &faces,
                                std::vector<TriangleType> &WrongTris,
                                int expVal=-1,
                                int numSub=3)
    {

        directions.clear();

        //not ct on border
        assert(!vPos.IsBorder());
        //the vertex should be the one on the pos
        assert(vPos.F()->V(vPos.E())==vPos.V());

        std::vector<TriangleType> SubFaces;
        std::vector<CoordType> Dir1,Dir2;
        std::vector<CoordType> Norms;
        std::vector<FaceType*> OriginalFaces;

        //find subfaces
        FindSubTriangles(vPos,SubFaces,OriginalFaces,numSub);

        //interpolate the cross field
        InterpolateCrossSubTriangles(SubFaces,OriginalFaces,Dir1,Norms);

        //then calculate the orthogonal
        for (size_t i=0;i<Dir1.size();i++)
        {
            CoordType TargetN=Norms[i];//vcg::Normal(SubFaces[i].P(0),SubFaces[i].P(1),SubFaces[i].P(2));
            CoordType CrossDir2=Dir1[i]^TargetN;
            CrossDir2.Normalize();
            Dir2.push_back(CrossDir2);
        }

        //then check the triangles
        CoordType CenterStar=vPos.V()->P();
        for (size_t i = 0; i < SubFaces.size(); ++i)
        {
            TriangleType t0=SubFaces[i];
            TriangleType t1=SubFaces[(i+1)%SubFaces.size()];
            CoordType N0=Norms[i];//vcg::Normal(t0.P(0),t0.P(1),t0.P(2));
            CoordType N1=Norms[(i+1)%Norms.size()];//vcg::Normal(t1.P(0),t1.P(1),t1.P(2));
            N0.Normalize();
            N1.Normalize();

            std::vector<CoordType> SubDEdges0;
            std::vector<CoordType> SubDEdges1;
            FindSubDir(t0,0,SubDEdges0,2);
            FindSubDir(t1,0,SubDEdges1,2);
            CoordType Bary0=SubDEdges0[1];
            CoordType Bary1=SubDEdges1[1];
            //then get the directions to the barycenter
            Bary0-=CenterStar;
            Bary1-=CenterStar;
            Bary0.Normalize();
            Bary1.Normalize();

            //then get the cross field of the 2 faces
            CoordType Dir1F0=Dir1[i];
            CoordType Dir2F0=Dir2[i];

            assert(fabs(Dir1F0*N0)<0.001);
            assert(fabs(Dir1F0*Dir2F0)<0.001);

            if ((Dir1F0*Bary0)<0)Dir1F0=-Dir1F0;
            if ((Dir2F0*Bary0)<0)Dir2F0=-Dir2F0;

            CoordType Dir1F1=Dir1[(i+1)%Dir1.size()];
            CoordType Dir2F1=Dir2[(i+1)%Dir2.size()];

            assert(fabs(Dir1F1*N1)<0.001);
            assert(fabs(Dir1F1*Dir2F1)<0.01);

            //find the most similar rotation of the cross field
            vcg::Matrix33<ScalarType> rotation=vcg::RotationMatrix(N0,N1);
            CoordType Dir1F0R=rotation*Dir1F0;
            CoordType Dir2F0R=rotation*Dir2F0;

            //then get the closest upf to K*PI/2 rotations
            Dir1F1=vcg::tri::CrossField<MeshType>::K_PI(Dir1F1,Dir1F0R,N1);
            Dir2F1=vcg::tri::CrossField<MeshType>::K_PI(Dir2F1,Dir2F0R,N1);

            //then check if cross the direction of the barycenter
            ScalarType versef0D1 = turn(N0, Bary0, Dir1F0);
            ScalarType versef1D1 = turn(N1, Bary1, Dir1F1);
            ScalarType versef0D2 = turn(N0, Bary0, Dir2F0);
            ScalarType versef1D2 = turn(N1, Bary1, Dir2F1);

            if ((Bary0*Bary1)<0)
            {

                WrongTris.push_back(t0);
                WrongTris.push_back(t1);
            }


            CoordType InterpDir;
            size_t tri_Index=-1;
            if ((versef0D1 * versef1D1) < ScalarType(0))
            {
                InterpolateDir(Dir1F0,Dir1F1,Bary0,Bary1,t0,t1,InterpDir,tri_Index);
            }
            else
            {
                if ((versef0D2 * versef1D2 )< ScalarType(0))
                InterpolateDir(Dir2F0,Dir2F1,Bary0,Bary1,t0,t1,InterpDir,tri_Index);
            }
            //no separatrix found continue
            if (tri_Index==-1)continue;

            //retrieve original face
            assert((tri_Index==0)||(tri_Index==1));
            int OrigFIndex=((i+tri_Index)%SubFaces.size())/numSub;
            assert(OrigFIndex>=0);
            assert(OrigFIndex<OriginalFaces.size());

            FaceType* currF=OriginalFaces[OrigFIndex];
            //add the data
            directions.push_back(InterpDir);
            faces.push_back(currF);
        }
        if (expVal==-1)return directions.size();
        if (directions.size()<=expVal)return directions.size();

        size_t sampledDir=directions.size();
        int to_erase=directions.size()-expVal;
        do
        {
            ReduceOneDirectionField(directions,faces);
            to_erase--;
        }while (to_erase!=0);
        return sampledDir;
    }

    static CoordType FollowDirection(const FaceType &f0,
                                     const  FaceType &f1,
                                     const  CoordType &dir0)
    {
        ///first it rotate dir to match with f1
        CoordType dirR=vcg::tri::CrossField<MeshType>::Rotate(f0,f1,dir0);
        ///then get the closest upf to K*PI/2 rotations
        CoordType dir1=f1.cPD1();
        CoordType ret=vcg::tri::CrossField<MeshType>::K_PI(dir1,dirR,f1.cN());
        return ret;
    }

    static int FollowDirectionI(const FaceType &f0,
                                const  FaceType &f1,
                                const  CoordType &dir0)
    {
        ///first it rotate dir to match with f1
        CoordType dirTarget=FollowDirection(f0,f1,dir0);
        CoordType dir[4];
        CrossVector(f1,dir);
        ScalarType best=-1;
        int ret=-1;
        for (int i=0;i<4;i++)
        {
            ScalarType dot=dir[i]*dirTarget;
            if (dot>best)
            {
                best=dot;
                ret=i;
            }
        }
        assert(ret!=-1);

        return ret;
    }

    static int FollowDirection(const FaceType &f0,
                               const  FaceType &f1,
                               int dir0)
    {
        ///first it rotate dir to match with f1
        CoordType dirS=CrossVector(f0,dir0);
        CoordType dirR=vcg::tri::CrossField<MeshType>::Rotate(f0,f1,dirS);
        ///then get the closest upf to K*PI/2 rotations
        CoordType dir1=f1.cPD1();
        //int ret=I_K_PI(dir1,dirR,f1.cN());
        CoordType dir[4];
        CrossVector(f1,dir);
        ScalarType best=-1;
        int ret=-1;
        for (int i=0;i<4;i++)
        {
            ScalarType dot=dir[i]*dirR;
            if (dot>best)
            {
                best=dot;
                ret=i;
            }
        }
        assert(ret!=-1);

        return ret;
    }

    static int FollowLineDirection(const FaceType &f0,
                                   const FaceType &f1,
                                   int dir)
    {
        ///first it rotate dir to match with f1
        CoordType dir0=CrossVector(f0,dir);
        CoordType dir0R=vcg::tri::CrossField<MeshType>::Rotate(f0,f1,dir0);
        ///then get the closest upf to K*PI/2 rotations
        CoordType dir1_test=CrossVector(f1,dir);
        CoordType dir2_test=-dir1_test;
        if ((dir1_test*dir0R)>(dir2_test*dir0R))
            return dir;
        return ((dir+2)%4);

    }

    ///fird a tranformation matrix to transform
    ///the 3D space to 2D tangent space specified
    ///by the cross field (where Z=0)
    static vcg::Matrix33<ScalarType> TransformationMatrix(const FaceType &f)
    {
        typedef typename FaceType::CoordType CoordType;
        typedef typename FaceType::ScalarType ScalarType;

        ///transform to 3d
        CoordType axis0=f.cPD1();
        CoordType axis1=f.cPD2();//axis0^f.cN();
        CoordType axis2=f.cN();

        return (vcg::TransformationMatrix(axis0,axis1,axis2));
    }


    ///transform a given angle in tangent space wrt X axis of
    ///tangest space will return the sponding 3D vector
    static CoordType TangentAngleToVect(const FaceType &f,const ScalarType &angle)
    {
        ///find 2D vector
        vcg::Point2<ScalarType> axis2D=vcg::Point2<ScalarType>(cos(angle),sin(angle));
        CoordType axis3D=CoordType(axis2D.X(),axis2D.Y(),0);
        vcg::Matrix33<ScalarType> Trans=TransformationMatrix(f);
        vcg::Matrix33<ScalarType> InvTrans=Inverse(Trans);
        ///then transform
        return (InvTrans*axis3D);
    }

    ///find an angle with respect to dirX on the plane perpendiculr to DirZ
    ///dirX and dirZ should be perpendicular
    static ScalarType TangentVectToAngle(const CoordType dirX,
                                         const CoordType dirZ,
                                         const CoordType &vect3D)
    {
        const CoordType dirY=dirX^dirZ;
        dirX.Normalize();
        dirY.Normalize();
        dirZ.Normalize();
        vcg::Matrix33<ScalarType> Trans=TransformationMatrix(dirX,dirY,dirZ);
        ///trensform the vector to the reference frame by rotating it
        CoordType vect_transf=Trans*vect3D;

        ///then put to zero to the Z coordinate
        vcg::Point2<ScalarType> axis2D=vcg::Point2<ScalarType>(vect_transf.X(),vect_transf.Y());
        axis2D.Normalize();

        ///then find the angle with respact to axis 0
        ScalarType alpha=atan2(axis2D.Y(),axis2D.X());	////to sum up M_PI?
        if (alpha<0)
            alpha=(2*M_PI+alpha);
        if (alpha<0)
            alpha=0;
        return alpha;
    }

    ///find an angle with respect to the tangent frame of given face
    static ScalarType VectToAngle(const FaceType &f,const CoordType &vect3D)
    {
        vcg::Matrix33<ScalarType> Trans=TransformationMatrix(f);

        ///trensform the vector to the reference frame by rotating it
        CoordType vect_transf=Trans*vect3D;

        ///then put to zero to the Z coordinate
        vcg::Point2<ScalarType> axis2D=vcg::Point2<ScalarType>(vect_transf.X(),vect_transf.Y());
        axis2D.Normalize();

        ///then find the angle with respact to axis 0
        ScalarType alpha=atan2(axis2D.Y(),axis2D.X());	////to sum up M_PI?
        if (alpha<0)
            alpha=(2*M_PI+alpha);
        if (alpha<0)
            alpha=0;
        return alpha;
    }

    ///transform a cross field into a couple of angles
    static void CrossFieldToAngles(const FaceType &f,
                                   ScalarType &alpha1,
                                   ScalarType &alpha2,
                                   int RefEdge=1)
    {
        CoordType axis0=f.cP1(RefEdge)-f.cP0(RefEdge);
        axis0.Normalize();
        CoordType axis2=f.cN();
        axis2.Normalize();
        CoordType axis1=axis2^axis0;
        axis1.Normalize();


        vcg::Matrix33<ScalarType> Trans=vcg::TransformationMatrix(axis0,axis1,axis2);

        //trensform the vector to the reference frame by rotating it
        CoordType trasfPD1=Trans*f.cPD1();
        CoordType trasfPD2=Trans*f.cPD2();

        //then find the angle with respact to axis 0
        alpha1=atan2(trasfPD1.Y(),trasfPD1.X());
        alpha2=atan2(trasfPD2.Y(),trasfPD2.X());
    }

    ///transform a cross field into a couple of angles
    static void AnglesToCrossField(FaceType &f,
                                   const ScalarType &alpha1,
                                   const ScalarType &alpha2,
                                   int RefEdge=1)
    {
          CoordType axis0=f.cP1(RefEdge)-f.cP0(RefEdge);
          axis0.Normalize();
          CoordType axis2=f.cN();
          axis2.Normalize();
          CoordType axis1=axis2^axis0;
          axis1.Normalize();

          vcg::Matrix33<ScalarType> Trans=vcg::TransformationMatrix(axis0,axis1,axis2);
          vcg::Matrix33<ScalarType> InvTrans=Inverse(Trans);

          CoordType PD1=CoordType(cos(alpha1),sin(alpha1),0);
          CoordType PD2=CoordType(cos(alpha2),sin(alpha2),0);

          //then transform and store in the face
          f.PD1()=(InvTrans*PD1);
          f.PD2()=(InvTrans*PD2);
    }

    ///return the 4 directiona of the cross field in 3D
    ///given a first direction as input
    static void CrossVector(const CoordType &dir0,
                            const CoordType &norm,
                            CoordType axis[4])
    {
        axis[0]=dir0;
        axis[1]=norm^axis[0];
        axis[2]=-axis[0];
        axis[3]=-axis[1];
    }

    ///return the 4 direction in 3D of
    ///the cross field of a given face
    static void CrossVector(const FaceType &f,
                            CoordType axis[4])
    {
        CoordType dir0=f.cPD1();
        CoordType dir1=f.cPD2();
        axis[0]=dir0;
        axis[1]=dir1;
        axis[2]=-dir0;
        axis[3]=-dir1;
    }

    ///return the 4 direction in 3D of
    ///the cross field of a given face
    static void CrossVector(const VertexType &v,
                            CoordType axis[4])
    {
        CoordType dir0=v.cPD1();
        CoordType dir1=v.cPD2();
        axis[0]=dir0;
        axis[1]=dir1;
        axis[2]=-dir0;
        axis[3]=-dir1;
    }

    ///return a specific direction given an integer 0..3
    ///considering the reference direction of the cross field
    static CoordType CrossVector(const FaceType &f,
                                 const int &index)
    {
        assert((index>=0)&&(index<4));
        CoordType axis[4];
        CrossVector(f,axis);
        return axis[index];
    }

    ///return a specific direction given an integer 0..3
    ///considering the reference direction of the cross field
    static CoordType CrossVector(const VertexType &v,
                                 const int &index)
    {
        assert((index>=0)&&(index<4));
        CoordType axis[4];
        CrossVector(v,axis);
        return axis[index];
    }

    ///set the cross field of a given face
    static void SetCrossVector(FaceType &f,
                               CoordType dir0,
                               CoordType dir1)
    {
        f.PD1()=dir0;
        f.PD2()=dir1;
    }

    ///set the face cross vector from vertex one
    static void SetFaceCrossVectorFromVert(FaceType &f)
    {
        const CoordType &t0=f.V(0)->PD1();
        const CoordType &t1=f.V(1)->PD1();
        const CoordType &t2=f.V(2)->PD1();
        const CoordType &N0=f.V(0)->N();
        const CoordType &N1=f.V(0)->N();
        const CoordType &N2=f.V(0)->N();
        const CoordType &NF=f.N();
        const CoordType bary=CoordType(0.33333,0.33333,0.33333);
        CoordType tF0,tF1;
        tF0=InterpolateCrossField(t0,t1,t2,N0,N1,N2,NF,bary);
        tF1=NF^tF0;
        tF0.Normalize();
        tF1.Normalize();
        SetCrossVector(f,tF0,tF1);

        //then set the magnitudo
        ScalarType mag1,mag2;
        for (int i=0;i<3;i++)
        {
           vcg::Matrix33<ScalarType> rotN=vcg::RotationMatrix(f.V(i)->N(),f.N());
           CoordType rotatedDir=rotN*f.V(i)->PD1();

           if (fabs(rotatedDir*tF0)>fabs(rotatedDir*tF1))
           {
               mag1+=fabs(f.V(i)->K1());
               mag2+=fabs(f.V(i)->K2());
           }
           else
           {
               mag1+=fabs(f.V(i)->K2());
               mag2+=fabs(f.V(i)->K1());
           }
        }

        f.K1()=mag1/(ScalarType)3;
        f.K2()=mag2/(ScalarType)3;

    }

    static void SetFaceCrossVectorFromVert(MeshType &mesh)
    {
        for (unsigned int i=0;i<mesh.face.size();i++)
        {
            FaceType *f=&mesh.face[i];
            if (f->IsD())continue;
            SetFaceCrossVectorFromVert(*f);
        }
    }

    ///set the face cross vector from vertex one
    static void SetVertCrossVectorFromFace(VertexType &v)
    {
        std::vector<FaceType *> faceVec;
        std::vector<int> index;
        vcg::face::VFStarVF(&v,faceVec,index);
        std::vector<CoordType> TangVect;
        std::vector<CoordType> Norms;
        for (unsigned int i=0;i<faceVec.size();i++)
        {
            TangVect.push_back(faceVec[i]->PD1());
            Norms.push_back(faceVec[i]->N());
        }
        std::vector<ScalarType> Weights(TangVect.size(),1.0/(ScalarType)TangVect.size());
        CoordType NRef=v.N();
        CoordType N0=faceVec[0]->N();
        CoordType DirRef=faceVec[0]->PD1();
        ///find the rotation matrix that maps between normals
        vcg::Matrix33<ScalarType> rotation=vcg::RotationMatrix(N0,NRef);
        DirRef=rotation*DirRef;


        CoordType tF1=vcg::tri::CrossField<MeshType>::InterpolateCrossField(TangVect,Weights,Norms,NRef);
        tF1.Normalize();
        CoordType tF2=NRef^tF1;
        tF2.Normalize();
        v.PD1()=tF1;
        v.PD2()=tF2;
    }

    static void SetVertCrossVectorFromFace(MeshType &mesh)
    {
        for (unsigned int i=0;i<mesh.vert.size();i++)
        {
            VertexType *v=&mesh.vert[i];
            if (v->IsD())continue;
            SetVertCrossVectorFromFace(*v);
        }
    }

    ///rotate a given vector from the tangent space
    ///of f0 to the tangent space of f1 by considering the difference of normals
    static CoordType Rotate(const FaceType &f0,const FaceType &f1,const CoordType &dir3D)
    {
        CoordType N0=f0.cN();
        CoordType N1=f1.cN();

        ///find the rotation matrix that maps between normals
        vcg::Matrix33<ScalarType> rotation=vcg::RotationMatrix(N0,N1);
        CoordType rotated=rotation*dir3D;
        return rotated;
    }

    // returns the 90 deg rotation of a (around n) most similar to target b
    /// a and b should be in the same plane orthogonal to N
    static CoordType K_PI(const CoordType &a, const CoordType &b, const CoordType &n)
    {
        CoordType c = (a^n).normalized();///POSSIBLE SOURCE OF BUG CHECK CROSS PRODUCT
        ScalarType scorea = a*b;
        ScalarType scorec = c*b;
        if (fabs(scorea)>=fabs(scorec)) return a*Sign(scorea); else return c*Sign(scorec);
    }

    // returns the 90 deg rotation of a (around n) most similar to target b
    /// a and b should be in the same plane orthogonal to N
    static int I_K_PI(const CoordType &a, const CoordType &b, const CoordType &n)
    {
        CoordType c = (n^a).normalized();
        ScalarType scorea = a*b;
        ScalarType scorec = c*b;
        if (fabs(scorea)>=fabs(scorec))///0 or 2
        {
            if (scorea>0)return 0;
            return 2;
        }else ///1 or 3
        {
            if (scorec>0)return 1;
            return 3;
        }
    }

    ///interpolate cross field with barycentric coordinates
    static CoordType InterpolateCrossField(const CoordType &t0,
                                           const CoordType &t1,
                                           const CoordType &t2,
                                           const CoordType &n0,
                                           const CoordType &n1,
                                           const CoordType &n2,
                                           const CoordType &target_n,
                                           const CoordType &bary)
    {
        std::vector<CoordType > V,Norm;
        std::vector<ScalarType > W;
        V.push_back(t0);
        V.push_back(t1);
        V.push_back(t2);
        Norm.push_back(n0);
        Norm.push_back(n1);
        Norm.push_back(n2);
        W.push_back(bary.X());
        W.push_back(bary.Y());
        W.push_back(bary.Z());

        CoordType sum=vcg::tri::InterpolateNRosy3D(V,Norm,W,4,target_n);
        return sum;
    }

    ///interpolate cross field with barycentric coordinates using normalized weights
    static  CoordType InterpolateCrossField(const std::vector<CoordType> &TangVect,
                                            const std::vector<ScalarType> &Weight,
                                            const std::vector<CoordType> &Norms,
                                            const CoordType &BaseNorm)
    {

        CoordType sum=InterpolateNRosy3D(TangVect,Norms,Weight,4,BaseNorm);
        return sum;
    }

    ///interpolate cross field with scalar weight
    static typename FaceType::CoordType InterpolateCrossFieldLine(const typename FaceType::CoordType &t0,
                                                                  const typename FaceType::CoordType &t1,
                                                                  const typename FaceType::CoordType &n0,
                                                                  const typename FaceType::CoordType &n1,
                                                                  const typename FaceType::CoordType &target_n,
                                                                  const typename FaceType::ScalarType &weight)
    {
        std::vector<CoordType > V,Norm;
        std::vector<ScalarType > W;
        V.push_back(t0);
        V.push_back(t1);
        Norm.push_back(n0);
        Norm.push_back(n1);
        W.push_back(weight);
        W.push_back(1-weight);
        InterpolateNRosy3D(V,Norm,&W,4,target_n);
    }


    ///return the difference of two cross field, values between [0,1]
    static typename FaceType::ScalarType DifferenceCrossField(const typename FaceType::CoordType &t0,
                                                              const typename FaceType::CoordType &t1,
                                                              const typename FaceType::CoordType &n)
    {
        CoordType trans0=t0;
        CoordType trans1=K_PI(t1,t0,n);
        ScalarType diff = vcg::AngleN(trans0,trans1)/(M_PI_4);
        return diff;
    }

    ///return the difference of two cross field, values between [0,1]
    static typename FaceType::ScalarType DifferenceLineField(const typename FaceType::CoordType &t0,
                                                             const typename FaceType::CoordType &t1,
                                                             const typename FaceType::CoordType &n)
    {
        CoordType trans0=t0;
        CoordType trans1=t1;
        if ((trans0*trans1)<0)trans1=-trans1;
        ScalarType angleD=vcg::Angle(trans0,trans1);
        assert(angleD>=0);
        assert(angleD<=M_PI_2);
        return (angleD/M_PI_2);
    }

    ///return the difference of two cross field, values between [0,1]
    static typename FaceType::ScalarType DifferenceCrossField(const typename vcg::Point2<ScalarType> &t0,
                                                              const typename vcg::Point2<ScalarType> &t1)
    {
        CoordType t03D=CoordType(t0.X(),t0.Y(),0);
        CoordType t13D=CoordType(t1.X(),t1.Y(),0);
        CoordType Norm=CoordType(0,0,1);
//        CoordType n=CoordType(0,0,1);
//        CoordType trans1=K_PI(t13D,t03D,n);
//        ScalarType diff=vcg::AngleN(trans0,trans1)/(M_PI_4);
        //ScalarType diff = 1-fabs(trans0*trans1);
        return DifferenceCrossField(t03D,t13D,Norm);
    }

    ///return the difference of two cross field, values between [0,1]
    static typename FaceType::ScalarType DifferenceLineField(const typename vcg::Point2<ScalarType> &t0,
                                                              const typename vcg::Point2<ScalarType> &t1)
    {
        CoordType t03D=CoordType(t0.X(),t0.Y(),0);
        CoordType t13D=CoordType(t1.X(),t1.Y(),0);
        CoordType Norm=CoordType(0,0,1);
//        CoordType n=CoordType(0,0,1);
//        CoordType trans1=K_PI(t13D,t03D,n);
//        ScalarType diff=vcg::AngleN(trans0,trans1)/(M_PI_4);
        //ScalarType diff = 1-fabs(trans0*trans1);
        return DifferenceLineField(t03D,t13D,Norm);
    }

    ///compute the mismatch between 2 directions
    ///each one si perpendicular to its own normal
    static int MissMatchByCross(const CoordType &dir0,
                                const CoordType &dir1,
                                const CoordType &N0,
                                const CoordType &N1)
    {
        CoordType dir0Rot=Rotate(dir0,N0,N1);
        CoordType dir1Rot=dir1;

        dir0Rot.Normalize();
        dir1Rot.Normalize();

        ScalarType angle_diff=VectToAngle(dir0Rot,N0,dir1Rot);

        ScalarType step=M_PI/2.0;
        int i=(int)floor((angle_diff/step)+0.5);
        int k=0;
        if (i>=0)
            k=i%4;
        else
            k=(-(3*i))%4;
        return k;
    }

    ///compute the mismatch between 2 faces
    static int MissMatchByCross(const FaceType &f0,
                                const FaceType &f1)
    {
        //CoordType dir0=CrossVector(f0,0);
        CoordType dir1=CrossVector(f1,0);

        CoordType dir1Rot=Rotate(f1,f0,dir1);
        dir1Rot.Normalize();

        ScalarType angle_diff=VectToAngle(f0,dir1Rot);

        ScalarType step=M_PI/2.0;
        int i=(int)floor((angle_diff/step)+0.5);
        int k=0;
        if (i>=0)
            k=i%4;
        else
            k=(-(3*i))%4;
        return k;
    }


    ///return true if a given vertex is singular,
    ///return also the missmatch
    static bool IsSingularByCross(const VertexType &v,int &missmatch)
    {
        typedef typename VertexType::FaceType FaceType;
        ///check that is on border..
        if (v.IsB())return false;

        std::vector<face::Pos<FaceType> > posVec;
        //SortedFaces(v,faces);
        face::Pos<FaceType> pos(v.cVFp(), v.cVFi());
        vcg::face::VFOrderedStarFF(pos, posVec);

        int curr_dir=0;
        for (unsigned int i=0;i<posVec.size();i++)
        {
            FaceType *curr_f=posVec[i].F();
            FaceType *next_f=posVec[(i+1)%posVec.size()].F();

            //find the current missmatch
            curr_dir=FollowDirection(*curr_f,*next_f,curr_dir);
        }
        missmatch=curr_dir;
        return(curr_dir!=0);
    }

    ///select singular vertices
    static void UpdateSingularByCross(MeshType &mesh)
    {
        bool hasSingular = vcg::tri::HasPerVertexAttribute(mesh,std::string("Singular"));
        bool hasSingularIndex = vcg::tri::HasPerVertexAttribute(mesh,std::string("SingularIndex"));

        typename MeshType::template PerVertexAttributeHandle<bool> Handle_Singular;
        typename MeshType::template PerVertexAttributeHandle<int> Handle_SingularIndex;

        if (hasSingular)
            Handle_Singular=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<bool>(mesh,std::string("Singular"));
        else
            Handle_Singular=vcg::tri::Allocator<MeshType>::template AddPerVertexAttribute<bool>(mesh,std::string("Singular"));

        if (hasSingularIndex)
            Handle_SingularIndex=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("SingularIndex"));
        else
            Handle_SingularIndex=vcg::tri::Allocator<MeshType>::template AddPerVertexAttribute<int>(mesh,std::string("SingularIndex"));

        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (mesh.vert[i].IsD())continue;
            //if (mesh.vert[i].IsB())continue;

            int missmatch;
            if (IsSingularByCross(mesh.vert[i],missmatch))
            {
                Handle_Singular[i]=true;
                Handle_SingularIndex[i]=missmatch;
            }
            else
            {
                Handle_Singular[i]=false;
                Handle_SingularIndex[i]=0;
            }
        }
    }


    static void GradientToCross(const FaceType &f,
                                const vcg::Point2<ScalarType> &UV0,
                                const vcg::Point2<ScalarType> &UV1,
                                const vcg::Point2<ScalarType> &UV2,
                                CoordType &dirU,
                                CoordType &dirV)
    {
        vcg::Point2<ScalarType> Origin2D=(UV0+UV1+UV2)/3;
        CoordType Origin3D=(f.cP(0)+f.cP(1)+f.cP(2))/3;

        vcg::Point2<ScalarType> UvT0=UV0-Origin2D;
        vcg::Point2<ScalarType> UvT1=UV1-Origin2D;
        vcg::Point2<ScalarType> UvT2=UV2-Origin2D;

        CoordType PosT0=f.cP(0)-Origin3D;
        CoordType PosT1=f.cP(1)-Origin3D;
        CoordType PosT2=f.cP(2)-Origin3D;

        CoordType Bary0,Bary1;
        vcg::InterpolationParameters2(UvT0,UvT1,UvT2,vcg::Point2<ScalarType>(1,0),Bary0);
        vcg::InterpolationParameters2(UvT0,UvT1,UvT2,vcg::Point2<ScalarType>(0,1),Bary1);

        //then transport to 3D
        dirU=PosT0*Bary0.X()+PosT1*Bary0.Y()+PosT2*Bary0.Z();
        dirV=PosT0*Bary1.X()+PosT1*Bary1.Y()+PosT2*Bary1.Z();

//        dirU-=Origin3D;
//        dirV-=Origin3D;
        dirU.Normalize();
        dirV.Normalize();
        //orient coherently
        CoordType Ntest=dirU^dirV;
        CoordType NTarget=vcg::Normal(f.cP(0),f.cP(1),f.cP(2));
        if ((Ntest*NTarget)<0)dirV=-dirV;

//        //then make them orthogonal
//        CoordType dirAvg=dirU^dirV;
        CoordType dirVTarget=NTarget^dirU;
        CoordType dirUTarget=NTarget^dirV;

         dirUTarget.Normalize();
         dirVTarget.Normalize();
         if ((dirUTarget*dirU)<0)dirUTarget=-dirUTarget;
         if ((dirVTarget*dirV)<0)dirVTarget=-dirVTarget;

         dirU=(dirU+dirUTarget)/2;
         dirV=(dirV+dirVTarget)/2;

         dirU.Normalize();
         dirV.Normalize();

//        ///compute non normalized normal
//        CoordType n  =  f.cN();

//        CoordType p0 =f.cP(1) - f.cP(0);
//        CoordType p1 =f.cP(2) - f.cP(1);
//        CoordType p2 =f.cP(0) - f.cP(2);

//        CoordType t[3];
//        t[0] =  -(p0 ^ n);
//        t[1] =  -(p1 ^ n);
//        t[2] =  -(p2 ^ n);

//        dirU = t[1]*UV0.X() + t[2]*UV1.X() + t[0]*UV2.X();
//        dirV = t[1]*UV0.Y() + t[2]*UV1.Y() + t[0]*UV2.Y();
    }

    static void MakeDirectionFaceCoherent(FaceType *f0,
                                          FaceType *f1)
    {
        CoordType dir0=f0->PD1();
        CoordType dir1=f1->PD1();

        CoordType dir0Rot=Rotate(*f0,*f1,dir0);
        dir0Rot.Normalize();

        CoordType targD=K_PI(dir1,dir0Rot,f1->N());

        f1->PD1()=targD;
        f1->PD2()=f1->N()^targD;
        f1->PD2().Normalize();
    }

    static void AdjustDirectionsOnTangentspace(MeshType &mesh)
    {
        for (size_t i=0;i<mesh.face.size();i++)
        {
            FaceType *f=&mesh.face[i];
            if (f->IsD())continue;
            CoordType Ntest=mesh.face[i].PD1()^mesh.face[i].PD2();
            Ntest.Normalize();
            CoordType Ntarget=mesh.face[i].N();
            if ((Ntest*Ntarget)>0.999)continue;

            //find the rotation matrix that maps between normals
            vcg::Matrix33<ScalarType> rotation=vcg::RotationMatrix(Ntest,Ntarget);
            mesh.face[i].PD1()=rotation*mesh.face[i].PD1();
            mesh.face[i].PD2()=rotation*mesh.face[i].PD2();
        }
    }

    static void OrientDirectionFaceCoherently(MeshType &mesh)
    {
        for (size_t i=0;i<mesh.face.size();i++)
        {
            FaceType *f=&mesh.face[i];
            if (f->IsD())continue;
            CoordType Ntest= CoordType::Construct( mesh.face[i].PD1()^mesh.face[i].PD2() );
            if ((Ntest*vcg::Normal(f->P(0),f->P(1),f->P(2)))<0)mesh.face[i].PD2()=-mesh.face[i].PD2();
        }
    }

    static void MakeDirectionFaceCoherent(MeshType &mesh,
                                          bool normal_diff=true)
    {
        vcg::tri::UpdateFlags<MeshType>::FaceClearV(mesh);
        vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);

        typedef typename std::pair<FaceType*,FaceType*> FacePair;
        std::vector<std::pair<ScalarType,FacePair> > heap;

        while (true)
        {
            bool found=false;
            for (int i=0; i<(int)mesh.face.size(); i++)
            {
                FaceType *f=&(mesh.face[i]);
                if (f->IsD())continue;
                if (f->IsV())continue;
                f->SetV();
                found=true;

                for (int i=0;i<f->VN();i++)
                {
                    FaceType *Fopp=f->FFp(i);
                    if (Fopp==f)continue;

                    FacePair entry(f,Fopp);

                    ScalarType val=0;
                    if (normal_diff)val=-(f->N()-Fopp->N()).Norm();

                    heap.push_back(std::pair<ScalarType,FacePair>(val,entry));
                }
                break;
            }

            if (!found)
            {

                vcg::tri::UpdateFlags<MeshType>::FaceClearV(mesh);
                return;///all faces has been visited
            }

            std::make_heap (heap.begin(),heap.end());

            while (!heap.empty())
            {
                std::pop_heap(heap.begin(), heap.end());
                FaceType *f0=heap.back().second.first;
                FaceType *f1=heap.back().second.second;
                assert(f0->IsV());
                heap.pop_back();

                MakeDirectionFaceCoherent(f0,f1);
                f1->SetV();
                for (int k=0; k<f1->VN(); k++)
                {
                    FaceType* f2 = f1->FFp(k);
                    if (f2->IsV())continue;
                    if (f2->IsD())continue;
                    if (f2==f1)continue;

                    FacePair entry(f1,f2);

                    ScalarType val=0;
                    if (normal_diff)val=-(f1->N()-f2->N()).Norm();

                    heap.push_back(std::pair<ScalarType,FacePair>(val,entry));
                    std::push_heap (heap.begin(),heap.end());
                }
            }
        }

    }

    ///transform curvature to UV space
    static vcg::Point2<ScalarType> CrossToUV(FaceType &f,int numD=0)
    {
        typedef typename FaceType::ScalarType ScalarType;
        typedef typename FaceType::CoordType CoordType;

        CoordType Curv=CrossVector(f,numD);
        Curv.Normalize();

        CoordType bary3d=(f.P(0)+f.P(1)+f.P(2))/3.0;
        vcg::Point2<ScalarType> Uv0=f.V(0)->T().P();
        vcg::Point2<ScalarType> Uv1=f.V(1)->T().P();
        vcg::Point2<ScalarType> Uv2=f.V(2)->T().P();
        vcg::Point2<ScalarType> baryUV=(Uv0+Uv1+Uv2)/3.0;
        CoordType direct3d=bary3d+Curv;
        CoordType baryCoordsUV;
        vcg::InterpolationParameters<FaceType,ScalarType>(f,direct3d,baryCoordsUV);
        vcg::Point2<ScalarType> curvUV=baryCoordsUV.X()*Uv0+
                baryCoordsUV.Y()*Uv1+
                baryCoordsUV.Z()*Uv2-baryUV;
        curvUV.Normalize();
        return curvUV;
    }

    static void InitDirFromWEdgeUV(MeshType &mesh)
    {
        for (size_t i=0;i<mesh.face.size();i++)
        {
            vcg::Point2<ScalarType> UV0=mesh.face[i].WT(0).P();
            vcg::Point2<ScalarType> UV1=mesh.face[i].WT(1).P();
            vcg::Point2<ScalarType> UV2=mesh.face[i].WT(2).P();
            GradientToCross(mesh.face[i],UV0,UV1,UV2,
                            CoordType::Construct(mesh.face[i].PD1()),
                            CoordType::Construct(mesh.face[i].PD2()) );
        }
        OrientDirectionFaceCoherently(mesh);
    }

};///end class
} //End Namespace Tri
} // End Namespace vcg
#endif
