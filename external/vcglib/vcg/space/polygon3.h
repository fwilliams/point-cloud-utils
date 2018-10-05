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

#ifndef POLYGON_H
#define POLYGON_H

#include <vcg/space/plane3.h>
#include <vcg/space/fitting3.h>
#include <vcg/space/point_matching.h>
#include <vcg/math/matrix33.h>

namespace vcg {

////return true if the
//template <class CoordType>
//bool CheckNormalizedCoords(CoordType dir)
//{
//    typedef typename CoordType::ScalarType ScalarType;
//    if(isnan(dir.X()))return false;
//    if(isnan(dir.Y()))return false;
//    if(isnan(dir.Z()))return false;
//    ScalarType Norm=dir.Norm();
//    if(fabs(Norm-1.f)>0.01f)return false;
//    return true;
//}

//return per vertex Normals of a polygonal face stored as a vector of coords
template <class CoordType>
void GetNormals(std::vector<CoordType> &Pos,
                std::vector<CoordType> &Norms)
{
    Norms.clear();
    int size=Pos.size();
    if (size<=2) return;
    for (int i=0;i<size;i++)
        Norms.push_back(Normal(Pos[i],Pos[(i+1)%size],Pos[(i+2)%size]).Normalize());
}

//return the normal of a polygonal face stored as a vector of coords
template <class CoordType>
CoordType Normal(std::vector<CoordType> &Pos)
{
    std::vector<CoordType> Norms;
    GetNormals(Pos,Norms);
    if (Norms.size()==0)
        return(CoordType(1,0,0));

    CoordType NSum=CoordType(0,0,0);
    for (size_t i=0;i<Norms.size();i++)
        NSum+=Norms[i];

    NSum.Normalize();
    return (NSum);
}

//return the area of a polygonal face stored as a vector of coords
template <class CoordType>
typename CoordType::ScalarType Area(const std::vector<CoordType> &Pos)
{
    typedef typename CoordType::ScalarType ScalarType;
    CoordType bary=CoordType(0,0,0);
    for (int i=0;i<Pos.size();i++)
        bary+=Pos[i];

    bary/=Pos.size();
    ScalarType Area=0;
    for (size_t i=0;i<Pos.size();i++)
    {
        CoordType p0=Pos[i];
        CoordType p1=Pos[(i+1)% Pos.size()];
        CoordType p2=bary;
        vcg::Triangle3<ScalarType> T(p0,p1,p2);
        Area+=(vcg::DoubleArea(T)/2);
    }
    return Area;
}

//return per vertex Normals of a polygonal face
template<class PolygonType>
void PolyNormals(const PolygonType &F,
                 std::vector<typename PolygonType::CoordType> &Norms)
{
    Norms.clear();
    if (F.VN()<=2) return;
    for (int i=0;i<F.VN();i++)
        Norms.push_back(Normal(F.cP0(i),F.cP1(i),F.cP2(i)).Normalize());
}

//return the barycenter of a polygonal face
template<class PolygonType>
typename PolygonType::CoordType PolyBarycenter(const PolygonType &F)
{
    typename PolygonType::CoordType bary(0,0,0);
    for (int i=0;i<F.VN();i++)
        bary+=F.cP(i);

    bary/=(typename PolygonType::ScalarType)F.VN();
    return bary;
}

//return the area of a polygonal face
template<class PolygonType>
typename PolygonType::ScalarType PolyArea(const PolygonType &F)
{
    typedef typename PolygonType::CoordType CoordType;
    typedef typename PolygonType::ScalarType ScalarType;

    CoordType bary=PolyBarycenter(F);
    ScalarType Area=0;
    for (size_t i=0;i<F.VN();i++)
    {
        CoordType p0=F.cP0(i);
        CoordType p1=F.cP1(i);
        CoordType p2=bary;
        vcg::Triangle3<ScalarType> T(p0,p1,p2);
        Area+=(vcg::DoubleArea(T)/2);
    }
    return Area;
}

//return the normal of a polygonal face
template<class PolygonType>
typename PolygonType::CoordType PolygonNormal(const PolygonType &F)
{
    typename PolygonType::CoordType n(0,0,0);

    for (int i=0;i<F.VN();i++)
        n+=Normal(F.cP0(i),F.cP1(i),F.cP2(i)).Normalize();

    return n.Normalize();
}

//return the perimeter of a polygonal face
template<class PolygonType>
typename PolygonType::ScalarType PolyPerimeter(const PolygonType &F)
{
    typedef typename PolygonType::ScalarType ScalarType;

    ScalarType SumL=0;
    for (int i=0;i<F.VN();i++)
    {
        ScalarType L=(F.cP0(i)-F.cP1(i)).Norm();
        SumL+=L;
    }
    return (SumL);
}

//return a Scalar value that encode the variance of the normals
//wrt the average one (1 means hight variance, 0 no variance)
template<class PolygonType>
typename PolygonType::ScalarType PolyNormDeviation(const PolygonType &F)
{
    typedef typename PolygonType::CoordType CoordType;
    typedef typename PolygonType::ScalarType ScalarType;

    std::vector<CoordType> Norms;
    PolyNormals(F,Norms);

    //calculate the Avg Normal
    CoordType AvgNorm(0,0,0);
    for (int i=0;i<Norms.size();i++)
        AvgNorm+=Norms[i];

    AvgNorm.Normalize();

    //if (!CheckNormalizedCoords(AvgNorm))return 1;

    ScalarType Dev=0;
    for (int i=0;i<Norms.size();i++)
        Dev+=pow((Norms[i]-AvgNorm).Norm()/2.0,2);

    Dev/=(ScalarType)Norms.size();
    Dev=sqrt(Dev);
    return Dev;
}

//return a Scalar value that encode the distance wrt ideal angle for each
//wrt the average one (1 correspond to hight variance, 0 no variance)
template<class PolygonType>
void PolyAngleDeviation(const PolygonType &F,
                        typename PolygonType::ScalarType &AvgDev,
                        typename PolygonType::ScalarType &MaxDev)
{
    typedef typename PolygonType::CoordType CoordType;
    typedef typename PolygonType::ScalarType ScalarType;
    assert(F.VN()>2);
    ScalarType IdealAngle=M_PI-(2*M_PI/(ScalarType)F.VN());
    assert(IdealAngle>0);

    //then compute the angle deviation
    MaxDev=0;
    AvgDev=0;

    for (int i=0;i<F.VN();i++)
    {
        CoordType dir0=F.cP0(i)-F.cP1(i);
        CoordType dir1=F.cP2(i)-F.cP1(i);

        ScalarType VAngle=vcg::Angle(dir0,dir1);
        assert(VAngle>=0);
        ScalarType VAngleDiff=fabs(VAngle-IdealAngle);

        if (VAngleDiff>MaxDev)MaxDev=VAngleDiff;

        AvgDev+=VAngleDiff;
    }
    AvgDev/=(ScalarType)F.VN();

    AvgDev/=(M_PI/2.0);
    MaxDev/=(M_PI/2.0);

    if (AvgDev>1)AvgDev=1;
    if (MaxDev>1)MaxDev=1;
}

//return the fitting plane of a polygonal face
template<class PolygonType>
vcg::Plane3<typename PolygonType::ScalarType> PolyFittingPlane(const PolygonType &F)
{
    typedef typename PolygonType::CoordType CoordType;
    typedef typename PolygonType::ScalarType ScalarType;
    vcg::Plane3<ScalarType> BestPL;
    assert(F.VN()>=3);
    std::vector<CoordType> pointVec;
    for (int i=0;i<F.VN();i++)
        pointVec.push_back(F.cP(i));

    vcg::FitPlaneToPointSet(pointVec,BestPL);
    return BestPL;
}

//return the flatness of a polygonal face as avg distance to the best fitting plane divided by half perimeter
template<class PolygonType>
typename PolygonType::ScalarType PolyFlatness(const PolygonType &F)
{
    typedef typename PolygonType::CoordType CoordType;
    typedef typename PolygonType::ScalarType ScalarType;

    if (F.VN()<=3)
        return 0;

    //average lenght
    ScalarType SumL=PolyPerimeter(F)/2.0;

    //diagonal distance
    vcg::Plane3<ScalarType> BestPL=PolyFittingPlane(F);

    //then project points on the plane
    ScalarType Flatness=0;
    for (int i=0;i<F.VN();i++)
    {
        CoordType pos=F.cP(i);
        CoordType proj=BestPL.Projection(pos);
        Flatness+=(pos-proj).Norm();
    }
    Flatness/=(ScalarType)F.VN();
    return((Flatness)/SumL);
}

//evaluate the PCA directions of a polygonal face
template<class PolygonType>
void PolyPCA(const PolygonType &F,
             typename PolygonType::CoordType PCA[])
{
    typedef typename PolygonType::CoordType CoordType;
    typedef typename PolygonType::ScalarType ScalarType;

    //compute the covariance matrix
    Eigen::Matrix3d EigenCovMat;
    //ComputeCovarianceMatrix(EigenCovMat);
    //compute covariance matrix
    ///compute the barycenter
    CoordType Barycenter=PolyBarycenter(F);

    // second cycle: compute the covariance matrix
    EigenCovMat.setZero();
    Eigen::Vector3d p;
    for (int i=0;i<F.VN();i++)
    {
        (F.cP(i)-Barycenter).ToEigenVector(p);
        EigenCovMat+= p*p.transpose(); // outer product
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d > eig(EigenCovMat);

    Eigen::Vector3d eval = eig.eigenvalues();
    Eigen::Matrix3d evec = eig.eigenvectors();

    eval = eval.cwiseAbs();
    int normInd,maxInd,minInd;

    ///get min and max coff ..
    ///the minumum is the Normal
    ///the other two the anisotropy directions
    eval.minCoeff(&normInd);
    eval.maxCoeff(&maxInd);
    minInd=(maxInd+1)%3;

    if (minInd==normInd)minInd=(normInd+1)%3;
    assert((minInd!=normInd)&&(minInd!=maxInd)&&(minInd!=maxInd));

    ///maximum direction of PCA
    PCA[0][0] = evec(0,maxInd);
    PCA[0][1] = evec(1,maxInd);
    PCA[0][2] = evec(2,maxInd);
    ///minimum direction of PCA
    PCA[1][0] = evec(0,minInd);
    PCA[1][1] = evec(1,minInd);
    PCA[1][2] = evec(2,minInd);
    ///Normal direction
    PCA[2][0] = evec(0,normInd);
    PCA[2][1] = evec(1,normInd);
    PCA[2][2] = evec(2,normInd);

    ScalarType LX=sqrt(eval[maxInd]);
    ScalarType LY=sqrt(eval[minInd]);
    //ScalarType LZ=sqrt(eval[normInd]);

    ///scale the directions
    PCA[0]*=LX;
    PCA[1]*=LY;
    //PCA[2]*=LZ;//.Normalize();
    PCA[2].Normalize();
}

//evaluate the PCA directions of a polygonal face
//scaled by the area of the face
template<class PolygonType>
void PolyScaledPCA(const PolygonType &F,
                   typename PolygonType::CoordType PCA[])
{
    typedef typename PolygonType::CoordType CoordType;
    typedef typename PolygonType::ScalarType ScalarType;

    std::vector<CoordType> SwapPos;

    ///compute the barycenter
    //CoordType Barycenter=PolyBarycenter(F);

    ///compute the Area
    ScalarType Area=PolyArea(F);

    PolyPCA(F,PCA);

    ScalarType Scale=sqrt(Area/(PCA[0].Norm()*PCA[1].Norm()));
    PCA[0]*=Scale;
    PCA[1]*=Scale;

}

//return the base template polygon as
//described by "Static Aware Grid Shells" by Pietroni et Al.
template<class CoordType>
void getBaseTemplatePolygon(int N,
                            std::vector<CoordType> &TemplatePos)
{
    typedef typename CoordType::ScalarType ScalarType;
    ///first find positions in the
    ///reference frame of the passed matrix
    ScalarType AngleInterval=2.0*M_PI/(ScalarType)N;
    ScalarType CurrAngle=0;
    TemplatePos.resize(N);
    for (size_t i=0;i<TemplatePos.size();i++)
    {
        ///find with trigonometric functions
        TemplatePos[i].X()=cos(CurrAngle);
        TemplatePos[i].Y()=sin(CurrAngle);
        TemplatePos[i].Z()=0;
        //            TemplatePos[i].Normalize();
        //            TemplatePos[i].X()*=Anisotropy;
        //            TemplatePos[i].Y()*=(1-Anisotropy);

        ///increment the angle
        CurrAngle+=AngleInterval;
    }
}

//return the rigidly aligned template polygon as
//described by "Static Aware Grid Shells" by Pietroni et Al.
template<class PolygonType>
void GetPolyTemplatePos(const PolygonType &F,
                        std::vector<typename PolygonType::CoordType> &TemplatePos,
                        bool force_isotropy=false)
{
    typedef typename PolygonType::CoordType CoordType;
    typedef typename PolygonType::ScalarType ScalarType;
    std::vector<CoordType>  UniformPos,UniformTempl;

    CoordType Barycenter=PolyBarycenter(F);

    getBaseTemplatePolygon(F.VN(),TemplatePos);

    CoordType PCA[3];
    PolyPCA(F,PCA);

    vcg::Matrix44<ScalarType> ToPCA,ToPCAInv;
    ToPCA.SetIdentity();

    CoordType dirX=PCA[0];
    CoordType dirY=PCA[1];
    CoordType dirZ=PCA[2];

    if (force_isotropy)
    {
        dirX.Normalize();
        dirY.Normalize();
        dirZ.Normalize();
//        CoordType dirXN=dirX;dirXN.Normalize();
//        CoordType dirYN=dirY;dirYN.Normalize();
//        CoordType dirZN=dirZ;dirZN.Normalize();

//        dirX=dirX*0.8+dirXN*0.2;
//        dirY=dirY*0.8+dirYN*0.2;
//        dirZ=dirZ*0.8+dirZN*0.2;
    }

    ///set the Rotation matrix
    ToPCA.SetColumn(0,dirX);
    ToPCA.SetColumn(1,dirY);
    ToPCA.SetColumn(2,dirZ);
    ToPCAInv=ToPCA;
    ToPCA=vcg::Inverse(ToPCA);

    ///then transform the polygon to PCA space
    for (int i=0;i<F.VN();i++)
    {
        ///translate
        CoordType Pos=F.cP(i)-Barycenter;
        ///rotate
        Pos=ToPCA*Pos;
        //retranslate
        UniformPos.push_back(Pos);
    }

    ///calculate the Area
    ScalarType AreaTemplate=Area(TemplatePos);
    ScalarType AreaUniform=Area(UniformPos);

//    if (TargetArea>0)
//    {
//        AreaUniform*=(AreaUniform/TargetArea);
//    }

    ScalarType Scale=sqrt(AreaTemplate/AreaUniform);

    for (size_t i=0;i<UniformPos.size();i++)
        UniformPos[i]*=Scale;

    ///check side
    CoordType N0=Normal(UniformPos);
    CoordType N1=Normal(TemplatePos);
    if ((N0*N1)<0)std::reverse(TemplatePos.begin(),TemplatePos.end());

    ///initialize
    std::vector<CoordType> FixPoints(UniformPos.begin(),UniformPos.end());
    std::vector<CoordType> MovPoints(TemplatePos.begin(),TemplatePos.end());

    ///add displacement along Z
    for (size_t i=0;i<FixPoints.size();i++)
    {
        FixPoints[i]+=CoordType(0,0,0.1);
        MovPoints[i]+=CoordType(0,0,0.1);
    }
    ///add original points
    FixPoints.insert(FixPoints.end(),UniformPos.begin(),UniformPos.end());
    MovPoints.insert(MovPoints.end(),TemplatePos.begin(),TemplatePos.end());

    ///then find the alignment
    vcg::Matrix44<ScalarType> Rigid;
    ///compute rigid match
    vcg::ComputeRigidMatchMatrix<ScalarType>(FixPoints,MovPoints,Rigid);

    ///then apply transformation
    UniformTempl.resize(TemplatePos.size(),CoordType(0,0,0));
    for (size_t i=0;i<TemplatePos.size();i++)
        UniformTempl[i]=Rigid*TemplatePos[i];

    ///then map back to 3D space
    for (size_t i=0;i<TemplatePos.size();i++)
    {
        TemplatePos[i]=UniformTempl[i];
        TemplatePos[i]*=1/Scale;
        TemplatePos[i]=ToPCAInv*TemplatePos[i];
    }

    for (size_t i=0;i<TemplatePos.size();i++)
        TemplatePos[i]+=Barycenter;

}

//compute the aspect ratio using the rigidly aligned template polygon as
//described by "Static Aware Grid Shells" by Pietroni et Al.
template<class PolygonType>
typename PolygonType::ScalarType PolyAspectRatio(const PolygonType &F,
                                                 bool isotropic=false)
{
    typedef typename PolygonType::CoordType CoordType;
    typedef typename PolygonType::ScalarType ScalarType;
    std::vector<CoordType> TemplatePos;

    GetPolyTemplatePos(F,TemplatePos,isotropic);

    ScalarType diff=0;
    assert((int)TemplatePos.size()==F.VN());

    ScalarType AreaP=PolyArea(F);
    for (size_t i=0;i<TemplatePos.size();i++)
        diff+=pow((TemplatePos[i]-F.cP(i)).Norm(),2)/AreaP;

    return(diff);
}

}
#endif // POLYGON_H
