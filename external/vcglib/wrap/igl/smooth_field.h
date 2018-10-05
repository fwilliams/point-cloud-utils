/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2014                                                \/)\/    *
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

#ifndef SMOOTHER_FIELD_H
#define SMOOTHER_FIELD_H

//eigen stuff
#include <eigenlib/Eigen/Sparse>

//vcg stuff
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>

//igl related stuff
#include <igl/n_polyvector.h>
#include <igl/principal_curvature.h>
#include <igl/igl_inline.h>
#include <igl/comiso/nrosy.h>

namespace vcg {
namespace tri {

enum SmoothMethod{SMMiq,SMNPoly};

template <class MeshType>
class FieldSmoother
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType CoordType;


    static void InitQualityByAnisotropyDir(MeshType &mesh)
    {
        std::vector<ScalarType> QVal;
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            ScalarType N1=fabs(mesh.vert[i].K1());
            ScalarType N2=fabs(mesh.vert[i].K2());
            QVal.push_back(N1);
            QVal.push_back(N2);
        }

        std::sort(QVal.begin(),QVal.end());
        int percUp=int(floor(QVal.size()*0.95+0.5));
        ScalarType trimUp=QVal[percUp];

        for (size_t i=0;i<mesh.vert.size();i++)
        {
            ScalarType N1=(mesh.vert[i].K1());
            ScalarType N2=(mesh.vert[i].K2());

           ScalarType NMax=std::max(N1,N2)/trimUp;
           ScalarType NMin=std::min(N1,N2)/trimUp;

           if (NMax>1)NMax=1;
           if (NMax<-1)NMax=-1;

           if (NMin>1)NMin=1;
           if (NMin<-1)NMin=-1;

           ScalarType CurvAni=(NMax-NMin)/2;
           mesh.vert[i].Q()=CurvAni;
        }
        vcg::tri::UpdateQuality<MeshType>::FaceFromVertex(mesh);
    }

    static void SetEdgeDirection(FaceType *f,int edge)
    {
        CoordType dir=f->P0(edge)-f->P1(edge);
        dir.Normalize();
        ScalarType prod1=fabs(dir*f->PD1());
        ScalarType prod2=fabs(dir*f->PD2());
        if (prod1>prod2)
        {
            f->PD1()=dir;
            f->PD2()=f->N()^dir;
        }else
        {
            f->PD2()=dir;
            f->PD1()=f->N()^dir;
        }
    }

    static void AddSharpEdgesConstraints(MeshType & mesh,
                                         const ScalarType &thr=0.2)
    {
        for (size_t i=0;i<mesh.face.size();i++)
            for (int j=0;j<mesh.face[i].VN();j++)
            {
                FaceType *f0=&mesh.face[i];
                FaceType *f1=f0->FFp(j);
                if (f0==f1)continue;
                CoordType N0=f0->N();
                CoordType N1=f1->N();
                if ((N0*N1)>thr)continue;
                SetEdgeDirection(f0,j);
                f0->SetS();
            }
    }

    static void AddBorderConstraints(MeshType & mesh)
    {
        for (size_t i=0;i<mesh.face.size();i++)
            for (int j=0;j<mesh.face[i].VN();j++)
            {
                FaceType *f0=&mesh.face[i];
                FaceType *f1=f0->FFp(j);
                assert(f1!=NULL);
                if (f0!=f1)continue;
                SetEdgeDirection(f0,j);
                f0->SetS();
            }
    }

    static void AddCurvatureConstraints(MeshType & mesh,const ScalarType &thr=0.5)
    {
        for (size_t i=0;i<mesh.face.size();i++)
            if (mesh.face[i].Q()>thr)mesh.face[i].SetS();
    }

    //hard constraints have selected face
    static void CollectHardConstraints( MeshType & mesh,
                                    Eigen::VectorXi &HardI,
                                    Eigen::MatrixXd &HardD,
                                    SmoothMethod SMethod,
                                    int Ndir)
    {
        //count number of hard constraints
        int numS=vcg::tri::UpdateSelection<MeshType>::FaceCount(mesh);
        HardI=Eigen::MatrixXi(numS,1);
        if ((Ndir==2)||(SMethod==SMMiq))
            HardD=Eigen::MatrixXd(numS,3);
        else
            HardD=Eigen::MatrixXd(numS,6);
        //then update them
        int curr_index=0;
        for (size_t i=0;i<mesh.face.size();i++)
        {
            if (!mesh.face[i].IsS())continue;

            CoordType dir=mesh.face[i].PD1();
            dir.Normalize();

            HardI(curr_index,0)=i;

            HardD(curr_index,0)=dir.X();
            HardD(curr_index,1)=dir.Y();
            HardD(curr_index,2)=dir.Z();
            if ((Ndir==4)&&(SMethod==SMNPoly))
            {
                dir=mesh.face[i].PD2();
                HardD(curr_index,3)=dir.X();
                HardD(curr_index,4)=dir.Y();
                HardD(curr_index,5)=dir.Z();
            }
            curr_index++;

        }

    }

    //hard constraints have selected face
    static void CollectSoftConstraints( MeshType & mesh,
                                        Eigen::VectorXi &SoftI,
                                        Eigen::MatrixXd &SoftD,
                                        Eigen::VectorXd &SoftW)
    {
        //count number of soft constraints
        int numS=vcg::tri::UpdateSelection<MeshType>::FaceCount(mesh);
        numS=mesh.fn-numS;
        //allocate eigen matrix
        SoftI=Eigen::MatrixXi(numS,1);
        SoftD=Eigen::MatrixXd(numS,3);
        SoftW=Eigen::MatrixXd(numS,1);

        //then update them
        int curr_index=0;
        for (size_t i=0;i<mesh.face.size();i++)
        {
            if (mesh.face[i].IsS())continue;

            CoordType dir=mesh.face[i].PD1();
            dir.Normalize();

            SoftI(curr_index,0)=i;

            SoftD(curr_index,0)=dir.X();
            SoftD(curr_index,1)=dir.Y();
            SoftD(curr_index,2)=dir.Z();

            SoftW(curr_index,0)=mesh.face[i].Q();

            curr_index++;

        }
    }

    static void SmoothMIQ(MeshType &mesh,
                          Eigen::VectorXi &HardI,   //hard constraints index
                          Eigen::MatrixXd &HardD,   //hard directions
                          Eigen::VectorXi &SoftI,   //soft constraints
                          Eigen::MatrixXd &SoftD,   //weight of soft constraints
                          Eigen::VectorXd &SoftW,   //soft directions
                          ScalarType alpha_soft,
                          int Ndir)
    {

        assert((Ndir==2)||(Ndir==4));
        Eigen::MatrixXi F;
        Eigen::MatrixXd V;

        MeshToMatrix<MeshType>::GetTriMeshData(mesh,F,V);

        Eigen::MatrixXd output_field;
        Eigen::VectorXd output_sing;

        igl::nrosy(V,F,HardI,HardD,SoftI,SoftW,SoftD,Ndir,alpha_soft,output_field,output_sing);

        //finally update the principal directions
        for (size_t i=0;i<mesh.face.size();i++)
        {
            CoordType dir1;
            dir1[0]=output_field(i,0);
            dir1[1]=output_field(i,1);
            dir1[2]=output_field(i,2);

            dir1.Normalize();
            CoordType dir2=mesh.face[i].N()^dir1;
            dir2.Normalize();

            ScalarType Norm1=mesh.face[i].PD1().Norm();
            ScalarType Norm2=mesh.face[i].PD2().Norm();

            mesh.face[i].PD1()=dir1*Norm1;
            mesh.face[i].PD2()=dir2*Norm2;
        }
    }

    static void SmoothNPoly(MeshType &mesh,
                            Eigen::VectorXi &HardI,   //hard constraints index
                            Eigen::MatrixXd &HardD,   //hard directions
                            int Ndir)
    {
        assert((Ndir==2)||(Ndir==4));

        Eigen::MatrixXi F;
        Eigen::MatrixXd V;

        MeshToMatrix<MeshType>::GetTriMeshData(mesh,F,V);

        Eigen::MatrixXd output_field;
        Eigen::VectorXd output_sing;

        igl::n_polyvector(V,F,HardI,HardD,output_field);

        //finally update the principal directions
        for (size_t i=0;i<mesh.face.size();i++)
        {
            CoordType dir1;
            dir1[0]=output_field(i,0);
            dir1[1]=output_field(i,1);
            dir1[2]=output_field(i,2);

            dir1.Normalize();
            CoordType dir2=mesh.face[i].N()^dir1;
            dir2.Normalize();

            ScalarType Norm1=mesh.face[i].PD1().Norm();
            ScalarType Norm2=mesh.face[i].PD2().Norm();

            mesh.face[i].PD1()=dir1*Norm1;
            mesh.face[i].PD2()=dir2*Norm2;
        }
    }

    static void PickRandomDir(MeshType &mesh,
                              int &indexF,
                              CoordType &Dir)
    {
        indexF=rand()%mesh.fn;
        FaceType *currF=&mesh.face[indexF];
        CoordType dirN=currF->N();
        dirN.Normalize();
        Dir=CoordType(1,0,0);
        if (fabs(Dir*dirN)>0.9)
            Dir=CoordType(0,1,0);
        if (fabs(Dir*dirN)>0.9)
            Dir=CoordType(0,0,1);

        Dir=dirN^Dir;
        Dir.Normalize();
    }

public:

    struct SmoothParam
    {
        //the 90° rotation independence while smoothing the direction field
        int Ndir;
        //the weight of curvature if doing the smoothing keeping the field close to the original one
        ScalarType alpha_curv;
        //align the field to border or not
        bool align_borders;
        //threshold to consider some edge as sharp feature and to use as hard constraint (0, not use)
        ScalarType sharp_thr;
        //threshold to consider some edge as high curvature anisotropyand to use as hard constraint (0, not use)
        ScalarType curv_thr;
        //the method used to smooth MIQ or "Designing N-PolyVector Fields with Complex Polynomials"
        SmoothMethod SmoothM;
        //the number of faces of the ring used ot esteem the curvature
        int curvRing;
        //this are additional hard constraints
        std::vector<std::pair<int,CoordType> > AddConstr;

        SmoothParam()
        {
            Ndir=4;
            curvRing=2;
            alpha_curv=0.0;

            align_borders=false;

            SmoothM=SMMiq;

            sharp_thr=0.0;
            curv_thr=0.4;
        }

    };

    static void SelectConstraints(MeshType &mesh,SmoothParam &SParam)
    {
        //clear all selected faces
        vcg::tri::UpdateFlags<MeshType>::FaceClear(mesh);

        //add curvature hard constraints
        //ScalarType Ratio=mesh.bbox.Diag()*0.01;

        if (SParam.curv_thr>0)
            AddCurvatureConstraints(mesh,SParam.curv_thr);///Ratio);

         //add alignment to sharp features
        if (SParam.sharp_thr>0)
            AddSharpEdgesConstraints(mesh,SParam.sharp_thr);

        //add border constraints
        if (SParam.align_borders)
            AddBorderConstraints(mesh);

        //aff final constraints
        for (int i=0;i<SParam.AddConstr.size();i++)
        {
            int indexF=SParam.AddConstr[i].first;
            CoordType dir=SParam.AddConstr[i].second;
            mesh.face[indexF].PD1()=dir;
            mesh.face[indexF].PD2()=mesh.face[indexF].N()^dir;
            mesh.face[indexF].PD1().Normalize();
            mesh.face[indexF].PD2().Normalize();
            mesh.face[indexF].SetS();
        }
    }

    static void GloballyOrient(MeshType &mesh)
    {
        vcg::tri::CrossField<MeshType>::MakeDirectionFaceCoherent(mesh,true);
    }


    static void InitByCurvature(MeshType & mesh,
                                int Nring,
                                bool UpdateFaces=true)
    {

        tri::RequirePerVertexCurvatureDir(mesh);

        Eigen::MatrixXi F;
        Eigen::MatrixXd V;

        Eigen::MatrixXd PD1,PD2,PV1,PV2;
        MeshToMatrix<MeshType>::GetTriMeshData(mesh,F,V);
        igl::principal_curvature(V,F,PD1,PD2,PV1,PV2,Nring);

        //then copy curvature per vertex
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            mesh.vert[i].PD1()=CoordType(PD1(i,0),PD1(i,1),PD1(i,2));
            mesh.vert[i].PD2()=CoordType(PD2(i,0),PD2(i,1),PD2(i,2));
            mesh.vert[i].K1()=PV1(i,0);
            mesh.vert[i].K2()=PV2(i,0);
        }
        if (!UpdateFaces)return;
        vcg::tri::CrossField<MeshType>::SetFaceCrossVectorFromVert(mesh);
        InitQualityByAnisotropyDir(mesh);
    }

    static void SmoothDirections(MeshType &mesh,
                                 int Ndir,
                                 SmoothMethod SMethod=SMNPoly,
                                 bool HardAsS=true,
                                 ScalarType alphaSoft=0)
    {

        Eigen::VectorXi HardI;   //hard constraints
        Eigen::MatrixXd HardD;   //hard directions
        Eigen::VectorXi SoftI;   //soft constraints
        Eigen::VectorXd SoftW;   //weight of soft constraints
        Eigen::MatrixXd SoftD;   //soft directions

        if (HardAsS)
            CollectHardConstraints(mesh,HardI,HardD,SMethod,Ndir);

        //collect soft constraints , miw only one that allows for soft constraints
        if ((alphaSoft>0)&&(SMethod==SMMiq))
          CollectSoftConstraints(mesh,SoftI,SoftD,SoftW);

        //add some hard constraints if are not present
        int numC=3;
        if ((SoftI.rows()==0)&&(HardI.rows()==0))
        {
            printf("Add Forced Constraints \n");
            fflush(stdout);
            HardI=Eigen::MatrixXi(numC,1);

            if ((Ndir==4)&&(SMethod==SMNPoly))
                HardD=Eigen::MatrixXd(numC,6);
            else
                HardD=Eigen::MatrixXd(numC,3);

            for (int i=0;i<numC;i++)
            {
                CoordType Dir;
                int indexF;
                PickRandomDir(mesh,indexF,Dir);

                HardI(i,0)=indexF;
                HardD(i,0)=Dir.X();
                HardD(i,1)=Dir.Y();
                HardD(i,2)=Dir.Z();

                if ((Ndir==4)&&(SMethod==SMNPoly))
                {
                    CoordType Dir1=mesh.face[indexF].N()^Dir;
                    Dir1.Normalize();
                    HardD(i,3)=Dir1.X();
                    HardD(i,4)=Dir1.Y();
                    HardD(i,5)=Dir1.Z();
                }
            }
        }

         //finally smooth
        if (SMethod==SMMiq)
         SmoothMIQ(mesh,HardI,HardD,SoftI,SoftD,SoftW,alphaSoft,Ndir);
        else
        {
            assert(SMethod==SMNPoly);
            SmoothNPoly(mesh,HardI,HardD,Ndir);
        }
    }


    static void SmoothDirections(MeshType &mesh,SmoothParam SParam)
    {
        //for the moment only cross and line field

        //initialize direction by curvature if needed
        if ((SParam.alpha_curv>0)||
             (SParam.sharp_thr>0)||
             (SParam.curv_thr>0))
            InitByCurvature(mesh,SParam.curvRing);

        SelectConstraints(mesh,SParam);
        //then do the actual smooth
        SmoothDirections(mesh,SParam.Ndir,SParam.SmoothM,true,SParam.alpha_curv);
    }

};

} // end namespace tri
} // end namespace vcg
#endif // SMOOTHER_FIELD_H
