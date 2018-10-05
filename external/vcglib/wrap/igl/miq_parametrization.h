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

#ifndef __MIQ_PARAMETRIZATION_H
#define __MIQ_PARAMETRIZATION_H

//igl stuff
#include <igl/cross_field_missmatch.h>
#include <igl/line_field_missmatch.h>
#include <igl/comb_line_field.h>
#include <igl/cut_mesh_from_singularities.h>
#include <igl/find_cross_field_singularities.h>
#include <igl/compute_frame_field_bisectors.h>
#include <igl/comiso/miq.h>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>

namespace vcg {
namespace tri {
template < class MeshType>         // Classe templatata su Tipo Mesh
class MiQParametrizer
{

    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::VertexType PolyVertexType;

public:
    struct MIQParameters
    {
        //the gradient of the parametrization 1 is the bb diagonal, as big is the gradient as small are the quads
        double gradient;
        //do the rounding or not across cuts... set to true to get a quadrangulation
        bool doRound;
        //do the round at once for each stiffness iteration or do it gradually.. gradually is more stable but much mor slow
        bool directRound;
        //the stiffness increment for ach iteration
        double stiffness;
        //the maximum number ofstiffness iteration to avoid folds
        int stiffness_iter;
        //local iteration to round integer variables ofr each stiffness round
        int local_iter;
        //this bool multiply the gradiens U or V separately times a constant to make hexagons
        bool hexaLine;
        //the threshold of normal dot product to consider a crease
        double crease_thr;
        //number of 90° rotation independence (4 for quad meshing, 2 to obtain a quad meshing for hexagonalization)
        int Ndir;
        //round or not the singularities
        bool round_singularities;
        //use the crease edges as feature or not
        bool crease_as_feature;
        //true if roound selected vert
        bool round_selected;

        MIQParameters()
        {
            gradient=80;
            doRound=true;
            directRound=true;
            round_singularities=true;
            crease_as_feature=false;
            round_selected=true;
            stiffness=5;
            stiffness_iter=10;
            local_iter=5;
            Ndir=4;
            crease_thr=0.2;
            hexaLine=false;
        }
    };


    static void GetFeatureLines(MeshType &trimesh,
                                std::vector<std::vector<int> > &feature_lines)
    {
        feature_lines.clear();

        for (size_t i=0;i<trimesh.face.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                //if (!trimesh.face[i].IsB(j))continue;
                if (!trimesh.face[i].IsCrease(j))continue;

                feature_lines.push_back(std::vector<int>());
                feature_lines.back().push_back(i);
                feature_lines.back().push_back(j);

            }
        }
    }

private:



    static void CrossFieldParam(MeshType &trimesh,
                                MIQParameters &MiqP)
    {

        Eigen::MatrixXi F;
        Eigen::MatrixXd V;
        vcg::tri::MeshToMatrix<MeshType>::GetTriMeshData(trimesh,F,V);

        //then get the principal directions
        Eigen::MatrixXd X1,X2;
        X1=Eigen::MatrixXd(trimesh.FN(), 3);
        for (size_t i=0;i<trimesh.face.size();i++)
        {
            CoordType Dir1=trimesh.face[i].PD1();
            Dir1.Normalize();
            for (int j=0;j<3;j++)
            {
                X1(i,j)=Dir1[j];
            }
        }

        Eigen::MatrixXd B1,B2,B3;
        igl::local_basis(V,F,B1,B2,B3);
        X2 = igl::rotate_vectors(X1, Eigen::VectorXd::Constant(1,M_PI/2), B1, B2);

        Eigen::MatrixXd UV;
        Eigen::MatrixXi FUV;

        //ScalarType gradsize=trimesh.bbox.Diag()*2;

        std::vector<std::vector<int> > hard_features;

        if (MiqP.crease_as_feature)
            GetFeatureLines(trimesh,hard_features);

        std::vector<int> extra_round;

        if (MiqP.round_selected)
        {
            for (int i=0;i<trimesh.vert.size();i++)
            {
                if (!trimesh.vert[i].IsS())continue;
                extra_round.push_back(i);
            }
        }

        igl::miq(V,F,X1,X2,UV,FUV,MiqP.gradient,MiqP.stiffness,MiqP.directRound,
                 MiqP.stiffness_iter,MiqP.local_iter,MiqP.doRound,MiqP.round_singularities,
                 extra_round,hard_features);

        // then copy UV
        for (size_t i=0;i<trimesh.face.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                int index=FUV(i,j);
                trimesh.face[i].WT(j).P()[0]=UV(index,0);
                trimesh.face[i].WT(j).P()[1]=UV(index,1);
            }
        }
    }

    static void LineFieldParam(MeshType &trimesh,
                               MIQParameters &MiqP)
    {

        Eigen::MatrixXi F;
        Eigen::MatrixXd V;
        vcg::tri::MeshToMatrix<MeshType>::GetTriMeshData(trimesh,F,V);

        //then get the principal directions
        Eigen::MatrixXd X1,X2;
        X1=Eigen::MatrixXd(trimesh.FN(), 3);
        for (size_t i=0;i<trimesh.face.size();i++)
        {
            CoordType Dir1=trimesh.face[i].PD1();
            Dir1.Normalize();
            for (int j=0;j<3;j++)
            {
                X1(i,j)=Dir1[j];
            }
        }

        Eigen::MatrixXd B1,B2,B3;
        igl::local_basis(V,F,B1,B2,B3);
        X2 = igl::rotate_vectors(X1, Eigen::VectorXd::Constant(1,M_PI/2), B1, B2);

        // Bisector field
        Eigen::MatrixXd BIS1, BIS2;

        // Combed bisector
        Eigen::MatrixXd BIS1_combed, BIS2_combed;

        // Per-corner, integer mismatches
        Eigen::MatrixXi MMatch;

        // Field singularities
        Eigen::VectorXi isSingularity, singularityIndex;

        // Per corner seams
        Eigen::MatrixXi Seams;

        // Combed field
        Eigen::MatrixXd X1_combed, X2_combed;


        // Global parametrization (with seams)
        Eigen::MatrixXd UV_seams;
        Eigen::MatrixXi FUV_seams;

        // Global parametrization
        Eigen::MatrixXd UV;
        Eigen::MatrixXi FUV;

        // Always work on the bisectors, it is more general
        igl::compute_frame_field_bisectors(V, F, X1, X2, BIS1, BIS2);

        // Comb the field, implicitly defining the seams
        igl::comb_line_field(V, F, BIS1, BIS1_combed);
        igl::local_basis(V,F,B1,B2,B3);
        BIS2_combed = igl::rotate_vectors(BIS1_combed, Eigen::VectorXd::Constant(1,M_PI/2), B1, B2);

        // Find the integer mismatches
        igl::line_field_missmatch(V, F, BIS1_combed, true, MMatch);

        // Find the singularities
        igl::find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex);

        // Cut the mesh, duplicating all vertices on the seams
        //igl::cut_mesh_from_singularities(V, F, MMatch, isSingularity, singularityIndex, Seams);
        igl::cut_mesh_from_singularities(V, F, MMatch,Seams);

        // Comb the frame-field accordingly
        igl::comb_frame_field(V, F, X1, X2, BIS1_combed, BIS2_combed, X1_combed, X2_combed);

        std::vector<std::vector<int> > hard_features;


        std::vector<int> extra_round;
        //collect extra vertex selected that need to be rounded
        if (MiqP.round_selected)
        {
            for (int i=0;i<trimesh.vert.size();i++)
            {
                if (!trimesh.vert[i].IsS())continue;
                extra_round.push_back(i);
            }
        }
        if (MiqP.crease_as_feature)
            GetFeatureLines(trimesh,hard_features);

        //scale gradient if needed
        ScalarType sqrt3=1.732050807568877;

        ScalarType GradX=0.5;
        ScalarType GradY=1;

        if (MiqP.hexaLine)
        {
            for (int i=0;i<X1_combed.rows();i++)
            {
                X1_combed(i)*=GradX;//*ScaleFact;
                X2_combed(i)*=GradY;//*ScaleFact;
            }
        }

//        igl::miq(V,F,X1_combed,X2_combed,BIS1_combed,BIS2_combed,
//                 MMatch,isSingularity,singularityIndex,Seams,
//                 UV,FUV,MiqP.gradient,MiqP.stiffness,MiqP.directRound,
//                 MiqP.stiffness_iter,MiqP.local_iter,MiqP.doRound,MiqP.round_singularities,extra_round,hard_features);
        igl::miq(V,F,X1_combed,X2_combed,
                 UV,FUV,MiqP.gradient,MiqP.stiffness,MiqP.directRound,
                 MiqP.stiffness_iter,MiqP.local_iter,MiqP.doRound,MiqP.round_singularities,extra_round,hard_features);

        // then copy UV
        for (size_t i=0;i<trimesh.face.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                int index=FUV(i,j);
                trimesh.face[i].WT(j).P()[0]=UV(index,0);//*4;
                trimesh.face[i].WT(j).P()[1]=UV(index,1);//*2;
            }
        }
    }

public:

    static void SetCreases(MeshType & mesh,
                           const ScalarType &thr=0.2,
                           bool setBorder=true)
    {
        for (size_t i=0;i<mesh.face.size();i++)
            for (int j=0;j<mesh.face[i].VN();j++)
            {
                FaceType *f0=&mesh.face[i];
                f0->ClearCrease(j);
            }


        for (size_t i=0;i<mesh.face.size();i++)
            for (int j=0;j<mesh.face[i].VN();j++)
            {
                FaceType *f0=&mesh.face[i];
                FaceType *f1=f0->FFp(j);

                if (f0==f1){f0->SetCrease(j);continue;}

                CoordType N0=f0->N();
                CoordType N1=f1->N();
                if ((N0*N1)>thr)continue;
                f0->SetCrease(j);

            }
    }

    static void MIQParametrize(MeshType &trimesh,
                               MIQParameters &MiqP)
    {
        if (MiqP.crease_as_feature)
            SetCreases(trimesh,MiqP.crease_thr);

        if (MiqP.Ndir==4)
            CrossFieldParam(trimesh,MiqP);
        else
            LineFieldParam(trimesh,MiqP);
    }
};

} // end namespace tri
} // end namespace vcg
#endif // VORO_CLUSTER_H
