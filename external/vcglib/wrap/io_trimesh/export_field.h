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
#ifndef __VCGLIB_EXPORTERFIELD
#define __VCGLIB_EXPORTERFIELD

#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>

namespace vcg {
namespace tri {
namespace io {

/** 
This class encapsulate a filter for saving field formats
*/
template <class MeshType>
class ExporterFIELD
{
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;

public:
    
    ///load a field on the mesh, it could be a vfield file (per vertex)
    ///or an ffield file (per face)
    static void SaveFaceFIELD(MeshType &mesh,
                              const char *path)
    {
        
        FILE *f = fopen(path,"wt");
        //if (!f)return false;
//            char word[512]; word[0]=0;
//            fscanf(f,"%s",word);
//            char c=0;
//            if (word[0]=='#') {
//                // skip comment line
//                while (fscanf(f,"%c",&c)!=EOF) if (c=='\n') break;
//            }
//            else
//            {
//                return false;
//            }
            int nf = mesh.fn;//-1;
            fprintf(f,"# frame generated with VCG \n");
            fprintf(f,"target frame \n");
            fprintf(f,"%d\n",nf);
//            
//            if( per_vertex && !HasPerVertexCurvatureDir(mesh)) throw vcg::MissingComponentException("PerVertexCurvatureDir");
//            if(!per_vertex && !HasPerFaceCurvatureDir(mesh))   throw vcg::MissingComponentException("PerFaceCurvatureDir");
            if (!HasPerFaceCurvatureDir(mesh))
                throw vcg::MissingComponentException("PerFaceCurvatureDir");
            
            fprintf(f,"k1	 k2	 k1v_x	 k1v_y	 k1v_z	 k2v_x	 k2v_y	 k2v_z\n");

            for (int i=0; i<nf; i++){
                vcg::Point3<float> u;
                u.Import(mesh.face[i].PD1());
                vcg::Point3<float> v;
                v.Import(mesh.face[i].PD2());
                
                fprintf(f,"1 1 %f %f %f %f %f %f\n",
                        (u.X()),(u.Y()),(u.Z()),
                        (v.X()),(v.Y()),(v.Z()));
            }
        fclose(f);
    }

    
    ///Save a 4 rosy format file as used by
    ///Interactive Visualization of Rotational Symmetry Fields on Surfaces
    ///Jonathan Palacios and Eugene Zhang
    static void Save4ROSY(MeshType &mesh,
                        const char *path)
    {
        FILE *f = fopen(path,"wt");
        fprintf(f,"%d\n",mesh.vn);
        fprintf(f,"4\n");
        for (unsigned int i=0;i<mesh.vert.size();i++)
        {
            float dirX=(float)mesh.vert[i].PD1().X();
            float dirY=(float)mesh.vert[i].PD1().Y();
            float dirZ=(float)mesh.vert[i].PD1().Z();
            fprintf(f,"%f %f %f \n",dirX,dirY,dirZ);

        }
        fclose(f);
    }

    //Save a 4 rosy format file as pair of angles
    static void Save2AngleFace(MeshType &mesh,
                              const char *path)
    {
        FILE *f = fopen(path,"wt");
        fprintf(f,"#%d param_field\n",mesh.fn);
        for (unsigned int i=0;i<mesh.face.size();i++)
        {
            ScalarType alpha1,alpha2;
            CoordType PD1Test=mesh.face[i].PD1();
            CoordType PD2Test=mesh.face[i].PD2();
            vcg::tri::CrossField<MeshType>::CrossFieldToAngles(mesh.face[i],alpha1,alpha2,1);
            fprintf(f,"%d %f %f \n",i,alpha1,alpha2);
        }
        fclose(f);
    }

}; // end class



} // end namespace tri
} // end namespace io
} // end namespace vcg

#endif

