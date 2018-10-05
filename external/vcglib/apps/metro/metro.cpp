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
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.23  2007/05/04 16:50:23  ganovelli
added plus types version (#ifdef _PLUS_TYPES_ to use it ).

Revision 1.22  2006/10/25 12:40:19  fiorin
Added possibility to use Octree as search structure:

Revision 1.21  2006/05/03 21:22:39  cignoni
added missing Include

Revision 1.20  2006/04/20 08:30:24  cignoni
small GCC compiling issues

Revision 1.19  2006/03/27 04:17:07  cignoni
moved to generic export.h

Revision 1.18  2006/01/10 13:20:40  cignoni
Changed ply::PlyMask to io::Mask

Revision 1.17  2005/10/02 23:11:00  cignoni
Version 4.06, Added possibility of using three different search structures UG Hash and AABB

Revision 1.16  2005/09/16 11:52:14  cignoni
removed wrong %v in vertex number printing

Revision 1.15  2005/04/04 10:36:36  cignoni
Release 4.05
Added saving of Error Histogram

Revision 1.14  2005/01/26 22:45:34  cignoni
Release 4.04
final updates for gcc compiling issues

Revision 1.13  2005/01/24 15:46:48  cignoni
Release 4.04
Moved to the library core the code for computing min distance froma a point to a mesh using a uniform grid.
Slightly faster.

Revision 1.12  2005/01/03 11:28:52  cignoni
Release 4.03
Better ply compatibility, and improved error reporting

Revision 1.11  2004/11/29 09:07:04  cignoni
Release 4.02
removed bug in printing Hausdorf distance,
removed bug in command line parsing,
upgraded import mesh library to support off format

Revision 1.10  2004/09/21 23:52:50  cignoni
Release 4.01

Revision 1.9  2004/09/20 16:29:08  ponchio
Minimal changes.

Revision 1.8  2004/09/20 15:17:28  cignoni
Removed bug in displays msec and better usage messages

Revision 1.7  2004/09/09 22:59:15  cignoni
Removed many small warnings

Revision 1.6  2004/07/15 00:15:16  cignoni
inflate -> offset

Revision 1.5  2004/06/24 09:08:31  cignoni
Official Release of Metro 4.00

Revision 1.4  2004/05/14 13:53:12  ganovelli
GPL  added

****************************************************************************/

// -----------------------------------------------------------------------------------------------

// standard libraries
#include <time.h>


#include <vcg/math/histogram.h>
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/component_ep.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/complex/algorithms/update/component_ep.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include "sampling.h"

using namespace std;
// project definitions.
// error messages

#define MSG_ERR_MESH_LOAD               "error loading the input meshes.\n"
#define MSG_ERR_INVALID_OPTION          "unable to parse option '%s'\n"
#define MSG_ERR_FILE_OPEN               "unable to open the output file.'n"
#define MSG_ERR_UNKNOWN_FORMAT          "unknown file format '%s'.\n"

// global constants
#define NO_SAMPLES_PER_FACE             10
#define N_SAMPLES_EDGE_TO_FACE_RATIO    0.1
#define BBOX_FACTOR                     0.1
#define INFLATE_PERCENTAGE			    0.02
#define MIN_SIZE					    125		/* 125 = 5^3 */
#define N_HIST_BINS                     256
#define PRINT_EVERY_N_ELEMENTS          1000


class CFace;
class CVertex;
struct UsedTypes:public vcg::UsedTypes< vcg::Use<CFace>::AsFaceType, vcg::Use<CVertex>::AsVertexType>{};
class CVertex   : public vcg::Vertex<UsedTypes,vcg::vertex::Coord3d,vcg::vertex::Qualityf,vcg::vertex::Normal3d,vcg::vertex::Color4b,vcg::vertex::BitFlags> {};
class CFace     : public vcg::Face< UsedTypes,vcg::face::VertexRef, vcg::face::Normal3d, vcg::face::EdgePlane,vcg::face::Color4b,vcg::face::Mark,vcg::face::BitFlags> {};
class CMesh     : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace> > {};


// -----------------------------------------------------------------------------------------------

using namespace vcg;


////////////////// Command line Flags and parameters
bool NumberOfSamples                = false;
bool SamplesPerAreaUnit             = false;
bool CleaningFlag=false;
// -----------------------------------------------------------------------------------------------

void Usage()
{
  printf("\nUsage:  "\
                                        "metro file1 file2 [opt]\n"\
                                        "Where opt can be:\n"\
                                        "  -v         disable vertex sampling\n"\
                                        "  -e         disable edge sampling\n"\
                                        "  -f         disable face sampling\n"\
                                        "  -u         ignore unreferred vertices\n"\
                                        "  -sx        set the face sampling mode\n"\
                                        "             where x can be:\n"\
                                        "              -s0  montecarlo sampling\n"\
                                        "              -s1  subdivision sampling\n"\
                                        "              -s2  similar triangles sampling (Default)\n"\
                                        "  -n#        set the required number of samples (overrides -A)\n"\
                                        "  -a#        set the required number of samples per area unit (overrides -N)\n"\
                                        "  -c         save a mesh with error as per-vertex colour and quality\n"\
                                        "  -C # #     Set the min/max values used for color mapping\n"\
                                        "  -L         Remove duplicated and unreferenced vertices before processing\n"\
                                        "  -h         write files with histograms of error distribution\n"\
                                        "  -G         Use a static Uniform Grid as Search Structure (default)\n"\
																				"  -O         Use an octree as a Search Structure\n"\
                                        "  -A         Use an AxisAligned Bounding Box Tree as Search Structure\n"\
                                        "  -H         Use an Hashed Uniform Grid as Search Structure\n"\
                                        "\n"
                                        "Default options are to sample vertexes, edge and faces by taking \n"
                                        "a number of samples that is approx. 10x the face number.\n"
                                        );
  exit(-1);
}


// simple aux function that compute the name for the file containing the stored computations
std::string SaveFileName(const std::string &filename)
{
 int pos=filename.find_last_of('.',filename.length());
 std::string fileout=filename.substr(0,pos)+"_metro.ply";
 return fileout;
}


// Open Mesh
void OpenMesh(const char *filename, CMesh &m)
{
  int err = tri::io::Importer<CMesh>::Open(m,filename);
  if(err) {
      printf("Error in reading %s: '%s'\n",filename,tri::io::Importer<CMesh>::ErrorMsg(err));
      if(tri::io::Importer<CMesh>::ErrorCritical(err)) exit(-1);
    }
  printf("read mesh `%s'\n", filename);
  if(CleaningFlag){
      int dup = tri::Clean<CMesh>::RemoveDuplicateVertex(m);
      int unref =  tri::Clean<CMesh>::RemoveUnreferencedVertex(m);
      printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,filename);
  }
}


int main(int argc, char**argv)
{
    CMesh                 S1, S2;
    float                ColorMin=0, ColorMax=0;
    double                dist1_max, dist2_max;
    unsigned long         n_samples_target, elapsed_time;
    double								n_samples_per_area_unit;
    int                   flags;

    // print program info
    printf("-------------------------------\n"
           "         Metro V.4.07 \n"
           "     http://vcg.isti.cnr.it\n"
           "   release date: "__DATE__"\n"
           "-------------------------------\n\n");

    if(argc <= 2)    Usage();
    // default parameters
    flags = SamplingFlags::VERTEX_SAMPLING |
          SamplingFlags::EDGE_SAMPLING |
          SamplingFlags::FACE_SAMPLING |
          SamplingFlags::SIMILAR_SAMPLING;

    // parse command line.
	  for(int i=3; i < argc;)
    {
      if(argv[i][0]=='-')
        switch(argv[i][1])
      {
        case 'h' : flags |= SamplingFlags::HIST; break;
        case 'v' : flags &= ~SamplingFlags::VERTEX_SAMPLING; break;
        case 'e' : flags &= ~SamplingFlags::EDGE_SAMPLING; break;
        case 'f' : flags &= ~SamplingFlags::FACE_SAMPLING; break;
        case 'u' : flags |= SamplingFlags::INCLUDE_UNREFERENCED_VERTICES; break;
        case 's'   :
          switch(argv[i][2])
          {
            case '0':  flags = (flags | SamplingFlags::MONTECARLO_SAMPLING  ) & (~ SamplingFlags::NO_SAMPLING );break;
            case '1':  flags = (flags | SamplingFlags::SUBDIVISION_SAMPLING ) & (~ SamplingFlags::NO_SAMPLING );break;
            case '2':  flags = (flags | SamplingFlags::SIMILAR_SAMPLING     ) & (~ SamplingFlags::NO_SAMPLING );break;
            default  :  printf(MSG_ERR_INVALID_OPTION, argv[i]);
              exit(0);
          }
          break;
        case 'n':  NumberOfSamples       = true;     n_samples_target        = (unsigned long) atoi(&(argv[i][2]));          break;
        case 'a':  SamplesPerAreaUnit    = true;     n_samples_per_area_unit = (unsigned long) atoi(&(argv[i][2])); break;
        case 'c':  flags |= SamplingFlags::SAVE_ERROR;   break;
        case 'L':  CleaningFlag=true; break;
        case 'C':  ColorMin=float(atof(argv[i+1])); ColorMax=float(atof(argv[i+2])); i+=2; break;
        case 'A':  flags |= SamplingFlags::USE_AABB_TREE;   printf("Using AABB Tree as search structure\n");           break;
        case 'G':  flags |= SamplingFlags::USE_STATIC_GRID; printf("Using static uniform grid as search structure\n"); break;
        case 'H':  flags |= SamplingFlags::USE_HASH_GRID;   printf("Using hashed uniform grid as search structure\n"); break;
				case 'O':  flags |= SamplingFlags::USE_OCTREE;      printf("Using octree as search structure\n");              break;
        default  :  printf(MSG_ERR_INVALID_OPTION, argv[i]);
          exit(0);
      }
      i++;
    }

		if(!(flags & SamplingFlags::USE_HASH_GRID) && !(flags & SamplingFlags::USE_AABB_TREE) && !(flags & SamplingFlags::USE_OCTREE))
       flags |= SamplingFlags::USE_STATIC_GRID;

    // load input meshes.
    OpenMesh(argv[1],S1);
    OpenMesh(argv[2],S2);

    string S1NewName=SaveFileName(argv[1]);
    string S2NewName=SaveFileName(argv[2]);

    if(!NumberOfSamples && !SamplesPerAreaUnit)
    {
        NumberOfSamples = true;
        n_samples_target = 10 * max(S1.fn,S2.fn);// take 10 samples per face
    }

    // compute face information
        tri::UpdateComponentEP<CMesh>::Set(S1);
        tri::UpdateComponentEP<CMesh>::Set(S2);

	// set bounding boxes for S1 and S2
		tri::UpdateBounding<CMesh>::Box(S1);
		tri::UpdateBounding<CMesh>::Box(S2);

    // set Bounding Box.
		Box3<CMesh::ScalarType>    bbox, tmp_bbox_M1=S1.bbox, tmp_bbox_M2=S2.bbox;
    bbox.Add(S1.bbox);
    bbox.Add(S2.bbox);
		bbox.Offset(bbox.Diag()*0.02);
	  S1.bbox = bbox;
	  S2.bbox = bbox;

    // initialize time info.
    int t0=clock();

    Sampling<CMesh> ForwardSampling(S1,S2);
    Sampling<CMesh> BackwardSampling(S2,S1);

    // print mesh info.
    printf("Mesh info:\n");
    printf(" M1: '%s'\n\tvertices  %7i\n\tfaces     %7i\n\tarea      %12.4f\n", argv[1], S1.vn, S1.fn, ForwardSampling.GetArea());
    printf("\tbbox (%7.4f %7.4f %7.4f)-(%7.4f %7.4f %7.4f)\n", tmp_bbox_M1.min[0], tmp_bbox_M1.min[1], tmp_bbox_M1.min[2], tmp_bbox_M1.max[0], tmp_bbox_M1.max[1], tmp_bbox_M1.max[2]);
    printf("\tbbox diagonal %f\n", (float)tmp_bbox_M1.Diag());
    printf(" M2: '%s'\n\tvertices  %7i\n\tfaces     %7i\n\tarea      %12.4f\n", argv[2], S2.vn, S2.fn, BackwardSampling.GetArea());
    printf("\tbbox (%7.4f %7.4f %7.4f)-(%7.4f %7.4f %7.4f)\n", tmp_bbox_M2.min[0], tmp_bbox_M2.min[1], tmp_bbox_M2.min[2], tmp_bbox_M2.max[0], tmp_bbox_M2.max[1], tmp_bbox_M2.max[2]);
    printf("\tbbox diagonal %f\n", (float)tmp_bbox_M2.Diag());

    // Forward distance.
    printf("\nForward distance (M1 -> M2):\n");
    ForwardSampling.SetFlags(flags);
    if(NumberOfSamples)
    {
        ForwardSampling.SetSamplesTarget(n_samples_target);
        n_samples_per_area_unit = ForwardSampling.GetNSamplesPerAreaUnit();
    }
    else
    {
        ForwardSampling.SetSamplesPerAreaUnit(n_samples_per_area_unit);
        n_samples_target = ForwardSampling.GetNSamplesTarget();
    }
    printf("target # samples      : %lu\ntarget # samples/area : %f\n", n_samples_target, n_samples_per_area_unit);
    ForwardSampling.Hausdorff();
    dist1_max  = ForwardSampling.GetDistMax();
    printf("\ndistances:\n  max  : %f (%f  wrt bounding box diagonal)\n", (float)dist1_max, (float)dist1_max/bbox.Diag());
    printf("  mean : %f\n", ForwardSampling.GetDistMean());
    printf("  RMS  : %f\n", ForwardSampling.GetDistRMS());
    printf("# vertex samples %9lu\n", ForwardSampling.GetNVertexSamples());
    printf("# edge samples   %9lu\n", ForwardSampling.GetNEdgeSamples());
    printf("# area samples   %9lu\n", ForwardSampling.GetNAreaSamples());
    printf("# total samples  %9lu\n", ForwardSampling.GetNSamples());
    printf("# samples per area unit: %f\n\n", ForwardSampling.GetNSamplesPerAreaUnit());

    // Backward distance.
    printf("\nBackward distance (M2 -> M1):\n");
    BackwardSampling.SetFlags(flags);
    if(NumberOfSamples)
    {
        BackwardSampling.SetSamplesTarget(n_samples_target);
        n_samples_per_area_unit = BackwardSampling.GetNSamplesPerAreaUnit();
    }
    else
    {
        BackwardSampling.SetSamplesPerAreaUnit(n_samples_per_area_unit);
        n_samples_target = BackwardSampling.GetNSamplesTarget();
    }
    printf("target # samples      : %lu\ntarget # samples/area : %f\n", n_samples_target, n_samples_per_area_unit);
    BackwardSampling.Hausdorff();
    dist2_max  = BackwardSampling.GetDistMax();
    printf("\ndistances:\n  max  : %f (%f  wrt bounding box diagonal)\n", (float)dist2_max, (float)dist2_max/bbox.Diag());
    printf("  mean : %f\n", BackwardSampling.GetDistMean());
    printf("  RMS  : %f\n", BackwardSampling.GetDistRMS());
    printf("# vertex samples %9lu\n", BackwardSampling.GetNVertexSamples());
    printf("# edge samples   %9lu\n", BackwardSampling.GetNEdgeSamples());
    printf("# area samples   %9lu\n", BackwardSampling.GetNAreaSamples());
    printf("# total samples  %9lu\n", BackwardSampling.GetNSamples());
    printf("# samples per area unit: %f\n\n", BackwardSampling.GetNSamplesPerAreaUnit());

    // compute time info.
    elapsed_time = clock() - t0;
    int n_total_sample=ForwardSampling.GetNSamples()+BackwardSampling.GetNSamples();
    double mesh_dist_max  = max(dist1_max , dist2_max);

    printf("\nHausdorff distance: %f (%f  wrt bounding box diagonal)\n",(float)mesh_dist_max,(float)mesh_dist_max/bbox.Diag());
    printf("  Computation time  : %d ms\n",(int)(1000.0*elapsed_time/CLOCKS_PER_SEC));
    printf("  # samples/second  : %f\n\n", (float)n_total_sample/((float)elapsed_time/CLOCKS_PER_SEC));

    // save error files.
    if(flags & SamplingFlags::SAVE_ERROR)
    {
      vcg::tri::io::PlyInfo p;
      p.mask|=vcg::tri::io::Mask::IOM_VERTCOLOR | vcg::tri::io::Mask::IOM_VERTQUALITY /* | vcg::ply::PLYMask::PM_VERTQUALITY*/ ;
      //p.mask|=vcg::ply::PLYMask::PM_VERTCOLOR|vcg::ply::PLYMask::PM_VERTQUALITY;
      if(ColorMax!=0 || ColorMin != 0){
        vcg::tri::UpdateColor<CMesh>::PerVertexQualityRamp(S1,ColorMin,ColorMax);
        vcg::tri::UpdateColor<CMesh>::PerVertexQualityRamp(S2,ColorMin,ColorMax);
      }
      tri::io::ExporterPLY<CMesh>::Save( S1,S1NewName.c_str(),true,p);
      tri::io::ExporterPLY<CMesh>::Save( S2,S2NewName.c_str(),true,p);
    }

    // save error files.
    if(flags & SamplingFlags::HIST)
    {
      ForwardSampling.GetHist().FileWrite("forward_result.csv");
      BackwardSampling.GetHist().FileWrite("backward_result.csv");
    }
   return 0;
}

// -----------------------------------------------------------------------------------------------
