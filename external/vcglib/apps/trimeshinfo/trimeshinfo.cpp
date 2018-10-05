/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
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
Revision 1.41  2007/03/08 22:49:35  cignoni
Changed the include order and the order in which cleaning filters are applied.

Revision 1.40  2006/05/16 21:55:28  cignoni
Updated to the new remove zero area syntax

Revision 1.39  2006/03/29 09:27:07  cignoni
Added managemnt of non critical errors

Revision 1.38  2006/03/29 08:17:56  corsini
fix typo

Revision 1.37  2006/03/27 07:15:59  cignoni
Added LoadMask

Revision 1.36  2006/02/28 14:54:10  corsini
Fix load mask initialization

Revision 1.35  2006/02/09 16:12:27  corsini
Add normal, color, texture information

Revision 1.34  2006/02/06 12:59:12  corsini
Fix mesh info structure initialization

Revision 1.33  2006/01/27 14:18:07  corsini
fix boolean entry type

Revision 1.32  2006/01/27 14:17:10  corsini
remove junk code

Revision 1.31  2006/01/27 13:35:51  corsini
fix signed/unsigned mismatch

Revision 1.30  2006/01/26 16:29:21  corsini
fix typo

Revision 1.29  2005/12/29 12:27:35  cignoni
Splitted IsComplexManifold in IsTwoManifoldFace and IsTwoManifoldVertex

Revision 1.28  2005/12/21 14:11:59  corsini
Out the number of self intersection

Revision 1.27  2005/12/21 13:26:58  corsini
Re-add save xml feature

Revision 1.26  2005/12/21 13:10:10  corsini
Move duplicated vertices routine
Modify genus computation call

Revision 1.25  2005/12/20 14:44:10  corsini
Modify html table

Revision 1.24  2005/12/20 13:29:41  corsini
Remove unuseful Clear function

Revision 1.23  2005/12/19 15:00:53  corsini
Disable xml output temporarily

Revision 1.22  2005/12/19 11:35:13  corsini
Add html output support

Revision 1.21  2005/12/16 11:16:36  corsini
Add manifold check to some properties

Revision 1.20  2005/12/15 11:20:00  corsini
Add vertex-face topology

Revision 1.19  2005/12/14 14:05:37  corsini
Adjust comments

Revision 1.18  2005/12/14 12:15:37  corsini
Re-add clean mesh saving feature

Revision 1.17  2005/12/13 15:46:30  corsini
Restructuring code

Revision 1.16  2005/12/12 12:09:08  cignoni
Changed names of clean function and tested inertia.h

Revision 1.15  2005/12/12 11:29:21  corsini
Minor changes

Revision 1.14  2005/12/12 10:48:16  corsini
Fix indentation

Revision 1.13  2005/11/17 00:42:03  cignoni
Changed include order
removed clean::initialize

Revision 1.12  2005/11/16 16:45:51  rita_borgo
Minor changes

Revision 1.11  2005/11/16 15:59:46  cignoni
Changed name of the component from <Flag> to <BitFlags>

Revision 1.10  2005/11/14 09:21:07  cignoni
Heavily restructured the code of Trimeshinfo.
Now divided the collecting part from the reporting one (xml and ascii)

Revision 1.9  2005/11/04 15:37:57  rita_borgo
Removed Debug option

Revision 1.8  2005/10/11 16:03:40  rita_borgo
Moved all the main functions inside clean.h

Revision 1.7	2005/09/30 15:48:46	 rita_borgo
Fixed	manifold Test

Revision 1.6 2005/09/30	14:13:01 rita_borgo
Problem: Text	not	aligned

Revision 1.5 2005/09/30	13:29:40 rita_borgo
Fixed	Manifold Test

Revision 1.4 2005/09/29	14:48:15 rita_borgo
Fixed	code related to	creation of	the	XML	file

Revision 1.3 2005/09/28	13:57:09 rita_borgo
Fixed	some printout	not	alligned

Revision 1.2 2005/09/28	10:46:04 rita_borgo
Added	possibility	of saving	File in	OFF	format

Revision 1.1 2005/09/20	10:15:27 rita_borgo
Changed	file name	to uniform with	other	solution projects,
before was main.cpp

Revision 1.8 2005/02/15	12:26:06 rita_borgo
Minor	changes	to self-intersection

Revision 1.7 2005/02/07	15:44:31 rita_borgo
Fixed	Color	and	Volume

Revision 1.6 2005/02/01	17:37:53 rita_borgo
Fixed	Volume and Color

Revision 1.5 2005/01/18	16:33:12 rita_borgo
Added	OFF	file Option

Revision 1.4 2005/01/17	18:19:00 rita_borgo
Added	new	routines.
Self-intersection	first	release

Revision 1.2 2005/01/03	16:13:09 rita_borgo
Added	Standard comments

****************************************************************************/

// Standard headers
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stack>

using namespace std;

// VCG headers
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/vertex/component.h>

#include <vcg/simplex/face/base.h>
#include <vcg/simplex/face/component.h>

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/edges.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/space/intersection/triangle_triangle3.h>

#include <vcg/math/histogram.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <vcg/simplex/face/pos.h> 
#include <vcg/complex/algorithms/inertia.h> 

#include "XMLTree.h"

#include <vcg/space/index/grid_static_ptr.h>
#include "defs.h"

using namespace vcg;


class CVertex;
class CFace;

struct MyTypes: public UsedTypes< Use<CVertex>::AsVertexType,Use<CFace>::AsFaceType>{};

class CVertex  : public Vertex< MyTypes, vertex::VFAdj, vertex::Coord3f,vertex::BitFlags, vertex::Normal3f > {};
class CFace    : public Face< MyTypes, face::FFAdj, face::VFAdj, face::VertexRef, face::Normal3f, face::BitFlags, face::Mark > {};
class CMesh    : public vcg::tri::TriMesh< vector<CVertex>, vector<CFace> > {};

typedef CMesh::VertexPointer VertexPointer;
typedef CMesh::VertexIterator VertexIterator;
typedef Point3<CMesh::ScalarType> Point3x;
typedef vector<Point3x> Hole;

typedef CMesh::VertexPointer VertexPointer;
typedef CMesh::VertexIterator VertexIterator;
typedef CMesh::FaceContainer FaceContainer;
typedef CMesh::ScalarType ScalarType;

struct MeshInfo
{
	string FileName;
	bool hasVNormal;
	bool hasFNormal;
	bool hasVColor;
	bool hasFColor;
	bool hasTexture;
	int vn,fn;
	bool VManifold;
	bool FManifold;
	int count_e,boundary_e,count_fd,count_uv,numholes;
	int BEdges;
	float Volume;
	int numcomponents,Genus;
	bool Regular,Semiregular;
	bool Orientable,Oriented;
	int dv;
	bool SelfIntersect;
	std::vector<CMesh::FaceType *> intersections;
};


static const int HTML_LINES = 31;
static const char * HTML_TABLE[HTML_LINES]=
{
"<html>",
"  <head>",
"    <meta http-equiv=\"content-type\" content=\"text/html; charset=ISO-8859-1\">",
"    <title>description.html</title>",
"  </head>",
"  <body>",
"    <span style=\"font-weight: bold;\"></span>",
"    <h2>TriMeshInfo V.1.2 - Results</h2>",
"    <hr style=\"width: 100%; height: 2px;\"><br>",
"    <table",
"      style=\"width: 100%; text-align: center; margin-left: auto; margin-right: auto;\" ",
"      border=\"1\" cellpadding=\"2\" cellspacing=\"2\">",
"      <tbody>",
"        <tr>",
"          <td>Name</td>",
"          <td>Vertices</td>",
"          <td>Faces</td>",
"          <td>Edges</td>",
"          <td>Holes/<br>Boundaries</td>",
"          <td>Connected<br>Components</td>",
"          <td>Isolated<br>Vertices</td>",
"          <td>Duplicated<br>Vertices</td>",
"          <td>Self<br>Interesection</td>",
"          <td>Manifold</td>",
"          <td>Orientable/<br>Oriented</td>",
"          <td>Genus</td>",
"        </tr>",
"      </tbody>",
"    </table>",
"  </body>",
"</html>"
};


int OpenMesh(const	char *filename, CMesh &m)
{
	printf("    Mesh loading...");

	int mask = 0;
	
	tri::io::Importer<CMesh>::LoadMask(filename,mask);
	int err = tri::io::Importer<CMesh>::Open(m, filename);

	if (err)
	{
   if(tri::io::Importer<CMesh>::ErrorCritical(err))
   {
		printf("\n    Error during loading %s: '%s'\n",filename, 
			tri::io::Importer<CMesh>::ErrorMsg(err));
		exit(-1);
   }
   else
   {
    printf("\n    Non Critical Troubles during loading %s: '%s'\n",filename, 
			tri::io::Importer<CMesh>::ErrorMsg(err));
   }
	}
	else
		printf(" done.\n\n");

	return mask;
}

void initMeshInfo(MeshInfo &mi)
{
	mi.vn = 0;
	mi.fn = 0;
	mi.hasVColor = false;
	mi.hasFColor = false;
	mi.hasFNormal = false;
	mi.hasVNormal = false;
	mi.hasTexture = false;
	mi.VManifold = false;
	mi.FManifold = false;
	mi.count_e = 0;
	mi.boundary_e = 0;
	mi.count_fd = 0;
	mi.count_uv = 0;
	mi.numholes = 0;
	mi.BEdges = 0;
	mi.Volume = 0;
	mi.numcomponents = 0;
	mi.Genus = 0;
	mi.Regular = false;
	mi.Semiregular = false;
	mi.Orientable = false;
	mi.Oriented = false;
	mi.dv = 0;
	mi.SelfIntersect = false;
	// mi.FileName has a constructor
	// mi.intersections has a constructor
}

void PrintMeshInfo(MeshInfo &mi)
{
	printf("    *** Mesh information ***\n\n");
	printf("    Mesh: '%s' \n", mi.FileName.c_str());
	printf("    Number of vertices: %d \n", mi.vn);
	printf("    Number of faces: %d \n",	mi.fn);
	printf("    Number of edges: %d \n", mi.count_e);
	printf("    Number of internal edges: %d \n", mi.count_e-mi.boundary_e);
	printf("    Number of boundary edges: %i \n", mi.boundary_e);
	printf("    Number of degenerated faces: %d\n",	mi.count_fd);
	printf("    Number of unreferenced vertices: %d\n",mi.count_uv);
	printf("    Number of duplicated vertices found: %d\n", mi.dv);
	printf("    Number of holes/boundaries: %d \n", mi.numholes);

	if (mi.hasVNormal)
		printf("    Has Per-Vertex Normal: YES\n");
	else
		printf("    Has Per-Vertex Normal: NO\n");

	if (mi.hasFNormal)
		printf("    Has Per-Face Normal: YES\n");
	else
		printf("    Has Per-Face Normal: NO\n");

	if (mi.hasVColor)
		printf("    Has Per-Vertex Color: YES\n");	
	else
		printf("    Has Per-Vertex Color: NO\n");
	
	if (mi.hasFColor)
		printf("    Has Per-Face Color: YES\n");
	else
		printf("    Has Per-Face Color: NO\n");

	if (mi.hasTexture)
		printf("    Has Texture information: YES\n");
	else
		printf("    Has Texture information: NO\n");

	if ((mi.VManifold && mi.FManifold )&&(mi.Oriented)&&(!mi.numholes))
		printf("    Volume: %f \n", mi.Volume);
	else 
		printf("    Volume: UNDEFINED  (a closed oriented manifold is required)\n");

	printf("    Number of connected components: %d\n",	mi.numcomponents);

	// Orientation
	if (!mi.VManifold && mi.FManifold)
	{
		printf("    Orientable Mesh: NO\n");
		printf("    Oriented Mesh: NO\n");
	}
	else
	{
		if (mi.Orientable)
      printf("    Orientable Mesh: YES\n");
		else 
			printf("    Orientable Mesh: NO\n");

		if (mi.Oriented)
			printf("    Oriented Mesh: YES\n");
		else 
			printf("    Oriented Mesh: NO\n");
	}

	// Manifold
	if (mi.VManifold && mi.FManifold ) 
		printf("    Manifold: YES\n");
	else
		printf("    Manifold: NO\n");

	// Genus
	if (mi.VManifold && mi.FManifold)
		printf("    Genus: %d \n", mi.Genus);
	else
		printf("    Genus: N/A \n");

	// Mesh Type
	if (mi.Regular) 
		printf("    Mesh Type: REGULAR\n");
	else if (mi.Semiregular)
		printf("    Mesh Type: SEMIREGULAR\n");
	else
		printf("    Mesh Type: IRREGULAR\n");

	// Further details
	if (mi.SelfIntersect)
		printf("    Self Intersection: %d\n", mi.intersections.size());
	else
		printf("    Self Intersection: NONE.\n");
}

void SaveXMLInfo(MeshInfo &mi)
{
	XMLTree	doc;
	doc.initializeMain();

	char s[256];
	sprintf(s,"%d",mi.vn);
	doc.addNode(s, VALUE_INTEGER, "Number of Vertices");
	sprintf(s,"%d",mi.fn);
	doc.addNode(s, VALUE_INTEGER, "Number of Faces");

	if(mi.VManifold && mi.FManifold)
		doc.addNode("false", VALUE_BOOL,"Manifold");
	else
		doc.addNode("true", VALUE_BOOL,"Manifold");

	sprintf(s,"%d",mi.count_e);
	doc.addNode(s, VALUE_INTEGER,"Number of Edges");
	sprintf(s,"%d",mi.count_fd);
	doc.addNode(s, VALUE_INTEGER,"Number of Degenerated Faces");

	sprintf(s,"%d",mi.count_uv);
	doc.addNode(s, VALUE_INTEGER,"Number of unreferenced vertices");
	sprintf(s,"%d",mi.numholes);
	doc.addNode(s, VALUE_INTEGER, "Number of Holes");
	sprintf(s,"%d",mi.BEdges);
	doc.addNode(s, VALUE_INTEGER, "Number of Border Edges");

	if (mi.hasVColor)
		doc.addNode("true", VALUE_BOOL, "Per-Vertex Color Information");
	else
		doc.addNode("false", VALUE_BOOL, "Per-Vertex Color Information");

	if (mi.hasFColor)
		doc.addNode("true", VALUE_BOOL, "Per-Face Color Information");
	else
		doc.addNode("false", VALUE_BOOL, "Per-Face Color Information");

	if (mi.hasVNormal)
		doc.addNode("true", VALUE_BOOL, "Per-Vertex Normal");
	else
		doc.addNode("false", VALUE_BOOL, "Per-Vertex Normal");

	if (mi.hasFNormal)
		doc.addNode("true", VALUE_BOOL, "Per-Face Normal");
	else
		doc.addNode("false", VALUE_BOOL, "Per-Face Normal");

	if (mi.hasTexture)
		doc.addNode("true", VALUE_BOOL, "Texture Information");
	else
		doc.addNode("false", VALUE_BOOL, "Texture Information");

	sprintf(s,"%f",mi.Volume);
	doc.addNode(s, VALUE_FLOAT,"Volume");
	sprintf(s,"%d",mi.numcomponents);
	doc.addNode(s, VALUE_INTEGER,"Number of Connected Components");
	sprintf(s,"%d",mi.Genus);
	doc.addNode(s, VALUE_INTEGER,"Genus");

	if (mi.Regular)
		doc.addNode("REGULAR", VALUE_STRING,"Type of Mesh");
	else if (mi.Semiregular)
		doc.addNode("SEMIREGULAR", VALUE_STRING,"Type of Mesh");
	else 
		doc.addNode("IRREGULAR", VALUE_STRING,"Type of Mesh");
	
	if (!mi.VManifold && mi.FManifold) 
	{
		doc.addNode("NO", VALUE_STRING,"Orientable Mesh");
		doc.addNode("NO", VALUE_STRING,"Oriented Mesh");
	}
	else
	{
		doc.addNode(mi.Orientable?"Yes":"No", VALUE_STRING,"Orientable Mesh");
		doc.addNode(mi.Oriented?"Yes":"No", VALUE_STRING,"Oriented Mesh");
	}

	sprintf(s,"%d",mi.dv);
	doc.addNode(s, VALUE_INTEGER,"Duplicated Vertices");
	doc.addNode(mi.SelfIntersect?"Yes":"No", VALUE_STRING,"Self Intersection");

	doc.finalizeMain();

	// save xml tree
	string filename = mi.FileName;

	int l = static_cast<int>(filename.size());
	int index = static_cast<int>(filename.find_first_of('.'));
	filename.erase(index, l - index);
	filename.append(".xml");

	doc.setName(filename.c_str());
	doc.printXMLTree();
}

void SaveMeshInfoHtmlTable(fstream &fout, MeshInfo &mi)
{
	fout << "        <tr>" << std::endl;
	fout << "          <td>" << mi.FileName << "</td>" << std::endl;
	fout << "          <td>" << mi.vn << "</td>" << std::endl;
	fout << "          <td>" << mi.fn << "</td>" << std::endl;
	fout << "          <td>" << mi.count_e << "</td>" << std::endl;
	
	if (mi.VManifold && mi.FManifold)
		fout << "          <td>" << mi.numholes << "</td>" << std::endl;
	else
		fout << "          <td>N/A</td>" << std::endl;


	fout << "          <td>" << mi.numcomponents << "</td>" << std::endl;
	fout << "          <td>" << mi.count_uv << "</td>" << std::endl;
	fout << "          <td>" << mi.dv << "</td>" << std::endl;

	if (mi.SelfIntersect)
		fout << "          <td>" << static_cast<unsigned int>(mi.intersections.size()) 
			<< "</td>" << std::endl;
	else
		fout << "          <td>None</td>" << std::endl;

	if(mi.VManifold && mi.FManifold)
		fout << "          <td>Yes</td>" << std::endl;
	else
		fout << "          <td>No</td>" << std::endl;

	if ((mi.Orientable)&&(mi.Oriented))
		fout << "          <td>Yes / Yes</td>" << std::endl;
	else if ((mi.Orientable)&&(!mi.Oriented))
		fout << "          <td>Yes / No</td>" << std::endl;
	else if (!mi.Orientable)
		fout << "          <td>No / No</td>" << std::endl;

	if (mi.VManifold && mi.FManifold)
		fout << "          <td>" << mi.Genus << "</td>" << std::endl;
	else
		fout << "          <td>N/A</td>" << std::endl;

	fout << "        </tr>" << std::endl;
}

void SaveHtmlInfo(MeshInfo &mi)
{
	char buff[1024];
	bool flagInsert = false;
	ifstream fin;
	fstream fout;
	
	// Try to open
	fin.open("result.html");
	long pos;
	if (fin.is_open())
	{
		while (!fin.eof())
		{
			pos = fin.tellg();
			fin.getline(buff, 1024);
			string str(buff);
			if (str == "      </tbody>")
				break;
		}
		flagInsert = true;
	}
	fin.close();

	if (flagInsert)
		fout.open("result.html", ios::in | ios::out);
	else
		fout.open("result.html", ios::out);

	if (!fout.is_open())
	{
		printf("\n  Impossible to write the HTML output file.\n");
	}
	else
	{
		if (flagInsert)
		{
			// Insert the mesh information into an existing table
			fout.seekp(pos, ios::beg);

			SaveMeshInfoHtmlTable(fout, mi);

			for (int i = HTML_LINES - 4; i < HTML_LINES; i++)
				fout << HTML_TABLE[i] << std::endl;
		}
		else
		{
			// Create a new table
			for (int i = 0; i < HTML_LINES - 4; i++)
				fout << HTML_TABLE[i] << std::endl;

			SaveMeshInfoHtmlTable(fout, mi);

			for (int i = HTML_LINES - 4; i < HTML_LINES; i++)
				fout << HTML_TABLE[i] << std::endl;
		}
	}

	fout.close();
}

int main(int argc, char ** argv)
{
	CMesh m;
	bool saveCleanMeshFlag = false;   // Save the clean mesh
	bool verboseFlag = true;          // Verbose mode on/off
	bool XmlFlag= false;              // XML output enabled/disabled
	bool HtmlFlag = false;            // HTML output enabled/disabled

	string SaveName;
	
	MeshInfo mi;
	initMeshInfo(mi);

	printf("\n  -------------------------------\n"
		"     TriMeshInfo V.1.23 \n"
		"     http://vcg.isti.cnr.it\n"
		"     release date: "__DATE__"\n"
		"  -------------------------------\n\n\n");


	// Parsing arguments
	///////////////////////////////////////////////
	
	if (argc <= 1)
	{
		printf(MSG_ERR_N_ARGS);
		exit(-1);
	}

	mi.FileName = argv[1];

	int i = 2;
	while (i < argc)
	{
		if (argv[i][0] == '-')
		{
			switch(argv[i][1])
			{
			case 'q' : 
				// Quiet mode, disable verbose mode
				verboseFlag = false;
			break;

			case 's':
				// Save the clean mesh with the name specified
				saveCleanMeshFlag = true;

				// Check clean mesh name (minimal check)
				if (i+1 >= argc)
				{
					printf("    Invalid output mesh name.\n\n");
					exit(-1);
				}
				else if (argv[i+1][0] != '-')
					SaveName = argv[i+1];
				else
				{
					printf("    Invalid output mesh name.\n\n");
					exit(-1);
				}

				i++;
			break;

			case 'x' :
				// Enable XML output
				XmlFlag = true;
			break;

			case 'h' :
				// Enable HTML output
				HtmlFlag = true;
			break;

			default:
				printf(MSG_ERR_INVALID_OPTION, argv[i]);
				exit(0);
			break;
			}
		}

		i++;
	};

	// Mesh loading
	//////////////////////////////////////////
	
	int load_mask = OpenMesh(mi.FileName.c_str(), m);

	if (load_mask & vcg::tri::io::Mask::IOM_VERTNORMAL)
		mi.hasVNormal = true;
	else
		mi.hasVNormal = false;

	if (load_mask & vcg::tri::io::Mask::IOM_FACENORMAL)
		mi.hasFNormal = true;
	else
		mi.hasFNormal = false;

	if (load_mask & vcg::tri::io::Mask::IOM_VERTCOLOR)
		mi.hasVColor = true;
	else
		mi.hasVColor = false;

	if (load_mask & vcg::tri::io::Mask::IOM_FACECOLOR)
		mi.hasFColor = true;
	else
		mi.hasFColor = false;

	if( (load_mask & vcg::tri::io::Mask::IOM_VERTTEXCOORD) ||
      (load_mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD) )
		mi.hasTexture = true;
	else
		mi.hasTexture = false;

	// Mesh processing
	//////////////////////////////////////////

	printf("    Mesh processing...\n\n");
	
	// Number of vertices
	mi.vn = m.vn;

	// Number of faces
	mi.fn = m.fn;

	// DUPLICATED VERTICES
	mi.dv = tri::Clean<CMesh>::RemoveDuplicateVertex(m);

	// DEGENERATED FACES => (faces with area zero)
	mi.count_fd = tri::Clean<CMesh>::RemoveDegenerateFace(m);
	mi.count_fd += tri::Clean<CMesh>::RemoveFaceOutOfRangeArea<false>(m,0);
	
	// UNREFERENCED VERTEX
	mi.count_uv = tri::Clean<CMesh>::RemoveUnreferencedVertex(m);
  
	// Update topology (face-to-face)
	tri::UpdateTopology<CMesh>::FaceFace(m);
	tri::UpdateTopology<CMesh>::VertexFace(m);

	// IS MANIFOLD?
	mi.VManifold = tri::Clean<CMesh>::CountNonManifoldVertexFF(m)>0;
	mi.FManifold = tri::Clean<CMesh>::CountNonManifoldEdgeFF(m)>0;
 
	// COUNT EDGES
	tri::Clean<CMesh>::CountEdges(m, mi.count_e, mi.boundary_e);

	// HOLES COUNT
	if(mi.VManifold && mi.FManifold)
	{
		mi.numholes = tri::Clean<CMesh>::CountHoles(m);
		tri::Clean<CMesh>::CountEdges(m, mi.BEdges,mi.numholes);
	}

	// CONNECTED COMPONENTS
	mi.numcomponents = tri::Clean<CMesh>::CountConnectedComponents(m);

	// ORIENTATION
	if (mi.VManifold && mi.FManifold)
		tri::Clean<CMesh>::IsOrientedMesh(m, mi.Oriented, mi.Orientable);
	else
	{
		mi.Oriented = false;
		mi.Orientable = false;
	}

	// Rebuild Vertex-Face topology
	tri::UpdateTopology<CMesh>::VertexFace(m);

	// VOLUME (require a closed oriented manifold)
	if ((mi.VManifold && mi.FManifold)&&(mi.Oriented)&&(!mi.numholes))
	{
		tri::Inertia<CMesh> mm;
		mm.Compute(m);
		mi.Volume = mm.Mass();

		// the sign of the volume depend on the mesh orientation
		if (mi.Volume < 0.0)
			mi.Volume = -mi.Volume;
	}

	// GENUS
	if(mi.VManifold && mi.FManifold) 
		mi.Genus = tri::Clean<CMesh>::MeshGenus(m, mi.numholes,
			mi.numcomponents, mi.count_e);
	
	// REGULARITY
	if (mi.VManifold && mi.FManifold)
		tri::Clean<CMesh>::IsRegularMesh(m, mi.Regular, mi.Semiregular);
	else
	{
		mi.Regular = false;
		mi.Semiregular = false;
	}

	// SELF INTERSECTION
	mi.SelfIntersect = tri::Clean<CMesh>::SelfIntersections(m, mi.intersections);

	// Mesh Information Output
	//////////////////////////////////////////

	// Print mesh information
	if(verboseFlag)
		PrintMeshInfo(mi);

	// Save mesh information in XML format
	if(XmlFlag)
		SaveXMLInfo(mi);

	// Save mesh information in HTML format
	if (HtmlFlag)
		SaveHtmlInfo(mi);

	// Save the clean mesh
	if (saveCleanMeshFlag)
	{
		printf("    Save the 'clean' mesh...");
		tri::io::Exporter<CMesh>::Save(m, SaveName.c_str());
		printf(" done.\n\n");
	}

	mi.intersections.clear();

	return 0;
}

