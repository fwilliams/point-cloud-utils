#include <stdio.h>
#include <vector>
#include<vcg/math/base.h>
#include<vcg/space/point3.h>
#include<vcg/space/point4.h>
#include<vcg/space/color4.h>
#include<vcg/math/matrix44.h>

#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/vertex/component.h>
#include <vcg/simplex/face/base.h>
#include <vcg/simplex/face/component.h>
#include <vcg/simplex/face/component_rt.h>

#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/update/flag.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/color.h>
#include<vcg/complex/algorithms/clean.h>

#include<wrap/io_trimesh/import_ply.h>
#include<wrap/io_trimesh/export_ply.h>

using namespace vcg;

class MyEdge;
class MyFaceC;
class MyFace;
class MyVertexC   : public VertexSimp2<MyVertexC,MyEdge,MyFaceC,vert::Coord3f,vert::Color4b,vert::Qualityf,vert::Normal3f,vert::BitFlags> {};
class MyFaceC     : public FaceSimp2< MyVertexC,MyEdge,MyFaceC,face::VertexRef, face::Normal3f,face::BitFlags> {};
class MyMeshC     : public tri::TriMesh< std::vector<MyVertexC>, std::vector<MyFaceC> > {};

class MyVertex   : public VertexSimp2<MyVertex,MyEdge,MyFace,vert::Coord3f,vert::Normal3f,vert::BitFlags> {};
class MyFace     : public FaceSimp2<  MyVertex,MyEdge,MyFace,face::VertexRef, face::Normal3f,face::BitFlags> {};
class MyMesh     : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

using namespace std;
using namespace tri;

/*
  class MyEdge;
  class MyFace;
  class MyEdgeC;
  class MyFaceC;

  class MyVertexC:public VertexVCVN<float,MyEdgeC,MyFace>{};
  class MyFaceC :public FaceFN<MyVertexC,MyEdgeC,MyFaceC>{};
  class MyMeshC: public tri::TriMesh< std::vector<MyVertexC>, std::vector<MyFaceC > >{};

  class MyVertex:public VertexVN<float,MyEdge,MyFace>{};
  class MyFace :public FaceFN<MyVertex,MyEdge,MyFace>{};
  class MyMesh: public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};
*/

//-------------------------------------------------------

  int			nummeshes;			// number of meshes extracted so far
  MyMesh		currentmesh;		// current mesh, read from stream and saved one completed
  MyMeshC		currentmeshC;		// current mesh, read from stream and saved one completed
  Matrix44f		currtrasf;
  float			angle;				// angle treshold for face deletion
  int			singlemap;			// single map mode, which map is to be saved. if -1 then all map are saved
  int			frommap;			// skip all maps BEFORE this index
  int			tomap;				// skip all maps AFTER this index
  bool			savecolor;			// if has color, save it on 3dmesh
  bool			hascolor;			// true if the current mesh has color
  bool          saveall;			// all elements are keeped (even invalids)
  bool          flipfaces;			// flip all faces
  int			todump;
  bool			dumpit;
  bool          unpack;
  bool			onlypoints;			// store only points
  bool			switchside;			// inverse triangulation order (swaping row->cols)




// read the current mesh from the stream
int readmesh(FILE* fp)
{
	int colnum;
	int rownum;
	int trinum;
	int numtokens;
	int rit,cit;
	char linebuf[256];
	int ii;
	float xx,yy,zz;	// position
	float rr,gg,bb;	// color
	float rf;		// reflectance
	MyMesh::FaceIterator fi;
	MyMesh::VertexIterator vi;
	MyMeshC::FaceIterator fiC;
	MyMeshC::VertexIterator viC;

	// cleaning mesh
	currentmesh.Clear();
	currentmeshC.Clear();

	// getting mesh size;
	fscanf(fp,"%i\n",&colnum);
	fscanf(fp,"%i\n",&rownum);

	// initial 4 lines [still don't know what is this :) :)]
	fscanf(fp,"%f %f %f\n", &xx, &yy, &zz);
	fscanf(fp,"%f %f %f\n", &xx, &yy, &zz);
	fscanf(fp,"%f %f %f\n", &xx, &yy, &zz);
	fscanf(fp,"%f %f %f\n", &xx, &yy, &zz);

	// now the transformation matrix
	fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(0,0)), &(currtrasf.ElementAt(0,1)), &(currtrasf.ElementAt(0,2)), &(currtrasf.ElementAt(0,3)));
	fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(1,0)), &(currtrasf.ElementAt(1,1)), &(currtrasf.ElementAt(1,2)), &(currtrasf.ElementAt(1,3)));
	fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(2,0)), &(currtrasf.ElementAt(2,1)), &(currtrasf.ElementAt(2,2)), &(currtrasf.ElementAt(2,3)));
	fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(3,0)), &(currtrasf.ElementAt(3,1)), &(currtrasf.ElementAt(3,2)), &(currtrasf.ElementAt(3,3)));

	// now the real data begins

	// first line, we should know if the format is
	// XX YY ZZ RF
	// or it is
	// XX YY ZZ RF RR GG BB

	// read the entire first line and then count the spaces. it's rude but it works :)
	ii=0;
	fread(&(linebuf[ii++]),1,1,fp);
	while(linebuf[ii-1] != '\n')
       fread(&(linebuf[ii++]),1,1,fp);
    linebuf[ii-1] = '\0'; // terminate the string

    numtokens=1;
	for(ii=0; ii<strlen(linebuf); ii++)
	{
		if(linebuf[ii] == ' ')
			numtokens++;
	}

	if(numtokens == 4)
	  hascolor = false;
	else if(numtokens == 7)
	  hascolor = true;
	else
	  return -1;

	// allocate all points
	printf("\n %i x %i \n", rownum, colnum);
	printf(" expect V %i F %i",(rownum*colnum),((rownum-1)*(colnum-1)*2));

	if(hascolor && savecolor)
	{
	 viC = Allocator<MyMeshC>::AddVertices(currentmeshC,(rownum*colnum));
	}
	else
	{
	 vi = Allocator<MyMesh>::AddVertices(currentmesh,(rownum*colnum));
	}

	// parse the first line....
	if(hascolor)
	{
	  printf("\n hascolor ");
	  sscanf(linebuf,"%f %f %f %f %f %f %f", &xx, &yy, &zz, &rf, &rr, &gg, &bb);
	}
	else
	{
	  printf("\n no color ");
	  sscanf(linebuf,"%f %f %f %f", &xx, &yy, &zz, &rf);
	}

	//addthefirstpoint
	if(hascolor && savecolor)
	{
	 (*viC).P()[0]=xx;
	 (*viC).P()[1]=yy;
	 (*viC).P()[2]=zz;
	 (*viC).Q()=rf;
	 (*viC).C()[0]=rr;
	 (*viC).C()[1]=gg;
	 (*viC).C()[2]=bb;
	 viC++;
	}
	else
	{
	 (*vi).P()[0]=xx;
	 (*vi).P()[1]=yy;
	 (*vi).P()[2]=zz;
	 vi++;
	}

	// now for each line until end of mesh (row*col)-1
	for(ii=0; ii<((rownum*colnum)-1); ii++)
	{
	 // read the stream
	 if(hascolor)
	   fscanf(fp,"%f %f %f %f %f %f %f", &xx, &yy, &zz, &rf, &rr, &gg, &bb);
	 else
	   fscanf(fp,"%f %f %f %f", &xx, &yy, &zz, &rf);

	 // add the point
	 if(hascolor && savecolor)
	 {
	  (*viC).P()[0]=xx;
	  (*viC).P()[1]=yy;
	  (*viC).P()[2]=zz;
	  (*viC).Q()=rf;
		(*viC).C()[0]=rr;
	  (*viC).C()[1]=gg;
	  (*viC).C()[2]=bb;
	  viC++;
 	 }
	 else
	 {
	  (*vi).P()[0]=xx;
	  (*vi).P()[1]=yy;
	  (*vi).P()[2]=zz;
	  vi++;
	 }

	}

	currentmesh.vn = currentmesh.vert.size();


	if(! onlypoints)
	{

		// now i can triangulate
		trinum = (rownum-1) * (colnum-1) * 2;

		if(hascolor && savecolor)
		{
		fiC= Allocator<MyMeshC>::AddFaces(currentmeshC,trinum);
		}
		else
		{
		fi= Allocator<MyMesh>::AddFaces(currentmesh,trinum);
		}


		currentmesh.fn = 0;
		currentmeshC.fn = 0;
		int v0i,v1i,v2i;
		for(rit=0; rit<rownum-1; rit++)
		for(cit=0; cit<colnum-1; cit++)
		if(hascolor && savecolor)
		{

			if(!switchside)
			{
				v0i = (rit  ) + ((cit  ) * rownum);
				v1i = (rit+1) + ((cit  ) * rownum);
				v2i = (rit  ) + ((cit+1) * rownum);
			}
			else
			{
				v0i = (cit  ) + ((rit  ) * colnum);
				v1i = (cit+1) + ((rit  ) * colnum);
				v2i = (cit  ) + ((rit+1) * colnum);
			}


			// upper tri
			(*fiC).V(2) = &(currentmeshC.vert[v0i]);
			(*fiC).V(1) = &(currentmeshC.vert[v1i]);
			(*fiC).V(0) = &(currentmeshC.vert[v2i]);

			if(flipfaces)
			{
			(*fiC).V(2) = &(currentmeshC.vert[v1i]);
			(*fiC).V(1) = &(currentmeshC.vert[v0i]);
			}

			currentmeshC.fn++;
			fiC++;

			if(!switchside)
			{
				v0i = (rit+1) + ((cit  ) * rownum);
				v1i = (rit+1) + ((cit+1) * rownum);
				v2i = (rit  ) + ((cit+1) * rownum);
			}
			else
			{
				v0i = (cit+1) + ((rit  ) * colnum);
				v1i = (cit+1) + ((rit+1) * colnum);
				v2i = (cit  ) + ((rit+1) * colnum);
			}

			// lower tri
			(*fiC).V(2) = &(currentmeshC.vert[v0i]);
			(*fiC).V(1) = &(currentmeshC.vert[v1i]);
			(*fiC).V(0) = &(currentmeshC.vert[v2i]);

			if(flipfaces)
			{
			(*fiC).V(2) = &(currentmeshC.vert[v1i]);
			(*fiC).V(1) = &(currentmeshC.vert[v0i]);
			}

			currentmeshC.fn++;
			fiC++;
		}
		else
		{
			// upper tri
			if(!switchside)
			{
				v0i = (rit  ) + ((cit  ) * rownum);
				v1i = (rit+1) + ((cit  ) * rownum);
				v2i = (rit  ) + ((cit+1) * rownum);
			}
			else
			{
				v0i = (cit  ) + ((rit  ) * colnum);
				v1i = (cit+1) + ((rit  ) * colnum);
				v2i = (cit  ) + ((rit+1) * colnum);
			}

			(*fi).V(2) = &(currentmesh.vert[v0i]);
			(*fi).V(1) = &(currentmesh.vert[v1i]);
			(*fi).V(0) = &(currentmesh.vert[v2i]);

			if(flipfaces)
			{
			(*fi).V(2) = &(currentmesh.vert[v1i]);
			(*fi).V(1) = &(currentmesh.vert[v0i]);
			}

			currentmesh.fn++;
			fi++;

			// lower tri
			if(!switchside)
			{
				v0i = (rit+1) + ((cit  ) * rownum);
				v1i = (rit+1) + ((cit+1) * rownum);
				v2i = (rit  ) + ((cit+1) * rownum);
			}
			else
			{
				v0i = (cit+1) + ((rit  ) * colnum);
				v1i = (cit+1) + ((rit+1) * colnum);
				v2i = (cit  ) + ((rit+1) * colnum);
			}

			(*fi).V(2) = &(currentmesh.vert[v0i]);
			(*fi).V(1) = &(currentmesh.vert[v1i]);
			(*fi).V(0) = &(currentmesh.vert[v2i]);

			if(flipfaces)
			{
			(*fi).V(2) = &(currentmesh.vert[v1i]);
			(*fi).V(1) = &(currentmesh.vert[v0i]);
			}

			currentmesh.fn++;
			fi++;
		}


		if(hascolor && savecolor)
		printf("\nV: %8i F: %8i \n", currentmeshC.vn, currentmeshC.fn);
		else
		printf("\nV: %8i F: %8i \n", currentmesh.vn, currentmesh.fn);
	}


	if(! saveall)
	{
	printf("deleting unsampled points \n");
	// now i delete all points in (0,0,0) that are unsampled points
	if(hascolor && savecolor)
	for(viC = currentmeshC.vert.begin(); viC != currentmeshC.vert.end(); viC++)
 	{
		if((*viC).P() == Point3f(0.0, 0.0, 0.0))
		{
			(*viC).SetD();
			currentmeshC.vn--;
		}
	}
	else
	for(vi = currentmesh.vert.begin(); vi != currentmesh.vert.end(); vi++)
	{
		if((*vi).P() == Point3f(0.0, 0.0, 0.0))
		{
			(*vi).SetD();
			currentmesh.vn--;
		}
	}
	}


	if(! onlypoints)
	{

		if(! saveall)
		{
		printf("deleting invalid faces \n");
		// and then i delete all faces with null vertices
		if(hascolor && savecolor)
		for(fiC = currentmeshC.face.begin(); fiC != currentmeshC.face.end(); fiC++)
		{
			if( ((*fiC).V(0)->IsD()) || ((*fiC).V(1)->IsD()) || ((*fiC).V(2)->IsD()) )
			{
				(*fiC).SetD();
				currentmeshC.fn--;
			}
		}
		else
		for(fi = currentmesh.face.begin(); fi != currentmesh.face.end(); fi++)
		{
			if( ((*fi).V(0)->IsD()) || ((*fi).V(1)->IsD()) || ((*fi).V(2)->IsD()) )
			{
				(*fi).SetD();
				currentmesh.fn--;
			}
		}
		}

		if(hascolor && savecolor)
		printf("V: %8i F: %8i \n", currentmeshC.vn, currentmeshC.fn);
		else
		printf("V: %8i F: %8i \n", currentmesh.vn, currentmesh.fn);


		// eliminate high angle triangles
		if((angle != 90)&&(!saveall))
		{
		printf(" culling by angle \n");
		float limit = cos( angle*3.14159265358979323846/180.0 );
		Point3f raggio;

		if(hascolor && savecolor)
		{
 		tri::UpdateNormals<MyMeshC>::PerFaceNormalized(currentmeshC);
		for(fiC = currentmeshC.face.begin(); fiC != currentmeshC.face.end(); fiC++)
		if(!(*fiC).IsD())
			{
				raggio = -((*fiC).V(0)->P() + (*fiC).V(1)->P() + (*fiC).V(2)->P()) / 3.0;
				raggio.Normalize();
				if(((*fiC).N() * raggio) < limit)
				{
				(*fiC).SetD();
				currentmeshC.fn--;
				}
			}
		}
		else
		{
 		vcg::tri::UpdateNormals<MyMesh>::PerFaceNormalized(currentmesh);
		for(fi = currentmesh.face.begin(); fi != currentmesh.face.end(); fi++)
		if(!(*fi).IsD())
			{
				raggio = -((*fi).V(0)->P() + (*fi).V(1)->P() + (*fi).V(2)->P()) / 3.0;
				raggio.Normalize();
				if((raggio * (*fi).N()) < limit)
				{
				(*fi).SetD();
				currentmesh.fn--;
				}
			}
		}

		}
	}

	currtrasf.transposeInPlace();

	// apply tranformation
	if(hascolor && savecolor)
	{
 	 for(viC = currentmeshC.vert.begin(); viC != currentmeshC.vert.end(); viC++)
	  if(!(*viC).IsD())
		{
		  (*viC).P() = currtrasf * (*viC).P();
		}
	}
	else
	{
 	 for(vi = currentmesh.vert.begin(); vi != currentmesh.vert.end(); vi++)
	  if(!(*vi).IsD())
		{
		  (*vi).P() = currtrasf * (*vi).P();
		}
	}

	if(hascolor && savecolor)
	{
     int dup = tri::Clean<MyMeshC>::RemoveDuplicateVertex(currentmeshC);
	 if(! onlypoints)
       int unref =  tri::Clean<MyMeshC>::RemoveUnreferencedVertex(currentmeshC);
	}
	else
	{
     int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(currentmesh);
	 if(! onlypoints)
		int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(currentmesh);
	}

	return 0;
}

// save each mesh in a separate file
void dounpack(FILE* fp)
{
 FILE* outf;
 char namef[128];
 int rnum;
 char linebuf[256];
 int ii;
 bool trovato;

 rnum=1;

	 trovato = false;

	 // search for the first integer
	 while(!trovato)
	 {
	  // read the entire first line and then count the spaces. it's rude but it works :)
	  ii=0;
	  fread(&(linebuf[ii++]),1,1,fp);
	  while(linebuf[ii-1] != '\n')
        fread(&(linebuf[ii++]),1,1,fp);
      linebuf[ii-1] = '\0'; // terminate the string

	  //check the string
	  if(strchr(linebuf,' ') == NULL)
		  trovato = true;
	 }

 while(!feof(fp))
 {

     sprintf(namef,"range%03i.ptx",rnum++);
	 outf = fopen(namef,"w");

	 // write first integer
	 fprintf(outf,"%s\n",linebuf);

	 // read and write next int
	 ii=0;
	 fread(&(linebuf[ii++]),1,1,fp);
	 while(linebuf[ii-1] != '\n')
       fread(&(linebuf[ii++]),1,1,fp);
     linebuf[ii-1] = '\0'; // terminate the string
	 fprintf(outf,"%s\n",linebuf);

	 // search for the next integer
	 while(!trovato)
	 {
	  // read the entire first line and then count the spaces. it's rude but it works :)
	  ii=0;
	  fread(&(linebuf[ii++]),1,1,fp);
	  while(linebuf[ii-1] != '\n')
        fread(&(linebuf[ii++]),1,1,fp);
      linebuf[ii-1] = '\0'; // terminate the string

	  //if not an integer then write it, otherwise close and remember for next step
	  if(strchr(linebuf,' ') == NULL)
	      trovato = true;
	  else
		  fprintf(outf,"%s\n",linebuf);
	 }

	 fclose(outf);
 }

}

// skip a mesh
int skipmesh(FILE* fp)
{
	int colnum;
	int rownum;
	int skiplines;
	char linebuf;

	if(feof(fp))
		return -1;

	// getting mesh size;
	fscanf(fp,"%i\n",&colnum);
	fscanf(fp,"%i\n",&rownum);

	printf("\n %i x %i \n", rownum, colnum);
	printf(" expect V %i F %i\n",(rownum*colnum),((rownum-1)*(colnum-1)*2));

	if(feof(fp))
		return -1;

	skiplines = (colnum * rownum) + 8; // have to skip (col * row) lines plus 8 lines for the header
	for(int ii=0; ii<skiplines; ii++)
	{
	 fread(&linebuf,1,1,fp);
	 while(linebuf != '\n')
       fread(&linebuf,1,1,fp);
	}

	return 0;
}


void parseparams(int argn, char** argvect)
{
 int pit;

 for(pit = 2; pit<argn; pit++)
 {
  if(argvect[pit][0] != '-')	// invalid param
  {
	  printf("invalid parameter\n");
  }
  else						// valid param
  {
	  if(argvect[pit][1] == 'a')	// angle
	  {
		angle = atof(&(argvect[pit][2]));
		printf("cutoff angle = %f \n",angle);
	  }
	  if(argvect[pit][1] == 'm')	// single map
	  {
		singlemap = atoi(&(argvect[pit][2]));
		frommap = 0;
		tomap = singlemap;
		printf("single map # %i \n",singlemap);
	  }
	  if(argvect[pit][1] == 'f')	// from map
	  {
		frommap = atoi(&(argvect[pit][2]));
		singlemap = -1;
		printf("start from map # %i \n",frommap);
	  }
	  if(argvect[pit][1] == 't')	// from map
	  {
		tomap = atoi(&(argvect[pit][2]));
		singlemap = -1;
		printf("end with map # %i \n",tomap);
	  }
	  if(argvect[pit][1] == 'd')	// dump to file
	  {
		todump = atoi(&(argvect[pit][2]));
		dumpit = true;
		printf("dumping # %i chars\n",todump);
	  }
	  if(argvect[pit][1] == 'u')	// unpack the file in different
	  {
		unpack = true;
		printf("UNPACKING \n");
	  }
	  if(argvect[pit][1] == 'c')	// save color if present
	  {
		savecolor = true;
		printf("saving color \n");
	  }
	  if(argvect[pit][1] == 'k')	// keep all
	  {
		saveall = true;
		printf("keeping all elements \n");
	  }
	  if(argvect[pit][1] == 'f')	// flip all tris
	  {
		flipfaces = true;
		printf("keeping all elements \n");
	  }
	  if(argvect[pit][1] == 'p')	// points only, do not triangulate
	  {
		onlypoints = true;
		printf("points only, do not triangulate \n");
	  }
	  if(argvect[pit][1] == 'r')	// points only, do not triangulate
	  {
		switchside = true;
		printf("swapped triangulation \n");
	  }


  }
 }

}

void printhelp()
{
  printf("-------------------------------------------------------------\n");
  printf("multiple PLY files will be extracted from the PTX \n");
  printf("\n");
  printf("USAGE:    ptx2ply filename.ptx [options]");
  printf("\n");
  printf("each map contained in the file will be saved in a PLY \n");
  printf("\n");
  printf("-aAA angle threshold for face elimination, glazing faces with\n");
  printf("     angle between normal and viewdirection > AA are removed \n");
  printf("     default is no cut.                                      \n");
  printf("     beware! this only works for range maps that still need  \n");
  printf("     to be tranformed in the final reference system          \n");
  printf("     (in this case the viewpoint is the origin)              \n");
  printf("\n");
  printf("-c   save color if present \n");
  printf("\n");
  printf("-mNN extract just the map NN, skip all the rest of the file  \n");
  printf("\n");
  printf("-fNN extract maps starting FROM index NN  \n");
  printf("\n");
  printf("-tNN extract maps UP TO index NN  \n");
  printf("\n");
  printf("-u   unpack the file generating a ptx for each map  \n");
  printf("\n");
  printf("-k   keep all elements (no points/tris discarded)  \n");
  printf("     only useful for debug  \n");
  printf("\n");
  printf("-f   flip all faces  \n");
  printf("\n");
  printf("-p   store points only  \n");
  printf("\n");
  printf("-r   during triangulation, swap rows->columns  \n");
  printf("\n");
  printf("MESH INDICES STARTS FROM 1                                   \n");
  printf("parameters -f and -t can be used together to specify an index\n");
  printf("range to be processed. \n");
  printf("-------------------------------------------------------------\n");
  exit(0);
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  FILE *fp;
  char filename[256];
  char modelname[128];


  //---------------------------------------------------------------
  printf("PTX to PLY conversion\n");
  //---------------------------------------------------------------

  //-- init params
  nummeshes = 1;
  currentmesh.Clear();
  angle = 90.0;
  singlemap = -1;
  frommap = 0;
  tomap = 99999;
  todump = 1024;
  dumpit = false;
  unpack = false;
  savecolor = false;
  saveall = false;
  flipfaces = false;
  onlypoints = false;
  switchside = false;


  if(argc < 2)
   printhelp();


  //--
  parseparams(argc, argv);

  strcpy(modelname,argv[1]);
  modelname[strlen(argv[1])-4] = '\0';

  fp = fopen(argv[1],"r");

  if(unpack)
	  dounpack(fp);

  while((!feof(fp)) && (nummeshes <= tomap))
  {
   printf("mesh %3i ",nummeshes);

   if((nummeshes >= frommap) && (nummeshes <= tomap) && ((singlemap == -1) || (singlemap == nummeshes)))
   {
     if(dumpit)
	 {
	  FILE* outf;
	  char cbuf;

	  outf = fopen("dump.txt","w");

	  for(int dit=0; dit<todump; dit++)
	  {
		  fread(&cbuf,1,1,fp);
		  fwrite(&cbuf,1,1,outf);
	  }

	  fclose(outf);
	  fclose(fp);
	  exit(0);
	 }

	 printf("reading ");
     readmesh(fp);

     sprintf(filename,"%s_%03i.ply",modelname,nummeshes);
     if(!feof(fp))
	 {
	  if(hascolor && savecolor)
	  {
	   int plyMask=tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY;
       tri::io::ExporterPLY<MyMeshC>::Save(currentmeshC,filename, plyMask);
	  }
	  else
       tri::io::ExporterPLY<MyMesh>::Save(currentmesh,filename);
	 }

     printf("ok! \n");

    }
	else
	{
     printf("skipping ");
	 skipmesh(fp);
     printf("ok! \n");
	}

   nummeshes++;
  }

  fclose(fp);

  return 0;
}

