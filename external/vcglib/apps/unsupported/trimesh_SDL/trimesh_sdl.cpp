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
Revision 1.4  2007/03/17 13:22:55  ganovelli
some more comment added

Revision 1.3  2006/11/21 17:44:06  cignoni
Update to the latest version of the trackball

Revision 1.2  2005/11/22 17:50:15  cignoni
Refactored the sample code.
Shortened a lot and removed all unused unnecessary stuff

Revision 1.1  2005/09/21 10:29:33  cignoni
Initial Relase

****************************************************************************/
#include <SDL.h>
#include <GL/glew.h>
#include <vector>

/*include the base definition for the vertex */
#include <vcg/simplex/vertex/base.h>

/*include the base definition for the face */
#include <vcg/simplex/face/base.h>

/*include the base definition for the trimesh*/
#include <vcg/complex/complex.h>

/*include the algorihm that update bounding box and normals*/
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>

/*include the importer from disk*/
#include <wrap/io_trimesh/import.h>

/*include wrapping of the trimesh towards opengl*/
#include <wrap/gl/trimesh.h>

/*include the trackball: NOTE, you the implementation of the trackball is found in the files:
wrap/gui/trackball.cpp and wrap/gui/trackmode.cpp. You should include these files in your solution
otherwise you'll get linking errors */
#include <wrap/gui/trackball.h>

using namespace vcg;
using namespace std;


class CFace;
class CVertex;
struct MyUsedTypes : public UsedTypes<	Use<CVertex>		::AsVertexType,
																				Use<CFace>			::AsFaceType>{};

/* define a vertex passing the attributes you want it to have. Each attributes has its own class.
Check vcg/simplex/vertex/component.h to find out the existing attributes. Note: then you could 
also personalized attributes */
class CVertex : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f >{};

/*same as for the vertes */ 
class CFace   : public Face<   MyUsedTypes, face::VertexRef, face::Normal3f > {};

/*the mesh is a container of vertices and a container of faces */ 
class CMesh   : public vcg::tri::TriMesh< vector<CVertex>, vector<CFace> > {};


////////////////////////////////////////////////////////////////////////////
// Globals: the mesh, the OpenGL wrapper to draw the mesh and the trackball.
CMesh mesh;
vcg::GlTrimesh<CMesh> glWrap;
vcg::Trackball track;
int drawMode;
int width =1024;
int height = 768;

////////////////////////////////////////////////////////////////////////////


void initGL()
{
  glClearColor(0, 0, 0, 0); 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
}


void myReshapeFunc(GLsizei w, GLsizei h)
{
  SDL_SetVideoMode(w, h, 0, SDL_OPENGL|SDL_RESIZABLE|SDL_DOUBLEBUF);
	width=w;
  height=h;
  glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
  initGL();
}

bool initSDL(const std::string &str) {
  /* Initialize SDL for video output */
  if ( SDL_Init(SDL_INIT_VIDEO) < 0 ) {
    fprintf(stderr, "Unable to initialize SDL: %s\n", SDL_GetError());
    exit(1);
  }
  if ( SDL_SetVideoMode(width, height, 0, SDL_OPENGL|SDL_RESIZABLE) == NULL ) {
    fprintf(stderr, "Unable to create OpenGL screen: %s\n", SDL_GetError());
    SDL_Quit();
    exit(2);
  }
  
  SDL_WM_SetCaption(str.c_str(), str.c_str());  
  myReshapeFunc(width, height);
  return true;
}

void display(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, width/(float)height, 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,5,   0,0,0,   0,1,0);    

    track.center=Point3f(0, 0, 0);
    track.radius= 1;

		track.GetView();
	track.Apply();
    glPushMatrix();
    float d=1.0f/mesh.bbox.Diag();
    glScale(d);
		glTranslate(-glWrap.m->bbox.Center());	

		// the trimesh drawing calls
		switch(drawMode)
		{
		  case 0: glWrap.Draw<vcg::GLW::DMSmooth,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  case 1: glWrap.Draw<vcg::GLW::DMPoints,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  case 2: glWrap.Draw<vcg::GLW::DMWire,     vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  case 3: glWrap.Draw<vcg::GLW::DMFlatWire, vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  case 4: glWrap.Draw<vcg::GLW::DMHidden,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  case 5: glWrap.Draw<vcg::GLW::DMFlat,     vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  default: break;
		}
    glPopMatrix();
    track.DrawPostApply();
    SDL_GL_SwapBuffers();
}

// The Event Loop Processor
int sdl_idle() {
  bool quit=false;
 	SDL_Event event;
	while( !quit ) {  
    SDL_WaitEvent(&event);
    switch( event.type ) {
      case SDL_QUIT:  quit = true; break; 
      case SDL_VIDEORESIZE : 			myReshapeFunc(event.resize.w,event.resize.h); 			break;
      case SDL_KEYDOWN:                                        
			  switch(event.key.keysym.sym) {
			    case SDLK_RCTRL:
			    case SDLK_LCTRL: track.ButtonDown(vcg::Trackball::KEY_CTRL); break;
			    case SDLK_q: exit(0); break;	
			    case SDLK_SPACE: drawMode=((drawMode+1)%6); printf("Current Mode %i\n",drawMode); break;	
			  }  break;
      case SDL_KEYUP: 
			  switch(event.key.keysym.sym) {
			    case SDLK_RCTRL:
			    case SDLK_LCTRL: track.ButtonUp(vcg::Trackball::KEY_CTRL); break;
			  }	break;
      case SDL_MOUSEBUTTONDOWN:   
	      switch(event.button.button) {
          case SDL_BUTTON_WHEELUP:    track.MouseWheel( 1); break;
          case SDL_BUTTON_WHEELDOWN:  track.MouseWheel(-1); break;
          case SDL_BUTTON_LEFT:	      track.MouseDown(event.button.x, (height - event.button.y), vcg::Trackball::BUTTON_LEFT); break;
          case SDL_BUTTON_RIGHT:	    track.MouseDown(event.button.x, (height - event.button.y), vcg::Trackball::BUTTON_RIGHT);break;
        } break;
      case SDL_MOUSEBUTTONUP:          
	      switch(event.button.button) {
          case SDL_BUTTON_LEFT:	      track.MouseUp(event.button.x, (height - event.button.y), vcg::Trackball::BUTTON_LEFT); break;
          case SDL_BUTTON_RIGHT:	    track.MouseUp(event.button.x, (height - event.button.y), vcg::Trackball::BUTTON_RIGHT);break;
        } break;
      case SDL_MOUSEMOTION: 
	      while(SDL_PeepEvents(&event, 1, SDL_GETEVENT, SDL_MOUSEMOTIONMASK));
	      track.MouseMove(event.button.x, (height - event.button.y));
	      break;  
      case SDL_VIDEOEXPOSE:
      default: break;
      }
		display();
	}

  SDL_Quit();
  return -1;
}



int main(int argc, char *argv[]) {	
	// Generic loading of the mesh from disk
	if(vcg::tri::io::Importer<CMesh>::Open(mesh,argv[1])!=0) {
      fprintf(stderr,"Error reading file %s\n",argv[1]);
			exit(0);
		}

  // Initialize the mesh itself
	vcg::tri::UpdateBounding<CMesh>::Box(mesh);      // update bounding box
  //vcg::tri::UpdateNormals<CMesh>::PerVertexNormalizePerFaceNormalized(mesh); // update Normals
  vcg::tri::UpdateNormals<CMesh>::PerVertexPerFace(mesh); // update Normals
	
	// Initialize the wrapper
  glWrap.m = &mesh;
  glWrap.Update();
	
  initSDL("SDL_minimal_viewer");
  initGL();
  sdl_idle();
	exit(0);
}



