/****************************************************************************
 * VCGLib                                                            o o     *
 * Visual and Computer Graphics Library                            o     o   *
 *                                                                _   O  _   *
 * Copyright(C) 2007                                                \/)\/    *
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

/**
 * Minimal   trimesh viewer made with AntTweakBar and freglut
 *
 * This sample shows how to use togheter: 
 * - the trimesh loading and initialization
 * - basic usage of the default manipulators (the "Trackball")
 */


#include <AntTweakBar.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <stdio.h>

/// vcg imports
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/create/platonic.h>

/// wrapper imports
#include <wrap/io_trimesh/import.h>
#include <wrap/gl/trimesh.h>
#include <wrap/gui/trackball.h>

/// declaring edge and face type

using namespace vcg;
class CFace;
class CVertex;
struct MyUsedTypes : public UsedTypes<	Use<CVertex>		::AsVertexType,
																				Use<CFace>			::AsFaceType>{};

/// compositing wanted proprieties
class CVertex : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags>{};
class CFace   : public vcg::Face<  MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::BitFlags > {};
class CMesh   : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace> > {};

/// the active mesh instance
CMesh mesh;

/// filename of the mesh to load
char * filename  = NULL;

/// the active mesh opengl wrapper
vcg::GlTrimesh<CMesh> glWrap;
/// the active manipulator
vcg::Trackball track;

/// window size
int width,height;

/// we choosed a subset of the avaible drawing modes
enum DrawMode{SMOOTH=0,PERPOINTS,WIRE,FLATWIRE,HIDDEN,FLAT};

/// the current drawmode
DrawMode drawmode;

/// Takes a GLUT MouseButton and returns the equivalent Trackball::Button
static vcg::Trackball::Button GLUT2VCG (int glut_button, int )
{
	int vcgbt = vcg::Trackball::BUTTON_NONE;

	switch(glut_button){
		case GLUT_LEFT_BUTTON: vcgbt |= vcg::Trackball::BUTTON_LEFT;	break;
		case GLUT_MIDDLE_BUTTON: vcgbt |= vcg::Trackball::BUTTON_RIGHT;	break;
		case GLUT_RIGHT_BUTTON: vcgbt |= vcg::Trackball::BUTTON_MIDDLE;	break;
	}

	int modifiers = glutGetModifiers();

	if (modifiers & GLUT_ACTIVE_SHIFT)	vcgbt |= vcg::Trackball::KEY_SHIFT;
	if (modifiers & GLUT_ACTIVE_CTRL)	vcgbt |= vcg::Trackball::KEY_CTRL;
	if (modifiers & GLUT_ACTIVE_ALT)	vcgbt |= vcg::Trackball::KEY_ALT;

	return vcg::Trackball::Button (vcgbt);
}

 




void Display(){
	glViewport(0, 0, width, height);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if(!mesh.face.empty()){
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(40, width /(float) height , 0.1, 100);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(0,0,5,   0,0,0,   0,1,0);
		track.center=vcg::Point3f(0, 0, 0);
		track.radius= 1;
		track.GetView();
		track.Apply();
		glPushMatrix();
		float d=1.0f/mesh.bbox.Diag();
		vcg::glScale(d);
		glTranslate(-glWrap.m->bbox.Center());	
		// the trimesh drawing calls
		switch(drawmode)
		{
		  case SMOOTH: 
	  		glWrap.Draw<vcg::GLW::DMSmooth,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();
	  		break;
		  case PERPOINTS: 
	  		glWrap.Draw<vcg::GLW::DMPoints,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();
	  		break;
		  case WIRE: 
			glWrap.Draw<vcg::GLW::DMWire,     vcg::GLW::CMNone,vcg::GLW::TMNone> ();
			break;
		  case FLATWIRE: 
	  		glWrap.Draw<vcg::GLW::DMFlatWire, vcg::GLW::CMNone,vcg::GLW::TMNone> ();
	  		break;
		  case HIDDEN: 
	  		glWrap.Draw<vcg::GLW::DMHidden,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();
	  		break;
		  case FLAT: 
	  		glWrap.Draw<vcg::GLW::DMFlat,     vcg::GLW::CMNone,vcg::GLW::TMNone> ();
	  		break;
		  default: 
	  		break;
		}

		glPopMatrix();
		track.DrawPostApply();
	}


    TwDraw();

    // Present frame buffer
    glutSwapBuffers();

    // Recall Display at next frame
    glutPostRedisplay();
}


void Reshape(int _width,int _height){
	width =  _width ;
	height = _height;
   TwWindowSize(width, height);
}
void Terminate(){}

void initMesh()
{
	// update bounding box
	vcg::tri::UpdateBounding<CMesh>::Box(mesh);
	// update Normals
	vcg::tri::UpdateNormals<CMesh>::PerVertexNormalizedPerFace(mesh);
	vcg::tri::UpdateNormals<CMesh>::PerFaceNormalized(mesh);
	// Initialize the opengl wrapper
 	glWrap.m = &mesh;
  	glWrap.Update();
}

void  TW_CALL loadMesh(void *)
{	
	if(filename==0) return;
   int err=vcg::tri::io::ImporterPLY<CMesh>::Open(mesh,(char*)filename);
	if(err!=0){
	  const char* errmsg=vcg::tri::io::ImporterPLY<CMesh>::ErrorMsg(err);
	}
	else
		initMesh();
}

void  TW_CALL loadTetrahedron(void *){
	vcg::tri::Tetrahedron(mesh);
	initMesh();
}

void TW_CALL loadDodecahedron(void * ){
	vcg::tri::Dodecahedron(mesh);
	initMesh();
}

void   keyReleaseEvent (unsigned char k,int x,int y)
{
	int modifiers = glutGetModifiers();
	if (modifiers & GLUT_ACTIVE_CTRL)
	  track.ButtonUp (Trackball::Button::KEY_CTRL);
	if (modifiers & GLUT_ACTIVE_SHIFT)
	  track.ButtonUp ( Trackball::Button::KEY_SHIFT);
	if (modifiers & GLUT_ACTIVE_ALT)
	  track.ButtonUp (Trackball::Button::KEY_ALT);
}
void   keyPressEvent (unsigned char k,int x,int  y)
{
	int modifiers = glutGetModifiers();
	if (modifiers & GLUT_ACTIVE_CTRL)
	  track.ButtonDown (Trackball::Button::KEY_CTRL);
	if (modifiers & GLUT_ACTIVE_SHIFT)
	  track.ButtonDown ( Trackball::Button::KEY_SHIFT);
	if (modifiers & GLUT_ACTIVE_ALT)
	  track.ButtonDown (Trackball::Button::KEY_ALT);

	TwEventKeyboardGLUT(k,x,y);
}

 void mousePressEvent(int bt,int state,int x,int y){
	 if(state == GLUT_DOWN)
		 track.MouseDown ( x , height   -  y  , GLUT2VCG (bt,state));
	 else
		 track.MouseUp ( x  , height   -  y  , GLUT2VCG (bt,state));

	TwEventMouseButtonGLUT(bt,state,x,y);

  };

void mouseMoveEvent (int x, int y )
{
	if(!TwEventMouseMotionGLUT(x,y))
	    track.MouseMove ( x  , height  -  y  );
}

  //void mouseReleaseEvent(QMouseEvent*e);
void wheelEvent(int wheel, int direction, int x, int y ){
	track.MouseWheel(wheel*direction);
}


void TW_CALL CopyCDStringToClient(char **destPtr, const char *src)
{
    size_t srcLen = (src!=NULL) ? strlen(src) : 0;
    size_t destLen = (*destPtr!=NULL) ? strlen(*destPtr) : 0;

    // Alloc or realloc dest memory block if needed
    if( *destPtr==NULL )
        *destPtr = (char *)malloc(srcLen+1);
    else if( srcLen>destLen )
        *destPtr = (char *)realloc(*destPtr, srcLen+1);

    // Copy src
    if( srcLen>0 )
        strncpy(*destPtr, src, srcLen);
    (*destPtr)[srcLen] = '\0'; // null-terminated string
}

int main(int argc, char *argv[])
{

	TwBar *bar; // Pointer to the tweak bar

    // Initialize AntTweakBar
    // (note that AntTweakBar could also be intialized after GLUT, no matter)
    if( !TwInit(TW_OPENGL, NULL) )
    {
        // A fatal error occured    
        fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
        return 1;
    }

	glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(640, 480);
    glutCreateWindow("AntTweakBar simple example using GLUT");
    glutCreateMenu(NULL);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

    // Set GLUT callbacks
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    atexit(Terminate);  // Called after glutMainLoop ends

	    // Set GLUT event callbacks
    // - Directly redirect GLUT mouse button events to AntTweakBar
	glutMouseFunc((GLUTmousebuttonfun)mousePressEvent);
    // - Directly redirect GLUT mouse motion events to AntTweakBar
    glutMotionFunc((GLUTmousemotionfun)mouseMoveEvent);
    // - Directly redirect GLUT mouse "passive" motion events to AntTweakBar (same as MouseMotion)
    glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
    // - Directly redirect GLUT key events to AntTweakBar
    glutKeyboardFunc((GLUTkeyboardfun)TwEventKeyboardGLUT);
    // - Directly redirect GLUT special key events to AntTweakBar
    glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);

	glutKeyboardFunc(keyPressEvent);
	glutKeyboardUpFunc(keyReleaseEvent);

	 	
	glutMouseWheelFunc(wheelEvent);
    bar = TwNewBar("TweakBar");

	TwCopyCDStringToClientFunc (CopyCDStringToClient);
	
	TwAddVarRW(bar,"Input",TW_TYPE_CDSTRING,&filename," label='Filepath' group=SetMesh help=` Name of the file to load` ");
	TwAddButton(bar,"Load from file",loadMesh,0,	" label='Load Mesh' group=SetMesh help=`load the mesh` ");
	TwAddButton(bar,"Use tetrahedron",loadTetrahedron,0,	" label='Make Tetrahedron' group=SetMesh help=`use tetrahedron.` ");
	TwAddButton(bar,"Use dodecahedron",loadDodecahedron,0,	" label='Make Dodecahedron' group=SetMesh help=`use dodecahedron.` ");


	// ShapeEV associates Shape enum values with labels that will be displayed instead of enum values
	TwEnumVal drawmodes[6] = { {SMOOTH, "Smooth"}, {PERPOINTS, "Per Points"}, {WIRE, "Wire"}, {FLATWIRE, "FlatWire"},{HIDDEN, "Hidden"},{FLAT, "Flat"}};
	// Create a type for the enum shapeEV
	TwType drawMode = TwDefineEnum("DrawMode", drawmodes, 6);
	// add 'g_CurrentShape' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [<] and [>].
	TwAddVarRW(bar, "Draw Mode", drawMode, &drawmode, " keyIncr='<' keyDecr='>' help='Change draw mode.' ");

	glutMainLoop();
}
