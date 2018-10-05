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
Revision 1.10  2005/11/12 06:48:47  cignoni
Version 1.0
Added management of point set, correct bug in printing on the screen,

Revision 1.9  2005/01/03 13:59:54  cignoni
Resolved min/max macro conflict

Revision 1.8  2004/09/30 00:57:42  ponchio
Reference to temporary fixed and indented.

Revision 1.7  2004/09/28 09:46:51  cignoni
Added MapFalseColor

Revision 1.6  2004/09/16 14:08:35  ponchio
gamma is a math function.

Revision 1.5  2004/09/10 14:02:20  cignoni
Added Cone directions

Revision 1.4  2004/09/09 22:34:38  cignoni
Integrated lost modifications...

Revision 1.3  2004/09/09 14:35:54  ponchio
Various changes for gcc compatibility

Revision 1.2  2004/07/11 22:13:30  cignoni
Added GPL comments


****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>
// this define is mandatory to avoid the conflicts due to the silly definition of 
// min and max macros in windows.h (included by glut...)
#define NOMINMAX 

#include <GL/glew.h>
#include <GL/glut.h>
#include <wrap/gl/space.h>


#include <wrap/callback.h>
#include <vcg/math/base.h>
#include <wrap/gui/trackball.h>
#include <vcg/simplex/vertex/with/vcvn.h>
#include <vcg/simplex/vertex/with/vcvn.h>
#include <vcg/simplex/face/with/fn.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/complex.h>
#include<wrap/io_trimesh/export_ply.h>
#include<wrap/io_trimesh/import_ply.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/bounding.h>
#include<vcg/complex/algorithms/update/color.h>

#include "visshader.h"
using namespace vcg;
using namespace std;




// Vertex, Face, Mesh and Grid definitions.
class MyEdge;
class AFace;
class AVertex   : public VertexVCVN< float ,MyEdge,AFace > {};
class AFace     : public FaceFN< AVertex,MyEdge,AFace > {};
class AMesh     : public tri::TriMesh< vector<AVertex>, vector<AFace> > {};

///////// Global ////////

int SampleNum=64;
int WindowRes=800;
unsigned int TexInd=0;
bool SwapFlag=false;
bool CullFlag=false;
bool ClosedFlag=false;
Point3f ConeDir(0,1,0);
float ConeAngleRad = math::ToRad(180.0f);

float lopass=0,hipass=1,gamma_correction=1;
float diff=.8;
float ambi=.2;

bool LightFlag=true;
bool ColorFlag=true;
bool FalseColorFlag=false;
bool ShowDirFlag=false;
int imgcnt=0;

Color4b BaseColor=Color4b::White;
Trackball QV;
Trackball QL;
Trackball *Q=&QV;

int ScreenH,ScreenW;
float ViewAngle=33;
vector<Point3f> ViewVector;

bool cb(const char *buf)
{
  printf(buf);
  return true;
}

void BuildOnePixelTexture(Color4b c, unsigned int &TexInd)
{
  if(TexInd==0) glGenTextures(1,&TexInd);

  glBindTexture(GL_TEXTURE_1D,TexInd);
  glTexImage1D(GL_TEXTURE_1D,0,GL_RGBA,1,0,GL_RGBA,GL_UNSIGNED_BYTE,&c);
  glEnable(GL_TEXTURE_1D);
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
}
void glutPrintf(int x, int y, const char * f, ... )
{
  glMatrixMode (GL_PROJECTION);   
  glPushMatrix();
  glLoadIdentity (); 
  glOrtho(0,ScreenW,0,ScreenH,-1,1);
  glMatrixMode (GL_MODELVIEW);    
  glPushMatrix();
  glLoadIdentity ();  
  
  int len, i;
  char buf[4096];
  va_list marker;
  va_start( marker, f );     

  int n = vsprintf(buf,f,marker);
  va_end( marker );              
  glPushAttrib(GL_ENABLE_BIT);
  glDisable(GL_LIGHTING);
  glColor3f(0,0,0);
  glRasterPos2f(x, y);
  len = (int) strlen(buf);
  for (i = 0; i < len; i++) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, buf[i]);
  }
  glMatrixMode (GL_PROJECTION);   
  glPopMatrix();
  glMatrixMode (GL_MODELVIEW);    
  glPopMatrix();
  glPopAttrib();  
}

// prototypes
void SaveTexturedGround();

void DrawViewVector()
{
  glDisable(GL_LIGHTING);  
  glColor3f(0,0,1);
  glBegin(GL_LINES);
  for(unsigned int i=0;i<ViewVector.size();++i)
    {
      glVertex3f(0,0,0);glVertex(ViewVector[i]);
    }
  glEnd();
}
void DrawLightVector()
{
  const int sz=5;
  glPushMatrix();
  QL.Apply();
  glDisable(GL_LIGHTING);  
  glBegin(GL_LINES);
  glColor3f(1,1,0);
  for(unsigned int i=0;i<=sz;++i)
    for(unsigned int j=0;j<=sz;++j)
      {
        glColor3f(1,1,0);
        glVertex3f(-1.0f+i*2.0/sz,-1.0f+j*2.0/sz,-1);
        glVertex3f(-1.0f+i*2.0/sz,-1.0f+j*2.0/sz, 1);
      }
  glEnd();
  glPopMatrix();
}

void Draw(AMesh &mm)
{
  if(mm.face.empty())
  { 
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    AMesh::VertexIterator vi;
    glBegin(GL_POINTS);
    for(vi=mm.vert.begin();vi!=mm.vert.end();++vi)
      {
        if(ColorFlag) glColor((*vi).C()); 
        glVertex((*vi).P());
      }
    glEnd();
    glPopAttrib();
  }
  else
  {
    AMesh::FaceIterator fi;
    glBegin(GL_TRIANGLES);
    for(fi=mm.face.begin();fi!=mm.face.end();++fi)
      {
        glNormal((*fi).V(0)->N()); if(ColorFlag) glColor((*fi).V(0)->C());  glVertex((*fi).V(0)->P());
        glNormal((*fi).V(1)->N()); if(ColorFlag) glColor((*fi).V(1)->C());  glVertex((*fi).V(1)->P());
        glNormal((*fi).V(2)->N()); if(ColorFlag) glColor((*fi).V(2)->C());  glVertex((*fi).V(2)->P());
      }
    glEnd();
  }
}

AMesh m;
VertexVisShader<AMesh> Vis(m);

string OutNameMsh;

	
/*  Called when the window is first opened and whenever 
 *  the window is reconfigured (moved or resized).
 */
void  ViewReshape(GLsizei w, GLsizei h)
{
  ScreenW=w; ScreenH=h;
  glViewport(0,0,w,h);
}

void  ViewDisplay (void)
{
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode (GL_PROJECTION);   
  glLoadIdentity (); 
  gluPerspective(ViewAngle,(float)ScreenW/ScreenH,1,7);
  glMatrixMode (GL_MODELVIEW);    
  glLoadIdentity ();  
  glPushMatrix();
  QL.Apply();
  glutPrintf(5,5,"Diffuse %04.2f   Ambient %04.2f "
    "   LowPass %04.2f    HiPass %04.2f    Gamma %04.2f    rgb = %03i:%03i:%03i",
	     diff,ambi,lopass,hipass,gamma_correction,BaseColor[0],BaseColor[1],BaseColor[2]);

  GLfloat light_position0[] = {0.0, 10.0, 300.0, 0.0};
  glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
  glPopMatrix();
  glTranslatef(0,0,-4);
  if(Q==&QL) DrawLightVector();	
  QL.GetView();
  QV.GetView();
  QV.Apply();
  if(ShowDirFlag) DrawViewVector();

  float d = 2.0/m.bbox.Diag();
  glScalef(d, d, d);
  glColor3f(diff,diff,diff); 
  glTranslate(-m.bbox.Center());
  if(LightFlag) glEnable(GL_LIGHTING);
  else glDisable(GL_LIGHTING);
  if(ColorFlag) glEnable(GL_COLOR_MATERIAL);
  else glDisable(GL_COLOR_MATERIAL);
  if(FalseColorFlag) glColorMaterial(GL_FRONT,GL_DIFFUSE);
  else  glColorMaterial(GL_FRONT,GL_AMBIENT);
  
  glMateriali(GL_FRONT,GL_SHININESS,0);
  float spec[4]={0,0,0,1};
  float ambientV[4]={ambi,ambi,ambi,1};
  float diffuseV[4]={diff,diff,diff,1};
  glMaterialfv(GL_FRONT,GL_SPECULAR,spec);
  glMaterialfv(GL_FRONT,GL_AMBIENT, ambientV);
  glMaterialfv(GL_FRONT,GL_DIFFUSE, diffuseV);
  glCullFace(GL_BACK);
    
  if(CullFlag) glEnable(GL_CULL_FACE);
  else        glDisable(GL_CULL_FACE);
 
  BuildOnePixelTexture(BaseColor,TexInd);
  Draw(m);
  glutSwapBuffers();
}
	
void ViewSpecialKey(int , int , int )
{
  glutPostRedisplay();
}
void Toggle(bool &flag) {flag = !flag;}
void UpdateVis()
{
  if( LightFlag && !FalseColorFlag) 
    Vis.MapVisibility(gamma_correction,lopass,hipass,ambi);
  if(!LightFlag && !FalseColorFlag) 
    Vis.MapVisibility(gamma_correction,lopass,hipass,1.0);
  if(FalseColorFlag)  
    Vis.MapFalseColor();
}
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

void ViewKey(unsigned char key, int , int )
{
  Point3f dir;
  switch (key) {
  case 27: exit(0);   	break;
  case 9: if(Q==&QV) Q=&QL;else Q=&QV;   	break;
  case 'l' :
    lopass=lopass+.05; 
    printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,gamma_correction); 
    UpdateVis(); break;
  case 'L' :
    lopass=lopass-.05; 
    printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,gamma_correction); 
    UpdateVis(); break;
  case 'h' :
    hipass=hipass-.05; 
    printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,gamma_correction); 
    UpdateVis(); break;
  case 'H' :
    hipass=hipass+.05; 
    printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,gamma_correction); 
    UpdateVis(); break;
  case 'd' :  diff+=.05; printf("Ambient %f Diffuse %f, \n",ambi,diff); 		UpdateVis(); break;
  case 'D' :  diff-=.05; printf("Ambient %f Diffuse %f, \n",ambi,diff); 		UpdateVis(); break;
  case 'a' :  ambi+=.05; printf("Ambient %f Diffuse %f, \n",ambi,diff); 		UpdateVis(); break;
  case 'A' :  ambi-=.05; printf("Ambient %f Diffuse %f, \n",ambi,diff); 		UpdateVis(); break;
    
  case 'e' :  ambi+=.05; diff-=.05; 
    printf("Ambient %f Diffuse %f, \n",ambi,diff); 
    UpdateVis(); break;
  case 'E' :  ambi-=.05; diff+=.05; 
    printf("Ambient %f Diffuse %f, \n",ambi,diff); 
    UpdateVis(); break;
  case 'p' :
    gamma_correction=gamma_correction-.05; 
    printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,gamma_correction); 
    UpdateVis(); break;
  case 'P' :
    gamma_correction=gamma_correction+.05; 
    printf("Lo %f, Hi %f Gamma %f\n",lopass,hipass,gamma_correction); 
    UpdateVis(); break;
  case 13 : 
    //Vis.ComputeUniform(SampleNum,ViewVector,cb); 
    Vis.ComputeUniformCone(SampleNum,ViewVector, ConeAngleRad,ConeDir,cb); 
    UpdateVis(); break;
  case ' ' : {
    Point3f dir = Q->camera.ViewPoint();
    printf("ViewPoint %f %f %f\n",dir[0],dir[1],dir[2]);
    dir.Normalize();
    dir=Inverse(Q->track.Matrix())*dir;
    printf("ViewPoint %f %f %f\n",dir[0],dir[1],dir[2]);
    dir.Normalize();
    Vis.ComputeSingle(dir,ViewVector,cb); 
    UpdateVis(); 
  } break;
  case 'r' : BaseColor[0]=min(255,BaseColor[0]+2);     break;
  case 'R' : BaseColor[0]=max(  0,BaseColor[0]-2);     break;
  case 'g' : BaseColor[1]=min(255,BaseColor[1]+2);     break;
  case 'G' : BaseColor[1]=max(  0,BaseColor[1]-2);     break;
  case 'b' : BaseColor[2]=min(255,BaseColor[2]+2);     break;
  case 'B' : BaseColor[2]=max(  0,BaseColor[2]-2);     break;
    
  case 'v' : Toggle(ShowDirFlag); break;
  case 'V' : 
    {
      SimplePic<Color4b> snapC;
      snapC.OpenGLSnap();
      char buf[128];
      sprintf(buf,"Snap%03i.ppm",imgcnt++);
      snapC.SavePPM(buf);
    }
  case 's' :
    Vis.SmoothVisibility();
    UpdateVis(); break;
  case 't' :
    Vis.SmoothVisibility(true);
    UpdateVis(); break;
  case 'S' :
    { 
      vcg::tri::io::PlyInfo p; 
      p.mask|=vcg::tri::io::Mask::IOM_VERTCOLOR  /* | vcg::ply::PLYMask::PM_VERTQUALITY*/ ;
      tri::io::ExporterPLY<AMesh>::Save(m,OutNameMsh.c_str(),false,p);
    }
    break;
  case 'C' : LightFlag = !LightFlag; 
    printf("Toggled Light %s\n",LightFlag?"on":"off"); 		
    UpdateVis(); break;
  case 'c' : ColorFlag = !ColorFlag; 
    printf("Toggled Color %s\n",ColorFlag?"on":"off"); break;
  case 'f' : FalseColorFlag = !FalseColorFlag; 
    printf("Toggled FalseColor %s\n",ColorFlag?"on":"off"); 
    UpdateVis(); break;
    
  case '1' : diff=0.80f; ambi=0.10f; gamma_correction=1.0; lopass=0.00f; hipass=1.00f; ColorFlag=false; UpdateVis(); break;
  case '2' : diff=0.65f; ambi=0.30f; gamma_correction=1.0; lopass=0.15f; hipass=0.80f; ColorFlag=true;  UpdateVis(); break;
  case '3' : diff=0.45f; ambi=0.50f; gamma_correction=1.0; lopass=0.20f; hipass=0.75f; ColorFlag=true;  UpdateVis(); break;
  case '4' : diff=0.35f; ambi=0.60f; gamma_correction=1.0; lopass=0.25f; hipass=0.70f; ColorFlag=true;  UpdateVis(); break;
  }
  glutPostRedisplay(); ;
} 
void ViewMenu(int val)
{
  ViewKey(val, 0, 0); 
}
/*********************************************************************/
// TrackBall Functions
/*********************************************************************/

int GW,GH; // Grandezza della finestra

void ViewMouse(int button, int state, int x, int y)
{
  static int KeyMod=0;
  static int glut_buttons=0;
  //printf("ViewMouse %i %i %i %i\n",x,y,button,state);
  int m_mask = 0;
  if(state == GLUT_DOWN) {
    KeyMod=glutGetModifiers();
    if(GLUT_ACTIVE_SHIFT & KeyMod)		m_mask |=  Trackball::KEY_SHIFT;
    if(GLUT_ACTIVE_ALT & KeyMod)			m_mask |=  Trackball::KEY_ALT;
    if(GLUT_ACTIVE_CTRL & KeyMod)			m_mask |=  Trackball::KEY_CTRL;

    glut_buttons |= (1<<button);
    Q->MouseDown(x, ScreenH-y, glut_buttons | m_mask);
  } else {
    if(GLUT_ACTIVE_SHIFT & KeyMod)		m_mask |=  Trackball::KEY_SHIFT;
    if(GLUT_ACTIVE_ALT & KeyMod)			m_mask |=  Trackball::KEY_ALT;
    if(GLUT_ACTIVE_CTRL & KeyMod)			m_mask |=  Trackball::KEY_CTRL;
    glut_buttons |= (1<<button);
    Q->MouseUp(x, ScreenH-y, glut_buttons | m_mask);
  }
}

void ViewMouseMotion(int x, int y)
{
  Q->MouseMove(x,ScreenH-y);
  glutPostRedisplay();
}

void SetLight()
{
  GLfloat light_ambient0[] = {0.0, 0.0, 0.0, 1.0};  
  GLfloat light_diffuse0[] = {1.0, 1.0, 1.0, 1.0};  
  GLfloat light_position0[] = {0.0, 10.0, 300.0, 0.0};
  glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light_diffuse0);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
  glEnable(GL_LIGHT0);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT,light_ambient0);
  
}
 
void  ViewInit (void) {
  SetLight();
  Q->Reset();
  Q->radius= 1;
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glClearColor (0.8, 0.8, 0.8, 0.0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING);
  
  //	glEnable(GL_BLEND);
  glShadeModel(GL_SMOOTH);
  //  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);     
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);	
}



int main(int argc, char** argv)
{
  if(argc<2) {
    printf(
	   "shadevis 1.0 \n"__DATE__"\n"
	   "Copyright 2003-2004 Visual Computing Lab I.S.T.I. C.N.R.\n"
	   "Paolo Cignoni (cignoni@isti.cnr.it)\n\n"
	   "Usage: shadevis file.ply [options]\n"
	   "Options:\n"
	   "     -w#      WindowResolution (default 600)\n"
	   "     -n#      Sample Directions (default 64)\n"
	   "     -z#      z offset (default 1e-3)\n"
	   "     -c       assume that the mesh is closed (slightly faster, default false)\n"
	   "     -f       Flip normal of the model\n"
	   "     -da #    Cone Direction Angle in degree (default 180)\n"
	   "     -dv # # # Cone Direction vector (default 0 0 1)\n"			
	   );
    
    return 1;
  }
  
  srand(time(0));
  int i=1;  
  while(i<argc 	&& (argv[i][0]=='-'))
    {
      switch(argv[i][1])
	{
	case 'd'  : 
	  if(argv[i][2] == 'a') { 
	    ConeAngleRad = math::ToRad(atof(argv[i+1])); ++i; break; 
	  }
	  if(argv[i][2] == 'v') { 
	    Point3f p(atof(argv[i+1]),atof(argv[i+2]),atof(argv[i+3]));
	    ConeDir = Normalize(p); 
	    i+=3; break; 
	  }
	  printf("Error unable to parse option '%s'\n",argv[i]); 
	  exit(0);
	  break;
	case 'n'  : SampleNum = atoi(argv[i]+2); break;
	case 'f'  : SwapFlag=false; break;
	case 'c'  : ClosedFlag=true; break;
        case 'w'  : WindowRes= atoi(argv[i]+2); 
	  printf("Set WindowRes to %i\n",WindowRes ); break;
        case 's'  : Vis.SplitNum= atoi(argv[i]+2); 
	  printf("Set SplitNum to %i\n",Vis.SplitNum ); break;
        case 'z'  : Vis.ZTWIST = atof(argv[i]+2); 
	  printf("Set ZTWIST to %f\n",Vis.ZTWIST ); break;
        default: {
	  printf("Error unable to parse option '%s'\n",argv[i]); 
	  exit(0);
	}
	}
      
      ++i;
    }
  
  
  string basename = argv[i];
  if(!(basename.substr(basename.length()-4)==".ply"))	{
    printf("Error: Unknown file extension %s\n",basename.c_str());
    return 1;
  }
  
  // loading original mesh
  int ret=tri::io::ImporterPLY<AMesh>::Open(m,argv[i]);
  if(ret) {printf("Error unable to open mesh %s : '%s' \n",argv[i],tri::io::ImporterPLY<AMesh>::ErrorMsg(ret));exit(-1);}
  tri::UpdateNormals<AMesh>::PerVertexNormalized(m);
  tri::UpdateBounding<AMesh>::Box(m);
  tri::UpdateColor<AMesh>::VertexConstant(m,Color4b::White);
  Vis.IsClosedFlag=ClosedFlag;
  Vis.Init();
  UpdateVis();

  printf("Mesh bbox (%f %f %f)-(%f %f %f)\n\n",
	 m.bbox.min[0],m.bbox.min[1],m.bbox.min[2],
	 m.bbox.max[0],m.bbox.max[1],m.bbox.max[2]);

  OutNameMsh=(string(argv[i]).substr(0,strlen(argv[i])-4));
  OutNameMsh+="_vis.ply";
  
  printf("Mesh       Output filename %s\n",OutNameMsh.c_str());
  
  printf("Mesh %iv %if bbox Diag %g\n",m.vn,m.fn,m.bbox.Diag());
  
  glutInit(&argc, argv);
  
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(WindowRes, WindowRes);
  glutInitWindowPosition (10,10);
  glutCreateWindow ("shadevis - Visual Computing Lab - vcg.isti.cnr.it ");
  glutDisplayFunc(ViewDisplay);
  glutReshapeFunc(ViewReshape); 
  glutKeyboardFunc(ViewKey);	
  glutSpecialFunc(ViewSpecialKey);	
  glutMouseFunc(ViewMouse);
  glutMotionFunc(ViewMouseMotion);
  
  ViewInit();	
  glewInit();	
  glutMainLoop();
  
  return(0);
}
