
#include <math.h>
#include <iostream>

#include <QtCore/QFile>
#include <QtGui/QMouseEvent>
#include <QtGui/QMessageBox>
#include <QtGui/QImage>
#include <QtGui/QFileDialog>

#include "glarea.h"
#include "cmesh.h"
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/flag.h>

#include <vcg/space/normal_extrapolation.h>
//#include "curvature.h"





using namespace std;
using namespace vcg;

GLArea::GLArea(QWidget *parent): QGLWidget(parent), pivot(NULL), smooth(false), radius(1.2) {
  
  tot = 1;

  setMouseTracking(true);  
//  init("uccello.ply", radius);
}       

bool GLArea::loadModel(const QString &file) {
    
     updateGL();
     return true;
}
void GLArea::open() {
  QString file = QFileDialog::getOpenFileName(this, "Select a ply file", "", "*.ply");
  if(!file.size()) return;
  init(file, radius);
}

void GLArea::init(QString file, float ballsize = 1.2) {

  int err = tri::io::Importer<CMesh>::Open(mesh, file.toAscii().data());  
  if(err) return;
  mesh.face.clear();
  mesh.fn = 0;

//  UpdateTopology<CMesh>::VertexFace(mesh);
//  UpdateTopology<CMesh>::FaceFace(mesh);
//  tri::UpdateFlags<CMesh>::FaceBorderFromFF(mesh);
//  tri::UpdateFlags<CMesh>::VertexBorderFromFace(mesh);

  //compute box;
  box = Box3f();  
  for(int i = 0; i < mesh.vert.size(); i++)
    box.Add(mesh.vert[i].P());
  
  
  float r = sqrt((box.Diag()*box.Diag())/mesh.vn);

//  mesh.face.clear();
//  mesh.fn = 0;      
  if(pivot) delete pivot;
  cout << "creating pibot\n";
  
  NormalExtrapolation<vector<CVertex> >::ExtrapolateNormals(mesh.vert.begin(), mesh.vert.end(), 10);
//  pivot = new Pivot<CMesh>(mesh, r*ballsize, 0.1, 0);
  pivot = new Pivot<CMesh>(mesh, 0, 0.1, 0);
//  pivot.build();
//  pivot->buildMesh();

}

void GLArea::save() {
    
     mesh.vn = mesh.vert.size();
     mesh.fn = mesh.face.size();
     tri::io::ExporterPLY<CMesh>::Save(mesh, "prova.ply");
     
}

void GLArea::addFace() { 
     pivot->addFace(); 
     updateGL(); 
/*     std::list<Hinge>::iterator li;
   for(li=pivot->front.begin();li!=pivot->front.end();++li)
     printf("(%d-%d-%d)",(*li).v0,(*li).v1,(*li).v2);
   printf("\n");*/
   
}

void GLArea::add10Faces() { 
     for(int i =0; i < 10; i++)
        if(-1 == pivot->addFace()) return;

     updateGL(); 
}


void GLArea::add100Faces() { 
     for(int i =0; i < 100; i++)
        if(-1 == pivot->addFace()) return;
     updateGL(); 
}

void GLArea::add1000Faces() { 
     for(int i =0; i < 1000; i++)
        if(-1 == pivot->addFace()) return;
     updateGL(); 
}

void GLArea::addAll() { 
  while(1) {
    for(int i = 0; i < 1000; i++) 
      if(0 > pivot->addFace()) return;
    updateGL(); 
  }
}

void GLArea::addTot() { 
  for(int i = 0; i < tot; i++) 
    if(0 > pivot->addFace()) return;
  updateGL(); 
}



void GLArea::viewSmooth(bool on) {
   smooth = on;
   updateGL();
}
void GLArea::initializeGL() { 
   glClearColor(1, 1, 1, 1); 
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_COLOR_MATERIAL);  

   glDisable(GL_BLEND);
   glEnable(GL_NORMALIZE);
   glDisable(GL_CULL_FACE);
   glCullFace(GL_BACK);
   glColor4f(1, 1, 1, 1);   

   glEnable(GL_LIGHTING);
   double st = 4; //1/sqrt(3);
   float lpos[4];
   lpos[0] = lpos[1] = lpos[2] = st;
   lpos[3] = 1;
   glLightfv(GL_LIGHT0, GL_POSITION, lpos);       
   
   float v[4] = {0.8, 0.8, 0.8, 0.0};
   glLightfv(GL_LIGHT0, GL_DIFFUSE, v);

   trackball.center=Point3f(0, 0, 0);
   trackball.radius= 1;

   glLoadIdentity();     
}


void GLArea::resizeGL(int w, int h) {   
  glViewport(0, 0, (GLint)w, (GLint)h);  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();  
  

  float r = w/(float)h;
  gluPerspective(60, r, 1, 4);

  glMatrixMode(GL_MODELVIEW);  
  
}


void GLArea::paintGL() {   
   
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);    
       
   glLoadIdentity();
   gluLookAt(0, 0, 3, 0, 0, 0, 0, 1, 0);        
   
   //Drawing the scene 
   glPushMatrix();
   trackball.GetView();
   trackball.Apply(true);

   Point3f c = -box.Center();
   float radius = 2.0f/box.Diag();
  
   if(mesh.face.size()>0) {
     CFace &face = mesh.face[0];
     CVertex *v[3];
     v[0] = face.V(0);
     c=-v[0]->P();
     radius=radius*5;
   }
   
   if(!pivot) return;

   glPushMatrix();
   glScalef(radius, radius, radius);
   glTranslatef(c[0], c[1], c[2]);
   
   if(pivot->front.size()>2)
   {
     glEnable(GL_LINE_SMOOTH);
     glColor4f(1, 0, 1, 0.1);
      glLineWidth(5);
      Pivot<CMesh>::Edgex &ee=pivot->front.front();
     int v0=ee.v0;
     int v1=ee.v1;
   glBegin(GL_LINES);
    glVertex3fv(mesh.vert[v0].P().V());
    glVertex3fv(mesh.vert[v1].P().V());
   glEnd();
  glLineWidth(1);
   }
   glEnable(GL_LIGHTING);
   glColor3f(0, 1, 0); 
        
        
           
   
   
   glBegin(GL_TRIANGLES);
  for(int i = 0; i < mesh.face.size(); i++) {
    CFace &face = mesh.face[i];
    CVertex *v[3];
    v[0] = face.V(0);
    v[1] = face.V(1);
    v[2] = face.V(2);

    Point3f n = (v[1]->P()- v[0]->P()) ^ (v[2]->P() - v[0]->P());
    //Point3f &n = face.N();
    glNormal3fv(&(n[0]));
  
    for(int k = 0; k < 3; k++) {      
      glVertex3fv((float *)&(v[k]->P()));
    }
  }
  glEnd();
  
  glEnable(GL_POLYGON_OFFSET_LINE);
  glPolygonOffset(-3, -3);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glDisable(GL_LIGHTING);  
  glColor3f(0, 0.5, 0);
  
  glPolygonOffset(-1, -1);
  glBegin(GL_TRIANGLES);
  for(int i = 0; i < mesh.face.size(); i++) {
    CFace &face = mesh.face[i];
    CVertex *v[3];
    v[0] = face.V(0);
    v[1] = face.V(1);
    v[2] = face.V(2);
    for(int k = 0; k < 3; k++) {      
      glVertex3fv((float *)&(v[k]->P()));
    }
  }
  glEnd();
  
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glDisable(GL_DEPTH_TEST);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
  for(list<Pivot<CMesh>::Edgex>::iterator k  = pivot->front.begin(); k != pivot->front.end(); k++) {
     glColor3f(1, 0, 0);                           
     Point3f &p0 = mesh.vert[(*k).v0].P();
     glVertex3fv(&(p0[0]));
     glColor3f(0, 0, 1);
     Point3f &p1 = mesh.vert[(*k).v1].P();
     glVertex3fv(&(p1[0]));
     
/*     glColor3f(1, 1, 0);
     Point3f middle = (mesh.vert[(*k).v0].P() + mesh.vert[(*k).v1].P())/2;
     glVertex3fv(&(middle[0]));
     glVertex3fv(&((*k).center[0]));*/
  }
  for(list<Pivot<CMesh>::Edgex>::iterator k  = pivot->deads.begin(); k != pivot->deads.end(); k++) {
     glColor3f(0, 0, 0);                           
     Point3f &p0 = mesh.vert[(*k).v0].P();
     glVertex3fv(&(p0[0]));
     Point3f &p1 = mesh.vert[(*k).v1].P();
     glVertex3fv(&(p1[0]));
  }
  glEnd();
  glEnable(GL_DEPTH_TEST);
  
  glPointSize(4.0f);
  glBegin(GL_POINTS);
  for(int i = 0; i < mesh.vert.size(); i++) {
    CVertex &v = mesh.vert[i];
    Point3f &p = v.P();
    if(v.IsD()) continue;
    if(v.IsV()) glColor3f(1, 0, 0);
    else if(v.IsB()) glColor3f(1, 1, 0);
    else continue;
    glVertex3f(p[0], p[1], p[2]);
  }
  glEnd();
  glColor3f(0, 0, 0);
  glPointSize(1.0f);
  
  glLineWidth(1.0f);
  glEnable(GL_LIGHTING);


 glBegin(GL_LINES);
  for(int i = 0; i < mesh.vert.size(); i++) {
    CVertex &v = mesh.vert[i];
    Point3f &p = v.P();
    if(v.IsD()) continue;
    glVertex3f(p[0], p[1], p[2]);
    Point3f q = p + v.N();
    glVertex3f(q[0], q[1], q[2]);
  }
  glEnd();
  
  
 


  glDisable(GL_POLYGON_OFFSET_LINE);
   
   
   
   

   glDisable(GL_LIGHTING);
   glPopMatrix();

   glScalef(radius, radius, radius);
   glTranslatef(c[0], c[1], c[2]); 
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//   cloud.draw();
   glColor3f(0.5, 0.5, 0.5);
   glPointSize(2);
   glBegin(GL_POINTS);
   for(int i = 0; i < mesh.vert.size(); i++) {
      CVertex &vert = mesh.vert[i];
      Point3f n = vert.N();
      Point3f p = vert.P();
      glNormal3f(n[0], n[1], n[2]);
      glVertex3f(p[0], p[1], p[2]);
   }
   glEnd();
   
   glDisable(GL_BLEND);
   glPopMatrix();
}

void GLArea::wheelEvent(QWheelEvent *e) {     
     if(e->delta() > 0)
       trackball.MouseWheel(1);
     else
       trackball.MouseWheel(-1);
     updateGL();
}

void GLArea::mouseMoveEvent(QMouseEvent *e) {
   trackball.MouseMove(e->x(), height() - e->y());
   updateGL();  
}
Trackball::Button QT2VCG(Qt::MouseButton qtbt,  Qt::KeyboardModifiers modifiers)
{
	int vcgbt=Trackball::BUTTON_NONE;
	if(qtbt & Qt::LeftButton		) vcgbt |= Trackball::BUTTON_LEFT;
	if(qtbt & Qt::RightButton		) vcgbt |= Trackball::BUTTON_RIGHT;
	if(qtbt & Qt::MidButton			) vcgbt |= Trackball::BUTTON_MIDDLE;
	if(modifiers & Qt::ShiftModifier		)	vcgbt |= Trackball::KEY_SHIFT;
	if(modifiers & Qt::ControlModifier ) vcgbt |= Trackball::KEY_CTRL;
	if(modifiers & Qt::AltModifier     ) vcgbt |= Trackball::KEY_ALT;
	return Trackball::Button(vcgbt);
}

void GLArea::keyReleaseEvent ( QKeyEvent * e )
{
      if(e->key()==Qt::Key_Control) trackball.MouseUp(0,0, QT2VCG(Qt::NoButton, Qt::ControlModifier ) );
      if(e->key()==Qt::Key_Shift) trackball.MouseUp(0,0, QT2VCG(Qt::NoButton, Qt::ShiftModifier ) );
      if(e->key()==Qt::Key_Alt) trackball.MouseUp(0,0, QT2VCG(Qt::NoButton, Qt::AltModifier ) );
}

void GLArea::mousePressEvent(QMouseEvent *e) {
  trackball.MouseDown(e->x(),height()-e->y(), QT2VCG(e->button(), e->modifiers() ) );
 //   if(e->button() == Qt::LeftButton)
       //trackball.MouseDown(e->x(), width() - e->y(), Trackball::BUTTON_LEFT);       
   // if(e->button() == Qt::RightButton)
     //   trackball.MouseDown(e->x(), width() - e->y(), Trackball::BUTTON_LEFT | Trackball::KEY_CTRL);       
updateGL();  
}

void GLArea::mouseReleaseEvent(QMouseEvent *e) { 
  trackball.MouseUp(e->x(),height()-e->y(), QT2VCG(e->button(), e->modifiers() ) );
  // if(e->button() == Qt::LeftButton)
      //trackball.MouseUp(e->x(), width() - e->y(), Trackball::BUTTON_LEFT);  
   //if(e->button() == Qt::RightButton)
    //  trackball.MouseUp(e->x(), width() - e->y(), Trackball::BUTTON_LEFT | Trackball::KEY_CTRL);                 
}
