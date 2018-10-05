/****************************************************************************
**
** Copyright (C) 2011 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/

#include "glwidget.h"
#include <wrap/qt/trackball.h>
#include <wrap/gl/picking.h>
#include <wrap/qt/anttweakbarMapper.h>

TwBar *bar;
char * filename;/// filename of the mesh to load
CMesh mesh;     /// the active mesh instance
vcg::GlTrimesh<CMesh> glWrap;    /// the active mesh opengl wrapper
vcg::Trackball track;     /// the active manipulator
GLW::DrawMode drawmode=GLW::DMFlatWire;     /// the current drawmode

void  TW_CALL loadTetrahedron(void *){
    vcg::tri::Tetrahedron(mesh);
    vcg::tri::UpdateBounding<CMesh>::Box(mesh);
    vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(mesh);
    vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(mesh);
    glWrap.m = &mesh;
    glWrap.Update();
}

void TW_CALL loadMesh(void *)
{
  if(filename==0) return;
  int err=vcg::tri::io::ImporterPLY<CMesh>::Open(mesh,(char*)filename);
  if(err==ply::E_NOERROR)
  {
    vcg::tri::UpdateBounding<CMesh>::Box(mesh);
    vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(mesh);
    vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(mesh);
    glWrap.m = &mesh;
    glWrap.Update();
  }
}

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
  filename=0;
  hasToPick=false;
  setWindowTitle(tr("Hello GL"));
  bar = TwNewBar("TweakBar");
  TwCopyCDStringToClientFunc (CopyCDStringToClient);

  TwAddVarRW(bar,"Input",TW_TYPE_CDSTRING, &filename," label='Filepath' group=SetMesh help=` Name of the file to load` ");
  TwAddButton(bar,"Load from file",loadMesh,0,	" label='Load Mesh' group=SetMesh help=`load the mesh` ");
  TwAddButton(bar,"Use tetrahedron",loadTetrahedron,0,	" label='Make Tetrahedron' group=SetMesh help=`use tetrahedron.` ");

  // ShapeEV associates Shape enum values with labels that will be displayed instead of enum values
  TwEnumVal drawmodes[6] = { {GLW::DMSmooth, "Smooth"}, {GLW::DMPoints, "Per Points"}, {GLW::DMWire, "Wire"}, {GLW::DMFlatWire, "FlatWire"},{GLW::DMHidden, "Hidden"},{GLW::DMFlat, "Flat"}};
  // Create a type for the enum shapeEV
  TwType drawMode = TwDefineEnum("DrawMode", drawmodes, 6);
  // add 'g_CurrentShape' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [<] and [>].
  TwAddVarRW(bar, "Draw Mode", drawMode, &drawmode, " keyIncr='<' keyDecr='>' help='Change draw mode.' ");
}

void GLWidget::initializeGL ()
{
  glewInit();
  glClearColor(0, 0, 0, 0);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
}

void GLWidget::resizeGL (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  TwWindowSize(w, h);
  initializeGL();
 }

void GLWidget::paintGL ()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, GLWidget::width()/(float)GLWidget::height(), 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3.5f,   0,0,0,   0,1,0);
    track.center=vcg::Point3f(0, 0, 0);
    track.radius= 1;
    track.GetView();
    glPushMatrix();
    track.Apply();
    glPushMatrix();
    if(mesh.vert.size()>0)
    {
      vcg::glScale(2.0f/mesh.bbox.Diag());
      glTranslate(-mesh.bbox.Center());
      glWrap.Draw(GLW::DrawMode(drawmode),GLW::CMNone,GLW::TMNone);
    }
    glPopMatrix();
    track.DrawPostApply();
    glPopMatrix();
    if(hasToPick)
    {
      hasToPick=false;
      Point3f pp;
      if(Pick<Point3f>(pointToPick[0],pointToPick[1],pp))
      {
        track.Translate(-pp);
        track.Scale(1.25f);
        QCursor::setPos(mapToGlobal(QPoint(width()/2+2,height()/2+2)));
      }
    }
    TwDraw();
}

void GLWidget::keyReleaseEvent (QKeyEvent * e)
{
  e->ignore ();
  if (e->key () == Qt::Key_Control)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
  if (e->key () == Qt::Key_Shift)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
  if (e->key () == Qt::Key_Alt) track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
  updateGL ();
}


void GLWidget::keyPressEvent (QKeyEvent * e)
{
  e->ignore ();
  if (e->key () == Qt::Key_Control) track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
  if (e->key () == Qt::Key_Shift)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
  if (e->key () == Qt::Key_Alt)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));

  TwKeyPressQt(e);
  updateGL ();
}

void GLWidget::mousePressEvent (QMouseEvent * e)
{
  if(!TwMousePressQt(e))
  {
  e->accept ();
  setFocus ();
  track.MouseDown (e->x (), height () - e->y (), QT2VCG (e->button (), e->modifiers ()));
  }
  updateGL ();
}

void GLWidget::mouseMoveEvent (QMouseEvent * e)
{
  if (e->buttons ()) {
    track.MouseMove (e->x (), height () - e->y ());
    updateGL ();
  }
  TwMouseMotion(e->x (), e->y ());
}

void GLWidget::mouseDoubleClickEvent (QMouseEvent * e)
{
  hasToPick=true;
  pointToPick=Point2i(e->x(),height()-e->y());
  updateGL ();
}


void GLWidget::mouseReleaseEvent (QMouseEvent * e)
{
  track.MouseUp (e->x (), height () - e->y (), QT2VCG (e->button (), e->modifiers ()));
  TwMouseReleaseQt(e);
  updateGL ();
}

void GLWidget::wheelEvent (QWheelEvent * e)
{
  const int WHEEL_STEP = 120;
  track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
  updateGL ();
}
