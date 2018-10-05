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
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.1  2007/10/18 08:52:06  benedetti
Initial release.


****************************************************************************/

#include "glarea.h"
#include <QMessageBox>
#include <QKeyEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#include <wrap/qt/trackball.h>

GLArea::GLArea (QWidget * parent)
          :QGLWidget (parent)
{
    drawmode= SMOOTH;
    GLArea::loadTetrahedron();
}

void GLArea::selectDrawMode(int mode){
    drawmode=DrawMode(mode);
    updateGL();
}

void GLArea::loadMesh(QString fileName)
{
   int err=vcg::tri::io::ImporterPLY<CMesh>::Open(mesh,(fileName.toStdString()).c_str());
    if(err!=0){
      const char* errmsg=vcg::tri::io::ImporterPLY<CMesh>::ErrorMsg(err);
          QMessageBox::warning(this,tr("Error Loading Mesh"),QString(errmsg));
    }
    initMesh("Loaded \""+fileName+"\".");
}

void GLArea::loadTetrahedron(){
    vcg::tri::Tetrahedron(mesh);
    initMesh(tr("Tethraedron [builtin]"));
}

void GLArea::loadDodecahedron(){
    vcg::tri::Dodecahedron(mesh);
    initMesh(tr("Dodecahedron [builtin]"));
}

void GLArea::initMesh(QString message)
{
    // update bounding box
    vcg::tri::UpdateBounding<CMesh>::Box(mesh);
    // update Normals
        vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(mesh);
        vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(mesh);
    // Initialize the opengl wrapper
    glWrap.m = &mesh;
    glWrap.Update();
    updateGL();
    emit setStatusBar(message);
}

void GLArea::initializeGL ()
{
  glClearColor(0, 0, 0, 0);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
}

void GLArea::resizeGL (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  initializeGL();
 }

void GLArea::paintGL ()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(25, GLArea::width()/(float)GLArea::height(), 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,5,   0,0,0,   0,1,0);
    track.center=vcg::Point3f(0, 0, 0);
    track.radius= 1;
    track.GetView();
    track.Apply();
    glPushMatrix();
    float d=2.0f/mesh.bbox.Diag();
    vcg::glScale(d);
    glTranslate(-glWrap.m->bbox.Center());
    // the trimesh drawing calls
    switch(drawmode)
    {
      case SMOOTH:
        glWrap.Draw<vcg::GLW::DMSmooth,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();
        break;
      case POINTS:
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

void GLArea::keyReleaseEvent (QKeyEvent * e)
{
  e->ignore ();
  if (e->key () == Qt::Key_Control)
    track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
  if (e->key () == Qt::Key_Shift)
    track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
  if (e->key () == Qt::Key_Alt)
    track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
  updateGL ();
}

void GLArea::keyPressEvent (QKeyEvent * e)
{
  e->ignore ();
  if (e->key () == Qt::Key_Control)
    track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
  if (e->key () == Qt::Key_Shift)
    track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
  if (e->key () == Qt::Key_Alt)
    track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));
  updateGL ();
}

void GLArea::mousePressEvent (QMouseEvent * e)
{
  e->accept ();
  setFocus ();
  track.MouseDown (QT2VCG_X(this,e), QT2VCG_Y(this,e), QT2VCG (e->button (), e->modifiers ()));
  updateGL ();
}

void GLArea::mouseMoveEvent (QMouseEvent * e)
{
  if (e->buttons ()) {
    track.MouseMove (QT2VCG_X(this,e), QT2VCG_Y(this,e));
    updateGL ();
  }
}

void GLArea::mouseReleaseEvent (QMouseEvent * e)
{
  track.MouseUp (QT2VCG_X(this,e), QT2VCG_Y(this,e), QT2VCG (e->button (), e->modifiers ()));
  updateGL ();
}

void GLArea::wheelEvent (QWheelEvent * e)
{
  const int WHEEL_STEP = 120;
  track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
  updateGL ();
}
