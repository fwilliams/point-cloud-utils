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
#include <QKeyEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#include <wrap/qt/trackball.h>
#include <cassert>
#include <wrap/gl/trimesh.h>
#include "mainwindow.h"

GLArea::GLArea (SharedDataOpenGLContext* sharedcontext,MainWindow* parent)
    :QGLWidget(NULL,sharedcontext),parwin(parent)
{
        
}

GLArea::~GLArea()
{

}

void GLArea::initializeGL()
{
	makeCurrent();
	glewExperimental=GL_TRUE;
	glewInit();
	glClearColor(0, 0, 0, 0);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

    GLenum err = glGetError();
    assert(err == GL_NO_ERROR);
    doneCurrent();
}

void GLArea::resizeGL (int w, int h)
{
	makeCurrent();
	glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    doneCurrent();
    //initializeGL();
}

void GLArea::paintGL ()
{
    if (parwin == NULL)
        return;
    SharedDataOpenGLContext::MultiViewManager* man = parwin->getMultiviewerManager();
    if (man == NULL)
        return;
    CMeshO& mesh = parwin->currentMesh();
    //glt.m = &mesh;
	makeCurrent();
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluPerspective(25, GLArea::width()/(float)GLArea::height(), 0.1, 100);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	gluLookAt(0,0,5,   0,0,0,   0,1,0);
	track.center=vcg::Point3f(0, 0, 0);
	track.radius= 1;
    track.GetView();
    track.Apply();
    glPushMatrix();
    if (mesh.VN() > 0)
    {
        float d=2.0f/mesh.bbox.Diag();
        vcg::glScale(d);
	    glTranslate(-mesh.bbox.Center());
        man->draw(context());
    }
	glPopMatrix();
	
	track.DrawPostApply();
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glPopAttrib();
    
	GLenum err = glGetError();
	assert(err == GL_NO_ERROR);
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
    makeCurrent();
	updateGL ();
    doneCurrent();
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
    makeCurrent();
	updateGL ();
    doneCurrent();
}

void GLArea::mousePressEvent (QMouseEvent * e)
{
	e->accept ();
	setFocus ();
	track.MouseDown (QT2VCG_X(this,e), QT2VCG_Y(this,e), QT2VCG (e->button (), e->modifiers ()));
    makeCurrent();
	updateGL ();
    doneCurrent();
}

void GLArea::mouseMoveEvent (QMouseEvent * e)
{
	if (e->buttons ()) {
		track.MouseMove (QT2VCG_X(this,e), QT2VCG_Y(this,e));
        makeCurrent();
		updateGL ();
        doneCurrent();
	}
}

void GLArea::mouseReleaseEvent (QMouseEvent * e)
{
	track.MouseUp (QT2VCG_X(this,e), QT2VCG_Y(this,e), QT2VCG (e->button (), e->modifiers ()));
    makeCurrent();
	updateGL ();
    doneCurrent();
}

void GLArea::wheelEvent (QWheelEvent * e)
{
	const int WHEEL_STEP = 120;
	track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
    makeCurrent();
	updateGL ();
    doneCurrent();
}

void GLArea::resetTrackBall()
{
	track.Reset();
	makeCurrent();
    updateGL();
    doneCurrent();
}



SharedDataOpenGLContext::SharedDataOpenGLContext(CMeshO& mesh,vcg::QtThreadSafeMemoryInfo& mi,QWidget* parent)
    :QGLWidget(parent),manager(mesh,mi,100000)
{
}

SharedDataOpenGLContext::~SharedDataOpenGLContext()
{
}

void SharedDataOpenGLContext::myInitGL()
{
    makeCurrent();
    glewInit();
    doneCurrent();
}

void SharedDataOpenGLContext::deAllocateBO()
{
    makeCurrent();
    manager.removeAllViewsAndDeallocateBO();
    doneCurrent();
}

void SharedDataOpenGLContext::setPerViewRendAtts( QGLContext* viewid,vcg::GLMeshAttributesInfo::PRIMITIVE_MODALITY_MASK mm,vcg::GLMeshAttributesInfo::RendAtts& atts )
{
    manager.setPerViewInfo(viewid,mm,atts);
}

void SharedDataOpenGLContext::manageBuffers()
{
    makeCurrent();
    manager.manageBuffers();
    doneCurrent();
}