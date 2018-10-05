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
Revision 1.4  2007/05/17 09:06:44  ganovelli
gestione double click

Revision 1.3  2006/12/10 23:29:57  ganovelli
added VFIterator (Pos is disabled in this version)

Revision 1.2  2006/12/10 22:17:18  ganovelli
cvs problem during frist committ. repeated

*/

#include <QtGui>
#include <GL/glew.h>
#include <QtOpenGL>

#include <math.h>
#include "glwidget.h"
#include <wrap/io_trimesh/import_PLY.h>
#include <wrap/gl/picking.h>
#include <wrap/gl/space.h>
#include <wrap/gl/pos.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/topology.h>


GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(parent)
{
    object = 0;

    trolltechGreen = QColor::fromCmykF(0.40, 0.0, 1.0, 0.0);
    trolltechPurple = QColor::fromCmykF(0.39, 0.39, 0.0, 0.0);
		track.SetPosition(vcg::Point3f(0.0,0.0,0.0));
		track.SetIdentity();
		track.radius = 0.4;
		pos.f=NULL;
		vfite.f = NULL;
		doPickVfIte = false;
		doPickPos = false;
}

GLWidget::~GLWidget()
{
    makeCurrent();
    glDeleteLists(object, 1);
}

QSize GLWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}

QSize GLWidget::sizeHint() const
{
    return QSize(400, 400);
}

void GLWidget::LoadTriMesh(QString &namefile)
{
	vcg::tri::io::ImporterPLY<MyStraightMesh>::Open(mesh,namefile.toAscii());
	vcg::tri::UpdateBounding<MyStraightMesh>::Box(mesh);
	vcg::tri::UpdateNormals<MyStraightMesh>::PerFace(mesh);
	vcg::tri::UpdateNormals<MyStraightMesh>::PerVertex(mesh);
	vcg::tri::UpdateTopology<MyStraightMesh>::FaceFace(mesh);
	vcg::tri::UpdateTopology<MyStraightMesh>::VertexFace(mesh);
	pos.f=0;
	vfite.f=NULL;
}

void GLWidget::OpenFile(){
	QStringList filters;
	

	QString	fileName = QFileDialog::getOpenFileName(this,tr("Open File"),".", filters.join("\n"));
	
	if (fileName.isEmpty())	return;
	else
	 LoadTriMesh( fileName );

}

void GLWidget::flipV( ){
	if(pos.f) pos.FlipV();
	repaint();
}
void GLWidget::flipE( ){
	if(pos.f) pos.FlipE();
	repaint();
}
void GLWidget::flipF( ){
	if(pos.f) pos.FlipF();
	repaint();
}
void GLWidget::nextE( ){
	if(pos.f) pos.NextE();
	repaint();
}
void GLWidget::nextB( ){
	if(pos.f) pos.NextB();
	repaint();
}

void GLWidget::nextVfite( ){
	if(vfite.F()) ++vfite;
	repaint();
}

void GLWidget::initializeGL()
{
    qglClearColor(trolltechPurple.dark());
    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
}

template <class VertexType>
void drawVertex(VertexType & v){
	glPushAttrib(0xffffffff);
	glPointSize(2.0);
	glBegin(GL_POINTS);
	glVertex(v.P());
	glEnd();
	glPopAttrib();
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,2,   0,0,0,   0,1,0);    
    
 		track.GetView();
    track.Apply();

		
		glScalef(1/glWrap.m->bbox.Diag(),1/glWrap.m->bbox.Diag(),1/glWrap.m->bbox.Diag());
		glTranslate(-glWrap.m->bbox.Center());	

		// to do some picking
		 MyStraightMesh::FaceType* fp=NULL;
		 MyStraightMesh::VertexType* vp=NULL;
		 	if(doPickPos)
		 	{
				std::vector<MyStraightMesh::FaceType*> res;
				int yes = vcg::Pick<MyStraightMesh::FaceContainer>(pic_x,ScreenH-pic_y+1,mesh.face,res,vcg::glTriangle3<MyStraightMesh::FaceType>,1,1);
				if(yes) 
					{fp = res[0];
						pos.Set(fp,0,fp->V(0));
					}
				doPickPos=false;
		  }else
		 	if(doPickVfIte)
			{
				std::vector<MyStraightMesh::VertexType*> res;
				int yes = vcg::Pick<MyStraightMesh::VertContainer>(pic_x,ScreenH-pic_y+1,mesh.vert,res,drawVertex<MyStraightMesh::VertexType>,3,3);
				if(yes) 
					{vp = res[0];
MyStraightMesh::FaceType* g  = vp->VFp();
						vfite=vcg::face::VFIterator<MyStraightMesh::FaceType>(vp);
					}

				doPickVfIte = false;
			}
		
		glWrap.Draw<vcg::GLW::DMFlatWire,vcg::GLW::CMNone,vcg::GLW::TMNone> ();
 
		if(pos.f!=NULL) {
			glPushAttrib(0xffffffff);
			glDisable(GL_LIGHTING);
			glColor3f(0.0,1.0,0.0);
			glDepthRange(0.0,0.999);
			vcg::GlPos<vcg::face::Pos<MyStraightMesh::FaceType> >::Draw(pos);
			glPopAttrib();
		}
		if(vfite.F()!=NULL) {
			glPushAttrib(0xffffffff);
			glDisable(GL_LIGHTING);
			glColor3f(0.0,1.0,0.0);
			glDepthRange(0.0,0.99);
			vcg::GlVfIterator<vcg::face::VFIterator<MyStraightMesh::FaceType> >::Draw(vfite);
			glPopAttrib();
		}

}





void GLWidget::mouseMoveEvent(QMouseEvent *e)
{
			track.MouseMove(e->x(),ScreenH-e->y()+1);
			repaint();

    //if (event->buttons() & Qt::LeftButton) {
    //    setXRotation(xRot + 8 * dy);
    //} else if (event->buttons() & Qt::RightButton) {
    //    setXRotation(xRot + 8 * dy);
    //}
 //   lastPos = event->pos();
}
 void GLWidget::keyPressEvent ( QKeyEvent * e ) {
		if((keypress == Qt::Key_Control)&&(e->key()==Qt::Key_Control))
			keypress = -1;
		else
			keypress = e->key();
 }
 void GLWidget::mouseDoubleClickEvent(QMouseEvent *e){
					if(e->button() == Qt::RightButton)
						doPickPos=true;

 }
 void GLWidget:: mousePressEvent(QMouseEvent *e)
{
		if( (keypress==Qt::Key_Control) && (e->button() == Qt::LeftButton) )
					track.MouseDown(e->x(),ScreenH-e->y()+1,vcg::Trackball::KEY_CTRL|vcg::Trackball::BUTTON_LEFT );
			else
				if(e->button() == Qt::LeftButton )
					track.MouseDown(e->x(),ScreenH-e->y()+1,vcg::Trackball::BUTTON_LEFT);
				else
					if(e->button() == Qt::RightButton)
					{
						doPickVfIte=true;
						pic_x = e->x();
						pic_y = e->y();
					}
			repaint();
		}

 void GLWidget::mouseReleaseEvent ( QMouseEvent * e ){
			
			if( (keypress==Qt::Key_Control) && (e->button() == Qt::LeftButton) )
					track.MouseUp(e->x(),ScreenH-e->y()+1,vcg::Trackball::KEY_CTRL );
			else
				if(e->button() == Qt::LeftButton )
					track.MouseUp(e->x(),ScreenH-e->y()+1,vcg::Trackball::BUTTON_LEFT);
			repaint();
		}


 void GLWidget::resizeGL(int w,int h){
			ScreenW=w; ScreenH=h;
			glViewport(0,0,w,h);
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(45,ScreenW/(float)ScreenH,0.01,5);
		}
 void GLWidget::wheelEvent ( QWheelEvent * e ){
			int v =  e->delta()/(float) 120;
			track.MouseWheel(v);
			repaint();
		}
