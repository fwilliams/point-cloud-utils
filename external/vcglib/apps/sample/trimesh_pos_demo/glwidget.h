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
#ifndef GLWIDGET_H_POS_DEMO
#define GLWIDGET_H_POS_DEMO

#include <GL/glew.h>
#include <QGLWidget>
#include "mesh_type.h"
#include <wrap/gui/trackball.h>
#include <wrap/gl/trimesh.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/topology.h>

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = 0);
    ~GLWidget();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    int xRotation() const { return xRot; }
		MyStraightMesh mesh;
		vcg::GlTrimesh<MyStraightMesh> glWrap;
		vcg::Trackball track;
		int ScreenH,ScreenW,pic_x,pic_y,keypress;
		bool doPickPos,doPickVfIte;
		vcg::face::Pos< MyStraightMesh::FaceType> pos;
		vcg::face::VFIterator< MyStraightMesh::FaceType> vfite;

public slots:
    void flipV( );
    void flipE( );
    void flipF( );
    void nextE( );
    void nextB( );
    void nextVfite( );

		void LoadTriMesh(QString& namefile);
		void OpenFile();

protected:
		void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mouseDoubleClickEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
		void wheelEvent ( QWheelEvent * e );
		void keyPressEvent ( QKeyEvent * e );

private:
    GLuint makeObject();
    void quad(GLdouble x1, GLdouble y1, GLdouble x2, GLdouble y2,
              GLdouble x3, GLdouble y3, GLdouble x4, GLdouble y4);
    void extrude(GLdouble x1, GLdouble y1, GLdouble x2, GLdouble y2);
    void normalizeAngle(int *angle);

    GLuint object;
    int xRot;
    QPoint lastPos;
    QColor trolltechGreen;
    QColor trolltechPurple;
};

#endif
