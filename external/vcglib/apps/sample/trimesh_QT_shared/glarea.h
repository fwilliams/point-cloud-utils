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

#ifndef GLAREA_H_
#define GLAREA_H_

/// Opengl related imports
#include <GL/glew.h>
#include <QGLWidget>

/// wrapper imports
#include <wrap/gui/trackball.h>
#include <wrap/qt/qt_thread_safe_memory_info.h>
#include <wrap/qt/qt_thread_safe_mesh_attributes_multi_viewer_bo_manager.h>
#include <wrap/gl/trimesh.h>
#include "mesh.h"

class MainWindow;

class SharedDataOpenGLContext : public QGLWidget
{
	Q_OBJECT
public:
    typedef vcg::QtThreadSafeGLMeshAttributesMultiViewerBOManager<CMeshO,QGLContext*> MultiViewManager;

	SharedDataOpenGLContext(CMeshO& mesh,vcg::QtThreadSafeMemoryInfo& mi,QWidget* parent = 0);
	~SharedDataOpenGLContext();

    void myInitGL();
	void deAllocateBO();
    void setPerViewRendAtts(QGLContext* view,vcg::GLMeshAttributesInfo::PRIMITIVE_MODALITY_MASK mm,vcg::GLMeshAttributesInfo::RendAtts& atts);
    void manageBuffers();
	MultiViewManager manager;
};


class GLArea:public QGLWidget
{
	Q_OBJECT 
public:
    GLArea (SharedDataOpenGLContext* sharedcontext,MainWindow* parent);
	~GLArea();
	void resetTrackBall();
    //unsigned int getId() const {return areaid;}
	/// we choosed a subset of the avaible drawing modes

signals:
		/// signal for setting the statusbar message
	void setStatusBar(QString message);
    void updateRenderModalityRequested(int);
protected:
	/// opengl initialization and drawing calls
	void initializeGL ();
	void resizeGL (int w, int h);
	void paintGL ();
	/// keyboard and mouse event callbacks
	void keyReleaseEvent(QKeyEvent * e);
	void keyPressEvent(QKeyEvent * e);
	void mousePressEvent(QMouseEvent*e);
	void mouseMoveEvent(QMouseEvent*e);
	void mouseReleaseEvent(QMouseEvent*e);
	void wheelEvent(QWheelEvent*e); 

private:
    MainWindow* parwin;
	/// the active manipulator
	vcg::Trackball track;
	/// mesh data structure initializer
	void initMesh(QString message);
    //unsigned int areaid;
    //vcg::GlTrimesh<CMeshO> glt;
};

#endif /*GLAREA_H_ */
