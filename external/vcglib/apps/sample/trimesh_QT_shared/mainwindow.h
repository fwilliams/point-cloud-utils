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

****************************************************************************/
#ifndef MAINWINDOW_H_
#define MAINWINDOW_H_

#include <GL/glew.h>
#include <QGLContext>
#include "ui_mainwindow.h"

#include "mesh.h"
#include <wrap/gl/gl_mesh_attributes_info.h>
#include <wrap/qt/qt_thread_safe_mesh_attributes_multi_viewer_bo_manager.h>
#include <wrap/qt/qt_thread_safe_memory_info.h>
#include <QProgressBar>
#include <QStatusBar>
#include <QComboBox>
#include "glarea.h"
enum MyDrawMode{MDM_SMOOTH=0,MDM_POINTS,MDM_WIRE,MDM_FLAT,MDM_QUAD_WIRE,MDM_QUAD_SMOOTH_WIRE};



class MainWindow:public QMainWindow
{
Q_OBJECT 
public:
   MainWindow(QWidget * parent = 0);
  ~MainWindow();

  CMeshO& currentMesh()  {return mesh;}
  SharedDataOpenGLContext::MultiViewManager* getMultiviewerManager();
  void updateRenderModality(GLArea* area,int clickedindex);
  static bool qCallBack(const int pos, const char * str);
public slots:
  void chooseMesh();
  void loadTetrahedron();
  void loadDodecahedron();
  void initMesh(QString& message);
  void updateRenderModality(int clickedindex);
signals:
  void loadMesh(QString newMesh);
  void updateRenderModalityRequested(int);
private:
  void initTable();
  Ui::mainWindow ui;
  GLArea* glar[2];
  SharedDataOpenGLContext* shared;
  vcg::QtThreadSafeMemoryInfo mi;
  QComboBox* rendbox[2];
  /// the active mesh instance
  CMeshO mesh;
  QMap<MyDrawMode,QString> stringrendtable;
  QMap<MyDrawMode,QPair<vcg::GLMeshAttributesInfo::PRIMITIVE_MODALITY_MASK,vcg::GLMeshAttributesInfo::RendAtts> > rendtable;
  static QProgressBar* qb;
  static QStatusBar* sb;
};

#endif /*MAINWINDOW_H_ */
