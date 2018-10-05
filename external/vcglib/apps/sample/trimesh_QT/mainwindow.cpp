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


#include "mainwindow.h"
#include <QFileDialog>

MainWindow::MainWindow (QWidget * parent):QMainWindow (parent)
{
  ui.setupUi (this);
  //connections

  //from toolFrame to glArea
  connect (ui.drawModeComboBox, SIGNAL (currentIndexChanged(int)),
           ui.glArea, SLOT (selectDrawMode(int)));

  connect (ui.loadTetrahedronPushButton, SIGNAL (clicked()),
           ui.glArea, SLOT (loadTetrahedron()));

  connect (ui.loadDodecahedronPushButton, SIGNAL (clicked()),
           ui.glArea, SLOT (loadDodecahedron()));

 //from toolFrame to glArea through mainwindow
  connect (ui.loadMeshPushButton, SIGNAL (clicked()),
           this, SLOT (chooseMesh()));
  connect (this, SIGNAL (loadMesh(QString)),
           ui.glArea, SLOT(loadMesh(QString)));

  //from glArea to statusbar
  connect (ui.glArea, SIGNAL (setStatusBar(QString)),
           ui.statusbar, SLOT (showMessage(QString)));
}

// mesh chooser file dialog
void MainWindow::chooseMesh()
{
  QString fileName = QFileDialog::getOpenFileName(this,
                     tr("Open Mesh"), QDir::currentPath(),
                     tr("Poly Model (*.ply)"));
  if(!fileName.isEmpty())
    emit loadMesh(fileName);
}
