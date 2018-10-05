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
#include "glarea.h"
#include <QGridLayout>
#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>
#include <QTime>

QProgressBar *MainWindow::qb;
QStatusBar* MainWindow::sb;

MainWindow::MainWindow (QWidget * parent)
	:QMainWindow (parent),mi(2000000000),mesh()
{
	ui.setupUi (this);
    initTable();
    QGridLayout* grid = new QGridLayout();
    //parent is set to NULL in order to avoid QT bug on MAC (business as usual...).
    //The QGLWidget are destroyed by hand in the MainWindow destructor...
    shared = new SharedDataOpenGLContext(mesh,mi,NULL);
	shared->setHidden(true);
    shared->myInitGL();

    for(unsigned int ii = 0;ii < 2;++ii)
	{
        rendbox[ii] = new QComboBox(this);
        for(QMap<MyDrawMode,QString>::iterator it = stringrendtable.begin();it != stringrendtable.end();++it)
            rendbox[ii]->addItem(it.value());
        rendbox[ii]->setCurrentIndex((int) MDM_SMOOTH);
        glar[ii] = new GLArea(shared,this);
        connect(rendbox[ii],SIGNAL(activated(int)),glar[ii],SIGNAL(updateRenderModalityRequested(int)));
        connect(glar[ii],SIGNAL(updateRenderModalityRequested(int)),this,SLOT(updateRenderModality(int)));
		//connect(shared,SIGNAL(dataReadyToBeRead(MyDrawMode,vcg::GLFeederInfo::ReqAtts&)),glar[ii], SLOT(updateRequested(MyDrawMode,vcg::GLFeederInfo::ReqAtts&)));
		grid->addWidget(rendbox[ii],0,ii,1,1,Qt::AlignRight);
        grid->addWidget(glar[ii],1,ii,1,1);
    }
    ui.glbox->setLayout(grid);

    connect(ui.actionLoad_Mesh,SIGNAL(triggered()),this,SLOT(chooseMesh()));
	connect (ui.actionLoad_Tetrahedron, SIGNAL (triggered()),this, SLOT (loadTetrahedron()));
    connect (ui.actionLoad_Dodecahedron, SIGNAL (triggered()),this, SLOT (loadDodecahedron()));

    sb = statusBar();
    qb=new QProgressBar(this);
    qb->setMaximum(100);
    qb->setMinimum(0);
    qb->reset();
    statusBar()->addPermanentWidget(qb,0);
}

// mesh chooser file dialog
void MainWindow::chooseMesh()
{
	mesh.Clear();
    QString plyext("ply");
    QString objext("obj");
    QString extoptions = QString("Poly Model (*.") + plyext + ");;OBJ Model (*." + objext + ")"; 
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open Mesh"), QDir::currentPath(),
		extoptions);
    QFileInfo fi(fileName);
    QTime loadingtime;
    loadingtime.start();
	int err=0;
    
    if (fi.suffix() == plyext)
        err=vcg::tri::io::ImporterPLY<CMeshO>::Open(mesh,(fileName.toStdString()).c_str(),qCallBack);
    else 
        if (fi.suffix() == objext)
        {
             int loadmask;
             err=vcg::tri::io::ImporterOBJ<CMeshO>::Open(mesh,(fileName.toStdString()).c_str(),loadmask,qCallBack);
        }
    int msec = loadingtime.elapsed();
	if(err!=0)
	{
		const char* errmsg=vcg::tri::io::ImporterPLY<CMeshO>::ErrorMsg(err);
		QMessageBox::warning(this,tr("Error Loading Mesh"),QString(errmsg));
	}
    QString msg = fileName + " vtx: " + QString::number(mesh.VN()) + " fcs: " + QString::number(mesh.FN()) + " loading time: " + QString::number(msec) + " msec"; 
	initMesh(msg);
}

void MainWindow::loadTetrahedron()
{
	mesh.Clear();
	vcg::tri::Tetrahedron(mesh);
	initMesh(tr("Tethraedron [builtin]"));
}

void MainWindow::loadDodecahedron()
{
	mesh.Clear();
	vcg::tri::Dodecahedron(mesh);
	initMesh(tr("Dodecahedron [builtin]"));
}

void MainWindow::initMesh(QString& message)
{
	if (shared != NULL)
		shared->deAllocateBO();
	// update bounding box
	vcg::tri::UpdateBounding<CMeshO>::Box(mesh);
	// update Normals
	vcg::tri::UpdateNormal<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh);
    QTime rdsetuptime;
    rdsetuptime.start();
    for(unsigned int ii = 0;ii < 2;++ii)
		if ((glar[ii] != NULL) && (rendbox[ii] != NULL))
        {
            glar[ii]->resetTrackBall();
            MyDrawMode mdm = (MyDrawMode) rendbox[ii]->currentIndex();
            QMap<MyDrawMode,QPair<vcg::GLMeshAttributesInfo::PRIMITIVE_MODALITY_MASK, vcg::GLMeshAttributesInfo::RendAtts> >::iterator it = rendtable.find(mdm);
            if (it == rendtable.end())
                return;
            shared->setPerViewRendAtts(glar[ii]->context(),it.value().first,it.value().second);
        }
    shared->manageBuffers();
    for(unsigned int ii = 0;ii < 2;++ii)
        glar[ii]->updateGL();
    qb->reset();
    int msec = rdsetuptime.elapsed();
    message += " bo creation: " + QString::number(msec) + " msec";
    statusBar()->showMessage(message);
    
}

void MainWindow::updateRenderModality(int clickedindex)
{
    GLArea* glasender = qobject_cast<GLArea*>(sender()); 
    if (glasender == NULL)
        return;
    updateRenderModality(glasender,clickedindex);
}

void MainWindow::updateRenderModality( GLArea* area,int clickedindex )
{
    if ((shared == NULL) || (area == NULL))
        return;
    MyDrawMode mdm = (MyDrawMode) clickedindex;
    QMap<MyDrawMode,QPair<vcg::GLMeshAttributesInfo::PRIMITIVE_MODALITY_MASK, vcg::GLMeshAttributesInfo::RendAtts> >::iterator it = rendtable.find(mdm);
    if (it == rendtable.end())
        return;
    shared->setPerViewRendAtts(area->context(),it.value().first,it.value().second);
    shared->manageBuffers();
    for(unsigned int ii = 0;ii < 2;++ii)
        glar[ii]->updateGL();
}

MainWindow::~MainWindow()
{
    for(int ii = 0;ii < 2;++ii)
		delete glar[ii];
	delete shared;
}

void MainWindow::initTable()
{
    stringrendtable[MDM_SMOOTH] = QString("Solid Smooth");
    stringrendtable[MDM_FLAT] = QString("Solid Flat");
    stringrendtable[MDM_WIRE] = QString("Wire");
    stringrendtable[MDM_POINTS] = QString("Points");
    stringrendtable[MDM_QUAD_WIRE] = QString("Wire Quad");
    stringrendtable[MDM_QUAD_SMOOTH_WIRE] = QString("Wire Solid Smooth Quad");

    vcg::GLMeshAttributesInfo::PRIMITIVE_MODALITY_MASK mmask = vcg::GLMeshAttributesInfo::PR_SOLID;
    vcg::GLMeshAttributesInfo::RendAtts ratts;
    ratts[vcg::GLMeshAttributesInfo::ATT_NAMES::ATT_VERTPOSITION] = true;
    ratts[vcg::GLMeshAttributesInfo::ATT_NAMES::ATT_VERTNORMAL] = true;
    rendtable[MDM_SMOOTH] = qMakePair(mmask,ratts);
    ratts.reset(true);
    ratts[vcg::GLMeshAttributesInfo::ATT_NAMES::ATT_FACENORMAL] = true;
    rendtable[MDM_FLAT] = qMakePair(mmask,ratts);
    ratts.reset(true);
    mmask = vcg::GLMeshAttributesInfo::PR_WIREFRAME_TRIANGLES;
    rendtable[MDM_WIRE] = qMakePair(mmask,ratts);
    ratts.reset(true);
    mmask = vcg::GLMeshAttributesInfo::PR_POINTS;
    ratts[vcg::GLMeshAttributesInfo::ATT_NAMES::ATT_VERTNORMAL] = true;
    rendtable[MDM_POINTS] = qMakePair(mmask,ratts);
    ratts.reset(true);
    mmask = vcg::GLMeshAttributesInfo::PR_WIREFRAME_EDGES;
    rendtable[MDM_QUAD_WIRE] = qMakePair(mmask,ratts);
    ratts.reset(true);
    mmask = vcg::GLMeshAttributesInfo::PR_WIREFRAME_EDGES | vcg::GLMeshAttributesInfo::PR_SOLID;
    ratts[vcg::GLMeshAttributesInfo::ATT_NAMES::ATT_VERTNORMAL] = true;
    rendtable[MDM_QUAD_SMOOTH_WIRE] = qMakePair(mmask,ratts);
}

bool MainWindow::qCallBack(const int pos, const char * str)
{
    int static lastPos=-1;
    if(pos==lastPos) return true;
    lastPos=pos;

    static QTime currTime = QTime::currentTime();
    if(currTime.elapsed()< 100) return true;
    currTime.start();
    sb->showMessage(str,5000);
    qb->show();
    qb->setEnabled(true);
    qb->setValue(pos);
    sb->update();
    qApp->processEvents();
    return true;
}

SharedDataOpenGLContext::MultiViewManager* MainWindow::getMultiviewerManager()
{
    if (shared == NULL)
        return NULL;
    return &(shared->manager);
}
