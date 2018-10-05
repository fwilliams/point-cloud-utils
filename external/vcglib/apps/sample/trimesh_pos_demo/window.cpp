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
Revision 1.2  2006/12/10 22:17:18  ganovelli
cvs problem during frist committ. repeated

*/
#include <QtGui>

#include "glwidget.h"
#include "window.h"
#include "mesh_type.h"


Window::Window()
{
    glWidget = new GLWidget;

		fvButton = createButton("FlipV()",SLOT(flipV( )));
		feButton = createButton("FlipE()",SLOT(flipE( )));
		ffButton = createButton("FlipF()",SLOT(flipF( )));
		neButton = createButton("NextE() {FlipE() + FlipF() }",SLOT(nextE( )));
		nbButton = createButton("NextB() ",SLOT(nextB( )));
		ldButton = createButton("Load TriMesh",SLOT(OpenFile( )));
		vfButton = createButton("++()",SLOT(nextVfite()));

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(glWidget);
    mainLayout->addWidget(fvButton);
    mainLayout->addWidget(feButton);
    mainLayout->addWidget(ffButton);
    mainLayout->addWidget(neButton);
    mainLayout->addWidget(nbButton);
    mainLayout->addWidget(vfButton);
    mainLayout->addWidget(ldButton);
    setLayout(mainLayout);

		glWidget->glWrap.m = &glWidget->mesh; 

    setWindowTitle(tr("TriMesh Pos Demo"));
}


QPushButton *Window::createButton(const char *text, const char *setterSlot)
{
    QPushButton *button = new QPushButton( text,0);
		button-> resize ( 50, 20 );

    if(!connect(button, SIGNAL(clicked()), glWidget, setterSlot))
			exit(0);
    return button;
}
