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
#ifndef WINDOW_H_POS_DEMO
#define WINDOW_H_POS_DEMO

#include <QWidget>
#include <QPushButton>

class QSlider;
class GLWidget;

class Window : public QWidget
{
    Q_OBJECT

public:
    Window();

private:
		QPushButton *createButton(const char *changedSignal, const char *setterSlot);

    GLWidget *glWidget;
		QPushButton * fvButton,*feButton,*ffButton,*neButton,*ldButton,*nbButton,*vfButton;
};												 
													 
#endif
