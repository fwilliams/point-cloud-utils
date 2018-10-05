#ifndef CG_RENDERAREA_H
#define CG_RENDERAREA_H

#include <iostream>

#include <QtOpenGL/QGLWidget>

#include <vcg/space/box3.h>
#include <wrap/gui/trackball.h>
#include <vcg/complex/algorithms/create/ball_pivoting.h>
#include "cmesh.h"


//#include "mls_surface.h"
//#include "advancing.h"

using namespace tri;
using namespace trimesh;
class GLArea : public QGLWidget {
    Q_OBJECT

public:               
                  
    GLArea(QWidget *parent = 0);
    bool smooth;
   
public slots:         
    
    bool loadModel(const QString &file);
    void addFace();
    void add10Faces();
    void add100Faces();
    void add1000Faces();
    void addAll();
    void setTot(int n) { tot = n; }
    void addTot();
    
    void open();
    void save();
    void setRadius(double _radius) {
      radius = _radius;
    }
    void viewSmooth(bool on);

protected:
    int tot;      
          
    float radius;
    vcg::Box3f box;
    CMesh mesh;
//    Pivot<CMesh> *pivot;
    Pivot<CMesh> *pivot;
    vcg::Trackball trackball;

    void init(QString file, float radius);          
    void draw();
             

    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    void wheelEvent(QWheelEvent *e);
    void mouseMoveEvent(QMouseEvent *e);
    void mousePressEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);        
    void keyReleaseEvent ( QKeyEvent * e );

};

#endif
