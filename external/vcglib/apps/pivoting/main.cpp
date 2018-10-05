#include <iostream>
using namespace std;

#include <QtGui/QApplication>
#include "mainwindow.h"

class QMyWindow: public QMainWindow {
   public:
     QMyWindow(QWidget *parent): QMainWindow(parent) {
       ui.setupUi(this);
       if(!connect(ui.add_face, SIGNAL(clicked()),
                   ui.area, SLOT(addFace()))) {
          cerr << "Could not connect addface\n";
       }

       if(!connect(ui.add_10_faces, SIGNAL(clicked()),
                   ui.area, SLOT(add10Faces()))) {
          cerr << "Could not connect addface\n";
       }

       if(!connect(ui.add_100_faces, SIGNAL(clicked()),
                   ui.area, SLOT(add100Faces()))) {
          cerr << "Could not connect addface\n";
       }

       if(!connect(ui.add_1000_faces, SIGNAL(clicked()),
                   ui.area, SLOT(add1000Faces()))) {
          cerr << "Could not connect addface\n";
       }
       if(!connect(ui.add_all, SIGNAL(clicked()),
                   ui.area, SLOT(addAll()))) {
          cerr << "Could not connect addface\n";
       }

       if(!connect(ui.add, SIGNAL(clicked()),
                   ui.area, SLOT(addTot()))) {
          cerr << "Could not connect addface\n";
       }


       if(!connect(ui.save, SIGNAL(clicked()),
                   ui.area, SLOT(save()))) {
          cerr << "Could not connect addface\n";
       }
         if(!connect(ui.smooth, SIGNAL(clicked(bool)),
                   ui.area, SLOT(viewSmooth(bool)))) {
          cerr << "Could not connect addface\n";
       }
       if(!connect(ui.open, SIGNAL(clicked()),
                   ui.area, SLOT(open()))) {
          cerr << "Could not connect addface\n";
       }
       if(!connect(ui.radius, SIGNAL(valueChanged(double)),
                   ui.area, SLOT(setRadius(double)))) {
          cerr << "Could not connect addface\n";
       }
       if(!connect(ui.tot, SIGNAL(valueChanged(int)),
                   ui.area, SLOT(setTot(int)))) {
          cerr << "Could not connect addface\n";
       }


  /*     connect(ui.world, SIGNAL(activated(const QString&)),
               ui.area, SLOT(loadWorld(const QString&)));

       connect(ui.shadow, SIGNAL(activated(int)),
               ui.area, SLOT(setShadowMode(int)));

       connect(ui.drawvolume, SIGNAL(clicked(bool)),
               ui.area, SLOT(setDrawVolume(bool)));

       connect(ui.track, SIGNAL(activated(int)),
               ui.area, SLOT(setTrackMode(int)));*/
     }   
   private:
     Ui::MainWindow ui;
};

int main(int argc, char *argv[]) {   
   QApplication app(argc, argv);
   QMyWindow *window = new QMyWindow(NULL);
   
   window->show();
   return app.exec();
}
