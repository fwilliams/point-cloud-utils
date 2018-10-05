/****************************************************************************
 * VCGLib                                                            o o     *
 * Visual and Computer Graphics Library                            o     o   *
 *                                                                _   O  _   *
 * Copyright(C) 2009                                                \/)\/    *
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

#include <QtCore>
#include "wrap/qt/img_qt.h"

void sample001_open_save_color(QString input_file, QString output_file)
{
  img::Image<> image;
  img::openQtRGB(input_file, image);
  img::saveQtRGB(image, output_file);
}

void sample002_open_save_grayscale(QString input_file, QString output_file)
{
  img::Image<1> image;
  img::openQtY(input_file, image);
  img::saveQtY(image, output_file);
}

void sample003_normalize(QString input_file, QString output_file)
{
  img::Image<> image;
  img::openQtRGB(input_file, image);
  img::saveQtRGB(img::getNormalized(image), output_file);
}


void sample004_boxfilter(QString input_file, QString output_file)
{
  img::Image<> image;
  img::openQtRGB(input_file, image);
  int radius=3;
  img::saveQtRGB(img::getBoxFiltered(image,radius), output_file);
}

void sample005_gaussiansmooth(QString input_file, QString output_file)
{
  img::Image<> image;
  img::openQtRGB(input_file, image);
  int radius=4;
  img::saveQtRGB(img::getGaussianSmoothed(image,radius), output_file);
}

void sample006_medianfilter(QString input_file, QString output_file)
{
  img::Image<> image;
  img::openQtRGB(input_file, image);
  int radius=5;
  img::saveQtRGB(img::getMedianFiltered(image,radius), output_file);
}

void sample007_unsharpmask(QString input_file, QString output_file)
{
  img::Image<> image;
  img::openQtRGB(input_file, image);

  int radius=4;
  double factor=0.6;

  img::saveQtRGB(img::getUnsharpMasked(image,radius,factor), output_file);
}


void sample008_laplacianfilter(QString input_file, QString output_file)
{
  img::Image<> image;
  img::openQtRGB(input_file, image);
  img::Image<> laplacianfiltered;
  img::LaplacianFilter(image,laplacianfiltered);
  img::saveQtRGB(img::getNormalized(laplacianfiltered), output_file);
}

void sample009_logfilter(QString input_file, QString output_file)
{
  img::Image<> image;
  img::openQtRGB(input_file, image);
  int radius=5;
  img::Image<> logfiltered;
  img::LoGFilter(image,logfiltered,radius);
  img::saveQtRGB(img::getNormalized(logfiltered), output_file);
}

void sample010_dogfilter(QString input_file, QString output_file)
{
  img::Image<> image;
  img::openQtRGB(input_file, image);
  // must be radius1 < radius2
  int radius1=2; 
  int radius2=4;
  img::Image<> dogfiltered;
  img::DoGFilter(image,dogfiltered,radius1,radius2);

  img::saveQtRGB(img::getNormalized(dogfiltered), output_file);
}

void sample011_general_convolutions(QString input_file, QString output_dir,QString output_suffix)
{
  img::Image<> image;
  img::openQtRGB(input_file, image);

  QVector< QPair< double*, QPair< QPair< int, int > , QString> > > mm;
  double *f;

  f=new double[9];
  f[0]= 0.0f; f[1]= 0.0f; f[2]= 0.0f;
  f[3]=-1.0f; f[4]= 1.0f; f[5]= 0.0f;
  f[6]= 0.0f, f[7]= 0.0f; f[8]= 0.0f;
  mm.push_back(qMakePair(f,qMakePair(qMakePair(3,3),QString("edge_enhance"))));    

  f=new double[9];
  f[0]= 2.0f; f[1]= 0.0f; f[2]= 0.0f;
  f[3]= 0.0f; f[4]=-1.0f; f[5]= 0.0f;
  f[6]= 0.0f, f[7]= 0.0f; f[8]=-1.0f;
  mm.push_back(qMakePair(f,qMakePair(qMakePair(3,3),QString("embross"))));    

  f=new double[15];
  f[0] = 1.0f; f[1] = 2.0f; f[2] =  0.0f; f[3] = 2.0f; f[4] = 1.0f;
  f[5] = 1.0f; f[6] = 2.0f; f[7] =-18.0f; f[8] = 2.0f; f[9] = 1.0f;
  f[10]= 1.0f; f[11]= 2.0f; f[12]=  0.0f; f[13]= 2.0f; f[14]= 1.0f;
  mm.push_back(qMakePair(f,qMakePair(qMakePair(5,3),QString("my_vert_edges"))));

  f=new double[15];
  f[0] = 1.0f; f[1] =  1.0f; f[2] = 1.0f; 
  f[3] = 2.0f; f[4] =  2.0f; f[5] = 2.0f; 
  f[6] = 0.0f; f[7] =-18.0f; f[8] = 0.0f; 
  f[9] = 2.0f; f[10]=  2.0f; f[11]= 2.0f; 
  f[12]= 1.0f; f[13]=  1.0f; f[14]= 1.0f;
  mm.push_back(qMakePair(f,qMakePair(qMakePair(3,5),QString("my_horiz_edges"))));

  QPair< double*, QPair< QPair< int, int > , QString> > m;

  foreach(m,mm){
    double* matrix=m.first;
    int matrix_width=((m.second).first).first;
    int matrix_height=((m.second).first).second;
    QString matrix_name=(m.second).second;

    img::Image<> convolved;
    img::convolution(image,convolved,matrix,matrix_width,matrix_height);

    delete [] matrix; 

    bool normalize=(img::minValue(convolved)<0.0f)||(img::maxValue(convolved)>=255.0f);
    QString output_file(output_dir+"/011_general_convolution_"+matrix_name+
                        "_"+(normalize?" normalized_":"")+output_suffix);

    if(normalize)
      img::saveQtRGB(img::getNormalized(convolved),output_file);
    else
      img::saveQtRGB(convolved,output_file);
  }
}

void img_filters(QString input_dir,QString image,QString output_dir)
{
  QString input_file(input_dir+"/"+image);
  sample001_open_save_color(input_file, output_dir+"/001-open_save_color_"+image);

  sample002_open_save_grayscale(input_file, output_dir+"/002-open_save_grayscale_"+image);

  sample003_normalize(input_file, output_dir+"/003_normalize_"+image);

  sample004_boxfilter(input_file, output_dir+"/004_boxfilter_"+image);

  sample005_gaussiansmooth(input_file, output_dir+"/005_gaussiansmooth_"+image);

  sample006_medianfilter(input_file, output_dir+"/006_medianfilter_"+image);

  sample007_unsharpmask(input_file, output_dir+"/007_unsharpmask_"+image);

  sample008_laplacianfilter(input_file, output_dir+"/008_laplacianfilter_normalized_"+image);

  sample009_logfilter(input_file, output_dir+"/009_logfilter_normalized_"+image);

  sample010_dogfilter(input_file, output_dir+"/010_dogfilter_normalized_"+image);

  sample011_general_convolutions(input_file, output_dir,image);
}

bool clean_dir(QDir dir); // utility, unrelated with the sample 

int main(int argc,char ** argv)
{
  if(argc<3) 
  {
    printf("Usage: img_filters <input_dir> <output_dir>\n");
    return 1;
  }
  printf("Executing img_filters over all images in \"%s\", ouput is in \"%s\"\n",  argv[1], argv[2]);

  QString input_dir(argv[1]);
  QString output_dir(argv[2]);

  QStringList readable_image_extensions = QStringList()
       << "*.bmp" << "*.gif" << "*.jpg" << "*.jpeg"
       << "*.png" << "*.pbm" << "*.pgm" << "*.ppm"
       << "*.tiff" << "*.xbm" << "*.xpm";

  QStringList image_list = QDir(input_dir).entryList(readable_image_extensions,QDir::Files|QDir::Readable,QDir::Name);
  assert(clean_dir(QDir(output_dir)));

  try {
    foreach(QString image, image_list)
      img_filters(input_dir,image,output_dir);  
  } catch (img::ImageException& e) {
    qDebug() << "caught ImageException, message:" << e.message;
  }
  return 0;
}

bool clean_dir(QDir dir){ // utility, unrelated with the sample
  if(!dir.exists()){
    qDebug() << QString("dir \"%1\" does not exists\n").arg(dir.path());
    return false;
  }
  foreach(QString e,dir.entryList(QDir::NoDotAndDotDot|QDir::Dirs|QDir::Files)){
    QFileInfo i(QString("%1/%2").arg(dir.path(),e));
    if(i.isDir()){
      if(!clean_dir(QDir(QString("%1/%2").arg(dir.path(),e)))){
        qDebug() << QString("cannot clean \"%1/%2\"\n").arg(dir.path(),e);
        return false;
      }
      if(!dir.rmdir(e)){
        qDebug() << QString("cannot remove \"%1/%2\"\n").arg(dir.path(),e);
        return false;
      }
    }else{
      if(!dir.remove(e)){
        qDebug() << QString("cannot remove \"%1/%2\"\n").arg(dir.path(),e);
       return false;
      }
    }
  }
  return true;
}
