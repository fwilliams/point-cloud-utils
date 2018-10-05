#ifndef IMG_FILTER_H_
#define IMG_FILTER_H_

// functions that operate on pixels
// in this file code readability is preferred over speed and memory optimizations

#include <math.h>
#include <vector>
#include <algorithm>

namespace img {

template<int Channels, typename ScalarType, bool Safe>
inline void normalize(Image<Channels,ScalarType,Safe> &image)
{      
  assert(image.isValid());
  assert(image.attributes.hasRange(0,255));
  if(Safe){
    if(!image.isValid())  throw ImageException("Invalid image");
    if(!image.attributes.hasRange(0,255)) throw ImageException("Invalid range attribute");
  }
  ScalarType max=maxValue(image);
  ScalarType min=minValue(image);  
  ScalarType scale=255.0f/(max-min);

  ScalarType* array = image.dataValues();
  int length =image.dataValuesSize();

  for(int offset=0;offset<length;offset++)
    array[offset] = (array[offset]-min)*scale;

}

template<int Channels, typename ScalarType, bool Safe>
inline Image<Channels,ScalarType,Safe> getNormalized(const Image<Channels,ScalarType,Safe> &image)
{
  Image<Channels,ScalarType,Safe> i(image);
  normalize(i);
  return i;
}

template<int Channels, typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void convolution(const Image<Channels,SrcScalarType,SrcSafe> &source,Image<Channels,DestScalarType,DestSafe> &destination,const DestScalarType *matrix,int matrix_width,int matrix_height)
{ 
  assert(source.isValid());
  assert(matrix != NULL);
  assert((matrix_width > 0) && ((matrix_width)%2 == 1));
  assert((matrix_height > 0) && ((matrix_height)%2 == 1));
  if(SrcSafe || DestSafe){
    if(!source.isValid()) throw ImageException("Invalid image");
    if( matrix == NULL) throw ImageException("NULL convolution matrix");
    if( !((matrix_width > 0) && ((matrix_width)%2 == 1)) ) throw ImageException("Matrix width must be a positive odd number");
    if( !((matrix_height > 0) && ((matrix_height)%2 == 1)) ) throw ImageException("Matrix height must be a positive odd number");
  }
  destination.setZero(source.width(),source.height()); 
  destination.attributes=source.attributes;
  
  int x_radius=(matrix_width-1)/2;
  int y_radius=(matrix_height-1)/2;

  // per ogni canale
  for(int channel=0; channel < Channels ;channel++){ // canali immagine

    // per ogni riga
    for (int y = 0; y < source.height(); y++){ // righe immagine

      // per ogni pixel nella riga
      for (int x = 0; x < source.width(); x++){ // colonne immagine
      DestScalarType sum=0.0f;
      int offset=0;
      
      for(int my = y-y_radius; my <= y+y_radius; my++) // righe matrice
        for(int mx = x-x_radius; mx <= x+x_radius; mx++) //colonne matrice
          sum += matrix[offset++] * static_cast<DestScalarType>(source.getValueAsClamped(mx,my,channel));
      
     destination.setValue(x,y,channel,sum);

      } // colonne immagine
    } // righe immagine
  } // canali immagine
}

template<int Channels, typename ScalarType, bool Safe> // get[filter]ed() functions constrain return type to be the same of the parameter 
inline Image<Channels,ScalarType,Safe> getConvolved(const Image<Channels,ScalarType,Safe> &image,const ScalarType *matrix,int matrix_width,int matrix_height)
{ 
  Image<Channels,ScalarType,Safe> i;
  convolution(image,i,matrix,matrix_width,matrix_height);
  return i;
}

template<int Channels, typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void boxFilter(const Image<Channels,SrcScalarType,SrcSafe> &source,Image<Channels,DestScalarType,DestSafe> &destination,int radius)
{ 
  assert(radius > 0);
  if(SrcSafe || DestSafe){
    if(radius <= 0) throw ImageException("Nonpositive radius");
  }
  int matrix_side=2*radius+1;
  int matrix_size=matrix_side*matrix_side;
  DestScalarType* matrix=new DestScalarType[matrix_size];
  DestScalarType val=1.0f/matrix_size;
  for(int i=0; i<matrix_size; i++)
    matrix[i]=val;
  convolution(source,destination,matrix,matrix_side,matrix_side);
  delete [] matrix; 
}

template<int Channels, typename ScalarType, bool Safe> // get[filter]ed() functions constrain return type to be the same of the parameter 
inline Image<Channels,ScalarType,Safe> getBoxFiltered(const Image<Channels,ScalarType,Safe> &image,int radius)
{ 
  Image<Channels,ScalarType,Safe> i;
  boxFilter(image,i,radius);
  return i;
}

template<typename ScalarType>
inline void _gaussian(const int &radius,const ScalarType &sigma,ScalarType * &matrix,int &matrix_side)
{
  matrix_side=2*radius+1;
  int matrix_size=matrix_side*matrix_side;
  matrix=new ScalarType[matrix_size];

  // calcolo la matrice
  int offset=0;
  const ScalarType c = -0.5f / (sigma*sigma);

  for(int y = -radius;y <= radius; y++)
    for(int x = -radius;x <= radius; x++) 
      matrix[offset++] = exp( (x*x+y*y)*c );

  // porto la matrice a somma 1
  ScalarType sum=0.0f;
  for(int i=0;i<matrix_size;i++)
    sum+=matrix[i];
  for(int i=0;i<matrix_size;i++)
    matrix[i]/=sum;
}

template<int Channels, typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void GaussianSmooth(const Image<Channels,SrcScalarType,SrcSafe> &source,Image<Channels,DestScalarType,DestSafe> &destination,const int radius)
{
  assert(radius > 0.0f);
  if(SrcSafe || DestSafe){
    if(radius <= 0.0f) throw ImageException("Nonpositive radius");
  }
  const DestScalarType sigma = radius/3.0f;
  DestScalarType *matrix=NULL;
  int matrix_side=0;
  _gaussian(radius,sigma,matrix,matrix_side);

  convolution(source,destination,matrix,matrix_side,matrix_side);
  delete [] matrix; 
}

template<int Channels, typename ScalarType, bool Safe> // get[filter]ed() functions constrain return type to be the same of the parameter 
inline Image<Channels,ScalarType,Safe> getGaussianSmoothed(const Image<Channels,ScalarType,Safe> &image,const int radius)
{ 
  Image<Channels,ScalarType,Safe> i;
  GaussianSmooth(image,i,radius);
  return i;
}

template<int Channels, typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void LaplacianFilter(const Image<Channels,SrcScalarType,SrcSafe> &source,Image<Channels,DestScalarType,DestSafe> &destination)
{
  DestScalarType *matrix=new DestScalarType[9];
  matrix[0]= 0.0f; matrix[1]= 1.0f; matrix[2]= 0.0f;
  matrix[3]= 1.0f; matrix[4]=-4.0f; matrix[5]= 1.0f;
  matrix[6]= 0.0f, matrix[7]= 1.0f; matrix[8]= 0.0f;

  convolution(source,destination,matrix,3,3);
  delete [] matrix; 
}

template<int Channels, typename ScalarType, bool Safe> // get[filter]ed() functions constrain return type to be the same of the parameter 
inline Image<Channels,ScalarType,Safe> getLaplacianFiltered(const Image<Channels,ScalarType,Safe> &image)
{ 
  Image<Channels,ScalarType,Safe> i;
  LaplacianFilter(image,i);
  return i;
}

template<typename ScalarType>
inline void _laplacian_of_gaussian(const int &radius,const ScalarType &sigma,ScalarType * &matrix,int &matrix_side)
{
  matrix_side=2*radius+1;
  int matrix_size=matrix_side*matrix_side;
  matrix=new ScalarType[matrix_size];
  
  // calcolo la matrice
  int offset=0;
    
  const ScalarType c1 = -0.5f/(sigma*sigma);
  ScalarType c;
  
  for(int y = -radius;y <= radius; y++)
    for(int x = -radius;x <= radius; x++){
      c = (x*x+y*y) * c1;
      matrix[offset++] = (1+c) * exp(c);
    }
      
  // porto la matrice a somma 1
  ScalarType sum=0.0f;
  for(int i=0;i<matrix_size;i++)
    sum+=matrix[i];
  for(int i=0;i<matrix_size;i++)
    matrix[i]/=sum;
}

template<int Channels, typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void LoGFilter(const Image<Channels,SrcScalarType,SrcSafe> &source,Image<Channels,DestScalarType,DestSafe> &destination,int radius)
{
  assert(radius > 0.0f);
  if(SrcSafe || DestSafe){
    if(radius <= 0.0f) throw ImageException("Nonpositive radius");
  }
  DestScalarType *matrix=NULL;
  int matrix_side=0;
  DestScalarType sigma = radius/3.0f;
  _laplacian_of_gaussian(radius,sigma,matrix,matrix_side);
  
  convolution(source,destination,matrix,matrix_side,matrix_side);
  delete [] matrix; 
}

template<int Channels, typename ScalarType, bool Safe> // get[filter]ed() functions constrain return type to be the same of the parameter 
inline Image<Channels,ScalarType,Safe> getLoGFiltered(const Image<Channels,ScalarType,Safe> &image,int radius)
{ 
  Image<Channels,ScalarType,Safe> i;
  LoGFilter(image,i,radius);
  return i;
}

template<int Channels, typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void DoGFilter(const Image<Channels,SrcScalarType,SrcSafe> &source,Image<Channels,DestScalarType,DestSafe> &destination,int radius1,int radius2)
{
  assert(radius1 > 0.0f);
  assert(radius2 > 0.0f);
  assert(radius2 > radius1);
  if(SrcSafe || DestSafe){
    if(radius1 <= 0.0f) throw ImageException("Nonpositive radius1");
    if(radius2 <= 0.0f) throw ImageException("Nonpositive radius2");
    if(radius2 <= radius1) throw ImageException("radius2 is less than radius1");
  }
  
  int matrix_side=0;

  DestScalarType *m1=NULL,*m2=NULL;
  const DestScalarType sigma1=radius1/3.0f;
  const DestScalarType sigma2=radius2/3.0f;
  
  // ottengo la prima gaussiana (col radius della seconda)
  _gaussian(radius2,sigma1,m1,matrix_side);
  // ottengo la seconda gaussiana
  _gaussian(radius2,sigma2,m2,matrix_side);
  
  int matrix_size=matrix_side*matrix_side;
  
  
  // sottraggo alla prima gaussiana la seconda
  for(int i=0;i<matrix_size;i++)
      m1[i] -= m2[i];

  // riporto la matrice a somma 1
  DestScalarType sum=0.0f;
  for(int i=0;i<matrix_size;i++)
    sum+=m1[i];
  for(int i=0;i<matrix_size;i++)
    m1[i]/=sum;
 
  convolution(source,destination,m1,matrix_side,matrix_side);
  delete [] m1; 
  delete [] m2; 
}

template<int Channels, typename ScalarType, bool Safe> // get[filter]ed() functions constrain return type to be the same of the parameter 
inline Image<Channels,ScalarType,Safe> getDoGFiltered(const Image<Channels,ScalarType,Safe> &image,int radius1,int radius2)
{ 
  Image<Channels,ScalarType,Safe> i;
  DoGFilter(image,i,radius1,radius2);
  return i;
}

template<int Channels, typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void UnsharpMask(const Image<Channels,SrcScalarType,SrcSafe> &source,Image<Channels,DestScalarType,DestSafe> &destination,int radius,float factor)
{
  assert(radius > 0);
  assert(factor > 0.0f); 
  if(SrcSafe || DestSafe){
    if(radius <= 0.0f) throw ImageException("Nonpositive radius");
    if(factor <= 0.0f) throw ImageException("Nonpositive factor");
  }
  // metto l'immagine smoothata in destination
  GaussianSmooth(source,destination,radius);
  
  DestScalarType* source_array = source.dataValues();
  DestScalarType* destination_array = destination.dataValues();
  int length = source.dataValuesSize();
 
  // unsharpo destination in loco
  for(int offset=0;offset<length;offset++)
      destination_array[offset] = source_array[offset]+factor*(source_array[offset]-destination_array[offset]);

}

template<int Channels, typename ScalarType, bool Safe> // get[filter]ed() functions constrain return type to be the same of the parameter 
inline Image<Channels,ScalarType,Safe> getUnsharpMasked(const Image<Channels,ScalarType,Safe> &image,int radius,float factor)
{ 
  Image<Channels,ScalarType,Safe> i;
  UnsharpMask(image,i,radius,factor);
  return i;
}

template<int Channels, typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void medianFilter(const Image<Channels,SrcScalarType,SrcSafe> &source,Image<Channels,DestScalarType,DestSafe> &destination,int radius)
{ 
  assert(source.isValid());
  assert(radius > 0);
  if(SrcSafe || DestSafe){
    if(!source.isValid()) throw ImageException("Invalid image");
    if(radius <= 0) throw ImageException("Nonpositive radius");
  }
  destination.setZero(source.width(),source.height()); 
  destination.attributes=source.attributes;
 
  // per ogni canale
  for(int channel=0; channel < Channels;channel++){ // canali immagine

    // per ogni riga
    for (int y = 0; y < source.height(); y++){ // righe immagine

      // per ogni pixel nella riga
      for (int x = 0; x < source.width(); x++){ // colonne immagine
        
        // memorizzo i valori dell'intorno
      std::vector<DestScalarType> v;
      for(int my = y-radius; my <= y+radius; my++) // righe intorno
        for(int mx = x-radius; mx <= x+radius; mx++) //colonne intorno
          if (source.isInside(mx,my))
          v.push_back(static_cast<DestScalarType>(source.getValue(mx,my,channel)));
         
      // ottengo la mediana 
      int s=v.size();
      assert(s>0);
      nth_element (v.begin(), v.begin()+(s/2), v.end());
      DestScalarType median = *(v.begin()+(s/2));
      if((s%2)==0) { // even s: mean of the 2 middle elements
      nth_element (v.begin(), v.begin()+(s/2)+1, v.end());
        median = ( *(v.begin()+(s/2)+1) + median ) /2.0f;
      }      
     // aggiorno il valore alla mediana dell'intorno
     destination.setValue(x,y,channel,median);

      } // colonne immagine
    } // righe immagine
  } // canali immagine
}

template<int Channels, typename ScalarType, bool Safe> // get[filter]ed() functions constrain return type to be the same of the parameter 
inline Image<Channels,ScalarType,Safe> getMedianFiltered(const Image<Channels,ScalarType,Safe> &image,int radius)
{ 
  Image<Channels,ScalarType,Safe> i;
  medianFilter(image,i,radius);
  return i;
}

template<int Channels, typename ScalarType1, bool Safe1,typename ScalarType2, bool Safe2> 
inline void channels_mean(const img::Image<Channels,ScalarType1,Safe1> &channels_image, img::Image<1,ScalarType2,Safe2> &mean_image)
{
  assert(channels_image.isValid());
  if(Safe1 || Safe2){
    if(!channels_image.isValid())  throw img::ImageException("channels_image rgb image");
  }
  mean_image.setZero(channels_image.width(),channels_image.height());
  mean_image.attributes=channels_image.attributes;
  
  for (int y_coord = 0; y_coord < channels_image.height(); ++y_coord)
    for (int x_coord = 0; x_coord < channels_image.width(); ++x_coord){
      ScalarType2 sum = ScalarType2(0.0);
      for (int channel=0; channel<Channels; ++channel)
        sum += static_cast<ScalarType2>(channels_image.getValue(x_coord, y_coord, channel));
      mean_image.setValue(x_coord, y_coord, 0, sum / ScalarType2(Channels) );
    } 
}

} //end namespace img
 
#endif /*IMG_FILTER_H_*/
