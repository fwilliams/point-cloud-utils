#ifndef IMG_CONVERT_H_
#define IMG_CONVERT_H_

// implementation of conversions between basic image types
// assumes linearity for both source and destination
// uses CIE 1931 luminance definition for default grayscale conversion:
// Y = 0.212671*R + 0.715160*G + 0.072169*B
// assumes 255 for default alpha_value parameter

#include "img/img_image.h"

namespace img {

// Y to RGB
template<typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void convert_Y_to_RGB(const Image<1,SrcScalarType,SrcSafe> &source, Image<3,DestScalarType,DestSafe> &destination)
{
  assert(source.isValid());
  if(SrcSafe || DestSafe){
    if(!source.isValid())  throw ImageException("Invalid source image");
  }
  destination.setZero(source.width(),source.height());
  for (int y_coord = 0; y_coord < source.height(); ++y_coord)
    for (int x_coord = 0; x_coord < source.width(); ++x_coord){
      DestScalarType Y = static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 0));
      destination.setValue(x_coord, y_coord, 0, Y);
      destination.setValue(x_coord, y_coord, 1, Y);
      destination.setValue(x_coord, y_coord, 2, Y);
    } 
  destination.attributes.setRange(source.attributes);
  destination.attributes.setGamma(source.attributes);
}

// Y to RGBA
template<typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void convert_Y_to_RGBA(const Image<1,SrcScalarType,SrcSafe> &source, Image<4,DestScalarType,DestSafe> &destination, DestScalarType alpha_value=255.0f )
{
  assert(source.isValid());
  if(SrcSafe || DestSafe){
    if(!source.isValid())  throw ImageException("Invalid source image");
  }
  destination.setZero(source.width(),source.height());
  for (int y_coord = 0; y_coord < source.height(); ++y_coord)
    for (int x_coord = 0; x_coord < source.width(); ++x_coord){
      DestScalarType Y = static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 0));
      destination.setValue(x_coord, y_coord, 0, Y);
      destination.setValue(x_coord, y_coord, 1, Y);
      destination.setValue(x_coord, y_coord, 2, Y);
      destination.setValue(x_coord, y_coord, 3, alpha_value);
    } 
  destination.attributes.setRange(source.attributes);
  destination.attributes.setGamma(source.attributes);
}

// RGB to Y
template<typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void convert_RGB_to_Y(const Image<3,SrcScalarType,SrcSafe> &source, Image<1,DestScalarType,DestSafe> &destination)
{
  assert(source.isValid());
  if(SrcSafe || DestSafe){
    if(!source.isValid())  throw ImageException("Invalid source image");
  }
  destination.setZero(source.width(),source.height());
  for (int y_coord = 0; y_coord < source.height(); ++y_coord)
    for (int x_coord = 0; x_coord < source.width(); ++x_coord){
      DestScalarType R = static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 0));
      DestScalarType G = static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 1));
      DestScalarType B = static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 2));
      destination.setValue(x_coord, y_coord, 0, 0.212671f*R + 0.715160f*G + 0.072169f*B);
    } 
  destination.attributes.setRange(source.attributes);
  destination.attributes.setGamma(source.attributes);
}

// RGB to RGBA
template<typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void convert_RGB_to_RGBA(const Image<3,SrcScalarType,SrcSafe> &source, Image<4,DestScalarType,DestSafe> &destination, DestScalarType alpha_value=255.0f )
{
  assert(source.isValid());
  if(SrcSafe || DestSafe){
    if(!source.isValid())  throw ImageException("Invalid source image");
  }
  destination.setZero(source.width(),source.height());
  for (int y_coord = 0; y_coord < source.height(); ++y_coord)
    for (int x_coord = 0; x_coord < source.width(); ++x_coord){
      destination.setValue(x_coord, y_coord, 0, static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 0)));
      destination.setValue(x_coord, y_coord, 1, static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 1)));
      destination.setValue(x_coord, y_coord, 2, static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 2)));
      destination.setValue(x_coord, y_coord, 3, alpha_value);
    } 
  destination.attributes.setRange(source.attributes);
  destination.attributes.setGamma(source.attributes);
}

// RGBA to Y
template<typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void convert_RGBA_to_Y(const Image<4,SrcScalarType,SrcSafe> &source, Image<1,DestScalarType,DestSafe> &destination)
{
  assert(source.isValid());
  if(SrcSafe || DestSafe){
    if(!source.isValid())  throw ImageException("Invalid source image");
  }
  destination.setZero(source.width(),source.height());
  for (int y_coord = 0; y_coord < source.height(); ++y_coord)
    for (int x_coord = 0; x_coord < source.width(); ++x_coord){
      DestScalarType R = static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 0));
      DestScalarType G = static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 1));
      DestScalarType B = static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 2));
      destination.setValue(x_coord, y_coord, 0, 0.212671f*R + 0.715160f*G + 0.072169f*B);
    } 
  destination.attributes.setRange(source.attributes);
  destination.attributes.setGamma(source.attributes);
}

// RGBA to RGB
template<typename SrcScalarType, bool SrcSafe,typename DestScalarType, bool DestSafe> 
inline void convert_RGBA_to_RGB(const Image<4,SrcScalarType,SrcSafe> &source, Image<3,DestScalarType,DestSafe> &destination)
{
  assert(source.isValid());
  if(SrcSafe || DestSafe){
    if(!source.isValid())  throw ImageException("Invalid source image");
  }
  destination.setZero(source.width(),source.height());
  for (int y_coord = 0; y_coord < source.height(); ++y_coord)
    for (int x_coord = 0; x_coord < source.width(); ++x_coord){
      destination.setValue(x_coord, y_coord, 0, static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 0)));
      destination.setValue(x_coord, y_coord, 1, static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 1)));
      destination.setValue(x_coord, y_coord, 2, static_cast<DestScalarType>(source.getValue(x_coord, y_coord, 2)));
    } 
  destination.attributes.setRange(source.attributes);
  destination.attributes.setGamma(source.attributes);
}

// range conversion
template<int Channels, typename SrcScalarType, bool SrcSafe, typename DestScalarType, bool DestSafe>
inline void convert_range_0_255_to_0_1(const Image<Channels,SrcScalarType,SrcSafe> &source, Image<Channels,DestScalarType,DestSafe> &destination)
{
  assert(source.isValid());
  assert(source.attributes.hasRange(0,255));
  if(SrcSafe || DestSafe){
    if(!source.isValid())  throw ImageException("Invalid source image");
    if(!source.attributes.hasRange(0,255)) throw ImageException("Invalid range attribute");
  }
  destination.setZero(source.width(),source.height());
  for(int offset=0;offset<source.dataValuesSize();offset++)
    destination.dataValues()[offset] = static_cast<DestScalarType>(source.dataValues()[offset]) / DestScalarType(255.0);
  destination.attributes=source.attributes;
  destination.attributes.setRange(0,1);
}

template<int Channels, typename SrcScalarType, bool SrcSafe, typename DestScalarType, bool DestSafe>
inline void convert_range_0_1_to_0_255(const Image<Channels,SrcScalarType,SrcSafe> &source, Image<Channels,DestScalarType,DestSafe> &destination)
{
  assert(source.isValid());
  assert(source.attributes.hasRange(0,1));
  if(SrcSafe || DestSafe){
    if(!source.isValid())  throw ImageException("Invalid source image");
    if(!source.attributes.hasRange(0,1)) throw ImageException("Invalid range attribute");
  }
  destination.setZero(source.width(),source.height());
  for(int offset=0;offset<source.dataValuesSize();offset++)
    destination.dataValues()[offset] = static_cast<DestScalarType>(source.dataValues()[offset]) * DestScalarType(255.0);
  destination.attributes=source.attributes;
  destination.attributes.setRange(0,255);
}

} //end namespace img

#endif /*IMG_CONVERT_H_*/
