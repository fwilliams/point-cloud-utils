#ifndef IMG_CS_BASE_H_
#define IMG_CS_BASE_H_

// warning: temporary name/location, will change in near future

#include "img/img.h"

namespace img {

// RGB: GAMMA -> LINEAR 
template<int Channels, typename ScalarType1, bool Safe1,typename ScalarType2, bool Safe2> 
inline void convert_gamma_precompensated_rgb_to_linear_rgb(const img::Image<Channels,ScalarType1,Safe1> &gamma_precompensated_rgb_image, img::Image<Channels,ScalarType2,Safe2> &linear_rgb_image)
{
  assert(gamma_precompensated_rgb_image.isValid());
  assert(!gamma_precompensated_rgb_image.attributes.hasColorspace(SRGB));
  assert(gamma_precompensated_rgb_image.attributes.hasRange(0.0,1.0));
  ScalarType2 old_gamma;
  gamma_precompensated_rgb_image.attributes.getGamma(old_gamma);
  assert((old_gamma>0.0)&&(old_gamma<1.0));
  if(Safe1 || Safe2){
    if(!gamma_precompensated_rgb_image.isValid())  throw ImageException("Invalid rgb image");
    if(gamma_precompensated_rgb_image.attributes.hasColorspace(SRGB)) throw ImageException("Invalid colorspace attribute");
    if(!gamma_precompensated_rgb_image.attributes.hasRange(0,1)) throw ImageException("Invalid range attribute");
    if(!((old_gamma>0.0)&&(old_gamma<1.0))) throw ImageException("Invalid gamma attribute");
  }
  linear_rgb_image.setZero(gamma_precompensated_rgb_image.width(),gamma_precompensated_rgb_image.height());
  linear_rgb_image.attributes=gamma_precompensated_rgb_image.attributes;
  
  for(int channel=0; channel<Channels; ++channel)
    for(int x_coord=0; x_coord<gamma_precompensated_rgb_image.width(); ++x_coord)
      for(int y_coord=0; y_coord<gamma_precompensated_rgb_image.height(); ++y_coord){
        linear_rgb_image.setValue(x_coord,y_coord,channel, pow(ScalarType2(gamma_precompensated_rgb_image.getValue(x_coord,y_coord,channel)),ScalarType2(1.0)/old_gamma));
      }
  linear_rgb_image.attributes.setGamma(ScalarType2(1.0));
}

// RGB: LINEAR -> GAMMA
// when converting from linear to gamma precompensated the gamma parameter mus be in range (0.0,1.0)
template<int Channels, typename ScalarType1, bool Safe1,typename ScalarType2, bool Safe2> 
inline void convert_linear_rgb_to_gamma_precompensated_rgb(const img::Image<Channels,ScalarType1,Safe1> &linear_rgb_image, img::Image<Channels,ScalarType2,Safe2> &gamma_precompensated_rgb_image, ScalarType2 gamma=ScalarType2(1.0/2.2))
{
  assert(linear_rgb_image.isValid());
  assert(!linear_rgb_image.attributes.hasColorspace(SRGB));
  assert(linear_rgb_image.attributes.hasRange(0.0,1.0));
  assert(linear_rgb_image.attributes.hasGamma(1.0));
  assert((gamma>0.0)&&(gamma<1.0));
  if(Safe1 || Safe2){
    if(!linear_rgb_image.isValid())  throw img::ImageException("Invalid rgb image");
    if(linear_rgb_image.attributes.hasColorspace(SRGB)) throw ImageException("Invalid colorspace attribute");
    if(!linear_rgb_image.attributes.hasRange(0.0,1.0)) throw ImageException("Invalid range attribute");
    if(!linear_rgb_image.attributes.hasGamma(1.0)) throw ImageException("Invalid gamma attribute");
    if(!((gamma>0.0)&&(gamma<1.0))) throw ImageException("Invalid gamma parameter");
  }
  gamma_precompensated_rgb_image.setZero(linear_rgb_image.width(),linear_rgb_image.height());
  gamma_precompensated_rgb_image.attributes=linear_rgb_image.attributes;
  
  for(int channel=0; channel<Channels; ++channel)
    for(int x_coord=0; x_coord<linear_rgb_image.width(); ++x_coord)
      for(int y_coord=0; y_coord<linear_rgb_image.height(); ++y_coord){
        gamma_precompensated_rgb_image.setValue(x_coord,y_coord,channel, pow(ScalarType2(linear_rgb_image.getValue(x_coord,y_coord,channel)),gamma));
      }
  gamma_precompensated_rgb_image.attributes.setGamma(gamma);
}

// SRGB: GAMMA -> LINEAR 
template<int Channels, typename ScalarType1, bool Safe1,typename ScalarType2, bool Safe2> 
inline void convert_gamma_precompensated_srgb_to_linear_srgb(const img::Image<Channels,ScalarType1,Safe1> &gamma_precompensated_srgb_image, img::Image<Channels,ScalarType2,Safe2> &linear_srgb_image)
{
  assert(gamma_precompensated_srgb_image.isValid());
  assert(gamma_precompensated_srgb_image.attributes.hasColorspace(SRGB));
  assert(gamma_precompensated_srgb_image.attributes.hasRange(0.0,1.0));
  assert(!gamma_precompensated_srgb_image.attributes.hasGamma(1.0));
  if(Safe1 || Safe2){
    if(!gamma_precompensated_srgb_image.isValid())  throw ImageException("Invalid rgb image");
    if(!gamma_precompensated_srgb_image.attributes.hasColorspace(SRGB)) throw ImageException("Invalid colorspace attribute");
    if(!gamma_precompensated_srgb_image.attributes.hasRange(0,1)) throw ImageException("Invalid range attribute");
    if(gamma_precompensated_srgb_image.attributes.hasGamma(1.0)) throw ImageException("Invalid gamma attribute");
  }
  linear_srgb_image.setZero(gamma_precompensated_srgb_image.width(),gamma_precompensated_srgb_image.height());
  linear_srgb_image.attributes=gamma_precompensated_srgb_image.attributes;
  
  for(int channel=0; channel<Channels; ++channel)
    for(int x_coord=0; x_coord<gamma_precompensated_srgb_image.width(); ++x_coord)
      for(int y_coord=0; y_coord<gamma_precompensated_srgb_image.height(); ++y_coord){
        ScalarType2 value = ScalarType2(gamma_precompensated_srgb_image.getValue(x_coord,y_coord,channel));
        if(value<=ScalarType2(0.04045))
          value = value / ScalarType2(12.92);
        else
          value = pow( (value+ScalarType2(0.055)) / ScalarType2(1.055), ScalarType2(2.4) );
        linear_srgb_image.setValue(x_coord, y_coord, channel, value);
      }
  linear_srgb_image.attributes.setGamma(ScalarType2(1.0));
}

// SRGB: LINEAR -> GAMMA
template<int Channels, typename ScalarType1, bool Safe1,typename ScalarType2, bool Safe2> 
inline void convert_linear_srgb_to_gamma_precompensated_srgb(const img::Image<Channels,ScalarType1,Safe1> &linear_srgb_image, img::Image<Channels,ScalarType2,Safe2> &gamma_precompensated_srgb_image)
{
  assert(linear_srgb_image.isValid());
  assert(linear_srgb_image.attributes.hasColorspace(SRGB));
  assert(linear_srgb_image.attributes.hasRange(0.0,1.0));
  assert(linear_srgb_image.attributes.hasGamma(1.0));
  if(Safe1 || Safe2){
    if(!linear_srgb_image.isValid())  throw img::ImageException("Invalid rgb image");
    if(!linear_srgb_image.attributes.hasColorspace(SRGB)) throw ImageException("Invalid colorspace attribute");
    if(!linear_srgb_image.attributes.hasRange(0.0,1.0)) throw ImageException("Invalid range attribute");
    if(!linear_srgb_image.attributes.hasGamma(1.0)) throw ImageException("Invalid gamma attribute");
  }
  gamma_precompensated_srgb_image.setZero(linear_srgb_image.width(),linear_srgb_image.height());
  gamma_precompensated_srgb_image.attributes=linear_srgb_image.attributes;
  
  for(int channel=0; channel<Channels; ++channel)
    for(int x_coord=0; x_coord<linear_srgb_image.width(); ++x_coord)
      for(int y_coord=0; y_coord<linear_srgb_image.height(); ++y_coord){
        ScalarType2 value = ScalarType2(linear_srgb_image.getValue(x_coord,y_coord,channel));
        if(value<=ScalarType2(0.0031308))
          value =  value * ScalarType2(12.92);
        else
          value = (ScalarType2(1.055) * pow(value, ScalarType2(1.0/2.4))) - ScalarType2(0.055);
        gamma_precompensated_srgb_image.setValue(x_coord, y_coord, channel, value);
      }
  gamma_precompensated_srgb_image.attributes.setGamma(ScalarType2(2.2)); // approssimation
}

// SRGB -> XYZ
template<typename ScalarType1, bool Safe1,typename ScalarType2, bool Safe2> 
inline void convert_srgb_to_xyz(const img::Image<3,ScalarType1,Safe1> &rgb_image, img::Image<3,ScalarType2,Safe2> &xyz_image)
{
  assert( rgb_image.isValid());
  assert( rgb_image.attributes.hasRange(0.0,1.0));
  assert( rgb_image.attributes.hasGamma(1.0));
  assert( rgb_image.attributes.hasColorspace(img::RGB) );
  if(Safe1 || Safe2){
    if(!rgb_image.isValid())  throw img::ImageException("Invalid rgb image");
    if(!rgb_image.attributes.hasRange(0.0,1.0))  throw img::ImageException("Invalid range attribute");
    if(!rgb_image.attributes.hasGamma(1.0))  throw img::ImageException("Invalid gamma attribute");
    if(!rgb_image.attributes.hasColorspace(img::RGB))  throw img::ImageException("Invalid colorspace attribute");
  }
  xyz_image.setZero(rgb_image.width(),rgb_image.height());
  xyz_image.attributes=rgb_image.attributes;

  for(int x_coord=0; x_coord<rgb_image.width(); ++x_coord)
    for(int y_coord=0; y_coord<rgb_image.height(); ++y_coord){
      const ScalarType2 R = static_cast<ScalarType2>(rgb_image.getValue(x_coord,y_coord,0));
      const ScalarType2 G = static_cast<ScalarType2>(rgb_image.getValue(x_coord,y_coord,1));
      const ScalarType2 B = static_cast<ScalarType2>(rgb_image.getValue(x_coord,y_coord,2));

      xyz_image.setValue(x_coord,y_coord,0, ( 0.412453f*R + 0.357580f*G + 0.180423f*B ) );
      xyz_image.setValue(x_coord,y_coord,1, ( 0.212671f*R + 0.715160f*G + 0.072169f*B ) );
      xyz_image.setValue(x_coord,y_coord,2, ( 0.019334f*R + 0.119193f*G + 0.950227f*B ) );
    }
    xyz_image.setColorspace(img::CIE_XYZ);
}

// XYZ -> SRGB
template<typename ScalarType1, bool Safe1,typename ScalarType2, bool Safe2> 
inline void convert_xyz_to_rgb(const img::Image<3,ScalarType1,Safe1> &xyz_image, img::Image<3,ScalarType2,Safe2> &rgb_image)
{
  assert( xyz_image.isValid());
  assert( xyz_image.attributes.hasRange(0.0,1.0));
  assert( xyz_image.attributes.hasGamma(1.0));
  assert( xyz_image.attributes.hasColorspace(img::CIE_XYZ) );
  
  if(Safe1 || Safe2){
    if(!xyz_image.isValid())  throw img::ImageException("Invalid xyz image");
    if(!xyz_image.attributes.hasRange(0.0,1.0))  throw img::ImageException("Invalid range attribute");
    if(!xyz_image.attributes.hasGamma(1.0))  throw img::ImageException("Invalid gamma attribute");
    if(!rgb_image.attributes.hasColorspace(img::CIE_XYZ))  throw img::ImageException("Invalid colorspace attribute");
  }
  rgb_image.setZero(xyz_image.width(),xyz_image.height());
  rgb_image.attributes=xyz_image.attributes;

  for(int x=0; x<xyz_image.width(); ++x)
    for(int y=0; y<xyz_image.height(); ++y){
      const ScalarType2 X = static_cast<ScalarType2>(xyz_image.getValue(x,y,0));
      const ScalarType2 Y = static_cast<ScalarType2>(xyz_image.getValue(x,y,1));
      const ScalarType2 Z = static_cast<ScalarType2>(xyz_image.getValue(x,y,2));
      
      const ScalarType2 R = 3.240479f*X + -1.537150f*Y + -0.498535f*Z;
      const ScalarType2 G = -0.969256f*X + 1.875992f*Y + 0.041556f*Z;
      const ScalarType2 B = 0.055648f*X + -0.204043f*Y + 1.057311f*Z;
      
      rgb_image.setValue(x,y,0, ( R<0.0f?0.0f:(R>1.0f?1.0f:R) ) );
      rgb_image.setValue(x,y,1, ( G<0.0f?0.0f:(G>1.0f?1.0f:G) ) );
      rgb_image.setValue(x,y,2, ( B<0.0f?0.0f:(B>1.0f?1.0f:B) ) );
    }
    rgb_image.setColorspace(img::RGB);
}

///ora nn ho tempo:
//// AltriRGB -> XYZ
//// XYZ -> AltriRGB
//// XYZ -> LAB
//// LAB -> XYZ
//template<typename ScalarType1, bool Safe1,typename ScalarType2, bool Safe2> 
//inline void convert_xyz_to_lab(const img::Image<3,ScalarType1,Safe1> &xyz_image, img::Image<3,ScalarType2,Safe2> &lab_image)
//{
//  assert(xyz_image.isValid());
//  if(Safe1 || Safe2){
//    if(!xyz_image.isValid())  throw img::ImageException("Invalid xyz image");
//  }
//  lab_image.setZero(xyz_image.width(),xyz_image.height());
//  
//  const ScalarType2 one_third=0.33333333333333333333333333333333333333333333333333333333333333f;
//  const ScalarType2 Xn = 0.9513f;
//  const ScalarType2 Yn = 1.000f;
//  const ScalarType2 Zn = 1.0886f;
//  
//  for(int x=0; x<xyz_image.width(); ++x)
//    for(int y=0; y<xyz_image.height(); ++y){
//      const ScalarType2 X_third = pow(static_cast<ScalarType2>(xyz_image.getValue(x,y,0))/Xn,one_third);
//      const ScalarType2 Y = static_cast<ScalarType2>(xyz_image.getValue(x,y,1))/Yn;
//      const ScalarType2 Y_third = pow(Y,one_third);
//      const ScalarType2 Z_third = pow(static_cast<ScalarType2>(xyz_image.getValue(x,y,2))/Zn,one_third);
//
//      lab_image.setValue(x,y,0, (Y > 0.008856f)?((116.0f*(Y_third)) - 16.0f):(903.3f*Y) );
//      lab_image.setValue(x,y,1, 500.0f * ((X_third) - (Y_third)) );
//      lab_image.setValue(x,y,2, 200.0f * ((Y_third) - (Z_third)) );
//    }
//}
//
//template<typename ScalarType1, bool Safe1,typename ScalarType2, bool Safe2> 
//inline void convert_lab_to_xyz(const img::Image<3,ScalarType1,Safe1> &lab_image, img::Image<3,ScalarType2,Safe2> &xyz_image)
//{
//  assert(0); // sistemare attributi
//  assert(lab_image.isValid());
//  if(Safe1 || Safe2){
//    if(!lab_image.isValid())  throw img::ImageException("Invalid lab image");
//  }
//  xyz_image.setZero(lab_image.width(),lab_image.height());
//
//  const ScalarType2 Xn = 0.9513f;
//  const ScalarType2 Yn = 1.000f;
//  const ScalarType2 Zn = 1.0886f;
//
//  for(int x=0; x<lab_image.width(); ++x)
//    for(int y=0; y<lab_image.height(); ++y){
//      const ScalarType2 P = (static_cast<ScalarType2>(lab_image.getValue(x,y,0))+16.0f)/116.0f;
//
//      xyz_image.setValue(x,y,0, Xn * pow(P + (lab_image.getValue(x,y,1)/500.0f), 3.0f) );
//      xyz_image.setValue(x,y,1, Yn * pow(P, 3.0f) );
//      xyz_image.setValue(x,y,2, Zn * pow(P - (lab_image.getValue(x,y,2)/200.0f), 3.0f) );
//    }
//}

} // end namespace img

#endif /*IMG_CS_BASE_H_*/
