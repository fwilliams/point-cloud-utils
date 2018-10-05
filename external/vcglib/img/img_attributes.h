#ifndef IMG_ATTRIBUTES_H_
#define IMG_ATTRIBUTES_H_

#include "img/img_base.h"

namespace img {

  enum COLORSPACE {
    UNDEFINED,
    RGB,
    SRGB,
    CIE_XYZ,
    CIE_LAB,
    CIE_LUV,
  };

template <typename ScalarType=double>
class ImgAttributes
{
public:

private:
  COLORSPACE _colorspace;
  
  ScalarType _gamma;
  
  ScalarType _range_min;
  ScalarType _range_max;
  
  ScalarType _reference_white_x;
  ScalarType _reference_white_y;
  ScalarType _reference_white_z;

public:
  ImgAttributes()
  : _colorspace(UNDEFINED), _gamma(ScalarType(0.0)),
    _range_min(ScalarType(0.0)), _range_max(ScalarType(0.0)),
    _reference_white_x(ScalarType(0.0)), _reference_white_y(ScalarType(0.0)), _reference_white_z(ScalarType(0.0))
  {
    STATIC_FLOAT_OR_DOUBLE_TYPECHECK( ScalarType );
  }
  
  ImgAttributes(const ImgAttributes<ScalarType>  &attributes)
  {
    STATIC_FLOAT_OR_DOUBLE_TYPECHECK( ScalarType );
    _colorspace = attributes._colorspace;
    _gamma = attributes._gamma;
    _range_min = attributes._range_min;
    _range_max = attributes._range_max;
    _reference_white_x = attributes._reference_white_x;
    _reference_white_y = attributes._reference_white_y;
    _reference_white_z = attributes._reference_white_z;
  }
  
  template<typename OtherScalarType>
  ImgAttributes(const ImgAttributes<OtherScalarType>  &attributes)
  {
    STATIC_FLOAT_OR_DOUBLE_TYPECHECK( ScalarType );
    _colorspace = attributes._colorspace;
    _gamma = static_cast<ScalarType>(attributes._gamma);
    _range_min = static_cast<ScalarType>(attributes._range_min);
    _range_max = static_cast<ScalarType>(attributes._range_max);
    _reference_white_x = static_cast<ScalarType>(attributes._reference_white_x);
    _reference_white_y = static_cast<ScalarType>(attributes._reference_white_y);
    _reference_white_z = static_cast<ScalarType>(attributes._reference_white_z);
  }
  
  template<typename OtherScalarType> 
  inline ImgAttributes< ScalarType> & operator =(const ImgAttributes<OtherScalarType> &attributes)
  {
    _colorspace = attributes._colorspace;
    _gamma = static_cast<ScalarType>(attributes._gamma);
    _range_min = static_cast<ScalarType>(attributes._range_min);
    _range_max = static_cast<ScalarType>(attributes._range_max);
    _reference_white_x = static_cast<ScalarType>(attributes._reference_white_x);
    _reference_white_y = static_cast<ScalarType>(attributes._reference_white_y);
    _reference_white_z = static_cast<ScalarType>(attributes._reference_white_z);
  }
  
  inline void reset()
  {
    _colorspace = UNDEFINED;
    _gamma = ScalarType(0.0);
    _range_min = ScalarType(0.0);
    _range_max = ScalarType(0.0);
    _reference_white_x = ScalarType(0.0);
    _reference_white_y = ScalarType(0.0);
    _reference_white_z = ScalarType(0.0);
  }

  // getters
  void getColorspace(COLORSPACE &ret_colorspace) const
  {
    ret_colorspace = _colorspace;
  }

  template<typename OtherScalarType> 
  void getGamma(OtherScalarType &ret_gamma) const
  {
    STATIC_FLOAT_OR_DOUBLE_TYPECHECK( OtherScalarType );
    ret_gamma = static_cast<OtherScalarType>(_gamma);
  }

  template<typename OtherScalarType> 
  void getRange(OtherScalarType &ret_range_min, OtherScalarType &ret_range_max) const
  {
    STATIC_FLOAT_OR_DOUBLE_TYPECHECK( OtherScalarType );
    ret_range_min = static_cast<OtherScalarType>(_range_min);
    ret_range_max = static_cast<OtherScalarType>(_range_max);
  }

  template<typename OtherScalarType> 
  void getReferenceWhite(OtherScalarType &ret_reference_white_x, OtherScalarType &ret_reference_white_y, OtherScalarType &ret_reference_white_z) const
  {
    STATIC_FLOAT_OR_DOUBLE_TYPECHECK( OtherScalarType );
    ret_reference_white_x = static_cast<OtherScalarType>(_reference_white_x);
    ret_reference_white_y = static_cast<OtherScalarType>(_reference_white_y);
    ret_reference_white_z = static_cast<OtherScalarType>(_reference_white_z);
  }

  // setters
  void setColorspace(COLORSPACE arg_colorspace)
  {
    _colorspace = arg_colorspace;
  }
  
  void setGamma(ScalarType arg_gamma)
  {
    _gamma = arg_gamma;
  }
  
  void setRange(ScalarType arg_range_min, ScalarType arg_range_max)
  {
    assert(arg_range_min<=arg_range_max);
    _range_min = arg_range_min;
    _range_max = arg_range_max;
  }
  
  void setReferenceWhite(ScalarType arg_reference_white_x, ScalarType arg_reference_white_y, ScalarType arg_reference_white_z)
  {
    _reference_white_x = arg_reference_white_x;
    _reference_white_y = arg_reference_white_y;
    _reference_white_z = arg_reference_white_z;
  }

  void setColorspace(const ImgAttributes<ScalarType> &attributes)
  {
    _colorspace = attributes._colorspace;
  }
  
  void setGamma(const ImgAttributes<ScalarType> &attributes)
  {
    _gamma = attributes._gamma;
  }
  
  void setRange(const ImgAttributes<ScalarType> &attributes)
  {
    _range_min = attributes._range_min;
    _range_max = attributes._range_max;
  }
  
  void setReferenceWhite(const ImgAttributes<ScalarType> &attributes)
  {
    _reference_white_x = attributes._reference_white_x;
    _reference_white_y = attributes._reference_white_y;
    _reference_white_z = attributes._reference_white_z;
  }

  // checks
  bool hasColorspace(COLORSPACE arg_colorspace) const
  {
    return _colorspace == arg_colorspace;
  }
  
  bool hasGamma(ScalarType arg_gamma) const
  {
    return _gamma == arg_gamma;
  }
  
  
  bool hasRange(ScalarType arg_range_min, ScalarType arg_range_max) const
  {
    return (_range_min == arg_range_min) &&
           (_range_max == arg_range_max);
  }
  
  bool hasReferenceWhite(ScalarType arg_reference_white_x, ScalarType arg_reference_white_y, ScalarType arg_reference_white_z) const
  {
    return (_reference_white_x == arg_reference_white_x) &&
           (_reference_white_y == arg_reference_white_y) &&
           (_reference_white_z == arg_reference_white_z);
  }
  
  bool hasColorspace(const ImgAttributes<ScalarType> &attributes) const
  {
    return _colorspace == attributes._colorspace;
  }
  
  bool hasGamma(const ImgAttributes<ScalarType> &attributes) const
  {
    return _gamma == attributes._gamma;
  }
  
  bool hasRange(const ImgAttributes<ScalarType> &attributes) const
  {
    return (_range_min == attributes._range_min) &&
           (_range_max == attributes._range_max);
  }
  
  bool hasReferenceWhite(const ImgAttributes<ScalarType> &attributes) const
  {
    return (_reference_white_x == attributes._reference_white_x) &&
           (_reference_white_y == attributes._reference_white_y) &&
           (_reference_white_z == attributes._reference_white_z);
  }
  
  
  bool operator==(const ImgAttributes<ScalarType> &attributes) const
  {
    return hasColorspace(attributes) && 
           hasGamma(attributes) &&
           hasRange(attributes) &&
           hasReferenceWhite(attributes);
  }

  static ImgAttributes createImgAttributes(COLORSPACE arg_colorspace = ScalarType(UNDEFINED),
                                   ScalarType arg_gamma = ScalarType(0.0),
                                   ScalarType arg_range_min = ScalarType(0.0),
                                   ScalarType arg_range_max = ScalarType(0.0),
                                   ScalarType arg_reference_white_x = ScalarType(0.0),
                                   ScalarType arg_reference_white_y = ScalarType(0.0),
                                   ScalarType arg_reference_white_z = ScalarType(0.0))
  {
    ImgAttributes attributes;
    attributes.setColorspace(arg_colorspace);
    attributes.setGamma(arg_gamma);
    attributes.setRange(arg_range_min, arg_range_max);
    attributes.setReferenceWhite(arg_reference_white_x, arg_reference_white_y, arg_reference_white_z);
    return attributes;
  }
  
};

} // end namespace img

#endif /*IMG_ATTRIBUTES_H_*/
