#ifndef IMG_IMAGE_H_
#define IMG_IMAGE_H_
/*! \file img_image.h
    \brief definition of the generic image class

    This header contains the image class definition.
*/

#include "img/img_base.h"
#include "img/img_attributes.h"

namespace img {

/*! \brief the generic image class

  The image class is templated over three parameters:
   - The number of the image channels: an image can have from one to an arbitrary large number of channels, i.e. values per pixel. The default  is to have 3 color channels.
   - The type of the values: the image stores the pixels' values in floating-point variables. The default is to have double floating point precision, the other option is to have single floating point precision. Static assertions are used to ensure that the value type of an image is a floating  point number.
   - The safeness: the image throws runtime exceptions when its member functions are called with "wrong" parameters. This behavior, active by default, can be disabled in order to speed up the computation. The parameter correctness is also independently checked with the dynamic assertions mechanism that can be disabled by compiling the code in "release" mode.

  The image data is packed by the pixel, to optimize multichannel processing efficiency. The image class interface is restricted to pixel and metadata access, processing and I/O are entirely delegated to external functions. The data is accessible in multiple ways: value and pixel wise, using float coordinates with nearest and bilinear interpolation, viewing the image with clamping, direct access to the data and so on.

  An auxiliary structured attribute contains all the metadata information, currently this data consists in the specification of the numeric range of the values, the image color space, the white point and the gamma compression of the image. The various functions checks these data before processing, failing if the image is not compatible with the operation in course.
*/
template<int Channels=3, typename ScalarType=double, bool Safe=true> 
class Image //TODO: Safe=false prima di mettere in vcg
{
private:
  // Data is private but controlled data wiews are
  // accessible from accessors

  /// width of the image
  int _width; 
  /// height of the image
  int _height; 
  /// data buffer
  ScalarType *_data; 

// constructors
public:
  /// the auxiliary structured attribute that contains all the metadata information.
  ImgAttributes<ScalarType> attributes;


  /*! \brief default image constructor

    Creates a 0 x 0 pixel image with no data and default attributes
  */
  Image() // creates a 0 x 0 image
  :_width(0),_height(0),_data(NULL),attributes()
  {
    STATIC_ASSERT( Channels>0 );
    STATIC_FLOAT_OR_DOUBLE_TYPECHECK( ScalarType );
  }

  /*! \brief (deep) copy constructor when all template parameters matches

    An explicit copy constructor is needed because when all template parameters matches the templated copy constructor is not consideredand and a wrong copy constructor is synthetized by compiler.

    \param image the image to be copied
  */
  Image(const Image<Channels,ScalarType,Safe> &image) // copy constructor (deep copy)
  {
    STATIC_ASSERT( Channels>0 );
    STATIC_FLOAT_OR_DOUBLE_TYPECHECK( ScalarType );

    assert(image._width > 0);
    assert(image._height > 0);
    assert(image._data != NULL);
    if(Safe) {
      if(image._width <= 0) throw ImageException("Image(Image): Nonpositive width");
      if(image._height <= 0) throw ImageException("Image(Image): Nonpositive height");
      if(image._data == NULL) throw ImageException("Image(Image): NULL data");     
    }    
    _width = image._width;
    _height = image._height;
    _data = new ScalarType[Channels * _width * _height];
    attributes=image.attributes;

    memcpy(_data, image._data, sizeof(ScalarType) * Channels * _width * _height);
  }

  /*! \brief (deep) copy constructor when some template parameters differs

    \param image the image to be copied
  */
  template<typename OtherScalarType, bool OtherSafe> 
  Image(const Image<Channels,OtherScalarType,OtherSafe> &image) // templated copy constructor (deep copy)
  :_width(0),_height(0),_data(NULL)
  {
    STATIC_ASSERT( Channels>0 );
    STATIC_FLOAT_OR_DOUBLE_TYPECHECK( ScalarType );

    assert(image._width > 0);
    assert(image._height > 0);
    assert(image._data != NULL);
    if(Safe || OtherSafe) {
      if(image._width <= 0) throw ImageException("Image(Image): Nonpositive width");
      if(image._height <= 0) throw ImageException("Image(Image): Nonpositive height");
      if(image._data == NULL) throw ImageException("Image(Image): NULL data");     
    }    
    _width = image._width;
    _height = image._height;
    _data = new ScalarType[Channels * _width * _height];
    attributes=image.attributes;

    if(typeid( ScalarType ) == typeid( OtherScalarType ))
      memcpy(_data, image._data, sizeof(ScalarType) * Channels * _width * _height);
    else
      for(int offset=0;offset< Channels * _width * _height; ++offset)
        _data[offset] = static_cast<ScalarType>(image._data[offset]);
  }

  /*! \brief blank image constructor

    Creates an arg_width x arg_height pixel image with 0-valued data and default attributes
    \param arg_width the width of the blank image
    \param arg_height the height of the blank image
  */
  Image(int arg_width,int arg_height)
  :_width(arg_width),_height(arg_height),_data(NULL),attributes()
  {
    STATIC_ASSERT( Channels>0 );
    STATIC_FLOAT_OR_DOUBLE_TYPECHECK( ScalarType );

    assert(arg_width>0);
    assert(arg_height>0);
    if(Safe) {
      if(arg_width <= 0) throw ImageException("Image(int,int): Nonpositive width");
      if(arg_height <= 0) throw ImageException("Image(int,int): Nonpositive height");
    }
    _data = new ScalarType[Channels * _width * _height];
    for(int offset=0;offset< Channels * _width * _height; ++offset)
      _data[offset] = 0.0f;
  }  

  /// the destructor
  ~Image()
  {
    attributes.reset();
    if(_data!=NULL){
      delete [] _data;
      _data=NULL;
    }
  }

// public functions
public:
  /*! \brief assignment operator (deep copy)

    \param image the image to be assigned to this instance
    \return a pointer to this instance
  */
  template<typename OtherScalarType, bool OtherSafe> 
  inline Image< Channels,ScalarType,Safe> & operator =(const Image<Channels,OtherScalarType,OtherSafe> &image)
  {
    assert(image._width > 0);
    assert(image._height > 0);
    assert(image._data != NULL);
    if(Safe || OtherSafe) {
      if(image._width <= 0) throw ImageException("operator =: Nonpositive width");
      if(image._height <= 0) throw ImageException("operator =: Nonpositive height");
      if(image._data == NULL) throw ImageException("operator =: NULL data");     
    }    
    _width = image._width;
    _height = image._height;
    _data = new ScalarType[Channels * _width * _height];
    attributes=image.attributes;

    if(typeid( ScalarType ) == typeid( OtherSafe ))
      memcpy(_data, image._data, sizeof(ScalarType) * Channels * _width * _height);
    else
      for(int offset=0;offset < Channels * _width * _height; ++offset)
        _data[offset] = static_cast<ScalarType>(image._data[offset]);
    return *this;
  }

  /*! \brief blanks and change the image dimensions

    Delete the current image data and create an arg_width x arg_height pixel image with 0-valued data and default attributes
    \param arg_width the width of the blank image
    \param arg_height the height of the blank image
  */
  inline void setZero(int arg_width, int arg_height)
  {
    assert(arg_width>0);
    assert(arg_height>0);
    if(Safe) {
      if(arg_width <= 0) throw ImageException("setZero: Nonpositive width");
      if(arg_height <= 0) throw ImageException("setZero: Nonpositive height");
    }
    if(_data!=NULL){
      delete [] _data;
      _data=NULL;
    }    
    _width = arg_width;
    _height = arg_height;
    _data = new ScalarType[Channels * _width * _height];
    attributes.reset();
    
    for(int offset=0;offset< Channels * _width * _height; ++offset)
      _data[offset] = 0.0f;
  }

  /*! \brief delete the image data

    Delete the current image data and create an 0 x 0 pixel image with no data.
  */
  inline void deleteData()
  {
    if(_data!=NULL){
      delete [] _data;
      _data=NULL;
    }
    _width = 0;
    _height = 0;    
  } 


  /*! \brief get all the values of a pixel

    \param x the horizontal coordinate of the pixel
    \param y the vertical coordinate of the pixel
    \param ret_pixel return parameter that is filled with the pixel values
  */
  inline void getPixel(int x, int y, ScalarType (& ret_pixel)[Channels]) const
  {
    assert( _data != NULL );
    assert( x >= 0 && x < _width );
    assert( y >= 0 && y < _height );
    if( Safe ){
      if ( _data == NULL ) throw ImageException("getPixel: NULL data");
      if ( !( x >= 0 && x < _width ) ) throw ImageException("getPixel: x out of bounds");
      if ( !( y >= 0 && y < _height ) ) throw ImageException("getPixel: y out of bounds");      
    }
    for (int channel=0;channel<Channels;++channel)
      ret_pixel[channel] = _data[ (x + y * _width) * Channels + channel ];
  }

  /*! \brief set all the values of a pixel

    \param x the horizontal coordinate of the pixel
    \param y the vertical coordinate of the pixel
    \param pixel the pixel values that are assigned to the pixel
  */
  inline void setPixel(int x, int y, const ScalarType (& pixel)[Channels])
  {
    assert( _data != NULL );
    assert( x >= 0 && x < _width );
    assert( y >= 0 && y < _height );
    if( Safe ){
      if ( _data == NULL ) throw ImageException("setPixel: NULL data");
      if ( !( x >= 0 && x < _width ) ) throw ImageException("setPixel: x out of bounds");
      if ( !( y >= 0 && y < _height ) ) throw ImageException("setPixel: y out of bounds");      
    }
    for (int channel=0;channel<Channels;++channel)
      _data[ (x + y * _width) * Channels + channel ] = pixel[channel];
  }

  /*! \brief get a single value of a pixel

    \param x the horizontal coordinate of the pixel
    \param y the vertical coordinate of the pixel
    \param channel the channel index
    \return the value of the channel at the pixel
  */
  inline ScalarType getValue(int x, int y, int channel) const
  {
    assert( _data != NULL );
    assert( x >= 0 && x < _width );
    assert( y >= 0 && y < _height );
    assert( channel >=0 && channel < Channels );
    if( Safe ){
      if ( _data == NULL ) throw ImageException("getFloat: NULL data");
      if ( !( x >= 0 && x < _width ) ) throw ImageException("getFloat: x out of bounds");
      if ( !( y >= 0 && y < _height ) ) throw ImageException("getFloat: y out of bounds");
      if ( !( channel >=0 && channel < Channels ) ) throw ImageException("channel out of bounds");
    }
    return _data[ (x + y * _width) * Channels + channel ];
  }

  /*! \brief set a single value of a pixel

    \param x the horizontal coordinate of the pixel
    \param y the vertical coordinate of the pixel
    \param channel the channel index
    \param value the value that is assigned to the channel at the pixel
  */
  inline void setValue(int x, int y, int channel, ScalarType value)
  {
    assert( _data != NULL );
    assert( x >= 0 && x < _width );
    assert( y >= 0 && y < _height );
    assert( channel >=0 && channel < Channels );
    if( Safe ){
      if ( _data == NULL ) throw ImageException("setFloat: NULL data");
      if ( !( x >= 0 && x < _width ) ) throw ImageException("setFloat: x out of bounds");
      if ( !( y >= 0 && y < _height ) ) throw ImageException("setFloat: y out of bounds");      
      if ( !( channel >=0 && channel < Channels ) ) throw ImageException("channel out of bounds");
    }
    _data[(x + y * _width) * Channels + channel] = value;
  }

  /*! \brief get all the values of a pixel, clamping if the coordinates are out of bounds

    \param x the horizontal coordinate of the pixel
    \param y the vertical coordinate of the pixel
    \param ret_pixel return parameter that is filled with the pixel values
  */
  inline void getPixelAsClamped(int x, int y, ScalarType (& ret_pixel)[Channels]) const
  {
    getPixel(x<0?0:(x<_width?x:_width-1), y<0?0:(y<_height?y:_height-1), ret_pixel);
  }

  /*! \brief get a single value of a pixel, clamping if the coordinates are out of bounds

    \param x the horizontal coordinate of the pixel
    \param y the vertical coordinate of the pixel
    \param channel the channel index
    \return the value of the channel at the pixel
  */
  inline float getValueAsClamped(int x, int y, int channel) const
  {
    return getValue(x<0?0:(x<_width?x:_width-1), y<0?0:(y<_height?y:_height-1),channel);
  }

  /*! \brief get all the values of a pixel, rounding the floating coordinates to the nearest pixel

    \param x the horizontal coordinate of the pixel
    \param y the vertical coordinate of the pixel
    \param ret_pixel return parameter that is filled with the pixel values
  */
  inline void nearestPixel(float x, float y, ScalarType (& ret_pixel)[Channels]) const
  {
  getPixel(static_cast<int>(floor(x+0.5f)),static_cast<int>(floor(y+0.5f)), ret_pixel);
  }

///// chi gli serve se la implementa, ora non c'ho tempo
//  inline void bilinearPixel(float x, float y, const ScalarType[] &pixel) const
//  {
//    assert( _data != NULL );
//    assert( x >= 0.0f && x <= _width-1.0f );
//    assert( y >= 0.0f && y <= _height-1.0f );
//    if( SAFE ){
//      if ( _data == NULL ) throw ImageException("bilinearPixel: NULL data");
//      if ( !( x >= 0.0f && x <= _width-1.0f ) ) throw ImageException("bilinearPixel: x out of bounds");
//      if ( !( y >= 0.0f && y <= _height-1.0f ) ) throw ImageException("bilinearPixel: y out of bounds");      
//    }    
//    
//  int x1 = static_cast<int>(floor(x));
//    int y1 = static_cast<int>(floor(y));
//  float a = 1.0f;
//  float b = 1.0f;
//  
//  if( x1 == _width )
//      x1--;
//  else
//      a = x - x1;
//  if( y1 == _height )
//      y1--;
//  else
//      b = y - y1;
//    
//  float a1 = (1.0f - a);
//    float b1 = (1.0f - b);
//    
//  int i = x1 + y1 * _width;
//
//  return _data[i] * a1 * b1
//         + _data[i + 1] * a * b1 
//         + _data[i + _width] * a1 * b + 
//         _data[i + _width + 1] * a * b;  
//  }

// accessors
public:
  /*! \brief get the width of the image

    \return the width of the image
  */
  inline int width() const
  {
    return _width;
  }

  /*! \brief get the height of the image

    \return the height of the image
  */
  inline int height() const
  {
    return _height;
  }

  /*! \brief set the width of the image

    \warning changing the image size without changing the data buffer accordingly can lead to memory errors
    \param width the new width of the image
  */
  inline void setWidth(int width)
  {
     assert(width > 0);
     if(Safe){
       if(width <= 0) throw ImageException("setWidth: Nonpositive width");
    }
    _width=width;
  }

  /*! \brief set the height of the image

    \warning changing the image size without changing the data buffer accordingly can lead to memory errors
    \param height the new height of the image
  */
   inline void setHeight(int height)
  {
    assert(height > 0);
    if(Safe){
       if(height <= 0) throw ImageException("setHeight: Nonpositive height");
    }
    _height=height;
  }

  /*! \brief get a const pointer to the image databuffer

    \warning this function exposes the internal structure of the image, mind your accesses.
    \return a const pointer to the image databuffer
  */
  inline ScalarType* dataValues() const
  {
    return _data;
  }

  /*! \brief get the size of the image databuffer

    \return the size of the image databuffer
  */
  inline int dataValuesSize() const
  {
    return _width * _height * Channels;
  }


  /*! \brief set the image databuffer

    \warning this function modifies the internal structure of the image!
    \param data a pointer to the new image databuffer
  */
  inline void setValues(ScalarType* data)
  {
    assert(data!=NULL);
    if(Safe){
      if(data == NULL) throw ImageException("setValues: NULL data");
    }
    _data = data;
  }

  /*! \brief checks if the given coordinates are inside the image bounds

    \param x the horizontal coordinate
    \param y the vertical coordinate
    \return true if the given coordinates are inside the image bounds, false otherwise
  */
  inline bool isInside(int x, int y) const
  {
    return x >= 0 && x < _width && y >= 0 && y < _height;
  }

  /*! \brief checks if the given float coordinates are inside the image bounds

    \param x the horizontal coordinate
    \param y the vertical coordinate
    \return true if the given coordinates are inside the image bounds, false otherwise
  */
  inline bool isInside(float x, float y) const
  {
    return x >= 0.0f && x <= _width-1.0f && y >= 0.0f && y <= _height-1.0f;
  }

  /*! \brief checks if the image has been initialized

    \return true if the data is not null and the minimun side is longer than 0 pixel, false otherwise
  */
  inline bool isValid() const
  {
    return _data != NULL && _width > 0 && _height > 0;
  }

  /*! \brief get the number of channels of the image

    \return the number of channels of the image
  */
  inline int channels() const
  {
    return Channels;
  }

  /*! \brief set a value in the image databuffer

    \warning this function modifies the internal structure of the image!
    \param i the index of the databuffer
    \param value the value to set in the i position of the databuffer
  */
  inline void setRawValue(int i, ScalarType value)
  {
     _data[i] = value;
  }
  
};

} //end namespace img

#endif /*IMG_IMAGE_H_*/
