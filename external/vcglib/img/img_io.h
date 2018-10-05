#ifndef IMG_IO_H_
#define IMG_IO_H_

// minimal input/output

#include <iostream>
#include <fstream>
#include <limits>
namespace img {

template<int Channels, typename ScalarType, bool Safe>
inline bool saveNativeFormat(const Image<Channels,ScalarType,Safe> &image, const char *filename)
{
  assert(image.isValid());
  if(Safe){
    if(!image.isValid()) throw ImageException("Invalid image");
  }
  using namespace std;
  ofstream output (filename, ios::out|ios::binary);
  
  if (output.is_open()) {
    int channels=Channels;
    output.write(reinterpret_cast<const char*>(& channels), sizeof(int));
    int scalartype=0;
    if(typeid(ScalarType) == typeid(float)) {
      scalartype=1;
    } else if(typeid(ScalarType) == typeid(double)) {
      scalartype=2;
    } else {
      assert(0);
    }
    output.write(reinterpret_cast<const char*>(& scalartype), sizeof(int));
    int width=image.width();
    output.write(reinterpret_cast<const char*>(& width), sizeof(int));
    int height=image.height();
    output.write(reinterpret_cast<const char*>(& height), sizeof(int));
    output.write(reinterpret_cast<const char*>(& image.attributes), sizeof(ImgAttributes<ScalarType>));    
    output.write(reinterpret_cast<const char*>(image.dataValues()), sizeof(ScalarType) * image.dataValuesSize());
    output.flush();
    output.close();
    return true;
  }
  if(Safe)
    throw ImageException("Unable to save file");
  return false;
}

template<int Channels, typename ScalarType, bool Safe>
inline bool openNativeFormat(const char *filename, Image<Channels,ScalarType,Safe> &image)
{

  using namespace std;
  ifstream input (filename, ios::in|ios::binary);
  if (input.is_open()) {    
    int channels;
    input.read(reinterpret_cast<char*>(&channels), sizeof(int));
    assert(channels==Channels);
    int scalartype;
    input.read(reinterpret_cast<char*>(&scalartype), sizeof(int));
    if(typeid(ScalarType) == typeid(float)) {
      assert(scalartype==1);
    } else if(typeid(ScalarType) == typeid(double)) {
      assert(scalartype==2);
    } else { 
      assert(0);
    }
    int width;
    input.read(reinterpret_cast<char*>(&width), sizeof(int));
    int height;
    input.read(reinterpret_cast<char*>(&height), sizeof(int));
    image.setZero(width,height);
    input.read(reinterpret_cast<char*>(& image.attributes), sizeof(ImgAttributes<ScalarType>));    
    input.read(reinterpret_cast<char*>(image.dataValues()), sizeof(ScalarType) * image.dataValuesSize());
    input.close();
    return true;
  }
  if(Safe)
    throw ImageException("Unable to open file");
  return false;
}

} //end namespace img

#endif /*IMG_IO_H_*/
