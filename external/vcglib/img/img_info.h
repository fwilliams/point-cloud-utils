#ifndef IMG_INFO_H_
#define IMG_INFO_H_

// functions that extrapolate information from read-only images

namespace img {

template<int Channels,typename ScalarType, bool Safe>
inline ScalarType minValue(const Image<Channels,ScalarType,Safe> &image)
{
  assert(image.isValid());
  if(Safe){
    if(!image.isValid())  throw ImageException("Invalid image");
  }
  ScalarType* array = image.dataValues();
  int length =image.dataValuesSize();

  ScalarType min = array[0];
  for(int offset=0;offset<length;offset++)
    if(min > array[offset])
      min = array[offset];
  return min;
}
  
template<int Channels,typename ScalarType, bool Safe>
inline ScalarType maxValue(const Image<Channels,ScalarType,Safe> &image)
{
  assert(image.isValid());
  if(Safe){
    if(!image.isValid())  throw ImageException("Invalid image");
  }
  ScalarType* array = image.dataValues();
  int length =image.dataValuesSize();

  ScalarType max = array[0];
  for(int offset=0;offset<length;offset++)
    if(max < array[offset])
      max = array[offset];
      
  return max;
}

} //end namespace img

#endif /*IMG_INFO_H_*/
