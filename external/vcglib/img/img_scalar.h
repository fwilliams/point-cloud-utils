#ifndef IMG_SCALAR_H_
#define IMG_SCALAR_H_
/*! \file img_scalar.h
    \brief functions that operate on the pixel values

    This header contains some utility functions for the scalar values of the image pixels.
*/

#include "img/img_base.h"

namespace img {

/*! \brief clamp a scalar value

  \param value the value that is to be clamped
  \param minval the minimum possible value after clamping
  \param maxval the the maximum possible value after clamping
  \return the clamped value
*/
template<typename ScalarType>
inline ScalarType clampValue(ScalarType value, ScalarType minval=ScalarType(0.0), ScalarType maxval=ScalarType(255.0))
{
  STATIC_FLOAT_OR_DOUBLE_TYPECHECK(ScalarType);
  if(value < minval) return minval;
  if(value > maxval) return maxval;
  return value;
}

/*! \brief convert a scalar value to the nearest integer

  \param value the value that is to be rounded
  \return the rounded value
*/
template<typename ScalarType>
inline int valueAsInt(ScalarType value)
{
  STATIC_FLOAT_OR_DOUBLE_TYPECHECK(ScalarType);
  return static_cast<int>(floor(static_cast<double>(value)+0.5f));
}

/*! \brief take the maximum between three scalar values

  \param a one of the three values
  \param b one of the three values
  \param c one of the three values
  \return the maximum value between a, b and c
*/
template<typename ScalarType>
inline ScalarType max3(ScalarType a, ScalarType b, ScalarType c)
{
  STATIC_FLOAT_OR_DOUBLE_TYPECHECK(ScalarType);
  if ( b<=a && c<=a ) return a;
  if ( a<=b && c<=b ) return b;
  return c;
}

/*! \brief take the minimum between three scalar values

  \param a one of the three values
  \param b one of the three values
  \param c one of the three values
  \return the minimum value between a, b and c
*/
template<typename ScalarType>
inline ScalarType min3(ScalarType a, ScalarType b, ScalarType c)
{
  STATIC_FLOAT_OR_DOUBLE_TYPECHECK(ScalarType);
  if ( a<=b && a<=c ) return a;
  if ( b<=a && b<=c ) return b;
  return c;
}

/*! \brief take the median between three scalar values

  \param a one of the three values
  \param b one of the three values
  \param c one of the three values
  \return the median value between a, b and c
*/
template<typename ScalarType>
inline ScalarType median3(ScalarType a, ScalarType b, ScalarType c)
{
  STATIC_FLOAT_OR_DOUBLE_TYPECHECK(ScalarType);
  if ( (b<=a && a<=c) || (c<=a && a<=b) ) return a;
  if ( (a<=b && b<=c) || (c<=b && b<=a) ) return b;
  return c;
}

/*! \brief compares two scalar values for (near) equality

  \param a one of the two values
  \param b one of the two values
  \param EPSILON the tolerance for considering a and b equal
  \return true if a differs from b for less than EPSILON, false otherwhise
*/
template<typename ScalarType>
inline bool almostEqual(ScalarType a, ScalarType b, ScalarType EPSILON=ScalarType(10e-5))
{
  STATIC_FLOAT_OR_DOUBLE_TYPECHECK(ScalarType);
  ScalarType diff=a-b;
  if (diff<0)
    diff *= ScalarType(-1);
  return diff < EPSILON;
}

} //end namespace img

#endif /*IMG_SCALAR_H_*/
