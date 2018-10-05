#ifndef IMG_BASE_H_
#define IMG_BASE_H_
/*! \file img_base.h
    \brief basic definitions for the img module

    This header contains the basic module definitions.
*/

/// base of the static assertion mechanism
template<bool> struct NON_TRUE_EXPR_CompileTimeError;
/// partial instantiation for the static assertion mechanism
template<> struct NON_TRUE_EXPR_CompileTimeError<true> {};

/// the static assertion mechanism
#define STATIC_ASSERT(exp) (NON_TRUE_EXPR_CompileTimeError< (exp) >())  

/// base of the static typecheck mechanism
template<typename> struct NON_FLOAT_OR_DOUBLE_TYPE_CompileTimeError;
/// partial instantiation for the static typecheck mechanism
template<> struct NON_FLOAT_OR_DOUBLE_TYPE_CompileTimeError<float> {};
/// partial instantiation for the static typecheck mechanism
template<> struct NON_FLOAT_OR_DOUBLE_TYPE_CompileTimeError<double> {};

/// the static typecheck mechanism
#define STATIC_FLOAT_OR_DOUBLE_TYPECHECK(type) (NON_FLOAT_OR_DOUBLE_TYPE_CompileTimeError< type >())

/// define NULL pointer value
#ifndef NULL
 #ifdef __cplusplus
  #define NULL 0
 #else
  #define NULL ((void *)0)
 #endif
#endif

#include <assert.h>
#include <math.h>
#include <exception>
#include <typeinfo>

/*! \brief the img module namespace

  this is the main image module namespace.
*/
namespace img {

/*! \brief the basic exception class

  this is the basic image exception class, it simply carries an error string to the console.
*/
class ImageException: public std::exception
{
public:
  /// the error string
  const char *message;
  /// default constructor
  ImageException():exception(),message("no message"){}
  /*! \brief message carrying constructor

    \param arg_message the error string
  */
  ImageException(const char *arg_message):exception(),message(arg_message){}
  /// the destructor
  virtual ~ImageException () throw (){}
};

} //end namespace img

#endif /*IMG_BASE_H_*/
