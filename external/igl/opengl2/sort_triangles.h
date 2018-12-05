// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OPENGL2_SORT_TRIANGLES_H
#define IGL_OPENGL2_SORT_TRIANGLES_H

#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  namespace opengl2
  {
    template <
      typename DerivedV,
      typename DerivedF,
      typename DerivedFF,
      typename DerivedI>
    IGL_INLINE void sort_triangles(
      const Eigen::MatrixBase<DerivedV> & V,
      const Eigen::MatrixBase<DerivedF> & F,
      Eigen::MatrixBase<DerivedFF> & FF,
      Eigen::MatrixBase<DerivedI> & I);
    template <
      typename DerivedV,
      typename DerivedF,
      typename DerivedFF,
      typename DerivedI>
    IGL_INLINE void sort_triangles_slow(
      const Eigen::MatrixBase<DerivedV> & V,
      const Eigen::MatrixBase<DerivedF> & F,
      Eigen::MatrixBase<DerivedFF> & FF,
      Eigen::MatrixBase<DerivedI> & I);
    //template <
    //  typename DerivedV,
    //  typename DerivedF,
    //  typename DerivedMV,
    //  typename DerivedP,
    //  typename DerivedFF,
    //  typename DerivedI>
    //IGL_INLINE void sort_triangles_robust(
    //  const Eigen::MatrixBase<DerivedV> & V,
    //  const Eigen::MatrixBase<DerivedF> & F,
    //  const Eigen::MatrixBase<DerivedMV> & MV,
    //  const Eigen::MatrixBase<DerivedP> & P,
    //  Eigen::MatrixBase<DerivedFF> & FF,
    //  Eigen::MatrixBase<DerivedI> & I);
    //template <
    //  typename DerivedV,
    //  typename DerivedF,
    //  typename DerivedFF,
    //  typename DerivedI>
    //IGL_INLINE void sort_triangles_robust(
    //  const Eigen::MatrixBase<DerivedV> & V,
    //  const Eigen::MatrixBase<DerivedF> & F,
    //  Eigen::MatrixBase<DerivedFF> & FF,
    //  Eigen::MatrixBase<DerivedI> & I);
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "sort_triangles.cpp"
#endif

#endif

