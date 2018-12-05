// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "model_proj_viewport.h"
#include "gl.h"

template <typename Derivedmodel, typename Derivedproj, typename Derivedviewport>
IGL_INLINE void igl::opengl2::model_proj_viewport(
    Eigen::MatrixBase<Derivedmodel> & model,
    Eigen::MatrixBase<Derivedproj> & proj,
    Eigen::MatrixBase<Derivedviewport> & viewport)
{
  Eigen::Matrix4d MV,P;
  Eigen::Vector4i VPi;
  glGetDoublev(GL_MODELVIEW_MATRIX,MV.data());
  glGetDoublev(GL_PROJECTION_MATRIX,P.data());
  glGetIntegerv(GL_VIEWPORT,VPi.data());
  viewport = VPi.cast<typename Derivedviewport::Scalar>();
  model = MV.cast<typename Derivedmodel::Scalar>();
  proj = P.cast<typename Derivedproj::Scalar>();
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::opengl2::model_proj_viewport<Eigen::Matrix<double, 4, 4, 0, 4, 4>, Eigen::Matrix<double, 4, 4, 0, 4, 4>, Eigen::Matrix<double, 4, 1, 0, 4, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, 4, 4, 0, 4, 4> >&, Eigen::MatrixBase<Eigen::Matrix<double, 4, 4, 0, 4, 4> >&, Eigen::MatrixBase<Eigen::Matrix<double, 4, 1, 0, 4, 1> >&);
template void igl::opengl2::model_proj_viewport<Eigen::Matrix<float, 4, 4, 0, 4, 4>, Eigen::Matrix<float, 4, 4, 0, 4, 4>, Eigen::Matrix<float, 4, 1, 0, 4, 1> >(Eigen::MatrixBase<Eigen::Matrix<float, 4, 4, 0, 4, 4> >&, Eigen::MatrixBase<Eigen::Matrix<float, 4, 4, 0, 4, 4> >&, Eigen::MatrixBase<Eigen::Matrix<float, 4, 1, 0, 4, 1> >&);
#endif
