// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_CGAL_INTERSECT_OTHER_H
#define IGL_COPYLEFT_CGAL_INTERSECT_OTHER_H
#include "../../igl_inline.h"
#include "RemeshSelfIntersectionsParam.h"

#include <Eigen/Dense>

namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      // INTERSECT_OTHER Given a triangle mesh (VA,FA) and another mesh (VB,FB)
      // find all pairs of intersecting faces. Note that self-intersections are
      // ignored.
      // 
      // Inputs:
      //   VA  #V by 3 list of vertex positions
      //   FA  #F by 3 list of triangle indices into VA
      //   VB  #V by 3 list of vertex positions
      //   FB  #F by 3 list of triangle indices into VB
      //   params   whether to detect only and then whether to only find first
      //     intersection
      // Outputs:
      //   IF  #intersecting face pairs by 2 list of intersecting face pairs,
      //     indexing FA and FB
      //   VVAB  #VVAB by 3 list of vertex positions
      //   FFAB  #FFAB by 3 list of triangle indices into VVA
      //   JAB  #FFAB list of indices into [FA;FB] denoting birth triangle
      //   IMAB  #VVAB list of indices stitching duplicates (resulting from
      //     mesh intersections) together
      template <
        typename DerivedVA,
        typename DerivedFA,
        typename DerivedVB,
        typename DerivedFB,
        typename DerivedIF,
        typename DerivedVVAB,
        typename DerivedFFAB,
        typename DerivedJAB,
        typename DerivedIMAB>
      IGL_INLINE bool intersect_other(
        const Eigen::MatrixBase<DerivedVA> & VA,
        const Eigen::MatrixBase<DerivedFA> & FA,
        const Eigen::MatrixBase<DerivedVB> & VB,
        const Eigen::MatrixBase<DerivedFB> & FB,
        const RemeshSelfIntersectionsParam & params,
        Eigen::MatrixBase<DerivedIF> & IF,
        Eigen::MatrixBase<DerivedVVAB> & VVAB,
        Eigen::MatrixBase<DerivedFFAB> & FFAB,
        Eigen::MatrixBase<DerivedJAB>  & JAB,
        Eigen::MatrixBase<DerivedIMAB> & IMAB);
      // Legacy wrapper for detect only using common types.
      //
      // Inputs:
      //   VA  #V by 3 list of vertex positions
      //   FA  #F by 3 list of triangle indices into VA
      //   VB  #V by 3 list of vertex positions
      //   FB  #F by 3 list of triangle indices into VB
      //   first_only  whether to only detect the first intersection.
      // Outputs:
      //   IF  #intersecting face pairs by 2 list of intersecting face pairs,
      //     indexing FA and FB
      // Returns true if any intersections were found
      //
      // See also: remesh_self_intersections
      IGL_INLINE bool intersect_other(
        const Eigen::MatrixXd & VA,
        const Eigen::MatrixXi & FA,
        const Eigen::MatrixXd & VB,
        const Eigen::MatrixXi & FB,
        const bool first_only,
        Eigen::MatrixXi & IF);
    }
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "intersect_other.cpp"
#endif
  
#endif

