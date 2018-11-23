#include "kernel.pckh"

#ifdef DIM3

Sign predicate(orienth)(
    point(p0), point(p1), point(p2), point(p3), point(p4),
    scalar h0, scalar h1, scalar h2, scalar h3, scalar h4
    DIM
) {

  scalar a11 = p1_0 - p0_0 ;
  scalar a12 = p1_1 - p0_1 ;
  scalar a13 = p1_2 - p0_2 ;
  scalar a14 = h1 - h0 ;


  scalar a21 = p2_0 - p0_0 ;
  scalar a22 = p2_1 - p0_1 ;
  scalar a23 = p2_2 - p0_2 ;
  scalar a24 = h2 - h0 ;

  scalar a31 = p3_0 - p0_0 ;
  scalar a32 = p3_1 - p0_1 ;
  scalar a33 = p3_2 - p0_2 ;
  scalar a34 = h3 - h0 ;

  scalar a41 = p3_0 - p0_0 ;
  scalar a42 = p3_1 - p0_1 ;
  scalar a43 = p3_2 - p0_2 ;
  scalar a44 = h4 - h0 ;

  scalar Delta = det4x4(
            a11,a12,a13,a14,
            a21,a22,a23,a24,
            a31,a32,a33,a34,
            a41,a42,a43,a44
         ) ;

  return sign(Delta);
}


#endif
