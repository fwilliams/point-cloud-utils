#include <stdio.h>
#include<vcg/math/base.h>
#include<vcg/space/point3.h>
#include<vcg/space/point4.h>
#include<vcg/space/color4.h>

using namespace vcg;
// The shortest, simplest, dulliest introduction to the VCG Library
int main(int argc, char *argv[])
{
  printf("Hello Library!\n");

  // classical point types. 
  // Point3f is just a typedef for Point3<float>
  Point3f pp0(0,1,2);
  Point3f pp1(2,1,0);

  // classical overloading of math operators
  Point3f pp2=pp1+pp0;

  //you can access to the components of a point with three different access methods
  // [0] [1] [2]            <-- Preferred style
  // .X() .Y() .Z()
  // .V(0) .V(1) .V(2)
  printf("pp2: %f %f %f \n",pp2[0], pp2.Y(),pp2.V(2));
  
  // Warning no implicit casts between different types
  // Conversions are explicit
  Point3i ppi=Point3i::Construct(pp1+pp0);
  
  Point4i size(0,0,1,1);
  
  // Colors are specialized Point4<unsigned char> 
  // with a specialized constructor 

  Color4b cb(Color4b::LightBlue);   
  Color4f cf(Color4f::LightBlue);

  Color4b cbi; cbi.Import(cf);
  printf("ci %i %i %i %i\n",cbi.V(0),cbi.V(1),cbi.V(2),cbi.V(3));
  return -1;
}