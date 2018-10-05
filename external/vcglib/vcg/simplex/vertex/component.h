/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#ifndef __VCG_MESH
#error "This file should not be included alone. It is automatically included by complex.h"
#endif
#ifndef __VCG_VERTEX_PLUS_COMPONENT
#define __VCG_VERTEX_PLUS_COMPONENT

namespace vcg {
namespace vertex {
/** \addtogroup VertexComponentGroup
  @{

*/
/*------------------------- Base Classes  -----------------------------------------*/

  template <class S>
  struct CurvatureDirBaseType{
          typedef Point3<S> VecType;
          typedef  S   ScalarType;
          CurvatureDirBaseType () {}
          Point3<S>max_dir,min_dir; // max and min curvature direction
          S k1,k2;// max and min curvature values
  };

/*------------------------- EMPTY CORE COMPONENTS -----------------------------------------*/

template <class TT> class EmptyCore: public TT {
public:
  typedef int FlagType;
  int &Flags()       { assert(0); static int dummyflags(0);  return dummyflags; }
  int cFlags() const { return 0; }
  static bool HasFlags()   { return false; }

  typedef vcg::Point3f CoordType;
  typedef CoordType::ScalarType      ScalarType;
  CoordType &P()       { assert(0); static CoordType coord(0, 0, 0); return coord; }
  CoordType cP() const { assert(0); static CoordType coord(0, 0, 0);  assert(0); return coord; }
  static bool HasCoord()   { return false; }
  inline bool IsCoordEnabled() const { return TT::VertexType::HasCoord();}

  typedef vcg::Point3s NormalType;
  NormalType &N()       { assert(0); static NormalType dummy_normal(0, 0, 0); return dummy_normal; }
  NormalType cN() const { assert(0); static NormalType dummy_normal(0, 0, 0); return dummy_normal; }
  static bool HasNormal()    { return false; }
  inline bool IsNormalEnabled() const { return TT::VertexType::HasNormal();}

  typedef float QualityType;
  QualityType &Q()       { assert(0); static QualityType dummyQuality(0); return dummyQuality; }
  QualityType cQ() const { assert(0); static QualityType dummyQuality(0); return dummyQuality; }
  static bool HasQuality()   { return false; }
  inline bool IsQualityEnabled() const { return TT::VertexType::HasQuality();}

  typedef vcg::Color4b ColorType;
  ColorType &C()       { static ColorType dumcolor(vcg::Color4b::White); assert(0); return dumcolor; }
  ColorType cC() const { static ColorType dumcolor(vcg::Color4b::White);  assert(0); return dumcolor; }
  static bool HasColor()   { return false; }
  inline bool IsColorEnabled() const { return TT::VertexType::HasColor();}

  typedef int  MarkType;
  void InitIMark()    {  }
  int cIMark()  const { assert(0); static int tmp=-1; return tmp;}
  int &IMark()        { assert(0); static int tmp=-1; return tmp;}
  static bool HasMark()   { return false; }
  inline bool IsMarkEnabled() const { return TT::VertexType::HasMark();}

  typedef ScalarType RadiusType;
  RadiusType &R()       { static ScalarType v = 0.0; assert(0 && "the radius component is not available"); return v; }
  RadiusType cR() const { static const ScalarType v = 0.0; assert(0 && "the radius component is not available"); return v; }
  static bool HasRadius()     { return false; }
  inline bool IsRadiusEnabled() const { return TT::VertexType::HasRadius();}

  typedef vcg::TexCoord2<float,1> TexCoordType;
  TexCoordType &T()       { static TexCoordType dummy_texcoord;  assert(0); return dummy_texcoord; }
  TexCoordType cT() const { static TexCoordType dummy_texcoord;  assert(0); return dummy_texcoord; }
  static bool HasTexCoord()   { return false; }
  inline bool IsTexCoordEnabled() const { return TT::VertexType::HasTexCoord();}

  typename TT::TetraPointer &VTp()       { static typename TT::TetraPointer tp = 0;  assert(0); return tp; }
  typename TT::TetraPointer cVTp() const  { static typename TT::TetraPointer tp = 0;  assert(0); return tp; }
  int &VTi() { static int z = 0; return z; }
  static bool HasVTAdjacency() { return false; }

  typename TT::FacePointer &VFp()       { static typename TT::FacePointer fp=0;  assert(0); return fp; }
  typename TT::FacePointer cVFp() const { static typename TT::FacePointer fp=0;  assert(0); return fp; }
  int &VFi()       { static int z=-1; assert(0); return z;}
  int cVFi() const { static int z=-1; assert(0); return z;}
  bool IsNull() const { return true; }
  static bool HasVFAdjacency()   { return false; }
  bool IsVFInitialized() const {return static_cast<const typename TT::VertexType *>(this)->cVFi()!=-1;}
  void VFClear() {
    if(IsVFInitialized()) {
      static_cast<typename TT::VertexPointer>(this)->VFp()=0;
      static_cast<typename TT::VertexPointer>(this)->VFi()=-1;
    }
  }

  typename TT::EdgePointer &VEp()       { static typename TT::EdgePointer ep=0;  assert(0); return ep; }
  typename TT::EdgePointer cVEp() const { static typename TT::EdgePointer ep=0;  assert(0); return ep; }
  int &VEi()       { static int z=-1; return z;}
  int cVEi() const { static int z=-1; return z;}
  static bool HasVEAdjacency()   {   return false; }
  bool IsVEInitialized() const {return static_cast<const typename TT::VertexType *>(this)->cVEi()!=-1;}
  void VEClear() {
    if(IsVEInitialized()) {
      static_cast<typename TT::VertexPointer>(this)->VEp()=0;
      static_cast<typename TT::VertexPointer>(this)->VEi()=-1;
    }
  }
  typename TT::HEdgePointer &VHp()       { static typename TT::HEdgePointer ep=0;  assert(0); return ep; }
  typename TT::HEdgePointer cVHp() const { static typename TT::HEdgePointer ep=0;  assert(0); return ep; }
  int &VHi()       { static int z=0; return z;}
  int cVHi() const { static int z=0; return z;}
  static bool HasVHAdjacency()   {   return false; }

  typedef float   CurScalarType;
  typedef float   ScalarTypeCur;
  typedef Point3f CurVecType;
  typedef Point2f CurvatureType;
  float &Kh()       { static float dummy = 0.f; assert(0);return dummy;}
  float &Kg()       { static float dummy = 0.f; assert(0);return dummy;}
  float cKh() const { static float dummy = 0.f; assert(0); return dummy;}
  float cKg() const { static float dummy = 0.f; assert(0); return dummy;}

  typedef CurvatureDirBaseType<float> CurvatureDirType;
  CurVecType &PD1()       {static CurVecType v(0,0,0); assert(0);return v;}
  CurVecType &PD2()       {static CurVecType v(0,0,0); assert(0);return v;}
  CurVecType cPD1() const {static CurVecType v(0,0,0); assert(0);return v;}
  CurVecType cPD2() const {static CurVecType v(0,0,0); assert(0);return v;}

  CurScalarType &K1()       { static ScalarType v = 0.0;assert(0);return v;}
  CurScalarType &K2()       { static ScalarType v = 0.0;assert(0);return v;}
  CurScalarType cK1() const {static ScalarType v = 0.0;assert(0);return v;}
  CurScalarType cK2() const  {static ScalarType v = 0.0;assert(0);return v;}

  static bool HasCurvature()			{ return false; }
  static bool HasCurvatureDir()			{ return false; }
  inline bool IsCurvatureEnabled() const { return TT::VertexType::HasCurvature();}
  inline bool IsCurvatureDirEnabled() const { return TT::VertexType::HasCurvatureDir();}

  template < class RightValueType>
  void ImportData(const RightValueType  & /*rVert*/ ) {
//			TT::ImportData( rVert);
  }
  static void Name(std::vector<std::string> & name){TT::Name(name);}
};

/*-------------------------- COORD ----------------------------------------*/
/*! \brief \em Generic Component: \b Geometric \b Position of the vertex
  Templated on the coordinate class. In practice you use one of the two specialized class Coord3f and Coord3d
  You can access to the coordinate of a vertex by mean of the P(),cP() member functions.
  */
template <class A, class T> class Coord: public T {
public:
  typedef A CoordType;
  typedef typename A::ScalarType      ScalarType;
  /// Return a const reference to the coordinate of the vertex
  inline const CoordType &P() const { return _coord; }
  /// Return a reference to the coordinate of the vertex
  inline       CoordType &P()       { return _coord; }
  /// Return a const reference to the coordinate of the vertex
  inline       CoordType cP() const { return _coord; }

  template < class RightValueType>
  void ImportData(const RightValueType  & rVert ) { if(rVert.IsCoordEnabled()) P().Import(rVert.cP()); T::ImportData( rVert); }
  static bool HasCoord()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Coord"));T::Name(name);}

private:
  CoordType _coord;
};
/// Specialized Coord Component in floating point precision.
template <class T> class Coord3f: public Coord<vcg::Point3f, T> {
public:	static void Name(std::vector<std::string> & name){name.push_back(std::string("Coord3f"));T::Name(name);}
};
/// Specialized Coord Component in double point precision.
template <class T> class Coord3d: public Coord<vcg::Point3d, T> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Coord3d"));T::Name(name);}
};

/*-------------------------- NORMAL ----------------------------------------*/
  /*! \brief \em Generic Component: \b %Normal of the vertex

    Templated on the Point3 class used to store the normal.
    In practice you use one of the two specialized class Normal3f and Normal3d.

    You can access to the normal of a vertex by mean of the N(),cN() member functions.

    \note Many algorithms assume that, for sake of precision coherence,
    the type of the normal is the same with respect to the type coord component.
    */

template <class A, class T> class Normal: public T {
public:
  typedef A NormalType;
  /// Return a const reference to the normal of the vertex
  inline const NormalType &N() const { return _norm; }
  /// Return a  reference to the normal of the vertex
  inline       NormalType &N()       { return _norm; }
  /// Return a const reference to the normal of the vertex
  inline       NormalType cN() const { return _norm; }
  template < class RightValueType>
  void ImportData(const RightValueType  & rVert ){
    if(rVert.IsNormalEnabled())  N().Import(rVert.cN());
    T::ImportData( rVert);
  }
  static bool HasNormal()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Normal"));T::Name(name);}

private:
  NormalType _norm;
};

template <class T> class Normal3s: public Normal<vcg::Point3s, T> {
public:static void Name(std::vector<std::string> & name){name.push_back(std::string("Normal3s"));T::Name(name);}
};
/// Specialized Normal component in floating point precision.
template <class T> class Normal3f: public Normal<vcg::Point3f, T> {
public:	static void Name(std::vector<std::string> & name){name.push_back(std::string("Normal3f"));T::Name(name);}
};
/// Specialized Normal component in double point precision.
template <class T> class Normal3d: public Normal<vcg::Point3d, T> {
public:	static void Name(std::vector<std::string> & name){name.push_back(std::string("Normal3d"));T::Name(name);}
};


/*-------------------------- INCREMENTAL MARK  ----------------------------------------*/
  /*! \brief Per vertex \b Incremental \b Mark

      It is just an int that allows to efficently un-mark the whole mesh. \sa UnmarkAll
      */

template <class T> class Mark: public T {
public:
  Mark():_imark(0){}
  /// Return a const reference to the incremental mark value
  inline const int &IMark() const { return _imark;}
  /// Return a reference to the incremental mark value
  inline       int &IMark()       { return _imark;}
  /// Return a const reference to the incremental mark value
  inline       int cIMark() const { return _imark;}
  static bool HasMark()      { return true; }
  inline void InitIMark()    { _imark = 0; }
  template < class RightValueType>
  void ImportData(const RightValueType  & rVert ) { if(rVert.IsMarkEnabled())  IMark() = rVert.cIMark(); T::ImportData( rVert); }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Mark"));T::Name(name);}

 private:
    int _imark;
};

/*-------------------------- TEXCOORD ----------------------------------------*/
  /*! \brief \em Generic Component: Per vertex \b Texture Coords

      Note that to have multiple different TexCoord for a single vertex
      (as it happens on atlas where a vertex can belong to two triangles
       mapped on different portionof the texture) you have two options:
      - explicit duplication of vertexes
      - use PerWedge Texture coords

      It is templated on the TextureCoord type. Usually you use the specialized classes TexCoord2f or TexCoord2d;
      See the TexCoord2 class to see how to access to texture coordinate values.

      */

template <class A, class TT> class TexCoord: public TT {
public:
  typedef A TexCoordType;

  /// Return a const reference to the Texture Coordinate
  const TexCoordType &T() const { return _t; }
        TexCoordType &T()       { return _t; }
        TexCoordType cT() const { return _t; }
    template < class RightValueType>
    void ImportData(const RightValueType  & rVert ) { if(rVert.IsTexCoordEnabled())  T() = rVert.cT(); TT::ImportData( rVert); }
  static bool HasTexCoord()   { return true; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("TexCoord"));TT::Name(name);}

private:
  TexCoordType _t;
};


template <class TT> class TexCoord2s: public TexCoord<TexCoord2<short,1>, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("TexCoord2s"));TT::Name(name);}
};
/// Specialized Texture component in floating point precision.
template <class TT> class TexCoord2f: public TexCoord<TexCoord2<float,1>, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("TexCoord2f"));TT::Name(name);}
};
/// Specialized Texture component in double precision.
template <class TT> class TexCoord2d: public TexCoord<TexCoord2<double,1>, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("TexCoord2d"));TT::Name(name);}
};

/*------------------------- FLAGS -----------------------------------------*/
  /*! \brief \em Component: Per vertex \b Flags

      This component stores a 32 bit array of bit flags.
      These bit flags are used for keeping track of selection, deletion, visiting etc. \sa \ref flags for more details on common uses of flags.
      */

template <class T> class BitFlags:  public T {
public:
  BitFlags(){_flags=0;}
  typedef int FlagType;
  inline const int &Flags() const {return _flags; }
  inline       int &Flags()       {return _flags; }
  inline       int cFlags() const {return _flags; }
  template < class RightValueType>
  void ImportData(const RightValueType  & rVert ) { if(RightValueType::HasFlags()) Flags() = rVert.cFlags(); T::ImportData( rVert); }
  static bool HasFlags()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("BitFlags"));T::Name(name);}

private:
  int  _flags;
};


/*-------------------------- Color  ----------------------------------*/
  /*! \brief \em Component: Per vertex \b Color
   *
   * Usually most of the library expects a color stored as 4 unsigned chars (so the component you use is a \c vertex::Color4b)
   * but you can also use float for the color components.
   */
template <class A, class T> class Color: public T {
public:
  Color():_color(vcg::Color4b::White) {}
  typedef A ColorType;
  inline const ColorType &C() const { return _color; }
  inline       ColorType &C()       { return _color; }
  inline       ColorType cC() const { return _color; }
  template < class RightValueType>
  void ImportData(const RightValueType  & rVert ) { if(rVert.IsColorEnabled()) C() = rVert.cC();  T::ImportData( rVert); }
  static bool HasColor()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Color"));T::Name(name);}

private:
  ColorType _color;
};

template <class TT> class Color4b: public Color<vcg::Color4b, TT> {
    public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Color4b"));TT::Name(name);}
};

/*-------------------------- Quality  ----------------------------------*/
  /*! \brief \em Component: Per vertex \b quality
The Quality Component is a generic place for storing a float. The term 'quality' is a bit misleading and it is due to its original storic meaning. You should intend it as a general purpose container.
\sa vcg::tri::UpdateColor for methods transforming quality into colors
\sa vcg::tri::UpdateQuality for methods to manage it

*/

template <class A, class TT> class Quality: public TT {
public:
  typedef A QualityType;
  Quality():_quality(0) {}

  inline const QualityType &Q() const { return _quality; }
  inline       QualityType &Q()       { return _quality; }
  inline       QualityType cQ() const {return _quality; }
  template < class RightValueType>
  void ImportData(const RightValueType  & rVert ) { if(rVert.IsQualityEnabled()) Q() = rVert.cQ(); TT::ImportData( rVert); }
  static bool HasQuality()   { return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("Quality"));TT::Name(name);}

private:
  QualityType _quality;
};

template <class TT> class Qualitys: public Quality<short, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualitys"));TT::Name(name);}
};
template <class TT> class Qualityf: public Quality<float, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualityf"));TT::Name(name);}
};
template <class TT> class Qualityd: public Quality<double, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Qualityd"));TT::Name(name);}
};

  /*-------------------------- Curvature   ----------------------------------*/

  /*! \brief \em Component: Per vertex basic \b curvature
    This component keeps the mean an gaussian curvature for a vertex. Used by some of the algorithms of vcg::tri::UpdateCurvature to store the computed curvatures.
      */
  template <class A, class TT> class Curvature: public TT {
  public:
    typedef Point2<A> CurvatureType;
    typedef typename CurvatureType::ScalarType ScalarTypeCur;
    const ScalarTypeCur &Kh() const { return _hk[0]; }
    const ScalarTypeCur &Kg() const { return _hk[1]; }
          ScalarTypeCur &Kh()       { return _hk[0]; }
          ScalarTypeCur &Kg()       { return _hk[1]; }
          ScalarTypeCur cKh() const { return _hk[0]; }
          ScalarTypeCur cKg() const { return _hk[1]; }

          template < class RightValueType>
          void ImportData(const RightValueType  & rVert ) {
            if(rVert.IsCurvatureEnabled()) {
              Kh() = rVert.cKh();
              Kg() = rVert.cKg();
            }
            TT::ImportData( rVert);
          }

    static bool HasCurvature()   { return true; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("Curvature"));TT::Name(name);}

  private:
    Point2<A> _hk;
  };


  template <class T> class Curvaturef: public Curvature< float, T> {
  public:	static void Name(std::vector<std::string> & name){name.push_back(std::string("Curvaturef"));T::Name(name);}
  };
  template <class T> class Curvatured: public Curvature<double , T> {
  public:	static void Name(std::vector<std::string> & name){name.push_back(std::string("Curvatured"));T::Name(name);}
  };

/*-------------------------- Curvature Direction ----------------------------------*/

  /*! \brief \em Component: Per vertex \b curvature \b directions
    This component keep the principal curvature directions. Used by some of the algorithms of vcg::tri::UpdateCurvature to store the computed curvatures.
      */

template <class A, class TT> class CurvatureDir: public TT {
public:
  typedef A CurvatureDirType;
    typedef typename CurvatureDirType::VecType CurVecType;
    typedef typename CurvatureDirType::ScalarType CurScalarType;

    CurVecType &PD1(){ return _curv.max_dir; }
    CurVecType &PD2(){ return _curv.min_dir; }
    const CurVecType &cPD1() const { return _curv.max_dir; }
    const CurVecType &cPD2() const { return _curv.min_dir; }

    CurScalarType &K1(){ return _curv.k1; }
    CurScalarType &K2(){ return _curv.k2; }
    const CurScalarType &cK1() const { return _curv.k1; }
    const CurScalarType &cK2() const { return _curv.k2; }
    template < class RightValueType>
    void ImportData(const RightValueType  & rVert ) {
      if(rVert.IsCurvatureDirEnabled()) {
        PD1().Import(rVert.cPD1());
        PD2().Import(rVert.cPD2());
        K1()  = rVert.cK1();  K2()  = rVert.cK2();
      }
      TT::ImportData( rVert);
    }

    static bool HasCurvatureDir()   { return true; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("CurvatureDir"));TT::Name(name);}

private:
  CurvatureDirType _curv;
};


template <class T> class CurvatureDirf: public CurvatureDir<CurvatureDirBaseType<float>, T> {
public:	static void Name(std::vector<std::string> & name){name.push_back(std::string("CurvatureDirf"));T::Name(name);}
};
template <class T> class CurvatureDird: public CurvatureDir<CurvatureDirBaseType<double>, T> {
public:	static void Name(std::vector<std::string> & name){name.push_back(std::string("CurvatureDird"));T::Name(name);}
};

/*-------------------------- Radius  ----------------------------------*/
  /*! \brief \em Component: Per vertex \b radius

    This component keep a floating point value meant to be the average distance from the surrounding vertices. Used in point clouds by some of the point splatting and MLS surface algorithms.
      */
  template <class A, class TT> class Radius: public TT {
  public:
    typedef A RadiusType;
    const RadiusType &R() const { return _radius; }
          RadiusType &R()       { return _radius; }
          RadiusType cR() const {return _radius; }
    template < class RightValueType>
    void ImportData(const RightValueType  & rVert ) { if(rVert.IsRadiusEnabled()) R() = rVert.cR(); TT::ImportData( rVert); }
    static bool HasRadius()   { return true; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("Radius"));TT::Name(name);}

  private:
    RadiusType _radius;
  };

template <class TT> class Radiusf: public Radius<float, TT> {
public: static void Name(std::vector<std::string> & name){name.push_back(std::string("Radiusf"));TT::Name(name);}
};


/*----------------------------- VEADJ ------------------------------*/
/*! \brief \em Component: Per vertex \b Vertex-Edge adjacency relation
    It stores a pointer to the first Edge of a list edges that is stored in a distributed way on the edges themselves.

\sa vcg::tri::UpdateTopology for functions that compute this relation
\sa iterators
*/

template <class T> class VEAdj: public T {
public:
  VEAdj(){_ep=0;_zp=-1;}
  typename T::EdgePointer &VEp()       {return _ep; }
  typename T::EdgePointer cVEp() const {return _ep; }
  int &VEi()       {return _zp; }
  int cVEi() const {return _zp; }
  template < class RightValueType>
  void ImportData(const RightValueType  & rVert ) {  T::ImportData( rVert); }
  static bool HasVEAdjacency()   {   return true; }
  static void Name(std::vector<std::string> & name){name.push_back(std::string("VEAdj"));T::Name(name);}

private:
  typename T::EdgePointer _ep ;
  int _zp ;
};

/*----------------------------- VFADJ ------------------------------*/
  /*! \brief \em Component: Per vertex \b Vertex-Face adjacency relation

It stores a pointer to the first face of a list of faces that is stored in a distributed way on the faces themselves.
Note that if you use this component it is expected that on the Face you use also the corresponding vcg::face::VFAdj component.

  \sa vcg::tri::UpdateTopology for functions that compute this relation
  \sa vcg::face::VFAdj
  \sa iterators
  */

  template <class T> class VFAdj: public T {
  public:
    VFAdj(){_fp=0;_zp=-1;}
    typename T::FacePointer &VFp()        { return _fp; }
    typename T::FacePointer cVFp() const  { return _fp; }
    int &VFi()       { return _zp; }
    int cVFi() const { return _zp; }
    bool IsNull() const { return _zp==-1;}
    template < class RightValueType>
    void ImportData(const RightValueType  & rVert ) { T::ImportData( rVert); }
    static bool HasVFAdjacency()   {   return true; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("VFAdj"));T::Name(name);}

  private:
    typename T::FacePointer _fp ;
    int _zp ;
  };

/*----------------------------- VHADJ ------------------------------*/

template <class T> class VHAdj: public T {
public:
    VHAdj(){_hp=0;_zp=-1;}
    typename T::HEdgePointer &VHp()       {return _hp; }
    typename T::HEdgePointer cVHp() const {return _hp; }
    int &VHi() {return _zp; }
    template < class RightValueType>
    void ImportData(const RightValueType  & rVert ) {  T::ImportData( rVert); }
    static bool HasVHAdjacency()   {   return true; }
    static void Name(std::vector<std::string> & name){name.push_back(std::string("VHAdj"));T::Name(name);}

private:
    typename T::HEdgePointer _hp ;
    int _zp ;
};

/*----------------------------- VTADJ ------------------------------*/

template <class T> class VTAdj: public T {
public:
    VTAdj() { _tp = 0; _zp=-1;}
    typename T::TetraPointer &VTp()       { return _tp; }
    typename T::TetraPointer cVTp() const { return _tp; }
    int &VTi() {return _zp; }
    static bool HasVTAdjacency() { return true; }
    static void Name( std::vector< std::string > & name ) { name.push_back( std::string("VTAdj") ); T::Name(name); }

private:
    typename T::TetraPointer _tp ;
    int _zp ;
};

  /** @} */   // End Doxygen VertexComponentGroup
  } // end namespace vert
}// end namespace vcg
#endif
