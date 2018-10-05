/****************************************************************************
* NanoPLY                                                                   *
* NanoPLY is a C++11 header-only library to read and write PLY file         *
*                                                                           *
* Copyright(C) 2014-2015                                                    *
* Visual Computing Lab                                                      *
* ISTI - Italian National Research Council                                  *
*                                                                           *
* This Source Code Form is subject to the terms of the Mozilla Public       *
* License, v. 2.0. If a copy of the MPL was not distributed with this       *
* file, You can obtain one at http://mozilla.org/MPL/2.0/.                  *
*                                                                           *
****************************************************************************/


#ifndef NANOPLY_HPP
#define NANOPLY_HPP

#include <vector>
#include <unordered_map>
#include <tuple>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <cstdio>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <string>


// Avoid conflicting declaration of min/max macros in windows headers
#if !defined(NOMINMAX) && (defined(_WIN32) || defined(_WIN32_)  || defined(WIN32) || defined(_WIN64))
# define NOMINMAX
# ifdef max
#  undef   max
#  undef   min
# endif
#endif

namespace nanoply
{

  /** Error Type.
  *	Error type returned by the open of a PLY file.
  */
  typedef enum NNP_ERROR
  {
    NNP_OK = 0x0000,  /**< No error. */
    NNP_UNABLE_TO_OPEN = 0x0001,  /**< The file cannot be opend. */
    NNP_MISSING_HEADER = 0x0002,  /**< The file does not contain a valid PLY header. */
    NNP_MISSING_FORMAT = 0x0004,  /**< The file has an invalid internal format. */
    NNP_INVALID_ELEMENT = 0x0008,  /**< The file has an invalid element. */
    NNP_INVALID_PROPERTY = 0x0010	  /**< The file has an invalid property. */
  } ErrorCode;


  /** PLY Element Entity.
  *	Element that can be saved in PLY file.
  */
  typedef enum NNP_ELEM
  {
    NNP_UNKNOWN_ELEM = 0x0, /**< Unknown element. */
    NNP_VERTEX_ELEM = 0x1, /**< Vertex element. */
    NNP_EDGE_ELEM = 0x2, /**< Edge element. */
    NNP_FACE_ELEM = 0x4  /**< Face element. */
  } PlyElemEntity;


  /** PLY Entity.
  *  Property that can be saved in a PLY file.
  */
  typedef enum NNP_ENTITY
  {
    NNP_UNKNOWN_ENTITY = 0x00000000,	/**< Unknown property. */
    NNP_PX = 0x00000001,	/**< Position x cordinate. */
    NNP_PY = 0x00000002,	/**< Position y cordinate. */
    NNP_PZ = 0x00000004,	/**< Position z cordinate. */
    NNP_PXYZ = 0x00000007,	/**< Position (x, y, z). */
    NNP_NX = 0x00000010,	/**< Normal x component. */
    NNP_NY = 0x00000020,	/**< Normal y component. */
    NNP_NZ = 0x00000040,	/**< Normal z component. */
    NNP_NXYZ = 0x00000070,	/**< Normal (x, y, z). */
    NNP_CR = 0x00000100,	/**< Color red component. */
    NNP_CG = 0x00000200,	/**< Color green component. */
    NNP_CB = 0x00000400,	/**< Color blue component. */
    NNP_CRGB = 0x00000700,	/**< Color RGB. */
    NNP_CA = 0x00000800,	/**< Color alpha component. */
    NNP_CRGBA = 0x00000F00,	/**< Color RGBA. */
    NNP_DENSITY = 0x00000008,	/**< Density or Radius property. */
    NNP_SCALE = 0x00000080,	/**< Scale property. */
    NNP_TEXTUREU = 0x00001000,	/**< Texture coordinate u. */
    NNP_TEXTUREV = 0x00002000,	/**< Texture coordinate v. */
    NNP_TEXTURE2D = 0x00003000,	/**< Texture coordinate 2D. */
    NNP_TEXTUREW = 0x00004000,	/**< Texture coordinate w. */
    NNP_TEXTURE3D = 0x00007000,	/**< Texture coordinate 3D. */
    NNP_TEXTUREINDEX = 0x00008000,	/**< Texture index. */
    NNP_QUALITY = 0x00010000,	/**< Quality property. */
    NNP_REFLECTANCE = 0x00020000,	/**< Reflectance property. */
    NNP_BITFLAG = 0x00040000,	/**< Bit flags. */
    NNP_K1 = 0x00080000,	/**< Main curvaure value k1. */
    NNP_K2 = 0x00100000,	/**< Main curvaure value k2. */
    NNP_KG = 0x00200000,	/**< Gaussian curvature value. */
    NNP_KH = 0x00400000, /**< Mean curvature value. */
    NNP_K1DIR = 0x00800000, /**< Curvature direction k1. */
    NNP_K2DIR = 0x01000000, /**< Curvature direction k2. */
    NNP_EDGE_V1 = 0x02000000,	/**< Index of the first vertex of the edge. */
    NNP_EDGE_V2 = 0x04000000,	/**< Index of the second vertex of the edge. */
    NNP_FACE_VERTEX_LIST = 0x08000000, /**< List of vertex indices for the face. */
    NNP_FACE_WEDGE_COLOR = 0x10000000, /**< List of colors for wedge. */
    NNP_FACE_WEDGE_NORMAL = 0x20000000, /**< List of normals for wedge. */
    NNP_FACE_WEDGE_TEX = 0x40000000  /**< List of texture coordinates for wedge. */
  } PlyEntity;


  /** PLY Type.
  *  Type of a PLY property.
  */
  typedef enum NNP_PLYTYPE
  {
    NNP_UNKNOWN_TYPE = 0x000000, /**< Unknown type. */
    NNP_FLOAT32 = 0x000001, /**< Float. */
    NNP_FLOAT64 = 0x000002, /**< Double. */
    NNP_INT8 = 0x000004, /**< Char. */
    NNP_INT16 = 0x000008, /**< Short. */
    NNP_INT32 = 0x000010, /**< Int. */
    NNP_UINT8 = 0x000020, /**< Unsigned Char. */
    NNP_UINT16 = 0x000040, /**< Unsigned Short. */
    NNP_UINT32 = 0x000080, /**< Unsigned Int. */
    NNP_LIST_UINT8_UINT32 = 0x000100, /**< List (size Unsigned Char) of Unsigned Int.  */
    NNP_LIST_INT8_UINT32 = 0x000200, /**< List (size Char) of Unsigned Int.  */
    NNP_LIST_UINT8_INT32 = 0x000400, /**< List (size Unsigned Char) of Int.  */
    NNP_LIST_INT8_INT32 = 0x000800, /**< List (size Char) of Int. */
    NNP_LIST_UINT8_FLOAT32 = 0x001000, /**< List (size Unsigned Char) of Float. */
    NNP_LIST_INT8_FLOAT32 = 0x002000, /**< List (size Char) of Float. */
    NNP_LIST_UINT8_FLOAT64 = 0x004000, /**< List (size Unsigned Char) of Double. */
    NNP_LIST_INT8_FLOAT64 = 0x008000, /**< List (size Char) of Double. */
    NNP_LIST_UINT8_UINT8 = 0x010000, /**< List (size Unsigned Char) of Unsigned Char.  */
    NNP_LIST_INT8_UINT8 = 0x020000, /**< List (size Char) of Unsigned Char.  */
    NNP_LIST_UINT8_INT8 = 0x040000, /**< List (size Unsigned Char) of Char.  */
    NNP_LIST_INT8_INT8 = 0x080000, /**< List (size Char) of Char. */
    NNP_LIST_UINT8_UINT16 = 0x100000, /**< List (size Unsigned Char) of Unsigned Short.  */
    NNP_LIST_INT8_UINT16 = 0x200000, /**< List (size Char) of Unsigned Short.  */
    NNP_LIST_UINT8_INT16 = 0x400000, /**< List (size Unsigned Char) of Short.  */
    NNP_LIST_INT8_INT16 = 0x800000  /**< List (size Char) of Short. */
  } PlyType;

  template <typename TEnum>
  inline std::size_t hashEnum(const TEnum & value)
  {
	  const std::hash<unsigned int> h{};
	  return h((unsigned int)value);
  }

} // end namespace nanoply

// Cast of an enum to an unsigned int (maybe only for Android)
namespace std
{
	template <> struct hash<nanoply::PlyElemEntity> { std::size_t operator () (const nanoply::PlyElemEntity & t) const { return nanoply::hashEnum(t); } };
	template <> struct hash<nanoply::PlyEntity    > { std::size_t operator () (const nanoply::PlyEntity     & t) const { return nanoply::hashEnum(t); } };
	template <> struct hash<nanoply::PlyType      > { std::size_t operator () (const nanoply::PlyType       & t) const { return nanoply::hashEnum(t); } };
} // end namespace std

namespace nanoply
{

  /**
  * @cond HIDDEN_SYMBOLS
  */
  template < size_t T> struct SizeT {};

  typedef std::vector<std::string> NameVector;
  typedef std::unordered_map<PlyType, NameVector> TypeMap;
  typedef std::unordered_map<PlyType, NameVector>::iterator TypeMapIterator;
  typedef std::unordered_map<PlyEntity, NameVector> EntityMap;
  typedef std::unordered_map<PlyEntity, NameVector>::iterator EntityMapIterator;
  typedef std::unordered_map<PlyElemEntity, NameVector> ElementMap;
  typedef std::unordered_map<PlyElemEntity, NameVector>::iterator ElementMapIterator;

  /* Names used for the PlyType */
  static TypeMap mapType({
      { PlyType::NNP_UNKNOWN_TYPE, NameVector({ "unknonw" }) },
      { PlyType::NNP_FLOAT32, NameVector({ "float", "float32" }) },
      { PlyType::NNP_FLOAT64, NameVector({ "double", "float64" }) },
      { PlyType::NNP_INT8, NameVector({ "char", "int8" }) },
      { PlyType::NNP_INT16, NameVector({ "short", "int16" }) },
      { PlyType::NNP_INT32, NameVector({ "int", "int32" }) },
      { PlyType::NNP_UINT8, NameVector({ "uchar", "uint8" }) },
      { PlyType::NNP_UINT16, NameVector({ "ushort", "uint16" }) },
      { PlyType::NNP_UINT32, NameVector({ "uint", "uint32" }) },
      { PlyType::NNP_LIST_UINT8_UINT32, NameVector({ "list uchar uint", "list uint8 uint32" }) },
      { PlyType::NNP_LIST_INT8_UINT32, NameVector({ "list char uint", "list int8 uint32" }) },
      { PlyType::NNP_LIST_UINT8_INT32, NameVector({ "list uchar int", "list uint8 int32" }) },
      { PlyType::NNP_LIST_INT8_INT32, NameVector({ "list char int", "list int8 int32" }) },
      { PlyType::NNP_LIST_UINT8_FLOAT32, NameVector({ "list uchar float", "list uint8 float32" }) },
      { PlyType::NNP_LIST_INT8_FLOAT32, NameVector({ "list char float", "list int8 float32" }) },
      { PlyType::NNP_LIST_UINT8_FLOAT64, NameVector({ "list uchar double", "list uint8 float64" }) },
      { PlyType::NNP_LIST_INT8_FLOAT64, NameVector({ "list char double", "list int8 float64" }) },
      { PlyType::NNP_LIST_UINT8_UINT8, NameVector({ "list uchar uchar", "list uint8 uint8" }) },
      { PlyType::NNP_LIST_INT8_UINT8, NameVector({ "list char uchar", "list int8 uint8" }) },
      { PlyType::NNP_LIST_UINT8_INT8, NameVector({ "list uchar char", "list uint8 int8" }) },
      { PlyType::NNP_LIST_INT8_INT8, NameVector({ "list char char", "list int8 int8" }) },
      { PlyType::NNP_LIST_UINT8_UINT16, NameVector({ "list uchar ushort", "list uint8 uint16" }) },
      { PlyType::NNP_LIST_INT8_UINT16, NameVector({ "list char ushort", "list int8 uint16" }) },
      { PlyType::NNP_LIST_UINT8_INT16, NameVector({ "list uchar short", "list uint8 int16" }) },
      { PlyType::NNP_LIST_INT8_INT16, NameVector({ "list char short", "list int8 int16" }) }
  });


  /* Names used for the PlyProperty */
  static EntityMap mapProp({
      { PlyEntity::NNP_UNKNOWN_ENTITY, NameVector({ "unknonw" }) },
      { PlyEntity::NNP_PX, NameVector({ "x", "px", "posx" }) },
      { PlyEntity::NNP_PY, NameVector({ "y", "py", "posy" }) },
      { PlyEntity::NNP_PZ, NameVector({ "z", "pz", "posz" }) },
      { PlyEntity::NNP_PXYZ, NameVector({ "x y z", "px py pz", "posx posy posz" }) },
      { PlyEntity::NNP_NX, NameVector({ "nx", "normalx" }) },
      { PlyEntity::NNP_NY, NameVector({ "ny", "normaly" }) },
      { PlyEntity::NNP_NZ, NameVector({ "nz", "normalz" }) },
      { PlyEntity::NNP_NXYZ, NameVector({ "nx ny nz", "normalx normaly normalz" }) },
      { PlyEntity::NNP_CR, NameVector({ "red", "diffuse_red", "r", "diffuse_r" }) },
      { PlyEntity::NNP_CG, NameVector({ "green", "diffuse_green", "g", "diffuse_g" }) },
      { PlyEntity::NNP_CB, NameVector({ "blue", "diffuse_blue", "b", "diffuse_b" }) },
      { PlyEntity::NNP_CA, NameVector({ "alpha", "diffuse_alpha", "a", "diffuse_a" }) },
      { PlyEntity::NNP_CRGB, NameVector({ "rgb", "diffuse_rgb" }) },
      { PlyEntity::NNP_CRGBA, NameVector({ "rgba", "diffuse_rgba" }) },
      { PlyEntity::NNP_DENSITY, NameVector({ "radius", "density" }) },
      { PlyEntity::NNP_SCALE, NameVector({ "scale", "value" }) },
      { PlyEntity::NNP_TEXTUREU, NameVector({ "texture_u", "u", "s" }) },
      { PlyEntity::NNP_TEXTUREV, NameVector({ "texture_v", "v", "t" }) },
      { PlyEntity::NNP_TEXTURE2D, NameVector({ "texture_uv", "uv" }) },
      { PlyEntity::NNP_TEXTUREW, NameVector({ "texture_w", "w" }) },
      { PlyEntity::NNP_TEXTURE3D, NameVector({ "texture_uvw", "uvw" }) },
      { PlyEntity::NNP_TEXTUREINDEX, NameVector({ "texnumber", "texid" }) },
      { PlyEntity::NNP_QUALITY, NameVector({ "quality", "confidence" }) },
      { PlyEntity::NNP_REFLECTANCE, NameVector({ "reflectance" }) },
      { PlyEntity::NNP_BITFLAG, NameVector({ "flags" }) },
      { PlyEntity::NNP_K1, NameVector({ "k1" }) },
      { PlyEntity::NNP_K2, NameVector({ "k2" }) },
      { PlyEntity::NNP_KG, NameVector({ "k" }) },
      { PlyEntity::NNP_KH, NameVector({ "h" }) },
      { PlyEntity::NNP_K1DIR, NameVector({ "k1dir" }) },
      { PlyEntity::NNP_K2DIR, NameVector({ "k2dir" }) },
      { PlyEntity::NNP_EDGE_V1, NameVector({ "vertex1", "v1" }) },
      { PlyEntity::NNP_EDGE_V2, NameVector({ "vertex2", "v2" }) },
      { PlyEntity::NNP_FACE_VERTEX_LIST, NameVector({ "vertex_index", "vertex_indices" }) },
      { PlyEntity::NNP_FACE_WEDGE_COLOR, NameVector({ "color" }) },
      { PlyEntity::NNP_FACE_WEDGE_NORMAL, NameVector({ "normal" }) },
      { PlyEntity::NNP_FACE_WEDGE_TEX, NameVector({ "texcoord" }) }
  });


  /* Names used for the PlyElement */
  static ElementMap mapElem({
      { PlyElemEntity::NNP_UNKNOWN_ELEM, NameVector({ "unknonw" }) },
      { PlyElemEntity::NNP_VERTEX_ELEM, NameVector({ "vertex" }) },
      { PlyElemEntity::NNP_EDGE_ELEM, NameVector({ "edge" }) },
      { PlyElemEntity::NNP_FACE_ELEM, NameVector({ "face" }) },
  });


  /* Returns the vector of possible name for the input PlyEntity */
  static inline const NameVector& PlyPropertyName(PlyEntity ent)
  {
    static NameVector emptyVec;
    if (mapProp.find(ent) != mapProp.end())
      return 	mapProp[ent];
    return  emptyVec;
  };


  /* Returns the vector of possible name for the input PlyType */
  static inline const NameVector& PlyTypeName(PlyType ent)
  {
    static NameVector emptyVec;
    if (mapType.find(ent) != mapType.end())
      return 	mapType[ent];
    return  emptyVec;
  };


  /* Returns the vector of possible name for the input PlyElement */
  static inline const NameVector& PlyElementName(PlyElemEntity ent)
  {
    static NameVector emptyVec;
    if (mapElem.find(ent) != mapElem.end())
      return 	mapElem[ent];
    return  emptyVec;
  };


  /* Return 1 for little endian, 0 for big endian*/
  inline int checkEndianness()
  {
    unsigned int x = 1;
    char *c = (char*)&x;
    return (int)*c;
  };

  /* Adjust the endianess of the input buffer*/
  inline void adjustEndianess(unsigned char* buffer, int typeSize, int count)
  {
    for (int i = 0; i < count; i++)
    {
      int offset = i*typeSize;
      for (int j = 0; j < typeSize / 2; j++)
      {
        unsigned char temp = buffer[offset + j];
        buffer[offset + j] = buffer[offset + typeSize - 1 - j];
        buffer[offset + typeSize - 1 - j] = temp;
      }
    }
  }
  /**
  * @endcond
  */



  /** PLY File.
  *  Class to manage the read and write of a PLY file using a memory buffer.
  */
  class PlyFile
  {
  private:
    std::fstream fileStream;	/**< Stream. */
    char mode;					/**< Mode of the stream (0 = read, 1 = write). */

    int64_t bufferSize;			/**< Size of the buffer. */
    int64_t bufferOffset;		/**< Offset for the next operation. */
    int64_t maxSize;			/**< Maximum size of the stream. */
    char* buffer;				/**< Buffer. */

  public:

    PlyFile();

    ~PlyFile();

    /**
    * Open the file in read mode.
    *
    * @param filename	name of the file.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    inline bool OpenFileToRead(const std::string &filename);

    /**
    * Open the file in write mode.
    *
    * @param filename	name of the file.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    inline bool OpenFileToWrite(const std::string &filename);

    /**
    * Read the next header line.
    *
    * @param line	read line.
    * @param last	true if the line is the last of header.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool NextHeaderLine(std::string &line, bool &last);

    /**
    * Write the line in the header.
    *
    * @param line	line to write.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool WriteHeaderLine(const std::string &line);

    /**
    * Read binary data from the file.
    *
    * @param dest	pointer where to copy the data.
    * @param nByte	number of byte to read.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool ReadBinaryData(char * & dest, int nByte);

    /**
    * Read ASCII data from the file.
    *
    * @param dest	reference to the container of the data.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    template <typename T>
    inline bool ReadAsciiData(T &dest);

    /**
    * Write binary data in the file.
    *
    * @param src	pointer to the data to write.
    * @param nByte	number of byte to write.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool WriteBinaryData(void *src, int nByte);

    /**
    * Write ASCII data in the file.
    *
    * @param src	reference to the container of the data.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    template <typename T>
    inline bool WriteAsciiData(const T &src);

    /**
    * Set the maximum size of the buffer.
    *
    * @return size	size of the buffer.
    */
    inline void SetBufferSize(int64_t size);

    /**
    * Force the write of the buffer in the file.
    */
    inline void Flush();
  };


  inline PlyFile::PlyFile()
  {
    buffer = NULL;
    maxSize = 10 * (1 << 20);
    mode = -1;
  }


  inline PlyFile::~PlyFile()
  {
    Flush();
    if (buffer != NULL)
      delete[] buffer;
    if (fileStream.is_open())
      fileStream.close();
  }


  inline bool PlyFile::OpenFileToRead(const std::string& filename)
  {
    if (fileStream.is_open())
      fileStream.close();
    mode = 0;
    fileStream.open(filename, std::fstream::in | std::fstream::binary);
    if (fileStream.fail())
      return false;
    bufferOffset = 0;
    return true;
  }


  inline bool PlyFile::OpenFileToWrite(const std::string& filename)
  {
    if (fileStream.is_open())
      fileStream.close();
    mode = 1;
    fileStream.open(filename, std::fstream::out | std::fstream::binary);
    if (fileStream.fail())
      return false;
    bufferOffset = 0;
    //fileStream.setf(std::ios::fixed, std::ios::floatfield);
    //fileStream.precision(7);
    return true;
  }


  inline bool PlyFile::NextHeaderLine(std::string& line, bool& last)
  {
    if (mode != 0)
      return false;
    std::getline(fileStream, line);
    std::transform(line.begin(), line.end(), line.begin(), ::tolower);
    last = false;
    if (line == "end_header")
      last = true;
    return true;
  }


  inline bool PlyFile::WriteHeaderLine(const std::string& line)
  {
    if (mode != 1)
      return false;
    fileStream << line;
    return true;
  }


  inline bool PlyFile::ReadBinaryData(char * & dest, int nByte)
  {
    if (mode != 0)
      return false;
    if (buffer == NULL)
    {
      buffer = new char[maxSize];
      fileStream.read(buffer, maxSize);
      bufferSize = fileStream.gcount();
      bufferOffset = 0;
    }
    else if (bufferOffset + nByte > bufferSize)
    {
      const size_t lastByte = bufferSize - bufferOffset;
      memcpy(buffer, &buffer[bufferOffset], lastByte);
      fileStream.read(&buffer[lastByte], maxSize - lastByte);
      bufferSize = fileStream.gcount() + lastByte;
      bufferOffset = 0;
    }
    //memcpy(dest, &buffer[bufferOffset], nByte);
		dest = buffer + bufferOffset;
    bufferOffset += nByte;
    return true;
  }


  template <typename T>
  inline bool PlyFile::ReadAsciiData(T& dest)
  {
    if (mode != 0)
      return false;
    fileStream >> dest;
    return true;
  }

  inline bool PlyFile::WriteBinaryData(void* src, int nByte)
  {
    if (mode != 1)
      return false;
    if (buffer == NULL)
    {
      buffer = new char[maxSize];
      bufferSize = maxSize;
      bufferOffset = 0;
    }
    else if (bufferOffset + nByte > bufferSize)
    {
      fileStream.write(buffer, bufferOffset);
      bufferOffset = 0;
    }
    memcpy(&buffer[bufferOffset], src, nByte);
    bufferOffset += nByte;
    return true;
  }


  template <typename T>
  inline bool PlyFile::WriteAsciiData(const T& src)
  {
    if (mode != 1)
      return false;
    fileStream << src;
    return true;
  }

  inline void PlyFile::SetBufferSize(int64_t size)
  {
    maxSize = size;
  }


  inline void PlyFile::Flush()
  {
    if (mode == 1)
      fileStream.write(buffer, bufferOffset);
  }




  /** PLY Property.
  *  Define a PLY property (entity and type).
  */
  class PlyProperty
  {
  public:

    std::string name;	/**< Property name. */
    PlyType type;		/**< Property type. */
    PlyEntity elem;		/**< Property entity. */
    bool validToWrite;  /**< Property validity (necessary to write the header). */

    /**
    * Constructor that sets the type and the entity of a standard PLY property.
    *
    * @param _t	Property type.
    * @param _e	Property entity.
    */
    inline PlyProperty(PlyType _t, PlyEntity _e) :type(_t), elem(_e), name(PlyPropertyName(_e)[0]), validToWrite(false){}

    /**
    * Constructor that sets the type, the entity and the name of a standard PLY property.
    *
    * @param _t		Property type.
    * @param _e		Property entity.
    * @param _n		Property name.
    */
    inline PlyProperty(PlyType _t, PlyEntity _e, std::string _n) :type(_t), elem(_e), name(_n), validToWrite(false){}

    /**
    * Constructor that sets the type and the name of a custom PLY property.
    *
    * @param _t		Property type.
    * @param _n		Property name.
    */
    inline PlyProperty(PlyType _t, std::string _n) :type(_t), elem(PlyEntity::NNP_UNKNOWN_ENTITY), name(_n), validToWrite(false){}

    /**
    * Get the description string of the property entity.
    *
    * @return	Description string of the property entity.
    */
    inline const char* EntityStr();

    /**
    * Get the name of the property entity.
    *
    * @return	Name of the property entity.
    */
    inline const char* EntityName();

    /**
    * Get the description string of the property type.
    *
    * @return	Description string of the property type.
    */
    inline const char* TypeStr();

    /**
    * Get the size in byte of the property type.
    *
    * @return	Size in byte of the property type.
    */
    inline int TypeSize();

    /**
    * Get the number of component of the property entity.
    *
    * @return	Number of component.
    */
    inline int CountValue();

    /**
    * Check if the property type is signed or unsigned.
    *
    * @return	true = signed, false = unsigned.
    */
    inline bool IsSigned();

    /**
    * Skip the property in an Ascii file.
    *
    * @param file	Opened file.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool SkipAsciiPropertyInFile(PlyFile &file);

    /**
    * Skip the property in a binary file.
    *
    * @param file	Opened file.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool SkipBinaryPropertyInFile(PlyFile &file);

    /**
    * Write the property string in the header of the PLY file.
    *
    * @param file	Opened file.
    * @return		  If successful returns true. Otherwise, it returns false.
    */
    inline bool WriteHeader(PlyFile &file);
  };


  inline const char* PlyProperty::EntityStr()
  {
    switch (this->elem)
    {
    case NNP_UNKNOWN_ENTITY:    return "NNP_UNKNOWN_ENTITY   ";
    case NNP_PX:    return "NNP_PX               ";
    case NNP_PY:    return "NNP_PY               ";
    case NNP_PZ:    return "NNP_PZ               ";
    case NNP_PXYZ:    return "NNP_PXYZ             ";
    case NNP_NX:    return "NNP_NX               ";
    case NNP_NY:    return "NNP_NY               ";
    case NNP_NZ:    return "NNP_NZ               ";
    case NNP_NXYZ:    return "NNP_NXYZ             ";
    case NNP_CR:    return "NNP_CR               ";
    case NNP_CG:    return "NNP_CG               ";
    case NNP_CB:    return "NNP_CB               ";
    case NNP_CRGB:    return "NNP_CRGB             ";
    case NNP_CA:    return "NNP_CA               ";
    case NNP_CRGBA:    return "NNP_CRGBA            ";
    case NNP_DENSITY:    return "NNP_DENSITY          ";
    case NNP_SCALE:    return "NNP_SCALE            ";
    case NNP_QUALITY:    return "NNP_QUALITY          ";
    case NNP_REFLECTANCE:    return "NNP_REFLECTANCE      ";
    case NNP_TEXTUREU:    return "NNP_TEXTUREU         ";
    case NNP_TEXTUREV:    return "NNP_TEXTUREV         ";
    case NNP_TEXTURE2D:    return "NNP_TEXTURE2D        ";
    case NNP_TEXTUREW:    return "NNP_TEXTUREW         ";
    case NNP_TEXTURE3D:    return "NNP_TEXTURE3D        ";
    case NNP_TEXTUREINDEX:    return "NNP_TEXTUREINDEX     ";
    case NNP_BITFLAG:    return "NNP_BITFLAG          ";
    case NNP_K1:    return "NNP_K1               ";
    case NNP_K2:    return "NNP_K2               ";
    case NNP_KG:    return "NNP_K                ";
    case NNP_KH:    return "NNP_H                ";
    case NNP_K1DIR:    return "NNP_K1DIR            ";
    case NNP_K2DIR:    return "NNP_K2DIR            ";
    case NNP_EDGE_V1:    return "NNP_EDGE_V1          ";
    case NNP_EDGE_V2:    return "NNP_EDGE_V2          ";
    case NNP_FACE_VERTEX_LIST:    return "NNP_FACE_VERTEX_LIST ";
    case NNP_FACE_WEDGE_COLOR:    return "NNP_FACE_WEDGE_COLOR ";
    case NNP_FACE_WEDGE_NORMAL:    return "NNP_FACE_WEDGE_NORMAL";
    case NNP_FACE_WEDGE_TEX:    return "NNP_FACE_WEDGE_TEX   ";

    default: assert(0);
      break;
    }
    return 0;
  }


  inline const char* PlyProperty::EntityName()
  {

    if (this->elem == PlyEntity::NNP_UNKNOWN_ENTITY)
      return (*this).name.c_str();
#ifdef USE_NOSTANDARDPLY_OUTPUT
    if (this->elem == PlyEntity::NNP_SCALE)
      return PlyPropertyName(this->elem)[1].c_str();
    if (this->elem == PlyEntity::NNP_QUALITY)
      return PlyPropertyName(this->elem)[1].c_str();
    if (this->elem == PlyEntity::NNP_CR)
      return PlyPropertyName(this->elem)[1].c_str();
    if (this->elem == PlyEntity::NNP_CB)
      return PlyPropertyName(this->elem)[1].c_str();
    if (this->elem == PlyEntity::NNP_CG)
      return PlyPropertyName(this->elem)[1].c_str();
    if (this->elem == PlyEntity::NNP_CA)
      return PlyPropertyName(this->elem)[1].c_str();
#endif
    const std::vector<std::string>& vect = PlyPropertyName(this->elem);
    if (vect.size() > 0)
      return vect[0].c_str();
    assert(0);
    return 0;
  }


  inline const char* PlyProperty::TypeStr()
  {
    switch (this->type)
    {
    case NNP_UNKNOWN_TYPE:    return "NNP_UNKNOWN_TYPE       ";
    case NNP_FLOAT32:    return "NNP_FLOAT32            ";
    case NNP_FLOAT64:    return "NNP_FLOAT64            ";
    case NNP_INT8:    return "NNP_INT8               ";
    case NNP_INT16:    return "NNP_INT16              ";
    case NNP_INT32:    return "NNP_INT32              ";
    case NNP_UINT8:    return "NNP_UINT8              ";
    case NNP_UINT16:    return "NNP_UINT16             ";
    case NNP_UINT32:    return "NNP_UINT32             ";
    case NNP_LIST_UINT8_UINT32:    return "NNP_LIST_UINT8_UINT32  ";
    case NNP_LIST_INT8_UINT32:    return "NNP_LIST_INT8_UINT32   ";
    case NNP_LIST_UINT8_INT32:    return "NNP_LIST_UINT8_INT32   ";
    case NNP_LIST_INT8_INT32:    return "NNP_LIST_INT8_INT32    ";
    case NNP_LIST_UINT8_FLOAT32:    return "NNP_LIST_UINT8_FLOAT32 ";
    case NNP_LIST_INT8_FLOAT32:    return "NNP_LIST_INT8_FLOAT32  ";
    case NNP_LIST_UINT8_FLOAT64:    return "NNP_LIST_UINT8_FLOAT64 ";
    case NNP_LIST_INT8_FLOAT64:    return "NNP_LIST_INT8_FLOAT64  ";
    case NNP_LIST_UINT8_UINT8:    return "NNP_LIST_UINT8_UINT8  ";
    case NNP_LIST_INT8_UINT8:    return "NNP_LIST_INT8_UINT8   ";
    case NNP_LIST_UINT8_INT8:    return "NNP_LIST_UINT8_INT8   ";
    case NNP_LIST_INT8_INT8:    return "NNP_LIST_INT8_INT8    ";
    case NNP_LIST_UINT8_UINT16:    return "NNP_LIST_UINT8_UINT16  ";
    case NNP_LIST_INT8_UINT16:    return "NNP_LIST_INT8_UINT16   ";
    case NNP_LIST_UINT8_INT16:    return "NNP_LIST_UINT8_INT16   ";
    case NNP_LIST_INT8_INT16:    return "NNP_LIST_INT8_INT16    ";
    default: assert(0);
      break;
    }
    return 0;
  }


  inline int PlyProperty::TypeSize()
  {
    switch (this->type)
    {
    case NNP_UNKNOWN_TYPE:
      return 0;
    case NNP_INT8:
    case NNP_UINT8:
      return sizeof(char);
    case NNP_INT16:
    case NNP_UINT16:
      return sizeof(short);
    case NNP_FLOAT32:
    case NNP_INT32:
    case NNP_UINT32:
      return sizeof(int);
    case NNP_FLOAT64:
      return sizeof(double);
    case NNP_LIST_UINT8_UINT32:
    case NNP_LIST_INT8_UINT32:
    case NNP_LIST_UINT8_INT32:
    case NNP_LIST_INT8_INT32:
      return sizeof(int);
    case NNP_LIST_UINT8_UINT16:
    case NNP_LIST_INT8_UINT16:
    case NNP_LIST_UINT8_INT16:
    case NNP_LIST_INT8_INT16:
      return sizeof(short);
    case NNP_LIST_UINT8_UINT8:
    case NNP_LIST_INT8_UINT8:
    case NNP_LIST_UINT8_INT8:
    case NNP_LIST_INT8_INT8:
      return sizeof(char);
    case NNP_LIST_UINT8_FLOAT32:
    case NNP_LIST_INT8_FLOAT32:
      return sizeof(float);
    case NNP_LIST_UINT8_FLOAT64:
    case NNP_LIST_INT8_FLOAT64:
      return sizeof(double);
    default: assert(0);
      break;
    }
    return 0;
  }


  inline bool PlyProperty::IsSigned()
  {
    switch (this->type)
    {
    case NNP_INT8:
    case NNP_INT16:
    case NNP_INT32:
    case NNP_FLOAT32:
    case NNP_FLOAT64:
    case NNP_LIST_INT8_UINT32:
    case NNP_LIST_INT8_INT32:
    case NNP_LIST_INT8_UINT16:
    case NNP_LIST_INT8_INT16:
    case NNP_LIST_INT8_UINT8:
    case NNP_LIST_INT8_INT8:
    case NNP_LIST_INT8_FLOAT32:
    case NNP_LIST_INT8_FLOAT64:
      return true;
    case NNP_UINT8:
    case NNP_UINT16:
    case NNP_UINT32:
    case NNP_LIST_UINT8_UINT32:
    case NNP_LIST_UINT8_INT32:
    case NNP_LIST_UINT8_UINT16:
    case NNP_LIST_UINT8_INT16:
    case NNP_LIST_UINT8_UINT8:
    case NNP_LIST_UINT8_INT8:
    case NNP_LIST_UINT8_FLOAT32:
    case NNP_LIST_UINT8_FLOAT64:
      return false;
    default:
      return false;
    }
  }


  inline int PlyProperty::CountValue()
  {
    if (this->elem == NNP_CRGB || this->elem == NNP_NXYZ || this->elem == NNP_PXYZ || this->elem == NNP_TEXTURE3D)
      return 3;
    else if (this->elem == NNP_CRGBA)
      return 4;
    else if (this->elem == NNP_TEXTURE2D)
      return 2;
    return 1;
  }


  inline bool PlyProperty::SkipAsciiPropertyInFile(PlyFile &file)
  {
    int count = CountValue();
    if (this->type >= NNP_LIST_UINT8_UINT32)
      file.ReadAsciiData(count);
    switch (type)
    {
    case NNP_INT8:
    case NNP_INT16:
    case NNP_INT32:
    case NNP_LIST_UINT8_INT32:
    case NNP_LIST_INT8_INT32:
    case NNP_LIST_UINT8_INT16:
    case NNP_LIST_INT8_INT16:
    case NNP_LIST_UINT8_INT8:
    case NNP_LIST_INT8_INT8:
    {
      int* temp = new int[count];
      for (int i = 0; i < count; i++)
        file.ReadAsciiData(temp[i]);
      delete[] temp;
      break;
    }
    case NNP_UINT8:
    case NNP_UINT16:
    case NNP_UINT32:
    case NNP_LIST_UINT8_UINT32:
    case NNP_LIST_INT8_UINT32:
    case NNP_LIST_UINT8_UINT16:
    case NNP_LIST_INT8_UINT16:
    case NNP_LIST_UINT8_UINT8:
    case NNP_LIST_INT8_UINT8:
    {
      unsigned int* temp = new unsigned int[count];
      for (int i = 0; i < count; i++)
        file.ReadAsciiData(temp[i]);
      delete[] temp;
      break;
    }
    case NNP_FLOAT32:
    case NNP_LIST_UINT8_FLOAT32:
    case NNP_LIST_INT8_FLOAT32:
    {
      float* temp = new float[count];
      for (int i = 0; i < count; i++)
        file.ReadAsciiData(temp[i]);
      delete[] temp;
      break;
    }
    case NNP_FLOAT64:
    case NNP_LIST_UINT8_FLOAT64:
    case NNP_LIST_INT8_FLOAT64:
    {
      double* temp = new double[count];
      for (int i = 0; i < count; i++)
        file.ReadAsciiData(temp[i]);
      delete[] temp;
      break;
    }
    }
    return true;
  }


  inline bool PlyProperty::SkipBinaryPropertyInFile(PlyFile& file)
  {
    char * temp = nullptr;
    int count = CountValue();
    if (this->type >= NNP_LIST_UINT8_UINT32)
    {
      int size;
      if (this->IsSigned())
      {
        file.ReadBinaryData(temp, sizeof(char));
        size = this->TypeSize() * int(*(reinterpret_cast<unsigned char *>(temp)));
      }
      else
      {
        file.ReadBinaryData(temp, sizeof(char));
        size = this->TypeSize() * int(*(reinterpret_cast<unsigned char *>(temp)));
      }
      file.ReadBinaryData(temp, size);
    }
    else
    {
      int size = this->TypeSize() * count;
      file.ReadBinaryData(temp, size);
    }
    return true;
  }


  inline bool PlyProperty::WriteHeader(PlyFile& file)
  {
    if (!validToWrite)
      return true;
    std::string name, type;
    type = PlyTypeName(this->type)[0];
    name = this->EntityName();
    std::vector<std::string> v;
    switch (this->elem)
    {
    case NNP_PXYZ:
    {
      v.push_back(PlyPropertyName(NNP_PX)[0]);
      v.push_back(PlyPropertyName(NNP_PY)[0]);
      v.push_back(PlyPropertyName(NNP_PZ)[0]);
      break;
    }
    case NNP_NXYZ:
    {
      v.push_back(PlyPropertyName(NNP_NX)[0]);
      v.push_back(PlyPropertyName(NNP_NY)[0]);
      v.push_back(PlyPropertyName(NNP_NZ)[0]);
      break;
    }
    case NNP_CRGB:
    {
      v.push_back(PlyPropertyName(NNP_CR)[0]);
      v.push_back(PlyPropertyName(NNP_CG)[0]);
      v.push_back(PlyPropertyName(NNP_CB)[0]);
      break;
    }
    case NNP_CRGBA:
    {
      v.push_back(PlyPropertyName(NNP_CR)[0]);
      v.push_back(PlyPropertyName(NNP_CG)[0]);
      v.push_back(PlyPropertyName(NNP_CB)[0]);
      v.push_back(PlyPropertyName(NNP_CA)[0]);
      break;
    }
    case NNP_TEXTURE2D:
    {
      v.push_back(PlyPropertyName(NNP_TEXTUREU)[0]);
      v.push_back(PlyPropertyName(NNP_TEXTUREV)[0]);
      break;
    }
    case NNP_TEXTURE3D:
    {
      v.push_back(PlyPropertyName(NNP_TEXTUREU)[0]);
      v.push_back(PlyPropertyName(NNP_TEXTUREV)[0]);
      v.push_back(PlyPropertyName(NNP_TEXTUREW)[0]);
      break;
    }
    default:
      v.push_back(name);

    }
    for (int i = 0; i < v.size(); i++)
    {
      std::stringstream s;
      s << "property " << type << " " << v[i] << "\n";
      if (!file.WriteHeaderLine(s.str()))
        return false;
    }
    return true;
  }




  /** PLY Element.
  *  Define a PLY Element as a collection of properties.
  */
  class PlyElement
  {
  public:

    std::string name;					/**< Name of the elment in the PLY header (for example "vertex", "face", ect.). */
    PlyElemEntity plyElem;				/**< Ply element entity. */
    size_t cnt;							/**< Number of instances of the elment in the PLY file. */
    std::vector<PlyProperty> propVec;   /**< Collection of properties that define the element. */
    bool validToWrite;					/**< Element validity (necessary to write the header). */

    /**
    * Default Constructor
    */
    inline PlyElement() :validToWrite(false){};

    /**
    * Constructor that sets the name, the properties and the number of instances of the element.
    *
    * @param _name		Name of the element.
    * @param prop		Vector of properties.
    * @param nElem		Number of instances.
    */
    inline PlyElement(std::string& _name, std::vector<PlyProperty> &prop, size_t nElem) :name(_name), cnt(nElem), propVec(prop), plyElem(PlyElemEntity::NNP_UNKNOWN_ELEM), validToWrite(false){};

    /**
    * Constructor that sets the entity, the properties and the number of instances of the element.
    *
    * @param ent		Element entity.
    * @param prop		Vector of properties.
    * @param nElem		Number of instances.
    */
    inline PlyElement(PlyElemEntity ent, std::vector<PlyProperty> &prop, size_t nElem) :name(PlyElementName(ent)[0]), cnt(nElem), propVec(prop), plyElem(ent), validToWrite(false){};

    /**
    * Parse the input line and add the properties to the element.
    * It assumes that the passed line has the folowing structure: "property PLYTYPE PLYELEMENT"
    *
    * @param line		Input line.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    inline bool AddProperty(std::string &line);

    /**
    * Initialize an element from the header line.
    *
    * @param elemStr	String with the element definition.
    * @param propStr	Strings with the property definition.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    inline bool InitFromHeader(std::string &elemStr, std::vector<std::string> &propStr);

    /**
    * Write the element descriport in the file header.
    *
    * @param file	Ply file.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool WriteHeader(PlyFile &file);

    /**
    * Skip the element in an Ascii file.
    *
    * @param file	Ply file
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool SkipAsciiElementsInFile(PlyFile &file);

    /**
    * Skip the element in a binary file.
    *
    * @param file	Ply file.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool SkipBinaryElementsInFile(PlyFile &file);

    /**
    * Check if the input entity is in the property of the element.
    *
    * @param entity		Input entity.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    inline bool Contains(NNP_ENTITY entity);
  };


  inline bool PlyElement::InitFromHeader(std::string &elemStr, std::vector<std::string> &propStr)
  {
    char* token;
    char* tempStr = &elemStr[0];
    token = strtok(tempStr, " \t");
    if (strstr(token, "element") == NULL)
      return false;
    token = strtok(0, " \t\n");
    name = std::string(token);
    plyElem = PlyElemEntity::NNP_UNKNOWN_ELEM;
    ElementMapIterator iter = mapElem.begin();
    bool found = false;
    while (iter != mapElem.end())
    {
      NameVector& v = (*iter).second;
      for (size_t i = 0; i < v.size(); i++)
      {
        if (v[i] == name)
        {
          found = true;
          break;
        }
      }
      if (found)
      {
        plyElem = (*iter).first;
        break;
      }
      iter++;
    }
    token = strtok(0, " \t\n");
    cnt = atoi(token);
    for (size_t i = 0; i < propStr.size(); i++)
      if (!AddProperty(propStr[i]))
        return false;
    unsigned int mask = 0;
    for (size_t i = 0; i < propVec.size(); i++)
      mask |= propVec[i].elem;
    std::vector<PlyProperty> compactPropVec;
    for (size_t i = 0; i < propVec.size(); i++)
    {
      switch (propVec[i].elem)
      {
      case NNP_NX:
      case NNP_NY:
        if ((mask & NNP_NXYZ) != NNP_NXYZ)
          compactPropVec.push_back(propVec[i]);
        break;
      case NNP_NZ:
        if ((mask & NNP_NXYZ) != NNP_NXYZ)
          compactPropVec.push_back(propVec[i]);
        else
          compactPropVec.push_back(PlyProperty(propVec[i].type, NNP_NXYZ));
        break;
      case NNP_PX:
      case NNP_PY:
        if ((mask & NNP_PXYZ) != NNP_PXYZ)
          compactPropVec.push_back(propVec[i]);
        break;
      case NNP_PZ:
        if ((mask & NNP_PXYZ) != NNP_PXYZ)
          compactPropVec.push_back(propVec[i]);
        else
          compactPropVec.push_back(PlyProperty(propVec[i].type, NNP_PXYZ));
        break;
      case NNP_CR:
      case NNP_CG:
        if (((mask & NNP_CRGB) != NNP_CRGB) && ((mask & NNP_CRGBA) != NNP_CRGBA))
          compactPropVec.push_back(propVec[i]);
        break;
      case NNP_CB:
        if (((mask & NNP_CRGB) != NNP_CRGB) & ((mask & NNP_CRGBA) != NNP_CRGBA))
          compactPropVec.push_back(propVec[i]);
        else if (((mask & NNP_CRGB) == NNP_CRGB) & ((mask & NNP_CRGBA) != NNP_CRGBA))
          compactPropVec.push_back(PlyProperty(propVec[i].type, NNP_CRGB));
        break;
      case NNP_CA:
        if (((mask & NNP_CRGB) != NNP_CRGB) & ((mask & NNP_CRGBA) != NNP_CRGBA))
          compactPropVec.push_back(propVec[i]);
        else if ((mask & NNP_CRGBA) == NNP_CRGBA)
          compactPropVec.push_back(PlyProperty(propVec[i].type, NNP_CRGBA));
        break;
      case NNP_TEXTUREU:
        if (((mask & NNP_TEXTURE2D) != NNP_TEXTURE2D) & ((mask & NNP_TEXTURE3D) != NNP_TEXTURE3D))
          compactPropVec.push_back(propVec[i]);
        break;
      case NNP_TEXTUREV:
        if (((mask & NNP_TEXTURE2D) != NNP_TEXTURE2D) & ((mask & NNP_TEXTURE3D) != NNP_TEXTURE3D))
          compactPropVec.push_back(propVec[i]);
        else if (((mask & NNP_TEXTURE2D) == NNP_TEXTURE2D) & ((mask & NNP_TEXTURE3D) != NNP_TEXTURE3D))
          compactPropVec.push_back(PlyProperty(propVec[i].type, NNP_TEXTURE2D));
        break;
      case NNP_TEXTUREW:
        if (((mask & NNP_TEXTURE2D) != NNP_TEXTURE2D) & ((mask & NNP_TEXTURE3D) != NNP_TEXTURE3D))
          compactPropVec.push_back(propVec[i]);
        else if ((mask & NNP_TEXTURE3D) == NNP_TEXTURE3D)
          compactPropVec.push_back(PlyProperty(propVec[i].type, NNP_TEXTURE3D));
        break;
      default:
        compactPropVec.push_back(propVec[i]);
        break;
      }
    }
    propVec.clear();
    propVec = compactPropVec;
    return true;
  }


  inline bool PlyElement::WriteHeader(PlyFile &file)
  {
    if (!validToWrite || cnt == 0)
      return true;
    bool ok = true;
    std::stringstream temp;
    temp << "element " << name << " " << cnt << "\n";
    if (file.WriteHeaderLine(temp.str()))
    {
      for (int i = 0; i < propVec.size(); i++)
        ok = propVec[i].WriteHeader(file);
    }
    else
      return false;
    return ok;
  }


  inline bool PlyElement::SkipAsciiElementsInFile(PlyFile &file)
  {
    for (int i = 0; i < this->cnt; ++i)
      for (int j = 0; j < this->propVec.size(); ++j)
        this->propVec[j].SkipAsciiPropertyInFile(file);
    return true;
  }


  inline bool PlyElement::SkipBinaryElementsInFile(PlyFile &file)
  {
    for (int i = 0; i < this->cnt; ++i)
      for (int j = 0; j < this->propVec.size(); ++j)
        this->propVec[j].SkipBinaryPropertyInFile(file);
    return true;
  }


  inline bool PlyElement::AddProperty(std::string &line)
  {
    char* token;
    char* tempStr = &line[0];
    token = strtok(tempStr, " \t");
    if (strstr(token, "property") == NULL)
      return false;
    char* typeStr = strtok(0, " \t\n");
    char *ty1, *ty2;
    std::string type;
    type.append(typeStr);
    if (strcmp(typeStr, "list") == 0)
    {
      ty1 = strtok(0, " \t\n");
      ty2 = strtok(0, " \t\n");
      type.append(" ");
      type.append(ty1);
      type.append(" ");
      type.append(ty2);
    }
    PlyType plyType = PlyType::NNP_UNKNOWN_TYPE;
    TypeMapIterator iterType = mapType.begin();
    bool found = false;
    while (iterType != mapType.end())
    {
      NameVector& v = (*iterType).second;
      for (size_t i = 0; i < v.size(); i++)
      {
        if (v[i] == type)
        {
          found = true;
          break;
        }
      }
      if (found)
      {
        plyType = (*iterType).first;
        break;
      }
      iterType++;
    }
    if (plyType == PlyType::NNP_UNKNOWN_TYPE)
      return false;
    char* nameStr = strtok(0, " \t\n");
    PlyEntity plyEntity = PlyEntity::NNP_UNKNOWN_ENTITY;
    EntityMapIterator iterEnt = mapProp.begin();
    found = false;
    while (iterEnt != mapProp.end())
    {
      NameVector& v = (*iterEnt).second;
      for (size_t i = 0; i < v.size(); i++)
      {
        if (v[i] == nameStr)
        {
          found = true;
          break;
        }
      }
      if (found)
      {
        plyEntity = (*iterEnt).first;
        break;
      }
      iterEnt++;
    }
    if (plyEntity != PlyEntity::NNP_UNKNOWN_ENTITY)
      propVec.push_back(PlyProperty(plyType, plyEntity, nameStr));
    else
      propVec.push_back(PlyProperty(plyType, nameStr));
    return true;
  }


  inline bool PlyElement::Contains(PlyEntity entity)
  {
    for (int i = 0; i < propVec.size(); i++)
    {
      if (propVec[i].elem == entity)
        return true;
    }
    return false;
  }




  /** PLY header info.
  *  Define the data of the PLY header
  */
  class Info
  {

  public:
    ErrorCode errInfo;						/**< Error code returned by the reading of a PLY file. */
    bool binary;							/**< Boolean about the file format (Binary = true, ASCII = false). */
    std::vector<PlyElement> elemVec;		/**< Elements defined in the header. */
    bool bigEndian;							/**< Endianess of the binary file. */
    std::string filename;					/**< Filename. */
    std::vector<std::string> textureFile;	/**< Texture file names. */

    /**
    * Default Constructor
    */
    inline Info();

    /**
    * Constructor that reads the header info from a file.
    *
    * @param filename	Path of the file to read.
    */
    inline Info(const std::string& filename);

    /**
    * Load the ply info from the header of the input filename.
    *
    * @param filename	Path of the file to read.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    inline bool LoadHeader(const std::string& filename);

    /**
    * Write the ply info in the header of the input file.
    *
    * @param file	File to write.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool WriteHeader(PlyFile& file);

    /**
    * Add the ply element to the header.
    *
    * @param pe	  Ply element to write in the header.
    */
    inline void AddPlyElement(PlyElement& pe);

    /**
    * Clear the error code.
    */
    inline void Clear() { errInfo = NNP_OK; }

    /**
    * Return the number of instances of the element with the input name
    *
    * @param name	Name of the element.
    * @return		The number of instances
    */
    inline size_t GetElementCount(std::string& name);

    /**
    * Return the number of instances of the element with the input element type
    *
    * @param e	Element type.
    * @return	The number of instances
    */
    inline size_t GetElementCount(PlyElemEntity e);

    /**
    * Return the number of vertex instances
    *
    * @return The number of vertex instances
    */
    inline size_t GetVertexCount();

    /**
    * Return the number of face instances
    *
    * @return The number of face instances
    */
    inline size_t GetFaceCount();

    /**
    * Return the number of edge instances
    *
    * @return The number of edge instances
    */
    inline size_t GetEdgeCount();

    /**
    * Return a reference to the element with a specific name
    *
    * @param name	Name of the element.
    * @return		The reference to the element
    */
    inline PlyElement* GetElement(std::string& name);

    /**
    * Return a reference to the element with a specific element type
    *
    * @param e	Element type.
    * @return	The reference to the element
    */
    inline PlyElement* GetElement(PlyElemEntity e);

    /**
    * Return a reference to the vertex element
    *
    * @return The reference to the vertex element
    */
    inline PlyElement* GetVertexElement();

    /**
    * Return a reference to the face element
    *
    * @return The reference to the face element
    */
    inline PlyElement* GetFaceElement();

    /**
    * Return a reference to the edge element
    *
    * @return The reference to the edge element
    */
    inline PlyElement* GetEdgeElement();

  };


  inline Info::Info()
  {
    this->binary = true;
    this->bigEndian = true;
    this->Clear();
  }


  inline Info::Info(const std::string& filename)
  {
    this->LoadHeader(filename);
  }


  inline bool Info::LoadHeader(const std::string& filename)
  {
    this->filename = filename;
    this->errInfo = NNP_OK;
    std::ifstream input(filename, std::ios::binary);
    if (!input.good())
    {
      this->errInfo = NNP_UNABLE_TO_OPEN;
      input.close();
      return false;
    }
    std::string buffer;
    std::getline(input, buffer);
    std::transform(buffer.begin(), buffer.end(), buffer.begin(), ::tolower);
    if (buffer != "ply")
    {
      this->errInfo = NNP_MISSING_HEADER;
      input.close();
      return false;
    }
    std::getline(input, buffer);
    std::transform(buffer.begin(), buffer.end(), buffer.begin(), ::tolower);
    std::size_t pos = buffer.find("format");
    if (pos == std::string::npos)
    {
      this->errInfo = NNP_MISSING_FORMAT;
      input.close();
      return false;
    }
    if (buffer.find("binary_lit") != std::string::npos)
    {
      this->binary = true;
      this->bigEndian = false;
    }
    else if (buffer.find("binary_big") != std::string::npos)
    {
      this->binary = true;
      this->bigEndian = true;
    }
    else if (buffer.find("ascii") != std::string::npos)
    {
      this->binary = false;
      this->bigEndian = false;
    }
    else
    {
      this->errInfo = NNP_MISSING_FORMAT;
      return false;
    }
    std::getline(input, buffer);
		std::string lowBuffer;
		lowBuffer.resize(buffer.size());		
		std::transform(buffer.begin(), buffer.end(), lowBuffer.begin(), ::tolower);
    while (lowBuffer != "end_header")
    {
      if (lowBuffer.find("element") != std::string::npos)
      {
        std::string elemStr = lowBuffer;
        std::vector<std::string> propStr;
        do
        {
          std::getline(input, buffer);
					lowBuffer.clear();
					lowBuffer.resize(buffer.size());
          std::transform(buffer.begin(), buffer.end(), lowBuffer.begin(), ::tolower);
          pos = lowBuffer.find("property");
          if (pos != std::string::npos)
            propStr.push_back(lowBuffer);
        } while (pos != std::string::npos);
        PlyElement pe;
        if (!pe.InitFromHeader(elemStr, propStr))
        {
          this->errInfo = NNP_INVALID_ELEMENT;
          input.close();
          return false;
        }
        elemVec.push_back(pe);
      }
      else
      {
        if (lowBuffer.find("comment texture") != std::string::npos)
          textureFile.push_back(buffer.substr(buffer.find(" ", 10) + 1));
        std::getline(input, buffer);
				lowBuffer.clear();
				lowBuffer.resize(buffer.size());
        std::transform(buffer.begin(), buffer.end(), lowBuffer.begin(), ::tolower);
      }
    }
    input.close();
    return true;
  }


  inline bool Info::WriteHeader(PlyFile& file)
  {
    bool ok = true;
    ok = file.WriteHeaderLine(std::string("ply\n"));
    if (this->binary)
      ok = file.WriteHeaderLine(std::string("format binary_little_endian 1.0\n"));
    else
      ok = file.WriteHeaderLine(std::string("format ascii 1.0\n"));
    ok = file.WriteHeaderLine(std::string("comment nanoply generated\n"));
    for (int i = 0; i < this->textureFile.size(); i++)
      ok = file.WriteHeaderLine(std::string("comment TextureFile ") + this->textureFile[i] + "\n");
    for (int i = 0; i < this->elemVec.size(); i++)
      ok = this->elemVec[i].WriteHeader(file);
    ok = file.WriteHeaderLine(std::string("end_header\n"));
    return ok;
  }


  inline void Info::AddPlyElement(PlyElement& pe)
  {
    elemVec.push_back(pe);
  }


  inline size_t Info::GetElementCount(std::string& name)
  {
    PlyElement* pe = GetElement(name);
    if (pe != NULL)
      return pe->cnt;
    return -1;
  }


  inline size_t Info::GetElementCount(PlyElemEntity e)
  {
    PlyElement* pe = GetElement(e);
    if (pe != NULL)
      return pe->cnt;
    return -1;
  }


  inline size_t Info::GetVertexCount()
  {
    return GetElementCount(PlyElemEntity::NNP_VERTEX_ELEM);
  }


  inline size_t Info::GetFaceCount()
  {
    return GetElementCount(PlyElemEntity::NNP_FACE_ELEM);
  }


  inline size_t Info::GetEdgeCount()
  {
    return GetElementCount(PlyElemEntity::NNP_EDGE_ELEM);
  }


  inline PlyElement* Info::GetElement(std::string& name)
  {
    for (int i = 0; i < elemVec.size(); i++)
    {
      if (elemVec[i].name == name)
        return &elemVec[i];
    }
    return NULL;
  }


  inline PlyElement* Info::GetElement(PlyElemEntity e)
  {
    for (int i = 0; i < elemVec.size(); i++)
    {
      if (elemVec[i].plyElem == e)
        return &elemVec[i];
    }
    return NULL;
  }


  inline PlyElement* Info::GetVertexElement()
  {
    return GetElement(PlyElemEntity::NNP_VERTEX_ELEM);
  }


  inline PlyElement* Info::GetFaceElement()
  {
    return GetElement(PlyElemEntity::NNP_FACE_ELEM);
  }


  inline PlyElement* Info::GetEdgeElement()
  {
    return GetElement(PlyElemEntity::NNP_EDGE_ELEM);
  }



  /** Abstract class for the descriptor of a Ply propertie.
  *	The class defines how a PlyProperty is saved in memory.
  */
  class DescriptorInterface
  {

  public:

    int64_t curPos;		/**< Position of the next property to read or to write. */

    void *base;			/**< Pointer to the memory location that contains the data of the property. */

    PlyEntity elem;		/**< Ply entity managed by the descriptor. */

    std::string name;	/**< Name of the PlyProperty*/

    /**
    * Constructor of the descriptor.
    *
    * @param _e	Ply entity managed by the descriptor.
    * @param _b	Pointer to the memory location that contains the data of the property.
    */
    inline DescriptorInterface(PlyEntity _e, void *_b) :curPos(0), elem(_e), base(_b), name(PlyPropertyName(_e)[0]){};

    /**
    * Constructor of the descriptor.
    *
    * @param _s		Name of the PlyProperty.
    * @param _b	Pointer to the memory location that contains the data of the property.
    */
    inline DescriptorInterface(std::string& _s, void *_b) :curPos(0), elem(PlyEntity::NNP_UNKNOWN_ENTITY), name(_s), base(_b){};

    /**
    * Restart the descriptor.
    */
    virtual void Restart() = 0;

    /**
    * Read the property data from the binary file.
    *
    * @param file		Input file.
    * @param prop		PLY property to read from the file.
    * @param fixEndian	If true the method adjust the endianess of the data.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    virtual bool ReadElemBinary(PlyFile &file, PlyProperty &prop, bool fixEndian) = 0;

    /**
    * Read the property data from the ascii file.
    *
    * @param file	Input file.
    * @param prop	PLY property to read from the file.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    virtual bool ReadElemAscii(PlyFile &file, PlyProperty &prop) = 0;

    /**
    * Write the property data in the binary file.
    *
    * @param file		Input file.
    * @param prop		PLY property to write in the file.
    * @param fixEndian	If true the method adjust the endianess of the data.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    virtual bool WriteElemBinary(PlyFile &file, PlyProperty &prop, bool fixEndian) = 0;

    /**
    * Write the property data in the ascii file.
    *
    * @param file	Input file.
    * @param prop	PLY property to write in the file.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    virtual bool WriteElemAscii(PlyFile &file, PlyProperty &prop) = 0;
  };



  /** Memory descriptor of a Ply element.
  *	The class defines how a PlyElement is saved in memory.
  */
  class ElementDescriptor
  {
  public:

    typedef std::vector<nanoply::DescriptorInterface*> PropertyDescriptor;

    std::string name;					/**< Name of the Ply element. */
    PlyElemEntity elem;					/**< PLY Element Entity. */
    PropertyDescriptor dataDescriptor;	/**< Vector of property descriptor. */

    /**
    * Constructor of the Ply element descriptor.
    *
    * @param _e		Ply Element entity managed by the descriptor.
    */
    inline ElementDescriptor(PlyElemEntity _e) : elem(_e), name(PlyElementName(_e)[0]){};

    /**
    * Constructor of the Ply element descriptor.
    *
    * @param _s		Name of the Ply element managed by the descriptor.
    */
    inline ElementDescriptor(std::string &_s) : elem(PlyElemEntity::NNP_UNKNOWN_ELEM), name(_s){};

    /**
    * Read all the properties of the element from the binary file.
    *
    * @param file		Input file.
    * @param elem		PLY element to read from the file.
    * @param fixEndian	If true the method adjust the endianess of the data.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    inline bool ReadElemBinary(PlyFile &file, PlyElement &elem, bool fixEndian);

    /**
    * Read all the property data of the element from the ascii file.
    *
    * @param file	Input file.
    * @param elem	PLY element to read from the file.
    * @return		If successful returns true. Otherwise, it returns false.
    */
    inline bool ReadElemAscii(PlyFile &file, PlyElement &elem);

    /**
    * Write all the property data of the element in the binary file.
    *
    * @param file		Input file.
    * @param elem		PLY element to write from the file.
    * @param fixEndian	If true the method adjust the endianess of the data.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    inline bool WriteElemBinary(PlyFile &file, PlyElement &elem, bool fixEndian);

    /**
    * Write all the property data of the element in the ascii file.
    *
    * @param file		Input file.
    * @param elem		PLY element to write from the file.
    * @return			If successful returns true. Otherwise, it returns false.
    */
    inline bool WriteElemAscii(PlyFile &file, PlyElement &elem);

    /**
    * Check if the properties defined in input element have a proper data descriport to write in the file.
    * It sets the variable "validToWrite" for all the Ply properties with a data descriport.
    *
    * @param elem		PLY element to write from the file.
    */
    inline void CheckDescriptor(PlyElement &elem);

  private:

    inline void ExtractDescriptor(PropertyDescriptor &descr, PlyElement &elem);

  };


  inline void ElementDescriptor::ExtractDescriptor(PropertyDescriptor& descr, PlyElement &elem)
  {
    for (int j = 0; j < elem.propVec.size(); j++)
    {
      PlyProperty& prop = elem.propVec[j];
      int i = 0;
      for (; i < dataDescriptor.size(); i++)
      {
        if (dataDescriptor[i]->elem == prop.elem)
        {
          if (prop.elem != PlyEntity::NNP_UNKNOWN_ENTITY)
          {
            descr.push_back(dataDescriptor[i]);
            break;
          }
          else //if (dataDescriptor[i]->name == prop.name)
          {
            std::string name1(dataDescriptor[i]->name);
            std::string name2(prop.name);
            std::transform(name1.begin(), name1.end(), name1.begin(), ::tolower);
            std::transform(name2.begin(), name2.end(), name2.begin(), ::tolower);
            if (name1 == name2)
            {
              descr.push_back(dataDescriptor[i]);
              break;
            }
          }
        }
      }
      if (i == dataDescriptor.size())
        descr.push_back(NULL);
    }
  }


  inline bool ElementDescriptor::ReadElemBinary(PlyFile &file, PlyElement &elem, bool fixEndian)
  {
    PropertyDescriptor descr;
    ExtractDescriptor(descr, elem);
    for (int i = 0; i < elem.cnt; i++)
    {
      for (int j = 0; j < elem.propVec.size(); j++)
      {
        PlyProperty& prop = elem.propVec[j];
        if (descr[j] != NULL)
          (*descr[j]).ReadElemBinary(file, prop, fixEndian);
        else
          prop.SkipBinaryPropertyInFile(file);
      }
    }
    return true;
  }


  inline bool ElementDescriptor::ReadElemAscii(PlyFile &file, PlyElement &elem)
  {
    PropertyDescriptor descr;
    ExtractDescriptor(descr, elem);
    for (int i = 0; i < elem.cnt; i++)
    {
      for (int j = 0; j < elem.propVec.size(); j++)
      {
        PlyProperty& prop = elem.propVec[j];
        if (descr[j] != NULL)
          (*descr[j]).ReadElemAscii(file, prop);
        else
          prop.SkipAsciiPropertyInFile(file);
      }
    }
    return true;
  }

  inline bool ElementDescriptor::WriteElemBinary(PlyFile &file, PlyElement &elem, bool fixEndian)
  {
    PropertyDescriptor descr;
    ExtractDescriptor(descr, elem);
    for (int i = 0; i < elem.cnt; i++)
    {
      for (int j = 0; j < elem.propVec.size(); j++)
      {
        if (descr[j] != NULL)
          (*descr[j]).WriteElemBinary(file, elem.propVec[j], fixEndian);
      }
    }
    return true;
  }

  inline bool ElementDescriptor::WriteElemAscii(PlyFile &file, PlyElement &elem)
  {
    PropertyDescriptor descr;
    ExtractDescriptor(descr, elem);
    for (int i = 0; i < elem.cnt; i++)
    {
      bool first = true;
      for (int j = 0; j < elem.propVec.size(); j++)
      {
        if (descr[j] != NULL)
        {
          if (!first)
            file.WriteAsciiData(std::string(" "));
          else
            first = false;
          (*descr[j]).WriteElemAscii(file, elem.propVec[j]);
        }
      }
      file.WriteAsciiData(std::string("\n"));
    }
    return true;
  }


  inline void ElementDescriptor::CheckDescriptor(PlyElement &elem)
  {
    if (elem.propVec.size() == 0)
    {
      elem.validToWrite = false;
      return;
    }
    elem.validToWrite = true;
    PropertyDescriptor descr;
    ExtractDescriptor(descr, elem);
    for (int j = 0; j < elem.propVec.size(); j++)
    {
      if (descr[j] != NULL)
        elem.propVec[j].validToWrite = true;
    }
  }



  /** Memory descriptor of a vector of properties.
  *	The class defines how a vector of PlyProperty is saved in memory.
  *
  *  @tparam CointainerType	Type of the container of the property
  *  @tparam VectorSize		Number of values stored in the property.
  *  @tparam ScalarType		Type of the values stored in the property.
  */
  template<class CointainerType, int VectorSize, typename ScalarType>
  class DataDescriptor : public DescriptorInterface
  {
  public:

    inline DataDescriptor();

    /**
    * Constructor of the descriptor.
    *
    * @param _e	Ply entity managed by the descriptor.
    * @param _b	Pointer to the memory location that contains the data of the property.
    */
    inline DataDescriptor(PlyEntity _e, void *_b) :DescriptorInterface(_e, _b){};

    /**
    * Constructor of the descriptor.
    *
    * @param _s	Name of the PlyProperty.
    * @param _b	Pointer to the memory location that contains the data of the property.
    */
    inline DataDescriptor(std::string& _s, void *_b) :DescriptorInterface(_s, _b){};

    inline void Restart();

    inline bool ReadElemBinary(PlyFile &file, PlyProperty &prop, bool fixEndian);

    inline bool ReadElemAscii(PlyFile &file, PlyProperty &prop);

    inline bool WriteElemBinary(PlyFile &file, PlyProperty &prop, bool fixEndian);

    inline bool WriteElemAscii(PlyFile &file, PlyProperty &prop);

  private:

    template<typename C>
    inline void ReadBinary(PlyFile &file, PlyProperty &prop, bool fixEndian);

    template<typename C>
    inline void ReadAscii(PlyFile &file, PlyProperty &prop);

    template<typename C>
    inline void WriteBinary(PlyFile &file, PlyProperty &prop, bool fixEndian);

    template<typename C>
    inline void WriteAscii(PlyFile &file, PlyProperty &prop);

  };


  template<class CointainerType, int VectorSize, typename ScalarType>
  inline void DataDescriptor<CointainerType, VectorSize, ScalarType>::Restart()
  {
    this->curPos = 0;
  }

  template<class ContainerType, int VectorSize, typename ScalarType>
  template<typename C>
  inline void DataDescriptor<ContainerType, VectorSize, ScalarType>::ReadBinary(PlyFile &file, PlyProperty &prop, bool fixEndian)
  {
    char * buffer = nullptr;
    int size;
    int count = prop.CountValue();
    int typeSize = prop.TypeSize();
    if (prop.type >= NNP_LIST_UINT8_UINT32)
    {
       file.ReadBinaryData(buffer, sizeof(char));
       const int cntList = int(*(reinterpret_cast<unsigned char *>(buffer)));
       size = typeSize * cntList;
       count = cntList;
    }
    else
      size = typeSize * count;
    file.ReadBinaryData(buffer, size);

    if (typeSize > 1 && fixEndian)
      adjustEndianess(reinterpret_cast<unsigned char *>(buffer), typeSize, count);

		unsigned char* baseProp = (unsigned char*)base + this->curPos*sizeof(ContainerType);
    C* temp = (C*)buffer;
    if ((prop.elem == NNP_CRGB || prop.elem == NNP_CRGBA))
    {
			float norm = 1.0f;
      if (std::is_same<ScalarType, float>::value && std::is_same<C, unsigned char>::value)
        norm = 1.0f / 255.0f;
      else if (std::is_same<ScalarType, unsigned char>::value && std::is_same<C, float>::value)
        norm = 255.0f;
			for (int i = 0; i < std::min(VectorSize, count); i++)
				*((ScalarType *)(baseProp + i*sizeof(ScalarType))) = ScalarType(temp[i] * norm);
    }
		else
		{
			for (int i = 0; i < std::min(VectorSize, count); i++)
				*((ScalarType *)(baseProp + i*sizeof(ScalarType))) = ScalarType(temp[i]);
		}
    ++(this->curPos);
  }


  template<class ContainerType, int VectorSize, typename ScalarType>
  inline bool DataDescriptor<ContainerType, VectorSize, ScalarType>::ReadElemBinary(PlyFile &file, PlyProperty &prop, bool fixEndian)
  {
    if (prop.elem != elem)
      return false;
    switch (prop.type)
    {
    case NNP_LIST_INT8_INT8:
    case NNP_LIST_UINT8_INT8:
    case NNP_INT8:				this->ReadBinary<char>(file, prop, fixEndian); break;
    case NNP_LIST_INT8_UINT8:
    case NNP_LIST_UINT8_UINT8:
    case NNP_UINT8:				this->ReadBinary<unsigned char>(file, prop, fixEndian); break;
    case NNP_LIST_INT8_INT16:
    case NNP_LIST_UINT8_INT16:
    case NNP_INT16:				this->ReadBinary<short>(file, prop, fixEndian); break;
    case NNP_LIST_INT8_UINT16:
    case NNP_LIST_UINT8_UINT16:
    case NNP_UINT16:			this->ReadBinary<unsigned short>(file, prop, fixEndian); break;
    case NNP_LIST_INT8_FLOAT32:
    case NNP_LIST_UINT8_FLOAT32:
    case NNP_FLOAT32:			this->ReadBinary<float>(file, prop, fixEndian); break;
    case NNP_LIST_UINT8_INT32:
    case NNP_LIST_INT8_INT32:
    case NNP_INT32:				this->ReadBinary<int>(file, prop, fixEndian); break;
    case NNP_LIST_UINT8_UINT32:
    case NNP_LIST_INT8_UINT32:
    case NNP_UINT32:			this->ReadBinary<unsigned int>(file, prop, fixEndian); break;
    case NNP_LIST_INT8_FLOAT64:
    case NNP_LIST_UINT8_FLOAT64:
    case NNP_FLOAT64:			this->ReadBinary<double>(file, prop, fixEndian); break;
    }
    return true;
  }



  template<class ContainerType, int VectorSize, typename ScalarType>
  template<typename C>
  inline void DataDescriptor<ContainerType, VectorSize, ScalarType>::ReadAscii(PlyFile &file, PlyProperty &prop)
  {
    int count = prop.CountValue();
    if (prop.type >= NNP_LIST_UINT8_UINT32)
      file.ReadAsciiData(count);

    C* temp = new C[count];
    for (int i = 0; i < count; i++)
      file.ReadAsciiData(temp[i]);

    unsigned char* baseProp = (unsigned char*)base + this->curPos*sizeof(ContainerType);
    if ((prop.elem == NNP_CRGB || prop.elem == NNP_CRGBA))
    {
			float norm = 1.0f;
      if (std::is_same<ScalarType, float>::value && prop.type == NNP_UINT8)
        norm = 1.0f / 255.0f;
      else if (std::is_same<ScalarType, unsigned char>::value && prop.type == NNP_FLOAT32)
        norm = 255.0f;
			for (int i = 0; i < std::min(VectorSize, count); i++)
				*((ScalarType *)(baseProp + i*sizeof(ScalarType))) = ScalarType(temp[i] * norm);
    }
		else
		{
			for (int i = 0; i < std::min(VectorSize, count); i++)
				*((ScalarType *)(baseProp + i*sizeof(ScalarType))) = ScalarType(temp[i]);
		}
    delete[] temp;
    ++(this->curPos);
  }



  template<class ContainerType, int VectorSize, typename ScalarType>
  inline bool DataDescriptor<ContainerType, VectorSize, ScalarType>::ReadElemAscii(PlyFile &file, PlyProperty &prop)
  {
    if (prop.elem != elem)
      return false;
    switch (prop.type)
    {
    case NNP_LIST_UINT8_INT8:
    case NNP_LIST_INT8_INT8:
    case NNP_INT8:				this->ReadAscii<int>(file, prop); break;
    case NNP_LIST_UINT8_UINT8:
    case NNP_LIST_INT8_UINT8:
    case NNP_UINT8:				this->ReadAscii<unsigned int>(file, prop); break;
    case NNP_LIST_UINT8_INT16:
    case NNP_LIST_INT8_INT16:
    case NNP_INT16:				this->ReadAscii<short>(file, prop); break;
    case NNP_LIST_UINT8_UINT16:
    case NNP_LIST_INT8_UINT16:
    case NNP_UINT16:			this->ReadAscii<unsigned short>(file, prop); break;
    case NNP_LIST_UINT8_FLOAT32:
    case NNP_LIST_INT8_FLOAT32:
    case NNP_FLOAT32:			this->ReadAscii<float>(file, prop); break;
    case NNP_LIST_UINT8_INT32:
    case NNP_LIST_INT8_INT32:
    case NNP_INT32:				this->ReadAscii<int>(file, prop); break;
    case NNP_LIST_UINT8_UINT32:
    case NNP_LIST_INT8_UINT32:
    case NNP_UINT32:			this->ReadAscii<unsigned int>(file, prop); break;
    case NNP_LIST_UINT8_FLOAT64:
    case NNP_LIST_INT8_FLOAT64:
    case NNP_FLOAT64:			this->ReadAscii<double>(file, prop); break;
    }
    return true;
  }



  template<class ContainerType, int VectorSize, typename ScalarType>
  template<typename C>
  inline void DataDescriptor<ContainerType, VectorSize, ScalarType>::WriteBinary(PlyFile &file, PlyProperty &prop, bool fixEndian)
  {
    (void)fixEndian;
    int count = prop.CountValue();
    C data[VectorSize];
    if (prop.type >= NNP_LIST_UINT8_UINT32)
    {
      if (prop.IsSigned())
      {
        char listSize = (char)VectorSize;
        file.WriteBinaryData(&listSize, 1);
        count = VectorSize;
      }
      else
      {
        unsigned char listSize = (unsigned char)VectorSize;
        file.WriteBinaryData(&listSize, 1);
        count = VectorSize;
      }
    }

    C temp = 0;
		unsigned char* baseProp = (unsigned char*)base + this->curPos*sizeof(ContainerType);
    if ((prop.elem == NNP_CRGB || prop.elem == NNP_CRGBA))
    {
			float norm = 1.0f;
      if (std::is_same<ScalarType, float>::value && std::is_same<C, unsigned char>::value)
        norm = 255.0f;
      else if (std::is_same<ScalarType, unsigned char>::value && std::is_same<C, float>::value)
        norm = 1.0f / 255.0f;
			for (int i = 0; i < std::min(VectorSize, count); i++)
				data[i] = (C)((*(ScalarType*)(baseProp + i*sizeof(ScalarType))) * norm);
    }
		else
		{
			for (int i = 0; i < std::min(VectorSize, count); i++)
				data[i] = (C)((*(ScalarType*)(baseProp + i*sizeof(ScalarType))));
		}

    if (sizeof(C) > 1 && fixEndian)
      adjustEndianess((unsigned char*)data, sizeof(C), std::min(VectorSize, count));

    file.WriteBinaryData(data, sizeof(C)*std::min(VectorSize, count));
    for (int i = 0; i < (count - VectorSize); i++)
      file.WriteBinaryData(&temp, sizeof(C));
    ++(this->curPos);
  }


  template<class ContainerType, int VectorSize, typename ScalarType>
  inline bool DataDescriptor<ContainerType, VectorSize, ScalarType>::WriteElemBinary(PlyFile &file, PlyProperty &prop, bool fixEndian)
  {
    if (prop.elem != elem)
      return false;
    switch (prop.type)
    {
    case NNP_LIST_INT8_INT8:
    case NNP_LIST_UINT8_INT8:
    case NNP_INT8:				this->WriteBinary<char>(file, prop, fixEndian); break;
    case NNP_LIST_INT8_UINT8:
    case NNP_LIST_UINT8_UINT8:
    case NNP_UINT8:				this->WriteBinary<unsigned char>(file, prop, fixEndian); break;
    case NNP_LIST_INT8_INT16:
    case NNP_LIST_UINT8_INT16:
    case NNP_INT16:				this->WriteBinary<short>(file, prop, fixEndian); break;
    case NNP_LIST_INT8_UINT16:
    case NNP_LIST_UINT8_UINT16:
    case NNP_UINT16:			this->WriteBinary<unsigned short>(file, prop, fixEndian); break;
    case NNP_LIST_INT8_FLOAT32:
    case NNP_LIST_UINT8_FLOAT32:
    case NNP_FLOAT32:			this->WriteBinary<float>(file, prop, fixEndian); break;
    case NNP_LIST_UINT8_INT32:
    case NNP_LIST_INT8_INT32:
    case NNP_INT32:				this->WriteBinary<int>(file, prop, fixEndian); break;
    case NNP_LIST_UINT8_UINT32:
    case NNP_LIST_INT8_UINT32:
    case NNP_UINT32:			this->WriteBinary<unsigned int>(file, prop, fixEndian); break;
    case NNP_LIST_INT8_FLOAT64:
    case NNP_LIST_UINT8_FLOAT64:
    case NNP_FLOAT64:			this->WriteBinary<double>(file, prop, fixEndian); break;
    }
    return true;
  }


  template<class ContainerType, int VectorSize, typename ScalarType>
  template<typename C>
  inline void DataDescriptor<ContainerType, VectorSize, ScalarType>::WriteAscii(PlyFile &file, PlyProperty &prop)
  {
    int count = prop.CountValue();
    if (prop.type >= NNP_LIST_UINT8_UINT32)
    {
      if (prop.IsSigned())
      {
        int listSize = (int)VectorSize;
        file.WriteAsciiData(listSize);
        count = VectorSize;
      }
      else
      {
        unsigned int listSize = (unsigned int)VectorSize;
        file.WriteAsciiData(listSize);
        count = VectorSize;
      }
      file.WriteAsciiData(std::string(" "));
    }

    C data[VectorSize];
		unsigned char* baseProp = (unsigned char*)base + this->curPos*sizeof(ContainerType);
    if ((prop.elem == NNP_CRGB || prop.elem == NNP_CRGBA))
    {
			float norm = 1.0;
      if (std::is_same<ScalarType, float>::value && prop.type == NNP_UINT8)
        norm = 255.0f;
      else if (std::is_same<ScalarType, unsigned char>::value && prop.type == NNP_FLOAT32)
        norm = 1.0f / 255.0f;
			for (int i = 0; i < std::min(VectorSize, count); i++)
				data[i] = (C)((*(ScalarType*)(baseProp + i*sizeof(ScalarType))) * norm);
    }
		else
		{
			for (int i = 0; i < std::min(VectorSize, count); i++)
				data[i] = (C)((*(ScalarType*)(baseProp + i*sizeof(ScalarType))));
		}
		   
   
    for (int i = 0; i < (count - VectorSize); i++)
      data[i] = 0;

    for (int i = 0; i < count; i++)
    {
      file.WriteAsciiData(data[i]);
      if (i < count - 1)
        file.WriteAsciiData(std::string(" "));
    }
    ++(this->curPos);
  }


  template<class ContainerType, int VectorSize, typename ScalarType>
  inline bool DataDescriptor<ContainerType, VectorSize, ScalarType>::WriteElemAscii(PlyFile &file, PlyProperty& prop)
  {
    if (prop.elem != elem)
      return false;
    if (prop.elem == PlyEntity::NNP_UNKNOWN_ENTITY && prop.name != name)
      return false;
    switch (prop.type)
    {
    case NNP_LIST_UINT8_INT8:
    case NNP_LIST_INT8_INT8:
    case NNP_INT8:				this->WriteAscii<int>(file, prop); break;
    case NNP_LIST_UINT8_UINT8:
    case NNP_LIST_INT8_UINT8:
    case NNP_UINT8:				this->WriteAscii<unsigned int>(file, prop); break;
    case NNP_LIST_UINT8_INT16:
    case NNP_LIST_INT8_INT16:
    case NNP_INT16:				this->WriteAscii<short>(file, prop); break;
    case NNP_LIST_UINT8_UINT16:
    case NNP_LIST_INT8_UINT16:
    case NNP_UINT16:			this->WriteAscii<unsigned short>(file, prop); break;
    case NNP_LIST_UINT8_FLOAT32:
    case NNP_LIST_INT8_FLOAT32:
    case NNP_FLOAT32:			this->WriteAscii<float>(file, prop); break;
    case NNP_LIST_UINT8_INT32:
    case NNP_LIST_INT8_INT32:
    case NNP_INT32:				this->WriteAscii<int>(file, prop); break;
    case NNP_LIST_UINT8_UINT32:
    case NNP_LIST_INT8_UINT32:
    case NNP_UINT32:			this->WriteAscii<unsigned int>(file, prop); break;
    case NNP_LIST_UINT8_FLOAT64:
    case NNP_LIST_INT8_FLOAT64:
    case NNP_FLOAT64:			this->WriteAscii<double>(file, prop); break;
    }
    return true;
  }



  template <size_t ActionType>
  inline bool ElemProcessing(ElementDescriptor& elemDescr, PlyElement &elem, PlyFile& file, bool fixEndian)
  {
    if ((elemDescr.elem != PlyElemEntity::NNP_UNKNOWN_ELEM && elemDescr.elem == elem.plyElem) ||
      (elemDescr.elem == PlyElemEntity::NNP_UNKNOWN_ELEM && elemDescr.name == elem.name))
    {
      if (ActionType == 0)
        elemDescr.ReadElemBinary(file, elem, fixEndian);
      else if (ActionType == 1)
        elemDescr.ReadElemAscii(file, elem);
      else if (ActionType == 2)
        elemDescr.WriteElemBinary(file, elem, fixEndian);
      else if (ActionType == 3)
        elemDescr.WriteElemAscii(file, elem);
      else if (ActionType == 4)
        elemDescr.CheckDescriptor(elem);
      return true;
    }
    return false;
  }

  typedef std::vector<ElementDescriptor*> MeshDescriptor;

  /**
  * Load a 3D model from a PLY file.
  *
  * @param meshElements			Vector that defines how to manage the ply element data in memory.
  * @param info					Info of the file to load.
  */
  inline bool OpenModel(Info& info, MeshDescriptor& meshElements)
  {
    PlyFile file;
    if (!file.OpenFileToRead(info.filename))
    {
      info.errInfo = NNP_UNABLE_TO_OPEN;
      return false;
    }
    bool last;
    std::string line;
    do
    {
      file.NextHeaderLine(line, last);
    } while (!last);
    bool fixEndian = false;
    if (checkEndianness() == 1)
    {
      if (info.bigEndian)
        fixEndian = true;
    }
    else
    {
      if (!info.bigEndian)
        fixEndian = true;
    }

    if (info.binary)
    {
      for (int i = 0; i < info.elemVec.size(); ++i)
      {
        PlyElement& pe = info.elemVec[i];
        int j = 0;
        for (; j < meshElements.size(); j++)
          if (ElemProcessing<0>(*meshElements[j], pe, file, fixEndian))
            break;
        if (j == meshElements.size())
          pe.SkipBinaryElementsInFile(file);
        //if (!TupleForEach(meshElements, pe, file, fixEndian, SizeT<0>()))
        //	pe.SkipBinaryElementsInFile(file);
      }
    }
    else
    {
      for (int i = 0; i < info.elemVec.size(); ++i)
      {
        PlyElement& pe = info.elemVec[i];
        int j = 0;
        for (; j < meshElements.size(); j++)
          if (ElemProcessing<1>(*meshElements[j], pe, file, false))
            break;
        if (j == meshElements.size())
          pe.SkipAsciiElementsInFile(file);
        //if (!TupleForEach(meshElements, pe, file, fixEndian, SizeT<1>()))
        //	pe.SkipAsciiElementsInFile(file);
      }
    }
    return true;
  }



  /**
  * Save a 3D model in a PLY file.
  *
  * @param filename				Path to the file to save.
  * @param meshElements			Vector that defines how to manage the ply element data in memory.
  * @param info					Info to saved in the PLY header.
  */
  inline bool SaveModel(std::string& filename, MeshDescriptor& meshElements, Info& info)
  {
    PlyFile file;
    if (!file.OpenFileToWrite(filename))
    {
      info.errInfo = NNP_UNABLE_TO_OPEN;
      return false;
    }
    for (int i = 0; i < info.elemVec.size(); ++i)
    {
      PlyElement& pe = info.elemVec[i];
      for (int j = 0; j < meshElements.size(); j++)
        if (ElemProcessing<4>(*meshElements[j], pe, file, false))
          break;
    }
    info.WriteHeader(file);
    bool fixEndian = false;
    if (checkEndianness() == 1)
    {
      if (info.bigEndian)
        fixEndian = true;
    }
    else
    {
      if (!info.bigEndian)
        fixEndian = true;
    }

    if (info.binary)
    {
      for (int i = 0; i < info.elemVec.size(); ++i)
      {
        PlyElement& pe = info.elemVec[i];
        if (pe.validToWrite)
        {
          for (int j = 0; j < meshElements.size(); j++)
            if (ElemProcessing<2>(*meshElements[j], pe, file, false))
              break;
        }
      }
    }
    else
    {
      for (int i = 0; i < info.elemVec.size(); ++i)
      {
        PlyElement& pe = info.elemVec[i];
        if (pe.validToWrite)
        {
          for (int j = 0; j < meshElements.size(); j++)
            if (ElemProcessing<3>(*meshElements[j], pe, file, false))
              break;
        }
      }
    }
    file.Flush();
    return true;
  }



  /**
  * @cond HIDDEN_SYMBOLS
  */
  template < typename TupleType, size_t ActionType>
  inline bool TupleForEach(TupleType &tuple, PlyElement &elem, PlyFile& file, bool fixEndian, SizeT<ActionType> a)
  {
    return TupleForEach(tuple, elem, file, fixEndian, SizeT<std::tuple_size<TupleType>::value>(), a);
  }

  template < typename TupleType, size_t ActionType>
  inline bool TupleForEach(TupleType &tuple, PlyElement &elem, PlyFile& file, bool fixEndian, SizeT<0> t, SizeT<ActionType> a) { return false; }

  template < typename TupleType, size_t N, size_t ActionType>
  inline bool TupleForEach(TupleType &tuple, PlyElement &elem, PlyFile& file, bool fixEndian, SizeT<N> t, SizeT<ActionType> a)
  {
    typename std::tuple_element<N - 1, TupleType>::type &elemDescr = std::get<N - 1>(tuple);
    if ((elemDescr.elem != PlyElemEntity::NNP_UNKNOWN_ELEM && elemDescr.elem == elem.plyElem) ||
      (elemDescr.elem == PlyElemEntity::NNP_UNKNOWN_ELEM && elemDescr.name == elem.name))
    {
      if (ActionType == 0)
        elemDescr.ReadElemBinary(file, elem, fixEndian);
      else if (ActionType == 1)
        elemDescr.ReadElemAscii(file, elem);
      else if (ActionType == 2)
        elemDescr.WriteElemBinary(file, elem, fixEndian);
      else if (ActionType == 3)
        elemDescr.WriteElemAscii(file, elem);
      else if (ActionType == 4)
        elemDescr.CheckDescriptor(elem);
      return true;
    }
    return TupleForEach(tuple, elem, file, fixEndian, SizeT<N - 1>(), a);
  }
  /**
  * @endcond
  */

} // end namespace nanoply
#endif // NANOPLY_HPP
