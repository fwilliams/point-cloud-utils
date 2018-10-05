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

#ifndef __VCGLIB_GRID_UTIL
#define __VCGLIB_GRID_UTIL

#include<vcg/space/index/base.h>
#include<vcg/space/box3.h>
#include <vcg/space/index/space_iterators.h>


#ifndef WIN32
#define __int64 long long
#define __cdecl
#endif

namespace vcg {

    /** BasicGrid

    Basic Class abstracting a gridded structure in a 3d space;
    Usueful for having coherent float to integer conversion in a unique place:
    Some Notes:
        - bbox is the real occupation of the box in the space;
        - siz is the number of cells for each side

    OBJTYPE:      Type of the indexed objects.
    SCALARTYPE:   Scalars type for structure's internal data (may differ from
                  object's scalar type).

    */

template <class SCALARTYPE>
class BasicGrid //:public SpatialIndex<SCALARTYPE>
{
public:

    typedef SCALARTYPE ScalarType;
    typedef Box3<ScalarType> Box3x;
    typedef Point3<ScalarType> CoordType;
    typedef BasicGrid<SCALARTYPE> GridType;

    Box3x bbox;

    CoordType dim;		/// Spatial Dimention (edge legth) of the bounding box
    Point3i siz;		/// Number of cells forming the grid
    CoordType voxel;	/// Dimensions of a single cell

    /*
     Derives the right values of Dim and voxel starting
     from the current values of siz and bbox
    */
    void ComputeDimAndVoxel()
        {
            this->dim  = this->bbox.max - this->bbox.min;
            this->voxel[0] = this->dim[0]/this->siz[0];
            this->voxel[1] = this->dim[1]/this->siz[1];
            this->voxel[2] = this->dim[2]/this->siz[2];
        }
    /* Given a 3D point, returns the coordinates of the cell where the point is
     * @param p is a 3D point
     * @return integer coordinates of the cell
     */
    inline Point3i GridP( const Point3<ScalarType> & p ) const
    {
        Point3i pi;
        PToIP(p, pi);
        return pi;
    }

    /* Given a 3D point p, returns the index of the corresponding cell
     * @param p is a 3D point in the space
     * @return integer coordinates pi of the cell
     */
    inline void PToIP(const CoordType & p, Point3i &pi ) const
    {
        CoordType t = p - bbox.min;
        pi[0] = int( t[0] / voxel[0] );
        pi[1] = int( t[1] / voxel[1] );
        pi[2] = int( t[2] / voxel[2] );
    }

    /* Given a cell index return the lower corner of the cell
     * @param integer coordinates pi of the cell
     * @return p is a 3D point representing the lower corner of the cell
     */
    template <class OtherScalarType>
    inline void IPiToPf(const Point3i & pi,  Point3<OtherScalarType> &p ) const
    {
        p[0] = bbox.min[0] + ((OtherScalarType)pi[0])*voxel[0];
        p[1] = bbox.min[1] + ((OtherScalarType)pi[1])*voxel[1];
        p[2] = bbox.min[2] + ((OtherScalarType)pi[2])*voxel[2];
    }
    /* Returns the matrix that applied to a point in grid space
     * transforms it in the original space.
     */
    inline Matrix44<ScalarType> IPtoPfMatrix() const
    {
      Matrix44<ScalarType> m; m.SetScale(voxel);
      Matrix44<ScalarType> t; t.SetTranslate(bbox.min);
      return t*m;
    }
    /* Given a cell index return the corresponding box
     * @param integer coordinates pi of the cell
     * @return b is the corresponding box in <ScalarType> coordinates
     */
    inline void IPiToBox(const Point3i & pi, Box3x & b ) const
    {
        CoordType p;
        p[0] = ((ScalarType)pi[0])*voxel[0];
        p[1] = ((ScalarType)pi[1])*voxel[1];
        p[2] = ((ScalarType)pi[2])*voxel[2];
        p += bbox.min;
        b.min = p;
        b.max = (p + voxel);
    }

    /* Given a cell index return the center of the cell itself
     * @param integer coordinates pi of the cell
     * @return b is the corresponding box in <ScalarType> coordinates
     */inline void IPiToBoxCenter(const Point3i & pi, CoordType & c ) const
    {
        CoordType p;
        IPiToPf(pi,p);
        c = p + voxel/ScalarType(2.0);
    }

    // Same of IPiToPf but for the case that you just want to transform
    // from a space to the other.
    template <class OtherScalarType>
    void IPfToPf(const Point3<OtherScalarType> & pi, Point3<OtherScalarType> &p ) const
    {
        p[0] = ((OtherScalarType)pi[0])*voxel[0] + bbox.min[0];
        p[1] = ((OtherScalarType)pi[1])*voxel[1] + bbox.min[1];
        p[2] = ((OtherScalarType)pi[2])*voxel[2] + bbox.min[2];
    }

    /* Given a cell in <ScalarType> coordinates, compute the corresponding cell in integer coordinates
     * @param b is the cell in <ScalarType> coordinates
     * @return ib is the correspondent box in integer coordinates
     */
  inline void BoxToIBox( const Box3x & b, Box3i & ib ) const
    {
    PToIP(b.min, ib.min);
    PToIP(b.max, ib.max);
        //assert(ib.max[0]>=0 && ib.max[1]>=0 && ib.max[2]>=0);
    }

    /* Given a cell in integer coordinates, compute the corresponding cell in <ScalarType> coordinates
     * @param ib is the cell in integer coordinates
     * @return b is the correspondent box in <ScalarType> coordinates
     */
    /// Dato un box in voxel ritorna gli estremi del box reale
    void IBoxToBox( const Box3i & ib, Box3x & b ) const
    {
        IPiToPf(ib.min,b.min);
        IPiToPf(ib.max+Point3i(1,1,1),b.max);
    }
};

template<class scalar_type>
void BestDim( const Box3<scalar_type> box, const scalar_type voxel_size, Point3i & dim )
{
    Point3<scalar_type> box_size = box.max-box.min;
    __int64 elem_num = (__int64)(box_size[0]/voxel_size +0.5) *( __int64)(box_size[1]/voxel_size +0.5) * (__int64)(box_size[2]/voxel_size +0.5);
    BestDim(elem_num,box_size,dim);
}
    /** Calcolo dimensioni griglia.
    Calcola la dimensione della griglia in funzione
    della ratio del bounding box e del numero di elementi
    */
    template<class scalar_type>
    void BestDim( const __int64 elems, const Point3<scalar_type> & size, Point3i & dim )
    {
        const __int64 mincells   = 1;		// Numero minimo di celle
        const double GFactor = 1;	// GridEntry = NumElem*GFactor
        double diag = size.Norm();	// Diagonale del box
        double eps  = diag*1e-4;		// Fattore di tolleranza

        assert(elems>0);
        assert(size[0]>=0.0);
        assert(size[1]>=0.0);
        assert(size[2]>=0.0);


        __int64 ncell = (__int64)(elems*GFactor);	// Calcolo numero di voxel
        if(ncell<mincells)
            ncell = mincells;

        dim[0] = 1;
        dim[1] = 1;
        dim[2] = 1;

        if(size[0]>eps)
        {
            if(size[1]>eps)
            {
                if(size[2]>eps)
                {
                    double k = pow((double)(ncell/(size[0]*size[1]*size[2])),double(1.0/3.f));
                    dim[0] = int(size[0] * k);
                    dim[1] = int(size[1] * k);
                    dim[2] = int(size[2] * k);
                }
                else
                {
                    dim[0] = int(::sqrt(ncell*size[0]/size[1]));
                    dim[1] = int(::sqrt(ncell*size[1]/size[0]));
                }
            }
            else
            {
                if(size[2]>eps)
                {
                    dim[0] = int(::sqrt(ncell*size[0]/size[2]));
                    dim[2] = int(::sqrt(ncell*size[2]/size[0]));
                }
                else
                    dim[0] = int(ncell);
            }
        }
        else
        {
            if(size[1]>eps)
            {
                if(size[2]>eps)
                {
                    dim[1] = int(::sqrt(ncell*size[1]/size[2]));
                    dim[2] = int(::sqrt(ncell*size[2]/size[1]));
                }
                else
                    dim[1] = int(ncell);
            }
            else if(size[2]>eps)
                dim[2] = int(ncell);
        }
    dim[0] = std::max(dim[0],1);
    dim[1] = std::max(dim[1],1);
    dim[2] = std::max(dim[2],1);
    }
}
#endif
