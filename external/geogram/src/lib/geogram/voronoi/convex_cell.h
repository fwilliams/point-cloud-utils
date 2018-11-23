/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef GEOGRAM_VORONOI_CONVEX_CELL
#define GEOGRAM_VORONOI_CONVEX_CELL

#ifndef STANDALONE_CONVEX_CELL
#  include <geogram/basic/common.h>
#  include <geogram/basic/memory.h>
#  include <geogram/basic/numeric.h>
#  include <geogram/basic/geometry.h>
#endif

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>



/**
 * \file geogram/voronoi/convex_cell.h
 * \brief Class to compute the intersection of a set of half-spaces in 3D.
 * \details Has its own types for points and vectors because it can be used
 *  independently from Geogram. In that case, define STANDALONE_CONVEX_CELL
 */


namespace VBW {

#ifdef STANDALONE_CONVEX_CELL
    using std::vector;
    typedef unsigned int index_t;
    typedef unsigned int global_index_t;
#   define vbw_assert(x) assert(x)
    struct vec2 {
	double x;
	double y;
    };
    struct vec3 {
	double x;
	double y;
	double z;
    };
    struct vec4 {
	double x;
	double y;
	double z;
	double w;
    };
#else    
    using GEO::vector;
    typedef unsigned int index_t;        // Always 32 bits
    typedef GEO::index_t global_index_t; // Possibly 64 bits in GARGANTUA mode
#   define vbw_assert(x) geo_debug_assert(x)
    using GEO::vec2;
    using GEO::vec3;
    using GEO::vec4;
#endif    
    
/******************************************************************************/


    /**
     * \brief Creates a vec2 from its components.
     * \param[in] x , y the components of the 
     *  vector.
     * \return the created vector.
     */
    inline vec2 make_vec2(
	double x, double y
    ) {
	vec2 result;
	result.x = x;
	result.y = y;
	return result;
    }
    

    /**
     * \brief Creates a vec3 from its components.
     * \param[in] x , y , z the components of the 
     *  vector.
     * \return the created vector.
     */
    inline vec3 make_vec3(
	double x, double y, double z
    ) {
	vec3 result;
	result.x = x;
	result.y = y;
	result.z = z;
	return result;
    }


    /**
     * \brief Computes the cross product between 
     *  two vectors.
     * \param[in] v1 , v2 the two vectors.
     * \return the cross product between \p v1 and 
     *  \p v2.
     */
    inline vec3 cross(vec3 v1, vec3 v2) {
	return make_vec3(
	    v1.y*v2.z - v1.z*v2.y,
	    v1.z*v2.x - v1.x*v2.z,
	    v1.x*v2.y - v1.y*v2.x
	);
    }

    /**
     * \brief Computes the dot product between 
     *  two vectors.
     * \param[in] v1 , v2 the two vectors.
     * \return the dot product between \p v1 and 
     *  \p v2.
     */
    inline double dot(vec3 v1, vec3 v2) {
	return (
	    v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
	);
    }

    /**
     * \brief Computes the squared length of a vector.
     * \param[in] v the vector.
     * \return the squared length of \p v.
     */
    inline double squared_length(vec3 v) {
	return (v.x*v.x + v.y*v.y + v.z*v.z);
    }

    /**
     * \brief Computes the squared distance between two points.
     * \param[in] v , w the two points.
     * \return the squared distance between \p v and \p w.
     */
    inline double squared_distance(vec3 v, vec3 w) {
	double dx = w.x-v.x;
	double dy = w.y-v.y;
	double dz = w.z-v.z;
	return (dx*dx+dy*dy+dz*dz);
    }

    /**
     * \brief Computes the length of a vector.
     * \param[in] v the vector.
     * \return the length of \p v.
     */
    inline double length(vec3 v) {
	return ::sqrt(squared_length(v));
    }

    /**
     * \brief Computes a normalized vector.
     * \param[in] v the vector.
     * \return a vector with the same direction
     *  as \v and unit length.
     */
    inline vec3 normalize(vec3 v) {
	double s = 1.0/length(v);
	return make_vec3(
	    s*v.x, s*v.y, s*v.z
	);
    }
    
    /**
     * \brief Creates a vec4 from its components.
     * \param[in] x , y , z , w the components of the 
     *  vector.
     * \return the created vector.
     */
    inline vec4 make_vec4(
	double x, double y, double z, double w
    ) {
	vec4 result;
	result.x = x;
	result.y = y;
	result.z = z;
	result.w = w;
	return result;
    }

    /**
     * \brief Computes the dot product between 
     *  two vectors.
     * \param[in] v1 , v2 the two vectors.
     * \return the dot product between \p v1 and 
     *  \p v2.
     */
    inline double dot(vec4 v1, vec4 v2) {
	return (
	    v1.x*v2.x + v1.y*v2.y +
	    v1.z*v2.z + v1.w*v2.w
	);
    }

    /**
     * \brief Computes the squared length of a vector.
     * \param[in] v the vector.
     * \return the squared length of \p v.
     */
    inline double squared_length(vec4 v) {
	return (
	    v.x*v.x + v.y*v.y +
	    v.z*v.z + v.w*v.w
	);
    }

    /**
     * \brief Computes the length of a vector.
     * \param[in] v the vector.
     * \return the length of \p v.
     */
    inline double length(vec4 v) {
	return ::sqrt(squared_length(v));
    }
    
    /**
     * \brief Some constants for the flags
     *  in TriangleWithFlags.
     * \see TriangleWithFlags.
     */
    enum {
	CONFLICT_MASK  = 32768, /**< \brief The mask for conflict triangles.  */
	END_OF_LIST    = 32767, /**< \brief Constant to indicate end of list. */
	VERTEX_AT_INFINITY = 0  /**< \brief Vertex at infinity.               */
    };


    /**
     * \brief Type for flags.
     */
    typedef unsigned char uchar;
    
    /**
     * \brief Type for local indices.
     * \details Valid values are between 0 and 32766. 
     *  Full range is not used due to bookkeeping reasons,
     *  \see TriangleWithFlags.
     */
    typedef unsigned short ushort;

    /**
     * \brief A triangle with the local indices of its three
     *  vertices.
     */
    struct Triangle {
	ushort i;
	ushort j;
	ushort k;
    };

    /**
     * \brief Creates a triangle from its three vertices.
     * \param[in] i , j , k the local indices of the three
     *  vertices.
     * \return The created triangle.
     */
    inline Triangle make_triangle(
	ushort i, ushort j, ushort k
    ) {
	Triangle result;
	result.i = i;
	result.j = j;
	result.k = k;
	return result;
    }

    /**
     * \brief A triangle with flags.
     * \details The flags are used for two purposes:
     *  - bits [0..15] are used to chain the triangles: 
     *   there are two lists of triangles, the valid triangles 
     *   and the free list. End of list is indicated by value 32767.
     *  - bit 16 (32768) is set if the triangle is in conflict.
     */
    struct TriangleWithFlags : public Triangle {
	ushort flags;
    };

    inline TriangleWithFlags make_triangle_with_flags(
	ushort i, ushort j, ushort k, ushort f
    ) {
	TriangleWithFlags result;
	result.i = i;
	result.j = j;
	result.k = k;
	result.flags = f;
	return result;
    }


/******************************************************************************/

    inline double det2x2(
	double a11, double a12,
	double a21, double a22
    ) {
	return a11*a22 - a12*a21;
    }
    
    inline double det3x3(
	double a11, double a12, double a13,
	double a21, double a22, double a23,
	double a31, double a32, double a33
    ) {
	return
	    a11*det2x2(a22,a23,a32,a33)
	    -a21*det2x2(a12,a13,a32,a33)
	    +a31*det2x2(a12,a13,a22,a23);
    }

    inline double det4x4(
	double a11, double a12, double a13, double a14,
	double a21, double a22, double a23, double a24,               
	double a31, double a32, double a33, double a34,  
	double a41, double a42, double a43, double a44  
    ) {
	double m12 = a21*a12 - a11*a22;
	double m13 = a31*a12 - a11*a32;
	double m14 = a41*a12 - a11*a42;
	double m23 = a31*a22 - a21*a32;
	double m24 = a41*a22 - a21*a42;
	double m34 = a41*a32 - a31*a42;
	
	double m123 = m23*a13 - m13*a23 + m12*a33;
	double m124 = m24*a13 - m14*a23 + m12*a43;
	double m134 = m34*a13 - m14*a33 + m13*a43;
	double m234 = m34*a23 - m24*a33 + m23*a43;
	
	return (m234*a14 - m134*a24 + m124*a34 - m123*a44);
    }   

/******************************************************************************/

    enum ConvexCellFlag {
	None        = 0,
	WithVGlobal = 1,
	WithTFlags  = 2
    };

    typedef index_t ConvexCellFlags;
    
    /**
     * \brief Computes the intersection between a set of halfplanes using
     *  Bowyer-Watson algorithm.
     * \details Implementation does not use exact predicates, and does not
     *  climb from a random vertex. Do not use with a large number of planes.
     */
    class GEOGRAM_API ConvexCell {
      public:

	/**
	 * \brief ConvexCell constructor.
	 * \param[in] flags a combination of WithVGlobal, WithTFlags
	 */
	ConvexCell(ConvexCellFlags flags = None);

	/**
	 * \brief Removes all vertices and triangles from this
	 *  ConvexCell.
	 * \details Keeps allocated memory for future use.
	 */
	void clear();
	
	/**
	 * \brief Initializes this ConvexCell to an axis-aligned
	 *  box.
	 * \details Previous contents of this ConvexCell are 
	 *  discarded.
	 * \param[in] xmin , ymin , zmin , xmax , ymax , zmax
	 *  the coordinates of the box.
	 */
	void init_with_box(
	    double xmin, double ymin, double zmin,
	    double xmax, double ymax, double zmax
	);
	
	/**
	 * \brief Saves the computed cell in alias wavefront
	 *  file format.
	 * \param[in] filename the name of the file where to
	 *  save the cell.
	 * \param[in] shrink shrinking factor to ease visualization.
	 */
	void save(const std::string& filename, double shrink=0.0) const;


	/**
	 * \brief Saves the computed cell in alias wavefront 
	 *  file format.
	 * \param[out] out a stream where to save the output.
	 * \param[in] v_offset offset applied to vertex indices.
	 * \param[in] shrink shrinking factor to ease visualization.
	 * \return the number of created vertices. 
	 */
	index_t save(
	    std::ostream& out, index_t v_offset=1, double shrink=0.0
	) const;
	
	/**
	 * \brief Clips this convex cell by a new plane.
	 * \details The positive side of the plane equation corresponds to
	 *  what is kept. In other words, the normal vector P.x, P.y, P.z 
	 *  points towards the interior of this ConvexCell.
	 * \param[in] P the plane equation.
	 */
	void clip_by_plane(vec4 P);

	/**
	 * \brief Clips this convex cell by a new plane and stores
	 *  the corresponding global index in the newly created vertex.
	 * \details The positive side of the plane equation corresponds to
	 *  what is kept. In other words, the normal vector P.x, P.y, P.z 
	 *  points towards the interior of this ConvexCell.
	 *  This function can only be called if global indices are stored.
	 * \param[in] P the plane equation.
	 * \param[in] j the global index of the plane.
	 */
	void clip_by_plane(vec4 P, global_index_t j);
	
	/**
	 * \brief Gets the number of triangles.
	 * \return the number of created triangles.
	 * \details The created triangles are not
	 *  necessarily valid ones. To get the valid triangles,
	 *  one needs to traverse the list from first_valid_.
	 */
	index_t nb_t() const {
	    return nb_t_;
	}

	/**
	 * \brief Gets the number of vertices.
	 * \return the number of vertices.
	 * \details Some vertices can be incident to no triangle.
	 *  The first six vertices correspond to the facets of the
	 *  initial axis aligned box passed to the constructor.
	 */
	index_t nb_v() const {
	    return nb_v_;
	}

	/**
	 * \brief Directly creates a new vertex.
	 * \param[in] P the plane equation attached to the vertex.
	 * \return the index of the newly created vertex.
	 */
	index_t create_vertex(vec4 P) {
	    if(nb_v_ == max_v_) {
		grow_v();
	    }
	    plane_eqn_[nb_v_] = P;
	    index_t result = nb_v_;
	    ++nb_v_;
	    return result;
	}

	/**
	 * \brief Directly creates a new triangle.
	 * \param[in] i , j, k the three vertices of the
	 *  triangle.
	 * \details The triangle is inserted into the list
	 *  of valid triangles.
	 * \return the index of the newly created triangle.
	 */
	index_t create_triangle(index_t i, index_t j, index_t k) {
	    vbw_assert(i < nb_v());
	    vbw_assert(j < nb_v());
	    vbw_assert(k < nb_v());
	    return new_triangle(i,j,k);
	}

	/**
	 * \brief Replaces a vertex with the vertex at infinity
	 *  in all facets.
	 * \param[in] v the vertex to be killed.
	 */
	void kill_vertex(index_t v);

	/**
	 * \brief Tests whether a vertex has a corresponding
	 *  facet in the cell.
	 * \details One needs to call compute_geometry() before
	 *  calling this function.
	 */
	bool vertex_is_contributing(index_t v) const {
	    geo_assert(!geometry_dirty_);
	    return v2t_[v] != END_OF_LIST;
	}

       /**
	* \brief Gets a triangle incident to a vertex.
	* \param[in] v vertex index.
	* \return a triangle incident to v.
	*/
	index_t vertex_triangle(index_t v) const {
	    geo_assert(!geometry_dirty_);	    
	    return v2t_[v];
	}
	
	/**
	 * \brief Computes the geometry and some cached information.
	 * \details Needs to be called before volume(),
	 *   facet_area() and barycenter().
	 */
	void compute_geometry();

	/**
	 * \brief Gets the dual facet area of a given vertex.
	 * \details compute_geometry() needs to be called before.
	 * \param[in] v the vertex.
	 * \return the dual facet area associated with v.
	 * \details terminate() needs to be called before 
	 *  calling this function.
	 */
	double facet_area(index_t v) const;

	/**
	 * \brief Computes the volume of this convex cell.
	 * \details compute_geometry() needs to be called before.
	 * \return the volume.
	 */
	double volume() const;

	/**
	 * \brief Computes the barycenter of this convex cell.
	 * \details compute_geometry() needs to be called before.
	 * \return the barycenter.
	 */
	vec3 barycenter() const;

	/**
	 * \brief Computes the squared radius with respect to a 
	 *  center.
	 * \return the maximum squared distance between center and
	 *  all the vertices of the cell.
	 */
	double squared_radius(vec3 center) const;

	/**
	 * \brief Tests whether this ConvexCell is empty.
	 * \details ConvexCell can be empty if everything was
	 *  clipped out.
	 * \retval true if this ConvexCell is empty.
	 * \retval false otherwise.
	 */
	bool empty() const {
	    return first_valid_ == END_OF_LIST;
	}

	/**
	 * \brief Gets the global vertex index from a local 
	 *  vertex index.
	 * \details Vertex indices correspond to planes (remember,
	 *  we are in dual form).
	 * \param[in] lv the local vertex index
	 * \return the global vertex index that corresponds to
	 *  lv.
	 */
	global_index_t v_global_index(index_t lv) const {
	    vbw_assert(has_vglobal_);
	    vbw_assert(lv < nb_v());
	    return vglobal_[lv];
	}

	/**
	 * \brief Tests whether a vertex with a given global index
	 *  exists in this ConvexCell.
	 * \param[in] v the global index.
	 * \retval true if there exists in this ConvexCell a vertex with
	 *  global index \p v.
	 * \retval false otherwise.
	 */
	bool has_v_global_index(global_index_t v) const;

	/**
	 * \brief Gets the first triangle.
	 * \return the index of the first triangle, or END_OF_LIST
	 *  if this ConvexCell is empty.
	 */
	ushort first_triangle() const {
	    return ushort(first_valid_);
	}

	/**
	 * \brief Gets the successor of a triangle.
	 * \param[in] t the index of a valid triangle.
	 * \return the index of the successor of \p t, or END_OF_LIST
	 *  if \p t is the last triangle.
	 */
	ushort next_triangle(ushort t) const {
	    return get_triangle_flags(t);
	}

	/**
	 * \brief Gets the point that corresponds to a triangle.
	 * \details If compute_geometry() was called, this gets
	 *  the previously computed point, else it is computed
	 *  and returned.
	 * \param[in] t the index of the triangle.
	 * \return the point that corresponds to triangle \p t.
	 */
	vec3 triangle_point(ushort t) const {
	    if(geometry_dirty_) {
		TriangleWithFlags T = get_triangle_and_flags(t);
		vec4 result = compute_triangle_point(T);
		vbw_assert(result.w != 0.0);
		return make_vec3(
		    result.x/result.w, result.y/result.w, result.z/result.w
		);
	    }
	    return triangle_point_[t];
	}

	global_index_t triangle_v_global_index(ushort t, index_t llv) const {
	    Triangle T = get_triangle(t);
	    ushort lv = ushort((llv==0)*T.i + (llv==1)*T.j + (llv==2)*T.k);
	    return v_global_index(lv);
	}

	index_t triangle_v_local_index(ushort t, index_t llv) const {
	    Triangle T = get_triangle(t);
	    return index_t((llv==0)*T.i + (llv==1)*T.j + (llv==2)*T.k);
	}
	
	bool triangle_is_user_marked(ushort t) {
	    vbw_assert(has_tflags_);
	    vbw_assert(t < max_t_);
	    return (tflags_[t] != 0);
	}

	void triangle_user_mark(ushort t) {
	    vbw_assert(has_tflags_);
	    vbw_assert(t < max_t_);
	    tflags_[t] = 1;
	}

	void triangle_user_unmark(ushort t) {
	    vbw_assert(has_tflags_);
	    vbw_assert(t < max_t_);
	    tflags_[t] = 0;
	}

	/**
	 * \brief Tests whether a cell has at least one vertex in conflict with
	 *  a halfspace.
	 * \param[in] P the equation of the halfspace.
	 * \retval true if there exists a triangle t such that 
	 *  triangle_is_in_conflict(P)
	 * \retval false otherwise.
	 */
	bool cell_has_conflict(const vec4& P) {
	    for(
		ushort t = first_triangle();
		t!=END_OF_LIST; t=next_triangle(t)
	    ) {
		TriangleWithFlags T = get_triangle_and_flags(t);
		if(triangle_is_in_conflict(T,P)) {
		    return true;
		}
	    }
	    return false;
	}

	/**
	 * \brief Tests whether a cell has all its vertices in conflict
	 *  with a plane.
	 * \param[in] P the equation of the halfspace.
	 * \retval true if all the triangles are in conflict with P.
	 * \retval false otherwise.
	 */
	bool cell_is_totally_in_conflict(const vec4& P) {
	    for(
		ushort t = first_triangle();
		t!=END_OF_LIST; t=next_triangle(t)
	    ) {
		TriangleWithFlags T = get_triangle_and_flags(t);
		if(!triangle_is_in_conflict(T,P)) {
		    return false;
		}
	    }
	    return true;
	}


	
	/**
	 * \brief Gets a triangle adjacent to another triangle by edge
	 *  local index.
	 * \param[in] t a triangle.
	 * \param[in] le local index of an edge of \p t (in 0..2).
	 * \return the triangle adjacent to \p t along \ p e.
	 */
	 index_t triangle_adjacent(index_t t, index_t le) const {
	     vbw_assert(t < max_t());
	     vbw_assert(le < 3);
	     Triangle T = get_triangle(t);
	     index_t v1 = index_t((le == 0)*T.j + (le == 1)*T.k + (le == 2)*T.i);
	     index_t v2 = index_t((le == 0)*T.k + (le == 1)*T.i + (le == 2)*T.j);
	     vbw_assert(vv2t(v1,v2) == t);
	     vbw_assert(vv2t(v2,v1) != END_OF_LIST);
	     return vv2t(v2,v1);
	 }

	/**
	 * \brief Gets the local index of a vertex in a triangle.
	 * \param[in] t a triangle.
	 * \param[in] v a vertex index.
	 * \return the local index of \p v in \p t (in 0..2).
	 */
	index_t triangle_find_vertex(index_t t, index_t v) const {
	    vbw_assert(t < max_t());
	    Triangle T = get_triangle(t);
	    index_t result = index_t((T.j == v) + 2*(T.k == v));
	    return result;
	}

	/**
	 * \brief Tests whether a triangle is infinite.
	 * \param[in] t the triangle
	 * \retval true if t is incident to the vertex at
	 *  infinity.
	 * \retval false otherwise.
	 */
	bool triangle_is_infinite(index_t t) const {
	    vbw_assert(t < max_t());
	    Triangle T = get_triangle(t);
	    return (
		T.i == VERTEX_AT_INFINITY ||
		T.j == VERTEX_AT_INFINITY ||
		T.k == VERTEX_AT_INFINITY
	    );
	}
	
	/**
	 * \brief Gets the equation of a plane associated with a vertex.
	 * \details The first six equations correspond to the six 
	 *  facets of a cube.
	 * \param[in] v the local index of the vertex.
	 */
        vec4 vertex_plane(index_t v) const {
	    vbw_assert(v < max_v());
	    return plane_eqn_[v];
	}

	/**
	 * \brief Gets the normal to the plane associated with a vertex.
	 * \details The first six equations correspond to the six 
	 *  facets of a cube.
	 * \param[in] v the local index of the vertex.
	 */
	vec3 vertex_plane_normal(index_t v) const {
	    vbw_assert(v != VERTEX_AT_INFINITY);
	    vbw_assert(v < max_v());
	    return make_vec3(
		plane_eqn_[v].x,
		plane_eqn_[v].y,
		plane_eqn_[v].z
	    );
	}
	
	/**
	 * \brief Tests whether a triangle is marked as conflict.
	 * \param[in] t a triangle.
	 * \retval true if \p t is marked as conflict.
	 * \retval false otherwise.
	 */
	bool triangle_is_marked_as_conflict(index_t t) const {
	    vbw_assert(t < max_t());
	    return (get_triangle_flags(t) & ushort(CONFLICT_MASK)) != 0;
	}
       
	/**
	 * \brief Tests whether a triangle is in conflict with a plane.
	 * \details A triangle is in conflict with a plane if feeding the point
	 *  associated with the triangle in the equation of the plane yields 
	 *  a negative number.
	 * \param[in] T a triangle.
	 * \param[in] eqn the four coefficients of the equation of the plane.
	 * \retval true if \p t is in conflict with \p eqn.
	 * \retval false otherwise.
	 */
        bool triangle_is_in_conflict(
	    TriangleWithFlags T, const vec4& eqn
	) const;

	/**
	 * \brief Creates a new triangle.
	 * \param[in] i , j , k the three vertices of the triangle.
	 */
	index_t new_triangle(index_t i, index_t j, index_t k) {
	    index_t result = first_free_;
	    if(result == END_OF_LIST) {
		result = nb_t_;
		++nb_t_;
		if(nb_t_ > max_t()) {
		    grow_t();
		}
	    } else {
		first_free_ = index_t(
		    get_triangle_flags(first_free_) & ~ushort(CONFLICT_MASK)
		);
	    }
	    vbw_assert(result < max_t());
	    t_[result] = make_triangle_with_flags(
		ushort(i), ushort(j), ushort(k), ushort(first_valid_)
	    );
	    set_vv2t(i, j, result);
	    set_vv2t(j, k, result);
	    set_vv2t(k, i, result);	    
	    first_valid_ = result;
	    if(has_tflags_) {
		tflags_[result] = 0;
	    }
	    return result;
	}
	
	/**
	 * \brief Creates a new triangle.
	 * \details Adjacency information is not used (kept for reference).
	 * \param[in] i , j , k the three vertices of the triangle.
	 * \param[in] adj0 , adj1 , adj2 the three adjacent triangles
	 *  (unused in this version).
	 * \return the index of the new triangle.
	 */
	index_t new_triangle(
	    index_t i, index_t j, index_t k,
	    index_t adj0, index_t adj1, index_t adj2
	) {
	    // Silence warnings
	    (void)(adj0);
	    (void)(adj1);
	    (void)(adj2);	    
	    return new_triangle(i, j, k);
	}

	/**
	 * \brief Computes the coordinates of the point
	 *  associated with a triangle.
	 * \param[in] T the triangle.
	 * \return the intersection between the three planes
	 *  associated with the three vertices of the triangle,
	 *  in homogeneous coordinates.
	 */
	vec4 compute_triangle_point(TriangleWithFlags T) const;

	/**
	 * \brief Gets the three vertices of a triangle.
	 * \param[in] t the triangle.
	 * \return a Triangle with the indices of the three vertices
	 *  of the triangle.
	 */
	Triangle get_triangle(index_t t) const {
	    vbw_assert(t < max_t());
	    return t_[t];
	}
	
	/**
	 * \brief Gets the flags associated with a triangle.
	 * \details Contains both the conflict flag and the
	 *  chaining.
	 * \param[in] t the triangle.
	 * \return the flags associated with \p t.
	 */
	ushort get_triangle_flags(index_t t) const {
	    vbw_assert(t < max_t());
	    return t_[t].flags;
	}

	/**
	 * \brief Sets the flags of a triangle.
	 * \param[in] t the triangle.
	 * \param[in] flags the flags to be set.
	 */
	void set_triangle_flags(index_t t, ushort flags) {
	    vbw_assert(t < max_t());
	    t_[t].flags = flags;
	}

	/**
	 * \brief Gets the three vertices of a triangle and its flags.
	 * \param[in] t the triangle.
	 * \return a TriangleWithFlags with the indices of the three vertices
	 *  of the triangle and the flags.
	 */
	TriangleWithFlags get_triangle_and_flags(index_t t) const {
	    vbw_assert(t < max_t());
	    return t_[t];
	}
	
	/**
	 * \brief Gets the triangle incident to an oriented edge.
	 * \param[in] v1 , v2 the two vertices of the oriented edge.
	 * \return the triangle incident to the oriented edge.
	 */
	index_t vv2t(index_t v1, index_t v2) const {
	    vbw_assert(v1 < max_v());
	    vbw_assert(v2 < max_v());
	    return index_t(vv2t_[max_v_*v1 + v2]);
	}

	/**
	 * \brief Sets the triangle incident to an oriented edge.
	 * \param[in] v1 , v2 the two vertices of the oriented edge.
	 * \param[in] t the triangle incident to the oriented edge.
	 */
	void set_vv2t(index_t v1, index_t v2, index_t t) {
	    vbw_assert(v1 < max_v());
	    vbw_assert(v2 < max_v());
	    vv2t_[max_v_*v1+v2] = ushort(t);
	}

	/**
	 * \brief Gets the maximum valid index for a triangle.
	 * \return the maximum valid index of a triangle.
	 */
	index_t max_t() const {
	    return max_t_;
	}

	/**
	 * \brief Gets the maximum valid index for a vertex.
	 * \return the maximum valid index of a vertex.
	 */
	index_t max_v() const {
	    return max_v_;
	}

	/**
	 * \brief Allocates more space for triangles.
	 * \details Makes max_t_ twice bigger.
	 */
	void grow_t();

	/**
	 * \brief Allocates more space for vertices.
	 * \details Makes max_v_ twice bigger.
	 */
	void grow_v();

      private:

	/** \brief number of allocated triangles */
	index_t max_t_;

	/** \brief number of allocated vertices */
	index_t max_v_;

	/** \brief indices of triangle vertices and flags */
	vector<TriangleWithFlags> t_;

	/**
	 * \brief vertex,vertex -> triangle adjacent to 
	 *   oriented edge.
	 */
	vector<ushort> vv2t_;

	/**
	 * \brief plane equation attached to each vertex,
	 *  as specified by clip_by_plane().
	 */
	vector<vec4> plane_eqn_;

	/** \brief number of used triangles. */
	index_t nb_t_;

	/** \brief number of used vertices. */	
	index_t nb_v_;

	/** \brief Head of the linked list of free triangles. */
	index_t first_free_;

	/** \brief Head of the linked list of valid triangles. */
	index_t first_valid_;

	/** 
	 * \brief true if triangle_point_ and t2v_ are 
	 *  not up to date. 
	 */
	bool geometry_dirty_;

	/**
	 * \brief dual vertex attached to each triangle.
	 */
	vector<vec3> triangle_point_;

	/** 
	 * \brief One triangle incident to each vertex, 
	 *  or END_OF_LIST if there is no such triangle.
	 */
	vector<ushort> v2t_;

	/**
	 * \brief Optional vector of gloval vertex indices.
	 */
	vector<global_index_t> vglobal_;
	
	/**
	 * \brief True if global vertex indices are stored.
	 */
	bool has_vglobal_;

	/**
	 * \brief Optional flags attached to the triangles.
	 */
	vector<uchar> tflags_;
	
	/**
	 * \brief True if triangle flags are stored.
	 */
	bool has_tflags_;
    };
}

namespace GEO {
    using VBW::ConvexCell;
}

#endif
