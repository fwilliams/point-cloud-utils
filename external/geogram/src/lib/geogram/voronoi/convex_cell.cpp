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

#include <geogram/voronoi/convex_cell.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

namespace VBW {

    ConvexCell::ConvexCell(ConvexCellFlags flags) :
	max_t_(64),
	max_v_(32),
	t_(max_t_),
	vv2t_(max_v_*max_v_),
	plane_eqn_(max_v_)
    {
	nb_t_ = 0;
	nb_v_ = 0;
	first_free_ = END_OF_LIST;
        first_valid_ = END_OF_LIST;
	geometry_dirty_ = true;
	has_vglobal_ = ((flags & WithVGlobal) != 0);
	if(has_vglobal_) {
	    vglobal_.assign(max_v_,index_t(-1));
	}
	has_tflags_  = ((flags & WithTFlags)  != 0);
	if(has_tflags_) {
	    tflags_.assign(max_t_,0);
	}
    }

    /***********************************************************************/

    void ConvexCell::clear() {
	nb_t_ = 0;
	nb_v_ = 0;
	first_free_ = END_OF_LIST;
        first_valid_ = END_OF_LIST;
	geometry_dirty_ = true;
#ifdef VBW_DEBUG
	// Initialize all triangle flags with something
	// different from VALID_TRIANGLE.
	for(index_t t=0; t<max_t(); ++t) {
	    set_triangle_flags(t, END_OF_LIST);
	}

	// Initializing edge table, just for
	// debugging purposes.
	/*
	for(index_t v1=0; v1<max_v(); ++v1) {
	    for(index_t v2=0; v2<max_v(); ++v2) {
		set_vv2t(v1, v2, END_OF_LIST);
	    }
	}*/
#endif
    }
    
    /***********************************************************************/

    void ConvexCell::init_with_box(
	double xmin, double ymin, double zmin,
	double xmax, double ymax, double zmax
    ) {
	clear();

	// The vertex at infinity.
 	plane_eqn_[0] = make_vec4(0,0,0,0);

	// Offset for the 6 bounding box plane equations.
	// Here they come first (offset is zero).
	index_t boff = 1;
	
	// The equations of the six faces of the bounding box.
	plane_eqn_[boff  ] = make_vec4( 1.0, 0.0, 0.0, -xmin);
	plane_eqn_[boff+1] = make_vec4(-1.0, 0.0, 0.0,  xmax);	
	plane_eqn_[boff+2] = make_vec4( 0.0, 1.0, 0.0, -ymin);
	plane_eqn_[boff+3] = make_vec4( 0.0,-1.0, 0.0,  ymax);	
	plane_eqn_[boff+4] = make_vec4( 0.0, 0.0, 1.0, -zmin);
	plane_eqn_[boff+5] = make_vec4( 0.0, 0.0,-1.0,  zmax);	

        //   Create the 8 triangles that correspond to the
        // 8 vertices of the bounding box.
 	//   (Unused) adjacency info. ----------------.
	//   Triangle vertices -.                     |
        //                      v                     v
	new_triangle(         boff+2,boff+5,boff+0, 1,4,2);
	new_triangle(         boff+5,boff+3,boff+0, 5,0,3);
	new_triangle(         boff+1,boff+5,boff+2, 0,6,3);
	new_triangle(         boff+5,boff+1,boff+3, 7,1,2);
	new_triangle(         boff+4,boff+2,boff+0, 0,5,6);
	new_triangle(         boff+4,boff+0,boff+3, 1,7,4);
	new_triangle(         boff+2,boff+4,boff+1, 7,2,4);
	new_triangle(         boff+4,boff+3,boff+1, 3,6,5);

	// We already created 6 vertices (for the 6 bounding box
	// plane equations) plus the vertex at infinity.
	nb_v_ = 7;

	geometry_dirty_ = true;
    }
    
    
    /***********************************************************************/

    void ConvexCell::save(const std::string& filename, double shrink) const {
	std::cerr << "====> Saving " << filename << std::endl;	
	std::ofstream out(filename.c_str());
	save(out, 1, shrink);
    }
    
    
    index_t ConvexCell::save(
	std::ostream& out, index_t v_offset, double shrink
    ) const {

	vec3 g = make_vec3(0.0, 0.0, 0.0);
	if(shrink != 0.0) {
	    const_cast<ConvexCell*>(this)->compute_geometry();
	    g = barycenter();
	}
	
	vector<index_t> v2t(nb_v(),index_t(-1));
	vector<index_t> t_index(nb_t(),index_t(-1));
	index_t nt=0;

	{
	    index_t t = first_valid_;
	    while(t != END_OF_LIST) { 
		TriangleWithFlags T = get_triangle_and_flags(t);
		vec4 p = compute_triangle_point(T);
		p.x /= p.w;
		p.y /= p.w;
		p.z /= p.w;
		p.w = 1.0;
		if(shrink != 0.0) {
		    p.x = shrink * g.x + (1.0 - shrink) * p.x;
		    p.y = shrink * g.y + (1.0 - shrink) * p.y;
		    p.z = shrink * g.z + (1.0 - shrink) * p.z;
		}
		out << "v " << p.x << " " << p.y << " " << p.z << std::endl;
		t_index[t] = nt;
		++nt;
		v2t[T.i] = t;
		v2t[T.j] = t;
		v2t[T.k] = t;
		t = index_t(T.flags);
	    }
	}
	
	for(index_t v=1; v<nb_v(); ++v) {
	    if(v2t[v] != index_t(-1)) {
		index_t t = v2t[v];
		out << "f ";
		do {
		    out << (t_index[t]+v_offset) << " ";
		    index_t lv = triangle_find_vertex(t,v);		   
		    t = triangle_adjacent(t, (lv + 1)%3);
		} while(t != v2t[v]);
		out << std::endl;
	    }
	}

	return nt;
    }

    /***********************************************************************/

    bool ConvexCell::has_v_global_index(global_index_t v) const {
	vbw_assert(has_vglobal_);
	for(index_t i=0; i<nb_v(); ++i) {
	    if(vglobal_[i] == v) {
		return true;
	    }
	}
	return false;
    }
    
    void ConvexCell::clip_by_plane(vec4 eqn, global_index_t j) {
	vbw_assert(has_vglobal_);
	clip_by_plane(eqn);
	vglobal_[nb_v()-1] = j;
    }
	
    void ConvexCell::clip_by_plane(vec4 eqn) {
	geometry_dirty_ = true;

	index_t lv = nb_v_;	
	if(lv == max_v()) {
	    grow_v();
	}
	plane_eqn_[lv] = eqn;
	vbw_assert(lv < max_v());
	++nb_v_;
	
	// Step 1: Find conflict zone and link conflicted triangles
	// (to recycle them in free list).

	index_t conflict_head = END_OF_LIST;
	index_t conflict_tail = END_OF_LIST;

        // Classify triangles, compute conflict list and valid list.
	// Note: This could be done by climbing from a random triangle,
	// but here we prefer complete linear scan for several reasons:
	//   - it is more robust to numerical errors (here we are not
	//     using exact predicates).
	//   - the code is simpler.
	//   - and more importantly, we got no more than a few tenths of
	//     vertices.

        index_t t = first_valid_;
        first_valid_ = END_OF_LIST;
        while(t != END_OF_LIST) { 
	    TriangleWithFlags T = get_triangle_and_flags(t);
	    if(triangle_is_in_conflict(T,eqn)) {
		set_triangle_flags(
		    t, ushort(conflict_head) | ushort(CONFLICT_MASK)
		);
		conflict_head = t;
		if(conflict_tail == END_OF_LIST) {
		    conflict_tail = t;
		}
	    } else {
		set_triangle_flags(t, ushort(first_valid_));
		first_valid_ = t;
	    }
	    t = index_t(T.flags);
	} 

	// Special case: no triangle in conflict.
	if(conflict_head == END_OF_LIST) {
	    return;
	}

	// Step 2: Triangulate cavity

        t = conflict_head;
        while(t != END_OF_LIST) { 
	    TriangleWithFlags T = get_triangle_and_flags(t);
	    
	    index_t adj1 = vv2t(T.j,T.i);
	    index_t adj2 = vv2t(T.k,T.j);
	    index_t adj3 = vv2t(T.i,T.k);	       
	    
	    ushort flg1 = get_triangle_flags(adj1);
	    ushort flg2 = get_triangle_flags(adj2);
	    ushort flg3 = get_triangle_flags(adj3);		
	    
	    if(	!(flg1 & ushort(CONFLICT_MASK)) ) {
		new_triangle(lv, T.i, T.j);
	    }

	    if(	!(flg2 & ushort(CONFLICT_MASK)) ) {
		new_triangle(lv, T.j, T.k);		
	    }

	    if(	!(flg3 & ushort(CONFLICT_MASK)) ) {
		new_triangle(lv, T.k, T.i);
	    }
	       
	    t = index_t(T.flags & ~ushort(CONFLICT_MASK));
	} 
	
	// Step 3: Recycle conflict list.
	
	set_triangle_flags(conflict_tail, ushort(first_free_));
	first_free_ = conflict_head;

    }

    /***********************************************************************/

    bool ConvexCell::triangle_is_in_conflict(
	TriangleWithFlags T, const vec4& eqn
    ) const {
	double det = 0.0;

	// If one of the vertices of the triangle is the
	// vertex at infinity, then the triangle is in conflict
	// with eqn if the oriented director vector of the intersection
	// between the two planes of the two other vertices
	// is opposite to the normal vector to eqn.

	if(T.i == VERTEX_AT_INFINITY) {
	    vec3 n2 = vertex_plane_normal(T.j);
	    vec3 n3 = vertex_plane_normal(T.k);
	    det = -det3x3(
		eqn.x, n2.x, n3.x, 
		eqn.y, n2.y, n3.y, 
		eqn.z, n2.z, n3.z
	    );
	} else if(T.j == VERTEX_AT_INFINITY) {
	    vec3 n3 = vertex_plane_normal(T.k);
	    vec3 n1 = vertex_plane_normal(T.i);
	    det = -det3x3(
		n1.x, eqn.x, n3.x,
		n1.y, eqn.y, n3.y,
		n1.z, eqn.z, n3.z
	    );
	} else if(T.k == VERTEX_AT_INFINITY) {
	    vec3 n1 = vertex_plane_normal(T.i);
	    vec3 n2 = vertex_plane_normal(T.j);
	    det = -det3x3(
		n1.x, n2.x, eqn.x,
		n1.y, n2.y, eqn.y,
		n1.z, n2.z, eqn.z
	    );
	} else {
	    
	    //   The triangle is in conflict with eqn if the 
	    // result of compute_triangle_point(t) injected in eqn 
	    // has a negative sign.
	    //   Examining the formula in compute_triangle_point(), this corresponds
	    // to (minus) the 4x4 determinant of the 4 plane equations
	    // developed w.r.t. the 4th column.
	    // (see Edelsbrunner - Simulation of Simplicity for similar examples
	    //  of computations).
	
	    vec4 p1 = vertex_plane(T.i);
	    vec4 p2 = vertex_plane(T.j);
	    vec4 p3 = vertex_plane(T.k);
	    det = det4x4(
		p1.x, p2.x, p3.x, eqn.x,
		p1.y, p2.y, p3.y, eqn.y,
		p1.z, p2.z, p3.z, eqn.z,
		p1.w, p2.w, p3.w, eqn.w
	    );
	}
	return (det > 0.0);
    }
    
    vec4 ConvexCell::compute_triangle_point(TriangleWithFlags T) const {

	double infinite_len = 16.0;
	
	// Special cases with one of the three vertices at infinity.
	if(T.i == VERTEX_AT_INFINITY) {
	    vec4 Pj = vertex_plane(T.j);
	    vec4 Pk = vertex_plane(T.k);
	    vec3 Njk = normalize(cross(
		make_vec3(Pj.x, Pj.y, Pj.z),
		make_vec3(Pk.x, Pk.y, Pk.z)
	    ));
	    index_t t_adj = vv2t(T.k, T.j);
	    vbw_assert(!triangle_is_infinite(t_adj));
	    TriangleWithFlags T_adj = t_[t_adj];
	    vec4 result = compute_triangle_point(T_adj);
	    result.x += result.w * Njk.x * infinite_len;
	    result.y += result.w * Njk.y * infinite_len;
	    result.z += result.w * Njk.z * infinite_len;
	    return result;
	} else if(T.j == VERTEX_AT_INFINITY) {
	    vec4 Pk = vertex_plane(T.k);
	    vec4 Pi = vertex_plane(T.i);
	    vec3 Nki = normalize(cross(
		make_vec3(Pk.x, Pk.y, Pk.z),
		make_vec3(Pi.x, Pi.y, Pi.z)
	    ));
	    index_t t_adj = vv2t(T.i, T.k);
	    vbw_assert(!triangle_is_infinite(t_adj));
	    TriangleWithFlags T_adj = t_[t_adj];
	    vec4 result = compute_triangle_point(T_adj);
	    result.x += result.w * Nki.x * infinite_len;
	    result.y += result.w * Nki.y * infinite_len;
	    result.z += result.w * Nki.z * infinite_len;
	    return result;
	} else if(T.k == VERTEX_AT_INFINITY) {
	    vec4 Pi = vertex_plane(T.i);
	    vec4 Pj = vertex_plane(T.j);
	    vec3 Nij = normalize(cross(
		make_vec3(Pi.x, Pi.y, Pi.z),
		make_vec3(Pj.x, Pj.y, Pj.z)
	    ));
	    index_t t_adj = vv2t(T.j, T.i);
	    vbw_assert(!triangle_is_infinite(t_adj));
	    TriangleWithFlags T_adj = t_[t_adj];
	    vec4 result = compute_triangle_point(T_adj);
	    result.x += result.w * Nij.x * infinite_len;
	    result.y += result.w * Nij.y * infinite_len;
	    result.z += result.w * Nij.z * infinite_len;
	    return result;
	}
	
        // Get the plane equations associated with each vertex of t
	vec4 pi1 = vertex_plane(T.i);
	vec4 pi2 = vertex_plane(T.j);
	vec4 pi3 = vertex_plane(T.k);

        // Find the intersection of the three planes using Kramer's formula.
	// (see Edelsbrunner - Simulation of Simplicity for other examples).
	// 
	// Kramer's formula: each component of the solution is obtained as
	//  the ratio of two determinants:
	//   - the determinant of the system where the ith column is replaced 
	//     with the rhs
	//   divided by:
	//   - the determinant of the system.
	// 
	// System of equations to be solved:
	// pi1.x * x + pi1.y * y + pi1.z * z = -pi1.w
	// pi2.x * x + pi2.y * y + pi2.z * z = -pi2.w
	// pi3.x * x + pi3.y * y + pi3.z * z = -pi3.w
	// 
	// Expression of the solution given by Cramer's formula:
	//     | -pi1.w   p1.y   pi1.z |   | pi1.x  pi1.y  pi1.z |
	// x = | -pi2.w   p2.y   pi2.z | / | pi2.x  pi2.y  pi2.z |
	//     | -pi3.w   p3.y   pi3.z |   | pi3.x  pi3.y  pi3.z |
       	// 
	//     |  pi1.x  -p1.w   pi1.z |   | pi1.x  pi1.y  pi1.z |
	// y = |  pi2.x  -p2.w   pi2.z | / | pi2.x  pi2.y  pi2.z |
	//     |  pi3.x  -p3.w   pi3.z |   | pi3.x  pi3.y  pi3.z |
	// 
	//     |  pi1.x   p1.y  -pi1.w |   | pi1.x  pi1.y  pi1.z |
	// z = |  pi2.x   p2.y  -pi2.w | / | pi2.x  pi2.y  pi2.z |
	//     |  pi3.x   p3.y  -pi3.w |   | pi3.x  pi3.y  pi3.z |
       
	vec4 result;

	result.x = -det3x3(
	    pi1.w, pi1.y, pi1.z,
	    pi2.w, pi2.y, pi2.z,
	    pi3.w, pi3.y, pi3.z
	);

	result.y = -det3x3(
	    pi1.x, pi1.w, pi1.z,
	    pi2.x, pi2.w, pi2.z,
	    pi3.x, pi3.w, pi3.z
	);
	
	result.z = -det3x3(
	    pi1.x, pi1.y, pi1.w,
	    pi2.x, pi2.y, pi2.w,
	    pi3.x, pi3.y, pi3.w
	);
	
	result.w = det3x3(
	    pi1.x, pi1.y, pi1.z,
	    pi2.x, pi2.y, pi2.z,
	    pi3.x, pi3.y, pi3.z
	);

	return result;
    }
    
    /***********************************************************************/

    void ConvexCell::grow_v() {
	vector<ushort> vv2t_new((max_v_*2)*(max_v_*2));
	for(index_t j=0; j<max_v_; ++j) {	
	    for(index_t i=0; i<max_v_; ++i) {
		vv2t_new[max_v_ * 2 * j + i] = vv2t_[max_v_ * j + i];
	    }
	}
	std::swap(vv2t_, vv2t_new);
	max_v_ *= 2;
	plane_eqn_.resize(max_v_);
	vglobal_.resize(max_v_, global_index_t(-1));
    }

    void ConvexCell::grow_t() {
	max_t_ *= 2;
	t_.resize(max_t_);
	if(has_tflags_) {
	    tflags_.resize(max_t_,0);
	}
    }

    /***********************************************************************/
    
    void ConvexCell::kill_vertex(index_t v) {
	for(index_t t=0; t<nb_t(); ++t) {
	    Triangle T = get_triangle(t);
	    if(T.i == v) {
		T.i = VERTEX_AT_INFINITY;
	    }
	    if(T.j == v) {
		T.j = VERTEX_AT_INFINITY;
	    }
	    if(T.k == v) {
		T.k = VERTEX_AT_INFINITY;
	    }
	    set_vv2t(T.i, T.j, t);
	    set_vv2t(T.j, T.k, t);
	    set_vv2t(T.k, T.i, t);
	    t_[t].i = T.i;
	    t_[t].j = T.j;
	    t_[t].k = T.k;	    
	}
    }

    /***********************************************************************/
    
    void ConvexCell::compute_geometry() {
	if(!geometry_dirty_) {
	    return;
	}

	triangle_point_.resize(nb_t());
	
	// Yes, need to do that with two
	// instructions: resize(nb_v(), END_OF_LIST)
	// does not reset the existing items to
	// END_OF LIST !! (had a hard time finding
	// this bug...)
	v2t_.resize(nb_v());
	v2t_.assign(nb_v(),END_OF_LIST);

        index_t t = first_valid_;
        while(t != END_OF_LIST) { 
	    TriangleWithFlags T = get_triangle_and_flags(t);
	    vec4 p = compute_triangle_point(T);
	    triangle_point_[t] = make_vec3(p.x/p.w, p.y/p.w, p.z/p.w);
	    v2t_[T.i] = ushort(t);
	    v2t_[T.j] = ushort(t);
	    v2t_[T.k] = ushort(t);
	    t = index_t(T.flags);
	}
	
	geometry_dirty_ = false;
    }

    inline double triangle_area(vec3 p1, vec3 p2, vec3 p3) {
	double Ux = p2.x - p1.x;
	double Uy = p2.y - p1.y;
	double Uz = p2.z - p1.z;
	double Vx = p3.x - p1.x;
	double Vy = p3.y - p1.y;
	double Vz = p3.z - p1.z;
	double Wx = Uy * Vz - Uz * Vy;
	double Wy = Uz * Vx - Ux * Vz;
	double Wz = Ux * Vy - Uy * Vx;
	return 0.5 * ::sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
    }
    
    double ConvexCell::facet_area(index_t v) const {
	vbw_assert(v < nb_v());
	vbw_assert(!geometry_dirty_);

	ushort t1t2[2];
	index_t cur=0;
	double result = 0.0;
	
	if(v2t_[v] != END_OF_LIST) {
	    index_t t = v2t_[v];
	    do {
		if(cur < 2) {
		    t1t2[cur] = ushort(t);
		} else {
		    result += triangle_area(
			triangle_point_[t1t2[0]],
			triangle_point_[t1t2[1]],
			triangle_point_[t]
		    );
		    t1t2[1] = ushort(t);
		}
		++cur;
		index_t lv = triangle_find_vertex(t,v);		   
		t = triangle_adjacent(t, (lv + 1)%3);
	    } while(t != v2t_[v]);
	}
	
	return result;
    }

    inline double tet_volume(vec3 p1, vec3 p2, vec3 p3, vec3 p4) {
	double Ux = p2.x - p1.x;
	double Uy = p2.y - p1.y;
	double Uz = p2.z - p1.z;
	
	double Vx = p3.x - p1.x;
	double Vy = p3.y - p1.y;
	double Vz = p3.z - p1.z;
	
	double Wx = p4.x - p1.x;
	double Wy = p4.y - p1.y;
	double Wz = p4.z - p1.z;

	double UVx = Uy * Vz - Uz * Vy;
	double UVy = Uz * Vx - Ux * Vz;
	double UVz = Ux * Vy - Uy * Vx;

	return ::fabs(
	    UVx * Wx + UVy * Wy + UVz * Wz
	) / 6.0;
    }
    
    double ConvexCell::volume() const {
	vbw_assert(!geometry_dirty_);
	double result = 0.0;

	ushort t_origin = END_OF_LIST;
	for(index_t v=0; v<nb_v_; ++v) {
	    if(v2t_[v] == END_OF_LIST) {
		continue;
	    }
	    if(t_origin == END_OF_LIST) {
		t_origin = v2t_[v];
		continue;
	    }
	    ushort t1t2[2];
	    index_t cur=0;
	    index_t t = v2t_[v];
	    do {
		if(cur < 2) {
		    t1t2[cur] = ushort(t);
		} else {
		    result += tet_volume(
			triangle_point_[t_origin],
			triangle_point_[t1t2[0]],
			triangle_point_[t1t2[1]],
			triangle_point_[t]
		    );
		    t1t2[1] = ushort(t);
		}
		++cur;
		index_t lv = triangle_find_vertex(t,v);		   
		t = triangle_adjacent(t, (lv + 1)%3);
	    } while(t != v2t_[v]);
	}
	return result;
    }

    vec3 ConvexCell::barycenter() const {
	vbw_assert(!geometry_dirty_);
	vec3 result = make_vec3(0.0, 0.0, 0.0);
	double m = 0.0;

	ushort t_origin = END_OF_LIST;
	for(index_t v=0; v<nb_v_; ++v) {
	    if(v2t_[v] == END_OF_LIST) {
		continue;
	    }
	    if(t_origin == END_OF_LIST) {
		t_origin = v2t_[v];
		continue;
	    }
	    ushort t1t2[2];
	    index_t cur=0;
	    index_t t = v2t_[v];
	    do {
		if(cur < 2) {
		    t1t2[cur] = ushort(t);
		} else {
		    vec3 p = triangle_point_[t_origin];
		    vec3 q = triangle_point_[t1t2[0]];
		    vec3 r = triangle_point_[t1t2[1]];
		    vec3 s = triangle_point_[t];
		    double cur_m = tet_volume(p,q,r,s);
		    m += cur_m;
		    result.x += cur_m*(p.x + q.x + r.x + s.x)/4.0;
		    result.y += cur_m*(p.y + q.y + r.y + s.y)/4.0;
		    result.z += cur_m*(p.z + q.z + r.z + s.z)/4.0;
		    t1t2[1] = ushort(t);
		}
		++cur;
		index_t lv = triangle_find_vertex(t,v);		   
		t = triangle_adjacent(t, (lv + 1)%3);
	    } while(t != v2t_[v]);
	}
	
	if(m != 0.0) {
	    result.x /= m;
	    result.y /= m;
	    result.z /= m;
	}
	return result;
    }
    
    /***********************************************************************/

    double ConvexCell::squared_radius(vec3 center) const {
	double result = 0.0;
        index_t t = first_valid_;
        while(t != END_OF_LIST) { 
	    TriangleWithFlags T = get_triangle_and_flags(t);
	    if(geometry_dirty_) {
		vec4 p4 = compute_triangle_point(T);
		vec3 p3 = make_vec3(
		    p4.x/p4.w, p4.y/p4.w, p4.z/p4.w
		);
		result = std::max(result, squared_distance(center,p3));		
	    } else {
		vec3 p = triangle_point_[t];
		result = std::max(result, squared_distance(center,p));
	    }
	    t = index_t(T.flags);
	}
	return result;
    }
    
}
