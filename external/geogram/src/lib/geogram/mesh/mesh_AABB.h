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

#ifndef GEOGRAM_MESH_MESH_AABB
#define GEOGRAM_MESH_MESH_AABB

/**
 * \file mesh_AABB.h
 * \brief Axis Aligned Bounding Box trees for accelerating
 *  geometric queries that operate on a Mesh.
 */

#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>
#include <geogram/basic/geometry.h>

namespace GEO {

    /**
     * \brief Axis Aligned Bounding Box tree of mesh facets.
     * \details Used to quickly compute facet intersection and
     *  to locate the nearest facet from 3d query points.
     */
    class GEOGRAM_API MeshFacetsAABB {
    public:
        /**
         * \brief Creates the Axis Aligned Bounding Boxes tree.
         * \param[in] M the input mesh. It can be modified,
         *  and will be triangulated (if
         *  not already a triangular mesh). The facets are
         *  re-ordered (using Morton's order, see mesh_reorder()).
         * \param[in] reorder if not set, Morton re-ordering is
         *  skipped (but it means that mesh_reorder() was previously
         *  called else the algorithm will be pretty unefficient).
         * \pre M.facets.are_simplices()
         */
        MeshFacetsAABB(Mesh& M, bool reorder = true);


	/**
	 * \brief Gets the mesh.
	 * \return a const reference to the mesh.
	 */
	const Mesh& mesh() const {
	    return mesh_;
	}
	
        /**
         * \brief Computes all the pairs of intersecting facets.
         * \param[in] action ACTION::operator(index_t,index_t) is
         *  invoked of all pairs of facets that have overlapping
         *  bounding boxes. triangles_intersection() needs to be
         *  called to detect the actual intersections.
         * \tparam ACTION user action class, that needs to define
         * operator(index_t,index_t), where the two indices are
         * the indices each pair of triangles that have intersecting
         * bounding boxes.
         */
        template <class ACTION>
        void compute_facet_bbox_intersections(
            ACTION& action
        ) const {
            intersect_recursive(
                action,
                1, 0, mesh_.facets.nb(),
                1, 0, mesh_.facets.nb()
            );
        }


        /**
         * \brief Computes all the intersections between a given
         *  box and the bounding boxes of all the facets.
         * \param[in] action ACTION::operator(index_t) is
         *  invoked for all facets that have a bounding
         *  box that intersects \p box_in.
         * \tparam ACTION user action class, that needs to define
         * operator(index_t), where the parameter is the index
         * of the triangle that has its bounding box intersecting 
         * \p box_in.
         */
        template< class ACTION >
        void compute_bbox_facet_bbox_intersections(
            const Box& box_in,
            ACTION& action
        ) const {
            bbox_intersect_recursive(
                action, box_in, 1, 0, mesh_.facets.nb()
            );
        }
        
        /**
         * \brief Finds the nearest facet from an arbitrary 3d query point.
         * \param[in] p query point
         * \param[out] nearest_point nearest point on the surface
         * \param[out] sq_dist squared distance between p and the surface.
         * \return the index of the facet nearest to point p.
         */
        index_t nearest_facet(
            const vec3& p, vec3& nearest_point, double& sq_dist
        ) const {
            index_t nearest_facet;
            get_nearest_facet_hint(p, nearest_facet, nearest_point, sq_dist);
            nearest_facet_recursive(
                p,
                nearest_facet, nearest_point, sq_dist,
                1, 0, mesh_.facets.nb()
            );
            return nearest_facet;
        }

        /**
         * \brief Computes the nearest point and nearest facet from
         * a query point, using user-specified hint.
         *
         * \details The hint is specified as reasonable initial values of
         * (nearest_facet, nearest_point, sq_dist). If multiple queries
         * are done on a set of points that has spatial locality,
         * the hint can be the result of the previous call.
         *
         * \param[in] p query point
         * \param[in,out] nearest_facet the nearest facet so far,
         *   or NO_FACET if not known yet
         * \param[in,out] nearest_point a point in nearest_facet
         * \param[in,out] sq_dist squared distance between p and
         *    nearest_point
         * \note On entry, \p sq_dist needs to be equal to the squared
         *   distance between \p p and \p nearest_point (it is easy to
         *   forget to update it when calling it within a loop).
         */
        void nearest_facet_with_hint(
            const vec3& p,
            index_t& nearest_facet, vec3& nearest_point, double& sq_dist
        ) const {
            if(nearest_facet == NO_FACET) {
                get_nearest_facet_hint(
                    p, nearest_facet, nearest_point, sq_dist
                );                
            }
            nearest_facet_recursive(
                p,
                nearest_facet, nearest_point, sq_dist,
                1, 0, mesh_.facets.nb()
            );
        }

        /**
         * \brief Computes the distance between an arbitrary 3d query
         *  point and the surface.
         * \param[in] p query point
         * \return the squared distance between \p p and the surface.
         */
        double squared_distance(const vec3& p) const {
            vec3 nearest_point;
            double result;
            nearest_facet(p, nearest_point, result);
            return result;
        }

	/**
	 * \brief Tests whether this surface mesh has an intersection
	 *  with a segment.
	 * \param[in] q1 , q2 the two extremities of the segment.
	 * \retval true if there exists an intersection between [q1 , q2]
	 *  and a facet of the mesh.
	 * \retval false otherwise.
	 */
	bool segment_intersection(const vec3& q1, const vec3& q2) const;


	/**
	 * \brief Finds the intersection between a segment and a surface that
	 *  is nearest to the first extremity of the segment.
	 * \param[in] q1 , q2 the two extremities of the segment.
	 * \param[out] t if there was an intersection, it is t*q2 + (1-t)*q1
	 * \param[out] f the intersected nearest facet or index_t(-1) if there
	 *  was no intersection.
	 * \retval true if there exists at least an intersection between [q1 , q2]
	 *  and a facet of the mesh.
	 * \retval false otherwise.
	 */
	bool segment_nearest_intersection(
	    const vec3& q1, const vec3& q2, double& t, index_t& f
	) const;

	
    protected:


        /**
         * \brief Computes all the facets that have a bbox that
         *  intersects a given bbox in a sub-tree of the AABB tree.
         *
         * Note that the tree structure is completely implicit,
         *  therefore the bounds of the (continuous) facet indices
         *  sequences that correspond to the facets contained
         *  in the two nodes are sent as well as the node indices.
         *
         * \param[in] action ACTION::operator(index_t) is
         *  invoked for all facet that has a bounding box that
         *  overlaps \p box.
         * \param[in] node index of the first node of the AABB tree
         * \param[in] b index of the first facet in \p node
         * \param[in] e one position past the index of the last
         *  facet in \p node
         */
        template <class ACTION>
        void bbox_intersect_recursive(
            ACTION& action,
            const Box& box,
            index_t node, index_t b, index_t e
        ) const {
            geo_debug_assert(e != b);

            // Prune sub-tree that does not have intersection
            if(!bboxes_overlap(box, bboxes_[node])) {
                return;
            }

            // Leaf case
            if(e == b+1) {
                action(b);
                return;
            }

            // Recursion
            index_t m = b + (e - b) / 2;
            index_t node_l = 2 * node;
            index_t node_r = 2 * node + 1;

            bbox_intersect_recursive(action, box, node_l, b, m);
            bbox_intersect_recursive(action, box, node_r, m, e);
        }
        
        /**
         * \brief Computes all the pairs of intersecting facets
         *  for two sub-trees of the AABB tree.
         *
         * Note that the tree structure is completely implicit,
         *  therefore the bounds of the (continuous) facet indices
         *  sequences that correspond to the facets contained
         *  in the two nodes are sent as well as the node indices.
         *
         * \param[in] action ACTION::operator(index_t,index_t) is
         *  invoked of all pairs of facets that have overlapping
         *  bounding boxes.
         * \param[in] node1 index of the first node of the AABB tree
         * \param[in] b1 index of the first facet in \p node1
         * \param[in] e1 one position past the index of the last
         *  facet in \p node1
         * \param[in] node2 index of the second node of the AABB tree
         * \param[in] b2 index of the first facet in \p node2
         * \param[in] e2 one position past the index of the second
         *  facet in \p node2
         */
        template <class ACTION>
        void intersect_recursive(
            ACTION& action,
            index_t node1, index_t b1, index_t e1,
            index_t node2, index_t b2, index_t e2
        ) const {
            geo_debug_assert(e1 != b1);
            geo_debug_assert(e2 != b2);

            // Since we are intersecting the AABBTree with *itself*,
            // we can prune half of the cases by skipping the test
            // whenever node2's facet index interval is greated than
            // node1's facet index interval.
            if(e2 <= b1) {
                return;
            }

            // The acceleration is here:
            if(!bboxes_overlap(bboxes_[node1], bboxes_[node2])) {
                return;
            }

            // Simple case: leaf - leaf intersection.
            if(b1 + 1 == e1 && b2 + 1 == e2) {
                action(b1, b2);
                return;
            }

            // If node2 has more facets than node1, then
            //   intersect node2's two children with node1
            // else
            //   intersect node1's two children with node2
            if(e2 - b2 > e1 - b1) {
                index_t m2 = b2 + (e2 - b2) / 2;
                index_t node2_l = 2 * node2;
                index_t node2_r = 2 * node2 + 1;
                intersect_recursive(action, node1, b1, e1, node2_l, b2, m2);
                intersect_recursive(action, node1, b1, e1, node2_r, m2, e2);
            } else {
                index_t m1 = b1 + (e1 - b1) / 2;
                index_t node1_l = 2 * node1;
                index_t node1_r = 2 * node1 + 1;
                intersect_recursive(action, node1_l, b1, m1, node2, b2, e2);
                intersect_recursive(action, node1_r, m1, e1, node2, b2, e2);
            }
        }

        /**
         * \brief Computes a reasonable initialization for
         *  nearest facet search.
         *
         * \details A good initialization makes the algorithm faster,
         *  by allowing early pruning of subtrees that provably
         *  do not contain the nearest neighbor.
         *
         * \param[in] p query point
         * \param[out] nearest_facet a facet reasonably near p
         * \param[out] nearest_point a point in nearest_facet
         * \param[out] sq_dist squared distance between p and nearest_point
         */
        void get_nearest_facet_hint(
            const vec3& p,
            index_t& nearest_facet, vec3& nearest_point, double& sq_dist
        ) const;

        /**
         * \brief The recursive function used by the implementation
         *  of nearest_facet().
         *
         * \details The first call may use get_nearest_facet_hint()
         * to initialize nearest_facet, nearest_point and sq_dist,
         * as done in nearest_facet().
         *
         * \param[in] p query point
         * \param[in,out] nearest_facet the nearest facet so far,
         * \param[in,out] nearest_point a point in nearest_facet
         * \param[in,out] sq_dist squared distance between p and nearest_point
         * \param[in] n index of the current node in the AABB tree
         * \param[in] b index of the first facet in the subtree under node \p n
         * \param[in] e one position past the index of the last facet in the
         *  subtree under node \p n
         */
        void nearest_facet_recursive(
            const vec3& p,
            index_t& nearest_facet, vec3& nearest_point, double& sq_dist,
            index_t n, index_t b, index_t e
        ) const;

        /**
         * \brief The recursive function used by the implementation
         *  of segment_intersection()
	 * \param[in] q1 , q2 the segment
	 * \param[in] dirinv precomputed 1/(q2.x-q1.x), 1/(q2.y-q1.y), 1/(q2.z-q1.z)
         * \param[in] n index of the current node in the AABB tree
         * \param[in] b index of the first facet in the subtree under node \p n
         * \param[in] e one position past the index of the last facet in the
         *  subtree under node \p n
	 * \retval true if their was an intersection
	 * \retval false otherwise
	 */
	bool segment_intersection_recursive(
	    const vec3& q1, const vec3& q2, const vec3& dirinv,
	    index_t n, index_t b, index_t e
	) const;

        /**
         * \brief The recursive function used by the implementation
         *  of segment_nearest_intersection()
	 * \param[in] q1 , q2 the segment
	 * \param[in] dirinv precomputed 1/(q2.x-q1.x), 1/(q2.y-q1.y), 1/(q2.z-q1.z)
         * \param[in] n index of the current node in the AABB tree
         * \param[in] b index of the first facet in the subtree under node \p n
         * \param[in] e one position past the index of the last facet in the
         *  subtree under node \p n
	 * \param[in,out] t the coordinate along [q1,q2] of the nearest 
	 *   intersection so-far.
	 * \param[in,out] f the nearest intersected facet so-far.
	 */
	void segment_nearest_intersection_recursive(
	    const vec3& q1, const vec3& q2, const vec3& dirinv,
	    index_t n, index_t b, index_t e,
	    double& t, index_t& f
	) const;

	
    protected:
        vector<Box> bboxes_;
        Mesh& mesh_;
    };

    /***********************************************************************/

    /**
     * \brief Axis Aligned Bounding Box tree of mesh tetrahedra.
     * \details Used to quickly find the tetrahedron that contains
     *  a given 3d point.
     */
    class GEOGRAM_API MeshCellsAABB {
    public:

        /**
         * \brief Symbolic constant for indicating that there
         *  is no containing tetrahedron.
         * \see containing_tet()
         */
        static const index_t NO_TET = index_t(-1);
        
        /**
         * \brief Creates the Axis Aligned Bounding Boxes tree.
         * \param[in] M the input mesh. It can be modified,
         *  The cells are re-ordered (using Morton's order, see mesh_reorder()).
         * \param[in] reorder if not set, Morton re-ordering is
         *  skipped (but it means that mesh_reorder() was previously
         *  called else the algorithm will be pretty unefficient).
         */
        MeshCellsAABB(Mesh& M, bool reorder = true);

	/**
	 * \brief Gets the mesh.
	 * \return a const reference to the mesh.
	 */
	const Mesh& mesh() const {
	    return mesh_;
	}
	
        /**
         * \brief Finds the index of a tetrahedron that contains a query point
         * \param[in] p a const reference to the query point
         * \param[in] exact specifies whether exact predicates should be used
         * \return the index of one of the tetrahedra that contains \p p or
         *  NO_TET if \p p is outside the mesh.
         * \note The input mesh needs to be tetrahedralized. If the mesh has
         *   arbitrary cells, then one may use instead containing_boxes().
         */
        index_t containing_tet(const vec3& p, bool exact =true) const {
            geo_debug_assert(mesh_.cells.are_simplices());
            return containing_tet_recursive(
                p, exact, 1, 0, mesh_.cells.nb()
            );
        }


        /**
         * \brief Computes all the intersections between a given
         *  box and the bounding boxes of all the cells.
         * \param[in] action ACTION::operator(index_t) is
         *  invoked for all cells that have a bounding
         *  box that intersects \p box_in.
         * \tparam ACTION user action class, that needs to define
         * operator(index_t), where the parameter is the index
         * of the cell that has its bounding box intersecting 
         * \p box_in.
         */
        template< class ACTION >
        void compute_bbox_cell_bbox_intersections(
            const Box& box_in,
            ACTION& action
        ) const {
            bbox_intersect_recursive(
                action, box_in, 1, 0, mesh_.cells.nb()
            );
        }



        /**
         * \brief Finds all the cells such that their bounding
         *  box contain a point.
         * \param[in] action ACTION::operator(index_t) is
         *  invoked for all cells that have a bounding
         *  box that contains \p p.
         * \tparam ACTION user action class, that needs to define
         * operator(index_t), where the parameter is the index
         * of the cell that has its bounding box containing
         * \p p.
         */
        template< class ACTION >
        void containing_boxes(
            const vec3& p,
            ACTION& action
        ) const {
            containing_boxes_recursive(
                action, p, 1, 0, mesh_.cells.nb()
            );
        }
        
    protected:

        /**
         * \brief Computes all the cells that have a bbox that
         *  intersects a given bbox in a sub-tree of the AABB tree.
         *
         * Note that the tree structure is completely implicit,
         *  therefore the bounds of the (continuous) facet indices
         *  sequences that correspond to the facets contained
         *  in the two nodes are sent as well as the node indices.
         *
         * \param[in] action ACTION::operator(index_t) is
         *  invoked for all cells that has a bounding box that
         *  overlaps \p box.
         * \param[in] box the query box
         * \param[in] node index of the first node of the AABB tree
         * \param[in] b index of the first facet in \p node
         * \param[in] e one position past the index of the last
         *  facet in \p node
         */
        template <class ACTION>
        void bbox_intersect_recursive(
            ACTION& action,
            const Box& box,
            index_t node, index_t b, index_t e
        ) const {
            geo_debug_assert(e != b);

            // Prune sub-tree that does not have intersection
            if(!bboxes_overlap(box, bboxes_[node])) {
                return;
            }

            // Leaf case
            if(e == b+1) {
                action(b);
                return;
            }

            // Recursion
            index_t m = b + (e - b) / 2;
            index_t node_l = 2 * node;
            index_t node_r = 2 * node + 1;

            bbox_intersect_recursive(action, box, node_l, b, m);
            bbox_intersect_recursive(action, box, node_r, m, e);
        }

        
        /**
         * \brief The recursive function used by the implementation
         *  of containing_tet().
         * \param[in] p a const reference to the query point
         * \param[in] exact specifies whether exact predicates should be used
         * \param[in] n index of the current node in the AABB tree
         * \param[in] b index of the first tet in the subtree under node \p n
         * \param[in] e one position past the index of the last tet in the
         *  subtree under node \p n
         * \return the index of one of the tetrahedra that contains \p p, or
         *  NO_TET if \p p is outside the mesh.
         */
        index_t containing_tet_recursive(
            const vec3& p, bool exact, 
            index_t n, index_t b, index_t e
        ) const;


        /**
         * \brief Computes all the cells that have a bbox that
         *  contain a given point in a sub-tree of the AABB tree.
         *
         * Note that the tree structure is completely implicit,
         *  therefore the bounds of the (continuous) facet indices
         *  sequences that correspond to the facets contained
         *  in the two nodes are sent as well as the node indices.
         *
         * \param[in] action ACTION::operator(index_t) is
         *  invoked for all cells that has a bounding box that
         *  contains \p p.
         * \param[in] p a const reference to the query point
         * \param[in] node index of the first node of the AABB tree
         * \param[in] b index of the first facet in \p node
         * \param[in] e one position past the index of the last
         *  facet in \p node
         */
        template <class ACTION>
        void containing_bboxes_recursive(
            ACTION& action,
            const vec3& p,
            index_t node, index_t b, index_t e
        ) const {
            geo_debug_assert(e != b);

            // Prune sub-tree that does not have intersection            
            if(!bboxes_[node].contains(p)) {
                return;
            }

            // Leaf case
            if(e == b+1) {
                action(b);
                return;
            }

            // Recursion
            index_t m = b + (e - b) / 2;
            index_t node_l = 2 * node;
            index_t node_r = 2 * node + 1;

            containing_bbox_recursive(action, p, node_l, b, m);
            containing_bbox_recursive(action, p, node_r, m, e);
        }
        
        vector<Box> bboxes_;
        Mesh& mesh_;
    };
    
}

#endif

