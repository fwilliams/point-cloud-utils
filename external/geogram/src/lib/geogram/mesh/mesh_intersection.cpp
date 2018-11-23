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

#include <geogram/mesh/mesh_intersection.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_fill_holes.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/triangle_intersection.h>
#include <geogram/basic/logger.h>

namespace {

    using namespace GEO;

    /**
     * \brief Computes the intersection between two triangular facets in
     *  a mesh
     * \param[in] M the mesh
     * \param[in] f1 index of the first facet
     * \param[in] f2 index of the second facet
     * \param[out] sym symbolic representation of the intersection (if any)
     * \return true if facets \p f1 and \p f2 have an intersection, false
     *  otherwise
     */
    bool triangles_intersect(
        const Mesh& M, index_t f1, index_t f2,
        vector<TriangleIsect>& sym
    ) {
        geo_debug_assert(M.facets.nb_vertices(f1) == 3);
        geo_debug_assert(M.facets.nb_vertices(f2) == 3);
        index_t c1 = M.facets.corners_begin(f1);
        const vec3& p1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c1));
        const vec3& p2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c1 + 1));
        const vec3& p3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c1 + 2));
        index_t c2 = M.facets.corners_begin(f2);
        const vec3& q1 = Geom::mesh_vertex(M, M.facet_corners.vertex(c2));
        const vec3& q2 = Geom::mesh_vertex(M, M.facet_corners.vertex(c2 + 1));
        const vec3& q3 = Geom::mesh_vertex(M, M.facet_corners.vertex(c2 + 2));
        return triangles_intersections(p1, p2, p3, q1, q2, q3, sym);
    }

    /**
     * \brief Tests whether two facets are adjacent
     * \details Two facets are adjacents if they share an edge
     * \param[in] M the mesh
     * \param[in] f1 index of the first facet
     * \param[in] f2 index of the second facet
     * \return true if facets \p f1 and \p f2 share an edge, false
     *  otherwise
     */
    bool facets_are_adjacent(const Mesh& M, index_t f1, index_t f2) {
        if(f1 == f2) {
            return true;
        }
        for(index_t c = M.facets.corners_begin(f1);
            c != M.facets.corners_end(f1); ++c) {
            if(M.facet_corners.adjacent_facet(c) == f2) {
                return true;
            }
        }
        return false;
    }

    /**
     * \brief Action class for storing intersections when traversing
     *  a AABBTree.
     */
    class StoreIntersections {
    public:
        /**
         * \brief Constructs the StoreIntersections
         * \param[in] M the mesh
         * \param[out] has_isect the flag that indicates for each facet
         *  whether it has intersections
         */
        StoreIntersections(
            const Mesh& M, vector<index_t>& has_isect
        ) :
            M_(M),
            has_intersection_(has_isect) {
            has_intersection_.assign(M_.facets.nb(), 0);
        }

        /**
         * \brief Determines the intersections between two facets
         * \details It is a callback for AABBTree traversal
         * \param[in] f1 index of the first facet
         * \param[in] f2 index of the second facet
         */
        void operator() (index_t f1, index_t f2) {
            // TODO: if facets are adjacents, test for
            // coplanarity.
            if(
                !facets_are_adjacent(M_, f1, f2) &&
                f1 != f2 && triangles_intersect(M_, f1, f2, sym_)
            ) {
                has_intersection_[f1] = 1;
                has_intersection_[f2] = 1;
            }
        }

    private:
        const Mesh& M_;
        vector<index_t>& has_intersection_;
        vector<TriangleIsect> sym_;
    };

    /**
     * \brief Deletes the intersecting facets from a mesh
     * \param[in] M the mesh
     * \param[in] nb_neigh number of rings of facets to delete around
     *  the intersecting facets
     * \return the number of facets that were deleted
     */
    index_t remove_intersecting_facets(Mesh& M, index_t nb_neigh = 0) {
        geo_assert(M.vertices.dimension() >= 3);
        mesh_repair(M, MESH_REPAIR_DEFAULT);  // it repairs and triangulates.

        vector<index_t> has_intersection;
        StoreIntersections action(M, has_intersection);
        MeshFacetsAABB AABB(M);
        AABB.compute_facet_bbox_intersections(action);

        for(index_t i = 1; i <= nb_neigh; i++) {
            for(index_t f = 0; f < M.facets.nb(); f++) {
                if(has_intersection[f] == 0) {
                    for(index_t c = M.facets.corners_begin(f);
                        c < M.facets.corners_end(f); c++
                     ) {
                        index_t f2 = M.facet_corners.adjacent_facet(c);
                        if(f2 != NO_FACET && has_intersection[f2] == i) {
                            has_intersection[f] = i + 1;
                            break;
                        }
                    }
                }
            }
        }

        index_t count = 0;
        for(index_t f = 0; f < has_intersection.size(); f++) {
            if(has_intersection[f] != 0) {
                count++;
            }
        }

        if(count != 0) {
            M.facets.delete_elements(has_intersection);
            Logger::out("Intersect")
                << "Removed " << count << " facets"
                << std::endl;
        } 
        return count;
    }
}

/****************************************************************************/

namespace GEO {

    void mesh_remove_intersections(
        Mesh& M, index_t max_iter
    ) {
        fill_holes(M, Numeric::max_float64());

        const double A = Geom::mesh_area(M);
        const double min_component_area = 5.0 * 0.001 * A;
        const index_t min_component_facets = 3;

        for(index_t i = 0; i < max_iter; i++) {
            index_t count = remove_intersecting_facets(M);
            if(count == 0) {
                return;
            }
            // Needs to be done before hole filling (to have clean hole borders)
            remove_small_connected_components(
                M, min_component_area, min_component_facets
            );
            fill_holes(M, Numeric::max_float64());
            // Needs to be done after hole filling
            // (removing bridges may generate small connected components)
            remove_small_connected_components(
                M, min_component_area, min_component_facets
            );
        }

        for(index_t neigh_size = 1; neigh_size < 4; neigh_size++) {
            for(index_t i = 0; i < max_iter; i++) {
                index_t count = remove_intersecting_facets(M, neigh_size);
                if(count == 0) {
                    return;
                }
                // Needs to be done before hole filling
                // (to have clean hole borders)
                remove_small_connected_components(
                    M, min_component_area, min_component_facets
                );
                fill_holes(M, Numeric::max_float64());
                // Needs to be done after hole filling
                // (removing bridges may generate small connected components)
                remove_small_connected_components(
                    M, min_component_area, min_component_facets
                );
            }
        }
    }
}

