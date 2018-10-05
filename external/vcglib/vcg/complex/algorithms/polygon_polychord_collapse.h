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
#ifndef POLYGON_POLYCHORD_COLLAPSE_H
#define POLYGON_POLYCHORD_COLLAPSE_H

#include <vector>
#include <list>
#if __cplusplus >= 201103L
  #include <unordered_map>
#else
  #include <map>
#endif
#include <set>
#include <algorithm>
#include <iterator>
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/jumping_pos.h>

namespace vcg {
namespace tri {
/** \addtogroup trimesh */

/**
* @brief The PolychordCollapse class provides methods to semplify a quad mesh, by collapsing the polychords.
*
* This class is an implementation of a method very similar to that for mesh semplification proposed
* by Daniels et al. in "Quadrilateral mesh simplification", see http://www.cs.utah.edu/~jdaniels/research/asia2008_qms.htm
* The main function is PolychordCollapse::CollapsePolychord() which deletes all the quadrilateral faces in a polychord.
* The polychords that can be collapsed in this case are those forming a closed loop (a ring) or that start and end to
* mesh borders. A way to preserve the structure of the singularities is also provided.
* The convenient method PolychordCollapse::CollapseAllPolychords() finds and collapses all the polychords on a mesh.
* The input mesh should be polygonal, i.e. it should have the vcg::face::PolyInfo component. Even though a generic
* triangle mesh can be given, actually the class does not perform any collapsing operation since it sees only triangles,
* in fact it does not consider faux edges.
*/
template < typename PolyMeshType >
class PolychordCollapse {
public:
  typedef typename PolyMeshType::CoordType      CoordType;
  typedef typename PolyMeshType::VertexType     VertexType;
  typedef typename PolyMeshType::VertexPointer  VertexPointer;
  typedef typename PolyMeshType::VertexIterator VertexIterator;
  typedef typename PolyMeshType::FaceType       FaceType;
  typedef typename PolyMeshType::FacePointer    FacePointer;
  typedef typename PolyMeshType::FaceIterator   FaceIterator;

  /**
  * @brief The PC_ResultCode enum codifies the result type of a polychord collapse operation.
  */
  enum PC_ResultCode {
    PC_SUCCESS          = 0x000,
    PC_NOTMANIF         = 0x001,
    PC_NOTQUAD          = 0x002,
    PC_NOLINKCOND       = 0x004,
    PC_SINGSIDEA        = 0x008,
    PC_SINGSIDEB        = 0x010,
    PC_SINGBOTH         = 0x020,
    PC_SELFINTERSECT    = 0x040,
    PC_NOMOREMANIF      = 0x080,
    PC_VOID             = 0x100,
    PC_OTHER            = 0x100
  };

  /**
  * @brief The PC_Chord struct identifies a coord of a polychord passing through a quad.
  */
  struct PC_Chord {
    unsigned long mark;
    PC_ResultCode q;
    PC_Chord * prev;
    PC_Chord * next;
    PC_Chord() : mark(std::numeric_limits<unsigned long>::max()), q(PC_VOID), prev(NULL), next(NULL) { }
    inline void Reset() {
      mark = std::numeric_limits<unsigned long>::max();
      q = PC_VOID;
      prev = next = NULL;
    }
  };

  /**
  * @brief The PC_Chords class gives efficient access to each coord (relative to a face).
  */
  class PC_Chords {
  public:
    /**
     * @brief PC_Chords constructor.
     * @note Since each face corresponds to two chords, the actual size of the vector of chords is 2*mesh.face.size().
     * @param mesh
     */
    PC_Chords (const PolyMeshType &mesh) : _Chords(2*mesh.face.size()), _currentChord(NULL) {
      Reset(mesh);
    }

    /**
     * @brief ResetMarks
     */
    void ResetMarks() {
      for (size_t i = 0; i < _Chords.size(); ++i)
        _Chords.at(i).mark = std::numeric_limits<unsigned long>::max();
    }

    /**
     * @brief Reset rearrages the container.
     * @note Since each face corresponds to two chords, the actual size of the vector of chords is 2*mesh.face.size().
     * @param mesh
     */
    void Reset(const PolyMeshType &mesh) {
      _Chords.resize(2*mesh.face.size());
      for (size_t j = 0; j < _Chords.size(); ++j)
        _Chords[j].Reset();
      _currentChord = NULL;

      PC_Chord *chord = NULL;
      long long j = 0;
      for (size_t i = 0; i < _Chords.size(); ++i) {
        // set the prev
        chord = NULL;
        if ((long long)i-1 >= 0) {
          chord = &_Chords[i-1];
          if (vcg::tri::HasPerFaceFlags(mesh)) {
            j = i-1;
            while (j >= 0 && mesh.face[j/2].IsD())
              --j;
            if (j >= 0)
              chord = &_Chords[j];
            else
              chord = NULL;
          }
        }
        _Chords[i].prev = chord;

        // set the next
        chord = NULL;
        if (i+1 < _Chords.size()) {
          chord = &_Chords[i+1];
          if (vcg::tri::HasPerFaceFlags(mesh)) {
            j = i+1;
            while (j < (long long)_Chords.size() && mesh.face[j/2].IsD())
              ++j;
            if (j < (long long)_Chords.size())
              chord = &_Chords[j];
            else
              chord = NULL;
          }
        }
        _Chords[i].next = chord;
      }
      if (mesh.face.size() > 0) {
        // set the current coord (first - not deleted - face)
        _currentChord = &_Chords[0];
        if (vcg::tri::HasPerFaceFlags(mesh) && mesh.face[0].IsD())
          _currentChord = _currentChord->next;
      }
    }

    /**
     * @brief operator [], given a face index and an offset, it returns (a reference to) its corresponding PC_Chord.
     * @param face_edge A std::pair<size_t, unsigned char>(face_index, offset). The offset should be 0 or 1.
     * @return A reference to the corresponding PC_Chord.
     */
    inline PC_Chord & operator[] (const std::pair<size_t, unsigned char> &face_edge) {
      assert(face_edge.first >= 0 && 2*face_edge.first+face_edge.second < _Chords.size());
      return _Chords[2*face_edge.first + face_edge.second];
    }
    /**
     * @brief operator [], given a face index and an offset, it returns (a const reference to) its corresponding PC_Chord.
     * @param face_edge A std::pair<size_t, unsigned char>(face_index, offset). The offset should be 0 or 1.
     * @return A reference to the corresponding PC_Chord.
     */
    inline const PC_Chord & operator[] (const std::pair<size_t, unsigned char> &face_edge) const {
      assert(face_edge.first >= 0 && 2*face_edge.first+face_edge.second < _Chords.size());
      return _Chords[2*face_edge.first + face_edge.second];
    }

    /**
     * @brief operator [], given a coord, it returns its corresponding face index and edge.
     * @param coord The coord pointer.
     * @return A std::pair <size_t, unsigned char>(face_index, offset) with offset being 0 or 1.
     */
    inline std::pair<size_t, unsigned char> operator[] (PC_Chord const * const coord) {
      assert(coord >= &_Chords[0] && coord < &_Chords[0]+_Chords.size());
      return std::pair<size_t, unsigned char>((coord - &_Chords[0])/2, (coord - &_Chords[0])%2);
    }

    /**
     * @brief UpdateCoord updates the coord information and links.
     * @param coord The coord to update.
     * @param mark The mark of the polychord.
     * @param resultCode The code for the type of the polychord.
     */
    inline void UpdateCoord (PC_Chord &coord, const unsigned long mark, const PC_ResultCode resultCode) {
      // update prev and next
      if (coord.q == PC_VOID) {
        if (coord.prev != NULL && &coord != _currentChord)
          coord.prev->next = coord.next;
        if (coord.next != NULL && &coord != _currentChord)
          coord.next->prev = coord.prev;
      }
      coord.mark = mark;
      coord.q = resultCode;
    }

    /**
     * @brief Next, if it's not at the end, it goes to the next coord.
     */
    inline void Next () {
      if (_currentChord != NULL)
        _currentChord = _currentChord->next;
    }

    /**
     * @brief GetCurrent returns the current FaceType pointer and edge.
     * @param face_edge A std::pair where to store the FaceType pointer and the edge index.
     */
    inline void GetCurrent (std::pair<size_t, unsigned char> &face_edge) {
      if (_currentChord != NULL) {
        face_edge.first = (_currentChord - &_Chords[0])/2;
        face_edge.second = (_currentChord - &_Chords[0])%2;
      } else {
        face_edge.first = std::numeric_limits<size_t>::max();
        face_edge.second = 0;
      }
    }

    /**
     * @brief End says if an end has been reached.
     * @return true if an end has been reached, false otherwise.
     */
    inline bool End () {
      return _currentChord == NULL;
    }

  private:
    std::vector<PC_Chord>   _Chords;
    PC_Chord                *_currentChord;
  };

  /**
   * @brief The LinkCondition class provides a tool to check if a polychord satisfies the link conditions.
   */
  class LinkConditions {
  private:
    typedef long int                LCVertexIndex;
    typedef std::set<LCVertexIndex> LCVertexStar;   ///< define the star of a vertex
    typedef long int                LCEdgeIndex;
    typedef std::set<LCEdgeIndex>   LCEdgeStar;     ///< define the set of edges whose star involves a vertex

    /**
     * @brief The LCVertex struct represents a vertex for the Link Conditions.
     */
    struct LCVertex {
      LCVertexStar star;  // vertex star
      LCEdgeStar edges;   // list of edges whose star involves this vertex
      LCVertex(){}        // default constructor
      LCVertex(const LCVertex &lcVertex) {  // copy constructor
        star = lcVertex.star;
        edges = lcVertex.edges;
      }
      LCVertex & operator=(const LCVertex &lcVertex) { // assignment operator
        star = lcVertex.star;
        edges = lcVertex.edges;
        return *this;
      }
      void reset() { star.clear(); edges.clear(); } // reset
    };

    /**
     * @brief The LCEdge struct represents an edge for the Link Conditions.
     */
    struct LCEdge {
      LCVertexIndex v1, v2;       // endpoints
      LCVertexStar star;          // edge star
      LCEdge() {v1 = v2 = -1;}    // default contructor
      LCEdge(const LCEdge &lcEdge) {  // copy constructor
        v1 = lcEdge.v1;
        v2 = lcEdge.v2;
        star = lcEdge.star;
      }
      LCEdge & operator=(const LCEdge &lcEdge) {  // assignment operator
        v1 = lcEdge.v1;
        v2 = lcEdge.v2;
        star = lcEdge.star;
        return *this;
      }
      void reset() {  // reset
        v1 = -1;
        v2 = -1;
        star.clear();
      }
    };

  public:
    /**
     * @brief LinkCondition constructor.
     * @param size The number of vertices of the mesh.
     */
    LinkConditions (const size_t size) : _lcVertices(size) { }

    /**
     * @brief Resize just resets the size of the container.
     * @param size
     */
    inline void Resize(const size_t size) {
      _lcVertices.resize(size);
      LC_ResetStars();
    }

    /**
     * @brief CheckLinkConditions checks if collapsing the polychord starting from startPos
     * satisfies the link conditions.
     * @warning The polychord starts from startPos and ends to itself (if it's a loop) or to a border. In the latter case,
     * call this method starting from the opposite border of the strip of quads.
     * @param mesh The mesh for getting the vertex index.
     * @param startPos The starting position of the polychord.
     * @return true if satisfied, false otherwise.
     */
    bool CheckLinkConditions (const PolyMeshType &mesh, const vcg::face::Pos<FaceType> &startPos) {
      assert(!startPos.IsNull());
      assert(mesh.vert.size() == _lcVertices.size());
      std::vector<LCEdge> lcEdges;
      LCVertexStar intersection;

      // reset the stars
      LC_ResetStars();

      // compute the stars
      LC_computeStars(mesh, startPos, lcEdges);

      // for each edge e = (v1,v2)
      // if intersection( star(v1) , star(v2) ) == star(e)
      //      then collapse e
      // else
      //      return false (i.e. link conditions not satisfied)
      for (size_t e = 0; e < lcEdges.size(); e++) {
        // compute the intersetion
        intersection.clear();
        std::set_intersection(_lcVertices[lcEdges[e].v1].star.begin(), _lcVertices[lcEdges[e].v1].star.end(),
                              _lcVertices[lcEdges[e].v2].star.begin(), _lcVertices[lcEdges[e].v2].star.end(),
                              std::inserter(intersection, intersection.end()));

        // if intersection( star(v1) , star(v2) ) != star(e) then return false
        if (intersection != lcEdges[e].star)
            return false;

        // else simulate the collapse
        LC_SimulateEdgeCollapse(lcEdges, e);
      }
      // at this point all collapses are possible, thus return true
      return true;
    }

  private:
    /**
     * @brief LC_ResetStars resets the stars on a polychord.
     */
    void LC_ResetStars() {
      for (size_t v = 0; v < _lcVertices.size(); ++v)
        _lcVertices[v].reset();
    }

    /**
     * @brief LC_computeStars computes the stars of edges and vertices of the polychord from the starting pos
     * either to itself (if it's a loop) or to the border edge.
     * @param mesh The mesh for getting the vertex index.
     * @param startPos Starting position.
     * @param lcEdges Vector of edge stars.
     */
    void LC_computeStars (const PolyMeshType &mesh, const vcg::face::Pos<FaceType> &startPos, std::vector<LCEdge> &lcEdges) {
      assert(!startPos.IsNull());
      assert(mesh.vert.size() == _lcVertices.size());
      vcg::face::Pos<FaceType> runPos = startPos;
      vcg::face::JumpingPos<FaceType> vStarPos;
      vcg::face::Pos<FaceType> eStarPos;
      LCEdgeIndex edgeInd = -1;
      size_t nEdges = 0;

      // count how many edges
      do {
        ++nEdges;
        // go on the next edge
        runPos.FlipE();
        runPos.FlipV();
        runPos.FlipE();
        runPos.FlipF();
      } while (runPos != startPos && !runPos.IsBorder());
      if (runPos.IsBorder())
        ++nEdges;

      // resize the vector of edges
      lcEdges.resize(nEdges);
      for (size_t e = 0; e < nEdges; ++e)
        lcEdges[e].reset();

      /// compute the star of all the vertices and edges seen from the polychord
      runPos = startPos;
      do {
        // access the next lcedge
        edgeInd++;
        // set lcvertices references
        lcEdges[edgeInd].v1 = vcg::tri::Index(mesh, runPos.V());
        lcEdges[edgeInd].v2 = vcg::tri::Index(mesh, runPos.VFlip());
        // add this edge to its vertices edge-stars
        _lcVertices[lcEdges[edgeInd].v1].edges.insert(edgeInd);
        _lcVertices[lcEdges[edgeInd].v2].edges.insert(edgeInd);
        // compute the star of this edge
        lcEdges[edgeInd].star.insert(lcEdges[edgeInd].v1);  // its endpoints, clearly
        lcEdges[edgeInd].star.insert(lcEdges[edgeInd].v2);  // its endpoints, clearly
        // navigate over the other vertices of this facet
        eStarPos = runPos;
        eStarPos.FlipE();
        eStarPos.FlipV();
        while (eStarPos.V() != runPos.VFlip()) {
          // add current vertex to the star of this edge
          lcEdges[edgeInd].star.insert(vcg::tri::Index(mesh, eStarPos.V()));
          // add this edge to the edge-star of the current vertex
          _lcVertices[vcg::tri::Index(mesh, eStarPos.V())].edges.insert(edgeInd);
          // go on
          eStarPos.FlipE();
          eStarPos.FlipV();
        }
        // go on the opposite facet
        if (!runPos.IsBorder()) {
          eStarPos = runPos;
          eStarPos.FlipF();
          eStarPos.FlipE();
          eStarPos.FlipV();
          while (eStarPos.V() != runPos.VFlip()) {
            // add current vertex to the star of this edge
            lcEdges[edgeInd].star.insert(vcg::tri::Index(mesh, eStarPos.V()));
            // add this edge to the edge-star of the current vertex
            _lcVertices[vcg::tri::Index(mesh, eStarPos.V())].edges.insert(edgeInd);
            // go on
            eStarPos.FlipE();
            eStarPos.FlipV();
          }
        }

        // compute the star of vertex v2
        runPos.FlipV();
        vStarPos.Set(runPos.F(), runPos.E(), runPos.V());
        // v2 is in its star
        _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].star.insert(vcg::tri::Index(mesh, vStarPos.V()));
        do {
          vStarPos.FlipV();
          vStarPos.FlipE();
          while (vStarPos.V() != runPos.V()) {
            // add the current vertex to the v2 star
            _lcVertices[vcg::tri::Index(mesh, runPos.V())].star.insert(vcg::tri::Index(mesh, vStarPos.V()));
            // add v2 to the star of the current vertex
            _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].star.insert(vcg::tri::Index(mesh, runPos.V()));
            vStarPos.FlipV();
            vStarPos.FlipE();
          }
          vStarPos.NextFE();
        } while (vStarPos != runPos);

        // compute the star of vertex v1
        runPos.FlipV();
        vStarPos.Set(runPos.F(), runPos.E(), runPos.V());
        // v1 is in its star
        _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].star.insert(vcg::tri::Index(mesh, vStarPos.V()));
        do {
          vStarPos.FlipV();
          vStarPos.FlipE();
          while (vStarPos.V() != runPos.V()) {
            // add the current vertex to the v2 star
            _lcVertices[vcg::tri::Index(mesh, runPos.V())].star.insert(vcg::tri::Index(mesh, vStarPos.V()));
            // add v2 to the star of the current vertex
            _lcVertices[vcg::tri::Index(mesh, vStarPos.V())].star.insert(vcg::tri::Index(mesh, runPos.V()));
            vStarPos.FlipV();
            vStarPos.FlipE();
          }
          vStarPos.NextFE();
        } while (vStarPos != runPos);

        // when arrive to a border, stop
        if (runPos != startPos && runPos.IsBorder())
          break;

        // go on the next edge
        runPos.FlipE();
        runPos.FlipV();
        runPos.FlipE();
        runPos.FlipF();
      } while (runPos != startPos);

      // check if the starting pos or the border has been reached
      assert(runPos == startPos || runPos.IsBorder());
    }

    /**
     * @brief LC_SimulateEdgeCollapse simulates an edge collapse by updating the stars involved.
     * @param lcEdges The vector of edges.
     * @param edgeInd The in dex of the edge to collapse.
     */
    void LC_SimulateEdgeCollapse (std::vector<LCEdge> &lcEdges, const LCEdgeIndex edgeInd) {
      // let v1 and v2 be the two end points
      LCVertexIndex v1 = lcEdges[edgeInd].v1;
      LCVertexIndex v2 = lcEdges[edgeInd].v2;
      LCVertexIndex v = -1;

      /// v2 merges into v1:
      // star(v1) = star(v1) U star(v2)
      _lcVertices[v1].star.insert(_lcVertices[v2].star.begin(), _lcVertices[v2].star.end());
      _lcVertices[v1].star.erase(v2);     // remove v2 from v1-star
      _lcVertices[v2].star.erase(v1);     // remove v1 from v2-star
      // foreach v | v2 \in star(v) [i.e. v \in star(v2)]
      //      star(v) = star(v) U {v1} \ {v2}
      for (typename LCVertexStar::iterator vIt = _lcVertices[v2].star.begin(); vIt != _lcVertices[v2].star.end(); ++vIt) {
        v = *vIt;
        if (v == v2)  // skip v2 itself
          continue;
        _lcVertices[v].star.insert(v1);
        _lcVertices[v].star.erase(v2);
      }
      /// update the star of the edges which include v1 and v2 in their star
      // foreach e | v1 \in star(e) ^ v2 \in star(e)
      //      star(e) = star(e) \ {v1,v2} U {v1}
      for (typename LCEdgeStar::iterator eIt = _lcVertices[v1].edges.begin(); eIt != _lcVertices[v1].edges.end(); ++eIt)
        lcEdges[*eIt].star.erase(v2);
      for (typename LCEdgeStar::iterator eIt = _lcVertices[v2].edges.begin(); eIt != _lcVertices[v2].edges.end(); ++eIt) {
        lcEdges[*eIt].star.erase(v2);
        lcEdges[*eIt].star.insert(v1);
      }
    }

    /**
     * @brief _lcVertices is a vector of vertex stars for the link conditions.
     */
    std::vector<LCVertex> _lcVertices;
  };


  // PolychordCollapse's methods begin here::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//  /**
//   * @brief CheckConsistent checks for consistency. ONLY FOR DEBUG.
//   * @param mesh
//   * @return
//   */
//  static bool CheckConsistent(PolyMeshType &mesh) {
//    vcg::tri::RequirePerFaceFlags(mesh);
//    vcg::tri::RequirePerFaceColor(mesh);
//    for (size_t f = 0; f < mesh.face.size(); ++f) {
//      if (!mesh.face[f].IsD()) {
//        for (int v = 0; v < mesh.face[f].VN(); ++v) {
//          if (!vcg::face::IsBorder(mesh.face[f], v)) {
//            if (mesh.face[f].FFp(v)->IsD()) {
//              mesh.face[f].C() = vcg::Color4b(vcg::Color4b::Magenta);
//              return false;
//            }
//            if (mesh.face[f].FFp(v)->FFp(mesh.face[f].FFi(v)) != &mesh.face[f]) {
//              mesh.face[f].C() = vcg::Color4b(vcg::Color4b::Yellow);
//              return false;
//            }
//          }
//        }
//      }
//    }
//    return true;
//  }

  /**
   * @brief MarkPolychords marks the chords of the polychord starting at startPos.
   * @param mesh The input mesh.
   * @param startPos The starting position.
   * @param chords The vector of chords.
   * @param mark The current mark, used to identify quads already visited.
   */
  static void MarkPolychords(const PolyMeshType &mesh,
                             const vcg::face::Pos<FaceType> &startPos,
                             PC_Chords &chords,
                             const unsigned long mark) {
    vcg::face::Pos<FaceType> runPos = startPos;
    std::pair<size_t, unsigned char> face_edge(std::numeric_limits<size_t>::max(), 0);
    do {
      assert(runPos.F()->VN() == 4);
      face_edge.first = vcg::tri::Index(mesh, runPos.F());
      face_edge.second = runPos.E() % 2;
      chords[face_edge].mark = mark;
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (runPos != startPos && !runPos.IsBorder());
  }

  /**
   * @brief CollapsePolychord performs all checks and then collapses the polychord.
   *
   * @warning This function deletes faces and vertices by calling
   * vcg::tri::Allocator<PolyMeshType>::DeleteFace() and
   * vcg::tri::Allocator<PolyMeshType>::DeleteVertex().
   * The object PC_Chords chords is used to track the polychords, and it has got
   * a size proportional to that of the mesh face container. If you actually
   * delete faces and vertices by calling vcg::tri::Allocator<PolyMeshType>::CompactFaceVector()
   * and vcg::tri::Allocator<PolyMeshType>::CompactVertexVector() after this function,
   * object PC_Chords chords then is not valid any more, so you MUST rearrange it
   * by calling PC_Chords.Reset(). For the same reason, you MUST rearrange LinkConditions linkConditions
   * by calling LinkConditions.Resize().
   * However, for efficiency, you SHOULD compact vertex and face containers at the end of all your
   * polychord collapsing operations, without having to rearrange chords and linkConditions.
   * The function CollapseAllPolychords() does this for you.
   *
   * @note Vertex flags, face flags, FF adjacency and FV adjacency are required. Not anything else.
   * Such components are automatically updated here. If the mesh has other components that may be
   * affected by this editing, you should update them later by yourself.
   *
   * @param mesh The polygonal mesh used for getting the face index and deleting the faces
   * (it SHOULD have the vcg::face::PolyInfo component).
   * @param pos Position of the polychord.
   * @param mark Mark for the current polychord.
   * @param chords Vector of chords.
   * @param linkConditions Link conditions checker.
   * @param checkSing true if singularities on both sides are not allowed.
   * @return A PC_ResultCode resulting from checks or PC_SUCCESS if the collapse has been performed.
   */
  static PC_ResultCode CollapsePolychord (PolyMeshType &mesh,
                                          const vcg::face::Pos<FaceType> &pos,
                                          const unsigned long mark,
                                          PC_Chords &chords,
                                          LinkConditions &linkConditions,
                                          const bool checkSing = true) {
    vcg::tri::RequireFFAdjacency(mesh);

    if (mesh.IsEmpty())
      return PC_VOID;

    if (pos.IsNull())
      return PC_VOID;

    vcg::face::Pos<FaceType> tempPos, startPos;

    // check if the sequence of facets is a polychord and find the starting coord
    PC_ResultCode resultCode = CheckPolychordFindStartPosition(pos, startPos, checkSing);
    // if not successful, visit the sequence for marking it and return
    if (resultCode != PC_SUCCESS && resultCode != PC_SINGSIDEA && resultCode != PC_SINGSIDEB) {
      // if not manifold, visit the entire polychord ending on the non-manifold edge
      if (resultCode == PC_NOTMANIF) {
        tempPos = pos;
        VisitPolychord(mesh, tempPos, chords, mark, resultCode);
        if (tempPos.IsManifold() && !tempPos.IsBorder()) {
          tempPos.FlipF();
          VisitPolychord(mesh, tempPos, chords, mark, resultCode);
        }
        return resultCode;
      }
      // if not quad, visit all the polychords passing through this coord
      if (resultCode == PC_NOTQUAD) {
        tempPos = startPos;
        do {
          if (!tempPos.IsBorder()) {
            tempPos.FlipF();
            VisitPolychord(mesh, tempPos, chords, mark, resultCode);
            tempPos.FlipF();
          }
          tempPos.FlipV();
          tempPos.FlipE();
        } while (tempPos != startPos);
        VisitPolychord(mesh, startPos, chords, mark, resultCode);
        return resultCode;
      }
      VisitPolychord(mesh, startPos, chords, mark, resultCode);
      return resultCode;
    }
    // check if the link conditions are satisfied
    // if not satisfied, visit the sequence for marking it and return
    if (!linkConditions.CheckLinkConditions(mesh, startPos)) {
      VisitPolychord(mesh, startPos, chords, mark, PC_NOLINKCOND);
      return PC_NOLINKCOND;
    }
    // mark the polychord's chords
    MarkPolychords(mesh, startPos, chords, mark);
    // check if the polychord does not intersect itself
    // if it self-intersects, visit the polychord for marking it and return
    if (IsPolychordSelfIntersecting(mesh, startPos, chords, mark)) {
      VisitPolychord(mesh, startPos, chords, mark, PC_SELFINTERSECT);
      return PC_SELFINTERSECT;
    }
    // check if manifoldness remains
    // if it will loose manifoldness, visit the sequence for marking it and return
    if (!WillPolychordBeManifold(mesh, startPos, chords, mark)) {
      VisitPolychord(mesh, startPos, chords, mark, PC_NOMOREMANIF);
      return PC_NOMOREMANIF;
    }
    // at this point the polychord is collapsable, visit it for marking
    VisitPolychord(mesh, startPos, chords, mark, PC_SUCCESS);

    // now collapse
    CoordType point;
//    int valenceA = 0, valenceB = 0;
    vcg::face::Pos<FaceType> runPos = startPos;
    vcg::face::JumpingPos<FaceType> tmpPos;
//    bool onSideA = false, onSideB = false;
    vcg::face::Pos<FaceType> sideA, sideB;
    typedef std::queue<VertexPointer *> FacesVertex;
    typedef std::pair<VertexPointer, FacesVertex> FacesVertexPair;
    typedef std::queue<FacesVertexPair> FacesVertexPairQueue;
    FacesVertexPairQueue vQueue;
    typedef std::pair<FacePointer *, FacePointer> FFpPair;
    typedef std::pair<char *, char> FFiPair;
    typedef std::pair<FFpPair, FFiPair> FFPair;
    typedef std::queue<FFPair> FFQueue;
    FFQueue ffQueue;
    std::queue<VertexPointer> verticesToDeleteQueue;
    std::queue<FacePointer> facesToDeleteQueue;

    runPos = startPos;
    do {
      // compute new vertex
      point = (runPos.V()->P() + runPos.VFlip()->P()) / 2.f;
      if (checkSing) {
        if (resultCode == PC_SINGSIDEA)
          point = runPos.V()->P();
        else if (resultCode == PC_SINGSIDEB)
          point = runPos.VFlip()->P();
      }
      runPos.V()->P() = point;
      // list the vertex pointer of the faces on the other side to be updated
      vQueue.push(FacesVertexPair());
      vQueue.back().first = runPos.V();
      tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
      tmpPos.FlipV();
      tmpPos.NextFE();    // go to next face
      while (tmpPos.F() != runPos.F()) {
        if (tmpPos.F() != runPos.FFlip())
          vQueue.back().second.push(&tmpPos.F()->V(tmpPos.VInd()));
        tmpPos.NextFE();    // go to next face
      }

      // enqueue to delete the other vertex
      verticesToDeleteQueue.push(runPos.VFlip());

      // list the adjacencies
      sideA = runPos;
      sideA.FlipE();
      sideA.FlipF();
      sideB = runPos;
      sideB.FlipV();
      sideB.FlipE();
      sideB.FlipF();
      // first side
      if (!sideA.IsBorder()) {
        ffQueue.push(FFPair(FFpPair(),FFiPair()));
        ffQueue.back().first.first = &sideA.F()->FFp(sideA.E());
        ffQueue.back().second.first = &sideA.F()->FFi(sideA.E());
        if (!sideB.IsBorder()) {
          ffQueue.back().first.second = sideB.F();
          ffQueue.back().second.second = sideB.E();
        } else {
          ffQueue.back().first.second = sideA.F();
          ffQueue.back().second.second = sideA.E();
        }
      }
      // second side
      if (!sideB.IsBorder()) {
        ffQueue.push(FFPair(FFpPair(),FFiPair()));
        ffQueue.back().first.first = &sideB.F()->FFp(sideB.E());
        ffQueue.back().second.first = &sideB.F()->FFi(sideB.E());
        if (!sideA.IsBorder()) {
          ffQueue.back().first.second = sideA.F();
          ffQueue.back().second.second = sideA.E();
        } else {
          ffQueue.back().first.second = sideB.F();
          ffQueue.back().second.second = sideB.E();
        }
      }

      // enqueue to delete the face
      facesToDeleteQueue.push(runPos.F());

      // go on next edge/face
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (runPos != startPos && !runPos.IsBorder());
    assert(runPos == startPos || vcg::face::IsBorder(*startPos.F(),startPos.E()));
    if (runPos.IsBorder()) {
      // compute new vertex on the last (border) edge
      point = (runPos.V()->P() + runPos.VFlip()->P()) / 2.f;
      if (checkSing) {
        if (resultCode == PC_SINGSIDEA)
          point = runPos.V()->P();
        else if (resultCode == PC_SINGSIDEB)
          point = runPos.VFlip()->P();
      }
      runPos.V()->P() = point;
      // list the vertex pointer of the faces on the other side to be updated
      vQueue.push(FacesVertexPair());
      vQueue.back().first = runPos.V();
      tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
      tmpPos.FlipV();
      tmpPos.NextFE();    // go to next face
      while (tmpPos.F() != runPos.F()) {
        vQueue.back().second.push(&tmpPos.F()->V(tmpPos.VInd()));
        tmpPos.NextFE();
      }

      // enqueue to delete the other vertex
      verticesToDeleteQueue.push(runPos.VFlip());
    }

    // update vertices
    while (!vQueue.empty()) {
      while (!vQueue.front().second.empty()) {
        *vQueue.front().second.front() = vQueue.front().first;
        vQueue.front().second.pop();
      }
      vQueue.pop();
    }

    // update adjacencies
    while (!ffQueue.empty()) {
      *ffQueue.front().first.first = ffQueue.front().first.second;
      *ffQueue.front().second.first = ffQueue.front().second.second;
      ffQueue.pop();
    }

    // delete faces
    while (!facesToDeleteQueue.empty()) {
      vcg::tri::Allocator<PolyMeshType>::DeleteFace(mesh, *facesToDeleteQueue.front());
      facesToDeleteQueue.pop();
    }

    // delete vertices
    while (!verticesToDeleteQueue.empty()) {
      vcg::tri::Allocator<PolyMeshType>::DeleteVertex(mesh, *verticesToDeleteQueue.front());
      verticesToDeleteQueue.pop();
    }

    return PC_SUCCESS;
  }

  /**
   * @brief CollapseAllPolychords finds and collapses all the polychords.
   * @param mesh The input polygonal mesh (it SHOULD have the vcg::face::PolyInfo component).
   * @param checkSing true if singularities on both sides of a polychord are not allowed.
   */
  static void CollapseAllPolychords (PolyMeshType &mesh, const bool checkSing = true) {
    vcg::tri::RequireFFAdjacency(mesh);

    if (mesh.IsEmpty())
      return;

    vcg::face::Pos<FaceType> pos;
    PC_ResultCode resultCode;
    std::pair<size_t, unsigned char> face_edge;
    // construct the link conditions checker
    LinkConditions linkConditions(mesh.vert.size());
    // construct the vector of chords
    PC_Chords chords(mesh);
    unsigned long mark = 0;

    // iterate over all the chords
    while (!chords.End()) {
      // get the current coord
      chords.GetCurrent(face_edge);
      resultCode = chords[face_edge].q;
      assert(resultCode == PC_VOID);
      // construct a pos on the face and edge of the current coord
      pos.Set(&mesh.face[face_edge.first], face_edge.second, mesh.face[face_edge.first].V(face_edge.second));
      // (try to) collapse the polychord
      resultCode = CollapsePolychord(mesh, pos, mark, chords, linkConditions, checkSing);
      // go to the next coord
      chords.Next();

      // increment the mark
      ++mark;
      if (mark == std::numeric_limits<unsigned long>::max()) {
        chords.ResetMarks();
        mark = 0;
      }
    }
  }

  /**
   * @brief FindPolychords lists all the valid polychords starting position of a mesh.
   * @param mesh The input mesh.
   * @param polychords The container of results.
   * @param loopsOnly true if closed polychords only must be listed, false for all polychords.
   */
  static void FindPolychords (PolyMeshType &mesh, std::deque< vcg::face::Pos<FaceType> > &polychords, const bool loopsOnly = false) {
    vcg::tri::RequireFFAdjacency(mesh);

    polychords.clear();

    if (mesh.IsEmpty())
      return;

    vcg::face::Pos<FaceType> pos, startPos;
    PC_ResultCode resultCode;
    std::pair<size_t, unsigned char> face_edge;
    // construct the vector of chords
    PC_Chords chords(mesh);
    unsigned long mark = 0;

    // iterate over all the chords
    while (!chords.End()) {
      // get the current coord
      chords.GetCurrent(face_edge);
      // construct a pos on the face and edge of the current coord
      pos.Set(&mesh.face[face_edge.first], face_edge.second, mesh.face[face_edge.first].V(face_edge.second));

      // check and find start pos
      resultCode = CheckPolychordFindStartPosition(pos, startPos, false);
      // visit the polychord
      if (resultCode == PC_SUCCESS || resultCode == PC_SINGBOTH || resultCode == PC_SINGSIDEA || resultCode == PC_SINGSIDEB) {
        VisitPolychord(mesh, startPos, chords, mark, PC_OTHER);
        // store a new polychord
        if (!loopsOnly)
          polychords.push_back(startPos);
        else if (!startPos.IsBorder())
          polychords.push_back(startPos);
      } else {
        if (resultCode == PC_NOTMANIF) {
          pos = startPos;
          VisitPolychord(mesh, pos, chords, mark, resultCode);
          if (pos.IsManifold() && !pos.IsBorder()) {
            pos.FlipF();
            VisitPolychord(mesh, pos, chords, mark, resultCode);
          }
        } else if (resultCode == PC_NOTQUAD) {
          // if not quad, visit all the polychords passing through this coord
          pos = startPos;
          do {
            if (!pos.IsBorder()) {
              pos.FlipF();
              VisitPolychord(mesh, pos, chords, mark, resultCode);
              pos.FlipF();
            }
            pos.FlipV();
            pos.FlipE();
          } while (pos != startPos);
          VisitPolychord(mesh, startPos, chords, mark, resultCode);
        }
        VisitPolychord(mesh, startPos, chords, mark, resultCode);
      }

      // go to the next coord
      chords.Next();
      // increment the mark
      ++mark;
      if (mark == std::numeric_limits<unsigned long>::max()) {
        chords.ResetMarks();
        mark = 0;
      }
    }
  }

  /**
   * @brief SplitPolychord splits a polychord into n polychords by inserting all the needed faces.
   * @param mesh is the input polygonal mesh.
   * @param pos is a position into the polychord (not necessarily the starting border). It will be updated with changes.
   * @param n is the number of polychords to replace the input one.
   * @param facesToUpdate is a vector of face pointers to be updated after re-allocation.
   * @param verticesToUpdate is a vector of vertex pointers to be updated after re-allocation.
   */
  static void SplitPolychord (PolyMeshType &mesh, vcg::face::Pos<FaceType> &pos, const size_t n,
                              std::vector<FacePointer *> &facesToUpdate, std::vector<VertexPointer *> &verticesToUpdate) {
    vcg::tri::RequireFFAdjacency(mesh);
    vcg::tri::RequirePerFaceFlags(mesh);

    if (mesh.IsEmpty())
      return;
    if (pos.IsNull())
      return;
    if (n <= 1)
      return;

    // remember which face vertex has pos, for later updating
    int posVInd = pos.VInd();

    vcg::face::Pos<FaceType> startPos, runPos;
    PC_ResultCode result = CheckPolychordFindStartPosition(pos, startPos, false);
    if (result != PC_SUCCESS && result != PC_SINGBOTH && result != PC_SINGSIDEA && result != PC_SINGSIDEB)
      return;

    // since every face has an orientation, ensure that the new polychords are inserted on the right of the starting pos
    startPos.FlipE();
    int e = startPos.E();
    startPos.FlipE();
    if (startPos.F()->Next(startPos.E()) == e)
      startPos.FlipV();

    // clear flags
    vcg::tri::UpdateFlags<PolyMeshType>::FaceClearV(mesh);
    vcg::tri::UpdateFlags<PolyMeshType>::FaceClearS(mesh);

    // count how many faces there are
    size_t fn1 = 0, fn2 = 0;
    runPos = startPos;
    do {
      // increase the number of faces
      if (runPos.F()->IsV()) {
        ++fn2;
        --fn1;
        runPos.F()->SetS();
      } else {
        ++fn1;
        runPos.F()->SetV();
      }

      // go onto the next face
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (!runPos.IsBorder() && runPos != startPos);

    // clear flags
    runPos = startPos;
    do {
      // clear visited
      runPos.F()->ClearV();
      // go onto the next face
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (!runPos.IsBorder() && runPos != startPos);

    // compute the number of faces and vertices that must be added to the mesh in order to insert the new polychords
    size_t FN = fn1 * (n - 1) + fn2 * (n * n - 1);
    size_t VN = fn1 * (n - 1) + fn2 * (n + 1) * (n - 1);
    if (startPos.IsBorder())
      VN += n - 1;

    // add the pos to update face and vertex pointer to the list of things to update after re-allocation
    facesToUpdate.push_back(&pos.F());
    verticesToUpdate.push_back(&pos.V());
    // add the starting position's face and vertex pointers to the list of things to update after re-allocation
    facesToUpdate.push_back(&startPos.F());
    verticesToUpdate.push_back(&startPos.V());

    runPos.SetNull();
    // add faces to the mesh
    FaceIterator firstAddedFaceIt = vcg::tri::Allocator<PolyMeshType>::AddFaces(mesh, FN, facesToUpdate);
    // add vertices to the mesh
    VertexIterator firstAddedVertexIt = vcg::tri::Allocator<PolyMeshType>::AddVertices(mesh, VN, verticesToUpdate);

    // delete the added starting position's face and vertex pointers
    facesToUpdate.pop_back();
    verticesToUpdate.pop_back();
    // delete the added pos to update face and vertex pointers
    facesToUpdate.pop_back();
    verticesToUpdate.pop_back();

    // allocate and initialize 4 vertices and ffAdj for each new face
    for (FaceIterator fIt = firstAddedFaceIt; fIt != mesh.face.end(); ++fIt) {
      fIt->Alloc(4);
      for (size_t j = 0; j < 4; ++j) {
        fIt->FFp(j) = &*fIt;
        fIt->FFi(j) = j;
      }
    }

    // two structures to store temporary face data and splitting information
    struct FaceData {
      FacePointer                 faceP;
      std::vector<FacePointer>    ffpAdj;
      std::vector<int>            ffiAdj;
      std::vector<VertexPointer>  fvpAdj;
      FaceData() : faceP(0), ffpAdj(4, 0), ffiAdj(4, 0), fvpAdj(4, 0) { }
    };
    struct FaceSubdivision {
      int                                   firstEdge;
      int                                   firstVertex;
      std::vector< std::vector<FaceData> >  subfaces;
    };
#if __cplusplus >= 201103L
    std::unordered_map<FacePointer,FaceSubdivision>                     faceSubdivisions;
    typename std::unordered_map<FacePointer,FaceSubdivision>::iterator  faceSubdivisionsIt;
    typename std::unordered_map<FacePointer,FaceSubdivision>::iterator  faceSubdivisionsPrevIt;
    typename std::unordered_map<FacePointer,FaceSubdivision>::iterator  faceSubdivisionsNeighbourIt;
#else
    std::map<FacePointer,FaceSubdivision                                faceSubdivisions;
    typename std::map<FacePointer,FaceSubdivision>::iterator            faceSubdivisionsIt;
    typename std::map<FacePointer,FaceSubdivision>::iterator            faceSubdivisionsPrevIt;
    typename std::map<FacePointer,FaceSubdivision>::iterator            faceSubdivisionsNeighbourIt;
#endif
    // structure to store temporary data to assign at external faces (close to polychord)
    struct ExternalFaceData {
      FacePointer faceTo;
      FacePointer faceFrom;
      int         edgeTo;
      int         edgeFrom;
      ExternalFaceData() : faceTo(0), faceFrom(0), edgeTo(0), edgeFrom(0) { }
      ExternalFaceData(const FacePointer &ft,
                       const FacePointer &ff,
                       const int et,
                       const int ef) : faceTo(ft), faceFrom(ff), edgeTo(et), edgeFrom(ef) { }
    };
    std::list<ExternalFaceData>                                         externalFaces;
    typename std::list<ExternalFaceData>::iterator                      externalFacesIt;

    int leftEdge, rightEdge, topEdge, bottomEdge, blVInd, brVInd, tlVInd, trVInd;
    int pleftEdge, prightEdge, ptopEdge, pbottomEdge, pblVInd, pbrVInd, ptlVInd, ptrVInd;
    CoordType fromPoint, toPoint;
    FacePointer faceP;
    // first pass: make subdivisions
    runPos = startPos;
    do {
      // create temporary data
      if (!runPos.F()->IsV()) {
        runPos.F()->SetV();
        faceSubdivisionsIt = faceSubdivisions.insert(std::make_pair(runPos.F(), FaceSubdivision())).first;
        faceSubdivisionsIt->second.firstEdge = runPos.E();
        faceSubdivisionsIt->second.firstVertex = runPos.VInd();
        if (runPos.F()->IsS())
          faceSubdivisionsIt->second.subfaces.resize(n, std::vector<FaceData>(n));
        else
          faceSubdivisionsIt->second.subfaces.resize(1, std::vector<FaceData>(n));
        // assign face pointers
        faceSubdivisionsIt->second.subfaces.at(0).at(0).faceP = runPos.F();
        for (size_t j = 1; j < n; ++j, ++firstAddedFaceIt)
          faceSubdivisionsIt->second.subfaces.at(0).at(j).faceP = &*firstAddedFaceIt;
        for (size_t i = 1; i < faceSubdivisionsIt->second.subfaces.size(); ++i)
          for (size_t j = 0; j < n; ++j, ++firstAddedFaceIt)
            faceSubdivisionsIt->second.subfaces.at(i).at(j).faceP = &*firstAddedFaceIt;
        // internal face pointers adj
        rightEdge = runPos.F()->Next(runPos.E());
        leftEdge = runPos.F()->Prev(runPos.E());
        for (size_t i = 0; i < faceSubdivisionsIt->second.subfaces.size(); ++i)
          for (size_t j = 0; j < n - 1; ++j) {
            faceSubdivisionsIt->second.subfaces.at(i).at(j).ffpAdj[rightEdge] = faceSubdivisionsIt->second.subfaces.at(i).at(j+1).faceP;
            faceSubdivisionsIt->second.subfaces.at(i).at(j).ffiAdj[rightEdge] = leftEdge;
            faceSubdivisionsIt->second.subfaces.at(i).at(j+1).ffpAdj[leftEdge] = faceSubdivisionsIt->second.subfaces.at(i).at(j).faceP;
            faceSubdivisionsIt->second.subfaces.at(i).at(j+1).ffiAdj[leftEdge] = rightEdge;
          }
        topEdge = runPos.F()->Next(rightEdge);
        bottomEdge = runPos.E();
        for (size_t i = 0; i < faceSubdivisionsIt->second.subfaces.size() - 1; ++i)
          for (size_t j = 0; j < n; ++j) {
            faceSubdivisionsIt->second.subfaces.at(i).at(j).ffpAdj[topEdge] = faceSubdivisionsIt->second.subfaces.at(i+1).at(j).faceP;
            faceSubdivisionsIt->second.subfaces.at(i).at(j).ffiAdj[topEdge] = bottomEdge;
            faceSubdivisionsIt->second.subfaces.at(i+1).at(j).ffpAdj[bottomEdge] = faceSubdivisionsIt->second.subfaces.at(i).at(j).faceP;
            faceSubdivisionsIt->second.subfaces.at(i+1).at(j).ffiAdj[bottomEdge] = topEdge;
          }
        // assign old vertex pointers
        blVInd = runPos.VInd();
        brVInd = runPos.F()->Next(blVInd);
        trVInd = runPos.F()->Next(brVInd);
        tlVInd = runPos.F()->Next(trVInd);
        faceSubdivisionsIt->second.subfaces.front().front().fvpAdj.at(blVInd) = runPos.F()->V(blVInd);
        faceSubdivisionsIt->second.subfaces.front().back().fvpAdj.at(brVInd) = runPos.F()->V(brVInd);
        faceSubdivisionsIt->second.subfaces.back().back().fvpAdj.at(trVInd) = runPos.F()->V(trVInd);
        faceSubdivisionsIt->second.subfaces.back().front().fvpAdj.at(tlVInd) = runPos.F()->V(tlVInd);
        // assign new internal vertex pointers
        for (size_t i = 0; i < faceSubdivisionsIt->second.subfaces.size() - 1; ++i)
          for (size_t j = 0; j < n - 1; ++j, ++firstAddedVertexIt) {
            faceSubdivisionsIt->second.subfaces.at(i).at(j).fvpAdj.at(trVInd) = &*firstAddedVertexIt;
            faceSubdivisionsIt->second.subfaces.at(i).at(j+1).fvpAdj.at(tlVInd) = &*firstAddedVertexIt;
            faceSubdivisionsIt->second.subfaces.at(i+1).at(j).fvpAdj.at(brVInd) = &*firstAddedVertexIt;
            faceSubdivisionsIt->second.subfaces.at(i+1).at(j+1).fvpAdj.at(blVInd) = &*firstAddedVertexIt;
          }
      }

      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (!runPos.IsBorder() && runPos != startPos);

    // update subdivision iterator
    if (runPos.IsBorder())
      faceSubdivisionsIt = faceSubdivisions.end();

    // second pass: assign edge vertices and subdivisions-to-subdivisions face-face adjacency
    runPos = startPos;
    do {
      // get current and previous subdivision
      faceSubdivisionsPrevIt = faceSubdivisionsIt;
      faceSubdivisionsIt = faceSubdivisions.find(runPos.F());

      // get original indices
      bottomEdge = faceSubdivisionsIt->second.firstEdge;
      rightEdge = runPos.F()->Next(bottomEdge);
      topEdge = runPos.F()->Next(rightEdge);
      leftEdge = runPos.F()->Next(topEdge);
      blVInd = faceSubdivisionsIt->second.firstVertex;
      brVInd = runPos.F()->Next(blVInd);
      trVInd = runPos.F()->Next(brVInd);
      tlVInd = runPos.F()->Next(trVInd);
      if (faceSubdivisionsPrevIt != faceSubdivisions.end()) {
        pbottomEdge = faceSubdivisionsPrevIt->second.firstEdge;
        prightEdge = runPos.FFlip()->Next(pbottomEdge);
        ptopEdge = runPos.FFlip()->Next(prightEdge);
        pleftEdge = runPos.FFlip()->Next(ptopEdge);
        pblVInd = faceSubdivisionsPrevIt->second.firstVertex;
        pbrVInd = runPos.F()->Next(pblVInd);
        ptrVInd = runPos.F()->Next(pbrVInd);
        ptlVInd = runPos.F()->Next(ptrVInd);
      }

      // assign bottom edge vertices (and vertex adjacency with the previous subdivision) and face-to-face adjacency
      if (runPos.E() == bottomEdge) {
        // get pre-existing coords
        fromPoint = faceSubdivisionsIt->second.subfaces.front().front().fvpAdj.at(blVInd)->P();
        toPoint = faceSubdivisionsIt->second.subfaces.front().back().fvpAdj.at(brVInd)->P();
        // assign new vertices
        for (size_t j = 0; j < n - 1; ++j, ++firstAddedVertexIt) {
          firstAddedVertexIt->P() = fromPoint + (toPoint - fromPoint) * (j + 1) / n;
          faceSubdivisionsIt->second.subfaces.front().at(j).fvpAdj.at(brVInd) = &*firstAddedVertexIt;
          faceSubdivisionsIt->second.subfaces.front().at(j+1).fvpAdj.at(blVInd) = &*firstAddedVertexIt;
        }
        if (faceSubdivisionsPrevIt != faceSubdivisions.end()) {
          if (runPos.F()->FFi(bottomEdge) == ptopEdge) {
            // update face-to-vertex adjacency
            for (size_t j = 0; j < n - 1; ++j) {
              faceSubdivisionsPrevIt->second.subfaces.back().at(j).fvpAdj.at(ptrVInd) = faceSubdivisionsIt->second.subfaces.front().at(j).fvpAdj.at(brVInd);
              faceSubdivisionsPrevIt->second.subfaces.back().at(j+1).fvpAdj.at(ptlVInd) = faceSubdivisionsIt->second.subfaces.front().at(j+1).fvpAdj.at(blVInd);
            }
            // update face-to-face adjacency
            for (size_t j = 0; j < n; ++j) {
              faceSubdivisionsIt->second.subfaces.front().at(j).ffpAdj.at(bottomEdge) = faceSubdivisionsPrevIt->second.subfaces.back().at(j).faceP;
              faceSubdivisionsIt->second.subfaces.front().at(j).ffiAdj.at(bottomEdge) = ptopEdge;
              faceSubdivisionsPrevIt->second.subfaces.back().at(j).ffpAdj.at(ptopEdge) = faceSubdivisionsIt->second.subfaces.front().at(j).faceP;
              faceSubdivisionsPrevIt->second.subfaces.back().at(j).ffiAdj.at(ptopEdge) = bottomEdge;
            }
          } else if (runPos.F()->FFi(bottomEdge) == prightEdge) {
            // update face-to-vertex adjacency
            for (size_t i = 0; i < n - 1; ++i) {
              faceSubdivisionsPrevIt->second.subfaces.at(i).back().fvpAdj.at(ptrVInd) = faceSubdivisionsIt->second.subfaces.front().at(n-i-1).fvpAdj.at(blVInd);
              faceSubdivisionsPrevIt->second.subfaces.at(i+1).back().fvpAdj.at(pbrVInd) = faceSubdivisionsIt->second.subfaces.front().at(n-i-2).fvpAdj.at(brVInd);
            }
            // update face-to-face adjacency
            for (size_t j = 0; j < n; ++j) {
              faceSubdivisionsIt->second.subfaces.front().at(j).ffpAdj.at(bottomEdge) = faceSubdivisionsPrevIt->second.subfaces.at(n-j-1).back().faceP;
              faceSubdivisionsIt->second.subfaces.front().at(j).ffiAdj.at(bottomEdge) = prightEdge;
              faceSubdivisionsPrevIt->second.subfaces.at(n-j-1).back().ffpAdj.at(prightEdge) = faceSubdivisionsIt->second.subfaces.front().at(j).faceP;
              faceSubdivisionsPrevIt->second.subfaces.at(n-j-1).back().ffiAdj.at(prightEdge) = bottomEdge;
            }
          } else {
            // must be pleftEdge
            // update face-to-vertex adjacency
            for (size_t i = 0; i < n - 1; ++i) {
              faceSubdivisionsPrevIt->second.subfaces.at(i).front().fvpAdj.at(ptlVInd) = faceSubdivisionsIt->second.subfaces.front().at(i).fvpAdj.at(brVInd);
              faceSubdivisionsPrevIt->second.subfaces.at(i+1).front().fvpAdj.at(pblVInd) = faceSubdivisionsIt->second.subfaces.front().at(i+1).fvpAdj.at(blVInd);
            }
            // update face-to-face adjacency
            for (size_t j = 0; j < n; ++j) {
              faceSubdivisionsIt->second.subfaces.front().at(j).ffpAdj.at(bottomEdge) = faceSubdivisionsPrevIt->second.subfaces.at(j).front().faceP;
              faceSubdivisionsIt->second.subfaces.front().at(j).ffiAdj.at(bottomEdge) = pleftEdge;
              faceSubdivisionsPrevIt->second.subfaces.at(j).front().ffpAdj.at(pleftEdge) = faceSubdivisionsIt->second.subfaces.front().at(j).faceP;
              faceSubdivisionsPrevIt->second.subfaces.at(j).front().ffiAdj.at(pleftEdge) = bottomEdge;
            }
          }
        } else {
          // must be on border
          // update face-to-face adjacency
          for (size_t j = 0; j < n; ++j) {
            faceSubdivisionsIt->second.subfaces.front().at(j).ffpAdj.at(bottomEdge) = faceSubdivisionsIt->second.subfaces.front().at(j).faceP;
            faceSubdivisionsIt->second.subfaces.front().at(j).ffiAdj.at(bottomEdge) = bottomEdge;
          }
        }
      } else if (runPos.E() == leftEdge) {
        // get pre-existing coords
        fromPoint = faceSubdivisionsIt->second.subfaces.front().front().fvpAdj.at(blVInd)->P();
        toPoint = faceSubdivisionsIt->second.subfaces.back().front().fvpAdj.at(tlVInd)->P();
        // assign new vertices
        for (size_t i = 0; i < n - 1; ++i, ++firstAddedVertexIt) {
          firstAddedVertexIt->P() = fromPoint + (toPoint - fromPoint) * (i + 1) / n;
          faceSubdivisionsIt->second.subfaces.at(i).front().fvpAdj.at(tlVInd) = &*firstAddedVertexIt;
          faceSubdivisionsIt->second.subfaces.at(i+1).front().fvpAdj.at(blVInd) = &*firstAddedVertexIt;
        }
        // can't be on border
        if (runPos.F()->FFi(leftEdge) == ptopEdge) {
          // update face-to-vertex adjacency
          for (size_t j = 0; j < n - 1; ++j) {
            faceSubdivisionsPrevIt->second.subfaces.back().at(j).fvpAdj.at(ptrVInd) = faceSubdivisionsIt->second.subfaces.at(n-j-1).front().fvpAdj.at(blVInd);
            faceSubdivisionsPrevIt->second.subfaces.back().at(j+1).fvpAdj.at(ptlVInd) = faceSubdivisionsIt->second.subfaces.at(n-j-2).front().fvpAdj.at(tlVInd);
          }
          // update face-to-face adjacency
          for (size_t i = 0; i < n; ++i) {
            faceSubdivisionsIt->second.subfaces.at(i).front().ffpAdj.at(leftEdge) = faceSubdivisionsPrevIt->second.subfaces.back().at(n-i-1).faceP;
            faceSubdivisionsIt->second.subfaces.at(i).front().ffiAdj.at(leftEdge) = ptopEdge;
            faceSubdivisionsPrevIt->second.subfaces.back().at(n-i-1).ffpAdj.at(ptopEdge) = faceSubdivisionsIt->second.subfaces.at(i).front().faceP;
            faceSubdivisionsPrevIt->second.subfaces.back().at(n-i-1).ffiAdj.at(ptopEdge) = leftEdge;
          }
        } else if (runPos.F()->FFi(leftEdge) == prightEdge) {
          // update face-to-vertex adjacency
          for (size_t i = 0; i < n - 1; ++i) {
            faceSubdivisionsPrevIt->second.subfaces.at(i).back().fvpAdj.at(ptrVInd) = faceSubdivisionsIt->second.subfaces.at(i).front().fvpAdj.at(tlVInd);
            faceSubdivisionsPrevIt->second.subfaces.at(i+1).back().fvpAdj.at(pbrVInd) = faceSubdivisionsIt->second.subfaces.at(i+1).front().fvpAdj.at(blVInd);
          }
          // update face-to-face adjacency
          for (size_t i = 0; i < n; ++i) {
            faceSubdivisionsIt->second.subfaces.at(i).front().ffpAdj.at(leftEdge) = faceSubdivisionsPrevIt->second.subfaces.at(i).back().faceP;
            faceSubdivisionsIt->second.subfaces.at(i).front().ffiAdj.at(leftEdge) = prightEdge;
            faceSubdivisionsPrevIt->second.subfaces.at(i).back().ffpAdj.at(prightEdge) = faceSubdivisionsIt->second.subfaces.at(i).front().faceP;
            faceSubdivisionsPrevIt->second.subfaces.at(i).back().ffiAdj.at(prightEdge) = leftEdge;
          }
        } else {
          // must be runPos.F()->FFi(leftEdge) == pleftEdge
          // update face-to-vertex adjacency
          for (size_t i = 0; i < n - 1; ++i) {
            faceSubdivisionsPrevIt->second.subfaces.at(i).front().fvpAdj.at(ptlVInd) = faceSubdivisionsIt->second.subfaces.front().at(n-i-1).fvpAdj.at(blVInd);
            faceSubdivisionsPrevIt->second.subfaces.at(i+1).front().fvpAdj.at(pblVInd) = faceSubdivisionsIt->second.subfaces.front().at(n-i-2).fvpAdj.at(tlVInd);
          }
          // update face-to-face adjacency
          for (size_t i = 0; i < n; ++i) {
            faceSubdivisionsIt->second.subfaces.at(i).front().ffpAdj.at(leftEdge) = faceSubdivisionsPrevIt->second.subfaces.at(n-i-1).front().faceP;
            faceSubdivisionsIt->second.subfaces.at(i).front().ffiAdj.at(leftEdge) = pleftEdge;
            faceSubdivisionsPrevIt->second.subfaces.at(n-i-1).front().ffpAdj.at(pleftEdge) = faceSubdivisionsIt->second.subfaces.at(i).front().faceP;
            faceSubdivisionsPrevIt->second.subfaces.at(n-i-1).front().ffiAdj.at(pleftEdge) = leftEdge;
          }
        }
      } else {
        // must be runPos.E() == rightEdge
        // get pre-existing coords
        fromPoint = faceSubdivisionsIt->second.subfaces.front().back().fvpAdj.at(brVInd)->P();
        toPoint = faceSubdivisionsIt->second.subfaces.back().back().fvpAdj.at(trVInd)->P();
        // assign new vertices
        for (size_t i = 0; i < n - 1; ++i, ++firstAddedVertexIt) {
          firstAddedVertexIt->P() = fromPoint + (toPoint - fromPoint) * (i + 1) / n;
          faceSubdivisionsIt->second.subfaces.at(i).back().fvpAdj.at(trVInd) = &*firstAddedVertexIt;
          faceSubdivisionsIt->second.subfaces.at(i+1).back().fvpAdj.at(brVInd) = &*firstAddedVertexIt;
        }
        // can't be on border
        if (runPos.F()->FFi(rightEdge) == ptopEdge) {
          // update face-to-vertex adjacency
          for (size_t j = 0; j < n - 1; ++j) {
            faceSubdivisionsPrevIt->second.subfaces.back().at(j).fvpAdj.at(ptrVInd) = faceSubdivisionsIt->second.subfaces.at(j).back().fvpAdj.at(trVInd);
            faceSubdivisionsPrevIt->second.subfaces.back().at(j+1).fvpAdj.at(ptlVInd) = faceSubdivisionsIt->second.subfaces.at(j+1).back().fvpAdj.at(brVInd);
          }
          // update face-to-face adjacency
          for (size_t i = 0; i < n; ++i) {
            faceSubdivisionsIt->second.subfaces.at(i).back().ffpAdj.at(rightEdge) = faceSubdivisionsPrevIt->second.subfaces.back().at(i).faceP;
            faceSubdivisionsIt->second.subfaces.at(i).back().ffiAdj.at(rightEdge) = ptopEdge;
            faceSubdivisionsPrevIt->second.subfaces.back().at(i).ffpAdj.at(ptopEdge) = faceSubdivisionsIt->second.subfaces.at(i).back().faceP;
            faceSubdivisionsPrevIt->second.subfaces.back().at(i).ffiAdj.at(ptopEdge) = rightEdge;
          }
        } else if (runPos.F()->FFi(rightEdge) == prightEdge) {
          // update face-to-vertex adjacency
          for (size_t i = 0; i < n - 1; ++i) {
            faceSubdivisionsPrevIt->second.subfaces.at(i).back().fvpAdj.at(ptrVInd) = faceSubdivisionsIt->second.subfaces.at(n-i-1).back().fvpAdj.at(brVInd);
            faceSubdivisionsPrevIt->second.subfaces.at(i+1).back().fvpAdj.at(pbrVInd) = faceSubdivisionsIt->second.subfaces.at(n-i-2).back().fvpAdj.at(trVInd);
          }
          // update face-to-face adjacency
          for (size_t i = 0; i < n; ++i) {
            faceSubdivisionsIt->second.subfaces.at(i).back().ffpAdj.at(rightEdge) = faceSubdivisionsPrevIt->second.subfaces.at(n-i-1).back().faceP;
            faceSubdivisionsIt->second.subfaces.at(i).back().ffiAdj.at(rightEdge) = prightEdge;
            faceSubdivisionsPrevIt->second.subfaces.at(n-i-1).back().ffpAdj.at(prightEdge) = faceSubdivisionsIt->second.subfaces.at(i).back().faceP;
            faceSubdivisionsPrevIt->second.subfaces.at(n-i-1).back().ffiAdj.at(prightEdge) = rightEdge;
          }
        } else {
          // must be runPos.F()->FFi(rightEdge) == pleftEdge
          // update face-to-vertex adjacency
          for (size_t i = 0; i < n - 1; ++i) {
            faceSubdivisionsPrevIt->second.subfaces.at(i).front().fvpAdj.at(ptlVInd) = faceSubdivisionsIt->second.subfaces.at(i).back().fvpAdj.at(trVInd);
            faceSubdivisionsPrevIt->second.subfaces.at(i+1).front().fvpAdj.at(pblVInd) = faceSubdivisionsIt->second.subfaces.at(i+1).back().fvpAdj.at(brVInd);
          }
          // update face-to-face adjacency
          for (size_t i = 0; i < n; ++i) {
            faceSubdivisionsIt->second.subfaces.at(i).back().ffpAdj.at(rightEdge) = faceSubdivisionsPrevIt->second.subfaces.at(i).front().faceP;
            faceSubdivisionsIt->second.subfaces.at(i).back().ffiAdj.at(rightEdge) = pleftEdge;
            faceSubdivisionsPrevIt->second.subfaces.at(i).front().ffpAdj.at(pleftEdge) = faceSubdivisionsIt->second.subfaces.at(i).back().faceP;
            faceSubdivisionsPrevIt->second.subfaces.at(i).front().ffiAdj.at(pleftEdge) = rightEdge;
          }
        }
      }

      // update subdivision's left and right sides face-to-face adjacency
      // go on left edge
      runPos.FlipE();
      runPos.FlipV();
      if (runPos.IsBorder()) {
        if (!runPos.F()->IsS()) {
          // must be runPos.E() == leftEdge and faceSubdivisionsIt->second.subfaces.size() == 1
          faceSubdivisionsIt->second.subfaces.front().front().ffpAdj.at(leftEdge) = faceSubdivisionsIt->second.subfaces.front().front().faceP;
          faceSubdivisionsIt->second.subfaces.front().front().ffiAdj.at(leftEdge) = leftEdge;
        }
      } else if (!runPos.FFlip()->IsV()) {
        // must be runPos.E() == leftEdge and faceSubdivisionsIt->second.subfaces.size() == 1
        assert(faceSubdivisionsIt->second.subfaces.size() == 1);
        faceSubdivisionsIt->second.subfaces.front().front().ffpAdj.at(leftEdge) = runPos.FFlip();
        faceSubdivisionsIt->second.subfaces.front().front().ffiAdj.at(leftEdge) = runPos.F()->FFi(leftEdge);
        externalFaces.push_back(ExternalFaceData(runPos.FFlip(),
                                                 faceSubdivisionsIt->second.subfaces.front().front().faceP,
                                                 runPos.F()->FFi(leftEdge),
                                                 leftEdge));
      } else if (!runPos.FFlip()->IsS() && !runPos.F()->IsS()) {
        // must be runPos.E() == leftEdge and faceSubdivisionsIt->second.subfaces.size() == 1
        faceSubdivisionsNeighbourIt = faceSubdivisions.find(runPos.FFlip());
        // must be faceSubdivisionsNeighbourIt != faceSubdivisions.end() and faceSubdivisionsNeighbourIt->second.subfaces.size() == 1
        pbottomEdge = faceSubdivisionsNeighbourIt->second.firstEdge;
        prightEdge = runPos.FFlip()->Next(pbottomEdge);
        ptopEdge = runPos.FFlip()->Next(prightEdge);
        pleftEdge = runPos.FFlip()->Next(ptopEdge);
        if (runPos.F()->FFi(leftEdge) == pleftEdge) {
          faceSubdivisionsIt->second.subfaces.front().front().ffpAdj.at(leftEdge) = faceSubdivisionsNeighbourIt->second.subfaces.front().front().faceP;
          faceSubdivisionsIt->second.subfaces.front().front().ffiAdj.at(leftEdge) = pleftEdge;
          faceSubdivisionsNeighbourIt->second.subfaces.front().front().ffpAdj.at(pleftEdge) = faceSubdivisionsIt->second.subfaces.front().front().faceP;
          faceSubdivisionsNeighbourIt->second.subfaces.front().front().ffiAdj.at(pleftEdge) = leftEdge;
        } else {
          // must be runPos.F()->FFi(leftEdge) == prightEdge
          faceSubdivisionsIt->second.subfaces.front().front().ffpAdj.at(leftEdge) = faceSubdivisionsNeighbourIt->second.subfaces.front().back().faceP;
          faceSubdivisionsIt->second.subfaces.front().front().ffiAdj.at(leftEdge) = prightEdge;
          faceSubdivisionsNeighbourIt->second.subfaces.front().back().ffpAdj.at(prightEdge) = faceSubdivisionsIt->second.subfaces.front().front().faceP;
          faceSubdivisionsNeighbourIt->second.subfaces.front().back().ffiAdj.at(prightEdge) = leftEdge;
        }
      }
      // go on right edge
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      if (runPos.IsBorder()) {
        if (!runPos.F()->IsS()) {
          // must be runPos.E() == rightEdge and faceSubdivisionsIt->second.subfaces.size() == 1
          faceSubdivisionsIt->second.subfaces.front().back().ffpAdj.at(rightEdge) = faceSubdivisionsIt->second.subfaces.front().back().faceP;
          faceSubdivisionsIt->second.subfaces.front().back().ffiAdj.at(rightEdge) = rightEdge;
        }
      } else if (!runPos.FFlip()->IsV()) {
        // must be runPos.E() == rightEdge and faceSubdivisionsIt->second.subfaces.size() == 1
        faceSubdivisionsIt->second.subfaces.front().back().ffpAdj.at(rightEdge) = runPos.FFlip();
        faceSubdivisionsIt->second.subfaces.front().back().ffiAdj.at(rightEdge) = runPos.F()->FFi(rightEdge);
        externalFaces.push_back(ExternalFaceData(runPos.FFlip(),
                                                 faceSubdivisionsIt->second.subfaces.front().back().faceP,
                                                 runPos.F()->FFi(rightEdge),
                                                 rightEdge));
      } else if (!runPos.FFlip()->IsS() && !runPos.F()->IsS()) {
        // must be runPos.E() == rightEdge and faceSubdivisionsIt->second.subfaces.size() == 1
        faceSubdivisionsNeighbourIt = faceSubdivisions.find(runPos.FFlip());
        // must be faceSubdivisionsNeighbourIt != faceSubdivisions.end() and faceSubdivisionsNeighbourIt->second.subfaces.size() == 1
        pbottomEdge = faceSubdivisionsNeighbourIt->second.firstEdge;
        prightEdge = runPos.FFlip()->Next(pbottomEdge);
        ptopEdge = runPos.FFlip()->Next(prightEdge);
        pleftEdge = runPos.FFlip()->Next(ptopEdge);
        if (runPos.F()->FFi(rightEdge) == pleftEdge) {
          faceSubdivisionsIt->second.subfaces.front().back().ffpAdj.at(rightEdge) = faceSubdivisionsNeighbourIt->second.subfaces.front().front().faceP;
          faceSubdivisionsIt->second.subfaces.front().back().ffiAdj.at(rightEdge) = pleftEdge;
          faceSubdivisionsNeighbourIt->second.subfaces.front().front().ffpAdj.at(pleftEdge) = faceSubdivisionsIt->second.subfaces.front().back().faceP;
          faceSubdivisionsNeighbourIt->second.subfaces.front().front().ffiAdj.at(pleftEdge) = rightEdge;
        } else {
          // must be runPos.F()->FFi(rightEdge) == prightEdge
          faceSubdivisionsIt->second.subfaces.front().back().ffpAdj.at(rightEdge) = faceSubdivisionsNeighbourIt->second.subfaces.front().back().faceP;
          faceSubdivisionsIt->second.subfaces.front().back().ffiAdj.at(rightEdge) = prightEdge;
          faceSubdivisionsNeighbourIt->second.subfaces.front().back().ffpAdj.at(prightEdge) = faceSubdivisionsIt->second.subfaces.front().back().faceP;
          faceSubdivisionsNeighbourIt->second.subfaces.front().back().ffiAdj.at(prightEdge) = rightEdge;
        }
      }

      // go on top edge
      runPos.FlipE();
      runPos.FlipV();
      if (runPos.IsBorder()) {
        // assign top edge vertices and face-to-border adjacency
        if (runPos.E() == topEdge) {
          // get pre-existing coords
          fromPoint = faceSubdivisionsIt->second.subfaces.back().front().fvpAdj.at(tlVInd)->P();
          toPoint = faceSubdivisionsIt->second.subfaces.back().back().fvpAdj.at(trVInd)->P();
          // assign new vertices
          for (size_t j = 0; j < n - 1; ++j, ++firstAddedVertexIt) {
            firstAddedVertexIt->P() = fromPoint + (toPoint - fromPoint) * (j + 1) / n;
            faceSubdivisionsIt->second.subfaces.back().at(j).fvpAdj.at(trVInd) = &*firstAddedVertexIt;
            faceSubdivisionsIt->second.subfaces.back().at(j+1).fvpAdj.at(tlVInd) = &*firstAddedVertexIt;
          }
          // update face-to-face adjacency
          for (size_t j = 0; j < n; ++j) {
            faceSubdivisionsIt->second.subfaces.back().at(j).ffpAdj.at(topEdge) = faceSubdivisionsIt->second.subfaces.back().at(j).faceP;
            faceSubdivisionsIt->second.subfaces.back().at(j).ffiAdj.at(topEdge) = topEdge;
          }
        } else if (runPos.E() == leftEdge) {
          // get pre-existing coords
          fromPoint = faceSubdivisionsIt->second.subfaces.front().front().fvpAdj.at(blVInd)->P();
          toPoint = faceSubdivisionsIt->second.subfaces.back().front().fvpAdj.at(tlVInd)->P();
          // assign new vertices
          for (size_t i = 0; i < n - 1; ++i, ++firstAddedVertexIt) {
            firstAddedVertexIt->P() = fromPoint + (toPoint - fromPoint) * (i + 1) / n;
            faceSubdivisionsIt->second.subfaces.at(i).front().fvpAdj.at(tlVInd) = &*firstAddedVertexIt;
            faceSubdivisionsIt->second.subfaces.at(i+1).front().fvpAdj.at(blVInd) = &*firstAddedVertexIt;
          }
          // update face-to-face adjacency
          for (size_t i = 0; i < n; ++i) {
            faceSubdivisionsIt->second.subfaces.at(i).front().ffpAdj.at(leftEdge) = faceSubdivisionsIt->second.subfaces.at(i).front().faceP;
            faceSubdivisionsIt->second.subfaces.at(i).front().ffiAdj.at(leftEdge) = leftEdge;
          }
        } else {
          // must be runPos.E() == rightEdge
          // get pre-existing coords
          fromPoint = faceSubdivisionsIt->second.subfaces.front().back().fvpAdj.at(brVInd)->P();
          toPoint = faceSubdivisionsIt->second.subfaces.back().back().fvpAdj.at(trVInd)->P();
          // assign new vertices
          for (size_t i = 0; i < n - 1; ++i, ++firstAddedVertexIt) {
            firstAddedVertexIt->P() = fromPoint + (toPoint - fromPoint) * (i + 1) / n;
            faceSubdivisionsIt->second.subfaces.at(i).back().fvpAdj.at(trVInd) = &*firstAddedVertexIt;
            faceSubdivisionsIt->second.subfaces.at(i+1).back().fvpAdj.at(brVInd) = &*firstAddedVertexIt;
          }
          // update face-to-face adjacency
          for (size_t i = 0; i < n; ++i) {
            faceSubdivisionsIt->second.subfaces.at(i).back().ffpAdj.at(rightEdge) = faceSubdivisionsIt->second.subfaces.at(i).back().faceP;
            faceSubdivisionsIt->second.subfaces.at(i).back().ffiAdj.at(rightEdge) = rightEdge;
          }
        }
      } else {
        // go onto the next face
        runPos.FlipF();
      }
    } while (!runPos.IsBorder() && runPos != startPos);

    // final pass: compute coords of new internal vertices, copy all data into mesh and clear flags
    for (faceSubdivisionsIt = faceSubdivisions.begin(); faceSubdivisionsIt != faceSubdivisions.end(); ++faceSubdivisionsIt) {
      // clear flags
      faceSubdivisionsIt->first->ClearV();
      if (faceSubdivisionsIt->first->IsS())
        faceSubdivisionsIt->first->ClearS();

      // get vertex indices
      blVInd = faceSubdivisionsIt->second.firstVertex;
      brVInd = faceSubdivisionsIt->first->Next(blVInd);

      // compute coords on bottom side and internal, horizontally
      for (size_t i = 1; i < faceSubdivisionsIt->second.subfaces.size(); ++i) {
        fromPoint = faceSubdivisionsIt->second.subfaces.at(i).front().fvpAdj.at(blVInd)->P();
        toPoint = faceSubdivisionsIt->second.subfaces.at(i).back().fvpAdj.at(brVInd)->P();
        for (size_t j = 1; j < n; ++j)
          faceSubdivisionsIt->second.subfaces.at(i).at(j).fvpAdj.at(blVInd)->P() = fromPoint + (toPoint - fromPoint) * j / n;
      }

      // finally, copy data into mesh
      for (size_t i = 0; i < faceSubdivisionsIt->second.subfaces.size(); ++i)
        for (size_t j = 0; j < n; ++j) {
          faceP = faceSubdivisionsIt->second.subfaces.at(i).at(j).faceP;
          for (size_t k = 0; k < faceSubdivisionsIt->second.subfaces.at(i).at(j).ffpAdj.size(); ++k)
            faceP->FFp(k) = faceSubdivisionsIt->second.subfaces.at(i).at(j).ffpAdj.at(k);
          for (size_t k = 0; k < faceSubdivisionsIt->second.subfaces.at(i).at(j).ffiAdj.size(); ++k)
            faceP->FFi(k) = faceSubdivisionsIt->second.subfaces.at(i).at(j).ffiAdj.at(k);
          for (size_t k = 0; k < faceSubdivisionsIt->second.subfaces.at(i).at(j).fvpAdj.size(); ++k)
            faceP->V(k) = faceSubdivisionsIt->second.subfaces.at(i).at(j).fvpAdj.at(k);
        }
    }
    // very last step: update external faces adjacency
    for (externalFacesIt = externalFaces.begin(); externalFacesIt != externalFaces.end(); ++externalFacesIt) {
      externalFacesIt->faceTo->FFp(externalFacesIt->edgeTo) = externalFacesIt->faceFrom;
      externalFacesIt->faceTo->FFi(externalFacesIt->edgeTo) = externalFacesIt->edgeFrom;
    }
    // very very last step: update pos
    pos.V() = pos.F()->V(posVInd);
  }

  /**
   * @brief SplitPolychord splits a polychord into n polychords by inserting all the needed faces.
   * @param mesh is the input polygonal mesh.
   * @param pos is a position into the polychord (not necessarily the starting border). It will be updated with changes.
   * @param n is the number of polychords to replace the input one.
   */
  static void SplitPolychord (PolyMeshType &mesh, vcg::face::Pos<FaceType> &pos, const size_t n) {
    std::vector<FacePointer *> facesToUpdate;
    std::vector<VertexPointer *> verticesToUpdate;
    SplitPolychord(mesh, pos, n, facesToUpdate, verticesToUpdate);
  }

  /**
   * @brief CheckPolychordFindStartPosition checks if it's a collapsable polychord.
   * @param pos Input The starting position.
   * @param startPos Output the new starting position (in case of borders).
   * @param checkSing true if singularities on both sides are not allowed.
   * @return PC_SUCCESS if it's a collapsable polychord, otherwise the code for the cause (startPos is on it).
   */
  static PC_ResultCode CheckPolychordFindStartPosition (const vcg::face::Pos<FaceType> &pos,
                                                        vcg::face::Pos<FaceType> &startPos,
                                                        const bool checkSing = true) {
    assert(!pos.IsNull());
    int valence = 0;
    bool singSideA = false, singSideB = false;
    bool borderSideA = false, borderSideB = false;
    bool polyBorderFound = false;
    vcg::face::JumpingPos<FaceType> jmpPos;

    startPos = pos;
    do {
      // check if it is a quad
      if (startPos.F()->VN() != 4)
        return PC_NOTQUAD;
      // check manifoldness
      if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
        return PC_NOTMANIF;
      startPos.FlipV();
      if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
        return PC_NOTMANIF;
      startPos.FlipV();

      // check if side A is on border
      startPos.FlipE();
      if (startPos.IsBorder())
        borderSideA = true;
      startPos.FlipE();
      // check if side B is on border
      startPos.FlipV();
      startPos.FlipE();
      if (startPos.IsBorder())
        borderSideB = true;
      startPos.FlipE();
      startPos.FlipV();

      // check if singularities are not in both sides
      if (checkSing) {
        // compute the valence of the vertex on side B
        startPos.FlipV();
        valence = startPos.NumberOfIncidentVertices();
        // if the vertex is on border increment its valence by 1 (virtually connect it to a dummy vertex)
        jmpPos.Set(startPos.F(), startPos.E(), startPos.V());
        if (jmpPos.FindBorder())
          ++valence;
        if (valence != 4)
          singSideB = true;
        // a 2-valence internl vertex cause a polychord to touch itself, producing non-2manifoldness
        // in that case, a 2-valence vertex is dealt as 2 singularities in both sides
        if (valence == 2 && !borderSideB)
          singSideA = true;
        // compute the valence of the vertex on side A
        startPos.FlipV();
        valence = startPos.NumberOfIncidentVertices();
        // if the vertex is on border increment its valence by 1 (virtually connect it to a dummy vertex)
        jmpPos.Set(startPos.F(), startPos.E(), startPos.V());
        if (jmpPos.FindBorder())
          ++valence;
        if (valence != 4)
          singSideA = true;
        // a 2-valence internal vertex cause a polychord to touch itself, producing non-2manifoldness
        // in that case, a 2-valence vertex is dealt as 2 singularities in both sides
        if (valence == 2 && !borderSideA)
          singSideB = true;
      }

      // if the first border has been reached, go on the other direction to find the other border
      if (startPos != pos && startPos.IsBorder() && !polyBorderFound) {
        startPos = pos;
        startPos.FlipF();
        polyBorderFound = true;
      }

      // if the other border has been reached, return
      if (polyBorderFound && startPos.IsBorder())
        break;

      // go to the next edge
      startPos.FlipE();
      startPos.FlipV();
      startPos.FlipE();
      // check manifoldness
      if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
        return PC_NOTMANIF;
      startPos.FlipV();
      if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
        return PC_NOTMANIF;
      startPos.FlipV();
      // go to the next face
      startPos.FlipF();
    } while (startPos != pos);

    // polychord with singularities on both sides can not collapse
    if ((singSideA && singSideB) ||
        (singSideA && borderSideB) ||
        (singSideB && borderSideA))
      return PC_SINGBOTH;

    // polychords that are rings and have borders on both sides can not collapse
    if (!polyBorderFound && borderSideA && borderSideB)
      return PC_SINGBOTH;

    // if there are singularities or borders on the side A, remember to keep coordinates on it
    if (singSideA || borderSideA)
      return PC_SINGSIDEA;
    // if there are singularities or borders on the side B, remember to keep coordinates on it
    if (singSideB || borderSideB)
      return PC_SINGSIDEB;

    return PC_SUCCESS;
  }

  /**
   * @brief VisitPolychord updates the information of a polychord.
   * @param mesh The mesh used for getting the face index.
   * @param startPos The starting position.
   * @param chords The vector of chords.
   * @param mark The mark.
   * @param q The visiting type.
   */
  static void VisitPolychord (const PolyMeshType &mesh,
                              const vcg::face::Pos<FaceType> &startPos,
                              PC_Chords &chords,
                              const unsigned long mark,
                              const PC_ResultCode q) {
    assert(!startPos.IsNull());
    vcg::face::Pos<FaceType> tmpPos, runPos = startPos;
    std::pair<size_t, unsigned char> face_edge(std::numeric_limits<size_t>::max(), 0);

    // follow the sequence of quads
    do {
      // check manifoldness
      tmpPos = runPos;
      do {
        if (runPos.F()->VN() != 4)  // non-quads are not visited
          return;
        if (!tmpPos.IsManifold()) {
          // update current coord
          face_edge.first = vcg::tri::Index(mesh, tmpPos.F());
          face_edge.second = tmpPos.E()%2;
          chords.UpdateCoord(chords[face_edge], mark, q);
          face_edge.second = (tmpPos.E()+1)%2;
          chords.UpdateCoord(chords[face_edge], mark, q);
          return;
        }
        tmpPos.FlipV();
        tmpPos.FlipE();
      } while (tmpPos != runPos);

      // update current coord
      face_edge.first = vcg::tri::Index(mesh, runPos.F());
      face_edge.second = runPos.E()%2;
      chords.UpdateCoord(chords[face_edge], mark, q);
      // if the polychord has to collapse, i.e. q == PC_SUCCESS, also visit the orthogonal coord
      if (q == PC_SUCCESS) {
        face_edge.second = (runPos.E()+1)%2;
        chords.UpdateCoord(chords[face_edge], mark, q);
      }

      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (runPos != startPos && !runPos.IsBorder() && runPos.F()->VN() == 4);
    assert(runPos == startPos || vcg::face::IsBorder(*startPos.F(),startPos.E())
           || runPos.F()->VN() != 4 || startPos.FFlip()->VN() != 4);
  }

  /**
   * @brief IsVertexAdjacentToAnyNonManifoldEdge checks if a vertex is adjacent to any non-manifold edge.
   * @param pos The starting position.
   * @return true if adjacent to non-manifold edges, false otherwise.
   */
  static bool IsVertexAdjacentToAnyNonManifoldEdge (const vcg::face::Pos<FaceType> &pos) {
    assert(!pos.IsNull());
    vcg::face::JumpingPos<FaceType> jmpPos;
    jmpPos.Set(pos.F(), pos.E(), pos.V());
    do {
      assert(!jmpPos.FFlip()->IsD());
      if (!jmpPos.IsManifold())
        return true;
      jmpPos.NextFE();
    } while (jmpPos != pos);
    return false;
  }

  /**
   * @brief IsPolychordSelfIntersecting checks if the input polychord intersects itself.
   * @warning Don't call this function without being sure that it's a polychord
   * (i.e. call CheckPolychordFindStartPoint() before calling IsPolychordSelfIntersecting().
   * @param mesh The mesh used for getting the face index.
   * @param startPos The starting position.
   * @param chords The vector of chords.
   * @param mark The current mark, used to identify quads already visited.
   * @return true if it intersects itself, false otherwise.
   */
  static bool IsPolychordSelfIntersecting (const PolyMeshType &mesh,
                                           const vcg::face::Pos<FaceType> &startPos,
                                           const PC_Chords &chords,
                                           const unsigned long mark) {
    assert(!startPos.IsNull());
    vcg::face::Pos<FaceType> runPos = startPos;
    vcg::face::Pos<FaceType> tmpPos;
    std::pair<size_t, unsigned char> face_edge(std::numeric_limits<size_t>::max(), 0);
    do {
      assert(runPos.F()->VN() == 4);
      // check if we've already crossed this face
      face_edge.first = vcg::tri::Index(mesh, runPos.F());
      face_edge.second = (runPos.E()+1)%2;
      if (chords[face_edge].mark == mark)
        return true;
      // if this coord is adjacent to another coord of the same polychord
      // i.e., this polychord touches itself without intersecting
      // it might cause a wrong collapse, producing holes and non-2manifoldness
      tmpPos = runPos;
      tmpPos.FlipE();
      if (!tmpPos.IsBorder()) {
        tmpPos.FlipF();
        face_edge.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge.second = (tmpPos.E()+1)%2;
        if (chords[face_edge].mark == mark)
          return true;
        // this should never hapen:
        face_edge.second = tmpPos.E()%2;
        if (chords[face_edge].mark == mark)
          return true;
      }
      tmpPos = runPos;
      tmpPos.FlipV();
      tmpPos.FlipE();
      if (!tmpPos.IsBorder()) {
        tmpPos.FlipF();
        face_edge.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge.second = (tmpPos.E()+1)%2;
        if (chords[face_edge].mark == mark)
          return true;
        // this should never hapen:
        face_edge.second = tmpPos.E()%2;
        if (chords[face_edge].mark == mark)
          return true;
      }
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (runPos != startPos && !runPos.IsBorder());

    return false;
  }

  /**
   * @brief WillPolychordBeManifold checks whether a polychord starting at startPos would cause non-manifoldness
   * if it was collapsed.
   * @note VisitPolychord() should be called before this method.
   * @param mesh The input mesh.
   * @param startPos The starting Pos.
   * @param chords The vector of chords.
   * @param mark The current mark, used to identify quads already visited.
   * @return true if manifoldness remains, false otherwise.
   */
  static bool WillPolychordBeManifold(const PolyMeshType &mesh,
                                      const vcg::face::Pos<FaceType> &startPos,
                                      PC_Chords &chords,
                                      const unsigned long mark) {
    assert(!startPos.IsNull());
    vcg::face::Pos<FaceType> runPos = startPos;
    vcg::face::JumpingPos<FaceType> tmpPos;
    std::pair<size_t, unsigned char> face_edge1(std::numeric_limits<size_t>::max(), 0);
    std::pair<size_t, unsigned char> face_edge2(std::numeric_limits<size_t>::max(), 0);
    bool in = true;
    unsigned int nTraversal = 0;

    // second step: check
    runPos = startPos;
    do {
      face_edge1.first = vcg::tri::Index(mesh, runPos.F());
      face_edge1.second = runPos.E() % 2;
      assert(chords[face_edge1].mark == mark);
      // check one vertex
      runPos.FlipV();
      in = true;
      nTraversal = 0;
      tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
      do {
        if (tmpPos.IsBorder() && in) {
          in = false;
          nTraversal++;
        }
        // go to next edge
        tmpPos.NextFE();
        // check if this face is already visited
        face_edge1.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge1.second = tmpPos.E() % 2;
        face_edge2.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge2.second = (tmpPos.E() + 1) % 2;
        if (in && chords[face_edge1].mark != mark && chords[face_edge2].mark != mark) {
          in = false;
          ++nTraversal;
        } else if (!in && (chords[face_edge1].mark == mark || chords[face_edge2].mark == mark)) {
          in = true;
          ++nTraversal;
        }
      } while (tmpPos != runPos);
      assert(in);
      assert(nTraversal % 2 == 0);
      if (nTraversal > 2)
        return false;

      // check other vertex
      runPos.FlipV();
      in = true;
      nTraversal = 0;
      tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
      do {
        if (tmpPos.IsBorder() && in) {
          in = false;
          ++nTraversal;
        }
        // go to next edge
        tmpPos.NextFE();
        // check if this face is already visited
        face_edge1.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge1.second = tmpPos.E() % 2;
        face_edge2.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge2.second = (tmpPos.E() + 1) % 2;
        if (in && chords[face_edge1].mark != mark && chords[face_edge2].mark != mark) {
          in = false;
          nTraversal++;
        } else if (!in && (chords[face_edge1].mark == mark || chords[face_edge2].mark == mark)) {
          in = true;
          ++nTraversal;
        }
      } while (tmpPos != runPos);
      assert(in);
      assert(nTraversal % 2 == 0);
      if (nTraversal > 2)
        return false;

      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (runPos != startPos && !runPos.IsBorder());
    if (runPos.IsBorder()) {
      // check one vertex
      runPos.FlipV();
      in = true;
      nTraversal = 0;
      tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
      do {
        if (tmpPos.IsBorder() && in) {
          in = false;
          ++nTraversal;
        }
        // go to next edge
        tmpPos.NextFE();
        // check if this face is already visited
        face_edge1.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge1.second = tmpPos.E() % 2;
        face_edge2.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge2.second = (tmpPos.E() + 1) % 2;
        if (in && chords[face_edge1].mark != mark && chords[face_edge2].mark != mark) {
          in = false;
          ++nTraversal;
        } else if (!in && (chords[face_edge1].mark == mark || chords[face_edge2].mark == mark)) {
          in = true;
          ++nTraversal;
        }
      } while (tmpPos != runPos);
      assert(in);
      assert(nTraversal % 2 == 0);
      if (nTraversal > 2)
        return false;

      // check other vertex
      runPos.FlipV();
      in = true;
      nTraversal = 0;
      tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
      do {
        if (tmpPos.IsBorder() && in) {
          in = false;
          ++nTraversal;
        }
        // go to next edge
        tmpPos.NextFE();
        // check if this face is already visited
        face_edge1.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge1.second = tmpPos.E() % 2;
        face_edge2.first = vcg::tri::Index(mesh, tmpPos.F());
        face_edge2.second = (tmpPos.E() + 1) % 2;
        if (in && chords[face_edge1].mark != mark && chords[face_edge2].mark != mark) {
          in = false;
          ++nTraversal;
        } else if (!in && (chords[face_edge1].mark == mark || chords[face_edge2].mark == mark)) {
          in = true;
          ++nTraversal;
        }
      } while (tmpPos != runPos);
      assert(in);
      assert(nTraversal % 2 == 0);
      if (nTraversal > 2)
        return false;
    }

    return true;
  }
};

}
}

#endif // POLYGON_Polychord_COLLAPSE_H
