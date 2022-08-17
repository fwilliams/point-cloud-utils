#include <npe.h>
#include <igl/adjacency_matrix.h>
#include <numeric>
#include <queue>

#include "common/common.h"


namespace {
  template < typename Atype, typename DerivedC, typename CountType>
  inline int connected_components(const Eigen::SparseMatrix<Atype> & adjacency_matrix,
                                  Eigen::PlainObjectBase<DerivedC> & out_vertex_components,
                                  std::vector<CountType> & out_component_counts) {
    using Index = typename Eigen::SparseMatrix<Atype>::Index;
    using SparseMatrixRowIter = typename Eigen::SparseMatrix<Atype>::InnerIterator;
    assert(adjacency_matrix.cols() == adjacency_matrix.rows() && "A should be square");
    
    const auto num_vertices = adjacency_matrix.rows();

    // A value of num_vertices  means not yet visited
    out_vertex_components.setConstant(num_vertices, 1, num_vertices);
    
    for(Eigen::Index f = 0; f < num_vertices; f += 1) {
      const CountType current_component = out_component_counts.size();
      CountType num_vertices_in_component = 0;

      // Enable Ctrl-C to cancel
      if (PyErr_CheckSignals() != 0) {
        throw pybind11::error_already_set();
      }

      // Already processed this vertex
      if(out_vertex_components(f) < num_vertices) {
        continue;
      }

      // Unprocessed vertex, run BFS to collect all vertices connected to it
      std::queue<Index> Q;
      Q.push(f);
      while(!Q.empty()) {
        // Enable Ctrl-C to cancel
        if (PyErr_CheckSignals() != 0) {
          throw pybind11::error_already_set();
        }

        const Index g = Q.front();
        Q.pop();

        if(out_vertex_components(g) < num_vertices) {  // already processed this vertex
          continue;
        }

        out_vertex_components(g) = current_component;
        num_vertices_in_component += 1;
        for(SparseMatrixRowIter it(adjacency_matrix, g); it; ++it) {
          const Index n = it.index();
          if(out_vertex_components(n) < num_vertices) {  // already processed this vertex
            continue;
          }
          Q.push(n);
        }
      }
      out_component_counts.push_back(num_vertices_in_component);
    }
    return out_component_counts.size();
  }
}

const char* connected_components_doc = R"igl_Qu8mg5v7(
Determine the connected components of a mesh

Parameters
----------
v: (#v, 3)-shaped array of mesh vertex positions (one vertex position per row)
f: (#f, 3)-shaped array of mesh face indexes into v (a row (fi, fj, fk) indicate the 3 vertices of a face)

Returns
-------
A tuple (cv, nv, cf, nf) where:
  - cv is a (#vertices,)-shaped array of integer indexes (starting from 0) indicating which
    component each vertex belongs to. i.e. cv[i] is the component of the vertex v[i].
  - nv is the number of vertices in each connected component. i.e. nv[j] is the number of vertices
    in component j
  - cf is a (#faces,)-shaped array of integer indexes (starting from 0) indicating which
    component each face belongs to. i.e. cf[i] is the component of the face f[i].
  - nf is the number of faces in each connected component. i.e. nf[j] is the number of faces in
    component j
)igl_Qu8mg5v7";
npe_function(connected_components)
npe_doc(connected_components_doc)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_long, dense_longlong)
npe_begin_code()
{
    Eigen::SparseMatrix<npe_Scalar_f> A;
    igl::adjacency_matrix(f, A);

    npe_Matrix_f c_v, c_f, v_counts, f_counts;
    std::vector<npe_Scalar_f> v_counts_vector;
    connected_components(A, c_v, v_counts_vector);
    v_counts.resize(v_counts_vector.size(), 1);
    for (int i = 0; i < v_counts_vector.size(); i += 1) {
      v_counts(i, 0) = v_counts_vector[i];
    }

    c_f.resize(f.rows(), 1);
    f_counts = npe_Matrix_f::Zero(v_counts.rows(), 1);
    for (int i = 0; i < c_f.rows(); i++) {
        npe_Scalar_f face_component = c_v(f(i, 0), 0);
        c_f(i, 0) = face_component;
        f_counts(face_component, 0) += 1;
    }

    return std::make_tuple(npe::move(c_v), npe::move(v_counts), 
                           npe::move(c_f), npe::move(f_counts));
}
npe_end_code()
