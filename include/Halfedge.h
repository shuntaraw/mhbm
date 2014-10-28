// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include "mesh.h"

namespace hbm {

    /// halfedge
    class Halfedge {
    public:
        Halfedge *pair() const;
        Halfedge *next() const;
        int vertex() const;
        int face() const;
        Halfedge *previous() const;
        void Invalidate();

    private:
        Halfedge *pair_ = 0; ///< pointer to the paired edge. 'nullptr' on border edges
        Halfedge *next_ = 0; ///< pointer to the next edge. nullprt if invalidated.
        int vertex_ = -1; ///< zero-based index of the destination vertex. -1 if invalidated.
        int face_ = -1; ///< zero-based index of the face connecting to this edge. -1 if invalidated.

        friend class HalfedgeMesh;
        friend class EditableHalfedgeMesh;
    };

    /// halfedge data structure
    class HalfedgeMesh {
    public:
        HalfedgeMesh() = default;
        HalfedgeMesh(const HalfedgeMesh&) = delete;
        HalfedgeMesh& operator=(const HalfedgeMesh&) = delete;

        int num_edges() const;
        int num_vertices() const;
        int num_faces() const;

        /// returns one of halfedges connecting to a vertex
        /// @return halfedge; the pair is null if the vertex is on the boundary
        Halfedge *vertex_edge(int vid) const;

        /// returns one of halfedges owned by a face
        Halfedge *face_edge(int fid) const;

        /// generate halfedge data structure
        void Construct(const CMesh& mesh);

    private:
        /// check if the halfedge structure is consistent (for debug)
        bool IsConsistent(const CMesh& mesh) const;

        /// generate a halfedge data structure
        void ConstructHalfedges(const CMesh& mesh);

        /// precompute the mapping from a face to a halfedge.
        void ConstructFaceMap(int nfaces);

        /// precompute the mapping from a vertex to a halfedge.
        void ConstructVertexMap(int nvertices);

    private:
        std::vector<Halfedge> halfedges_; ///< halfedges
        std::vector<Halfedge *> vertex_edges_; ///< mapping from vertex to halfedge
        std::vector<Halfedge *> face_edges_; ///< mapping from face to halfedge
    };


    std::vector<int> GetAdjacentVertices(const HalfedgeMesh& hds, int vid);

    std::vector<int> GetAdjacentFaces(const HalfedgeMesh& hds, int vid);

    int GetDegree(const HalfedgeMesh& hds, int vid);

    bool IsBorder(const HalfedgeMesh& hds, int vid);


} // namespace slib
