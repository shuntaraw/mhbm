// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <cassert>
#include <vector>
#include "exception_util.h"

namespace slib {
namespace mesh {

/// halfedge
class Halfedge {
public:
    Halfedge *pair()const {
        return pair_;
    }

    Halfedge *next()const {
        return next_;
    }

    int vertex()const {
        return vertex_;
    }

    int face()const {
        return face_;
    }

    Halfedge *previous()const {
        auto prev = next();
        while (prev->next() != this) {
            prev = prev->next();
        }
        return prev;
    }

    void Invalidate() {
        pair_ = 0;
        next_ = 0;
        vertex_ = -1;
        face_ = -1;
    }

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

    int num_edges() const {
        return halfedges_.size();
    }

    int num_vertices() const {
        return vertex_edges_.size();
    }

    int num_faces() const {
        return face_edges_.size();
    }

    /// returns one of halfedges connecting to a vertex
    /// @return halfedge; the pair is null if the vertex is on the boundary
    Halfedge *vertex_edge(int vid) const {
        return vertex_edges_[vid];
    }

    /// returns one of halfedges owned by a face
    Halfedge *face_edge(int fid) const {
        return face_edges_[fid];
    }

    /// generate halfedge data structure
    /// @tparam Mesh polygon mesh
    template <typename Mesh>
    void Construct(const Mesh& mesh) {
        ConstructHalfedges(mesh);
        assert(IsConsistent(mesh));
        ConstructVertexMap(mesh.vertices.size());
        ConstructFaceMap(mesh.faces.size());
    }

private:
    /// check if the halfedge structure is consistent (for debug)
    /// @tparam Mesh polygon mesh
    template <typename Mesh>
    bool IsConsistent(const Mesh& mesh) const {
        int num_unused = 0;
        for (auto& edge : halfedges_) {
            // pointers
            if (edge.face() == -1 || !edge.next() || edge.vertex() == -1) {
                if (edge.face() == -1 && !edge.next() && edge.vertex() == -1) {
                    num_unused++;
                } else {
                    ThrowRuntimeError("uninitialized edge");
                }
            }

            // face
            int fid = edge.face();
            const Halfedge *current = &edge;
            for (int e = 0; e < mesh.faces[fid].index.size(); e++) {
                if (current->face() != fid) {
                    ThrowRuntimeError("face inconsistent");
                }
                current = current->next();
            }
            if (current != &edge) {
                ThrowLogicError("face not closed");
            }

            // pair
            if (edge.pair()) {
                if (edge.pair()->pair() != &edge) {
                    ThrowRuntimeError("inconsistent halfedge pair");
                }
            }
        }
        if (num_unused) {
            std::clog << num_unused << " unused halfedges" << std::endl;
        }
        return true;
    }

    /// generate a halfedge data structure
    /// @tparam Mesh polygon mesh
    template <typename Mesh>
    void ConstructHalfedges(
        const Mesh& mesh ///< [in] polygon mesh
    ) {
        // reserve memory
        int nhalfedges = 0;
        for (auto& face : mesh.faces) {
            nhalfedges += face.index.size(); // number of halfedges
        }
        halfedges_.clear();
        halfedges_.reserve(nhalfedges);

        std::vector<std::tuple<int, int, Halfedge *>> edges;
        edges.reserve(nhalfedges);

        // insert halfedges without pairs
        for (int fid = 0; fid < mesh.faces.size(); fid++) {
            auto& face = mesh.faces[fid];
            Halfedge *first, *prev = 0;
            for (int vid : face.index) {
                halfedges_.push_back(Halfedge());
                Halfedge *current = &halfedges_.back();

                // face
                current->face_ = fid;
                // vertex
                current->vertex_ = vid;
                // next edge
                if (prev) {
                    prev->next_ = current;
                    edges.push_back(std::make_tuple(prev->vertex(), current->vertex(), current));
                } else {
                    first = current;
                }
                prev = current;
            }

            // the next of the last is the first
            prev->next_ = first;
            edges.push_back(std::make_tuple(prev->vertex(), first->vertex(), first));
        }

        std::sort(edges.begin(), edges.end());

        // pair
        for (auto it = edges.begin(), end = edges.end(); it != end; ++it) {
            auto current = std::get<2>(*it);
            if (!current->pair()) {
                auto prev = current;
                while (prev->next() != current) {
                    prev = prev->next_;
                }
                auto upper = std::upper_bound(
                                 it + 1,
                                 end,
                                 std::make_tuple(current->vertex(), prev->vertex(), (Halfedge *)0));
                if (upper == edges.end() || std::get<0>(*upper) != current->vertex() || std::get<1>(*upper) != prev->vertex()) {
                    // border edge
                } else {
                    auto pair = std::get<2>(*upper);
                    current->pair_ = pair;
                    if (pair->pair()) {
                        ThrowRuntimeError("duplicate edge (%d-%d)", pair->vertex(), current->vertex());
                    }
                    pair->pair_ = current;
                    auto upper2 = upper + 1;
                    if (upper2 != end && std::get<0>(*upper2) == current->vertex() && std::get<1>(*upper2) == prev->vertex()) {
                        ThrowRuntimeError("duplicate edge (%d-%d)", current->vertex(), prev->vertex());
                    }
                }
            }
        }
    }

    /// precompute the mapping from a face to a halfedge.
    void ConstructFaceMap(int nfaces) {
        face_edges_.resize(nfaces, nullptr);
        for (auto& edge : halfedges_) {
            face_edges_[edge.face()] = &edge;
        }
    }

    /// precompute the mapping from a vertex to a halfedge.
    void ConstructVertexMap(int nvertices) {
        vertex_edges_.resize(nvertices, nullptr);
        std::vector<bool> border(vertex_edges_.size(), false);
        int nerrors = 0;
        for (auto& edge : halfedges_) {
            int vid = edge.vertex();
            if (!edge.pair() || vertex_edges_[vid] == nullptr) { // prioritize boundary edges
                vertex_edges_[vid] = &edge;
            }
            if (!edge.pair()) { // border
                if (border[vid]) { // already visited
                    nerrors++;
                }
                border[vid] = true;
            }
        }
        if (nerrors) {
            std::clog << "warning: " << nerrors << " dangling vertices" << std::endl;
        }
        nerrors = 0;
        for (auto e : vertex_edges_) {
            if (!e) {
                nerrors++;
            }
        }
        if (nerrors) {
            std::clog << "warning: " << nerrors << " unused vertices" << std::endl;
        }
    }

private:
    std::vector<Halfedge> halfedges_; ///< halfedges
    std::vector<Halfedge *> vertex_edges_; ///< mapping from vertex to halfedge
    std::vector<Halfedge *> face_edges_; ///< mapping from face to halfedge
};


inline std::vector<int> GetAdjacentVertices(const HalfedgeMesh& hds, int vid)  {
    std::vector<int> out;
    auto e = hds.vertex_edge(vid);
    auto prev = e->previous();
    out.push_back(prev->vertex());
    auto first = e;
    while (e && e->next()->pair() != first) {
        out.push_back(e->next()->vertex());
        e = e->next()->pair();
    }
    return out;
}

inline std::vector<int> GetAdjacentFaces(const HalfedgeMesh& hds, int vid)  {
    std::vector<int> out;
    auto e = hds.vertex_edge(vid);
    auto first = e;
    do {
        out.push_back(e->face());
        e = e->next()->pair();
    } while (e != first && e);
    return out;
}

inline int GetDegree(const HalfedgeMesh& hds, int vid)  {
    auto e = hds.vertex_edge(vid);
    if (!e) {
        return 0;
    }
    int degree = 1;
    auto first = e;
    while (e) {
        e = e->next()->pair();
        if (e == first) {
            return degree;
        }
        degree++;
    }
    return degree;
}

inline bool IsBorder(const HalfedgeMesh& hds, int vid)  {
    return !hds.vertex_edge(vid)->pair();
}

/// halfedge with split/collapse functionality
class EditableHalfedgeMesh {
public:
    EditableHalfedgeMesh() = default;
    EditableHalfedgeMesh(const EditableHalfedgeMesh&) = delete;
    EditableHalfedgeMesh& operator=(const EditableHalfedgeMesh&) = delete;
    ~EditableHalfedgeMesh() {
        for (auto he : halfedges_) {
            delete he;
        }
    }

    int num_edges() const {
        return halfedges_.size();
    }

    int num_vertices() const {
        return vertex_edges_.size();
    }

    int num_faces() const {
        return face_edges_.size();
    }

    /// returns one of halfedges connecting to a vertex
    /// @return halfedge; the pair is null if the vertex is on the boundary
    Halfedge *vertex_edge(int vid) const {
        return vertex_edges_[vid];
    }

    /// returns one of halfedges owned by a face
    Halfedge *face_edge(int fid) const {
        return face_edges_[fid];
    }

    Halfedge *halfedge(int start, int end) const {
        Halfedge *e = vertex_edges_[start]->next();
        auto first = e;
        while (e->vertex() != end) {
            e = e->pair();
            if (!e || (e = e->next()) == first) {
                ThrowRuntimeError("non-existent edge");
            }
        }
        return e;
    }

    /// generate halfedge data structure
    /// @tparam Mesh polygon mesh
    template <typename Mesh>
    void Construct(const Mesh& mesh) {
        GenerateGraph(mesh);
        assert(IsConsistent(mesh));
        GenerateVertexMap(mesh.vertices.size());
        GenerateFaceMap(mesh.faces.size());
    }

    // split an edge into two, resulting in one additional vertex, and one or two additional faces for respectively border or non-border edges.
    void SplitEdge(Halfedge *he) {
        // -->[]-2->
        //    |^1 /
        //   0|| /
        // ...v|L
        halfedges_.push_back(new Halfedge);
        auto he0 = halfedges_.back();
        halfedges_.push_back(new Halfedge);
        auto he1 = halfedges_.back();
        halfedges_.push_back(new Halfedge);
        auto he2 = halfedges_.back();

        if (vertex_edges_[he->vertex()] == he) {
            vertex_edges_[he->vertex()] = he2;
        }
        if (face_edges_[he->face()] == he->next()) {
            face_edges_[he->face()] = he;
        }

        int newvid = num_vertices();
        int newfid = num_faces();
        face_edges_.push_back(he1);

        he0->next_ = he->next()->next();
        he0->pair_ = he1;
        he0->vertex_ = he->next()->vertex();
        he0->face_ = he->face();

        he1->next_ = he2;
        he1->pair_ = he0;
        he1->vertex_ = newvid;
        he1->face_ = newfid;

        he2->next_ = he->next();
        he2->vertex_ = he->vertex();
        he2->face_ = newfid;

        he->next()->next_ = he1;
        he->next()->face_ = newfid;

        he->next_ = he0;
        he->vertex_ = newvid;

        if (he->pair()) {
            he->pair()->pair_ = 0; // temporarily
            SplitEdge(he->pair()); // vertex_edges_ will be updated.
            he->pair()->pair_ = he2;
            he2->pair_ = he->pair();

            he->pair()->next()->pair()->next()->pair_ = he;
            he->pair_ = he->pair()->next()->pair()->next();

            assert(he->previous()->vertex() == he->pair()->vertex());
            assert(he0->previous()->vertex() == he0->pair()->vertex());
            assert(he1->previous()->vertex() == he1->pair()->vertex());
            assert(he2->previous()->vertex() == he2->pair()->vertex());
        } else {
            he2->pair_ = 0;
            vertex_edges_.push_back(he);
        }
    }

#if 0
    void CollapseEdge(const Halfedge *e) {
        // end vertex is replaced with start vertex
        int start_vid = e->prev()->vertex();
        auto ve = vertex_edges_[e->vertex()];
        auto first = ve;
        do {
            ve->vertex_ = start_vid;
            ve = ve->next()->pair_;
        } while (ve && ve != first);
        vertex_edges_[e->vertex()] = 0; // TODO: assuming non-dangling

        auto collapse_face = [&](Halfedge * p) {
            if (p->next()->next()->next() == p) {
                // triangular face is deleted
                face_edges_[p->face()] = 0;
                // reconnect edges
                auto n = p->next_;
                auto nn = n->next_;
                if (n->pair()) {
                    n->pair_->pair_ = nn->pair_;
                }
                if (nn->pair()) {
                    nn->pair_->pair_ = n->pair_;
                }
                // invalidate edges
                p->Invalidate();
                n->Invalidate();
                nn->Invalidate();
            } else {
                if (face_edges_[p->face()] == p) {
                    face_edges_[p->face()] = p->next_;
                }
                auto previous = p->next_;
                while (previous->next() != p) {
                    previous = previous->next_;
                }
                previous->next_ = p->next_;
                p->Invalidate();
            }
        };

        // paired face
        Halfedge *m = const_cast<Halfedge *>(e);
        if (m->pair()) {
            collapse_face(m->pair_);
        }
        // current face
        collapse_face(m);
    }
#endif
private:
    /// check if the halfedge structure is consistent (for debug)
    /// @tparam Mesh polygon mesh
    template <typename Mesh>
    bool IsConsistent(const Mesh& mesh) const {
        int num_unused = 0;
        for (auto he : halfedges_) {
            // pointers
            if (he->face() == -1 || !he->next() || he->vertex() == -1) {
                if (he->face() == -1 && !he->next() && he->vertex() == -1) {
                    num_unused++;
                } else {
                    ThrowRuntimeError("uninitialized edge");
                }
            }

            // face
            int fid = he->face();
            auto current = he;
            for (int e = 0; e < mesh.faces[fid].index.size(); e++) {
                if (current->face() != fid) {
                    ThrowRuntimeError("face inconsistent");
                }
                current = current->next();
            }
            if (current != he) {
                ThrowLogicError("face not closed");
            }

            // pair
            if (he->pair()) {
                if (he->pair()->pair() != he) {
                    ThrowRuntimeError("inconsistent halfedge pair");
                }
            }
        }
        if (num_unused) {
            std::clog << num_unused << " unused halfedges" << std::endl;
        }
        return true;
    }

    /// generate a halfedge data structure
    /// @tparam Mesh polygon mesh
    template <typename Mesh>
    void GenerateGraph(const Mesh& mesh) {
        // reserve memory
        int nhalfedges = 0;
        for (auto& face : mesh.faces) {
            nhalfedges += face.index.size(); // number of halfedges
        }
        for (auto he : halfedges_) {
            delete he;
        }
        halfedges_.clear();
        halfedges_.reserve(nhalfedges);

        std::vector<std::tuple<int, int, Halfedge *>> edges;
        edges.reserve(nhalfedges);

        // insert halfedges without pairs
        for (int fid = 0; fid < mesh.faces.size(); fid++) {
            auto& face = mesh.faces[fid];
            Halfedge *first, *prev = 0;
            for (int vid : face.index) {
                halfedges_.push_back(new Halfedge);
                auto current = halfedges_.back();

                // face
                current->face_ = fid;
                // vertex
                current->vertex_ = vid;
                // next edge
                if (prev) {
                    prev->next_ = current;
                    edges.push_back(std::make_tuple(prev->vertex(), current->vertex(), current));
                } else {
                    first = current;
                }
                prev = current;
            }

            // the next of the last is the first
            prev->next_ = first;
            edges.push_back(std::make_tuple(prev->vertex(), first->vertex(), first));
        }

        std::sort(edges.begin(), edges.end());

        // pair
        for (auto it = edges.begin(), end = edges.end(); it != end; ++it) {
            auto current = std::get<2>(*it);
            if (!current->pair()) {
                auto prev = current;
                while (prev->next() != current) {
                    prev = prev->next_;
                }
                auto upper = std::upper_bound(
                                 it + 1,
                                 end,
                                 std::make_tuple(current->vertex(), prev->vertex(), (Halfedge *)0));
                if (upper == edges.end() || std::get<0>(*upper) != current->vertex() || std::get<1>(*upper) != prev->vertex()) {
                    // border edge
                } else {
                    auto pair = std::get<2>(*upper);
                    current->pair_ = pair;
                    if (pair->pair()) {
                        ThrowRuntimeError("duplicate edge (%d-%d)", pair->vertex(), current->vertex());
                    }
                    pair->pair_ = current;
                    auto upper2 = upper + 1;
                    if (upper2 != end && std::get<0>(*upper2) == current->vertex() && std::get<1>(*upper2) == prev->vertex()) {
                        ThrowRuntimeError("duplicate edge (%d-%d)", current->vertex(), prev->vertex());
                    }
                }
            }
        }
    }

    /// precompute the mapping from a face to a halfedge.
    void GenerateFaceMap(int nfaces) {
        face_edges_.resize(nfaces, 0);
        for (auto he : halfedges_) {
            face_edges_[he->face()] = he; // last occurence
        }
    }

    /// precompute the mapping from a vertex to a halfedge.
    void GenerateVertexMap(int nvertices) {
        vertex_edges_.resize(nvertices, 0);
        std::vector<bool> border(vertex_edges_.size(), false);
        int nerrors = 0;
        for (auto he : halfedges_) {
            int vid = he->vertex();
            if (!he->pair() || !vertex_edges_[vid]) { // prioritize boundary edges
                vertex_edges_[vid] = he;
            }
            if (!he->pair()) { // border
                if (border[vid]) { // already visited
                    nerrors++;
                }
                border[vid] = true;
            }
        }
        if (nerrors) {
            std::clog << "warning: " << nerrors << " dangling vertices" << std::endl;
        }
        nerrors = 0;
        for (auto e : vertex_edges_) {
            if (!e) {
                nerrors++;
            }
        }
        if (nerrors) {
            std::clog << "warning: " << nerrors << " unused vertices" << std::endl;
        }
    }

private:
    std::vector<Halfedge *> halfedges_; // halfedges (stable vector)
    std::vector<Halfedge *> vertex_edges_; // mapping from vertex to halfedge
    std::vector<Halfedge *> face_edges_; // mapping from face to halfedge
};

} // namespace mesh
} // namespace slib
