// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#include "precompile.h"
#include "Halfedge.h"

namespace hbm {

Halfedge *Halfedge::pair()const {
    return pair_;
}

Halfedge *Halfedge::next()const {
    return next_;
}

int Halfedge::vertex()const {
    return vertex_;
}

int Halfedge::face()const {
    return face_;
}

Halfedge *Halfedge::previous()const {
    auto prev = next();
    while (prev->next() != this) {
        prev = prev->next();
    }
    return prev;
}

void Halfedge::Invalidate() {
    pair_ = 0;
    next_ = 0;
    vertex_ = -1;
    face_ = -1;
}



int HalfedgeMesh::num_edges() const {
    return halfedges_.size();
}

int HalfedgeMesh::num_vertices() const {
    return vertex_edges_.size();
}

int HalfedgeMesh::num_faces() const {
    return face_edges_.size();
}

/// returns one of halfedges connecting to a vertex
/// @return halfedge; the pair is null if the vertex is on the boundary
Halfedge *HalfedgeMesh::vertex_edge(int vid) const {
    return vertex_edges_[vid];
}

/// returns one of halfedges owned by a face
Halfedge *HalfedgeMesh::face_edge(int fid) const {
    return face_edges_[fid];
}

/// generate halfedge data structure
/// @tparam Mesh polygon mesh
void HalfedgeMesh::Construct(const CMesh& mesh) {
    ConstructHalfedges(mesh);
    assert(IsConsistent(mesh));
    ConstructVertexMap(mesh.vertices.size());
    ConstructFaceMap(mesh.faces.size());
}

bool HalfedgeMesh::IsConsistent(const CMesh& mesh) const {
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
void HalfedgeMesh::ConstructHalfedges(
    const CMesh& mesh ///< [in] polygon mesh
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
void HalfedgeMesh::ConstructFaceMap(int nfaces) {
    face_edges_.resize(nfaces, nullptr);
    for (auto& edge : halfedges_) {
        face_edges_[edge.face()] = &edge;
    }
}

/// precompute the mapping from a vertex to a halfedge.
void HalfedgeMesh::ConstructVertexMap(int nvertices) {
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



std::vector<int> GetAdjacentVertices(const HalfedgeMesh& hds, int vid)  {
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

std::vector<int> GetAdjacentFaces(const HalfedgeMesh& hds, int vid)  {
    std::vector<int> out;
    auto e = hds.vertex_edge(vid);
    auto first = e;
    do {
        out.push_back(e->face());
        e = e->next()->pair();
    } while (e != first && e);
    return out;
}

int GetDegree(const HalfedgeMesh& hds, int vid)  {
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

bool IsBorder(const HalfedgeMesh& hds, int vid)  {
    return !hds.vertex_edge(vid)->pair();
}

} // hbm
