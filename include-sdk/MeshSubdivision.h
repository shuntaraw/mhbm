// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <map>
#include "MeshBase.h"
#include "Halfedge.h"

namespace slib {
namespace mesh {

/// modified-butterfly subdivision surface
namespace mb {

/// helper function for modified-butterfly subdivision surface
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
CVector<float, 3> InsertBoundary(const Mesh& mesh, const hbm::Halfedge *edge) {
    // -b--a<-a--b-
    //  _\/_\/_\/_
    float a =  9.0 / 16;
    float b = -1.0 / 16;
    CVector<float, 3> pos {0, 0, 0};
    pos += a * mesh.vertices[edge->vertex()].position;

    auto e = edge->next();
    while (e->pair()) {
        e = e->pair()->next();
    }
    pos += b * mesh.vertices[e->vertex()].position;

    e = edge->next()->next();
    pos += a * mesh.vertices[e->vertex()].position;

    e = edge->next();
    while (e->next()->pair()) {
        e = e->next()->pair()->next();
    }
    pos += b * mesh.vertices[e->vertex()].position;
    return pos;
}

/// @tparam Mesh polygon mesh
template <typename Mesh>
static
CVector<float, 3> InsertBorder(
    const Mesh& mesh,
    const hbm::Halfedge *edge) {
    // '0' is boundary
    //   9...2...3
    //  .'. / \ .'.
    // 8...0<--1...4
    //  '.' \ / '.'
    //   7...6...5
    float w = 0;
    float a = 1.0 / 2 - w; // 0,1
    float b = 1.0 / 8 + 2 * w; // 2,6
    float c = -1.0 / 16 - w; // 3,5,7,9
    float d = w; // 4,8
    CVector<float, 3> pos {0, 0, 0};
    pos += a * mesh.vertices[edge->vertex()].position; // 0
    pos += a * mesh.vertices[edge->pair()->vertex()].position; // 1
    pos += b * mesh.vertices[edge->pair()->next()->vertex()].position; // 2
    pos += b * mesh.vertices[edge->next()->vertex()].position; // 6

    auto reflect = [](const CVector<float, 3>& a, const CVector<float, 3>& b, const CVector<float, 3>& c) {
        //    a
        // c  |  c'
        //    b
        CVector<float, 3> ab = normalized_of(b - a);
        return 2 * (dot(ab, c - a) * ab + a) - c;
    };


    if (edge->pair()->next()->pair()) {
        // '3' exists
        pos += c * mesh.vertices[edge->pair()->next()->pair()->next()->vertex()].position;
    } else {
        // generate '3'
        pos += c * reflect(
                   mesh.vertices[edge->pair()->vertex()].position,
                   mesh.vertices[edge->pair()->next()->vertex()].position,
                   mesh.vertices[edge->vertex()].position);
    }
    if (edge->next()->next()->pair()) {
        // '5' exists
        pos += c * mesh.vertices[edge->next()->next()->pair()->next()->vertex()].position;
    } else {
        // generate '5'
        pos += c * reflect(
                   mesh.vertices[edge->next()->vertex()].position,
                   mesh.vertices[edge->next()->next()->vertex()].position,
                   mesh.vertices[edge->vertex()].position);
    }
    if (edge->next()->pair()) {
        // '7' exists
        pos += c * mesh.vertices[edge->next()->pair()->next()->vertex()].position;
    } else {
        // generate '7'
        pos += c * reflect(
                   mesh.vertices[edge->vertex()].position,
                   mesh.vertices[edge->next()->vertex()].position,
                   mesh.vertices[edge->pair()->vertex()].position);
    }
    if (edge->pair()->next()->next()->pair()) {
        // '9' exists
        pos += c * mesh.vertices[edge->pair()->next()->next()->pair()->next()->vertex()].position;
    } else {
        // generate '9'
        pos += c * reflect(
                   mesh.vertices[edge->vertex()].position,
                   mesh.vertices[edge->pair()->next()->vertex()].position,
                   mesh.vertices[edge->pair()->vertex()].position);
    }
    return pos;
}

/// helper
/// @tparam Mesh polygon mesh
template <typename Mesh>
static
CVector<float, 3> InsertOrdinary(
    const Mesh& mesh,
    const hbm::Halfedge *edge) {
    //   9---2---3
    //  / \ / \ / \
    // 8---0<--1---4
    //  \ / \ / \ /
    //   7---6---5
    float w = 0;
    float a = 1.0 / 2 - w; // 0,1
    float b = 1.0 / 8 + 2 * w; // 2,6
    float c = -1.0 / 16 - w; // 3,5,7,9
    float d = w; // 4,8
    CVector<float, 3> pos {0, 0, 0};
    pos += a * mesh.vertices[edge->vertex()].position;
    pos += a * mesh.vertices[edge->pair()->vertex()].position;
    edge = edge->next();
    pos += b * mesh.vertices[edge->vertex()].position;
    edge = edge->pair()->next();
    pos += c * mesh.vertices[edge->vertex()].position;
    edge = edge->pair()->next();
    pos += d * mesh.vertices[edge->vertex()].position;
    edge = edge->pair()->next();
    pos += c * mesh.vertices[edge->vertex()].position;
    edge = edge->pair()->next();
    pos += b * mesh.vertices[edge->vertex()].position;
    edge = edge->pair()->next();
    pos += c * mesh.vertices[edge->vertex()].position;
    edge = edge->pair()->next();
    pos += d * mesh.vertices[edge->vertex()].position;
    edge = edge->pair()->next();
    pos += c * mesh.vertices[edge->vertex()].position;
    return pos;
}

/// helper
/// @tparam Mesh polygon mesh
template <typename Mesh>
static
CVector<float, 3> InsertExtraordinary(
    const Mesh& mesh, // polygon mesh
    const hbm::Halfedge *edge, // an edge on which a new vertex is generated
    int degree // degree of the extraordinary vertex
) {
    CVector<float, 3> pos = 0.75f * mesh.vertices[edge->vertex()].position;
    switch (degree) {
    case 3:
        //    1
        //  / | \
        // (  <--0
        //  \ | /
        //    1
        pos +=  5.0 / 12 * mesh.vertices[edge->pair()->vertex()].position;
        pos += -1.0 / 12 * mesh.vertices[edge->next()->vertex()].position;
        pos += -1.0 / 12 * mesh.vertices[edge->pair()->next()->vertex()].position;
        return pos;
    case 4:
        //    1
        //  / | \
        // 2--<--0
        //  \ | /
        //    1
        pos +=  3.0 / 8 * mesh.vertices[edge->pair()->vertex()].position;
        pos += -1.0 / 8 * mesh.vertices[edge->next()->pair()->next()->vertex()].position;
        return pos;
    default:
        //  k2---k1
        //  : \ / \
        //  :--<---0
        //  : / \ /
        //   2---1
        assert(degree > 4);
        const hbm::Halfedge *first = edge;
        for (int j = 0; j < degree; j++) {
            float s = (0.25 + cos(2 * M_PI * j / degree) + 0.5 * cos(4 * M_PI * j / degree)) / degree;
            pos += s * mesh.vertices[edge->pair()->vertex()].position;
            edge = edge->next()->pair();
        }
        assert(edge == first);
        return pos;
    }
}

/// helper
/// @tparam Mesh polygon mesh
template <typename Mesh>
static
CVector<float, 3> InsertVertex(
    const Mesh& mesh, // polygon mesh
    const hbm::Halfedge *edge // an edge on which a new vertex is generated
) {
    // boundary edge
    if (!edge->pair()) {
        return InsertBoundary(mesh, edge);
    }

    // destination vertex is on boundary
    int degree0 = 1;
    for (auto he = edge->next()->pair(); he != edge; he = he->next()->pair()) {
        degree0++;
        if (!he) {
            return InsertBorder(mesh, edge);
        }
    }

    // source vertex is on boundary
    int degree1 = 1;
    for (auto he = edge->pair()->next()->pair(); he != edge->pair(); he = he->next()->pair()) {
        degree1++;
        if (!he) {
            return InsertBorder(mesh, edge->pair());
        }
    }

    if (degree0 == 6) {
        // both are ordinary
        if (degree1 == 6) {
            return InsertOrdinary(mesh, edge);
        } else {
            // source vertex is extraordinary
            return InsertExtraordinary(mesh, edge->pair(), degree1);
        }
    } else {
        if (degree1 == 6) {
            // destination vertex is extraordinary
            return InsertExtraordinary(mesh,  edge, degree0);
        } else {
            // both are extraordinary
            return InsertExtraordinary(mesh,  edge, degree0) / 2 + InsertExtraordinary(mesh,  edge->pair(), degree1) / 2;
        }
    }
}

}

/// loop's subdivision surfaces
namespace loop {

/// helper
/// @tparam Mesh polygon mesh
template <typename Mesh>
static
CVector<float, 3> InsertVertex(
    const Mesh& mesh,
    const hbm::Halfedge *edge) {
    if (!edge->pair()) {
        // boundary
        return mesh.vertices[edge->vertex()].position / 2 + mesh.vertices[edge->next()->next()->vertex()].position / 2;
    } else {
        // non-boundary
        CVector<float, 3> pos {0, 0, 0};
        pos += 3.0 / 8 * mesh.vertices[edge->vertex()].position;
        pos += 3.0 / 8 * mesh.vertices[edge->pair()->vertex()].position;
        pos += 1.0 / 8 * mesh.vertices[edge->next()->vertex()].position;
        pos += 1.0 / 8 * mesh.vertices[edge->pair()->next()->vertex()].position;
        return pos;
    }
}

/// helper
/// @tparam Mesh polygon mesh
template <typename Mesh>
static
CVector<float, 3> UpdateVertex(
    const Mesh& mesh,
    const hbm::Halfedge *edge) {
    CVector<float, 3> pos {0, 0, 0};
    int degree = 1;
    for (auto he = edge->next()->pair(); he != edge; he = he->next()->pair()) {
        degree++;
        if (!he) {
            // boundary vertex
            pos += 3.0 / 4 * mesh.vertices[edge->vertex()].position;
            auto first = edge;
            while (first->pair()) {
                first = first->pair()->next()->next();
            }
            pos += 1.0 / 8 * mesh.vertices[first->next()->next()->vertex()].position;
            auto last = edge;
            while (last->next()->pair()) {
                last = last->next()->pair();
            }
            pos += 1.0 / 8 * mesh.vertices[last->next()->vertex()].position;
            return pos;
        }
    }
    // non-boundary vertex
    assert(degree > 2);
    float s;
    if (degree == 3) {
        s = 3.0 / 16;
    } else {
        s = (5.0 / 8 - pow(3.0 / 8 + 1.0 / 4 * cos(2 * M_PI / degree), 2)) / degree;
    }
    pos += (1 - s * degree) * mesh.vertices[edge->vertex()].position;
    for (int i = 0; i < degree; i++) {
        pos += s * mesh.vertices[edge->next()->vertex()].position;
        edge = edge->next()->pair();
    }
    return pos;
}

}

/// generate suvdivision surface by the modified butterfly algorithm.
/// the existing vertices remain unchanged.
/// the number of faces is increased by a factor of four.
/// the resultant mesh tends to suffer from high-frequency noise.
/// @see http://mrl.nyu.edu/~dzorin/papers/zorin1996ism.pdf
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
Mesh SubdivideModifiedButterfly(
    const Mesh& mesh // polygon mesh
) {
    hbm::HalfedgeMesh hds;
    hds.Construct(mesh);
    auto vertices = mesh.vertices;
    std::vector<typename Mesh::CFace> faces;
    faces.reserve(mesh.faces.size() * 4);
    // TODO: std::map is slow
    std::map<std::pair<int, int>, int> new_vid; // (lo_vid,hi_vid)->new_vid
    for (int fid = 0; fid < mesh.faces.size(); fid++) {
        auto& face = mesh.faces[fid];
        assert(face.index.size() == 3);
        auto he = hds.face_edge(fid);
        //      v0
        //     / \ 
        //    v3--v4
        //   / \ / \ 
        //  v2--v5--v1
        int vids[6] = {he->vertex(), he->next()->vertex(), he->next()->next()->vertex(), };
        for (int fv = 0; fv < 3; fv++, he = he->next()) {
            auto value = std::make_pair(
                             std::make_pair(he->vertex(), he->previous()->vertex()),
                             vertices.size());
            if (value.first.first > value.first.second) {
                std::swap(value.first.first, value.first.second);
            }
            auto m = new_vid.insert(value);
            if (m.second) {
                typename Mesh::CVertex vtx;
                vtx.position = mb::InsertVertex(mesh, he);
                vertices.push_back(std::move(vtx));
            }
            vids[3 + fv] =  m.first->second;
        }
        int triangle[4][3] = {
            {0, 4, 3},
            {1, 5, 4},
            {3, 4, 5},
            {2, 3, 5}
        };
        for (int t = 0; t < 4; t++) {
            typename Mesh::CFace face;
            face.resize(3);
            face.index[0] = vids[triangle[t][0]];
            face.index[1] = vids[triangle[t][1]];
            face.index[2] = vids[triangle[t][2]];
            faces.push_back(std::move(face));
        }
    }
    Mesh out;
    out.vertices = std::move(vertices);
    out.faces = std::move(faces);
    return out;
}

/// generate suvdivision surface by Loop's scheme.
/// the order of existing vertices remains unchanged.
/// the number of faces is increased by a factor of four.
/// @tparam Mesh polygon mesh
/// @see http://research.microsoft.com/en-us/um/people/cloop/thesis.pdf
template <typename Mesh>
inline
Mesh SubdivideLoop(
    const Mesh& mesh // polygon mesh
) {
    hbm::HalfedgeMesh hds;
    hds.Construct(mesh);
    auto vertices = mesh.vertices;

    // update existing coordinates
    for (int vid = 0; vid < mesh.vertices.size(); vid++) {
        auto he = hds.vertex_edge(vid);
        if (!he) {
            continue;    // unused vertex
        }
        vertices[vid].position = loop::UpdateVertex(mesh, he);
    }

    // generate vertices and faces
    std::vector<typename Mesh::CFace> faces;
    faces.reserve(mesh.faces.size() * 4);
    // TODO: std::map is slow
    std::map<std::pair<int, int>, int> new_vid; // (lo_vid,hi_vid)->new_vid
    for (int fid = 0; fid < mesh.faces.size(); fid++) {
        auto& face = mesh.faces[fid];
        assert(face.index.size() == 3);
        auto he = hds.face_edge(fid);
        //      v0
        //     / \ 
        //    v3--v4
        //   / \ / \ 
        //  v2--v5--v1
        int vids[6] = {he->vertex(), he->next()->vertex(), he->next()->next()->vertex(), };
        for (int fv = 0; fv < 3; fv++, he = he->next()) {
            auto value = std::make_pair(
                             std::make_pair(he->vertex(), he->previous()->vertex()),
                             vertices.size());
            if (value.first.first > value.first.second) {
                std::swap(value.first.first, value.first.second);
            }
            auto m = new_vid.insert(value);
            if (m.second) {
                typename Mesh::CVertex vtx;
                vtx.position = loop::InsertVertex(mesh, he);
                vertices.push_back(std::move(vtx));
            }
            vids[3 + fv] = m.first->second;
        }

        int triangle[4][3] = {
            {0, 4, 3},
            {1, 5, 4},
            {3, 4, 5},
            {2, 3, 5}
        };
        for (int t = 0; t < 4; t++) {
            typename Mesh::CFace face;
            face.resize(3);
            face.index[0] = vids[triangle[t][0]];
            face.index[1] = vids[triangle[t][1]];
            face.index[2] = vids[triangle[t][2]];
            faces.push_back(std::move(face));
        }
    }
    Mesh out;
    out.vertices = std::move(vertices);
    out.faces = std::move(faces);
    return out;
}

} // namespace mesh
} // namespace slib
