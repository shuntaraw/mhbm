// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <boost/range.hpp>
#include <vector>
#include <valarray>
#include <array>
#include "MatrixBase.h"

/// shared utilities
namespace slib {
/// MeshVertexmesh manipulation
namespace mesh {

/// vertex coordinate
struct CVertexPosition {
    CVector<float, 3> position = CVector<float, 3> {0, 0, 0} ; ///< 3D coordinate
};

/// vertex normal
struct CVertexNormal {
    CVector<float, 3> normal = CVector<float, 3> {0, 0, 0} ; ///< normal
};

/// vertex color
struct CVertexColor {
    CVector<unsigned char, 4> color = CVector<unsigned char, 4> {192, 192, 192, 255} ; ///< RGB color
};

/// texture coordinate
struct CVertexTexCoord {
    CVector<float, 2> texcoord = CVector<float, 2> {0, 0} ; ///< texture coordinate in [0,1]
};

/// triangle
struct CTriangleIndex {
    std::array<int, 3> index; ///< indices of vertices
    /// change the number of vertices (do nothing)
    void resize(size_t n) {
        if (n != 3) {
            ThrowRuntimeError("not a triangle");
        }
        return;
    }
};

/// polygon
struct CFaceIndex {
    std::valarray<int> index; ///< indices of vertices
    /// change the number of vertices (do nothing)
    void resize(int n) {
        index.resize(n);
    }
};

/// face normal
struct CFaceNormal {
    CVector<float, 3> normal; ///< normal vector
};

/// polygon mesh
/// @tparam CVertexType vertex
/// @tparam  CFaceType face
template <typename CVertexType = CVertexPosition, typename CFaceType = CTriangleIndex>
class CPolygonMesh {
public:
    typedef CVertexType CVertex;
    typedef CFaceType CFace;

    CPolygonMesh() {}

    CPolygonMesh(const CPolygonMesh&   m) :
        vertices(m.vertices),
        faces(m.faces) {}

    CPolygonMesh(CPolygonMesh&&  m) :
        vertices(std::move(m.vertices)),
        faces(std::move(m.faces)) {}

    CPolygonMesh& operator=(const CPolygonMesh& m) {
        if (this != &m) {
            vertices =  m.vertices;
            faces =  m.faces;
        }
        return *this;
    }

    CPolygonMesh& operator=(CPolygonMesh && m) {
        vertices = std::move(m.vertices);
        faces = std::move(m.faces);
        return *this;
    }

    std::vector<CVertexType> vertices;///< vertices
    std::vector<CFaceType> faces;///< faces
};

} // namespace mesh
} // namespace slib
