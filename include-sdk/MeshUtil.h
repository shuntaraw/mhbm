// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <random>
#include <map>
#include "MatrixUtil.h"
#include "MeshBase.h"
#include "MeshLaplacian.h"
#include "Halfedge.h"

namespace slib {
namespace mesh {

/// mesh normal calculation
namespace normal {

/// update face normal
template <typename Mesh>
inline
void UpdateFaces(
    Mesh& mesh, ///< polygon mesh
    typename std::enable_if < !std::is_base_of<CFaceNormal, typename Mesh::CFace>::value >::type * = 0 ///< (dummy)
) {}

/// update face normal
template <typename Mesh>
inline
void UpdateFaces(
    Mesh& mesh, ///< polygon mesh
    typename std::enable_if<std::is_base_of<CFaceNormal, typename Mesh::CFace>::value>::type * = 0 ///< (dummy)
) {
    const int nvertices = mesh.vertices.size();
    int nfaces = mesh.faces.size();
    for (int fid = 0; fid < nfaces; fid++) {
        auto& face = mesh.faces[fid];
        assert(face.index.size() > 2);

        // average
        slib::CVector<float, 3> normal;
        normal.fill_with(0);
        for (int fvid = 0; fvid < face.index.size() - 2; fvid++) {
            auto v01 = mesh.vertices[face.index[fvid + 1]].position - mesh.vertices[face.index[0]].position;
            auto v02 = mesh.vertices[face.index[fvid + 2]].position - mesh.vertices[face.index[0]].position;
            auto n = cross(v01, v02);
            if (dot(n, n)) {
                normal += n; // weighted by area
            }
        }

        // normalize
        if (dot(normal, normal)) {
            normal = normalized_of(normal);
        }

        face.normal = normal;
    }
}

/// update vertex normal
template <typename Mesh>
inline
void UpdateVertices(
    Mesh& mesh, ///< polygon mesh
    typename std::enable_if < !std::is_base_of<CVertexNormal, typename Mesh::CVertex>::value >::type * = 0 ///< (dummy)
) {}

/// update vertex normal by the cotangent weight formula.
/// @see http://multires.caltech.edu/pubs/diffGeoOps.pdf
template <typename Mesh>
inline
void UpdateVertices(
    Mesh& mesh, ///< polygon mesh
    typename std::enable_if<std::is_base_of<CVertexNormal, typename Mesh::CVertex>::value>::type * = 0 ///< (dummy)
) {
    int nvertices = mesh.vertices.size();
    for (auto& v : mesh.vertices) {
        v.normal = {0, 0, 0};
    }

    // average face normals
    int nfaces = mesh.faces.size();
    if (!nfaces) {
        // TODO: to calculate normal from point clouds
        std::clog << "warning: no face" << std::endl;
        return;
    }
    int num_collapsed = 0;
    for (auto& face : mesh.faces) {
        int nv = face.index.size();
        int v0 = face.index[0];
        for (int vid = 0; vid < nv - 2; vid++) {
            int v1 = face.index[vid + 1];
            int v2 = face.index[vid + 2];
            auto v01 = mesh.vertices[v1].position - mesh.vertices[v0].position;
            auto v12 = mesh.vertices[v2].position - mesh.vertices[v1].position;
            auto v20 = mesh.vertices[v0].position - mesh.vertices[v2].position;

            auto n = cross(v01, v12);

            float cos0 = -dot(v01, v20);
            float sin0 = norm2_of(cross(v01, v20));
            float cos1 = -dot(v01, v12);
            float sin1 = norm2_of(cross(v01, v12));
            float cos2 = -dot(v12, v20);
            float sin2 = norm2_of(cross(v12, v20));
            if (sin0 == 0 || sin1 == 0 || sin2 == 0) {
                num_collapsed++;
                continue;
            }
            float cot0 = cos0 / sin0;
            float cot1 = cos1 / sin1;
            float cot2 = cos2 / sin2;

            float triangle = norm2_of(n) / 2;
            float a0, a1, a2; // Voronoi area
            if (cos0 < 0) {
                a0 = triangle / 2;
                a1 = triangle / 4;
                a2 = triangle / 4;
            } else if (cos1 < 0) {
                a0 = triangle / 4;
                a1 = triangle / 2;
                a2 = triangle / 4;
            } else if (cos2 < 0) {
                a0 = triangle / 4;
                a1 = triangle / 4;
                a2 = triangle / 2;
            } else {
                a0 = (cot1 * dot(v20, v20) + cot2 * dot(v01, v01)) / 8;
                a1 = (cot2 * dot(v01, v01) + cot0 * dot(v12, v12)) / 8;
                a2 = (cot0 * dot(v12, v12) + cot1 * dot(v20, v20)) / 8;
            }
            mesh.vertices[v0].normal += n * a0;
            mesh.vertices[v1].normal += n * a1;
            mesh.vertices[v2].normal += n * a2;
        }
    }

    if (num_collapsed) {
        std::cerr <<  "warning: " << num_collapsed << " triangles are collapsed and ignored." << std::endl;
    }

    // normalize
    for (auto& v : mesh.vertices) {
        auto norm2 = norm2_of(v.normal);
        if (norm2) {
            v.normal /= norm2;
        }
    }
}

}

/// update normals
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void UpdateNormal(Mesh& mesh) {
    normal::UpdateFaces<Mesh>(mesh);
    normal::UpdateVertices<Mesh>(mesh);
}

////////////////////////////////////////////////////////////////////// bounding box

/// compute upper and lower boundaries of a mesh
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void CalculateBoundingBox(
    const Mesh& mesh, ///< polygon mesh
    CVector<float, 3>& lower, ///< upper bound
    CVector<float, 3>& upper ///< lower bound
) {
    lower.fill_with(std::numeric_limits<float>::max());
    upper.fill_with(-std::numeric_limits<float>::max());
    for (auto& v : mesh.vertices) {
        for (int d = 0; d < 3; d++) {
            lower[d] = std::min(lower[d], v.position[d]);
            upper[d] = std::max(upper[d], v.position[d]);
        }
    }
}

/// compute the bounding sphere by Ritter's algorithm.
/// complexity = O(nd).
/// @see http://en.wikipedia.org/wiki/Bounding_sphere
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void CalculateBoundingSphere(
    const Mesh& mesh, ///< polygon mesh
    CVector<float, 3>& center, ///< center of bounding sphere
    float& radius ///< radius of bounding sphere
) {
    if (mesh.vertices.empty()) {
        std::clog << "warning: empty mesh" << std::endl;
        center = {0, 0, 0};
        radius = 0;
        return;
    }
    // x = random point
    auto& x = mesh.vertices[0].position;
    int nvertices = mesh.vertices.size();

    // y = farthest point from x
    int iy = 0;
    float maxd2 = 0;
    for (int i = 1; i < nvertices; i++) {
        float d2 = norm2_squared_of(mesh.vertices[i].position - x);
        if (maxd2 < d2) {
            maxd2 = d2;
            iy = i;
        }
    }
    auto& y = mesh.vertices[iy].position;

    // z = farthest point from y
    int iz = 0;
    for (int i = 1; i < nvertices; i++) {
        float d2 = norm2_squared_of(mesh.vertices[i].position - y);
        if (maxd2 < d2) {
            maxd2 = d2;
            iz = i;
        }
    }
    auto& z = mesh.vertices[iz].position;

    // initial estimates
    center = (y + z) / 2;
    radius = norm2_of(y - z) / 2;

    while (1) {
        // p = farthest point from current sphere
        int ip = -1;
        float tolerance = 1.001; // to avoid floating point error
        maxd2 = radius * radius * tolerance;
        for (int i = 0; i < nvertices; i++) {
            float d2 = norm2_squared_of(center - mesh.vertices[i].position);
            if (maxd2 < d2) {
                maxd2 = d2;
                ip = i;
            }
        }
        if (ip == -1) {
            break;
        } else {
            // p-------c'+---c---+
            //     d-r     r   r
            // c'= p + (c-p)(d+r)/d/2
            //   = c*(d+r)/(2d) + p*(d-r)/(2d)
            // r'= d-r
            float dist = sqrt(maxd2);
            center = ((dist - radius) * mesh.vertices[ip].position + (dist + radius) * center) / (2 * dist);
            radius = (dist + radius) / 2;
        }
    }
}

////////////////////////////////////////////////////////////////////// clean-up

/// eliminate faces.
/// vertices remain unchanged.
template <typename Mesh>
inline
void DeleteFaces(
    Mesh& mesh, ///< polygon mesh
    const std::vector<int>& fid_to_delete ///< list of face ids.
) {
    int noldfaces = mesh.faces.size();
    std::vector<bool> face_usage(noldfaces, true);
    for (int fid : fid_to_delete) {
        face_usage[fid] = false;
    }
    std::vector<Mesh::CFace> newfaces;
    int nnewfaces = std::count(face_usage.begin(), face_usage.end(), true);
    newfaces.reserve(nnewfaces);
    for (int fid = 0; fid < noldfaces; fid++) {
        if (face_usage[fid]) {
            newfaces.push_back(mesh.faces[fid]);
        }
    }
    if (mesh.faces.size() - newfaces.size()) {
        std::clog << mesh.faces.size() - newfaces.size() << " faces removed." << std::endl;
    }
    mesh.faces = std::move(newfaces);
}

/// eliminate vertices.
/// connected faces are also eliminated.
template <typename Mesh>
inline
void DeleteVertices(
    Mesh& mesh, ///< polygon mesh
    const std::vector<int>& vid_to_delete///< list of face ids.
) {
    int noldvertices = mesh.vertices.size();
    std::vector<bool> usage(noldvertices, true);
    for (auto vid : vid_to_delete) {
        usage[vid] = false;
    }
    std::vector<int> newvid(noldvertices);
    std::vector<typename Mesh::CVertex> newvertices;
    newvertices.reserve(noldvertices - vid_to_delete.size());
    for (int i = 0; i < noldvertices; i++) {
        if (usage[i]) {
            newvid[i] = newvertices.size();
            newvertices.push_back(mesh.vertices[i]);
        }
    }
    if (mesh.faces.size()) {
        std::vector<int> fid_to_delete;
        int noldfaces = mesh.faces.size();
        for (int fid = 0; fid < noldfaces; fid++) {
            auto& face = mesh.faces[fid];
            for (int vid : face.index) {
                if (!usage[vid]) {
                    fid_to_delete.push_back(fid);
                    break;
                }
            }
        }
        DeleteFaces(mesh, fid_to_delete);
        for (auto& face : mesh.faces) {
            for (int& vid : face.index) {
                vid = newvid[vid];
            }
        }
    }
    if (mesh.vertices.size() - newvertices.size()) {
        std::clog << mesh.vertices.size() - newvertices.size() << " vertices removed." << std::endl;
    }
    mesh.vertices = std::move(newvertices);
}

/// select unused vertices.
template <typename Mesh>
inline std::vector<int> SelectUnusedVertices(const Mesh& mesh) {
    int nvertices = mesh.vertices.size();
    std::vector<bool> vertex_usage(nvertices, false);
    for (auto& face : mesh.faces) {
        for (int vid : face.index) {
            vertex_usage[vid] = true;
        }
    }
    std::vector<int> selected;
    for (int vid = 0; vid < nvertices; vid++) {
        if (!vertex_usage[vid]) {
            selected.push_back(vid);
        }
    }
    return selected;
}

/// select small faces.
/// surface connectivity is preserved.
template <typename Mesh>
inline std::vector<int> SelectCollapsedFaces(
    const Mesh& mesh, ///< polygon mesh
    float min_area = 0,
    float min_edge_length = 0
) {
    int nfaces = mesh.faces.size();
    std::vector<int> selected;
    for (int fid = 0; fid < nfaces; fid++) {
        auto& face = mesh.faces[fid];
        // calculate area
        float area = 0;
        auto& v0 = mesh.vertices[face.index[0]].position;
        if (norm2_of(v0 - mesh.vertices[face.index[1]].position) < min_edge_length) {
            selected.push_back(fid);
            continue;
        }
        for (int fv = 0; fv < face.index.size() - 2; fv++) {
            auto& v1 = mesh.vertices[face.index[fv + 1]].position;
            auto& v2 = mesh.vertices[face.index[fv + 2]].position;
            if (norm2_of(v1 - v2) <= min_edge_length) {
                selected.push_back(fid);
                break;
            }
            auto cp = cross(v1 - v0, v2 - v1);
            area += norm2_of(cp);
        }
        // delete
        if (area <= min_area) {
            selected.push_back(fid);
        }
    }
    return selected;
}

template <typename Mesh>
inline std::vector<int> SelectNonManifoldVertices(const Mesh& mesh, bool allow_boundary = true) {
    int nfaces = mesh.faces.size();
    std::vector<int> selected;

    std::vector<std::tuple<int, int, int>> edges; // (vid, vid)->fid
    edges.reserve(3 * nfaces);
    for (int fid = 0; fid < nfaces; fid++) {
        auto& face = mesh.faces[fid];
        for (int fv = 0; fv < face.index.size(); fv++) {
            edges.push_back(
                std::make_tuple(face.index[fv], face.index[(fv + 1) % face.index.size()], fid)
            );
        }
    }

    auto comp = [](const std::tuple<int, int, int>& t1, const std::tuple<int, int, int>& t2) {
        if (std::get<0>(t1) < std::get<0>(t2)) {
            return true;
        } else if (std::get<0>(t1) > std::get<0>(t2)) {
            return false;
        } else {
            if (std::get<1>(t1) < std::get<1>(t2)) {
                return true;
            } else {
                return false;
            }
        }
    };

    std::sort(edges.begin(), edges.end(), comp);

    for (auto it1 = edges.begin(), it2 = it1 + 1, end = edges.end(); it2 != end; ++it1, ++it2) {
        if (std::get<0>(*it1) == std::get<1>(*it2) &&
            std::get<1>(*it1) == std::get<0>(*it2)) {
            selected.push_back(std::get<0>(*it1));
            selected.push_back(std::get<1>(*it1));
        }
    }

    if (!allow_boundary) {
        for (auto& e : edges) {
            if (std::find_if(edges.begin(), edges.end(), [&](const std::tuple<int, int, int>& t) {
            return
                std::get<0>(e) == std::get<1>(t) &&
                    std::get<1>(e) == std::get<0>(t);
            }) == edges.end()) {
                selected.push_back(std::get<0>(e));
                selected.push_back(std::get<1>(e));
            }
        }
    }
    return selected;
}

#if 0
template <typename Mesh>
inline std::vector<int> SelectNonManifoldFaces(const Mesh& mesh, bool allow_boundary = true) {
    std::map<std::pair<int, int>, int> edge; // (vid,vid)->fid
    int nfaces = mesh.faces.size();
    std::vector<int> selected;
    for (int fid = 0; fid < nfaces; fid++) {
        auto& face = mesh.faces[fid];
        for (int fv = 0; fv < face.index.size(); fv++) {
            auto e = std::make_pair(std::make_pair(face.index[fv], face.index[(fv + 1) % face.index.size()]), fid);
            auto it = edge.insert(e);
            if (!it.second) {
                // duplicated edge
                selected.push_back(fid);
                selected.push_back(it.first->second);
            }
        }
    }
    if (allow_boundary) {
        return selected;
    }
    for (auto& e : edge) {
        auto pair = std::make_pair(e.first.second, e.first.first);
        if (edge.find(pair) == edge.end()) {
            // boundary edge
            selected.push_back(e.second);
        }
    }
    return selected;
}
#endif

/// eliminate duplicated faces
/// @return true if deleted
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
std::vector<int> SelectDuplicatedFaces(Mesh& mesh) {
    std::vector<std::pair<std::tuple<int, int, int>, int>> faces; // v1,v2,v3 -> fid
    for (int fid = 0; fid < mesh.faces.size(); fid++) {
        auto& face = mesh.faces[fid];
        int begin = 0;
        int minvid = face.index[0];
        int nfv = face.index.size();
        for (int i = 1; i < nfv; i++) {
            if (minvid > face.index[i]) {
                minvid = face.index[i];
                begin = i;
            }
        }
        auto key = std::make_tuple(face.index[begin], face.index[(begin + 1) % nfv], 0); // 'v3' is dummy
        auto val = std::make_pair(key, fid);
        for (int i = 2; i < nfv; i++) {
            std::get<2>(val.first) = face.index[(begin + i) % nfv]; // update 'v3'
            faces.push_back(val);
        }
    }
    std::sort(faces.begin(), faces.end());
    std::vector<int> selected;
    for (int i = 0; i < faces.size() - 1; i++) {
        if (faces[i].first == faces[i + 1].first) {
            selected.push_back(faces[i].second);
        }
    }
    return selected;
}

/// find connected components
/// @tparam Mesh polygon mesh
/// @return {{vid,...},{vid,...},...}
template <typename Mesh>
inline
std::vector<std::vector<int>> FindComponents(const Mesh& mesh) {
    const int nvertices = mesh.vertices.size();
    const int nfaces = mesh.faces.size();

    std::vector<int> min_connected_vid(nvertices);
    for (int vid = 0; vid < nvertices; vid++) {
        min_connected_vid[vid] = vid;
    }

    std::function<int(int)> get_min_connected_vid = [&](int vid) {
        if (min_connected_vid[vid] < vid) {
            return get_min_connected_vid(min_connected_vid[vid]);
        } else {
            return vid;
        }
    };

    for (auto& face : mesh.faces) {
        int min_vid = nvertices;
        for (int vid : face.index) {
            min_vid = std::min(min_vid, get_min_connected_vid(vid));
        }
        for (int vid : face.index) {
            int old_min_vid = get_min_connected_vid(vid);
            if (old_min_vid != vid) {
                min_connected_vid[old_min_vid] = min_vid;
            }
            min_connected_vid[vid] = min_vid;
        }
    }

    // classify the vertices
    // TODO: std::map is slow
    std::map<int, std::vector<int>> components; // min_vid,{vid,...}
    for (int vid = 0; vid < nvertices; vid++) {
        int min_vid = get_min_connected_vid(vid);
        components[min_vid].push_back(vid);
    }

    std::vector<std::vector<int>> vertex_components;
    for (auto it = components.begin(), end = components.end(); it != end; ++it) {
        vertex_components.push_back(std::move(it->second));
#ifdef _DEBUG
        if (vertex_components.back().size() > 1) {
            std::clog << "component = " << vertex_components.size() << ": " << vertex_components.back().size() << " vertices" << std::endl;
        }
#endif
    }
    return vertex_components;
}

/// select small connected components
template <typename Mesh>
inline
std::vector<int> SelectSmallComponentFaces(
    const Mesh& mesh, ///< polygon mesh
    int max_num_components, ///< maximum number of components
    int min_num_vertices ///< minimum number of vertices in a component
) {
    auto vertex_component = FindComponents(mesh);
    std::vector<int> selected;
    if (vertex_component.size() == 1) {
        return selected;
    }

    // erase small components
    auto remove_end =
        std::remove_if(
            vertex_component.begin(),
            vertex_component.end(),
            [min_num_vertices](const std::vector<int>& v1)->bool { return v1.size() >= min_num_vertices; });

    int toremove = std::distance(remove_end, vertex_component.end()) - max_num_components;
    if (toremove > 0) {
        // sort in ascending order
        std::sort(
            remove_end,
            vertex_component.end(),
            [](const std::vector<int>& v1, const std::vector<int>& v2)->bool { return v1.size() < v2.size(); });
        remove_end += toremove;
    }
    int nremove = std::distance(vertex_component.begin(), remove_end);
    if (nremove < 1) {
        return selected;
    }
    std::clog << "removed " << nremove << " components" << std::endl;

    // list of faces to be removed
    std::vector<bool> removed_vertices(mesh.vertices.size(), false);
    for (auto it = vertex_component.begin(); it != remove_end; ++it) {
        for (auto v = it->begin(); v != it->end(); ++v) {
            removed_vertices[*v] = true;
        }
    }

    for (int fid = 0; fid < mesh.faces.size(); fid++) {
        for (int vid : mesh.faces[fid].index) {
            if (removed_vertices[vid]) {
                selected.push_back(fid);
                break;
            }
        }
    }
    return selected;
}

/// invert the orientation of faces
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void FlipFaces(Mesh& mesh) {
    int nfaces = mesh.faces.size();
    for (int fid = 0; fid < nfaces; fid++) {
        auto& face = mesh.faces[fid];
        int nfv = face.vertices.size();
        for (int i = 0; i < nfv / 2; i++) {
            std::swap(face.index[i], face.index[nfv - 1 - i]);
        }
    }
}

/// mesh transformation
namespace transform {

/// update vertex
template <typename Vertex, int NRows>
inline
void UpdateVertex(
    Vertex& vertex, ///<
    const CMatrix<float, NRows, 4>& mat, ///< affine transformation
    typename std::enable_if<std::is_base_of<CVertexNormal, Vertex>::value>::type * = 0
) {
    vertex.position = AffineTransform(mat, vertex.position);
    vertex.normal = RotateVector(mat, vertex.normal);
}

/// update vertex
template <typename Vertex, int NRows>
inline
void UpdateVertex(
    Vertex& vertex, ///<
    const CMatrix<float, NRows, 4>& mat, ///< affine transformation
    typename std::enable_if < !std::is_base_of<CVertexNormal, Vertex>::value >::type * = 0
) {
    vertex.position = AffineTransform(mat, vertex.position);
}

template <typename Face, int NRows>
inline
void UpdateFace(
    Face& face, ///<
    const CMatrix<float, NRows, 4>& mat, ///< affine transformation
    typename std::enable_if<std::is_base_of<CFaceNormal, Face>::value>::type * = 0
) {
    face.normal = RotateVector(mat, face.normal);
}

template <typename Face, int NRows>
inline
void UpdateFace(
    Face& face, ///<
    const CMatrix<float, NRows, 4>& mat, ///< affine transformation
    typename std::enable_if < !std::is_base_of<CFaceNormal, Face>::value >::type * = 0
) {}
}

/// affine transform mesh
/// @tparam Mesh polygon mesh
/// @tparam NRows number of rows of transformation matrix
template <typename Mesh, int NRows>
inline
void AffineTransform(
    Mesh& mesh, ///< polygon mesh
    const CMatrix<float, NRows, 4>& mat ///< affine transformation
) {
    //transform::TransformVertices(mesh, mat);
    //transform::TransformFaces(mesh, mat);
    for (auto& v : mesh.vertices) {
        transform::UpdateVertex(v, mat);
    }
    for (auto& f : mesh.faces) {
        transform::UpdateFace(f, mat);
    }
}

////////////////////////////////////////////////////////////////////// edit

/// append two meshes
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void Append(
    Mesh& dst, ///< base polygon mesh
    const Mesh& src ///< additional polygon mesh
) {
    int noldvertices = dst.vertices.size();
    int noldfaces = dst.faces.size();
    dst.vertices.insert(dst.vertices.end(), src.vertices.begin(), src.vertices.end());
    dst.faces.insert(dst.faces.end(), src.faces.begin(), src.faces.end());

    for (int fid = noldfaces; fid < dst.faces.size(); fid++) {
        for (int& vid : dst.faces[fid].index) {
            vid += noldvertices;
        }
    }
}

/// compute the volume of a polygon mesh
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
float CalcVolume(const Mesh& mesh) {
    std::vector<std::tuple<int, int>> edges;
    for (auto& face : mesh.faces) {
        for (int fvid = 0; fvid < face.size(); fvid++) {
            auto e = std::make_tuple(face.index[fvid], face.index[(fvid + 1) % face.size()]);
            if (std::get<0>(e) > std::get<1>(e)) {
                std::swap(std::get<0>(e), std::get<1>(e));
            }
            edges.push_back(e);
        }
    }
    std::sort(edges.begin(), edges.end());
    if (edges.size() % 2) {
        ThrowRuntimeError("odd number of halfedges");
    }
    for (int eid = 0; eid < edges.size(); eid += 2) {
        if (edges[eid] != edges[eid + 1]) {
            ThrowRuntimeError("invalid halfedge pair");
        }
    }
    float volume = 0;
    for (auto& face : mesh.faces) {
        auto& p0 = mesh.vertices[face.index[0]].position;
        for (int fvid = 0; fvid < face.size() - 2; fvid++) {
            auto& p1 = mesh.vertices[face.index[fvid + 1]].position;
            auto& p2 = mesh.vertices[face.index[fvid + 2]].position;
            auto mat =  make_matrix_from_column_vectors(p0, p1, p2);
            volume += determinant_of(mat) / 6;
        }
    }
    return volume;
}

template <typename Mesh>
inline
void AddGaussianNoise(Mesh& mesh, float stddev) {
    CMatrix<float> coordinate;
    ConvertMeshToCoordinate(mesh, coordinate);
    slib::CSparseMatrix<float> L;
    GetMeanCurvatureNormalLaplacian(mesh, L, false);
    auto normal = L.MultiplyTo('N', coordinate);
    std::mt19937 rng;
    rng.seed(0);
    std::normal_distribution<float> gen(0, stddev);
    for (int vid = 0; vid < mesh.vertices.size(); vid++) {
        mesh.vertices[vid].position += normalized_of(make_vector_from_row(normal, vid)) * gen(rng) ;
    }
}

/// Taubin mesh smoothing
///
/// @see http://www.cs.jhu.edu/~misha/Fall07/Papers/Taubin95.pdf
template <typename Mesh>
inline
void TaubinSmooth(Mesh& mesh, ///<
                  float lambda, ///< smoothing strengh in a single iteration in [0,1] (lambda)
                  int num_iterations ///< number of iteration
                 ) {
    assert(lambda > 0 && lambda < 1);
    float kpb = 0.1; // pass-band frequency (see the paper above)
    float mu = 1.0 / (kpb - 1 / lambda);
    CMatrix<float> coordinate, curvature_normal;
    ConvertMeshToCoordinate(mesh, coordinate);
    CSparseMatrix<float> L;
    for (int i = 0; i < num_iterations; i++) {
        GetMeanValueCoordinateLaplacian(mesh, L);
        curvature_normal = L.MultiplyTo('N', coordinate);
        coordinate -= lambda * curvature_normal;
        ConvertCoordinateToMesh(coordinate, mesh);

        // un-shrink
        GetMeanValueCoordinateLaplacian(mesh, L);
        curvature_normal = L.MultiplyTo('N', coordinate);
        coordinate -= mu * curvature_normal;
        ConvertCoordinateToMesh(coordinate, mesh);
    }
}

/// convert mesh coordinates to an mx3 coordinate matrix
/// @tparam Mesh polygon mesh
/// @tparam T type of each differencial coordinate
/// @return  coordinate matrix
template <typename Mesh, typename T /*= float*/>
inline
void ConvertMeshToCoordinate(
    const Mesh& mesh, ///< polygon mesh
    CMatrix<T>& coordinates  ///< mx3 coordinate matrix
) {
    coordinates.resize(mesh.vertices.size(), 3);
    for (int r = 0; r < mesh.vertices.size(); r++) {
        for (int c = 0; c < 3; c++) {
            coordinates(r, c) = mesh.vertices[r].position[c];
        }
    }
}

/// copy an mx3 coordinate matrix to mesh coordinates
/// @tparam Mesh polygon mesh
/// @tparam T coordinate
template <typename Mesh, typename T>
inline
void ConvertCoordinateToMesh(
    const CMatrix<T>& coordinates, ///< mx3 coordinate matrix
    Mesh& mesh ///< polygon mesh. vertex coordinates will be updated.
) {
    assert(coordinates.num_rows() == mesh.vertices.size());
    assert(coordinates.num_cols() == 3);
    for (int r = 0; r < coordinates.num_rows(); r++) {
        mesh.vertices[r].position = make_vector_from_row(coordinates, r) ;
    }
}

template <typename Mesh1, typename Mesh2>
inline
bool IsSameTopology(
    const Mesh1& m1,
    const Mesh2& m2
) {
    if (m1.vertices.size() != m2.vertices.size()) {
        return false;
    }
    int nfaces = m1.faces.size();
    if (nfaces != m2.faces.size()) {
        return false;
    }
    for (int fid = 0; fid < nfaces; fid++) {
        auto& f1 = m1.faces[fid];
        auto& f2 = m2.faces[fid];
        int nv = f1.index.size();
        if (nv != f2.index.size()) {
            return false;
        }
        for (int v = 0; v < nv; v++) {
            if (f1.index[v] != f2.index[v]) {
                return false;
            }
        }
    }
    return true;
}

} // namespace mesh
} // namespace slib
