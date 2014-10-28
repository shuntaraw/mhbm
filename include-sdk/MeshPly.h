// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include "MeshBase.h"
#include "rply.h"

namespace slib {
namespace mesh {
/// loading ply meshes
namespace ply {

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex>
inline
int ReadVertexPosition(
    p_ply_argument argument ///<  rply argument
) {
    std::vector<CVertex> *vertices;
    long index;
    ply_get_argument_user_data(argument, reinterpret_cast<void **>(&vertices), &index);
    long instance_index;
    ply_get_argument_element(argument, 0, &instance_index);
    (*vertices)[instance_index].position[index] = ply_get_argument_value(argument);
    return 1;
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex>
inline
int ReadVertexNormal(
    p_ply_argument argument ///<  rply argument
) {
    std::vector<CVertex> *vertices;
    long index;
    ply_get_argument_user_data(argument, reinterpret_cast<void **>(&vertices), &index);
    long instance_index;
    ply_get_argument_element(argument, 0, &instance_index);
    (*vertices)[instance_index].normal[index] = ply_get_argument_value(argument);
    return 1;
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex>
inline
int ReadVertexColor(
    p_ply_argument argument ///<  rply argument
) {
    std::vector<CVertex> *vertices;
    long index;
    ply_get_argument_user_data(argument, reinterpret_cast<void **>(&vertices), &index);
    long instance_index;
    ply_get_argument_element(argument, 0, &instance_index);
    (*vertices)[instance_index].color[index] = ply_get_argument_value(argument);
    return 1;
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex>
inline
int ReadVertexTexCoord(
    p_ply_argument argument ///<  rply argument
) {
    std::vector<CVertex> *vertices;
    long index;
    ply_get_argument_user_data(argument, reinterpret_cast<void **>(&vertices), &index);
    long instance_index;
    ply_get_argument_element(argument, 0, &instance_index);
    (*vertices)[instance_index].texcoord[index] = ply_get_argument_value(argument);
    return 1;
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CFace>
inline
int ReadFaceIndex(
    p_ply_argument argument ///<  rply argument
) {
    std::vector<CFace> *faces;
    ply_get_argument_user_data(argument, reinterpret_cast<void **>(&faces), 0);
    long instance_index;
    ply_get_argument_element(argument, 0, &instance_index);
    p_ply_property property;
    long length, value_index;
    ply_get_argument_property(argument, &property, &length, &value_index);
    if (value_index == -1) {
        if (std::is_base_of<CTriangleIndex, CFace>::value) {
            if (length != 3) {
                ThrowRuntimeError("not a triangle");
            }
        } else {
            (*faces)[instance_index].resize(length);
        }
    } else {
        (*faces)[instance_index].index[value_index] = ply_get_argument_value(argument);
    }
    return 1;
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex, typename CFace>
inline
void SetReadVertexPositionFace(
    p_ply ply, ///<  rply argument
    std::vector<CVertex> *vertices, ///< [out] mesh
    std::vector<CFace> *faces ///< [out] mesh
) {
    int nvertices =
        ply_set_read_cb(ply, "vertex", "x", ReadVertexPosition<CVertex>, vertices, 0);
    ply_set_read_cb(ply, "vertex", "y", ReadVertexPosition<CVertex>, vertices, 1);
    ply_set_read_cb(ply, "vertex", "z", ReadVertexPosition<CVertex>, vertices, 2);
    vertices->resize(nvertices);
    int nfaces = ply_set_read_cb(ply, "face", "vertex_indices", ply::ReadFaceIndex<CFace>, faces, 0);
    faces->resize(nfaces);
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex>
inline
void SetReadVertexNormal(
    p_ply ply, ///<  rply argument
    std::vector<CVertex> *vertices, ///< [out] mesh
    typename std::enable_if < !std::is_base_of<CVertexNormal, CVertex>::value >::type * = 0
) {
    // do nothing
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex>
inline
void SetReadVertexNormal(
    p_ply ply, ///<  rply argument
    std::vector<CVertex> *vertices, ///< [out] mesh
    typename std::enable_if<std::is_base_of<CVertexNormal, CVertex>::value>::type * = 0
) {
    ply_set_read_cb(ply, "vertex", "nx", ReadVertexNormal<CVertex>, vertices, 0);
    ply_set_read_cb(ply, "vertex", "ny", ReadVertexNormal<CVertex>, vertices, 1);
    ply_set_read_cb(ply, "vertex", "nz", ReadVertexNormal<CVertex>, vertices, 2);
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex>
inline
void SetReadVertexColor(
    p_ply ply, ///<  rply argument
    std::vector<CVertex> *vertices, ///< [out] mesh
    typename std::enable_if < !std::is_base_of<slib::mesh::CVertexColor, CVertex>::value >::type * = 0
) {
    // do nothing
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex>
inline
void SetReadVertexColor(
    p_ply ply, ///<  rply argument
    std::vector<CVertex> *vertices, ///< [out] mesh
    typename std::enable_if<std::is_base_of<CVertexColor, CVertex>::value>::type * = 0
) {
    // VisualSFM
    ply_set_read_cb(ply, "vertex", "diffuse_red", ReadVertexColor<CVertex>, vertices, 0);
    ply_set_read_cb(ply, "vertex", "diffuse_green", ReadVertexColor<CVertex>, vertices, 1);
    ply_set_read_cb(ply, "vertex", "diffuse_blue", ReadVertexColor<CVertex>, vertices, 2);
    //
    ply_set_read_cb(ply, "vertex", "red", ReadVertexColor<CVertex>, vertices, 0);
    ply_set_read_cb(ply, "vertex", "green", ReadVertexColor<CVertex>, vertices, 1);
    ply_set_read_cb(ply, "vertex", "blue", ReadVertexColor<CVertex>, vertices, 2);
    ply_set_read_cb(ply, "vertex", "alpha", ReadVertexColor<CVertex>, vertices, 3);
    //
    ply_set_read_cb(ply, "vertex", "r", ReadVertexColor<CVertex>, vertices, 0);
    ply_set_read_cb(ply, "vertex", "g", ReadVertexColor<CVertex>, vertices, 1);
    ply_set_read_cb(ply, "vertex", "b", ReadVertexColor<CVertex>, vertices, 2);
    ply_set_read_cb(ply, "vertex", "a", ReadVertexColor<CVertex>, vertices, 3);
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex>
inline
void SetReadVertexTexCoord(
    p_ply ply, ///<  rply argument
    std::vector<CVertex> *vertices, ///< [out] mesh
    typename std::enable_if < !std::is_base_of<CVertexTexCoord, CVertex>::value >::type * = 0
) {
    // do nothing
}

/// callback for rply
/// @tparam Mesh polygon mesh
template <typename CVertex>
inline
void SetReadVertexTexCoord(
    p_ply ply, ///<  rply argument
    std::vector<CVertex> *vertices, ///< [out] mesh
    typename std::enable_if<std::is_base_of<CVertexTexCoord, CVertex>::value>::type * = 0
) {
    ply_set_read_cb(ply, "vertex", "u", ReadVertexTexCoord<CVertex>, vertices, 0);
    ply_set_read_cb(ply, "vertex", "v", ReadVertexTexCoord<CVertex>, vertices, 1);
}

/// callback for rply
/// @tparam CVertex vertex
template <typename CVertex>
inline
void ExportNormal(
    std::ofstream& ofs,
    const CVertex& v,
    typename std::enable_if < !std::is_base_of<CVertexNormal, CVertex>::value >::type * = 0) {}

/// callback for rply
/// @tparam CVertex vertex
template <typename CVertex>
inline
void ExportNormal(
    std::ofstream& ofs,
    const CVertex& v,
    typename std::enable_if<std::is_base_of<CVertexNormal, CVertex>::value>::type * = 0) {
    ofs.write((const char *)v.normal.data(), 12);
}

/// callback for rply
/// @tparam CVertex vertex
template <typename CVertex>
inline
void ExportColor(
    std::ofstream& ofs,
    const CVertex& v,
    typename std::enable_if < !std::is_base_of<CVertexColor, CVertex>::value >::type * = 0) {}

/// callback for rply
/// @tparam CVertex vertex
template <typename CVertex>
inline
void ExportColor(
    std::ofstream& ofs,
    const CVertex& v,
    typename std::enable_if<std::is_base_of<CVertexColor, CVertex>::value>::type * = 0) {
    ofs.write((const char *)v.color.data(), 4);
}
/// callback for rply
template <typename CVertex>
inline
void ExportTexCoord(
    std::ofstream& ofs,
    const CVertex& v,
    typename std::enable_if < !std::is_base_of<CVertexTexCoord, CVertex>::value >::type * = 0) {}

/// callback for rply
template <typename CVertex>
inline
void ExportTexCoord(
    std::ofstream& ofs,
    const CVertex& v,
    typename std::enable_if<std::is_base_of<CVertexTexCoord, CVertex>::value>::type * = 0) {
    ofs.write((const char *)v.texcoord().data(), 8);
}

} // namespace

/// import a mesh file in ply format
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void ReadPly(
    Mesh& mesh, ///< [out] mesh
    const std::string& filename ///<  filename
) {
    std::clog << "mesh <= " << filename << std::endl;
    p_ply ply = ply_open(filename.c_str(), 0, 0, 0);
    if (!ply) {
        ThrowRuntimeError("failed to open " + filename);
    }
    try {
        if (!ply_read_header(ply)) {
            ThrowRuntimeError("failed to parse " + filename);
        }

        mesh.vertices.clear();
        mesh.faces.clear();
        ply::SetReadVertexPositionFace(ply, &mesh.vertices, &mesh.faces);

        // x,y,z components of normal (if present)
        ply::SetReadVertexNormal(ply, &mesh.vertices);

        // vertex color (if present)
        ply::SetReadVertexColor(ply, &mesh.vertices);

        if (!ply_read(ply)) {
            ThrowRuntimeError("failed to read " + filename);
        }
    } catch (const std::exception&) {
        ply_close(ply);
        mesh.vertices.clear();
        mesh.faces.clear();
        throw;
    }
    ply_close(ply);
}

/// export a mesh file in ply format
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void WritePly(
    const Mesh& mesh, ///<  polygon mesh
    const std::string& filename ///<  filename
) {
    std::clog << "mesh => " << filename << std::endl;
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }

    ofs << "ply\n"
        << "format binary_little_endian 1.0\n"
        << "element vertex " << mesh.vertices.size() << "\n"
        << "property float x\n"
        << "property float y\n"
        << "property float z\n";

    if (std::is_base_of<CVertexNormal, typename Mesh::CVertex>::value) {
        ofs << "property float nx\n"
            << "property float ny\n"
            << "property float nz\n";
    }
    if (std::is_base_of<CVertexColor, typename Mesh::CVertex>::value) {
        ofs << "property uchar red\n"
            << "property uchar green\n"
            << "property uchar blue\n"
            << "property uchar alpha\n";
    }
    if (std::is_base_of<CVertexTexCoord, typename Mesh::CVertex>::value) {
        ofs << "property float u\n"
            << "property float v\n";
    }

    ofs << "element face " << mesh.faces.size() << "\n"
        << "property list uchar int vertex_indices\n"
        << "end_header\n";

    for (auto& v : mesh.vertices) {
        ofs.write((const char *)v.position.data(), sizeof(float) * 3);
        ply::ExportNormal(ofs, v);
        ply::ExportColor(ofs, v);
        ply::ExportTexCoord(ofs, v);
    }

    for (int i = 0; i < mesh.faces.size(); i++) {
        ofs.put(mesh.faces[i].index.size());
        for (int v = 0; v < mesh.faces[i].index.size(); v++) {
            int id = mesh.faces[i].index[v];
            ofs.write((const char *)&id, sizeof(int));
        }
    }
}
}
}
