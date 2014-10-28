// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <sstream>
#include <boost/timer/timer.hpp>
#include <boost/algorithm/string/replace.hpp>
#include "MatrixUtil.h"
#include "MeshBase.h"
#include "MeshPly.h"
#include "MeshObj.h"

namespace slib {
namespace mesh {

/// import a mesh file in off format
/// @tparam CMesh polygon mesh
/// @see http://www.geom.uiuc.edu/software/geomview/geomview_6.html#SEC36
template <typename CMesh>
inline
void ReadOff(CMesh& mesh, const std::string& filename) {
    std::clog << "mesh <= " << filename << std::endl;
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }

    char buf[BUFSIZ];
    auto read_word = [&]() {
        while (1) {
            ifs >> buf;
            if (ifs.fail()) {
                ThrowRuntimeError("premature EOF: " + filename);
            }
            if (buf[0] == '#') {
                ifs.getline(buf, BUFSIZ);
            } else {
                return buf;
            }
        }
    };

    // header
    bool has_ST = false;
    bool has_C = false;
    bool has_N = false;
    bool has_n = false;
    int dimension = 3;
    char *p = read_word();
    if (p[0] == 'S' && p[1] == 'T') {
        has_ST = true;
        p += 2;
    }
    if (p[0] == 'C') {
        has_C = true;
        p++;
    }
    if (p[0] == 'N') {
        has_N = true;
        p++;
    }
    if (strncmp(p, "4OFF", 4) == 0) {
        dimension = 4;
        read_word(); // nvertices
    } else if (strncmp(p, "nOFF", 4) == 0) {
        dimension = atoi(read_word());
        read_word(); // nvertices
    } else if (strncmp(p, "4nOFF", 5) == 0) {
        dimension = atoi(read_word()) + 1;
        read_word(); // nvertices
    } else if (strncmp(p, "OFF", 3) == 0) {
        read_word(); // nvertices
    }
    if (dimension != 3 && dimension != 4) {
        ThrowRuntimeError("unsupported dimension: " + filename);
    }

    // number of elements
    int nvertices = atoi(buf);
    if (!nvertices) {
        ThrowRuntimeError("no vertex found: " + filename);
    }
    int nfaces = atoi(read_word());
    int nedges = atoi(read_word()); // unused

    mesh.vertices.assign(nvertices, typename CMesh::CVertex());
    mesh.faces.assign(nfaces, typename CMesh::CFace());

    // vertices
    for (int i = 0; i < nvertices; i++) {
        // coordinates
        for (int d = 0; d < 3; d++) {
            mesh.vertices[i].position[d] = atof(read_word());
        }
        if (dimension == 4) {
            float f = atof(read_word());
            if (f == 0) {
                ThrowRuntimeError("infine coordinate: " + filename);
            }
            mesh.vertices[i].position /= f;
        }
        // discard others
        if (has_N) {
            read_word();
            read_word();
            read_word();
        }
        if (has_C) {
            read_word();
            read_word();
            read_word();
        }
        if (has_ST) {
            read_word();
            read_word();
        }
    }

    for (int i = 0; i < nfaces; i++) {
        int num = atoi(read_word());
        mesh.faces[i].resize(num);
        for (int j = 0; j < num; j++) {
            mesh.faces[i].index[j] = atof(read_word());
        }
        ifs.getline(buf, BUFSIZ); // discard face color
    }
    if (ifs.fail()) {
        ThrowRuntimeError("premature EOF: " + filename);
    }
}

/// export a mesh file in off format
/// @tparam CVertex vertex
/// @tparam CFace face
/// @see http://www.geom.uiuc.edu/software/geomview/geomview_6.html#SEC36
template <typename CVertex, typename CFace>
inline
void WriteOff(const CPolygonMesh<CVertex, CFace>& mesh, const std::string& filename) {
    std::clog << "mesh => " << filename << std::endl;
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }
    ofs << "OFF\n";
    ofs << mesh.vertices.size() << " " << mesh.faces.size() << " 0\n";
    for (auto& v : mesh.vertices) {
        ofs << v.position[0] << " " << v.position[1] << " " << v.position[2] << "\n";
    }
    for (int i = 0; i < mesh.faces.size(); i++) {
        ofs << mesh.faces[i].index.size() << " ";
        for (int j = 0; j < mesh.faces[i].index.size(); j++) {
            ofs << mesh.faces[i].index[j] << " ";
        }
        ofs << "\n";
    }
}

/// import a mesh file in geo format
/// @tparam CMesh polygon mesh
/// @see http://lc.cray.com/doc/movie/
template <typename CMesh>
inline
void ReadGeo(CMesh& mesh, const std::string& filename) {
    std::clog <<  "mesh <= " << filename << std::endl;
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }
    int nparts, nvertices, nelements, nedges;
    ifs >> nparts >> nvertices >> nelements >> nedges;
    if (ifs.fail() || !nvertices) {
        ThrowRuntimeError("invalid header: " + filename);
    }
    int ebegin, eend;
    ifs >> ebegin >> eend;
    if (eend != nelements + 1) {
        ThrowRuntimeError("multiple parts not supported");
    }
    mesh.vertices.assign(nvertices, typename CMesh::CVertex());
    for (int i = 0; i < nvertices; i++) {
        ifs >> mesh.vertices[i].position[0] >> mesh.vertices[i].position[1] >> mesh.vertices[i].position[2];
    }
    mesh.faces.clear();
    while (nedges > 0) {
        typename CMesh::CFace face;
        face.resize(3);
        ifs >> face.index[0] >> face.index[1] >> face.index[2];
        nedges -= 3;
        face.index[0]--;
        face.index[1]--;
        face.index[2]--;
        if (ifs.fail()) {
            ThrowRuntimeError("premature EOF: " + filename);
        }
        while (1) {
            if (face.index[2] < 0) {
                face.index[2] = -face.index[2] - 2;
                mesh.faces.push_back(std::move(face));
                break;
            } else {
                mesh.faces.push_back(face);
                std::clog <<  "splitting a non-triangle face." << std::endl;
                face.index[1] = face.index[2];
                ifs >> face.index[2];
                nedges--;
                face.index[2]--;
                if (ifs.fail()) {
                    ThrowRuntimeError("premature EOF: " + filename);
                }
            }
        }
    }
}

/// export a mesh file in geo format
/// @tparam CMesh polygon mesh
/// @see http://lc.cray.com/doc/movie/
template <typename CMesh>
inline
void WriteGeo(const CMesh& mesh, const std::string& filename) {
    std::clog << "mesh => " << filename << std::endl;
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }
    int nedges = 0;
    for (auto& f : mesh.faces) {
        nedges += f.index.size();
    }
    ofs << "1 " << mesh.vertices.size() << " " << mesh.faces.size() << " " << nedges << std::endl;
    ofs << "1 " << mesh.faces.size() + 1 << std::endl;
    for (int vid = 0; vid < mesh.vertices.size(); vid++) {
        ofs << mesh.vertices[vid].position[0] << " " << mesh.vertices[vid].position[1] << " " << mesh.vertices[vid].position[2] << " ";
        if (vid % 2) {
            ofs << std::endl;
        }
    }
    int nvalues = 0;
    for (auto& f : mesh.faces) {
        int nv = f.index.size();
        for (int v = 0; v < nv; v++) {
            if (v == nv - 1) {
                ofs << -f.index[v] - 1 << " ";
            } else {
                ofs << f.index[v] + 1 << " ";
            }
            nvalues++;
            if (nvalues == 10) {
                nvalues = 0;
                ofs << std::endl;
            }
        }
    }
}

/// import a mesh file in ppd format
/// @tparam CMesh polygon mesh
template <typename CMesh>
inline
void ReadPpd(CMesh& mesh, const std::string& filename) {
    std::clog << "mesh <= " << filename << std::endl;
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }
    int nvertices = 0, ntriangles = 0;
    std::string label;
    ifs >> label; // header
    if (label != "header") {
        ThrowRuntimeError("invalid format: " + filename);
    }
    while (1) {
        ifs >> label;
        if (ifs.fail()) {
            ThrowRuntimeError("invalid format: " + filename);
        } else if (label == "vertex") {
            ifs >> nvertices;
        } else if (label == "face") {
            ifs >> ntriangles;
        } else if (label == "end") {
            break;
        } else {
            int dummy;
            ifs >> dummy;
        }
    }
    if (!nvertices) {
        ThrowRuntimeError("invalid format: " + filename);
    }

    while (1) {
        ifs >> label;
        if (ifs.fail()) {
            ThrowRuntimeError("invalid format: " + filename);
        } else if (label == "vertex") {
            break;
        }
    }

    mesh.vertices.clear();
    mesh.vertices.resize(nvertices);
    for (int i = 0; i < nvertices; i++) {
        int id;
        CVector<float, 3> pos;
        ifs >> id >> pos[0] >> pos[1] >> pos[2];
        mesh.vertices[id - 1].position = pos;
    }
    ifs >> label;
    if (label != "end") {
        ThrowRuntimeError("invalid format: " + filename);
    }

    mesh.faces.resize(ntriangles);
    while (1) {
        ifs >> label;
        if (ifs.fail()) {
            ThrowRuntimeError("invalid format: " + filename);
        } else if (label == "face") {
            break;
        }
    }
    for (int i = 0; i < ntriangles; i++) {
        int id;
        mesh.faces[i].resize(3);
        int index[3];
        ifs >> id >> index[0] >> index[1] >> index[2];
        mesh.faces[i].index[0] = index[0] - 1;
        mesh.faces[i].index[1] = index[1] - 1;
        mesh.faces[i].index[2] = index[2] - 1;
    }
    ifs >> label;
    if (label != "end") {
        ThrowRuntimeError("invalid format: " + filename);
    }
}

/// import a mesh file in txt format
/// @tparam CMesh polygon mesh
template <typename CMesh>
inline
void ReadTxt(
    CMesh& mesh, // mesh
    const std::string& filename// filename
) {
    std::clog <<  "mesh <= " << filename << std::endl;
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }
    mesh.vertices.clear();
    mesh.faces.clear();
    while (1) {
        typename CMesh::CVertex vtx;
        ifs >> vtx.position[0] >> vtx.position[1] >> vtx.position[2];
        if (ifs.fail()) {
            break;
        }
        mesh.vertices.push_back(vtx);
    }
}

template <typename CMesh>
inline
void WriteTxt(
    const CMesh& mesh, // mesh
    const std::string& filename// filename
) {
    std::clog << "mesh => " << filename << std::endl;
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }
    for (auto& v : mesh.vertices) {
        ofs << v.position[0] << " " << v.position[1] << " " << v.position[2] << std::endl;
    }
}

template <typename Mesh>
inline
void WriteStl(const Mesh& mesh, const std::string& filename) {
    std::clog << "mesh => " << filename << std::endl;
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }
    for (int i = 0; i < 80; i++) {
        ofs.put(0);
    }
    unsigned int ntriangles = 0;
    for (auto& f : mesh.faces) {
        ntriangles += f.index.size() - 2;
    }
    ofs.write((const char *)&ntriangles, 4);
    auto write_normal = [](std::ofstream & ofs, const CVector<float, 3>& v0, const CVector<float, 3>& v1, const CVector<float, 3>& v2) {
        CVector<float, 3> normal = cross(v1 - v0, v2 - v0);
        if (dot(normal, normal)) {
            normal = normalized_of(normal);
        }
        ofs.write((const char *)normal.data(), 12);
    };
    for (auto& f : mesh.faces) {
        auto& v0 = mesh.vertices[f.index[0]].position;
        for (int fvid = 0; fvid < f.index.size() - 2; fvid++) {
            auto& v1 = mesh.vertices[f.index[fvid + 1]].position;
            auto& v2 = mesh.vertices[f.index[fvid + 2]].position;
            write_normal(ofs, v0, v1, v2);
            ofs.write((const char *)v0.data(), 12);
            ofs.write((const char *)v1.data(), 12);
            ofs.write((const char *)v2.data(), 12);
            ofs.put(0);
            ofs.put(0);
        }
    }
}

template <typename CVertexType, typename CFaceType>
inline
void Read(CPolygonMesh<CVertexType, CFaceType>& mesh, const std::string& filename) {
    auto ext = strrchr(filename.c_str(), '.');
    if (!ext) {
        ThrowRuntimeError("unsupported format: " + filename);
    } else if (!stricmp(ext, ".ply")) {
        ReadPly(mesh, filename);
    } else if (!stricmp(ext, ".obj") || !stricmp(ext, ".smf")) {
        ReadObj(mesh, filename);
    } else if (!stricmp(ext, ".off")) {
        ReadOff(mesh, filename);
    } else if (!stricmp(ext, ".csv") || !stricmp(ext, ".txt")) {
        ReadTxt(mesh, filename);
    } else if (!stricmp(ext, ".geo")) {
        ReadGeo(mesh, filename);
    } else if (!stricmp(ext, ".ppd")) {
        ReadPpd(mesh, filename);
    } else {
        ThrowRuntimeError("unsupported format:" + filename);
    }
    std::clog << "#vertices = " << mesh.vertices.size() << ", #faces = " << mesh.faces.size() << std::endl;
}

template <typename CVertexType, typename CFaceType>
inline
void Write(const CPolygonMesh<CVertexType, CFaceType>& mesh, const std::string& filename) {
    auto ext = strrchr(filename.c_str(), '.');
    if (!ext) {
        ThrowRuntimeError("unsupported format: " + filename);
    } else if (!stricmp(ext, ".ply")) {
        WritePly(mesh, filename);
    } else if (!stricmp(ext, ".obj")) {
        WriteObj(mesh, filename);
    } else if (!stricmp(ext, ".off")) {
        WriteOff(mesh, filename);
    } else if (!stricmp(ext, ".csv") || !stricmp(ext, ".txt")) {
        WriteTxt(mesh, filename);
    } else if (!stricmp(ext, ".stl")) {
        WriteStl(mesh, filename);
    } else if (!stricmp(ext, ".geo")) {
        WriteGeo(mesh, filename);
    } else {
        ThrowRuntimeError("unsupported format:" + filename);
    }
}

} // namespace mesh
} // namespace slib
