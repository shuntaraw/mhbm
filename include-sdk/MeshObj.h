#pragma once
#include "MeshBase.h"

namespace slib {
namespace mesh {
/// loading obj meshes
namespace obj {

/// callback for obj
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void ExportTexCoord(
    std::ofstream& ofs,
    const Mesh& mesh,
    typename std::enable_if < !std::is_base_of<CVertexTexCoord, typename Mesh::CVertex>::value >::type * = 0) {}

/// callback for obj
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void ExportTexCoord(
    std::ofstream& ofs,
    const Mesh& mesh,
    typename std::enable_if<std::is_base_of<CVertexTexCoord, typename Mesh::CVertex>::value>::type * = 0) {
    for (auto& v : mesh.vertices) {
        ofs << "vt " << v.texcoord(0) << " " << v.texcoord(1) << "\n";
    }
}

} // namespace obj

/// import a mesh file in obj format
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void ReadObj(
    Mesh& mesh, ///< [out] mesh
    const std::string& filename ///<  filename
) {
    std::cout << "mesh <= " << filename << std::endl;
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }

    mesh.vertices.clear();
    mesh.faces.clear();
    while (1) {
        char buf[BUFSIZ];
        ifs.getline(buf, BUFSIZ);

        if (ifs.fail()) {
            break;
        }

        std::istringstream str(buf);
        char token[BUFSIZ];
        str >> token;

        if (str.fail()) {
            // empty line
        } else if (!strcmp(token, "v")) {
            mesh.vertices.push_back(typename Mesh::CVertex());
            str >> mesh.vertices.back().position[0] >> mesh.vertices.back().position[1] >> mesh.vertices.back().position[2];
        } else if (!strcmp(token, "f")) {
            typename Mesh::CFace face;
            face.resize(3);
            for (int i = 0; i < 3; i++) {
                str >> face.index[i];
                face.index[i]--;
                if (str.fail()) {
                    ThrowRuntimeError("premature end of face indices");
                }
                if (str.peek() == '/') {
                    char unused[16];
                    str >> unused;
                }
            }
            while (1) {
                mesh.faces.push_back(face);

                face.index[1] = face.index[2];
                str >> face.index[2];
                face.index[2]--;
                if (str.fail()) {
                    break; // end of indices
                }
                std::clog << "triangulating face #" << mesh.faces.size() << std::endl;
                if (str.peek() == '/') {
                    char unused[16];
                    str >> unused;
                }
            }
        }
    }
}

/// export a mesh file in obj format
/// @tparam Mesh polygon mesh
template <typename Mesh>
inline
void WriteObj(
    const Mesh& mesh, ///<  polygon mesh
    const std::string& filename ///<  filename
) {
    std::clog << "mesh => " << filename << std::endl;
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        ThrowRuntimeError("failed to open " + filename);
    }

    if (std::is_base_of<CVertexTexCoord, typename Mesh::CVertex>::value) {
        auto mtl = boost::replace_all_copy(filename, ".obj", ".mtl");
        auto png = boost::replace_all_copy(filename, ".obj", ".bmp");
        std::ofstream ofs2(mtl);
        ofs2 << "newmtl material_0\n"
             "Ka 1 1 1\n"
             "Kd 1 1 1\n"
             "Ks 1 1 1\n"
             "map_Kd " << png << "\n";
        ofs << "mtllib " << mtl << "\n"
            << "usemtl material_0\n";
        std::clog << "material => " << mtl << std::endl;
    }
    for (auto& v : mesh.vertices) {
        ofs << "v " << v.position[0] << " " << v.position[1] << " " << v.position[2] << "\n";
    }

    if (std::is_base_of<CVertexTexCoord, typename Mesh::CVertex>::value) {
        obj::ExportTexCoord(ofs, mesh);
    }

    for (int i = 0; i < mesh.faces.size(); i++) {
        ofs << "f ";
        for (int v = 0; v < mesh.faces[i].index.size(); v++) {
            int vid = mesh.faces[i].index[v] + 1;
            ofs << vid;
            if (std::is_base_of<CVertexTexCoord, typename Mesh::CVertex>::value) {
                ofs << "/" << vid;
            }
            ofs << " ";
        }
        ofs << "\n";
    }
}
}
}
