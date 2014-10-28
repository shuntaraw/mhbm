// Copyright (c) 2012-2013 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once
#include "mkl_csr_driver.h"
#include "MatrixUtil.h"
#include "Halfedge.h"

namespace slib {
namespace mesh {

/// compute voronoi areas
/// @tparam Mesh polygon mesh
/// @tparam T area
/// @return area
template <typename Mesh, typename T>
inline
void CalculateVoronoiArea(
    const Mesh& mesh, ///< polygon mesh
    std::vector<T>& area ///< voronoi area
) {
    int nvertices = mesh.vertices.size();
    area.resize(nvertices, 0);
    int nerrors = 0;
    for (auto& face : mesh.faces) {
        int nv = face.index.size();
        int v0 = face.index[0];
        for (int fvid = 0; fvid < nv - 2; fvid++) {
            int v1 = face.index[fvid + 1];
            int v2 = face.index[fvid + 2];

            CVector<T, 3> e01 = mesh.vertices[v1].position - mesh.vertices[v0].position;
            CVector<T, 3> e12 = mesh.vertices[v2].position - mesh.vertices[v1].position;
            CVector<T, 3> e20 = mesh.vertices[v0].position - mesh.vertices[v2].position;

            T cos0 = -dot(e01, e20);
            T sin0 = norm2_of(cross(e01, e20));
            T cos1 = -dot(e01, e12);
            T sin1 = norm2_of(cross(e01, e12));
            T cos2 = -dot(e12, e20);
            T sin2 = norm2_of(cross(e12, e20));
            if (sin0 == 0 || sin1 == 0 || sin2 == 0) {
                nerrors++;
                continue;
            }
            T cot0 = cos0 / sin0;
            T cot1 = cos1 / sin1;
            T cot2 = cos2 / sin2;

            T triangle = sin0 / 2;
            T a0, a1, a2; // Voronoi area
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
                a0 = (cot1 * dot(e20, e20) + cot2 * dot(e01, e01)) / 8;
                a1 = (cot2 * dot(e01, e01) + cot0 * dot(e12, e12)) / 8;
                a2 = (cot0 * dot(e12, e12) + cot1 * dot(e20, e20)) / 8;
            }

            area[v0] += a0;
            area[v1] += a1;
            area[v2] += a2;
        }
    }
    if (nerrors) {
        std::cerr <<  "warning: " << nerrors << " collapsed triangle(s)" << std::endl;
    }
}

/// construct a laplacian matrix of mean curvature normal operator (aka cotangent weight)
/// @tparam Mesh polygon mesh
/// @tparam T type of each differencial coordinate
/// @return differential coordinate operator
/// @see http://multires.caltech.edu/pubs/diffGeoOps.pdf
template <typename Mesh, typename T>
inline
void GetMeanCurvatureNormalLaplacian(
    const Mesh& mesh, ///< polygon mesh
    CSparseMatrix<T>& L,
    bool divide_by_area ///< normalize differential coordinaes by voronoi area
) {
    int nerrors = 0;
    int nvertices = mesh.vertices.size();
    MatrixGenerator<T> genL;
    std::vector<T> area(nvertices, 0);
    for (auto& face : mesh.faces) {
        if (face.index.size() != 3) {
            ThrowRuntimeError("not triangular mesh");
        }

        int v0 = face.index[0];
        int v1 = face.index[1];
        int v2 = face.index[2];

        auto e01 = mesh.vertices[v1].position - mesh.vertices[v0].position;
        auto e12 = mesh.vertices[v2].position - mesh.vertices[v1].position;
        auto e20 = mesh.vertices[v0].position - mesh.vertices[v2].position;

        T cos0 = -dot(e01, e20);
        T sin0 = norm2_of(cross(e01, e20));
        T cos1 = -dot(e01, e12);
        T sin1 = norm2_of(cross(e01, e12));
        T cos2 = -dot(e12, e20);
        T sin2 = norm2_of(cross(e12, e20));

        if (sin0 == 0 || sin1 == 0 || sin2 == 0) {
            nerrors++;
            continue;
        }

        T cot0 = cos0 / sin0;
        T cot1 = cos1 / sin1;
        T cot2 = cos2 / sin2;

        genL.Add(v0, v0, (cot1 + cot2) / 2);
        genL.Add(v0, v1, -cot2 / 2);
        genL.Add(v0, v2, -cot1 / 2);

        genL.Add(v1, v1, (cot2 + cot0) / 2);
        genL.Add(v1, v2, -cot0 / 2);
        genL.Add(v1, v0, -cot2 / 2);

        genL.Add(v2, v2, (cot0 + cot1) / 2);
        genL.Add(v2, v0, -cot1 / 2);
        genL.Add(v2, v1, -cot0 / 2);

        if (divide_by_area) {
            T triangle = sin0 / 2;
            T a0, a1, a2; // Voronoi area
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
                a0 = (cot1 * dot(e20, e20) + cot2 * dot(e01, e01)) / 8;
                a1 = (cot2 * dot(e01, e01) + cot0 * dot(e12, e12)) / 8;
                a2 = (cot0 * dot(e12, e12) + cot1 * dot(e20, e20)) / 8;
            }
            area[v0] += a0;
            area[v1] += a1;
            area[v2] += a2;
        }
    }
    if (nerrors) {
        std::cerr <<  "warning: " << nerrors << " collapsed triangle(s)" << std::endl;
    }

    if (divide_by_area) {
        MatrixGenerator<T> genS;
        for (int r = 0; r < nvertices; r++) {
            assert(area[r]);
            genS.Add(r, r, 1 / area[r]);
        }
        L = genS.GenerateSparse(nvertices, nvertices).MultiplyTo('N', genL.GenerateSparse(nvertices, nvertices));
    } else {
        L = genL.GenerateSparse(nvertices, nvertices);
    }
}

/// construct a laplacian matrix of mean value coordinate (asymmetric).
/// @see http://cs.brown.edu/courses/cs224/papers/mean_value.pdf
/// @tparam Mesh polygon mesh
/// @tparam T type of each differencial coordinate
/// @return differential coordinate operator
template <typename Mesh , typename T>
inline
void GetMeanValueCoordinateLaplacian(const Mesh& mesh, CSparseMatrix<T>& L) {
    int nerrors = 0;
    const int nvertices = mesh.vertices.size();
    const int nfaces = mesh.faces.size();
    MatrixGenerator<T> genL;
    for (auto& face : mesh.faces) {
        if (face.index.size() != 3) {
            ThrowRuntimeError("not triangular mesh");
        }
        int v0 = face.index[0];
        const int v1 = face.index[ 1];
        const int v2 = face.index[ 2];

        auto& p0 = mesh.vertices[v0].position;
        auto& p1 = mesh.vertices[v1].position;
        auto& p2 = mesh.vertices[v2].position;

        auto n01 = normalized_of(p1 - p0);
        auto n12 = normalized_of(p2 - p1);
        auto n20 = normalized_of(p0 - p2);

        T l01 = 1 / norm2_of(p1 - p0);
        T l12 = 1 / norm2_of(p2 - p1);
        T l20 = 1 / norm2_of(p0 - p2);

        T cos0 = -dot(n20, n01);
        T cos1 = -dot(n01, n12);
        T cos2 = -dot(n12, n20);

        if (cos0 == -1 || cos1 == -1 || cos2 == -1) {
            ThrowRuntimeError("collapsed triangle");
        }

        T tan0 = sqrt((1 - cos0) / (1 + cos0));
        T tan1 = sqrt((1 - cos1) / (1 + cos1));
        T tan2 = sqrt((1 - cos2) / (1 + cos2));

        genL.Add(v0, v0, tan0 * (l01 + l20));
        genL.Add(v0, v1, -tan0 * l01);
        genL.Add(v0, v2, -tan0 * l20);

        genL.Add(v1, v1, tan1 * (l12 + l01));
        genL.Add(v1, v2, -tan1 * l12);
        genL.Add(v1, v0, -tan1 * l01);

        genL.Add(v2, v2, tan2 * (l20 + l12));
        genL.Add(v2, v0, -tan2 * l20);
        genL.Add(v2, v1, -tan2 * l12);
    }

    L = genL.GenerateSparse(nvertices, nvertices);
    MatrixGenerator<T> genS;
    for (int i = 0; i < nvertices; i++) {
        genS.Add(i, i, 1 / L(i, i));
    }
    L = genS.GenerateSparse(nvertices, nvertices).MultiplyTo('N', L);
}

/// construct a laplacian matrix of topological laplacian (area weighted)
/// @tparam Mesh polygon mesh
/// @tparam T type of each differencial coordinate
/// @return differential coordinate operator
template <typename Mesh, typename T>
inline
void GetAdjacencyMatrix(const Mesh& mesh, CSparseMatrix<T>& A) {
    std::vector<std::tuple<int, int>> edges;
    for (auto& face : mesh.faces) {
        for (int fvid = 0; fvid < face.index.size(); fvid++) {
            int v0 = face.index[fvid];
            int v1 = face.index[(fvid + 1) % face.index.size()];
            if (v0 > v1) {
                std::swap(v0, v1);
            }
            edges.push_back(std::make_tuple(v0, v1));
        }
    }
    std::sort(edges.begin(), edges.end());
    auto end = std::unique(edges.begin(), edges.end());
    int nvertices = mesh.vertices.size();
    MatrixGenerator<T> genA;
    for (auto edge = edges.begin(); edge != end; ++edge) {
        int v0 = std::get<0>(*edge);
        int v1 = std::get<1>(*edge);
        genA.Add(v0, v1, 1);
        genA.Add(v1, v0, 1);
    }
    A = genA.GenerateSparse(nvertices, nvertices);
}

/// construct a laplacian matrix of topological laplacian (area weighted)
/// @tparam Mesh polygon mesh
/// @tparam T type of each differencial coordinate
template <typename Mesh, typename T>
inline
void GetTopologicalGradient(
    const Mesh& mesh,///< polygon mesh
    CSparseMatrix<T>& L
) {
    MatrixGenerator<T> gen;
    int nedges = 0;
    for (auto& f : mesh.faces) {
        int nv = f.index.size();
        for (int v = 0; v < nv; v++) {
            gen.Add(nedges, f.index[v], -1);
            gen.Add(nedges, f.index[(v + 1) % nv], 1);
            nedges++;
        }
    }
    L = gen.GenerateSparse(nedges, mesh.vertices.size()/*,false,true*/);
}

/// construct a laplacian matrix of topological laplacian (area weighted)
/// @tparam Mesh polygon mesh
/// @tparam T type of each differencial coordinate
template <typename Mesh, typename T>
inline
void GetTopologicalLaplacian(
    const Mesh& mesh,///< polygon mesh
    CSparseMatrix<T>& L
) {
    CSparseMatrix<T> A;
    GetAdjacencyMatrix(mesh, A);
    slib::CVector<T> degrees(mesh.vertices.size(), 1);
    degrees.fill_with(1);
    degrees = A.MultiplyTo('N', degrees);
    auto D = slib::make_sparse_diagonal_matrix(degrees);
    L = D.AddTo('N', -1, A);
}

}
}
