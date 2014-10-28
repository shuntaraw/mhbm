// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include "MeshBase.h"
#include "MeshFile.h"
#include "MeshUtil.h"
#include "MeshLaplacian.h"

namespace hbm {

/// vertex of a mesh
struct CCustomVertex :
        public slib::mesh::CVertexPosition,
        public slib::mesh::CVertexColor,
        public slib::mesh::CVertexNormal {
    float elasticity = 1; ///< elasticity
    float area = 1; ///<  voronoi area
};

/// face of a mesh
struct CCustomFace :
    public slib::mesh::CTriangleIndex {};

/// triangular mesh
class CMesh : public slib::mesh::CPolygonMesh<CCustomVertex, CCustomFace> {
public:
    CMesh() = default;
    CMesh(const CMesh&) = default;
    CMesh& operator=(const CMesh& m) = default;
    CMesh(CMesh&& m) : slib::mesh::CPolygonMesh<CCustomVertex, CCustomFace>(std::move(m)) {}
    CMesh& operator=(CMesh && m) {
        slib::mesh::CPolygonMesh<CCustomVertex, CCustomFace>::operator=(std::move(m));
        return *this;
    }

    /// update vertex property
    void UpdateGeometry() {
        if (faces.empty()) {
            return;
        }
        // normal
        slib::mesh::UpdateNormal(*this);
        // area
        std::vector<float> area;
        slib::mesh::CalculateVoronoiArea(*this, area);
        int nvertices = this->vertices.size();
        for (int i = 0; i < nvertices; i++) {
            this->vertices[i].area = area[i];
        }
    }
};

/// compute a laplacian matrix to generate differential coordinate
/// @return laplacian matrix
inline
slib::CSparseMatrix<double> GenerateLaplacianMatrix(CMesh mesh, const slib::CMatrix<double>& org_pos
                                                   ) {
    slib::mesh::ConvertCoordinateToMesh(org_pos, mesh);
    slib::CSparseMatrix<double> L;
    slib::mesh::GetMeanCurvatureNormalLaplacian(mesh, L, true);
    return L;
}

inline
void Dump(const hbm::CMesh& src) {
    static int id = 0;
    char filename[_MAX_PATH];
    sprintf(filename, "dump-%05d.ply", id++);
    slib::mesh::Write(src, filename);
}

/// affine transform 3D coordinate matrix
/**
@verbatim
[x,y,z] = [R(x,y,z)+(tx,ty,tz)]
[  :  ]   [        :          ]
@endverbatim
*/
template <typename T1, typename T2>
inline
void AffineTransformCoordinateMatrix(
    const slib::CMatrix<T1, 4, 4>& trans, // affine transformation
    slib::CMatrix<T2>& coordinates // mx3 coordinate matrix
) {
    const int nrows = coordinates.num_rows();
    for (int r = 0; r < nrows; r++) {
        auto p =  AffineTransform(trans, make_vector_from_row(coordinates, r));
        for (int c = 0; c < 3; c++) {
            coordinates(r, c) = p[c];
        }
    }
}

/// reorder 3D coordinate matrix from mx3 to 3mx1
/**
@verbatim
[x y z] -> [x]
[: : :]    [y]
[     ]    [z]
[     ]    [:]
@endverbatim
*/
/// @return 3mx1 coordinate matrix
inline
slib::CMatrix<double> VectorizeCoordinates(const slib::CMatrix<double>& mat) {
    assert(mat.num_cols() == 3);
    slib::CMatrix<double> ret(mat.num_rows() * 3, 1);
    for (int c = 0; c < 3; c++) {
        for (int r = 0; r < mat.num_rows(); r++) {
            ret(3 * r + c, 0) = mat(r, c);
        }
    }
    return ret;
}

/// reorder 3D coordinate matrix from mx3 to 3mx1
/**
@verbatim
[x] -> [x y z]
[y]    [  :  ]
[z]
[:]
@endverbatim
*/
/// @return 3mx1 coordinate matrix
inline
slib::CMatrix<double> StackCoordinates(const slib::CMatrix<double>& mat) {
    assert(mat.num_cols() == 1);
    assert(mat.num_rows() % 3 == 0);
    int nrows = mat.num_rows() / 3;
    slib::CMatrix<double> ret(nrows, 3);
    for (int c = 0; c < 3; c++) {
        for (int r = 0; r < nrows; r++) {
            ret(r, c) = mat(3 * r + c, 0);
        }
    }
    return ret;
}

/// expand a matrix from nxn to 3nx3n block diagonal
/**
@verbatim
[e  ]   [e   ]
[ \ ] = [ e  ]
[   ]   [  e ]
[   \]
@endverbatim
*/
/// @return 3nx3n block diagonal matrix
inline
slib::CSparseMatrix<double> ExpandLaplacianMatrix(const slib::CSparseMatrix<double>& mat) {
    //const int offset = mat.IsOneBased() ? 1 : 0;
    slib::MatrixGenerator<double> gen;
    for (int r = 0; r < mat.num_rows(); r++) {
        for (int idx = mat.row_index_ptr()[r] - 1; idx < mat.row_index_ptr()[r + 1] - 1; idx++) {
            int c = mat.col_index_ptr()[idx] - 1;
            double v = mat.element_ptr()[idx];
            gen.Add(3 * r + 0, 3 * c + 0, v);
            gen.Add(3 * r + 1, 3 * c + 1, v);
            gen.Add(3 * r + 2, 3 * c + 2, v);
        }
    }
    return gen.GenerateSparse(mat.num_rows() * 3, mat.num_cols() * 3);
}

}// namespace hbm
