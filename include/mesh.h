// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include "MeshBase.h"
#include "mkl_csr_driver.h"

namespace hbm {

/// vertex of a mesh
struct MeshVertex :
        public slib::mesh::CVertexPosition,
        public slib::mesh::CVertexColor,
        public slib::mesh::CVertexNormal {
    float elasticity = 1; ///< elasticity
    float area = 1; ///<  voronoi area
};

/// face of a mesh
typedef slib::mesh::CTriangleIndex  MeshFace;

/// triangular mesh
class MHBMLIB_API CMesh : public slib::mesh::CPolygonMesh<MeshVertex, MeshFace> {
public:
    CMesh() = default;
    CMesh(const CMesh&) = default;
    CMesh& operator=(const CMesh& m) = default;
    CMesh(CMesh&& m) : slib::mesh::CPolygonMesh<MeshVertex, MeshFace>(std::move(m)) {}
    CMesh& operator=(CMesh && m) {
        slib::mesh::CPolygonMesh<MeshVertex, MeshFace>::operator=(std::move(m));
        return *this;
    }

    /// update vertex property
    void UpdateGeometry();

    void Dump()const;
    void Read(const std::string& filename);
    void Write(const std::string& filename)const;

    /// compute a laplacian matrix to generate differential coordinate
    /// @return laplacian matrix
    slib::CSparseMatrix<double> GenerateLaplacianMatrix(const slib::CMatrix<double>& org_pos)const;
    void ConvertMeshToCoordinate(slib::CMatrix<double>& x)const ;
    void ConvertCoordinateToMesh(const slib::CMatrix<double>& x);
    void AffineTransform(const slib::CMatrix<float, 4, 4>& mat);
    std::vector<int> SelectDuplicatedFaces()const;
    std::vector<int> SelectNonManifoldVertices()const;
    std::vector<int> SelectCollapsedFaces()const;
    std::vector<int> SelectUnusedVertices()const;
    void DeleteVertices(const std::vector<int>& vid_to_delete);
    void DeleteFaces(const std::vector<int>& fid_to_delete) ;
    CMesh SubdivideModifiedButterfly() const;
    CMesh SubdivideLoop() const;
};



}// namespace hbm
