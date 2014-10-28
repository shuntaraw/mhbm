// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#include "constant.h"

#include "MeshFile.h"
#include "MeshUtil.h"
#include "MeshLaplacian.h"
#include "mesh.h"
#include "MeshSubdivision.h"

namespace hbm {


/// update vertex property
void CMesh::UpdateGeometry() {
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

void CMesh::Dump() const {
    static int id = 0;
    char filename[_MAX_PATH];
    sprintf(filename, "dump-%05d.ply", id++);
    Write(filename);
}

void CMesh::Read(const std::string& filename) {
    slib::mesh::Read(*this, filename);
}

void CMesh::Write(const std::string& filename)const {
    slib::mesh::Write(*this, filename);
}

/// compute a laplacian matrix to generate differential coordinate
/// @return laplacian matrix
slib::CSparseMatrix<double> CMesh::GenerateLaplacianMatrix(const slib::CMatrix<double>& org_pos)const {
    CMesh work(*this);
    work.ConvertCoordinateToMesh(org_pos );
    slib::CSparseMatrix<double> L;
    slib::mesh::GetMeanCurvatureNormalLaplacian(work, L ,false );
//    slib::mesh::GetMeanValueCoordinateLaplacian(work, L );
    return L;
}

void CMesh::ConvertMeshToCoordinate(slib::CMatrix<double>& x) const {
    slib::mesh::ConvertMeshToCoordinate(*this, x);
}

void CMesh::ConvertCoordinateToMesh(const slib::CMatrix<double>& x) {
    slib::mesh::ConvertCoordinateToMesh(x, *this);
}
void CMesh::AffineTransform(const slib::CMatrix<float, 4, 4>& mat) {
    slib::mesh::AffineTransform(*this, mat);
}

std::vector<int> CMesh::SelectDuplicatedFaces()const {
    return slib::mesh::SelectDuplicatedFaces(*this);
}

std::vector<int> CMesh::SelectNonManifoldVertices()const {
    return slib::mesh::SelectNonManifoldVertices(*this);
}

std::vector<int> CMesh::SelectCollapsedFaces()const {
    return slib::mesh::SelectCollapsedFaces(*this);
}

std::vector<int> CMesh::SelectUnusedVertices()const {
    return slib::mesh::SelectUnusedVertices(*this);
}

void CMesh::DeleteVertices(const std::vector<int>& vid_to_delete) {
    slib::mesh::DeleteVertices(*this, vid_to_delete);
}

void CMesh::DeleteFaces(const std::vector<int>& fid_to_delete) {
    slib::mesh::DeleteFaces(*this, fid_to_delete);
}

CMesh CMesh::SubdivideModifiedButterfly() const {
    return slib::mesh::SubdivideModifiedButterfly(*this);
}

CMesh CMesh::SubdivideLoop() const {
    return slib::mesh::SubdivideLoop(*this);
}

}// namespace hbm
