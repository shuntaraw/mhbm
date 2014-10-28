// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include "constant.h"
#include "CMesh.h"
#include "correspondence.h"

namespace hbm {

    /// deformable mesh
    class MHBMLIB_API DeformableMesh {
public:
    /// provide a range to mesh vertices which automatically updates related data structure after updating the coordinates.
    class VertexRange : public boost::iterator_range<std::vector<CMesh::CVertex>::iterator> {
    public:
        /// default constructor
        VertexRange(CMesh *m) : Base(m->vertices.begin(), m->vertices.end()), mesh_(m) {}
        /// destructor. related data structure is updated.
        ~VertexRange() {
            mesh_->UpdateGeometry();
        }
    private:
        typedef boost::iterator_range<std::vector<CMesh::CVertex>::iterator> Base;
        CMesh *mesh_;
    };

    DeformableMesh() {}

    /// set a surface model
    void LoadMesh(CMesh m);

    /// revert deformation.
    void ResetCoordinates();

    /// deform current and original meshes globally by affine transformation.
    void ApplyTransformation(const slib::CMatrix<float, 4, 4>& mat);

    /// replace vertices with those of a deformed mesh.
    template <typename T>
    void SetCoordinates(const slib::CMatrix<T>& pos) {
        mesh_.ConvertCoordinateToMesh(pos);
        mesh_.UpdateGeometry();
    }

    /// generate subdivision surface by Loop's algorithm.
    void SubdivideApproximatingSubdivision();

    /// generate subdivision surface by the Modified Butterfly algorithm.
    void SubdivideInterpolatingSubdivision();

    /// get mesh
    const CMesh& mesh() const {
        return mesh_;
    }

    const slib::CMatrix<double>& org_pos() const {
        return org_pos_;
    }

    /// get the range to mesh vertices. the internal data structure of mesh will be updated when the range is destroyed.
    VertexRange vertex_range() {
        return VertexRange(&mesh_);
    }

    /// set landmarks in the matrix format.
    void set_landmark(slib::CSparseMatrix<double> landmark) {
        landmark_ = std::move(landmark);
    }

    /// get landmarks in the matrix format.
    const slib::CSparseMatrix<double>& landmark() const {
        return landmark_;
    }

private:
    /// update landmark matrices by adding extra columns after subdivision.
    void SubdivideCorrespondence();

private:
    CMesh mesh_; // pointer to source mesh
    slib::CMatrix<double> org_pos_; // mx3 matrix of the original 3D coordinates of source mesh
    slib::CSparseMatrix<double> landmark_; // source landmark
};

enum class TRANSFORMATION {
    TRANSLATION,
    RIGID,
    SIMILARITY,
    size,
};

/// deform current mesh by non-rigid ICP method.
MHBMLIB_API slib::CMatrix<double> DeformNonrigid(
    const CMesh& src_mesh, ///< source mesh
    const CMesh& dst_mesh, ///< target mesh
    const slib::CSparseMatrix<double> *landmarkC,///< landmark of source; 'nullptr' if unused
    const slib::CSparseMatrix<double> *landmarkD,///< landmark of target; 'nullptr' if unused
    const slib::CMatrix<double>& src_org,
    float correspondence_weight, ///< weight of nearest-neighbor correspondence
    float landmark_weight, ///< weight of landmark correspondence
    //float rigidness_weight,
    bool point_plane_distance, ///<   true if point-to-plane constraint is used, false if point-to-point.
    TRANSFORMATION deformation_model, ///<   local deformation model = { translation, rigidity, similarity }
    float max_distance2, ///< distance threshold
    float min_cosangle, ///< angle threshold
    bool allow_border, ///< if the correspondence to mesh boundaries is allowed
    hbm::SEARCH_DIRECTION direction ///< if the correspondence is from static- to deformable- mesh.
);

/// deform current and original meshes globally by affine transformation.
MHBMLIB_API slib::CMatrix<float, 4, 4> EstimateAffine(
    const CMesh& src_mesh, ///< source mesh
    const CMesh& dst_mesh, ///< target mesh
    const slib::CSparseMatrix<double> *landmarkC, ///< landmark of source; 'nullptr' if unused
    const slib::CSparseMatrix<double> *landmarkD, ///< landmark of target; 'nullptr' if unused
    float landmark_weight, ///< weight of landmark correspondence
    TRANSFORMATION deformation_model, ///< global deformation model = { translation, rigid, similarity, pre-scaling }
    float max_distance2, ///< distance threshold
    float min_cosangle, ///< angle threshold
    bool allow_border, ///< if the correspondence to mesh boundaries is allowed
    SEARCH_DIRECTION direction, ///< if the correspondence is from static- to deformable- mesh.
    bool enalble_icp  ///< true if ICP
);

} // namespace hbm
