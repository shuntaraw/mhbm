// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include "constant.h"
#include "CMesh.h"
#include "KdTree.h"
#include "Halfedge.h"

namespace hbm {

/// for registration
struct PositionPair {
    slib::CVector<float, 3> src;///< source point
    slib::CVector<float, 3> dst;///< destination point
    float weight;///< weight
};

/// vertex or point on a triangle
struct MeshCoordinate {
    union {
        size_t vid; ///< target vertex id
        size_t fid; ///< target triangle id
    };
    float t1; ///< target barycentric coordinate wrt e01
    float t2; ///< target barycentric coordinate wrt e02

    /// set the target to a vertex, not a point on a triangle
    void set_as_vertex() {
        t1 = -1;
    }

    /// check if the target type is a vertex, not a point on a triangle
    bool is_vertex() const {
        return t1 == -1;
    }

    /// interpolate coordinates
    slib::CVector<float, 3> get_position(const CMesh *mesh///< mesh
                                        ) const {
        auto& index = mesh->faces[fid].index;
        return
            (1 - t1 - t2) * mesh->vertices[index[0]].position +
            t1 * mesh->vertices[index[1]].position +
            t2 * mesh->vertices[index[2]].position;
    }

    /// interpolate normals
    slib::CVector<float, 3> get_normal(const CMesh *mesh///<mesh
                                      ) const {
        auto& index = mesh->faces[fid].index;
        return
            normalized_of((1 - t1 - t2) * mesh->vertices[index[0]].normal +
                          t1 * mesh->vertices[index[1]].normal +
                          t2 * mesh->vertices[index[2]].normal);
    }

    /// convert to coordinates
    void ToCoordinate(const hbm::CMesh *mesh,///<mesh
                      slib::CVector<float, 3>& pos///<point
                     ) const {
        if (is_vertex()) {
            pos = mesh->vertices[vid].position;
        } else {
            pos = get_position(mesh);
        }
    }

    /// reformat to normal matrix
    void ToNormalMatrix(const hbm::CMesh *mesh, ///<mesh
                        const slib::CVector<float, 3>& normal, ///<normal
                        int row, ///< row
                        slib::MatrixGenerator<double>& gen///< matrix generator
                       ) const {
        if (is_vertex()) {
            gen.Add(row, 3 * vid + 0, normal[0]);
            gen.Add(row, 3 * vid + 1, normal[1]);
            gen.Add(row, 3 * vid + 2, normal[2]);
        } else {
            auto& face = mesh->faces[fid];
            gen.Add(row, 3 * face.index[0] + 0, normal[0] * (1 - t1 - t2));
            gen.Add(row, 3 * face.index[0] + 1, normal[1] * (1 - t1 - t2));
            gen.Add(row, 3 * face.index[0] + 2, normal[2] * (1 - t1 - t2));
            gen.Add(row, 3 * face.index[1] + 0, normal[0] * t1);
            gen.Add(row, 3 * face.index[1] + 1, normal[1] * t1);
            gen.Add(row, 3 * face.index[1] + 2, normal[2] * t1);
            gen.Add(row, 3 * face.index[2] + 0, normal[0] * t2);
            gen.Add(row, 3 * face.index[2] + 1, normal[1] * t2);
            gen.Add(row, 3 * face.index[2] + 2, normal[2] * t2);
        }
    }

    /// reformat to coordinate matrix
    void ToPositionMatrix(const hbm::CMesh *mesh, ///<mesh
                          int row, ///<row
                          slib::MatrixGenerator<double>& gen///< matrix generator
                         ) const {
        if (is_vertex()) {
            gen.Add(row, vid, 1);
        } else {
            auto& face = mesh->faces[fid];
            gen.Add(row, face.index[0], 1 - t1 - t2);
            gen.Add(row, face.index[1], t1);
            gen.Add(row, face.index[2], t2);
        }
    }

};

/// correspondence represented in vertex/face index and barycentric coordinates
struct MeshCoordinatePair {
    MeshCoordinate src;///<source coordinate
    MeshCoordinate dst;///<destination coordinate
    float weight;///<weight
};

/// set of correspondences represented in sparse matrices.
/// let 's' and 'd' be nx3 coordinate matrices of source and target meshes,
/// then the correspondences between two meshes are given by 'WCs = WDs'.
struct MatrixPair {
    MatrixPair() = default;
    MatrixPair(const MatrixPair&) = delete;
    MatrixPair& operator=(const MatrixPair&) = delete;
    MatrixPair(MatrixPair&&  p)
        : C(std::move(p.C)),
          D(std::move(p.D)),
          W(std::move(p.W)) {}
    MatrixPair& operator=(MatrixPair && p) {
        C = std::move(p.C);
        D = std::move(p.D);
        W = std::move(p.W);
        return *this;
    }
    /// release internal memory
    void Clear() {
        C.Clear();
        D.Clear();
        W.Clear();
    }
    slib::CSparseMatrix<double> C; ///< correspondence of deformable_mesh
    slib::CSparseMatrix<double> D; ///< correspondence of static_mesh
    slib::CSparseMatrix<double> W; ///< weight
};

/// index-based correspondences from vertices to a mesh.
class MHBMLIB_API MeshCorrespondence {
public:
    MeshCorrespondence() = default;
    MeshCorrespondence(const MeshCorrespondence&) = delete;
    MeshCorrespondence& operator=(const MeshCorrespondence& cor) = delete;
    MeshCorrespondence(MeshCorrespondence&&  cor)
        : src_mesh_(cor.src_mesh_),
          dst_mesh_(cor.dst_mesh_),
          data_(std::move(cor.data_)) {}
    MeshCorrespondence& operator=(MeshCorrespondence && cor) {
        src_mesh_ = cor.src_mesh_;
        dst_mesh_ = cor.dst_mesh_;
        data_ = std::move(cor.data_);
        cor.src_mesh_ = 0;
        cor.dst_mesh_ = 0;
        return *this;
    }

    /// set mesh pointers
    void SetMesh(const CMesh *src , const CMesh *dst) {
        src_mesh_ = src ;
        dst_mesh_ = dst ;
        data_.clear();
    }

    /// convert to pairs of 3D coordinates
    /// 'src' and 'dst' are swapped if 'reverse_' is true.
    /// @return correspondence represented in 3D corrdinates
    std::vector<PositionPair> ToCoordinate() const;

    /// convert to sparse matrices with position constraints
    /// 'C' and 'D' are swapped if 'reverse_' is true.
    /// @return correspondence represented in sparse matrices
    MatrixPair ToPositionMatrix() const;

    /// convert to sparse matrices with normal constraints
    /// 'C' and 'D' are swapped if 'reverse_' is true.
    /// @return correspondence represented in sparse matrices
    MatrixPair ToNormalMatrix() const;

    /// add a correspondence
    void Append(const MeshCoordinatePair& p) {
        data_.push_back(p);
    }

    /// reserve memory
    void Reserve(size_t n) {
        data_.reserve(n);
    }

    /// add correspondences
    void Append(const MeshCorrespondence& p) {
        if (src_mesh_ != p.src_mesh_ || dst_mesh_ != p.dst_mesh_) {
            ThrowRuntimeError("incompatible meshes");
        }
        data_.insert(data_.end(), p.data_.begin(), p.data_.end());
    }

    /// reformat to pairs of coordinate
    const std::vector<MeshCoordinatePair>& ToMeshCoordinatePair() const {
        return data_;
    }

    /// invert correspondence
    void Invert();

    void ScaleWeight(float scale)  {
        for (auto& c : data_) {
            c.weight *= scale;
        }
    }
private:
    const CMesh *src_mesh_; ///< source mesh, from which correspondences are calculated.
    const CMesh *dst_mesh_; ///< destination mesh, to which correspondences are calculated.
    std::vector<MeshCoordinatePair> data_; ///< correspondences
};

/// find correspondence from 'src' to 'dst'
class MHBMLIB_API ClosestPointSearch {
public:
    /// set internal pointers
    void SetParameters(
        const CMesh *src ,
        const CMesh *dst ,
        float max_distance2,
        float min_cosangle,
        bool allow_border);

    /// callback called when dst mesh is updated
    void UpdateDstMesh();

    /// return correspondences
    MeshCorrespondence Find() const;

    /// search for mesh correspondences
    bool FindNearestVertex(
        const CMesh::CVertex& query,
        size_t& nearest,
        float& distance2) const;

    /// search for adjacent faces
    bool SearchAdjacentFaces(
        const CMesh::CVertex& query,
        size_t nearest_vid,
        float& distance2,
        MeshCoordinatePair& correspondence) const;

private:
    /// check if the point is at the border of a mesh
    bool IsBorder(const MeshCoordinatePair& c) const;

private:
    const CMesh *src_mesh_; ///< source mesh, from which correspondences are calculated.
    const CMesh *dst_mesh_; ///< destination mesh, to which correspondences are calculated.
    float max_distance2_; ///< distance threshold
    float min_cosangle_; ///< angle threshold
    bool allow_border_; ///< if the correspondence to mesh boundaries is allowed
    KdTree/*<CMesh::CVertex>*/ dst_kdtree_; ///< KD-tree for nearest point search. defined as mutable as it behave as a cache.
    hbm::HalfedgeMesh dst_halfedge_; ///< half edge data structure for adjacent face traversal
};

/// tags for search direction
enum class SEARCH_DIRECTION {
    FORWARD,///< forward
    BACKWARD,///< backward
    BIDIRECTIONAL,///< bidirectional
};

/// find correspondences between two meshes
class MHBMLIB_API MultidirectionalClosestPointSearch {
public:
    /// construct from a mesh.
    void SetParameters(const CMesh *src , ///< triangular mesh to be deformed
                       const CMesh *dst , ///< triangular mesh or point cloud
                       float max_distance2,///< max squared distance
                       float min_cosangle,///< min cosine of normal angle
                       bool allow_border,///< wheather correspondence to boarder is allowed
                       SEARCH_DIRECTION direction ///< direction of corresondence search
                      );

    /// find correspondences for point-plane correspondence. this may be called in many times in the ICP loop.
    MeshCorrespondence Find() const;

private:
    ClosestPointSearch forward_; ///< forward correspondence
    mutable ClosestPointSearch backward_;///< backward correspondence
    SEARCH_DIRECTION direction_;///< direction of correspondence search
};

} // hbm
