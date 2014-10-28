// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#include "constant.h"

#include "correspondence.h"

template <typename T>
T clamp(T value, T low, T high) {
    return (value < low) ? low : (value > high ? high : value);
}

/// calculate Tukey's biweight function.
/// @see http://research.microsoft.com/en-us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html
/// @return weight
float CalculateWeight(float distance2, float max_distance2) {
    assert(distance2 >= 0);
    assert(distance2 <= max_distance2);
    return 1 - distance2 / max_distance2;
}

/// find the closest point to a triangle in 3D.
///
/// @see http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
///  ^t
/// p2
///  |\ [query]
///  | \   
/// p0--p1>s
void FindNearestPointToTriangle(const slib::CVector<float, 3>& query, // 3D coordinate of a query point
                                const slib::CVector<float, 3>& p0, // 3D coordinate of vertex 0
                                const slib::CVector<float, 3>& p1, // 3D coordinate of vertex 1
                                const slib::CVector<float, 3>& p2, // 3D coordinate of vertex 2
                                float& s, // barycentric coordinate on edge 01
                                float& t, // barycentric coordinate on edge 02
                                float& distance2 // squared distance to the nearest point
                               ) {
    auto BP = p0 - query;
    auto E0 = p1 - p0;
    auto E1 = p2 - p0;

    float a = dot(E0, E0);
    float b = dot(E0, E1);
    float c = dot(E1, E1);
    float d = dot(E0, BP);
    float e = dot(E1, BP);
    float f = dot(BP, BP);

    float det = a * c - b * b;
    s = b * e - c * d;
    t = b * d - a * e;

    if (det == 0) { // degenerated
        if (a == 0) { // p01---p2 (|E0|=0)
            s = 0;
            // Q(t)=ctt+2et+f
            // Q'(t)=0 when t=-e/c
            if (c == 0) {
                t = 0;
            } else {
                t = clamp(-e / c, 0.f, 1.f);
            }
        } else if (c == 0) { // p02---p1 (|E1|=0)
            t = 0;
            // Q(s,t)=ass+2ds+f
            // Q'(s,t)=0 when s=-d/a    (a!=0)
            s = clamp(-d / a, 0.f, 1.f);
        } else if (b > 0) {
            if (a > c) { // p0---p2---p1
                t = 0;
                // Q(s)=ass+2ds+f
                // Q'(s)/2=0 when s=-d/a    (a!=0)
                s = clamp(-d / a, 0.f, 1.f);
            } else { // p0---p1---p2
                s = 0;
                // Q(s,t)=ctt+2et+f
                // Q'(s,t)=0 when t=-e/c    (c!=0)
                t = clamp(-e / c, 0.f, 1.f);
            }
        } else { // p1---p0---p2
            // Q(s)=ass+2b(1-s)s+c(1-s)(1-s)+2ds+2e(1-s)+f
            // Q'(s)/2=as+b-2bs-c(1-s)+d-e
            //        =(a-2b+c)s+(b-c+d-e)
            //        =0 when s=-(b-c+d-e)/(a-2b+c)    (a-2b+c!=0)
            s = clamp((c + e - b - d) / (a + c - 2 * b), 0.f, 1.f);
            t = 1 - s;
        }
    } else {
        // 2 ^t
        //  \|
        //  p2
        // 3 |\ 1
        //   |0\
        // -p0--p1->s
        // 4 | 5 \ 6
        int region;
        if (s + t <= det) {
            if (s < 0) {
                if (t < 0) {
                    region = 4;
                } else {
                    region = 3;
                }
            } else if (t < 0) {
                region = 5;
            } else {
                region = 0;
            }
        } else {
            if (s < 0) {
                region = 2;
            } else if (t < 0) {
                region = 6;
            } else {
                region = 1;
            }
        }

        switch (region) {
        case 0: { // inside
            assert(det);
            float invdet = 1 / det;
            s *= invdet;
            t *= invdet;
            break;
        }
        case 1: { // on s+t=1
            // F(s) = Q(s,1-s) = (a-2b+c)s^2 + 2(b-c+d-e)s + (c+2e+f)
            // F'(s)/2 = (a-2b+c)s + (b-c+d-e)
            // F'(S) = 0 when S = (c+e-b-d)/(a-2b+c)
            // a-2b+c = |E0-E1|^2 > 0, so only sign of c+e-b-d need be considered
            float numer = c + e - b - d;
            if (numer <= 0) {
                s = 0;
            } else {
                float denom = a - 2 * b + c; // positive quantity
                s = (numer >= denom ? 1 : numer / denom);
            }
            t = 1 - s;
            break;
        }
        case 2:
            // Grad(Q) = 2(as+bt+d,bs+ct+e)
            // (0,-1)*Grad(Q(0,1)) = (0,-1)*(b+d,c+e) = -(c+e)
            // (1,-1)*Grad(Q(0,1)) = (1,-1)*(b+d,c+e) = (b+d)-(c+e)
            // min on edge s+t=1 if (1,-1)*Grad(Q(0,1)) < 0)
            // min on edge s=0 otherwise
            if (c + e > b + d) { // minimum on edge s+t=1
                float numer = c + e - b - d; // > 0
                float denom = a - 2 * b + c; // > 0
                s = (numer >= denom ? 1 : numer / denom);
                t = 1 - s;
            } else { // minimum on edge s=0
                // c>=0
                s = 0;
                t = ((c + e) <= 0 ? 1 : (e >= 0 ? 0 : -e / c));
            }
            break;
        case 3: // on s=0
            // F(t) = Q(0,t) = ct^2 + 2et + f
            // F'(t)/2 = ct+e
            // F'(T) = 0 when T = -e/c
            s = 0;
            // c>=0
            t = (e >= 0 ? 0 : (-e >= c ? 1 : -e / c));
            break;
        case 4:
            // Grad(Q) = 2(as+bt+d,bs+ct+e)
            // (1,0)*Grad(Q(0,0)) = (1,0)*(d,e) = d
            // (0,1)*Grad(Q(0,0)) = (0,1)*(d,e) = e
            // min on edge s=0 if (1,0)*Grad(Q(0,0)) < 0)
            // min on edge t=0 otherwise
            if (d < 0) { // minimum on edge s=0
                s = 0;
                // c>=0
                t = (e >= 0 ? 0 : (-e >= c ? 1 : -e / c));
            } else { // minimum on edge t=0
                // a>=0
                s = (d >= 0 ? 0 : (-d >= a ? 1 : -d / a));
                t = 0;
            }
            break;
        case 5: // on t=0
            // F(s) = Q(s,0) = as^2 + 2ds + f
            // F'(s)/2 = at+d
            // F'(s) = 0 when s = -d/a
            s = (d >= 0 ? 0 : (-d >= a ? 1 : -d / a));
            t = 0;
            break;
        case 6:
            // Grad(Q) = 2(as+bt+d,bs+ct+e)
            // (-1,0)*Grad(Q(1,0)) = (-1,0)*(a+d,b+e) = -(a+d)
            // (-1,1)*Grad(Q(1,0)) = (-1,1)*(a+d,b+e) = (b+e)-(a+d)
            // min on edge s+t=1 if (-1,1)*Grad(Q(1,0)) < 0)
            // min on edge t=0 otherwise
            if ((a + d) > (b + e)) { // minimum on edge s+t=1
                // a-2b+c = |E0-E1|^2 > 0, so only sign of c+e-b-d need be considered
                float numer = c + e - b - d;
                if (numer < 0) {
                    s = 0;
                } else {
                    float denom = a - 2 * b + c;
                    s = (numer >= denom ? 1 : numer / denom);
                }
                t = 1 - s;
            } else { // minimum on edge t=0
                s = 0;
                // c>=0
                t = ((c + e) <= 0 ? 1 : (e >= 0 ? 0 : -e / c));
            }
            break;
        default:
            ThrowLogicError("undefined label");
        }
    }
    distance2 = std::max(0.0f, a * s * s + 2 * b * s * t + c * t * t + 2 * d * s + 2 * e * t + f);
}

namespace hbm {

slib::CVector<float, 3> MeshCoordinate::get_position(const CMesh *mesh///< mesh
                                                    ) const {
    auto& index = mesh->faces[fid].index;
    return
        (1 - t1 - t2) * mesh->vertices[index[0]].position +
        t1 * mesh->vertices[index[1]].position +
        t2 * mesh->vertices[index[2]].position;
}

slib::CVector<float, 3> MeshCoordinate::get_normal(const CMesh *mesh///<mesh
                                                  ) const {
    auto& index = mesh->faces[fid].index;
    return
        normalized_of((1 - t1 - t2) * mesh->vertices[index[0]].normal +
                      t1 * mesh->vertices[index[1]].normal +
                      t2 * mesh->vertices[index[2]].normal);
}

void MeshCoordinate::ToCoordinate(const hbm::CMesh *mesh,///<mesh
                                  slib::CVector<float, 3>& pos///<point
                                 ) const {
    if (is_vertex()) {
        pos = mesh->vertices[vid].position;
    } else {
        pos = get_position(mesh);
    }
}

void MeshCoordinate::ToNormalMatrix(const hbm::CMesh *mesh, ///<mesh
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

void MeshCoordinate::ToPositionMatrix(const hbm::CMesh *mesh, ///<mesh
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

void ClosestPointSearch::SetParameters(const CMesh *src ,
                                       const CMesh *dst ,
                                       float max_distance2,
                                       float min_cosangle,
                                       bool allow_border) {
    src_mesh_ = src ;
    dst_mesh_ = dst ;
    max_distance2_ = max_distance2;
    min_cosangle_ = min_cosangle;
    allow_border_ = allow_border;

    dst_kdtree_.Construct(
        dst_mesh_->vertices.begin(),
        dst_mesh_->vertices.end()); // initial update
    dst_halfedge_.Construct(*dst_mesh_);
}

void ClosestPointSearch::UpdateDstMesh() {
    dst_kdtree_.Reconstruct();
}

MeshCorrespondence ClosestPointSearch::Find() const {
    MeshCorrespondence correspondence;
    correspondence.SetMesh(src_mesh_, dst_mesh_);

    int nvertices = src_mesh_->vertices.size();
    correspondence.Reserve(nvertices);
    #pragma omp parallel for
    for (int src_vid = 0; src_vid < nvertices; src_vid++) {
        // find the nearest vertex
        MeshCoordinatePair c;
        float distance2;
        if (!FindNearestVertex(src_mesh_->vertices[src_vid], c.dst.vid, distance2)) {
            // not found
            continue;
        }

        // find the nearest point on triangles
        if (dst_mesh_->faces.size()) {
            if (!SearchAdjacentFaces(src_mesh_->vertices[src_vid], c.dst.vid,  distance2, c)) {
                // isolated point in meshes
                continue;
            }
            if (!allow_border_ && IsBorder(c)) {
                // on border
                continue;
            }
        } else {
            c.dst.set_as_vertex();
        }

        c.src.vid = src_vid;
        c.src.set_as_vertex();
        c.weight = CalculateWeight(distance2, max_distance2_) * src_mesh_->vertices[src_vid].area;

        #pragma omp critical
        {
            correspondence.Append(c);
        }
    }
    return correspondence;
}

bool ClosestPointSearch::FindNearestVertex(const CMesh::CVertex& query,
        size_t& nearest_vid,
        float& distance2) const {
    if (!dst_kdtree_.GetKNearest(query.position.data(), max_distance2_,  1, &nearest_vid, &distance2) ||
        dot(query.normal, dst_mesh_->vertices[nearest_vid].normal) < min_cosangle_) {
        return false;
    }
    return true;
}

bool ClosestPointSearch::SearchAdjacentFaces(const CMesh::CVertex& query,
        size_t nearest_vid,
        float& distance2,
        MeshCoordinatePair& correspondence
                                            ) const {
    auto first = dst_halfedge_.vertex_edge(nearest_vid);
    if (!first) {
        // isolated point
        return false;
    }

    // nearest vertex
    correspondence.dst.fid = first->face();
    if (nearest_vid == dst_mesh_->faces[first->face()].index[0]) {
        correspondence.dst.t1 = 0;
        correspondence.dst.t2 = 0;
    } else if (nearest_vid == dst_mesh_->faces[first->face()].index[1]) {
        correspondence.dst.t1 = 1;
        correspondence.dst.t2 = 0;
    } else {
        correspondence.dst.t1 = 0;
        correspondence.dst.t2 = 1;
    }

    auto edge = first;
    do {
        auto& face = dst_mesh_->faces[edge->face()];

        //// check face normal
        //if (dot(face.normal, query.normal) < min_cosangle_) {
        //    edge = edge->next()->pair();
        //    continue;
        //}

        // search for the nearest
        float s, t, d2;
        FindNearestPointToTriangle(
            query.position,
            dst_mesh_->vertices[face.index[0]].position,
            dst_mesh_->vertices[face.index[1]].position,
            dst_mesh_->vertices[face.index[2]].position,
            s,
            t,
            d2);
        if (distance2 > d2) {
            // interpolated normal on the face
            MeshCoordinate dst;
            dst.fid = edge->face();
            dst.t1 = s;
            dst.t2 = t;
            auto normal = dst.get_normal(dst_mesh_);
            //normalized_of(
            //(1-s-t)* dst_mesh_->vertices[face.index[0]].normal+
            //s*dst_mesh_->vertices[face.index[1]].normal+
            //t*dst_mesh_->vertices[face.index[2]].normal
            //);
            // check face normal
            if (dot(normal, query.normal) >= min_cosangle_) {
                // update nearest
                correspondence.dst = dst;
                //correspondence.dst.fid = edge->face();
                //correspondence.dst.t1 = s;
                //correspondence.dst.t2 = t;
                distance2 = d2;
            }
        }
        edge = edge->next()->pair();
    } while (edge && edge != first);
    return true;
}

bool ClosestPointSearch::IsBorder(const MeshCoordinatePair& correspondence) const {
    auto halfedge = dst_halfedge_.face_edge(correspondence.dst.fid);
    auto& face = dst_mesh_->faces[correspondence.dst.fid];

    // nearest on e20
    if (correspondence.dst.t1 == 0) {
        while (halfedge->vertex() != face.index[0]) {
            halfedge = halfedge->next();
        }
        if (!halfedge->pair()) {
            return true;
        }
    }

    // nearest is on e01
    if (correspondence.dst.t2 == 0) {
        while (halfedge->vertex() != face.index[1]) {
            halfedge = halfedge->next();
        }
        if (!halfedge->pair()) {
            return true;
        }
    }

    // nearest is on e12
    if (correspondence.dst.t1 + correspondence.dst.t2 == 1) {
        while (halfedge->vertex() != face.index[2]) {
            halfedge = halfedge->next();
        }
        if (!halfedge->pair()) {
            return true;
        }
    }
    return false; // not on border
}

void MultidirectionalClosestPointSearch::SetParameters(const CMesh *src , // triangular mesh
        const CMesh *dst , // triangular mesh or point cloud
        float max_distance2,
        float min_cosangle,
        bool allow_border,
        SEARCH_DIRECTION direction
                                                      ) {
    direction_ = direction;
    switch (direction) {
    case SEARCH_DIRECTION::FORWARD:
        forward_.SetParameters(src , dst , max_distance2, min_cosangle, allow_border);
        break;
    case SEARCH_DIRECTION::BACKWARD:
        backward_.SetParameters(dst , src , max_distance2, min_cosangle, allow_border);
        break;
    case SEARCH_DIRECTION::BIDIRECTIONAL:
        forward_.SetParameters(src , dst , max_distance2, min_cosangle, allow_border);
        backward_.SetParameters(dst , src , max_distance2, min_cosangle, allow_border);
        break;
    default:
        ThrowLogicError("undefined label");
    }
}

MeshCorrespondence MultidirectionalClosestPointSearch::Find() const {
    switch (direction_) {
    case SEARCH_DIRECTION::FORWARD:
        return forward_.Find();
    case SEARCH_DIRECTION::BACKWARD: {
        backward_.UpdateDstMesh(); // mutable
        auto bwd = backward_.Find();
        bwd.Invert();
        return bwd;
    }
    case SEARCH_DIRECTION::BIDIRECTIONAL: {
        auto fwd = forward_.Find();
        backward_.UpdateDstMesh(); // mutable
        auto bwd = backward_.Find();
        bwd.Invert();
        fwd.Append(bwd);
        fwd.ScaleWeight(0.5);
        return fwd;
    }
    default:
        ThrowLogicError("undefined label");
    }
}

void MeshCorrespondence::Append(const MeshCorrespondence& p) {
    if (src_mesh_ != p.src_mesh_ || dst_mesh_ != p.dst_mesh_) {
        ThrowRuntimeError("incompatible meshes");
    }
    data_.insert(data_.end(), p.data_.begin(), p.data_.end());
}

void MeshCorrespondence::Invert() {
    std::swap(src_mesh_, dst_mesh_);
    for (auto& c : data_) {
        std::swap(c.src, c.dst);
    }
}

void MeshCorrespondence::ScaleWeight(float scale)  {
    for (auto& c : data_) {
        c.weight *= scale;
    }
}

std::vector<PositionPair> MeshCorrespondence::ToCoordinate() const {
    auto size = data_.size();
    std::vector<PositionPair> out(size);
    for (int i = 0; i < size; i++) {
        auto& c = data_[i];
        c.src.ToCoordinate(src_mesh_, out[i].src);
        c.dst.ToCoordinate(dst_mesh_, out[i].dst);
        out[i].weight = c.weight;
    }
    return out;
}

MatrixPair MeshCorrespondence::ToPositionMatrix() const {
    slib::MatrixGenerator<double> genC, genD, genW;
    auto size = data_.size();
    for (int i = 0; i < size; i++) {
        auto& c = data_[i];
        c.src.ToPositionMatrix(src_mesh_, i, genC);
        c.dst.ToPositionMatrix(dst_mesh_, i, genD);
        genW.Add(i, i, sqrt(c.weight));
    }
    MatrixPair ret;
    ret.C = genC.GenerateSparse(size, src_mesh_->vertices.size());
    ret.D = genD.GenerateSparse(size, dst_mesh_->vertices.size());
    ret.W = genW.GenerateSparse(size, size);
    return ret;
}

MatrixPair MeshCorrespondence::ToNormalMatrix() const {
    slib::MatrixGenerator<double> genC, genD, genW;
    auto size = data_.size();
    for (int i = 0; i < size; i++) {
        auto& c = data_[i];
        slib::CVector<float, 3> normal;
        if (c.src.is_vertex()) {
            if (c.dst.is_vertex()) {
                normal = dst_mesh_->vertices[c.dst.vid].normal;
            } else {
                auto& index = dst_mesh_->faces[c.dst.fid].index;
                normal = c.dst.get_normal(dst_mesh_);
            }
        } else {
            auto& index = src_mesh_->faces[c.src.fid].index;
            normal = c.src.get_normal(src_mesh_);
        }
        c.src.ToNormalMatrix(src_mesh_, normal, i, genC);
        c.dst.ToNormalMatrix(dst_mesh_, normal, i, genD);
        genW.Add(i, i, sqrt(c.weight));
    }
    MatrixPair ret;
    ret.C = genC.GenerateSparse(size, 3 * src_mesh_->vertices.size());
    ret.D = genD.GenerateSparse(size, 3 * dst_mesh_->vertices.size());
    ret.W = genW.GenerateSparse(size, size);
    return ret;
}

} // namespace hbm
