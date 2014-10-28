// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#include "precompile.h"

#include "constant.h"
#include "correspondence.h"
#include "deformation.h"
#include "landmark.h"

template <typename T>
T clamp(T value, T low, T high) {
    return (value < low) ? low : (value > high ? high : value);
}

namespace hbm {


/// compute rotation matrix from a covariance
/// @return maximized rRr
double SolveRotationSVD(
    const slib::CMatrix<double, 3, 3>& covar,///< sum_rTr
    slib::CMatrix<double, 3, 3>& rot ///< rotation
) {
    slib::CMatrix<double, 3, 3> U, Vt;
    slib::CVector<double, 3> s;
    const char jobu = 'A';
    const char jobvt = 'A';
    const int m = 3;
    const int n = 3;
    const int lda = 3;
    const int ldu = 3;
    const int ldvt = 3;
    double work[15];
    const int lwork = 15;
    int info;
    double mat[9];
    std::copy_n(covar.data(), 9, mat);
    DGESVD(&jobu, &jobvt, &m, &n, mat, &lda, s.data(), U.data(), &ldu, Vt.data(), &ldvt, work, &lwork, &info);
    rot = transpose_of(U * Vt);
    double det = determinant_of(rot);
    if (det < 0) {
        rot = transpose_of(U * slib::make_diagonal_matrix<double>(1, 1, -1) * Vt);
    }
    return s[0] + s[1] + s[2];
}


/// solve similarity from point correspondences.
/// dst = s.R.src + t:
/// @return squared residual |s.R.src + t - dst|^2
/// @see http://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf
double SolveSimilarity(
    const std::vector<PositionPair>& correspondence,
    double *scale,///<[in,out] scale (fixed to 1 if nullptr)
    slib::CMatrix<double, 3, 3> *rot,///<[in,out] rotation (fixed to identity if nullptr)
    slib::CVector<double, 3> *trans ///<[in,out] translation (cannot be nullptr)
) {
    if (correspondence.size() < 3) {
        ThrowRuntimeError("too few correspondence");
    }
    // use double instead of float; each weight can be very small
    double sum_weight = 0;
    slib::CVector<double, 3> src_center ;  // zero vector
    slib::CVector<double, 3> dst_center ;  // zero vector
    src_center.fill_with(0);
    dst_center.fill_with(0);
    for (auto& c : correspondence) {
        src_center += c.weight * c.src;
        dst_center += c.weight * c.dst;
        sum_weight += c.weight;
    }
    if (sum_weight == 0) {
        ThrowRuntimeError("too small weights");
    }
    src_center /= sum_weight;
    dst_center /= sum_weight;
    if (rot == 0) {
        assert(scale == 0);
        *trans = dst_center - src_center;
        double residual2 = 0;
        for (auto& c : correspondence) {
            residual2 += c.weight * norm2_squared_of(c.src + *trans - c.dst);
        }
        return residual2 / sum_weight;
    }

    // use double instead of float; each weight can be very small
    slib::CMatrix<double, 3, 3> covar ;
    covar.fill_with(0);
    double sum_src2 = 0;
    double sum_dst2 = 0;
    for (auto& c : correspondence) {
        auto s = c.src - src_center;
        auto d = c.dst - dst_center;
        covar += s * transpose_of(d) * c.weight;
        sum_src2 += dot(s, s) * c.weight;
        sum_dst2 += dot(d, d) * c.weight;
    }

    slib::CMatrix<double, 3, 3> drot;
    double consistency = SolveRotationSVD(covar, drot);
    *rot = drot;
    if (scale) {
        *scale = sqrt(sum_dst2 / sum_src2);
        *trans = dst_center - *scale **rot * src_center;
        return (*scale **scale * sum_src2 + sum_dst2 - 2 * *scale * consistency) / sum_weight;
    } else {
        *trans = dst_center - *rot * src_center;
        return (sum_src2 + sum_dst2 - 2 * consistency) / sum_weight;
    }
}

/// partial matrix assignment.
/// @return m1 = s.m2 (partial)
void CopyPartialMatrix(
    slib::CMatrix<double>& m1, // matrix
    const slib::CMatrix<double>& m2, // matrix
    double scale = 1, // scale
    int row = 0, // start row
    int col = 0 // start column
) {
    assert(m2.num_rows() + row <= m1.num_rows());
    assert(m2.num_cols() + col <= m1.num_cols());
    for (int mc = 0; mc < m2.num_cols(); mc++) {
        for (int mr = 0; mr < m2.num_rows(); mr++) {
            m1(row + mr, col + mc) = scale * m2(mr, mc);
        }
    }
}

/// as-rigid-as-possible coordinates
/// @return nx3 or 3nx1 differential coordinate matrix
slib::CMatrix<double> CalculateDifferentialCoordinate(
    const slib::CMatrix<double>& src_org, // original (nx3)
    const slib::CMatrix<double>& src_pos, // current (nx3)
    const slib::CSparseMatrix<double>& src_laplacian, // laplacian (nxn)
    bool allow_scaling // scale
) {
    const int nvertices = src_org.num_rows();
    assert(nvertices == src_pos.num_rows());
    assert(nvertices == src_laplacian.num_rows());
    assert(nvertices == src_laplacian.num_cols());

    //const int offset = src_laplacian.IsOneBased() ? 1 : 0;
    std::vector<slib::CMatrix<double, 3, 3>> rotation(nvertices);

    // solve rotation (and scale) assuming known coordinates
    #pragma omp parallel for
    for (int r = 0; r < nvertices; r++) {
        // use double instead of float. weights can be very small.
        slib::CMatrix<double, 3, 3> covariance ; // zero vector
        covariance.fill_with(0);
        double orgsum2 = 0;
        double newsum2 = 0;
        for (int i = src_laplacian.row_index_ptr()[r] - 1; i < src_laplacian.row_index_ptr()[r + 1] - 1; i++) {
            int c = src_laplacian.col_index_ptr()[i] - 1;
            if (r == c) {
                // diagonal element
                continue;
            }
            // unweighted. cotangent weights can be negative, resulting in invalid results.
            auto orgedge =  make_vector_from_row(src_org, r) - make_vector_from_row(src_org, c) ;
            auto newedge =  make_vector_from_row(src_pos, r) - make_vector_from_row(src_pos, c) ;
            covariance += orgedge * transpose_of(newedge);
            if (allow_scaling) {
                orgsum2 += dot(orgedge, orgedge);
                newsum2 += dot(newedge, newedge);
            }
        }

        // solve rotation
        SolveRotationSVD(covariance, rotation[r]);

        // solve scale
        if (allow_scaling) {
            double scale = sqrt(newsum2 / orgsum2);
            scale =  clamp<double> (scale, MIN_SIMILARITY_SCALE , MAX_SIMILARITY_SCALE);
            rotation[r] *= scale;
        }
    }

    slib::CMatrix<double> LRv(nvertices, 3);

    #pragma omp parallel for
    for (int r = 0; r < nvertices; r++) {
        slib::CVector<double, 3> p;
        p.fill_with(0);
        for (int i = src_laplacian.row_index_ptr()[r] - 1; i < src_laplacian.row_index_ptr()[r + 1] - 1; i++) {
            int c = src_laplacian.col_index_ptr()[i] - 1;
            if (c == r) {
                continue; // diagonal
            }
            double w_ij = -src_laplacian.element_ptr()[i];
            double w_ji = -src_laplacian(c, r);
            slib::CMatrix<double, 3, 3> RR = w_ij * rotation[r] +  w_ji * rotation[c ];
            slib::CVector<double, 3> orgedge =  make_vector_from_row(src_org, r) -  make_vector_from_row(src_org, c) ;
            p += RR * orgedge;
        }
        LRv(r, 0) = p[0];
        LRv(r, 1) = p[1];
        LRv(r, 2) = p[2];
    }
    return LRv;
}

slib::CSparseMatrix<double> GenerateSymmetricLaplacian(const slib::CSparseMatrix<double>& L) {
    slib::MatrixGenerator<double> gen;
    int num_rows = L.num_rows();
    for (int r = 0; r < num_rows; r++) {
        for (int i = L.row_index_ptr()[r] - 1; i < L.row_index_ptr()[r + 1] - 1; i++) {
            int c = L.col_index_ptr()[i] - 1;
            if (r == c) {
                continue;
            }
            double v = L.element_ptr()[i];
            gen.Add(r, c, v);
            gen.Add(c, r, v);
            gen.Add(r, r, -v);
            gen.Add(c, c, -v);
        }
    }
    return gen.GenerateSparse(num_rows, num_rows);
}

/// deform in the local coordinate of 'src':
double SolveLaplacianDeformationConstrainedPosition(
    const slib::CSparseMatrix<double>& src_laplacian, // (nxn)
    slib::CMatrix<double>& src_pos, // [in,out] (nx3)
    const slib::CMatrix<double>& src_org, // transformed (nx3)
    const slib::CMatrix<double>& wDu, // (mx3)
    const slib::CSparseMatrix<double>& wC, // (mxn)
    TRANSFORMATION model // local deformation model = { translation, rigid, similarity }
) {
    const int nvertices = src_laplacian.num_rows();
    const int ncorrespondences = wDu.num_rows();

    assert(src_laplacian.num_rows() == nvertices);
    assert(src_laplacian.num_cols() == nvertices);
    assert(src_pos.num_rows() == nvertices);
    assert(src_pos.num_cols() == 3);
    assert(src_org.num_rows() == nvertices);
    assert(src_org.num_cols() == 3);
    assert(wDu.num_rows() == ncorrespondences);
    assert(wDu.num_cols() == 3);
    assert(wC.num_rows() == ncorrespondences);
    assert(wC.num_cols() == nvertices);


    // A = [ L]
    //     [wC]
    slib::CSparseMatrix<double> A;
    slib::PardisoDriver<double> driver;

    // b = [ * ]
    //     [wDu]
    slib::CMatrix<double> b(nvertices + ncorrespondences, 3);
    CopyPartialMatrix(b, wDu, 1, nvertices, 0);

    double residual = 0;
    switch (model) {
    case TRANSFORMATION::TRANSLATION:
        A = src_laplacian;
        A.AppendRows(1, wC);
        driver.Factorize(1, 'N', A.MultiplyTo('T', A), false); // AtA

        // b = [LRv]
        //     [wDu]
        CopyPartialMatrix(b, src_laplacian.MultiplyTo('N', src_org));

        driver.Solve(A.MultiplyTo('T', b), src_pos, false);

        // |b - Ax|^2
        b -= A.MultiplyTo('N', src_pos);
        residual = norm2_squared_of(b); // Fourobenius norm
        break;

    case TRANSFORMATION::RIGID:
    case TRANSFORMATION::SIMILARITY: {
        A = GenerateSymmetricLaplacian(src_laplacian);
        A.AppendRows(1, wC);
        driver.Factorize(1, 'N', A.MultiplyTo('T', A), false); // AtA

        double prev_residual = std::numeric_limits<float>::max();
        for (int niters = 0;  ; niters++) {
            // b = [LRv]
            //     [wDu]
            bool allow_scaling = (model == TRANSFORMATION::SIMILARITY);
            auto LRv = CalculateDifferentialCoordinate(src_org, src_pos, src_laplacian, allow_scaling);  // LRv
            CopyPartialMatrix(b, LRv);

            driver.Solve(A.MultiplyTo('T', b), src_pos, false);

            //prev_residual = residual;
            residual = norm2_squared_of(b - A.MultiplyTo('N', src_pos)); // Fourobenius norm

            std::clog << "\r|r| = " << sqrt(residual);

            if ((niters % MIN_ICP_ITERATIONS) == 0) {
                if (prev_residual / residual <= 1 + MIN_ICP_IMPROVEMENT) {
                    break;
                } else {
                    prev_residual = residual;
                }
            }
        }
        std::clog << "\r" ;
    }
    break;
    default:
        ThrowLogicError("undefined label");
    }
    return residual;
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

/// deform in the local coordinate of 'src'
double SolveLaplacianDeformationConstrainedNormal(
    const slib::CSparseMatrix<double>& src_laplacian, // nxn laplacian matrix
    slib::CMatrix<double>& src_pos, // [in,out] nx3 coordinate matrix
    const slib::CMatrix<double>& src_org, // nx3 original coordinate matrix
    const slib::CMatrix<double>& wDu, // mx1 matrix of weighted anchor points
    const slib::CSparseMatrix<double>& wC, // mx3n weighted correspondence matrix
    TRANSFORMATION model // local deformation model = { translation, rigid, similarity }
) {
    int nvertices = src_laplacian.num_rows();
    int ncorrespondences = wDu.num_rows();

    assert(src_laplacian.num_rows() == nvertices);
    assert(src_laplacian.num_cols() == nvertices);
    assert(src_pos.num_rows() == nvertices);
    assert(src_pos.num_cols() == 3);
    assert(src_org.num_rows() == nvertices);
    assert(src_org.num_cols() == 3);
    assert(wDu.num_rows() == ncorrespondences);
    assert(wDu.num_cols() == 1);
    assert(wC.num_rows() == ncorrespondences);
    assert(wC.num_cols() == 3 * nvertices);

    // A = [ L]
    //     [wC]
    slib::CSparseMatrix<double> A;
    slib::PardisoDriver<double> driver;

    // b = [ * ]
    //     [wDu]
    slib::CMatrix<double> b(3 * src_pos.num_rows() + wDu.num_rows(), 1);
    CopyPartialMatrix(b, wDu, 1, 3 * src_pos.num_rows(), 0);

    double residual = 0;
    switch (model) {
    case TRANSFORMATION::TRANSLATION:
        A = ExpandLaplacianMatrix(src_laplacian);
        A.AppendRows(1, wC);
        driver.Factorize(1, 'N', A.MultiplyTo('T', A), false);

        // b = [LRv]
        //     [wDu]
        CopyPartialMatrix(b, VectorizeCoordinates(src_laplacian.MultiplyTo('N', src_org)));

        driver.Solve(A.MultiplyTo('T', b), src_pos, false);

        b -= A.MultiplyTo('N', src_pos);
        residual = norm2_squared_of(b);

        src_pos = StackCoordinates(src_pos);
        break;

    case TRANSFORMATION::RIGID:
    case TRANSFORMATION::SIMILARITY: {
        A = ExpandLaplacianMatrix(GenerateSymmetricLaplacian(src_laplacian));
        A.AppendRows(1, wC);
        driver.Factorize(1, 'N', A.MultiplyTo('T', A), false);

        double prev_residual = std::numeric_limits<float>::max();
        for (int niters = 0;  ; niters++) {
            // b = [LRv]
            //     [wDu]
            slib::CMatrix<double> LRv = CalculateDifferentialCoordinate(
                                            src_org,
                                            src_pos,
                                            src_laplacian,
                                            model == TRANSFORMATION::SIMILARITY); // LRv
            CopyPartialMatrix(b, VectorizeCoordinates(LRv));

            driver.Solve(A.MultiplyTo('T', b), src_pos, false);

            residual = norm2_squared_of(b - A.MultiplyTo('N', src_pos));
            std::clog << "\r|r| = " << sqrt(residual);
            src_pos = StackCoordinates(src_pos);

            if ((niters % MIN_ICP_ITERATIONS) == 0) {
                if (prev_residual / residual <= 1 + MIN_ICP_IMPROVEMENT) {
                    break;
                } else {
                    prev_residual = residual;
                }
            }
        }
        std::clog << "\r";
    }
    break;
    default:
        ThrowLogicError("undefined label");
    }
    return residual;
}

void DeformableMesh::LoadMesh(CMesh m) {
    mesh_ = std::move(m);
    mesh_.UpdateGeometry();
    mesh_.ConvertMeshToCoordinate(org_pos_);
}

void DeformableMesh::ResetCoordinates() {
    int nvertices = mesh_.vertices.size();
    for (int r = 0; r < nvertices; r++) {
        mesh_.vertices[r].position  =  make_vector_from_row(org_pos_, r) ;
    }
    mesh_.UpdateGeometry();
}

void DeformableMesh::ApplyTransformation(const slib::CMatrix<float, 4, 4>& mat) {
    mesh_.AffineTransform(mat);
    mesh_.UpdateGeometry();
}

std::vector<PositionPair> ConvertLandmarkToCoordinates(const CMesh& src, const slib::CSparseMatrix<double>& C, const CMesh& dst, const slib::CSparseMatrix<double>& D, float weight) {
    slib::CMatrix<double> coordinates;
    src.ConvertMeshToCoordinate(coordinates);
    auto src_lm = C.MultiplyTo('N', coordinates);
    dst.ConvertMeshToCoordinate(coordinates);
    auto dst_lm = D.MultiplyTo('N', coordinates);
    int ncors = src_lm.num_rows();
    if (ncors != dst_lm.num_rows()) {
        ThrowRuntimeError("invalid number of landmarks");
    }
    std::vector<PositionPair> ret(ncors);
    for (int i = 0; i < ncors; i++) {
        ret[i].src =  make_vector_from_row(src_lm, i) ;
        ret[i].dst =  make_vector_from_row(dst_lm, i) ;
        ret[i].weight = weight;
    }
    return ret;
}

/// deform current and original meshes globally by affine transformation.
MHBMLIB_API slib::CMatrix<float, 4, 4> EstimateAffine(
    const CMesh& src_mesh,
    const CMesh& dst_mesh,
    const slib::CSparseMatrix<double> *landmarkC,
    const slib::CSparseMatrix<double> *landmarkD,
    float landmark_weight,
    TRANSFORMATION model,
    float max_distance2,
    float min_cosangle,
    bool allow_border,
    SEARCH_DIRECTION direction,
    bool enalble_icp
) {
    if (landmarkC) {
        assert(landmarkD);
        assert(landmarkC->num_rows() == landmarkD->num_rows());
        assert(landmarkC->num_cols() == src_mesh.vertices.size());
        assert(landmarkD->num_cols() == dst_mesh.vertices.size());
    } else {
        assert(!landmarkD);
        assert(enalble_icp);
    }

    // accumulated transformation
    float scale = 1;
    slib::CMatrix<float, 4, 4> rigid = slib::make_diagonal_matrix<float>(1, 1, 1, 1);

    // current working mesh
    CMesh cur_mesh = src_mesh;

    // nearest-neighbor search
    MultidirectionalClosestPointSearch mcp;
    if (enalble_icp) {
        mcp.SetParameters(&cur_mesh, &dst_mesh, max_distance2, min_cosangle, allow_border, direction);
    }

    double residual = 0;
    double prev_residual = std::numeric_limits<float>::max();
    for (int niters = 0;  ; niters++) {
        std::vector<PositionPair> corpos;
        if (enalble_icp) {
            corpos = mcp.Find().ToCoordinate();
        }
        if (landmarkC) {
            auto lmpos = ConvertLandmarkToCoordinates(cur_mesh, *landmarkC, dst_mesh, *landmarkD, landmark_weight);
            corpos.insert(corpos.end(), lmpos.begin(), lmpos.end());
        }
        int ncors = corpos.size();
        switch (model) {
        case TRANSFORMATION::TRANSLATION: {
            // cur = src+t
            // dst <-> cur+t'
            slib::CVector<double, 3> t;
            residual = SolveSimilarity(corpos, 0, 0, &t);
            rigid = make_translation_matrix(t) * rigid;
            cur_mesh.vertices = src_mesh.vertices;
            cur_mesh.AffineTransform(rigid);
        }
        break;
        case TRANSFORMATION::RIGID: {
            // cur = [Rt]src
            // dst <-> {Rt}cur
            slib::CMatrix<double, 3, 3> R;
            slib::CVector<double, 3> t;
            residual = SolveSimilarity(corpos, 0, &R, &t);
            rigid =  make_affine_matrix(R, t) * rigid;
            slib::NormalizeRotation(rigid);
            cur_mesh.vertices = src_mesh.vertices;
            cur_mesh.AffineTransform(rigid);
        }
        break;
        case TRANSFORMATION::SIMILARITY: {
            // cur = [Rt][s]src
            // [Rt]'cur = [s]src
            // [Rt]'dst <-> {Rt}{s}[s].src
            auto invRt = inverse_rigidity_of(rigid);
            for (auto& c : corpos) {
                c.src =  AffineTransform(invRt, c.src);
                c.dst =  AffineTransform(invRt, c.dst);
            }
            double s;
            slib::CMatrix<double, 3, 3> R;
            slib::CVector<double, 3> t;
            residual = SolveSimilarity(corpos, &s, &R, &t);
            scale *= s;
            rigid *=  make_affine_matrix(R, t);
            slib::NormalizeRotation(rigid);
            cur_mesh.vertices = src_mesh.vertices;
            cur_mesh.AffineTransform(rigid * slib::make_diagonal_matrix(scale, scale, scale, 1));
            cur_mesh.UpdateGeometry();
        }
        break;
        default:
            ThrowLogicError("undefined label");
        }

        std::clog << niters << ": #points=" << ncors << ", |r|=" << std::fixed << sqrt(residual) << std::endl;
        if ((niters % MIN_ICP_ITERATIONS) == 0) {
            if (prev_residual / residual <= 1 + MIN_ICP_IMPROVEMENT) {
                break;
            } else {
                prev_residual = residual;
            }
        }
        if (!enalble_icp) {
            break;
        }
    }

    return rigid * slib::make_diagonal_matrix(scale, scale, scale, 1);
}

void DeformableMesh::SubdivideCorrespondence() {
    int ncols = mesh().vertices.size() - landmark_.num_cols();
    landmark_.AppendColumns(ncols);
}

void DeformableMesh::SubdivideApproximatingSubdivision() {
    // subdivide original
    CMesh org_mesh = mesh();
    org_mesh.ConvertCoordinateToMesh(org_pos());
    org_mesh = org_mesh.SubdivideLoop();
    org_mesh.ConvertMeshToCoordinate(org_pos_);

    // subdivide current
    int noldvertices = mesh().vertices.size();
    mesh_ = mesh().SubdivideLoop();

    // update elasticity
    auto vertices = vertex_range();
    for (int vid = noldvertices; vid < vertices.size(); vid++) {
        vertices[vid].elasticity = 1;
    }
    for (auto& face : mesh().faces) {
        float elasticity = -1;
        for (int vid : face.index) {
            if (vid < noldvertices) {
                // old vertex, at most one per face
                elasticity = mesh().vertices[vid].elasticity;
                break;
            }
        }
        if (elasticity != -1) {
            for (int vid : face.index) {
                if (mesh().vertices[vid].elasticity > elasticity) {
                    vertices[vid].elasticity = elasticity;    // set to the lowest
                }
            }
        }
    }
    if (landmark_.num_rows()) {
        SubdivideCorrespondence();
    }
}

void DeformableMesh::SubdivideInterpolatingSubdivision() {
    // subdivide original
    CMesh org_mesh = mesh();
    org_mesh.ConvertCoordinateToMesh(org_pos());
    org_mesh = org_mesh.SubdivideModifiedButterfly();
    org_mesh.ConvertMeshToCoordinate(org_pos_);

    // subdivide current
    int noldvertices = mesh().vertices.size();
    mesh_ = mesh().SubdivideModifiedButterfly();

    // update elasticity
    auto vertices = vertex_range();
    for (int vid = noldvertices; vid < vertices.size(); vid++) {
        vertices[vid].elasticity = 1;
    }
    for (auto& face : mesh().faces) {
        float elasticity = -1;
        for (int vid : face.index) {
            if (vid < noldvertices) {
                // old vertex, at most one per face
                elasticity = mesh().vertices[vid].elasticity;
                break;
            }
        }
        if (elasticity != -1) {
            for (int vid : face.index) {
                if (mesh().vertices[vid].elasticity > elasticity) {
                    vertices[vid].elasticity = elasticity;    // set to the lowest
                }
            }
        }
    }
    if (landmark_.num_rows()) {
        SubdivideCorrespondence();
    }
}

/// affine transform 3D coordinate matrix
/**
@verbatim
[x,y,z] = [R(x,y,z)+(tx,ty,tz)]
[  :  ]   [        :          ]
@endverbatim
*/
void AffineTransformCoordinateMatrix(
    const slib::CMatrix<double, 4, 4>& trans, // affine transformation
    slib::CMatrix<double>& coordinates // mx3 coordinate matrix
) {
    const int nrows = coordinates.num_rows();
    for (int r = 0; r < nrows; r++) {
        auto p = AffineTransform(trans, make_vector_from_row(coordinates, r));
        for (int c = 0; c < 3; c++) {
            coordinates(r, c) = p[c];
        }
    }
}

void TransformOrg(const CMesh& mesh, slib::CMatrix<double>& org) {
    std::vector<PositionPair> pair(mesh.vertices.size());
    for (int vid = 0; vid < mesh.vertices.size(); vid++) {
        pair[vid].src = make_vector_from_row(org, vid);
        pair[vid].dst = mesh.vertices[vid].position;
        pair[vid].weight = 1;
    }
    double scale;
    slib::CMatrix<double, 3, 3> rot;
    slib::CVector<double, 3> trans;
    SolveSimilarity(pair,
                    &scale,
                    &rot,
                    &trans);
    AffineTransformCoordinateMatrix(make_affine_matrix(scale * rot, trans), org);
}

MHBMLIB_API slib::CMatrix<double> DeformNonrigid(
    const CMesh& src_mesh,
    const CMesh& dst_mesh,
    const slib::CSparseMatrix<double> *landmarkC,
    const slib::CSparseMatrix<double> *landmarkD,
    const slib::CMatrix<double>& src_org,
    float correspondence_weight,
    float landmark_weight,
    //float rigidness_weight,
    bool point_plane_distance,
    TRANSFORMATION model,
    float max_distance2,
    float min_cosangle,
    bool allow_border,
    SEARCH_DIRECTION direction
) {
    // laplacian
    auto L = src_mesh.GenerateLaplacianMatrix(src_org/*, rigidness_weight*/);

    // current working mesh
    auto cur_mesh = src_mesh;
    slib::CMatrix<double> cur_pos, dst_pos;
    src_mesh.ConvertMeshToCoordinate(cur_pos);
    dst_mesh.ConvertMeshToCoordinate(dst_pos);
    if (point_plane_distance) {
        dst_pos = VectorizeCoordinates(dst_pos);
    }

    // if model == transation, src_org should be rotated properly.
    auto src_rotated_org = src_org;
    TransformOrg(src_mesh, src_rotated_org);

    MatrixPair mat;
    MultidirectionalClosestPointSearch mcp;
    if (correspondence_weight > 0) {
        // ICP correspondence
        mcp.SetParameters(&cur_mesh, &dst_mesh, max_distance2, min_cosangle, allow_border, direction);
    } else if (landmark_weight && landmarkC && landmarkD) {
        // landmark
        if (point_plane_distance) {
            mat.C = ExpandLaplacianMatrix(*landmarkC);
            mat.D = ExpandLaplacianMatrix(*landmarkD);
        } else {
            mat.C = *landmarkC;
            mat.D = *landmarkD;
        }
        mat.C.Scale(landmark_weight);
        mat.D.Scale(landmark_weight);
    } else {
        ThrowRuntimeError("invalid parameter");
    }

    double residual = 0;
    double prev_residual = std::numeric_limits<float>::max();
    for (int niters = 0;  ; niters++) {
        // nearest neighbors
        if (correspondence_weight > 0) {
            if (point_plane_distance) {
                mat = mcp.Find().ToNormalMatrix();
            } else {
                mat = mcp.Find().ToPositionMatrix();
            }
            mat.C = mat.W.MultiplyTo('N', mat.C);
            mat.D = mat.W.MultiplyTo('N', mat.D);
            mat.C.Scale(correspondence_weight);
            mat.D.Scale(correspondence_weight);

            if (landmark_weight && landmarkC && landmarkD) {
                if (point_plane_distance) {
                    mat.C.AppendRows(landmark_weight, ExpandLaplacianMatrix(*landmarkC));
                    mat.D.AppendRows(landmark_weight, ExpandLaplacianMatrix(*landmarkD));
                } else {
                    mat.C.AppendRows(landmark_weight, *landmarkC);
                    mat.D.AppendRows(landmark_weight, *landmarkD);
                }
            }
        }

        if (mat.C.num_rows() < 3) {
            ThrowRuntimeError("insufficient correspondence");
        }

        if (point_plane_distance) {
            residual = SolveLaplacianDeformationConstrainedNormal(
                           L, cur_pos, src_rotated_org,
                           mat.D.MultiplyTo('N', dst_pos), mat.C, model);
        } else {
            residual = SolveLaplacianDeformationConstrainedPosition(
                           L, cur_pos, src_rotated_org,
                           mat.D.MultiplyTo('N', dst_pos), mat.C, model);
        }

        std::clog << "|r| = " << sqrt(residual) << ", #pairs = " << mat.C.num_rows() << std::endl;
        cur_mesh.ConvertCoordinateToMesh(cur_pos);
        cur_mesh.UpdateGeometry();

        if ((niters % MIN_ICP_ITERATIONS) == 0) {
            if (prev_residual / residual <= 1 + MIN_ICP_IMPROVEMENT) {
                break;
            } else {
                prev_residual = residual;
            }
        }
        if (correspondence_weight == 0) {
            break;
        }
    }
    return cur_pos;
}

} // namespace hbm