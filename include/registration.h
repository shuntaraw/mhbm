#include "MatrixUtil.h"
#include "mkl_lapack_driver.h"

namespace slib {
namespace mesh {

template <typename T>
struct Correspondence {
    CVector<T, 3> src;
    CVector<T, 3> dst;
    T weight;
};

/// solve similarity from point correspondences.
/// dst = s.R.src + t:
/// @return SSD |s.R.src + t - dst|^2
/// @see http://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf
template <typename T>
inline
T SolveSimilarity(
    const std::vector<Correspondence<T>>& correspondence,
    T *scale, ///<[in,out] scale (nullptr for translation/rigidity)
    CMatrix<T, 3, 3> *rot,///<[in,out] rotation (nullptr for translation)
    CVector<T, 3> *trans ///<[in,out] translation
) {
    if (correspondence.size() < 3) {
        ThrowRuntimeError("insufficient number of correspondence");
    }
    assert(trans);

    // solve translation
    double sum_weight = 0;
    CVector<double, 3> src_center, dst_center;
    src_center.fill_with(0);
    dst_center.fill_with(0);
    for (auto& c : correspondence) {
        src_center += c.weight * c.src;
        dst_center += c.weight * c.dst;
        sum_weight += c.weight;
    }
    if (sum_weight == 0) {
        ThrowRuntimeError("small weights");
    }
    src_center /= sum_weight;
    dst_center /= sum_weight;

    if (rot == 0) {
        // translation
        assert(scale == 0);
        *trans = dst_center - src_center;
        double residual2 = 0;
        for (auto& c : correspondence) {
            residual2 += c.weight * norm2_squared_of(c.src + *trans - c.dst);
        }
        return residual2 / sum_weight;
    }

    // solve rotation
    CMatrix<double, 3, 3> covar;
    covar.fill_with(0);
    double sum_src2 = 0;
    double sum_dst2 = 0;
    for (auto& c : correspondence) {
        CVector<double, 3> s = c.src - src_center;
        CVector<double, 3> d = c.dst - dst_center;
        covar += s * transpose_of(d) * c.weight;
        sum_src2 += dot(s, s) * c.weight;
        sum_dst2 += dot(d, d) * c.weight;
    }
    CMatrix<double, 3, 3> U, Vt;
    CVector<double, 3> s;
    LAPACKE_gesvd<T,3,3>(covar, U, s, Vt);
    *rot = transpose_of(U * Vt);
    if (determinant_of(*rot) < 0) {
        *rot = transpose_of(U * make_diagonal_matrix(1., 1., -1.) * Vt);
    }
    double consistency = s[0] + s[1] + s[2];

    if (scale) {
        // similarity
        *scale = sqrt(sum_dst2 / sum_src2);
        *trans = dst_center - *scale **rot * src_center;
        return (*scale **scale * sum_src2 + sum_dst2 - 2 * *scale * consistency) / sum_weight;
    } else {
        // rigidity
        *trans = dst_center - *rot * src_center;
        return (sum_src2 + sum_dst2 - 2 * consistency) / sum_weight;
    }
}

}
}
