#include "MatrixUtil.h"
#include "mkl_lapack_driver.h"

namespace hbm {

template <typename T>
struct Correspondence {
    slib::CVector<T, 3> src;
    slib::CVector<T, 3> dst;
    T weight;
};

/// solve similarity from point correspondences.
/// dst = s.R.src + t:
/// @return SSD |s.R.src + t - dst|^2
/// @see http://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf
double SolveSimilarity(
    const std::vector<Correspondence<double>>& correspondence,
    double *scale, ///<[in,out] scale (nullptr for translation/rigidity)
    slib::CMatrix<double, 3, 3> *rot,///<[in,out] rotation (nullptr for translation)
    slib::CVector<double, 3> *trans ///<[in,out] translation
    );

}
