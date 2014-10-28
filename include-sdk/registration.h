#include "MatrixUtil.h"
#include "mkl_lapack_driver.h"

namespace hbm {

/// weighted point correspondence
template <typename T>
struct Correspondence {
    slib::CVector<T, 3> src;///< source
    slib::CVector<T, 3> dst;///< destination
    T weight;///< weight
};

/// solve similarity from point correspondences.
/// dst = s.R.src + t:
/// @return SSD |s.R.src + t - dst|^2
/// @see http://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf
double SolveSimilarity(
    const std::vector<Correspondence<double>>& correspondence,///< vector of correspondence
    double *scale, ///< scale (nullptr for translation/rigidity)
    slib::CMatrix<double, 3, 3> *rot,///< rotation (nullptr for translation)
    slib::CVector<double, 3> *trans ///< translation
);

}
