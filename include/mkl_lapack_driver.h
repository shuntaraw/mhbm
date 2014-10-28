// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <boost/timer/timer.hpp>
#include <mkl_lapack.h>
#include "MatrixBase.h"
#include "exception_util.h"
#include "mkl_driver.h"

namespace slib {

/// Cholesky factorization of a symmetric (Hermitian) positive-definite matrix
inline int LAPACKE_potrf(int layout, char uplo, int n, float *a, int lda) {
    LAPACKE_spotrf(layout, uplo, n, a, lda);
}
inline int LAPACKE_potrf(int layout, char uplo, int n, double *a, int lda) {
    LAPACKE_dpotrf(layout, uplo, n, a, lda);
}
template <typename T, int NR, int NC>
inline void LAPACKE_potrf(
    CMatrix<T, NR, NC>& mat ///< lower triangle
) {
    assert(mat.num_rows() == mat.num_cols());
    int info = LAPACKE_potrf(LAPACK_COL_MAJOR, 'L', mat.num_rows(), mat.data(), mat.num_rows());
    if (info < 0) {
        ThrowRuntimeError("%d-th parameter had an illegal value in DPFTRF()", -info);
    } else if (info > 0) {
        ThrowRuntimeError("the leading minor of order %d is not positive-definite in DPFTRF()", info);
    }
}

/// solve linear least squares (LLS) problems by QR or LQ factorization.
inline int LAPACKE_gels(int layout, char trans, int m,
                        int n, int nrhs, float *a,
                        int lda, float *b, int ldb) {
    return LAPACKE_sgels(layout, trans, m, n, nrhs, a, lda, b, ldb);
}
inline int LAPACKE_gels(int layout, char trans, int m,
                        int n, int nrhs, double *a,
                        int lda, double *b, int ldb) {
    return LAPACKE_dgels(layout, trans, m, n, nrhs, a, lda, b, ldb);
}
template <typename T, int NR, int NC, int NRhs>
inline
void LAPACKE_gels(
    const CMatrix<T, NR, NC>& mat,
    CMatrix<T, NR, NRhs>& vec ///< vector, solution
) {
    assert(vec.num_rows() == std::max(mat.num_rows(), mat.num_cols()));
    auto A = mat;
    int info = LAPACKE_gels(LAPACK_COL_MAJOR, 'N', A.num_rows(), A.num_cols(), vec.num_cols(), A.data(), A.num_rows(), vec.data(), vec.num_rows());
    if (info < 0) {
        ThrowRuntimeError("%d-th parameter had an illegal value in DGELS()", -info);
    } else if (info > 0) {
        ThrowRuntimeError("%d-th diagonal element of the triangular factor is zero in DGELS()", info);
    }
    vec.resize(mat.num_cols(), vec.num_cols());
}

/// Singular Value Decomposition.
inline int LAPACKE_gesvd(int layout, char jobu, char jobvt,
                         int m, int n, float *a, int lda,
                         float *s, float *u, int ldu, float *vt,
                         int ldvt, float *superb) {
    return LAPACKE_sgesvd(layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}
inline int LAPACKE_gesvd(int layout, char jobu, char jobvt,
                         int m, int n, double *a,
                         int lda, double *s, double *u, int ldu,
                         double *vt, int ldvt, double *superb) {
    return LAPACKE_dgesvd(layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}
template <typename T, int NRows, int NCols>
inline
void LAPACKE_gesvd(
    const CMatrix<T, NRows, NCols>& mat, ///<
    CMatrix<T, NRows, NRows>& matU, ///<
    CVector < T, (NRows < NCols ? NRows : NCols) > & vecW, ///<
    CMatrix<T, NCols, NCols>& matVt ///<
) {
    assert(mat.num_rows() == matU.num_rows());
    assert(mat.num_rows() == matU.num_cols());
    assert(mat.num_cols() == matVt.num_rows());
    assert(mat.num_cols() == matVt.num_cols());
    assert(vecW.num_rows() == std::min(mat.num_rows(), mat.num_cols()));
    auto A = mat; // copy
    slib::CVector<T> superb(std::min(mat.num_rows(), mat.num_cols()) - 1);
    int info = LAPACKE_gesvd(LAPACK_COL_MAJOR, 'A', 'A', A.num_rows(), A.num_cols(), A.data(), A.num_rows(), vecW.data(), matU.data(), matU.num_rows(), matVt.data(), matVt.num_rows(), superb.data());
    if (info < 0) {
        ThrowRuntimeError("%d-th parameter had an illegal value in ?GESVD()", -info);
    } else if (info > 0) {
        ThrowRuntimeError("%d superdiagonals of the intermediate bidiagonal form did not converge to zero in ?GESVD()", info);
    }
}

/// Computes all eigenvalues and, optionally, eigenvectors of a real symmetric / Hermitian matrix.
inline int LAPACKE_syev(int layout, char jobz, char uplo, int n,
                        float *a, int lda, float *w) {
    return LAPACKE_ssyev(layout, jobz, uplo, n, a, lda, w);
}
inline int LAPACKE_syev(int layout, char jobz, char uplo, int n,
                        double *a, int lda, double *w) {
    return LAPACKE_dsyev(layout, jobz, uplo, n, a, lda, w);
}
template <typename T, int NRows>
inline void LAPACKE_syev(
    CMatrix<T, NRows, NRows>& mat, ///< upper triangle
    CVector<T, NRows>& value  ///<
) {
    assert(mat.num_rows() == mat.num_cols());
    assert(mat.num_cols() == value.num_rows());
    int info = LAPACKE_syev(LAPACK_COL_MAJOR, 'V', 'U', mat.num_rows(), mat.data(), mat.num_rows(), value.data());
    if (info < 0) {
        ThrowRuntimeError("%d-th parameter of ?SYEV() is invalid" , -info);
    } else if (info > 0) {
        ThrowRuntimeError("?SYEV() did not converge");
    }
}

/// Computes all eigenvalues and, optionally, eigenvectors of a real symmetric / Hermitian matrix using the Relatively Robust Representations.
inline int LAPACKE_syevr(int  layout, char jobz, char range, char uplo,
                         int n, float *a, int lda, float vl,
                         float vu, int il, int iu, float abstol,
                         int *m, float *w, float *z, int ldz,
                         int *isuppz) {
    return LAPACKE_ssyevr(layout, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz);
}
inline int LAPACKE_syevr(int  layout, char jobz, char range, char uplo,
                         int n, double *a, int lda, double vl,
                         double vu, int il, int iu,
                         double abstol, int *m, double *w, double *z,
                         int ldz, int *isuppz) {
    return LAPACKE_dsyevr(layout, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz);
}
template <typename T, int NR1, int NC1, int NR2>
inline void LAPACKE_syevr(
    CMatrix<T, NR1, NC1>& mat, ///< upper triangle, eigenvectors in the first columns
    CVector<T, NR2>& value,  ///< eigenvectors
    char small_or_large ///< 'L' larger / 'S' smaller eigenvalues
) {
    assert(mat.num_rows() == mat.num_cols());
    int n = mat.num_rows();
    int m = value.num_rows();
    int il, iu;
    switch (small_or_large) {
    case 'L':
        // greater
        il = n - m + 1;
        iu = n;
        break;
    case 'S':
        // smaller
        il = 1;
        iu = m;
        break;
    default:
        ThrowLogicError("undefined label");
    }
    CVector<T> w(n);
    CMatrix<T> z(n, m);
    std::vector<int> isuppz(2 * m);
    int info = LAPACKE_syevr(LAPACK_COL_MAJOR, 'V', 'I', 'U', n, mat.data(), n, 0, 0, il, iu, 0, &m, w.data(), z.data(), n, isuppz.data());
    if (info < 0) {
        ThrowRuntimeError("%d-th parameter of ?SYEVR() is invalid", -info);
    } else if (info > 0) {
        ThrowRuntimeError("?SYEVR() did not converge");
    }
    mat.copy_from(z.data());
    value.copy_from(w.data());
}

/// Computes the LU factorization of a general m-by-n matrix.
inline int LAPACKE_getrf(int layout, int m, int n,
                         float *a, int lda, int *ipiv) {
    return LAPACKE_sgetrf(layout, m, n, a, lda, ipiv);
}
inline int LAPACKE_getrf(int layout, int m, int n,
                         double *a, int lda, int *ipiv) {
    return LAPACKE_dgetrf(layout, m, n, a, lda, ipiv);
}
template <typename T, int NRows, int NCols, int NPivs>
inline void LAPACKE_getrf(
    CMatrix<T, NRows, NCols>& mat, ///< matrix, L (except unit diagonal) and U.
    CVector<int, NPivs>& pivot  ///< pivot indices
) {
    assert(pivot.num_rows() == std::min(mat.num_rows(), mat.num_cols()));
    int info = LAPACKE_getrf(LAPACK_COL_MAJOR, mat.num_rows(), mat.num_cols(), mat.data(), mat.num_rows(), pivot.data());
    if (info < 0) {
        ThrowRuntimeError("%d-th parameter had an illegal value in ?GETRF()", -info);
    } else if (info > 0) {
        std::clog << "warining: singular matrix in ?GETRF()" << std::endl;
    }
}

/// Solves a system of linear equations with an LU-factored square matrix, with multiple right-hand sides.
inline int LAPACKE_getrs(int layout, char trans, int n,
                         int nrhs, const float *a, int lda,
                         const int *ipiv, float *b, int ldb) {
    return LAPACKE_sgetrs(layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
}
inline int LAPACKE_getrs(int layout, char trans, int n,
                         int nrhs, const double *a, int lda,
                         const int *ipiv, double *b, int ldb) {
    return LAPACKE_dgetrs(layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
}
template <typename T, int NRows, int NRhs>
inline void LAPACKE_getrs(
    const CMatrix<T, NRows, NRows>& mat, ///< L (except unit diagonal) and U
    const CVector<int, NRows>& pivot,  ///< pivot indices
    CMatrix<T, NRows, NRhs>& vec///<
) {
    assert(mat.num_rows() == mat.num_cols());
    assert(pivot.num_rows() == mat.num_rows());
    assert(vec.num_rows() == mat.num_rows());
    int info = LAPACKE_getrs(LAPACK_COL_MAJOR, 'N', mat.num_rows(), vec.num_cols(), mat.data(), mat.num_rows(), pivot.data(), vec.data(), vec.num_rows());
    if (info < 0) {
        ThrowRuntimeError("%d-th parameter had an illegal value in ?GETRS()", -info);
    }
}

/// Computes the inverse of an LU-factored general matrix.
inline int LAPACKE_sgetri(int layout, int n, float *a,
                          int lda, const int *ipiv) {
    LAPACKE_sgetri(layout, n, a, lda, ipiv);
}
inline int LAPACKE_dgetri(int layout, int n, double *a,
                          int lda, const int *ipiv) {
    LAPACKE_dgetri(layout, n, a, lda, ipiv);
}
template <typename T, int NRows>
inline void LAPACKE_getri(
    CMatrix<T, NRows, NRows>& mat, ///< L (except unit diagonal) and U, inverse
    const CVector<int, NRows>& pivot  ///< pivot indices
) {
    assert(mat.num_rows() == mat.num_cols());
    assert(pivot.num_rows() == mat.num_rows());
    int info = LAPACKE_getri(LAPACK_COL_MAJOR, mat.num_rows(), mat.data(), mat.num_rows(), pivot.data());
    if (info < 0) {
        ThrowRuntimeError("%d-th parameter had an illegal value in ?GETRI()", -info);
    } else if (info > 0) {
        ThrowRuntimeError("%d-th diagonal element of the factor U is zero in ?GETRI()", info);
    }
}

// diagonal scaling on a vector
// x = diag(d)*x
inline void LASCL2(int m, int n, const float *d, float *x, int incx) {
    SLASCL2(&m, &n, d, x, &incx);
}
inline void LASCL2(int m, int n, const double *d, double *x, int incx) {
    DLASCL2(&m, &n, d, x, &incx);
}
template <typename T, int NR1, int NR2, int NC2>
inline void LASCL2(const CVector<T, NR1>& d, CMatrix<T, NR2, NC2>& x) {
    assert(d.num_rows() == x.num_rows());
    int m = d.num_rows();
    int n = 1;
    int incx = 1;
    for (int c = 0; c < x.num_cols(); c++) {
        LASCL2(m, n, d.data(), x.col_data(c), incx);
    }
}

// reciprocal diagonal scaling on a vector
inline void LARSCL2(int m, int n, const float *D, float *x, int incx) {
    SLARSCL2(&m, &n, D, x, &incx);
}
inline void LARSCL2(int m, int n, const double *D, double *x, int incx) {
    DLARSCL2(&m, &n, D, x, &incx);
}
template <typename T, int NR1, int NR2, int NC2>
inline void LARSCL2(const CVector<T, NR1>& D, CMatrix<T, NR2, NC2>& x) {
    assert(D.num_rows() == x.num_rows());
    int m = D.num_rows();
    int n = 1;
    int incx = 1;
    for (int c = 0; c < x.num_cols(); c++) {
        LARSCL2(m, n, D.data(), x.col_data(c), incx);
    }
}

/// compute the determinant of a matrix
/// @tparam T  element
/// @tparam NRows number of rows
template <typename T, int NRows>
inline T determinant_of(
    const CMatrix<T, NRows, NRows>& a
) {
    auto mat = a;
    CVector<int, NRows> pivot(a.num_rows());
    LAPACKE_getrf(mat, pivot);
    double det = 1;
    for (int i = 0; i < a.num_rows(); i++) {
        if (pivot[i] != i + 1) {
            det *= -mat(i, i);
        } else {
            det *= mat(i, i);
        }
    }
    return det;
}

/// compute an inverse matrix
/// @tparam T  element
/// @tparam NRows number of rows
template <typename T, int NRows>
inline CMatrix<T, NRows, NRows> inverse_of(const CMatrix<T, NRows, NRows>& mat) {
    auto A = mat;
    CVector<int, NRows> pivot(A.num_rows());
    LAPACKE_getrf(A, pivot);
    LAPACKE_getri(A, pivot);
    return A;
}

/// find a rotation matrix closest to a given matrix wrt the Frobenius norm
/// @tparam T  element
/// @tparam NRows number of rows
template <typename T, int NRows>
inline
void NormalizeRotation(CMatrix<T, NRows, NRows>& mat) {
    assert((mat.num_rows() == 3 && mat.num_cols() == 3) || (mat.num_rows() == 4 && mat.num_cols() == 4));
    CMatrix<T, 3, 3> R;
    for (int c = 0; c < 3; c++) {
        for (int r = 0; r < 3; r++) {
            R(r, c) = mat(r, c);
        }
    }
    CMatrix<T, 3, 3> U, Vt;
    CVector<T, 3> W;
    LAPACKE_gesvd<T, 3, 3>(R, U, W, Vt);
    R = U * Vt;
    for (int c = 0; c < 3; c++) {
        for (int r = 0; r < 3; r++) {
            mat(r, c) = R(r, c);
        }
    }
}

} // namespace slib
