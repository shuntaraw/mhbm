// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include "mkl_driver.h"
#include <fstream>
#include <mkl_blas.h>
#include "MatrixUtil.h"

namespace slib {

// CBLAS drivers. scalar is typed as 'double' to resolve ambiguity in overloading

inline void cblas_axpy(int n, float alpha, const float *x, int incx, float *y, int incy) {
    assert(x != y);
    cblas_saxpy(n, alpha, x, incx, y, incy);
}
inline void cblas_axpy(int n, double alpha, const double *x, int incx, double *y, int incy) {
    assert(x != y);
    cblas_daxpy(n, alpha, x, incx, y, incy);
}
template <typename T, int NR1, int NC1, int NR2, int NC2>
inline void cblas_axpy(double alpha, const CMatrix<T, NR1, NC1>& x, CMatrix<T, NR2, NC2>& y) {
    assert(x.num_rows() == y.num_rows());
    assert(x.num_cols() == y.num_cols());
    cblas_axpy(x.num_rows()*x.num_cols(), T(alpha), x.data(), 1, y.data(), 1);
}

inline void cblas_axpby(int n, float alpha, const float *x, int incx, float beta, float *y, int incy) {
    assert(x != y);
    cblas_saxpby(n, alpha, x, incx, beta, y, incy);
}
inline void cblas_axpby(int n, double alpha, const double *x, int incx, double beta, double *y, int incy) {
    assert(x != y);
    cblas_daxpby(n, alpha, x, incx, beta, y, incy);
}
template <typename T, int NR1, int NC1, int NR2, int NC2>
inline void cblas_axpby(double a, CMatrix<T, NR1, NC1>& x, double beta, const CMatrix<T, NR2, NC2>& y) {
    assert(x.num_rows() == y.num_rows());
    assert(x.num_cols() == y.num_cols());
    cblas_axpby(x.num_rows()*x.num_cols(), T(alpha), x.data(), 1, beta, y.data(), 1);
}

inline float cblas_dot(int n, const float  *x, int incx, const float  *y, int incy) {
    return cblas_sdot(n, x, incx, y, incy);
}
inline double cblas_dot(int n, const double  *x, int incx, const double  *y, int incy) {
    return cblas_ddot(n, x, incx, y, incy);
}
template <typename T, int NR1, int NC1, int NR2, int NC2>
inline T cblas_dot(const CMatrix<T, NR1, NC1>& x, const CMatrix<T, NR2, NC2>& y) {
    assert(x.num_rows() == y.num_rows());
    assert(x.num_cols() == y.num_cols());
    return cblas_dot(x.num_rows() * x.num_cols(), x.data(), 1, y.data(), 1);
}

inline float cblas_asum(int n, const float *x, int incx) {
    return cblas_sasum(n, x, incx);
}
inline double cblas_asum(int n, const double *x, int incx) {
    return cblas_dasum(n, x, incx);
}
template <typename T, int NR, int NC>
inline T cblas_asum(const CMatrix<T, NR, NC>& x) {
    return cblas_asum(x.num_rows() * x.num_cols(), x.data(), 1);
}

inline void cblas_gemv(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transa, int m, int n, float alpha, const float *A, int lda, const float *x, int incx, float beta, float *y, int incy) {
    assert(A != y);
    assert(x != y);
    cblas_sgemv(layout, transa, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
inline void cblas_gemv(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transa, int m, int n, double alpha, const double *A, int lda, const double *x, int incx, double beta, double *y, int incy) {
    assert(A != y);
    assert(x != y);
    cblas_dgemv(layout, transa, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
template <typename T, int NR1, int NC1, int NR2, int NC2, int NR3, int NC3>
inline void cblas_gemv(CBLAS_TRANSPOSE transa, double alpha, const CMatrix<T, NR1, NC1>& A, const CMatrix<T, NR2, NC2>& x, double beta, CMatrix<T, NR3, NC3>& y) {
    assert(y.num_rows() == (transa == CblasTrans ? A.num_cols() : A.num_rows()));
    assert(x.num_rows() == (transa == CblasTrans ? A.num_rows() : A.num_cols()));
    assert(x.num_cols() == 1);
    assert(y.num_cols() == 1);
    return cblas_gemv(CblasColMajor, transa, A.num_rows(), A.num_cols(), T(alpha), A.data(), A.num_rows(), x.data(), 1, beta, y.data(), 1);
}

inline void cblas_gemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb, int m, int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc) {
    assert(A != C);
    assert(B != C);
    cblas_sgemm(layout, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
inline void cblas_gemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb, int m, int n, int k, double alpha, const double *A, int lda, const double *B, int ldb, double beta, double *C, int ldc) {
    assert(A != C);
    assert(B != C);
    cblas_dgemm(layout, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
template <typename T, int NR1, int NC1, int NR2, int NC2, int NR3, int NC3>
inline void cblas_gemm(CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb, double alpha, const CMatrix<T, NR1, NC1>& A, const CMatrix<T, NR2, NC2>& B, double beta,  CMatrix<T, NR3, NC3>& C) {
    int m = (transa == CblasTrans) ? A.num_cols() : A.num_rows();
    int n = (transb == CblasTrans) ? B.num_rows() : B.num_cols();
    int k = (transa == CblasTrans) ? A.num_rows() : A.num_cols();
    assert(m == C.num_rows());
    assert(n == C.num_cols());
    assert(k == ((transb == CblasTrans) ? B.num_cols() : B.num_rows()));
    cblas_gemm(CblasColMajor, transa, transb, m, n, k, T(alpha), A.data(), A.num_rows(), B.data(), B.num_rows(), T(beta), C.data(), C.num_rows());
}

inline void cblas_scal(int n, float alpha, float *x, int incx) {
    return cblas_sscal(n, alpha, x, incx);
}
inline void cblas_scal(int n, double alpha,  double *x, int incx) {
    cblas_dscal(n, alpha, x, incx);
}
template <typename T, int NR, int NC>
inline void cblas_scal(double alpha, const CMatrix<T, NR, NC>& x) {
    cblas_scal(x.num_rows()*x.num_cols(), T(alpha), x.data(), 1);
}

inline float cblas_nrm2(int n, const float *x, int incx) {
    return cblas_snrm2(n, x, incx);
}
inline double cblas_nrm2(int n, const double *x, int incx) {
    return cblas_dnrm2(n, x, incx);
}
template <typename T, int NR, int NC>
inline T cblas_nrm2(const CMatrix<T, NR, NC>& x) {
    return cblas_nrm2(x.num_rows() * x.num_cols(), x.data(), 1);
}

// y = alpha*A*x+beta*y where A is diagonal
inline void cblas_sbmv(int n, float alpha, const float *a, const float *x, float beta, float *y) {
    cblas_ssbmv(CblasColMajor, CblasUpper, n, 0, alpha, a, 1, x, 1, beta, y, 1);
}
inline void cblas_sbmv(int n, double alpha, const double *a, const double *x, double beta, double *y) {
    cblas_dsbmv(CblasColMajor, CblasUpper, n, 0, alpha, a, 1, x, 1, beta, y, 1);
}
template <typename T, int NR1, int NR2, int NR3>
inline void cblas_sbmv(double alpha, const CVector<T, NR1>& A, const CVector<T, NR2>& x, double beta, CVector<T, NR3>& y) {
    assert(A.num_rows() == x.num_rows());
    assert(A.num_rows() == y.num_rows());
    return cblas_sbmv(A.num_rows(), T(alpha), A.data(), x.data(), T(beta), y.data());
}

// Y = alpha*A*X+beta*Y where A is diagonal
inline void cblas_sbmm(int m, int n, float alpha, const float *a, const float *x, float beta, float *y) {
    for (int c = 0; c < n; c++) {
        cblas_ssbmv(CblasColMajor, CblasUpper, m, 0, alpha, a, 1, x + c * m, 1, beta, y + c * m, 1);
    }
}
inline void cblas_sbmm(int m, int n, double alpha, const double *a, const double *x, double beta, double *y) {
    for (int c = 0; c < n; c++) {
        cblas_dsbmv(CblasColMajor, CblasUpper, m, 0, alpha, a, 1, x + c * m, 1, beta, y + c * m, 1);
    }
}
template <typename T, int NR1, int NR2, int NC2, int NR3, int NC3>
inline void cblas_sbmm(double alpha, const CVector<T, NR1>& A, const CMatrix<T, NR2, NC2>& X, double beta, CMatrix<T, NR3, NC3>& Y) {
    assert(A.num_rows() == X.num_rows());
    assert(A.num_rows() == Y.num_rows());
    assert(X.num_cols() == Y.num_cols());
    return cblas_sbmm(A.num_rows(), X.num_cols(), T(alpha), A.data(), X.data(), T(beta), Y.data());
}

// utilities for CMatrix

template <int NRows, int NCols, typename Scalar>
inline void ScaleMatrixInPlace(CMatrix<float, NRows, NCols>& m, Scalar s) {
    mkl_simatcopy('C', 'N', m.num_rows(), m.num_cols(), float(s), m.data(), m.num_rows(), m.num_rows());
}

template <int NRows, int NCols, typename Scalar>
inline void ScaleMatrixInPlace(CMatrix<double, NRows, NCols>& m, Scalar s) {
    mkl_dimatcopy('C', 'N', m.num_rows(), m.num_cols(), double(s), m.data(), m.num_rows(), m.num_rows());
}

template <int NRows1, int NCols1, int NRows2, int NCols2, typename Scalar>
inline void ScaleAddMatrixInPlace(CMatrix<float, NRows1, NCols1>& m1, Scalar s, const CMatrix<float, NRows2, NCols2>& m2) {
    cblas_axpy(float(s), m2, m1);
}

template <int NRows1, int NCols1, int NRows2, int NCols2, typename Scalar>
inline void ScaleAddMatrixInPlace(CMatrix<double, NRows1, NCols1>& m1, Scalar s, const CMatrix<double, NRows2, NCols2>& m2) {
    cblas_axpy(double(s), m2, m1);
}

template <int NRows1, int NCols1, int NRows2, int NCols2, int NRows3, int NCols3>
inline void MultiplyMatrix(
    CMatrix<float, NRows1, NCols1>& m1,
    const CMatrix<float, NRows2, NCols2>& m2,
    const CMatrix<float, NRows3, NCols3>& m3
) {
    cblas_gemm(CblasNoTrans, CblasNoTrans, 1, m2, m3, 0, m1);
}

template <int NRows1, int NCols1, int NRows2, int NCols2, int NRows3, int NCols3>
inline void MultiplyMatrix(
    CMatrix<double, NRows1, NCols1>& m1,
    const CMatrix<double, NRows2, NCols2>& m2,
    const CMatrix<double, NRows3, NCols3>& m3
) {
    cblas_gemm(CblasNoTrans, CblasNoTrans, 1, m2, m3, 0, m1);
}

template <int NRows, int NCols>
inline float norm2_of(const CMatrix<float, NRows, NCols>& mat) {
    return cblas_nrm2(mat);
}

template <int NRows, int NCols>
inline double norm2_of(const CMatrix<double, NRows, NCols>& mat) {
    return cblas_nrm2(mat);
}

template <int NRows, int NCols>
inline float norm2_squared_of(const CMatrix<float, NRows, NCols>& mat) {
    return cblas_dot(mat, mat);
}

template <int NRows, int NCols>
inline double norm2_squared_of(const CMatrix<double, NRows, NCols>& mat) {
    return cblas_dot(mat, mat);
}

} // namespace slib
