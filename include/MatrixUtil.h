// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include "MatrixBase.h"
#include "exception_util.h"

namespace slib {

/// compute the infinity norm
template <typename T, int NRows, int NCols>
inline
T norm_infty_of(const CMatrix<T, NRows, NCols>& mat) {
    T norm = 0;
    for (int i = 0; i < mat.num_rows()*mat.num_cols(); i++) {
        norm = std::max(norm, abs(mat[i]));
    }
    return norm;
}

/// compute the 1-norm
template <typename T, int NRows, int NCols>
inline
T norm1_of(const CMatrix<T, NRows, NCols>& mat) {
    T norm = 0;
    for (int i = 0; i < mat.num_rows()*mat.num_cols(); i++) {
        norm += abs(mat[i]);
    }
    return norm;
}

/// compute the square of 2-norm
template <typename T, int NRows, int NCols>
inline
T norm2_squared_of(const CMatrix<T, NRows, NCols>& mat) {
    T norm = 0;
    for (int i = 0; i < mat.num_rows()*mat.num_cols(); i++) {
        norm += mat[i] * mat[i];
    }
    return norm;
}

/// compute the 2-norm
template <typename T, int NRows, int NCols>
inline
decltype(std::sqrt(T())) norm2_of(
    const CMatrix<T, NRows, NCols>& mat //
) {
    return std::sqrt(norm2_squared_of(mat));
}

/// compute a unit vector
template <typename T, int NRows, int NCols>
inline
CMatrix < decltype(T() / std::sqrt(T())), NRows, NCols > normalized_of(const CMatrix<T, NRows, NCols>& vec) {
    return CMatrix<T, NRows, NCols>(vec) /= norm2_of(vec);
}

/// compute a dot product
/// @tparam T element
/// @tparam NRows number of rows
/// @tparam NCols number of columns
template <typename T, int NRows, int NCols>
inline
T dot(
    const CMatrix<T, NRows, NCols>& v1, //
    const CMatrix<T, NRows, NCols>& v2 //
) {
    assert(v1.num_rows() == v2.num_rows());
    assert(v1.num_cols() == 1 && v2.num_cols() == 1);
    T ret = 0;
    for (int i = 0; i < v1.num_rows(); i++) {
        ret += v1[i] * v2[i];
    }
    return ret;
}

/// compute a cross product
/// @tparam T element
template <typename T>
inline
CVector<T, 3> cross(const CVector<T, 3>& v1, const CVector<T, 3>& v2) {
    return CVector <T, 3> {
        v1(1, 0) *v2(2, 0) - v1(2, 0) *v2(1, 0),
        v1(2, 0) *v2(0, 0) - v1(0, 0) *v2(2, 0),
        v1(0, 0) *v2(1, 0) - v1(1, 0) *v2(0, 0)
    };
}

/// compute a cross product
/// @tparam T element
template <typename T>
inline
T cross(const CVector<T, 2>& v1, const CVector<T, 2>& v2) {
    return v1(0, 0) * v2(1, 0) - v1(1, 0) * v2(0, 0);
}

/// construct a skew-symmetric matrix
template <typename T>
inline
CMatrix<T, 3, 3> skew_symmetric_of(CVector<T, 3>& vec) {
    return make_matrix<T, 3, 3>(
               0, -vec(2, 0), vec(1, 0),
               vec(2, 0), 0, -vec(0, 0),
               -vec(1, 0), vec(0, 0), 0
           );
}

/// construct  a homogeneous vector (a vector with an extra element of 1)
/// @tparam T element
/// @tparam NRows number of rows
/// @return homogeneous vector
template <typename T, int NRows>
inline
CMatrix < T, NRows + 1, 1 > homogeneous_of(
    const CVector<T, NRows>& vec // euclidean vector
) {
    CMatrix < T, NRows + 1, 1 > iret;
    for (int i = 0; i < NRows; i++) {
        iret[i] = vec[i];
    }
    iret[NRows] = 1;
    return iret;
}

/// construct  a euclidean vector (a vector divided by an additional last element)
/// @tparam T element
/// @tparam NRows number of rows
/// @return euclidean vector
template <typename T, int NRows>
inline
CMatrix < T, NRows - 1, 1 > euclidean_of(
    const CVector<T, NRows>& vec // homogeneous vector
) {
    CMatrix < T, NRows - 1, 1 > iret;
    for (int i = 0; i < NRows - 1; i++) {
        iret[i] = vec[i] / vec[NRows - 1];
    }
    return iret;
}

/// compute the trace of a matrix
template <typename T, int NRows>
inline
T trace_of(
    const CMatrix<T, NRows, NRows>& mat // matrix
) {
    assert(mat.num_rows() == mat.num_cols());
    T t = T(0);
    for (int i = 0; i < mat.num_rows(); i++) {
        t += mat(i, i);
    }
    return t;
}

/// compute a transposed matrix
template <typename T, int NRows, int NCols>
inline
CMatrix<T, NCols, NRows> transpose_of(
    const CMatrix<T, NRows, NCols>& mat
) {
    CMatrix<T, NCols, NRows> iret(mat.num_cols(), mat.num_rows());
    for (int c = 0; c < mat.num_cols(); c++) {
        for (int r = 0; r < mat.num_rows(); r++) {
            iret(c, r) = mat(r, c);
        }
    }
    return iret;
}

/// compute the determinant of a matrix
/// @tparam T element
template <typename T>
inline
T determinant_of(const CMatrix<T, 2, 2>& a) {
    return a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
}

/// compute the determinant of a matrix
/// @tparam T element
template <typename T>
inline
T determinant_of(const CMatrix<T, 3, 3>& a) {
    return a(0, 0) * a(1, 1) * a(2, 2) +
           a(0, 1) * a(1, 2) * a(2, 0) +
           a(0, 2) * a(1, 0) * a(2, 1) -
           a(0, 0) * a(1, 2) * a(2, 1) -
           a(0, 1) * a(1, 0) * a(2, 2) -
           a(0, 2) * a(1, 1) * a(2, 0);
}

/// compute the determinant of a matrix
/// @tparam T element
template <typename T>
inline
T determinant_of(const CMatrix<T, 4, 4>& a) {
    return ((a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * a(2, 2) + (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)) * a(2, 1) + (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * a(2, 0)) * a(3, 3) +
           ((a(0, 1) * a(1, 0) - a(0, 0) * a(1, 1)) * a(2, 3) + (a(0, 0) * a(1, 3) - a(0, 3) * a(1, 0)) * a(2, 1) + (a(0, 3) * a(1, 1) - a(0, 1) * a(1, 3)) * a(2, 0)) * a(3, 2) +
           ((a(0, 0) * a(1, 2) - a(0, 2) * a(1, 0)) * a(2, 3) + (a(0, 3) * a(1, 0) - a(0, 0) * a(1, 3)) * a(2, 2) + (a(0, 2) * a(1, 3) - a(0, 3) * a(1, 2)) * a(2, 0)) * a(3, 1) +
           ((a(0, 2) * a(1, 1) - a(0, 1) * a(1, 2)) * a(2, 3) + (a(0, 1) * a(1, 3) - a(0, 3) * a(1, 1)) * a(2, 2) + (a(0, 3) * a(1, 2) - a(0, 2) * a(1, 3)) * a(2, 1)) * a(3, 0);
}

/// compute an inverse matrix
/// @tparam T element
template <typename T>
inline
CMatrix<T, 2, 2> inverse_of(const CMatrix<T, 2, 2>& a) {
    T det = determinant_of(a);
    if (det == 0) {
        ThrowRuntimeError("singular matrix");
    }
    return make_vector(
               a(1, 1), -a(1, 0),
               -a(0, 1), a(0, 0)
           ) / det;
}

/// compute an inverse matrix
/// @tparam T element
template <typename T>
inline
CMatrix<T, 3, 3> inverse_of(const CMatrix<T, 3, 3>& a) {
    T det = determinant_of(a);
    if (det == 0) {
        ThrowRuntimeError("singular matrix");
    }
    return make_matrix(
               (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)),
               (a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2)),
               (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)),

               (a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2)),
               (a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0)),
               (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)),

               (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0)),
               (a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1)),
               (a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0))
           ) / det;
}

/// compute an inverse matrix
/// @tparam T element
template <typename T>
inline
CMatrix<T, 4, 4> inverse_of(const CMatrix<T, 4, 4>& a) {
    T det = determinant_of(a);
    if (det == 0) {
        ThrowRuntimeError("singular matrix");
    }
    return make_matrix(
               ((a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) * a(3, 3) + (a(1, 3) * a(2, 1) - a(1, 1) * a(2, 3)) * a(3, 2) + (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)) * a(3, 1)),
               -((a(0, 1) * a(2, 2) - a(0, 2) * a(2, 1)) * a(3, 3) + (a(0, 3) * a(2, 1) - a(0, 1) * a(2, 3)) * a(3, 2) + (a(0, 2) * a(2, 3) - a(0, 3) * a(2, 2)) * a(3, 1)),
               ((a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * a(3, 3) + (a(0, 3) * a(1, 1) - a(0, 1) * a(1, 3)) * a(3, 2) + (a(0, 2) * a(1, 3) - a(0, 3) * a(1, 2)) * a(3, 1)),
               -((a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * a(2, 3) + (a(0, 3) * a(1, 1) - a(0, 1) * a(1, 3)) * a(2, 2) + (a(0, 2) * a(1, 3) - a(0, 3) * a(1, 2)) * a(2, 1)),

               -((a(1, 0) * a(2, 2) - a(1, 2) * a(2, 0)) * a(3, 3) + (a(1, 3) * a(2, 0) - a(1, 0) * a(2, 3)) * a(3, 2) + (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)) * a(3, 0)),
               ((a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0)) * a(3, 3) + (a(0, 3) * a(2, 0) - a(0, 0) * a(2, 3)) * a(3, 2) + (a(0, 2) * a(2, 3) - a(0, 3) * a(2, 2)) * a(3, 0)),
               -((a(0, 0) * a(1, 2) - a(0, 2) * a(1, 0)) * a(3, 3) + (a(0, 3) * a(1, 0) - a(0, 0) * a(1, 3)) * a(3, 2) + (a(0, 2) * a(1, 3) - a(0, 3) * a(1, 2)) * a(3, 0)),
               ((a(0, 0) * a(1, 2) - a(0, 2) * a(1, 0)) * a(2, 3) + (a(0, 3) * a(1, 0) - a(0, 0) * a(1, 3)) * a(2, 2) + (a(0, 2) * a(1, 3) - a(0, 3) * a(1, 2)) * a(2, 0)),

               ((a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0)) * a(3, 3) + (a(1, 3) * a(2, 0) - a(1, 0) * a(2, 3)) * a(3, 1) + (a(1, 1) * a(2, 3) - a(1, 3) * a(2, 1)) * a(3, 0)),
               -((a(0, 0) * a(2, 1) - a(0, 1) * a(2, 0)) * a(3, 3) + (a(0, 3) * a(2, 0) - a(0, 0) * a(2, 3)) * a(3, 1) + (a(0, 1) * a(2, 3) - a(0, 3) * a(2, 1)) * a(3, 0)),
               ((a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * a(3, 3) + (a(0, 3) * a(1, 0) - a(0, 0) * a(1, 3)) * a(3, 1) + (a(0, 1) * a(1, 3) - a(0, 3) * a(1, 1)) * a(3, 0)),
               -((a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * a(2, 3) + (a(0, 3) * a(1, 0) - a(0, 0) * a(1, 3)) * a(2, 1) + (a(0, 1) * a(1, 3) - a(0, 3) * a(1, 1)) * a(2, 0)),

               -((a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0)) * a(3, 2) + (a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2)) * a(3, 1) + (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) * a(3, 0)),
               ((a(0, 0) * a(2, 1) - a(0, 1) * a(2, 0)) * a(3, 2) + (a(0, 2) * a(2, 0) - a(0, 0) * a(2, 2)) * a(3, 1) + (a(0, 1) * a(2, 2) - a(0, 2) * a(2, 1)) * a(3, 0)),
               -((a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * a(3, 2) + (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)) * a(3, 1) + (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * a(3, 0)),
               ((a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * a(2, 2) + (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)) * a(2, 1) + (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * a(2, 0))
           ) / det;
}

/// compute a pseudo inverse matrix
template <typename T, int NRows, int NCols>
inline
CMatrix<T, NCols, NRows> pseudo_inverse_of(const CMatrix<T, NRows, NCols>& mat) {
    auto transposed = transpose_of(mat);
    if (mat.num_cols() < mat.num_rows()) {
        return inverse_of(transposed * mat) * transposed;
    } else {
        return transposed * inverse_of(mat * transposed);
    }
}

// float/double 4x4 affine matrix
// float/double 3x3 or 4x4 rotation matrix
// float/double 3x1 or 4x1 translation vector

/// construct a rigid transformation
template <typename T>
inline
CMatrix<T, 4, 4> inverse_rigidity_of(const CMatrix<T, 4, 4>& a) {
    assert(a(3, 0) == 0 && a(3, 1) == 0 && a(3, 2) == 0 && a(3, 3) == 1);
    CMatrix<T, 4, 4> inv = make_matrix<T, 4, 4>(
                               a(0, 0), a(1, 0), a(2, 0), 0,
                               a(0, 1), a(1, 1), a(2, 1), 0,
                               a(0, 2), a(1, 2), a(2, 2), 0,
                               0, 0, 0, 1);
    auto Rt = AffineTransform(inv, CVector <T, 3> { a(0, 3), a(1, 3), a(2, 3)});
    inv(0, 3) = -Rt[0];
    inv(1, 3) = -Rt[1];
    inv(2, 3) = -Rt[2];
    return inv;
}

/// construct a similarity transformation
/// s*R*x+t -> 1/s*R'*(x-t) = 1/s*R'*x-1/s*R'*t
template <typename T>
inline
CMatrix<T, 4, 4> inverse_similarity_of(const CMatrix<T, 4, 4>& a) {
    assert(a(3, 0) == 0 && a(3, 1) == 0 && a(3, 2) == 0 && a(3, 3) == 1);
    float s2 = norm2_squared_of(make_vector_from_column(a, 0));
    if (s2 == 0) {
        ThrowRuntimeError("non-invertible transformation");
    }
    CMatrix<T, 4, 4> inv = make_matrix<T, 4, 4>(
                               a(0, 0) / s2, a(1, 0) / s2, a(2, 0) / s2, 0,
                               a(0, 1) / s2, a(1, 1) / s2, a(2, 1) / s2, 0,
                               a(0, 2) / s2, a(1, 2) / s2, a(2, 2) / s2, 0,
                               0, 0, 0, 1);
    auto Rt = AffineTransform(inv, make_vector(a(0, 3), a(1, 3), a(2, 3)));
    inv(0, 3) = -Rt[0];
    inv(1, 3) = -Rt[1];
    inv(2, 3) = -Rt[2];
    return inv;
}

/// construct a translation matrix
template <typename T>
inline
CMatrix<T, 4, 4> make_translation_matrix(T tx, T ty, T tz) {
    CMatrix<T, 4, 4> ret = make_diagonal_matrix(1, 1, 1, 1);
    ret(0, 3) = tx;
    ret(1, 3) = ty;
    ret(2, 3) = tz;
    return ret;
}

/// extract a translation vector
template <typename T, int NRows, int NCols>
inline
CVector<T, 3> make_translation_vector(const CMatrix<T, NRows, NCols>& mat) {
    assert(mat.num_rows() == 3 || mat.num_rows() == 4);
    assert(mat.num_cols() == 3 || mat.num_cols() == 4);
    return make_vector(mat(0, 3), mat(1, 3), mat(2, 3));
}

/// construct a translation matrix
template <typename T, int NRows, int NCols>
inline
CMatrix<T, 4, 4> make_translation_matrix(const CMatrix<T, NRows, NCols>& vec) {
    assert(vec.num_rows() == 3 || vec.num_rows() == 4);
    assert(vec.num_cols() == 1);
    return make_translation_matrix(vec[0], vec[1], vec[2]);
}

/// construct a translation matrix
template <typename T, int NRows>
inline
CMatrix<T, 4, 4> make_rotation_matrix(const CMatrix<T, NRows, NRows>& rot) {
    static_assert(NRows == 3 || NRows == 4, "invalid size");
    CMatrix<T, 4, 4> ret;
    ret.fill_with(0);
    for (int c = 0; c < 3; c++) {
        for (int r = 0; r < 3; r++) {
            ret(r, c) = rot(r, c);
        }
    }
    ret(3, 3) = 1;
    return ret;
}

/// construct a rotation matrix
/// @tparam T element
/// @tparam NRows number of rows
template <typename T , int NRows>
inline
void ConvertToEulerAngle(
    const CMatrix<T, NRows, NRows>& mat, //
    T& a, //
    T& b, //
    T& c //
) {
    static_assert(NRows == 3 || NRows == 4, "invalid size");
    // assume b>= 0
    if (abs(mat(2, 2)) != 1) {
        a = atan2(mat(1, 2), mat(0, 2));
        b = atan2(sqrt(mat(0, 2) * mat(0, 2) + mat(1, 2) * mat(1, 2)), mat(2, 2));
        c = atan2(mat(2, 1), -mat(2, 0));
    } else {
        // assume c = 0
        a = atan2(mat(1, 0), mat(1, 1));
        b = 0;
        c = 0;
    }
}

/// defompose a rotation matrix to roll-pitch-yaw angles
/// @tparam T element
/// @tparam NRows number of rows
/// @see http://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles
template <typename T, int NRows, int NCols>
inline
void ConvertToRollPitchYawAngle(
    const CMatrix<T, NRows, NCols>& mat, // rot(roll)*rot(pitch)*rot(yaw)
    T& roll , // angle
    T& pitch , // angle, in (-ƒÎ/2,ƒÎ/2)
    T& yaw  // angle
) {
    static_assert(NRows == 0 || NRows == 3 || NRows == 4, "invalid size");
    static_assert(NCols == 0 || NCols == 3 || NCols == 4, "invalid size");
    // assume pitch in (-ƒÎ/2,ƒÎ/2)
    T cp = sqrt(mat(0, 0) * mat(0, 0) + mat(1, 0) * mat(1, 0));
    if (cp) {
        roll = atan2(mat(1, 0), mat(0, 0)); // in (-ƒÎ, ƒÎ)
        pitch = atan2(-mat(2, 0), cp); // in (-ƒÎ/2, ƒÎ/2)
        yaw = atan2(mat(2, 1), mat(2, 2)); // in (-ƒÎ, ƒÎ)
    } else {
        // cp=0, sp=?}1
        // assume sy=0, cy=sp
        // [0   -sr*cy  cr  ]
        // [0   cr*cy   sr  ]
        // [-sp 0       0   ]
        roll = atan2(-mat(0, 1), mat(1, 1)); // (-ƒÎ, ƒÎ)
        pitch = -mat(2, 0) * M_PI / 2; // -ƒÎ/2 or ƒÎ/2
        yaw = 0; // assumption for disambiguation
    }
}

/// construct a rotation matrix from Tait-Bryan angles
/// @return rot(roll)*rot(pitch)*rot(yaw) or rot(z)*rot(y)*rot(x)
/// @tparam T element
/// @see http://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles
inline
CMatrix<float, 4, 4> make_rotation_matrix_from_roll_pitch_yaw(float roll, float pitch, float yaw) {
    float cr = cos(roll);
    float sr = sin(roll);
    float cp = cos(pitch);
    float sp = sin(pitch);
    float cy = cos(yaw);
    float sy = sin(yaw);
    return make_matrix<float, 4, 4>(
               cr * cp, cr * sp * sy - sr * cy, cr * sp * cy + sr * sy, 0,
               sr * cp, sr * sp * sy + cr * cy, sr * sp * cy - cr * sy, 0,
               -sp, cp * sy, cp * cy, 0,
               0, 0, 0, 1);
}

/// construct a rotation matrix
/// @tparam T element
inline
CMatrix<float, 4, 4> make_rotation_matrix_from_euler_angle(float a, float b, float c) {
    float ca = cos(a);
    float sa = sin(a);
    float cb = cos(b);
    float sb = sin(b);
    float cc = cos(c);
    float sc = sin(c);
    return make_matrix<float, 4, 4>(
               ca * cb * cc - sa * sc, -ca * cb * sc - sa * cc, ca * sb, 0,
               sa * cb * cc + ca * sc, -sa * cb * sc + ca * cc, sa * sb, 0,
               -sb * cc, sb * sc, cb, 0,
               0, 0, 0, 1);
}

/// construct a rotation matrix
/// @tparam T element
template <typename T>
inline
CMatrix<T, 4, 4> make_rotation_matrix(const CVector<T, 3>& axis, T cosangle, T sinangle) {
    if (dot(axis, axis) == 0) {
        return make_diagonal_matrix(1, 1, 1, 1);
    }

    auto u = normalized_of(axis);
    return make_matrix<float, 4, 4>(
               u[0] * u[0] + (1 - u[0] * u[0]) * cosangle,
               u[0] * u[1] * (1 - cosangle) - u[2] * sinangle,
               u[0] * u[2] * (1 - cosangle) + u[1] * sinangle,
               0,
               u[0] * u[1] * (1 - cosangle) + u[2] * sinangle,
               u[1] * u[1] + (1 - u[1] * u[1]) * cosangle,
               u[1] * u[2] * (1 - cosangle) - u[0] * sinangle,
               0,
               u[0] * u[2] * (1 - cosangle) - u[1] * sinangle,
               u[1] * u[2] * (1 - cosangle) + u[0] * sinangle,
               u[2] * u[2] + (1 - u[2] * u[2]) * cosangle,
               0,
               0, 0, 0, 1);
}

/// construct a rotation matrix
/// @tparam T element
template <typename T>
inline
CMatrix<T, 4, 4> make_rotation_matrix(const CVector<T, 3>& axis, T angle) {
    return make_rotation_matrix<T>(axis, cos(angle), sin(angle));
}

/// construct a rotation matrix
/// @tparam T element
template <typename T>
inline
CMatrix<T, 4, 4> make_rotation_matrix(const CVector<T, 3>& rotvec) {
    return make_rotation_matrix<T>(normalized_of(rotvec), norm2_of(rotvec));
}

/// construct a rotation matrix
/// @tparam T element
template <typename T>
inline
CMatrix<T, 4, 4> make_rotation_matrix(const CVector<T, 3>& src, const CVector<T, 3>& dst) {
    auto x = cross(normalized_of(src), normalized_of(dst));
    auto c = dot(normalized_of(src), normalized_of(dst));
    auto s = norm2_of(x);
    if (s) {
        return make_rotation_matrix(normalized_of(x), c, s);
    } else {
        return make_diagonal_matrix(1, 1, 1, 1);
    }
}

/// construct a rotation matrix
/// @tparam T element
inline
CMatrix<float, 4, 4> make_rotation_matrix_from_mouse_drag(float rx, float ry, float dx, float dy) {
    CVector<float, 3> v1 = normalized_of<float, 3, 1>(CVector <float, 3> {rx, ry, 1.f});
    CVector<float, 3> v2 = normalized_of<float, 3, 1>(CVector <float, 3> {rx + dx, ry + dy, 1.f});
    CVector<float, 3> cr = cross(v1, v2);
    return make_rotation_matrix(cr, dot(v1, v2), norm2_of(cr));
}

#if defined(MK_LBUTTON)
/// construct a rotation matrix
/// @tparam T element
inline
CMatrix<float, 4, 4> make_rigid_matrix_from_mouse_drag(float sx, float sy, float dx, float dy, int nFlags) {
    // left button
    if ((nFlags & MK_LBUTTON) &&
        !(nFlags & MK_MBUTTON) &&
        !(nFlags & MK_RBUTTON)) {
        return make_translation_matrix<float>(dx, -dy, 0);
    }
    // middle button
    else if ((!(nFlags & MK_LBUTTON) &&
              (nFlags & MK_MBUTTON) &&
              !(nFlags & MK_RBUTTON)) ||
             ((nFlags & MK_LBUTTON) &&
              !(nFlags & MK_MBUTTON) &&
              (nFlags & MK_RBUTTON))) {
        return make_translation_matrix<float>(0, 0, dy);
    }
    // right button
    else if (!(nFlags & MK_LBUTTON) &&
             !(nFlags & MK_MBUTTON) &&
             (nFlags & MK_RBUTTON)) {
        return make_rotation_matrix_from_mouse_drag(sx, -sy, dx, -dy);
    } else {
        ThrowLogicError("undefined mouse button");
    }
}
#endif

/// construct a rigid transformation matrix
/// @tparam T element
/// @tparam nMatRows number of rows
/// @tparam nVecRows number of rows
template <typename T, int MatDim, int VecDim>
inline
CMatrix<T, 4, 4> make_affine_matrix(const CMatrix<T, MatDim, MatDim>& rot, const CMatrix<T, VecDim, 1>& trans) {
    assert(rot.num_rows() == rot.num_cols());
    assert(rot.num_rows() == 3 || rot.num_rows() == 4);
    return {
        rot(0, 0), rot(1, 0), rot(2, 0), 0,
        rot(0, 1), rot(1, 1), rot(2, 1), 0,
        rot(0, 2), rot(1, 2), rot(2, 2), 0,
        trans(0, 0), trans(1, 0), trans(2, 0), 1,
    };
}

/// affine transform a vector
/// @tparam T1 element
/// @tparam T2 element
/// @return 3x1 vector
template <typename T1, int MatRows, int MatCols, typename T2, int VecRows>
inline
CMatrix<decltype(T1()*T2()), 3, 1> AffineTransform(
    const CMatrix<T1, MatRows, MatCols>& mat, // 3x4 or 4x4
    const CMatrix<T2, VecRows, 1>& vec
) {
    assert(mat.num_rows() == 3 || mat.num_rows() == 4);
    assert(mat.num_cols() == 4);
    assert(vec.num_rows() == 3);
    decltype(T1()*T2()) z = 1;
    if (mat.num_rows() == 4) {
        z = mat(3, 0) * vec[0] + mat(3, 1) * vec[1] + mat(3, 2) * vec[2] + mat(3, 3);
    }
    return CVector <decltype(T1()*T2()), 3> {
        (mat(0, 0) * vec[0] + mat(0, 1) * vec[1] + mat(0, 2) * vec[2] + mat(0, 3)) / z,
        (mat(1, 0) * vec[0] + mat(1, 1) * vec[1] + mat(1, 2) * vec[2] + mat(1, 3)) / z,
        (mat(2, 0) * vec[0] + mat(2, 1) * vec[1] + mat(2, 2) * vec[2] + mat(2, 3)) / z
    };
}

/// rotate a vector
/// @tparam T element
/// @return 3x1 vector
template <typename T>
inline
CVector<T, 3> RotateVector(
    const CMatrix<T, 4, 4>& mat, // 3x3, 3x4, or 4x4 matrix; the 4th column is ignored if 3x4 or 4x4
    const CVector<T, 3>& vec // 3x1 vector
) {
    return {
        mat(0, 0) *vec[0] + mat(0, 1) *vec[1] + mat(0, 2) *vec[2],
        mat(1, 0) *vec[0] + mat(1, 1) *vec[1] + mat(1, 2) *vec[2],
        mat(2, 0) *vec[0] + mat(2, 1) *vec[1] + mat(2, 2) *vec[2]
    };
}

template <typename T>
inline
void DecomposeProjection(const CMatrix<T, 3, 4>& projection,
                         CMatrix<T, 3, 4>& intrinsic,
                         CMatrix<T, 4, 4>& extrinsic) {
    CVector<T, 3> p1 = make_vector(projection(0, 0), projection(0, 1), projection(0, 2));
    CVector<T, 3> p2 = make_vector(projection(1, 0), projection(1, 1), projection(1, 2));
    CVector<T, 3> p3 = make_vector(projection(2, 0), projection(2, 1), projection(2, 2));
    T p03 = projection(0, 3);
    T p13 = projection(1, 3);
    T p23 = projection(2, 3);

    // intrinsics
    T theta = acos(-dot(cross(p1, p3), cross(p2, p3)) / (norm2_of(cross(p1, p3)) * norm2_of(cross(p2, p3))));
    T au = norm2_of(cross(p1, p3)) * sin(theta);
    T av = norm2_of(cross(p2, p3)) * sin(theta);
    T u0 = dot(p1, p3);
    T v0 = dot(p2, p3);

    // extrinsics
    CVector<T, 3> r1 = (p1 + (p2 - p3 * v0) * au / av * cos(theta) - p3 * u0) / au;
    CVector<T, 3> r2 = (p2 - p3 * v0) * sin(theta) / au;
    CVector<T, 3> r3 = p3;
    CVector<T, 3> t = make_vector(p03 / au + (p13 - v0 * p23) / av * cos(theta) - u0 * p23 / au,
                                  sin(theta) / av * (p13 - v0 * p23),
                                  p23);
    intrinsic = {
        au, 0, 0,
        -au / tan(theta), av / sin(theta), 0,
        u0, v0, 1,
        0, 0, 0
    };
    extrinsic = {
        r1[0], r2[0], r3[0], 0,
        r1[1], r2[1], r3[1], 0,
        r1[2], r2[2], r3[2], 0,
        t[0], t[1], t[2], 1
    };
//#ifdef _DEBUG
//    CMatrix<T, 3, 4> error = intrinsic * extrinsic  - projection;
//    std::clog << "error=";
//    Dump(error);
//#endif
}

} // namespace slib
