// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <cassert>

#include <algorithm>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <vector>
#include <sstream>
#include <type_traits>

#include "exception_util.h"

namespace slib {

template <typename T, int NRows, int NCols>
class CMatrix;

/// for operator*=()
template <typename T, int NRows, int NCols, typename S>
void ScaleMatrixInPlace(CMatrix<T, NRows, NCols>& m, S s);

/// for operator+=() and operator-=()
template <typename T, int NRows, int NCols, typename T2, int NRows2, int NCols2, typename S>
void ScaleAddMatrixInPlace(CMatrix<T, NRows, NCols>& m1, S s, const CMatrix<T2, NRows2, NCols2>& m2);

/// for operator*()
template <typename T1, int NRows1, int NCols1, typename T2, int NRows2, int NCols2, typename T3, int NRows3, int NCols3>
void MultiplyMatrix(CMatrix<T1, NRows1, NCols1>& m1, const CMatrix<T2, NRows2, NCols2>& m2, const CMatrix<T3, NRows3, NCols3>& m3);

/// matrix
///
/// set NRows and NCols to non-zero for a fixed-size matrix
template <typename T, int NRows = 0, int NCols = 0>
class CMatrix {
public:
    CMatrix() = default;
    CMatrix(const CMatrix&) = default;
    CMatrix& operator=(const CMatrix&) = default;

    CMatrix(CMatrix&& mat) : elements_(std::move(mat.elements_)) {}

    CMatrix& operator=(CMatrix && mat) {
        elements_ = std::move(mat.elements_);
        return *this;
    }

    template <typename T2, int NRows2, int NCols2>
    CMatrix(const CMatrix<T2, NRows2, NCols2>& mat) {
        resize(mat.num_rows(), mat.num_cols());
        copy_from(mat.data());
    }

    template <typename T2, int NRows2, int NCols2>
    CMatrix& operator=(const CMatrix<T2, NRows2, NCols2>& mat) {
        resize(mat.num_rows(), mat.num_cols());
        copy_from(mat.data());
        return *this;
    }

    CMatrix(std::initializer_list<T> list) {
        assert(NRows * NCols == list.size());
        std::copy(list.begin(), list.end(), elements_.data());
    }

    /// constructor with size and elements.
    explicit CMatrix(int nrows, int ncols = 1) {
        resize(nrows, ncols);
    }

    explicit CMatrix(const T *p, int nrows = NRows, int ncols = NCols) {
        resize(nrows, ncols);
        assert(nrows * ncols);
        std::copy_n(p, nrows * ncols, elements_.data());
    }

    /// @return an element
    T& operator()(int r, int c) {
        assert(0 <= r && r < num_rows() && 0 <= c && c < num_cols());
        return elements_.data()[r + num_rows() * c];
    }

    /// @return an element
    T const& operator()(int r, int c) const {
        assert(0 <= r && r < num_rows() && 0 <= c && c < num_cols());
        return elements_.data()[r + num_rows() * c];
    }

    /// @return an element
    T& operator[](int i) {
        assert(0 <= i && i < num_rows() * num_cols());
        return elements_.data()[i];
    }

    /// @return an element
    T const& operator[](int i) const {
        assert(0 <= i && i < num_rows() * num_cols());
        return elements_.data()[i];
    }

    /// addition assignment
    template <typename T2, int NRows2, int NCols2>
    CMatrix& operator+=(const CMatrix<T2, NRows2, NCols2>& mat) {
        ScaleAddMatrixInPlace(*this, 1, mat);
        return*this;
    }

    /// subtraction assignment
    template <typename T2, int NRows2, int NCols2>
    CMatrix& operator-=(const CMatrix<T2, NRows2, NCols2>& mat) {
        ScaleAddMatrixInPlace(*this, -1, mat);
        return*this;
    }

    /// multiplication assignment.
    /// size may be changed for variable-size matrix.
    template <typename T2, int NRows2, int NCols2>
    CMatrix& operator*=(const CMatrix<T2, NRows2, NCols2>& mat) {
        return operator=(*this * mat);
    }

    /// multiplication assignment
    template <typename S, typename /*Reserved*/ = typename std::enable_if<std::is_arithmetic<S>::value>::type>
    CMatrix& operator*=(S scalar) {
        ScaleMatrixInPlace(*this, scalar);
        return*this;
    }

    /// division assignment
    template <typename S>
    CMatrix& operator/=(S scalar) {
        if (scalar == 0) {
            ThrowRuntimeError("division by zero");
        }
        return *this *= 1.0 / scalar; // divide by float or double
    }

    /// unary minus
    CMatrix < decltype(T() * -1), NRows, NCols > operator-() const {
        return CMatrix < decltype(T() * -1), NRows, NCols > (*this) *= -1;
    }

    /// addition
    template <typename T2, int NRows2, int NCols2>
    CMatrix < decltype(T() + T2()), NRows, NCols > operator+(const CMatrix<T2, NRows2, NCols2>& mat) const {
        return CMatrix < decltype(T() + T2()), NRows, NCols > (*this) += mat;
    }

    /// subtraction
    template <typename T2, int NRows2, int NCols2>
    CMatrix < decltype(T() - T2()), NRows, NCols > operator-(const CMatrix<T2, NRows2, NCols2>& mat) const {
        return CMatrix < decltype(T() - T2()), NRows, NCols > (*this) -= mat;
    }

    /// multiplication
    template <typename T2, int NRows2, int NCols2>
    CMatrix<decltype(T() * T2()), NRows, NCols2> operator*(const CMatrix<T2, NRows2, NCols2>& mat) const {
        CMatrix<decltype(T() * T2()), NRows, NCols2> ret(num_rows(), mat.num_cols());
        MultiplyMatrix(ret, *this, mat);
        return ret;
    }

    /// multiplication
    template <typename S, typename /*Reserved*/ = typename std::enable_if<std::is_arithmetic<S>::value>::type>
    CMatrix<decltype(S()*T()), NRows, NCols> operator*(S scalar) const {
        return CMatrix<decltype(S()*T()), NRows, NCols>(*this) *= scalar;
    }

    /// division
    template <typename S, typename /*Reserved*/ = typename std::enable_if<std::is_arithmetic<S>::value>::type>
    CMatrix < decltype(T() / S()), NRows, NCols > operator/(S scalar) const {
        return CMatrix < decltype(T() / S()), NRows, NCols > (*this) /= scalar;
    }

    /// multiplication
    template <typename S, typename /*Reserved*/ = typename std::enable_if<std::is_arithmetic<S>::value>::type>
    friend CMatrix<decltype(S()*T()), NRows, NCols> operator*(S scalar, const CMatrix& mat) {
        return mat * scalar;
    }

    /// equality
    template <typename T2, int NRows2, int NCols2>
    bool operator==(const CMatrix<T2, NRows2, NCols2>& mat) const {
        if (num_rows() != mat.num_rows() || num_cols() != mat.num_cols()) {
            return false;
        }
        return std::equal(elements_.data(), elements_.data() + num_rows() * num_cols(), mat.data());
    }

    /// inequality
    template <typename T2, int NRows2, int NCols2>
    bool operator!=(const CMatrix<T2, NRows2, NCols2>& mat2) const {
        return !(*this == mat2);
    }

    bool empty() const {
        return num_rows()*num_cols()==0;
    }

    /// copy elements from a pointer to elements
    template <typename T2>
    CMatrix& copy_from(const T2 *elements) {
        assert(elements);
        std::copy_n(elements, num_rows() * num_cols(), elements_.data());
        return *this;
    }

    /// change the size of matrix. the elements are not initialized.
    void resize(int nrows, int ncols = 1) {
        return elements_.resize(nrows, ncols);
    }

    /// change the size of matrix. the elements are not initialized.
    void reshape(int nrows, int ncols) {
        return elements_.reshape(nrows, ncols);
    }

    /// initialize elements with a value
    void fill_with(T val) {
        std::fill_n(elements_.data(), num_rows()*num_cols(), val);
    }

    void append_rows(int n) {
        elements_.append_rows(n);
    }

    template <typename T2, int NR2, int NC2>
    void append_rows(const CMatrix<T2, NR2, NC2>& rows) {
        assert(num_cols() == rows.num_cols());
        int r = num_rows();
        elements_.append_rows(rows.num_rows());
        for (int c = 0; c < rows.num_cols(); c++) {
            std::copy_n(rows.col_data(c), rows.num_rows(), &operator()(r, c));
        }
    }

    void append_cols(int n) {
        elements_.append_cols(n);
    }

    template <typename T2, int NR2, int NC2>
    void append_cols(const CMatrix<T2, NR2, NC2>& cols) {
        assert(num_rows() == cols.num_rows());
        int c = num_cols();
        elements_.append_cols(cols.num_cols());
        std::copy_n(cols.data(), cols.num_rows()*cols.num_cols(), col_data(c));
    }

    /// @return a pointer to elements
    const T *data() const {
        return elements_.data();
    }

    /// @return a pointer to elements
    T *data() {
        return elements_.data();
    }

    /// @return a pointer to elements
    const T *col_data(int c) const {
        return elements_.data() + c * num_rows();
    }

    /// @return a pointer to elements
    T *col_data(int c) {
        return elements_.data() + c * num_rows();
    }

    /// @return the number of rows
    int num_rows() const {
        return elements_.num_rows();
    }

    /// @return the number of columns
    int num_cols() const {
        return elements_.num_cols();
    }

private:
    class VariableArray {
    public:
        VariableArray() : data_(0), nrows_(NRows), ncols_(NCols) {}
        VariableArray(const VariableArray& v) {
            assert(!NRows || v.nrows_ == NRows);
            assert(!NCols || v.ncols_ == NCols);
            data_ = new T[v.nrows_ * v.ncols_];
            nrows_ = v.nrows_;
            ncols_ = v.ncols_;
            std::copy_n(v.data_, nrows_ * ncols_, data_);
        }
        VariableArray& operator=(const VariableArray& v) {
            assert(!NRows || v.nrows_ == NRows);
            assert(!NCols || v.ncols_ == NCols);
            if (nrows_ * ncols_ < v.nrows_ * v.ncols_) {
                delete[] data_;
                data_ = new T[v.nrows_ * v.ncols_];
            }
            nrows_ = v.nrows_;
            ncols_ = v.ncols_;
            std::copy_n(v.data_, nrows_ * ncols_, data_);
            return *this;
        }
        VariableArray(VariableArray&& v) : data_(v.data_), nrows_(v.nrows_), ncols_(v.ncols_) {
            assert(!NRows || v.nrows_ == NRows);
            assert(!NCols || v.ncols_ == NCols);
            v.data_ = 0;
            v.nrows_ = 0;
            v.ncols_ = 0;
        }
        VariableArray& operator=(VariableArray && v) {
            assert(!NRows || v.nrows_ == NRows);
            assert(!NCols || v.ncols_ == NCols);
            data_ = v.data_;
            nrows_ = v.nrows_;
            ncols_ = v.ncols_;
            v.data_ = 0;
            v.nrows_ = 0;
            v.ncols_ = 0;
            return *this;
        }
        ~VariableArray() {
            delete[]data_;
        }
        void resize(int nrows, int ncols = 1) {
            assert(!NRows || nrows == NRows);
            assert(!NCols || ncols == NCols);
            int nr = std::min(nrows_, nrows);
            int nc = std::min(ncols_, ncols);
            if (nrows_ * ncols_ < nrows * ncols) {
                // expand
                T *backup = data_;
                data_ = new T[nrows * ncols];
                for (int c = 0; c < nc; c++) {
                    for (int r = 0; r < nr; r++) {
                        data_[r + nrows * c] = backup[r + nrows_ * c];
                    }
                }
                delete[] backup;
            } else {
                // shrink
                if (nrows < nrows_) {
                    // [04]    [02]4
                    // [15] -> [13]5
                    //  26
                    //  37
                    for (int c = 1; c < nc; c++) {
                        for (int r = 0; r < nr; r++) {
                            data_[r + nrows * c] = data_[r + nrows_ * c];
                        }
                    }
                } else if (nrows > nrows_) {
                    // [02]46    [03]
                    // [13]57 -> [14]
                    //            25
                    for (int c = nc - 1; c > 0; c--) {
                        for (int r = nr - 1; r >= 0; r--) {
                            data_[r + nrows * c] = data_[r + nrows_ * c];
                        }
                    }
                }
            }
            nrows_ = nrows;
            ncols_ = ncols;
        }
        void reshape(int nrows, int ncols) {
            assert(!NRows || NRows == nrows);
            assert(!NCols || NCols == ncols);
            assert(nrows_ * ncols_ == nrows * ncols);
            nrows_ = nrows;
            ncols_ = ncols;
        }
        void append_rows(int n) {
            assert(!NRows);
            int old_nrows = nrows_;
            resize(nrows_ + n, ncols_);
            for (int c = 0; c < ncols_; c++) {
                for (int r = old_nrows; r < nrows_ ; r++) {
                    data_[r + c * nrows_ ] = 0;
                }
            }
        }
        void append_cols(int n) {
            assert(!NCols);
            int old_ncols = ncols_;
            resize(nrows_, ncols_ + n);
            std::fill_n(data_ + nrows_ * old_ncols, nrows_ * n, 0);
        }
        int num_rows() const {
            return nrows_;
        }
        int num_cols() const {
            return ncols_;
        }
        T *data() {
            return data_;
        }
        const T *data() const {
            return data_;
        }
    private:
        T *data_;
        int nrows_;
        int ncols_;
    };

    class FixedArray {
    public:
        FixedArray() = default;
        FixedArray(const FixedArray&) = default;
        FixedArray& operator=(const FixedArray&) = default;
        void resize(int nrows, int ncols = 1) {
            assert(NRows == nrows && NCols == ncols);
        }
        void reshape(int nrows, int ncols) {
            assert(NRows == nrows && NCols == ncols);
        }
        void append_rows(int n) {
            assert(n == 0);
        }
        void append_cols(int n) {
            assert(n == 0);
        }
        int num_rows() const {
            return NRows;
        }
        int num_cols() const {
            return NCols;
        }
        T *data() {
            return data_;
        }
        const T *data() const {
            return data_;
        }
    private:
        T data_[NRows *NCols];
    };

#if 0
    template <typename T, int NRows, int NCols>
    class RefArray {
    public:
        RefArray() = delete;
        RefArray(const RefArray&) = default;
        RefArray& operator=(const RefArray&) = default;
        RefArray(T *data, int dimension) : data_(data), dimension_(dimension) {}
        void resize(int nrows, int ncols = 1) {
            assert(NRows == nrows && NCols == ncols);
        }
        void reshape(int nrows, int ncols) {
            assert(NRows == nrows && NCols == ncols);
        }
        void append_rows(int n) {
            assert(n == 0);
        }
        void append_cols(int n) {
            assert(n == 0);
        }
        int num_rows() const {
            return NRows;
        }
        int num_cols() const {
            return NCols;
        }
        T *data() {
            ThrowLogicError("invalid operation");
        }
        const T *data() const {
            ThrowLogicError("invalid operation");
        }
    private:
        T *const data_;
        int dimension_; // leading dimension
    };
#endif

    typename std::conditional < NRows *NCols != 0, FixedArray , VariableArray >::type elements_; // storage for fixed-size or variable-size elements
};

/// import from a file
template <typename T, int NRows, int NCols>
inline
void Read(CMatrix<T, NRows, NCols>& mat, const std::string& filename) {
    std::clog << "matrix <= " << filename << std::endl;
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        ThrowRuntimeError("failed to open file: " + filename);
    }
    // read the first line to determine with number of cols
    std::string line;
    getline(ifs, line);
    std::istringstream str(line);
    std::vector<T> elements;
    int ncols = 0;
    while (1) {
        T e;
        str >> e;
        if (str.fail()) {
            break;
        }
        ncols++;
        elements.push_back(e);
    }

    while (1) {
        T e;
        ifs >> e;
        if (ifs.fail()) {
            break;
        }
        elements.push_back(e);
    }
    int nrows = elements.size() / ncols;
    mat.resize(nrows, ncols);
    for (int r = 0; r < nrows; r++) {
        for (int c = 0; c < ncols; c++) {
            mat(r, c) = elements[r * ncols + c];
        }
    }
}

/// export to a file
template <typename T, int NRows, int NCols>
inline
void Write(const CMatrix<T, NRows, NCols>& mat, const std::string& filename) {
    std::clog << "matrix => " << filename << std::endl;
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        ThrowRuntimeError("failed to open file: " + filename);
    }
    ofs << std::scientific;
    for (int r = 0; r < mat.num_rows(); r++) {
        for (int c = 0; c < mat.num_cols(); c++) {
            ofs << mat(r, c) << " ";
        }
        ofs << std::endl;
    }
}

#define DUMP(x) { std::clog << #x << ":" << std::endl; Dump(x, true); }

/// print elements (for debug)
template <typename T, int NRows, int NCols>
inline void Dump(
    const CMatrix<T, NRows, NCols>& mat,
    bool force = false, // dump even if the matrix is huge
    std::ostream& ost = std::cout  // stream
) {
    const int max_size = 12;
    ost << "matrix = " << mat.num_rows() << "x" << mat.num_cols() << std::endl;
    if (force || (mat.num_rows() < max_size && mat.num_cols() < max_size)) {
        ost << std::scientific;
        for (int r = 0; r < mat.num_rows(); r++) {
            for (int c = 0; c < mat.num_cols(); c++) {
                ost << (float)mat(r, c) << " ";
            }
            ost << std::endl;
        }
    }
}

/// vector
template <typename T, int NRows = 0>
using CVector = CMatrix<T, NRows, 1>;

/*
/// construct a vector
template <typename T, int NRows, typename T2>
inline
void make_vector_impl(CVector<T, NRows>& vec, int i, T2 e) {
    assert(i == NRows - 1);
    vec[i] = e;
}
template <typename T, int NRows, typename T2, typename...Ts>
inline
void make_vector_impl(CVector<T, NRows>& vec, int i, T2 e, Ts...rest) {
    vec[i] = e;
    make_vector_impl(vec, i + 1, rest...);
}
template <typename T, typename...Ts>
inline
CVector < T, sizeof...(Ts) + 1 > make_vector(T e, Ts...rest) {
    CVector < T, sizeof...(Ts) + 1 > vec;
    vec[0] = e;
    make_vector_impl(vec, 1, rest...);
    return vec;
}
*/
/// construct a matrix
template<typename T, int NRows, int NCols, typename T2>
inline
void impl_make_matrix(CMatrix<T, NRows, NCols>& mat, int i, T2 e) {
    mat(NRows - 1, NCols - 1) = T(e);
}
template<typename T, int NRows, int NCols, typename T2, typename...Ts>
inline
void impl_make_matrix(CMatrix<T, NRows, NCols>& mat, int i, T2 e, Ts...rest) {
    int r = i / NCols;
    int c = i % NCols;
    mat(r, c) = T(e);
    impl_make_matrix(mat, i + 1, rest...);
}
template<typename T, int NRows, int NCols, typename...Ts>
inline
CMatrix<T, NRows, NCols> make_matrix(Ts...elements) {
    static_assert(NRows * NCols == sizeof...(Ts), "invalid arguments");
    CMatrix<T, NRows, NCols> mat;
    impl_make_matrix(mat, 0, elements...);
    return mat;
}

/// construct a diagonal matrix
template <typename T, int NDims>
inline
CMatrix<T, NDims, NDims> make_diagonal_matrix(const CVector<T, NDims>& v) {
    CMatrix<T, NDims, NDims> ret(v.num_rows(), v.num_rows());
    ret.fill_with(0);
    for (int i = 0; i < v.num_rows(); i++) {
        ret(i, i) = v(i, 0);
    }
    return ret;
}

template<typename T, int NDims, typename T2>
inline
void impl_make_diagonal_matrix(CMatrix<T, NDims, NDims>& mat, int i, T2 e) {
    assert(i == NDims - 1);
    mat(i, i) = e;
}
template<typename T, int NDims, typename T2, typename ...Ts>
inline
void impl_make_diagonal_matrix(CMatrix<T, NDims, NDims>& mat, int i, T2 e, Ts...rest) {
    mat(i, i) = e;
    impl_make_diagonal_matrix(mat, i + 1, rest...);
}
template<typename T, typename ...Ts>
inline
CMatrix < T, sizeof...(Ts) + 1, sizeof...(Ts) + 1 > make_diagonal_matrix(T e00, Ts...rest) {
    CMatrix < T, sizeof...(Ts) + 1, sizeof...(Ts) + 1 > mat;
    mat.fill_with(0);
    mat(0, 0) = e00;
    impl_make_diagonal_matrix(mat, 1, rest...);
    return mat;
}

template<typename T, int NRows, int NCols>
inline
void impl_make_matrix_from_column_vectors(CMatrix<T, NRows, NCols>& mat, int c, const CVector<T, NRows>& v) {
    assert(c == NCols - 1);
    std::copy_n(v.data(), NRows, &mat(0, c));
}
template<typename T, int NRows, int NCols, typename ...Ts>
inline
void impl_make_matrix_from_column_vectors(CMatrix<T, NRows, NCols>& mat, int c, const CVector<T, NRows>& v, Ts...rest) {
    std::copy_n(v.data(), NRows, &mat(0, c));
    impl_make_matrix_from_column_vectors(mat, c + 1, rest...);
}
template<typename T, int NRows, typename ...Ts>
inline
CMatrix < T, NRows, sizeof...(Ts) + 1 > make_matrix_from_column_vectors(const CVector<T, NRows>& v, Ts...rest) {
    CMatrix < T, NRows, sizeof...(Ts) + 1 > mat;
    mat.fill_with(0);
    std::copy_n(v.data(), NRows, &mat(0, 0));
    impl_make_matrix_from_column_vectors(mat, 1, rest...);
    return mat;
}
/*
/// construct a matrix from a set of row vectors
template <typename T, int NRows, int NCols>
inline
CMatrix<T, NRows, NCols> make_matrix_from_column_vectors(std::initializer_list<CVector<T, NRows>> list) {
auto it = list.begin();
CMatrix<T, NRows, NCols> mat(it->num_rows(), list.size());
for (int c = 0; c < mat.num_cols(); c++, ++it) {
for (int r = 0; r < mat.num_rows(); r++) {
mat(r, c) = (*it)[r];
}
}
return mat;
}
*/

/// construct a matrix from a set of row vectors
template <typename T, int NRows, int NCols>
inline
CMatrix<T> make_matrix_from_matrix(const CMatrix<T, NRows, NCols>& mat, int sr, int sc, int nr, int nc) {
    assert(sr >= 0 && sr + nr <= mat.num_rows());
    assert(sc >= 0 && sc + nc <= mat.num_cols());
    CMatrix<T> ret(nr, nc);
    for (int c = 0; c < nc; c++) {
        for (int r = 0; r < nr; r++) {
            ret(r, c) = mat(sr + r, sc + c);
        }
    }
    return ret;
}

/// construct a row vector
template <typename T, int NRows, int NCols>
inline
CVector<T, NCols> make_vector_from_row(const CMatrix<T, NRows, NCols>& mat, int row) {
    CVector<T, NCols> vec(mat.num_cols(), 1);
    for (int c = 0; c < mat.num_cols(); c++) {
        vec[c] = mat(row, c);
    }
    return vec;
}

/// construct a column vector
template <typename T, int NRows, int NCols>
inline
CVector<T, NRows> make_vector_from_column(const CMatrix<T, NRows, NCols>& mat, int col) {
    CVector<T, NRows> vec(mat.num_rows(), 1);
    vec.copy_from(mat.data() + col * mat.num_rows());
    return vec;
}

/// scale a matrix.
/// @return m*=s
template <typename T, int NRows, int NCols, typename S>
inline
void ScaleMatrixInPlace(CMatrix<T, NRows, NCols>& m, S s) {
    static_assert(std::is_arithmetic<S>::value, "invalid scalar");
    assert(m.num_rows());
    std::transform(
        m.data(),
        m.data() + m.num_rows()*m.num_cols(),
        m.data(),
    [s](T e) {
        return e * s;
    });
}

/// scale and add a matrix.
/// @return m1+=s*m2
template <typename T, int NRows, int NCols, typename T2, int NRows2, int NCols2, typename S>
inline
void ScaleAddMatrixInPlace(CMatrix<T, NRows, NCols>& m1, S s, const CMatrix<T2, NRows2, NCols2>& m2) {
    assert(m1.num_rows() == m2.num_rows());
    assert(m1.num_cols() == m2.num_cols());
    int num = m1.num_rows() * m1.num_cols();
    for (int i = 0; i < num; i++) {
        m1[i] += s * m2[i];
    }
}

/// multiply matrices.
/// @return m1=m2*m3
template <typename T1, int NRows1, int NCols1, typename T2, int NRows2, int NCols2, typename T3, int NRows3, int NCols3>
inline
void MultiplyMatrix(CMatrix<T1, NRows1, NCols1>& m1, const CMatrix<T2, NRows2, NCols2>& m2, const CMatrix<T3, NRows3, NCols3>& m3) {
    assert(m2.num_rows() == m1.num_rows());
    assert(m2.num_cols() == m3.num_rows());
    assert(m3.num_cols() == m1.num_cols());
    assert((void *)&m2 != (void *)&m1);
    assert((void *)&m3 != (void *)&m1);
    for (int c = 0; c < m1.num_cols(); c++) {
        for (int r = 0; r < m1.num_rows(); r++) {
            T1 v = T1(0);
            for (int i = 0; i < m2.num_cols(); i++) {
                v += m2(r, i) * m3(i, c);
            }
            m1(r, c) = v;
        }
    }
}

} // namespace slib
