// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <mkl.h>
#include <mkl_spblas.h>
#include <vector>
#include <tuple>
#include <algorithm>
#include <string>
#include "MatrixBase.h"
#include "exception_util.h"

namespace slib {

template <typename T>
class MatrixGenerator;

/// sparse matrix in CSR, 3-array variation, one-based indexing.
/// @see http://software.intel.com/sites/products/documentation/hpc/mkl/mklman/GUID-9FCEB1C4-670D-4738-81D2-F378013412B0.htm
template <typename T>
class CSparseMatrix {
public:
    CSparseMatrix() : ncols_(0)  {
        row_index_.push_back(1); // one-based
    }

    CSparseMatrix(const CSparseMatrix& mat) = default;

    CSparseMatrix& operator=(const CSparseMatrix& mat) = default;

    CSparseMatrix(CSparseMatrix&& mat) :
        element_(std::move(mat.element_)),
        row_index_(std::move(mat.row_index_)),
        col_index_(std::move(mat.col_index_)),
        ncols_(mat.ncols_) {}

    CSparseMatrix& operator=(CSparseMatrix && mat) {
        element_ = std::move(mat.element_);
        row_index_ = std::move(mat.row_index_);
        col_index_ = std::move(mat.col_index_);
        ncols_ = mat.ncols_;
        return *this;
    }

    void Clear() {
        element_.clear();
        row_index_.assign(1, 1); // one-based
        col_index_.clear();
        ncols_ = 0;
    }

    void Resize(int nrows, int ncols) {
        if (row_index_.size() > nrows + 1) {
            // remove some rows
            row_index_.resize(nrows + 1);
            col_index_.erase(col_index_.begin() + row_index_.back(), col_index_.end());
            element_.erase(element_.begin() + row_index_.back(), element_.end());
        } else {
            // append empty rows if necessary
            while (row_index_.size() < nrows + 1) {
                row_index_.push_back(row_index_.back());
            }
        }
        if (ncols_ <= ncols) {
            // append empty cols
            ncols_ = ncols;
        } else {
            for (int r = 0; r < num_rows(); r++) {
                for (int i = row_index_[r ] - 1; i < row_index_[r + 1] - 1;) {
                    if (col_index_[i] >= ncols) {
                        col_index_.erase(col_index_.begin() + i);
                        element_.erase(element_.begin() + i);
                        for (int ir = r + 1; ir < row_index_.size(); ir++) {
                            row_index_[ir]--;
                        }
                    } else {
                        i++;
                    }
                }
            }
            ncols_ = ncols;
        }
    }

    void AppendRows(int nrows) {
        assert(nrows > 0);
        for (int i = 0; i < nrows; i++) {
            row_index_.push_back(row_index_.back());
        }
    }

    void AppendRows(T scale, const CSparseMatrix& rows) {
        int diff = row_index_.back() - 1;
        for (auto it = rows.row_index_.begin() + 1, end = rows.row_index_.end(); it != end; ++it) {
            row_index_.push_back(*it + diff);
        }
        if (scale) {
            col_index_.insert(col_index_.end(), rows.col_index_.begin(), rows.col_index_.end());
            for (auto it = rows.element_.begin(), end = rows.element_.end(); it != end; ++it) {
                element_.push_back(*it * scale);
            }
        }
    }

    void AppendColumns(int ncols) {
        assert(ncols > 0);
        ncols_ += ncols;
    }

    void AppendColumns(T scale, const CSparseMatrix& cols) {
        int oldncols = ncols_;
        ncols_ += cols.num_cols();
        if (scale) {
            for (int r = 0, nrows = num_rows(); r < nrows; r++) {
                for (int idx = cols.row_index_[r] - 1, end = cols.row_index_[r + 1] - 1; idx < end; idx++) {
                    // TODO: insert multiple elements at a time
                    int c = oldncols + cols.col_index_[idx] - 1;
                    T e = cols.element_[idx];
                    Add(r, c, e);
                }
            }
        }
    }

    void Add(int r, int c, T v) {
        assert(r >= 0 && r < num_rows());
        assert(c >= 0 && c < num_cols());
        if (!v) {
            return;
        }

        int idx = row_index_[r] - 1;
        for (int end = row_index_[r + 1] - 1; idx < end; idx++) {
            int col = col_index_[idx] - 1;
            if (c == col) {
                // add to an existing element
                element_[idx] += v;
                return;
            } else if (c < col) {
                // insert a new element
                std::clog << "warning: inserting element to sparse matrix" << std::endl;
                break;
            }
        }
        // insert or append a new element
        element_.insert(element_.begin() + idx, v);
        col_index_.insert(col_index_.begin() + idx, c + 1);
        for (r++; r < row_index_.size(); r++) {
            row_index_[r]++;
        }
    }

    /// count the number of non-zero elements
    int num_nonzeros() const {
        return row_index_.back() - 1;
    }

    /// return the number of rows
    int num_rows() const {
        return row_index_.size() - 1;
    }

    /// return the number of columns
    int num_cols() const {
        return ncols_;
    }

    /// return an element
    T operator()(int r, int c) const {
        assert(r >= 0 && r < num_rows());
        assert(c >= 0 && c < num_cols());
        for (int idx = row_index_[r] - 1, end = row_index_[r + 1] - 1; idx < end; idx++) {
            int col = col_index_[idx] - 1;
            if (col < c) {
                continue;
            } else if (col == c) {
                return element_[idx];
            } else {
                return 0;
            }
        }
        return 0;
    }

    /// return the pointer to non-zero elements
    const T *element_ptr() const {
        if (element_.size()) {
            return element_.data();
        } else {
            return 0;
        }
    }

    /// return the pointer to row indices
    const int *row_index_ptr() const {
        return row_index_.data();
    }

    /// return the pointer to column indices
    const int *col_index_ptr() const {
        if (col_index_.size()) {
            return col_index_.data();
        } else {
            return 0;
        }
    }

    /// C := A+beta*op(B)
    CSparseMatrix AddTo(char transa, T beta, const CSparseMatrix& B) const {
        switch (transa) {
        case 'N':
            assert(num_rows() == B.num_rows() && num_cols() == B.num_cols());
            break;
        case 'T':
            assert(num_rows() == B.num_cols() && num_cols() == B.num_rows());
            break;
        default:
            ThrowLogicError("undefined label");
        }

        CSparseMatrix C;
        C.ncols_ = ncols_;
        C.row_index_.resize(num_rows() + 1);
        int job = 1;
        int sort = 0;
        int m = num_rows();
        int n = num_cols();
        int ierr;
        if (std::is_same<T, float>::value) {
            mkl_scsradd(&transa, &job, &sort, &m, &n, (float *)element_ptr(), (int *)col_index_ptr(), (int *)row_index_ptr(), (float *)&beta, (float *)B.element_ptr(), (int *)B.col_index_ptr(), (int *)B.row_index_ptr(), 0, 0, C.row_index_.data(), 0, &ierr);
        } else {
            mkl_dcsradd(&transa, &job, &sort, &m, &n, (double *)element_ptr(), (int *)col_index_ptr(), (int *)row_index_ptr(), (double *)&beta, (double *)B.element_ptr(), (int *)B.col_index_ptr(), (int *)B.row_index_ptr(), 0, 0, C.row_index_.data(), 0, &ierr);
        }
        if (ierr) {
            ThrowRuntimeError("mkl_?csradd() returned error = %d", ierr);
        }
        int nzmax = C.row_index_.back() - 1;
        C.element_.resize(nzmax);
        C.col_index_.resize(nzmax);
        nzmax++; // BUG: ???
        job = 0;
        if (std::is_same<T, float>::value) {
            mkl_scsradd(&transa, &job, &sort, &m, &n, (float *)element_ptr(), (int *)col_index_ptr(), (int *)row_index_ptr(), (float *)&beta, (float *)B.element_ptr(), (int *)B.col_index_ptr(), (int *)B.row_index_ptr(), (float *)C.element_ptr(), (int *)C.col_index_ptr(), (int *)C.row_index_ptr(), &nzmax, &ierr);
        } else {
            mkl_dcsradd(&transa, &job, &sort, &m, &n, (double *)element_ptr(), (int *)col_index_ptr(), (int *)row_index_ptr(), (double *)&beta, (double *)B.element_ptr(), (int *)B.col_index_ptr(), (int *)B.row_index_ptr(), (double *)C.element_ptr(), (int *)C.col_index_ptr(), (int *)C.row_index_ptr(), &nzmax, &ierr);
        }
        if (ierr) {
            ThrowRuntimeError("mkl_?csradd() returned error = %d", ierr);
        }
        return C;
    }

    /// scale and assign.
    /// C = op(A)*B
    CSparseMatrix MultiplyTo(char transa, const CSparseMatrix& B) const {
        CSparseMatrix C;
        switch (transa) {
        case 'N':
            assert(num_cols() == B.num_rows());
            C.row_index_.resize(num_rows() + 1);
            break;
        case 'T':
            assert(num_rows() == B.num_rows());
            C.row_index_.resize(num_cols() + 1);
            break;
        default:
            ThrowLogicError("undefined label");
        }

        C.ncols_ = B.ncols_;
        int job = 1;
        int sort = 0;
        int m = num_rows();
        int n = num_cols();
        int ierr;
        if (std::is_same<T, float>::value) {
            mkl_scsrmultcsr(&transa, &job, &sort, &m, &n, (int *)&B.ncols_, (float *)element_ptr(), (int *)col_index_ptr(), (int *)row_index_ptr(), (float *)B.element_ptr(), (int *)B.col_index_ptr(), (int *)B.row_index_ptr(), 0, 0, C.row_index_.data(), 0, &ierr);
        } else {
            mkl_dcsrmultcsr(&transa, &job, &sort, &m, &n, (int *)&B.ncols_, (double *)element_ptr(), (int *)col_index_ptr(), (int *)row_index_ptr(), (double *)B.element_ptr(), (int *)B.col_index_ptr(), (int *)B.row_index_ptr(), 0, 0, C.row_index_.data(), 0, &ierr);
        }
        if (ierr) {
            ThrowRuntimeError("mkl_?csrmultcsr() returned error = %d", ierr);
        }
        int nzmax = C.row_index_.back() - 1;
        C.element_.resize(nzmax);
        C.col_index_.resize(nzmax);
        job = 0;
        if (std::is_same<T, float>::value) {
            mkl_scsrmultcsr(&transa, &job, &sort, &m, &n, (int *)&B.ncols_, (float *)element_ptr(), (int *)col_index_ptr(), (int *)row_index_ptr(), (float *)B.element_ptr(), (int *)B.col_index_ptr(), (int *)B.row_index_ptr(), (float *)C.element_ptr(), (int *)C.col_index_ptr(), (int *)C.row_index_ptr(), &nzmax, &ierr);
        } else {
            mkl_dcsrmultcsr(&transa, &job, &sort, &m, &n, (int *)&B.ncols_, (double *)element_ptr(), (int *)col_index_ptr(), (int *)row_index_ptr(), (double *)B.element_ptr(), (int *)B.col_index_ptr(), (int *)B.row_index_ptr(), (double *)C.element_ptr(), (int *)C.col_index_ptr(), (int *)C.row_index_ptr(), &nzmax, &ierr);
        }
        if (ierr) {
            ThrowRuntimeError("mkl_?csrmultcsr() returned error = %d", ierr);
        }
        return C;
    }

    /// scale
    CSparseMatrix& Scale(T scale) {
        for (int i = 0; i < num_nonzeros(); i++) {
            element_[i] *= scale;
        }
        return *this;
    }

    /// multiply two matrices.
    /// dense <= op(sparse) x dense.
    CMatrix<T> MultiplyTo(char transa, const CMatrix<T>& B) const {
        CMatrix<T> C;
        switch (transa) {
        case 'N':
            assert(num_cols() == B.num_rows());
            C.resize(num_rows(), B.num_cols());
            break;
        case 'T':
            assert(num_rows() == B.num_rows());
            C.resize(num_cols(), B.num_cols());
            break;
        default:
            ThrowLogicError("undefined label");
        }

        int m = num_rows();
        int n = C.num_cols();
        int k = num_cols();
        T alpha = 1;
        char matdescra[] = {
            'G', // general
            0, // lower triangular
            0, // non-unit
            'F', // one-based indexing
            0, // reserved
            0, // reserved
        };
        int ldb = B.num_rows();
        T beta = 0;
        int ldc = C.num_rows();
        if (std::is_same<T, float>::value) {
            mkl_scsrmm(&transa, &m, &n, &k, (float *)&alpha, matdescra, (float *)element_ptr(), (int *)col_index_ptr(), (int *)row_index_ptr(), (int *)row_index_ptr() + 1, (float *)B.data(), &ldb, (float *)&beta, (float *)C.data(), &ldc);
        } else {
            mkl_dcsrmm(&transa, &m, &n, &k, (double *)&alpha, matdescra, (double *)element_ptr(), (int *)col_index_ptr(), (int *)row_index_ptr(), (int *)row_index_ptr() + 1, (double *)B.data(), &ldb, (double *)&beta, (double *)C.data(), &ldc);
        }
        return C;
    }

    /// import from a file
    void Read(const std::string& filename) {
        std::clog << "matrix <= " << filename << std::endl;
        std::ifstream ifs(filename);
        if (!ifs.is_open()) {
            ThrowRuntimeError("failed to open: " + filename);
        }
        int nrows;
        ifs >> nrows >> ncols_;
        row_index_.assign(1, 1);
        col_index_.clear();
        element_.clear();
        int row = 0;
        while (1) {
            int r, c;
            double e;
            ifs >> r >> c >> e;
            if (ifs.fail()) {
                break;
            }
            if (r < row) {
                ThrowRuntimeError("out of order");
            }
            for (; row < r; row++) {
                row_index_.push_back(element_.size() + 1);
            }
            col_index_.push_back(c + 1);
            element_.push_back(e);
        }
        for (; row < nrows + 1; row++) {
            row_index_.push_back(element_.size() + 1);
        }
    }

    /// export to a file
    void Write(const std::string& filename) const {
        std::clog << "matrix => " << filename << std::endl;
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs.is_open()) {
            ThrowRuntimeError("failed to open: " + filename);
        }
        ofs << num_rows() << " " << num_cols() << "\n";
        ofs << std::scientific;
        for (int r = 0, nrows = num_rows(); r < nrows; r++) {
            for (int i = row_index_[r] - 1, end = row_index_[r + 1] - 1; i < end; i++) {
                ofs << r << " " << col_index_[i] - 1 << " " << element_[i] << "\n";
            }
        }
    }

    /// convert the matrix format
    CMatrix<T> ToDense() const {
        CMatrix<T> mat(num_rows(), num_cols());
        mat.fill_with(0);
        for (int r = 0; r < num_rows(); r++) {
            for (int idx = row_index_[r] - 1, end = row_index_[r + 1] - 1; idx < end; idx++) {
                mat(r, col_index_[idx ] - 1) = element_[idx  ];
            }
        }
        return mat;
    }

    // discard below diagonal
    CSparseMatrix ToUpperTriangle() const {
        CSparseMatrix ret;
        int nrows = num_rows();
        int nnz = num_nonzeros();
        ret.row_index_.reserve(num_rows() + 1);
        ret.col_index_.reserve(nnz);
        ret.element_.reserve(nnz);
        ret.ncols_ = num_cols();
        ret.row_index_.push_back(1);
        for (int r = 0; r < nrows; r++) {
            for (int idx = row_index_[r] - 1; idx < row_index_[r + 1] - 1; idx++) {
                int c = col_index_ptr()[idx] - 1;
                if (r <= c) {
                    ret.col_index_.push_back(c + 1);
                    ret.element_.push_back(element_[idx]);
                }
            }
            ret.row_index_.push_back(ret.element_.size() + 1);
        }
        return ret;
    }

    bool IsUpperTriangle() const {
        int nrows = num_rows();
        for (int r = 0; r < nrows; r++) {
            for (int idx = row_index_[r] - 1; idx < row_index_[r + 1] - 1; idx++) {
                int c = col_index_[idx] - 1;
                if (r > c) {
                    return false;
                }
            }
        }
        return true;
    }

    bool IsStructuallySymmetric() const {
        auto has_element = [this](int r, int c) {
            for (int idx = row_index_[r] - 1; idx < row_index_[r + 1] - 1; idx++) {
                int col = col_index_[idx] - 1;
                if (col == c) {
                    return true;
                } else if (col > c) {
                    return false;
                }
            }
            return false;
        };
        int nrows = num_rows();
        for (int r = 0; r < nrows; r++) {
            for (int idx = row_index_[r] - 1; idx < row_index_[r + 1] - 1; idx++) {
                int c = col_index_[idx] - 1;
                if (!has_element(c, r)) {
                    return false;
                }
            }
            if (!has_element(r, r)) {
                return false;
            }
        }
        return true;
    }

    bool IsValid() const {
        int nrows = num_rows();
        int nnz = num_nonzeros();
        if (col_index_.size() != nnz || element_.size() != nnz) {
            return false;
        }
        for (int r = 0; r < nrows; r++) {
            if (row_index_[r] > row_index_[r + 1]) {
                return false;
            }
            for (int idx = row_index_[r]; idx < row_index_[r + 1] - 1; idx++) {
                if (col_index_[idx - 1] >= col_index_[idx]) {
                    return false;
                }
            }
        }
        return true;
    }

    friend class MatrixGenerator<T>;

private:
    std::vector<T> element_; // values of non-zero elements
    std::vector<int> row_index_; // row indices
    std::vector<int> col_index_; // column indices of non-zero elements
    int ncols_; // number of columns
};

/// helper class to construct a CSparseMatrix object
/// @tparam T element
template <typename T>
class MatrixGenerator {
public:
    MatrixGenerator() : sorted_(true) {}

    /// clear old data
    void Clear() {
        data_.clear();
        sorted_ = true;
    }

    /// add an element
    void Add(int r, int c, T v) {
        if (v) {
            // check if the elements are in order
            if (sorted_ && !data_.empty()) {
                int pr = std::get<0>(data_.back());
                int pc = std::get<1>(data_.back());
                if (pr > r) {
                    sorted_ = false;
                } else if (pr == r && pc > c) {
                    sorted_ = false;
                }
            }
            data_.push_back(std::make_tuple(r, c, v));
        }
    }

    /// add a submatrix
    template <typename T2>
    void AddBlock(int row, int col, T scale, const CSparseMatrix<T2>& mat) {
        if (scale) {
            for (int r = 0; r < mat.num_rows(); r++) {
                for (int idx = mat.row_index_ptr()[r] - 1; idx < mat.row_index_ptr()[r + 1] - 1; idx++) {
                    int c = mat.col_index_ptr()[idx] - 1;
                    Add(r + row, c + col, scale * mat.element_ptr()[idx]);
                }
            }
        }
    }

    /// add a submatrix
    template <typename T2, int NRows, int NCols>
    void AddBlock(int row, int col, T scale, const CMatrix<T2, NRows, NCols>& mat) {
        if (scale) {
            for (int c = 0; c < mat.num_cols(); c++) {
                for (int r = 0; r < mat.num_rows(); r++) {
                    if (mat(r, c)) {
                        Add(r + row, c + col, scale * mat(r, c));
                    }
                }
            }
        }
    }

    /// construct a CMatrix object
    CMatrix<T> GenerateDense(int nrows, int ncols) const {
        CMatrix<T> mat(nrows, ncols);
        mat.fill_with(0);
        for (auto& e : data_) {
            mat(std::get<0>(e), std::get<1>(e)) += std::get<2>(e);
        }
        return mat;
    }

    /// construct a CSparseMatrix object
    CSparseMatrix<T> GenerateSparse(int nrows, int ncols) {
        CSparseMatrix<T> mat;
        mat.ncols_ = ncols;
        mat.row_index_.clear();
        mat.col_index_.clear();
        mat.element_.clear();
        if (!sorted_) {
            std::sort(data_.begin(), data_.end());
            sorted_ = true;
        }
        for (auto& e : data_) {
            int r = std::get<0>(e);
            int c = std::get<1>(e);
            T v = std::get<2>(e);
            if (r < 0 || r >= nrows || c < 0 || c >= ncols) {
                ThrowRuntimeError("matrix index out of range");
            }
            if (mat.row_index_.size() == r + 1 && mat.col_index_.back() == c + 1) {
                // add
                mat.element_.back() += v;
            } else {
                // append
                while (mat.row_index_.size() < r + 1) {
                    mat.row_index_.push_back(mat.element_.size() + 1);
                }
                mat.col_index_.push_back(c + 1);
                mat.element_.push_back(v);
            }
        }
        while (mat.row_index_.size() < nrows + 1) {
            mat.row_index_.push_back(mat.element_.size() + 1);
        }
        return mat;
    }

private:
    std::vector<std::tuple<int, int, T>> data_; // indices and values of non-zero elements
    bool sorted_;
};


/// construct a CSparseMatrix object
/// @tparam T element
/// @return sparse matrix
template <typename T>
inline
CSparseMatrix<T> make_sparse_identity_matrix(int nrows, T scale) {
    MatrixGenerator<T> gen;
    for (int r = 0; r < nrows; r++) {
        gen.Add(r, r, scale);
    }
    return gen.Generate(nrows, nrows);
}

/// construct a CSparseMatrix object
/// @tparam T element
/// @return sparse matrix
template <typename T>
inline
CSparseMatrix<T> make_sparse_diagonal_matrix(const CVector<T>& diagonal) {
    MatrixGenerator<T> gen;
    for (int r = 0; r < diagonal.num_rows(); r++) {
        gen.Add(r, r, diagonal[r]);
    }
    return gen.GenerateSparse(diagonal.num_rows(), diagonal.num_rows());
}

/// construct a CSparseMatrix object
/// @tparam T element
/// @return sparse matrix
template <typename T>
inline
CSparseMatrix<T> make_sparse_matrix(const CMatrix<T>& mat) {
    MatrixGenerator<T> gen;
    for (int r = 0; r < mat.num_rows(); r++) {
        for (int c = 0; c < mat.num_cols(); c++) {
            int idx = r + c * mat.num_rows();
            T e = mat[idx];
            if (e) {
                gen.Add(r, c, e);
            }
        }
    }
    return gen.GenerateSparse(mat.num_rows(), mat.num_cols());
}

/// print the debug information of  matrix
template <typename T>
inline void Dump(
    CSparseMatrix<T>& mat,
    bool force = false,
    std::ostream& ost = std::cout  // stream
) {
    float rate = 100.0 * mat.num_nonzeros();
    if (rate) {
        rate /= mat.num_rows() * mat.num_cols();
    }
    ost << "matrix = " << mat.num_rows() << "x" << mat.num_cols() << ", " << mat.num_nonzeros() << " non-zeros, " << rate << "% occupied" << std::endl;
    if (force || (mat.num_rows() < 12 && mat.num_cols() < 12)) {
        for (int r = 0; r < mat.num_rows(); r++) {
            std::cout << r << ": ";
            for (int idx = mat.row_index_ptr()[r] - 1; idx < mat.row_index_ptr()[r + 1] - 1; idx++) {
                ost << "(" << r << "," << (mat.col_index_ptr()[idx  ] - 1) << ")" << mat.element_ptr()[idx  ] << "\t";
            }
            ost << std::endl;
        }
    }
}

} // namespace slib
