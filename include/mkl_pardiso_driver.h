// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <string>
#include <iostream>
#include <boost/filesystem.hpp>
#include <mkl_pardiso.h>
#include "mkl_driver.h"
#include "mkl_csr_driver.h"
#include "exception_util.h"

#ifdef _WIN32
#define VC_EXTRALEAN
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h> // GlobalMemoryStatusEx(),GlobalMemoryStatusEx(),Get/SetEnvironmentVariable()
#else
#include <unistd.h> // sysconf()
#endif

#ifdef _WIN32
inline size_t GetTotalSystemMemory() {
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
}
inline std::string GetEnv(const std::string& name) {
    char var[1024] = {0};
    GetEnvironmentVariable(name.c_str(), var, sizeof(var));
    return var;
}
inline void SetEnv(const std::string& name, const std::string& var) {
    SetEnvironmentVariable(name.c_str(), var.c_str());
    std::clog << name << " = " << var << std::endl;
}
#else
inline size_t GetTotalSystemMemory() {
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}
inline std::string GetEnv(const std::string& name) {
    return getenv(name.c_str());
}
inline void SetEnv(const std::string& name, const std::string& var) {
    std::string arg = name + "=" + var;
    setenv(arg.c_str());
    std::clog << name << " = " << var << std::endl;
}
#endif

namespace slib {

inline void SetPardisoEnv() {
    auto env = GetEnv("MKL_PARDISO_OOC_MAX_CORE_SIZE");
    if (env.empty()) {
        SetEnv("MKL_PARDISO_OOC_MAX_CORE_SIZE", std::to_string(GetTotalSystemMemory() / 1024 / 1024));
        env = GetEnv("TEMP");
        if (!boost::filesystem::exists(env)) {
            env = GetEnv("TMP");
        }
        if (!boost::filesystem::exists(env)) {
            env = ".";
        }
        env += "\\ooc_temp";
        SetEnv("MKL_PARDISO_OOC_PATH", env);
    }
}

inline void ThrowPardisoError(int error) {
    std::string message =  "PARDISO error (" + std::to_string(error) + "): " ;
    switch (error) {
    case -1:
        message += "input inconsistent";
        break;
    case -2:
        message += "not enough memory";
        break;
    case -3:
        message += "reordering problem";
        break;
    case -4:
        message += "zero pivot, numerical factorization or iterative refinement problem";
        break;
    case -5:
        message += "unclassified (internal) error";
        break;
    case -6:
        message += "reordering failed (matrix types 11 and 13 only)";
        break;
    case -7:
        message += "diagonal matrix is singular";
        break;
    case -8:
        message += "32-bit integer overflow problem";
        break;
    case -9:
        message += "not enough memory for OOC";
        break;
    case -10:
        message += "problems with opening OOC temporary files";
        break;
    case -11:
        message += "read/write problems with the OOC data file";
        break;
    default:
        message += "unknown error";
        break;
    }
    ThrowRuntimeError(message);
}

// @param mtype:
//  1 = Real and structurally symmetric
//  2 = Real and symmetric positive definite
// -2 = Real and symmetric indefinite
// 11 = Real and nonsymmetric matrix
template <typename T>
inline void PARDISO_driver(int mtype, char transa, const CSparseMatrix<T>& A, const CMatrix<T>& B, CMatrix<T>& X, bool debug) {
    static_assert(std::is_floating_point<T>::value, "T must be floating point");
    assert(mtype == 1 || mtype == 2 || mtype == -2 || mtype == 11);
    assert(transa == 'N' || transa == 'T');
    assert(A.num_rows() == A.num_cols());
    assert(B.num_rows() == A.num_rows());
    assert(&B != &X);
    X.resize(A.num_rows(), B.num_cols());

    // set environment variables
    SetPardisoEnv();

    // initialize
    void *pt[64]; // parameters for pardiso
    MKL_INT iparm[64]; // parameters for pardiso
    PARDISOINIT(pt, &mtype, iparm);

    // set parameters
    iparm[11] = (transa == 'T' ? 2 : 0);
    iparm[26] = (debug ? 1 : 0); // matrix checker
    iparm[27] = std::is_same<T, float>::value ? 1 : 0; // precision
    iparm[59] = 1; // use an out-of-core algorithm when the problem is too large

    // solve
    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT phase;
    MKL_INT n = A.num_rows();
    void *a = const_cast<T *>(A.element_ptr());
    int *ia = const_cast<int *>(A.row_index_ptr());
    int *ja = const_cast<int *>(A.col_index_ptr());
    MKL_INT nrhs = B.num_cols();
    MKL_INT msglvl = debug ? 1 : 0; // message
    void *b = const_cast<T *>(B.data());
    void *x = X.data();
    MKL_INT error;
    if (debug) {
        // analysis
        phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, 0, &nrhs, iparm, &msglvl, 0, 0, &error);
        float mem_in_gb = std::max(iparm[14], iparm[15] + iparm[16]) / 1024. / 1024.;
        std::clog << "total peak memory = " << std::fixed << mem_in_gb << " GB" << std::endl;
        if (error != 0) {
            ThrowPardisoError(error);
        }
        // factorization
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, 0, &nrhs, iparm, &msglvl, 0, 0, &error);
        if (error != 0) {
            ThrowPardisoError(error);
        }
        // solve
        phase = 33;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, 0, &nrhs, iparm, &msglvl, b, x, &error);
    } else {
        // analysis, factorization, solve
        MKL_INT phase = 13;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, 0, &nrhs, iparm, &msglvl, b, x, &error);
    }
    if (error != 0) {
        ThrowPardisoError(error);
    }
    // clear
    phase = -1;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, 0, &nrhs, iparm, &msglvl, b, x, &error);
}

///  driver for PARDISO
/// @tparam T  element
template <typename T>
class PardisoDriver {
    static_assert(std::is_floating_point<T>::value, "T must be floating point");
public:
    PardisoDriver() {
        SetPardisoEnv();
    }
    PardisoDriver(const PardisoDriver&) = delete;
    PardisoDriver& operator=(const PardisoDriver&) = delete;
    PardisoDriver(PardisoDriver&&) = delete;
    PardisoDriver& operator=(PardisoDriver &&) = delete;

    /// destructor
    ~PardisoDriver() {
        Clear();
    }

    // factorize
    // mtype =
    //  1 : Real and structurally symmetric
    //  2 : Real and symmetric positive definite
    // -2 : Real and symmetric indefinite
    // 11 : Real and nonsymmetric matrix
    void Factorize(int mtype, char transa, const CSparseMatrix<T>& A, bool debug) {
        Clear();
        A_ = A;
#ifdef _DEBUG
        if (mtype == 1) {
            assert(A_.IsStructuallySymmetric());
        } else if (mtype == 2 || mtype == -2) {
            assert(A_.IsUpperTriangle());
        }
#endif
        // initialize
        mtype_ = mtype; // real and asymmetric
        PARDISOINIT(pt_, &mtype_, iparm_);

        iparm_[11] = (transa == 'T' ? 2 : 0);
        iparm_[26] = (debug ? 1 : 0); // matrix checker
        iparm_[27] = (std::is_same<T, float>::value ? 1 : 0);		// precision
        iparm_[59] = 1; // use an out-of-core algorithm when the problem is too large

        Phase12(debug);
    }

    /// solve a linear system
    void Solve(const CMatrix<T>& B, CMatrix<T>& X, bool debug) {
        assert(A_.num_cols() == B.num_rows());
        assert(&B != &X);
        X.resize(A_.num_cols(), B.num_cols());

        // solve
        MKL_INT maxfct = 1;
        MKL_INT mnum = 1;
        MKL_INT phase = 33; // solve
        MKL_INT n = A_.num_rows();
        void *a = (void *)A_.element_ptr();
        int *ia = (int *)A_.row_index_ptr();
        int *ja = (int *)A_.col_index_ptr();
        MKL_INT nrhs = B.num_cols();
        MKL_INT msglvl = (debug ? 1 : 0); // message
        void *b = (void *)B.data();
        void *x = X.data();
        MKL_INT error;
        PARDISO(pt_, &maxfct, &mnum, &mtype_, &phase, &n, a, ia, ja, 0, &nrhs, iparm_, &msglvl, b, x, &error);

        if (error != 0) {
            ThrowPardisoError(error);
        }
    }

private:
    /// clear
    void Clear() {
        if (A_.num_rows()) {
            MKL_INT phase = -1;
            MKL_INT n = 0;
            MKL_INT msglvl = 0; // message
            MKL_INT error;
            PARDISO(pt_, 0, 0, &mtype_, &phase, &n, 0, 0, 0, 0, 0, iparm_, &msglvl, 0, 0, &error);
            A_.Clear();
        }
    }

    void Phase12(bool debug) {
        MKL_INT maxfct = 1;
        MKL_INT mnum = 1;
        MKL_INT phase;
        MKL_INT n = A_.num_rows();
        T *a = (T *)A_.element_ptr();
        int *ia = (int *)A_.row_index_ptr();
        int *ja = (int *)A_.col_index_ptr();
        MKL_INT nrhs = 0;
        MKL_INT msglvl = (debug ? 1 : 0); // message
        MKL_INT error;
        if (debug) {
            // analysis
            phase = 11;
            PARDISO(pt_, &maxfct, &mnum, &mtype_, &phase, &n, a, ia, ja, 0, &nrhs, iparm_, &msglvl, 0, 0, &error);
            if (error != 0) {
                ThrowPardisoError(error);
            }
            std::clog << "total peak memory = " << std::fixed << (std::max(iparm_[14], iparm_[15] + iparm_[16]) / 1024.0 / 1024.0) << " GB" << std::endl;
            //  factorization
            phase = 22;
            PARDISO(pt_, &maxfct, &mnum, &mtype_, &phase, &n, a, ia, ja, 0, &nrhs, iparm_, &msglvl, 0, 0, &error);
        } else {
            // analysis and factorization
            phase = 12;
            PARDISO(pt_, &maxfct, &mnum, &mtype_, &phase, &n, a, ia, ja, 0, &nrhs, iparm_, &msglvl, 0, 0, &error);
        }
        if (error != 0) {
            ThrowPardisoError(error);
        }
    }

private:
    CSparseMatrix<T> A_; // system matrix in CSR format
    void *pt_[64]; // parameters for pardiso
    MKL_INT iparm_[64]; // parameters for pardiso
    MKL_INT mtype_; // parameters for pardiso
};

} // namespace slib
