// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

// utility
#define TRACE printf
#include "MatrixBase.h"
#include "MatrixUtil.h"
#include "mkl_blas_driver.h"
#include "mkl_lapack_driver.h"
#include "mkl_csr_driver.h"
#include "mkl_pardiso_driver.h"

#ifdef MHBMLIB_EXPORTS
#define MHBMLIB_API __declspec(dllexport)
#elif defined(MHBMLIB_IMPORTS)
#define MHBMLIB_API __declspec(dllimport)
#else
#define MHBMLIB_API
#endif

namespace hbm {

/// version string
    extern MHBMLIB_API const std::string VERSION;

/// threshold of iterations for ICP
    extern MHBMLIB_API int MIN_ICP_ITERATIONS;

/// threshold of improvement for ICP
    extern MHBMLIB_API float MIN_ICP_IMPROVEMENT;

/// lower bound for scale parameters
    extern MHBMLIB_API float MIN_SIMILARITY_SCALE;

/// upper bound for scale parameters
    extern MHBMLIB_API float MAX_SIMILARITY_SCALE;


    /// tags for search direction
    enum class SEARCH_DIRECTION {
        FORWARD,///< forward
        BACKWARD,///< backward
        BIDIRECTIONAL,///< bidirectional
    };

    enum class TRANSFORMATION {
        TRANSLATION,
        RIGID,
        SIMILARITY,
        size,
    };

}