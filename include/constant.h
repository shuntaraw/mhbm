// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

/// preprocessor definitions for visual c++
#define _CRT_SECURE_NO_WARNINGS
/// preprocessor definitions for visual c++
#define _USE_MATH_DEFINES

/// print debug message to console
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
/// internal macro for DLL
///
/// define MHBMLIB_IMPORTS to import symbols from sdk.dll
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
    size
};

/// tags for registration and local deformation
enum class TRANSFORMATION {
    TRANSLATION,///< translation
    RIGID,///< rigid transformation
    SIMILARITY,///< conformal transformation (similarity transformation)
    size
};

}