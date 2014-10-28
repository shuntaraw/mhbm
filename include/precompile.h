// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

// C
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

// C++

// extra C
#include <omp.h> // omp_set_num_threads()

// extra C++
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string_regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

// utility
#define TRACE printf
#include "exception_util.h"
#include "Halfedge.h"
#include "KdTree.h"
#include "MatrixUtil.h"
#include "MeshFile.h"
#include "MeshUtil.h"
#include "MeshSubdivision.h"
#include "mkl_blas_driver.h"
#include "mkl_lapack_driver.h"
#include "mkl_csr_driver.h"
#include "mkl_pardiso_driver.h"
#include "numeric_util.h"