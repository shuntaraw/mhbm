// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#include "constant.h"

#include <boost/algorithm/string.hpp>

namespace hbm {

/// version string
MHBMLIB_API const std::string VERSION = boost::replace_all_copy<std::string>("1.0.0 ($Revision: 8625 $)", "$", "");

/// threshold of iterations for ICP
MHBMLIB_API int MIN_ICP_ITERATIONS = 5;

/// threshold of improvement for ICP
MHBMLIB_API float MIN_ICP_IMPROVEMENT = 1e-6;

/// lower bound for scale parameters
MHBMLIB_API float MIN_SIMILARITY_SCALE = 0.1;

/// upper bound for scale parameters
MHBMLIB_API float MAX_SIMILARITY_SCALE = 1.5;

}