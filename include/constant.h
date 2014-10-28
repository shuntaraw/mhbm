// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#ifdef MHBMLIB_EXPORTS
#define MHBMLIB_API __declspec(dllexport)
#elif defined(MHBMLIB_IMPORTS)
#define MHBMLIB_API __declspec(dllimport)
#else
#define MHBMLIB_API 
#endif

namespace hbm {

    /// version string
    static const std::string VERSION = boost::replace_all_copy<std::string>("1.0.0 ($Revision: 8625 $)", "$", "");

    /// threshold of iterations for ICP 
    static int MIN_ICP_ITERATIONS = 5;

    /// threshold of improvement for ICP 
    static float MIN_ICP_IMPROVEMENT = 1e-6;

    /// lower bound for scale parameters
    static float MIN_SIMILARITY_SCALE = 0.1;

    /// upper bound for scale parameters
    static float MAX_SIMILARITY_SCALE = 1.5;

}