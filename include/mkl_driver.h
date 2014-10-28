// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

/// @see http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
#ifdef _DLL // /MD
#pragma comment(lib, "mkl_rt.lib")
#else // /MT
#if defined(_M_IX86)
#pragma comment(lib, "mkl_intel_c.lib")
#elif defined(_M_X64)
#pragma comment(lib, "mkl_intel_lp64.lib")
#endif
#pragma comment(lib, "mkl_intel_thread.lib")
#pragma comment(lib, "mkl_core.lib")
#pragma comment(lib, "libiomp5md.lib")
#endif

#include <mkl.h>
