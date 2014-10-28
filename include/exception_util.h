// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <cstdarg>
#include <string>

/// throw runtime_error() with filename and line number
#define ThrowRuntimeError(...) slib::throw_error<std::runtime_error>(__FILE__, __LINE__, __VA_ARGS__ )

/// throw runtime_error() with filename and line number
#define ThrowLogicError(...) slib::throw_error<std::logic_error>(__FILE__, __LINE__, __VA_ARGS__ )

namespace slib {

/// throw error
template <typename Error>
inline
void throw_error(const char *file, int line, const std::string& message) {
    throw Error(
#ifndef NDEBUG
        std::string(file) + "(" + std::to_string(line) + "): " +
#endif
        message);
}

/// construct runtime error
template <typename Error>
inline
void throw_error(const char *file, int line, const char *fmt, ...) {
    std::string message;
    va_list arguments;
    va_start(arguments, fmt);
    int length = _vscprintf(fmt, arguments);
    message.assign(length, 0);
    vsprintf(&message[0], fmt, arguments);
    va_end(arguments);
    throw_error<Error>(file, line, message);
}

}
