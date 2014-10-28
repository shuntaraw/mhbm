// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#include <cstdarg>
#include <string>

namespace slib {

/// construct a string using the printf format
inline
std::string string_format(const char *fmt, ...) {
    std::string buffer;
    va_list arguments;
    va_start(arguments, fmt);
    int length = _vscprintf(fmt, arguments);
    buffer.assign(length, 0);
    vsprintf(&buffer[0], fmt, arguments);
    va_end(arguments);
    return buffer;
}

/// @see http://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
/// @deprecated use <boost/algorithm/string/replace.hpp> instead.
inline
std::string string_replace(const std::string& s, const std::string& from, const std::string& to) {
    if (from.empty()) {
        return s;
    }
    std::string str = s;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
    return str;
}

}
