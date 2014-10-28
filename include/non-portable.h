// Copyright (c) 2012-2014 Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
// and the National Institute of Advanced Industrial Science and Technology

#pragma once

#ifdef _WIN32
#define VC_EXTRALEAN
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h> // GlobalMemoryStatusEx(),GlobalMemoryStatusEx(),Get/SetEnvironmentVariable()
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
