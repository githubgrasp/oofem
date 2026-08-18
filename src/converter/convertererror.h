#pragma once
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace converter {


  // Simple message (no formatting)
[[noreturn]] inline void error(const std::string &msg) {
    std::fputs(msg.c_str(), stderr);
    std::fputc('\n', stderr);
    std::exit(EXIT_FAILURE);
}

// printf-style formatting: %s, %d, %f, etc.
[[noreturn]] inline void errorf(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::vfprintf(stderr, fmt, args);
    va_end(args);
    std::fputc('\n', stderr);
    std::exit(EXIT_FAILURE);
}

// Convenience for the old “error4” pattern
[[noreturn]] inline void error4(const char *fmt, int a, int b, int c) {
    errorf(fmt, a, b, c);
}

inline FILE* fopen_or_die(const std::string& name, const char* mode) {
    if (FILE* f = std::fopen(name.c_str(), mode)) return f;
    errorf("Cannot open output file %s", name.c_str());
    return nullptr; // unreachable
}

  
} // namespace converter
