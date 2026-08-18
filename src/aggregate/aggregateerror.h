#pragma once
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace aggregate {

/// Print `msg` to stderr and exit with `EXIT_FAILURE`. Use for fatal,
/// non-recoverable problems detected anywhere in the aggregate module.
[[noreturn]] inline void error(const std::string &msg) {
    std::fputs(msg.c_str(), stderr);
    std::fputc('\n', stderr);
    std::exit(EXIT_FAILURE);
}

/// printf-style fatal: format with `fmt` and varargs, write to stderr,
/// then exit with `EXIT_FAILURE`. Avoids string concatenation at call
/// sites that report numeric context (file paths, indices, sizes).
[[noreturn]] inline void errorf(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::vfprintf(stderr, fmt, args);
    va_end(args);
    std::fputc('\n', stderr);
    std::exit(EXIT_FAILURE);
}

} // namespace aggregate
