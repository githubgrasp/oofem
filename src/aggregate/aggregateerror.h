#pragma once
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace aggregate {

[[noreturn]] inline void error(const std::string &msg) {
    std::fputs(msg.c_str(), stderr);
    std::fputc('\n', stderr);
    std::exit(EXIT_FAILURE);
}

[[noreturn]] inline void errorf(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::vfprintf(stderr, fmt, args);
    va_end(args);
    std::fputc('\n', stderr);
    std::exit(EXIT_FAILURE);
}

} // namespace aggregate
