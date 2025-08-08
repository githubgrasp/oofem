#pragma once
#include <vector>
#include <iostream>
#include <cstdlib>

namespace generator {

template <class T>
inline bool includes1(const std::vector<T*> &v, int i) {
    return i >= 1 && static_cast<size_t>(i) <= v.size() && v[i - 1] != nullptr;
}

template <class T>
inline void ensure_size1(std::vector<T*> &v, int i) {
    if (i < 1) return;
    if (static_cast<size_t>(i) > v.size()) v.resize(i, nullptr);
}

template <class T>
inline void put1(std::vector<T*> &v, int i, T* ptr) {
    ensure_size1(v, i);
    v[i - 1] = ptr;
}

template <class T>
inline void put1_replace(std::vector<T*> &v, int i, T* ptr) {
    if (i < 1) {
        std::cerr << "put1_replace: index must be >= 1 (got " << i << ")\n";
        std::exit(EXIT_FAILURE);
    }
    ensure_size1(v, i);
    if (v[i - 1] && v[i - 1] != ptr) delete v[i - 1];
    v[i - 1] = ptr;
}

template <class T>
inline T* at1(const std::vector<T*> &v, int i) {
    if (!includes1(v, i)) return nullptr;
    return v[i - 1];
}

template <class T>
inline int size1(const std::vector<T*> &v) {
    return static_cast<int>(v.size());
}
 
} // namespace generator
