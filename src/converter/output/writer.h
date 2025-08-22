#pragma once
#include <string>
class Grid;
struct Writer {
  virtual ~Writer() = default;
  virtual void write(const Grid& grid, const std::string& baseName) = 0;
};
