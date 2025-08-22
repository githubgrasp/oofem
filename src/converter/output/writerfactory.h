#pragma once
#include <memory>
#include <string>
struct Writer;
std::unique_ptr<Writer> makeWriter(const std::string& gridTypeLower);
