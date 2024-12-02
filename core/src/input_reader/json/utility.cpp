#include "input_reader/json/utility.h"

namespace InputReader {
namespace Json {
Object::Object(const std::string &name, const std::string &path)
    : name(name), path(path) {}
Attribute::Attribute(const std::string &name) : name(name) {}
} // namespace Json
} // namespace InputReader
