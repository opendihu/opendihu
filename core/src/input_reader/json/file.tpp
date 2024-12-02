#include "input_reader/json/file.h"

namespace InputReader {
namespace Json {
template <typename T> bool File::readAttr(const char *name, T &out) const {
  out = content_["__attributes"][name].template get<T>();
  return true;
}
} // namespace Json
} // namespace InputReader
