#include "input_reader/hdf5/file.h"

namespace InputReader {
namespace HDF5 {
template <typename T, std::enable_if_t<std::is_same<T, int32_t>::value, bool>>
bool File::readAttr(const char *name, T &out) const {
  return this->readAttribute(name, (void *)&out) >= 0;
}

template <typename T, std::enable_if_t<std::is_same<T, double>::value, bool>>
bool File::readAttr(const char *name, T &out) const {
  return this->readAttribute(name, (void *)&out) >= 0;
}

template <typename T,
          std::enable_if_t<std::is_same<T, std::string>::value, bool>>
bool File::readAttr(const char *name, T &out) const {
  bool found = false;
  for (const auto &e : attributes_) {
    if (e.name == name) {
      out.resize(e.size + 1, '\0');
      found = true;
      break;
    }
  }
  assert(found);
  return this->readAttribute(name, (void *)out.data()) >= 0;
}
} // namespace HDF5
} // namespace InputReader
