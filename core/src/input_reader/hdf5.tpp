#include "input_reader/hdf5.h"

namespace InputReader {
template <int nComponents>
bool HDF5::readNestedDoubleVector(
    const char *name, std::vector<std::array<double, nComponents>> &out) const {
  std::vector<double> o;
  if (!this->readDoubleVector(name, o)) {
    return false;
  }
  out.resize(o.size() / nComponents);
  for (size_t i = 0; i < (o.size() / nComponents); i++) {
    for (size_t j = 0; j < nComponents; j++) {
      out[i][j] = o[(i * nComponents) + j];
    }
  }

  return true;
}
} // namespace InputReader
