#include "input_reader/generic.h"

namespace InputReader {
template <int D>
bool Generic::readDoubleVecD(const char *name,
                             std::array<std::vector<double>, D> &out,
                             const std::string &groupName) const {
  std::vector<double> values;
  bool ret = this->readDoubleVector(name, values, groupName);
  if (!ret) {
    return ret;
  }
  for (size_t i = 0; i < D; i++) {
    out[i] = std::vector<double>();
    out[i].reserve(values.size() / D);
  }
  for (size_t i = 0; i < values.size(); i++) {
    out[i % D].push_back(values[i]);
  }

  return true;
}
} // namespace InputReader
