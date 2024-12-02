#include "input_reader/hdf5/full_dataset.h"

namespace InputReader {
namespace HDF5 {
FullDataset::FullDataset(const char *file) : File(file) {}

bool FullDataset::readIntVector(const char *name, std::vector<int32_t> &out,
                                const std::string &groupName) const {
  return this->readDataset(name, groupName, H5T_NATIVE_INT, out);
}

bool FullDataset::readDoubleVector(const char *name, std::vector<double> &out,
                                   const std::string &groupName) const {
  return this->readDataset(name, groupName, H5T_NATIVE_DOUBLE, out);
}
} // namespace HDF5
} // namespace InputReader
