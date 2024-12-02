#include "input_reader/hdf5/partial.h"

namespace InputReader {
namespace HDF5 {
Partial::Partial(const char *file) : File(file) {
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize_);
  MPI_Comm_rank(MPI_COMM_WORLD, &ownRank_);
}

bool Partial::readIntVector(const char *name, std::vector<int32_t> &out,
                            const std::string &groupName) const {
  return this->readDataset(name, groupName, H5T_NATIVE_INT, out);
}

bool Partial::readDoubleVector(const char *name, std::vector<double> &out,
                               const std::string &groupName) const {
  return this->readDataset(name, groupName, H5T_NATIVE_DOUBLE, out);
}
} // namespace HDF5
} // namespace InputReader
