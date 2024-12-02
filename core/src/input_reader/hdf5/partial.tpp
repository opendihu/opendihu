#include "input_reader/hdf5/partial.h"

namespace InputReader {
namespace HDF5 {
template <typename T>
bool Partial::readDataset(const char *name, const std::string &groupName,
                          hid_t memTypeId, std::vector<T> &out) const {
  const herr_t RANK = 1;

  const std::string *fullName = getFullDatasetName(name, groupName);
  if (!fullName) {
    LOG(ERROR) << "Name: " << name << " not found";
    return false;
  }
  herr_t err;
  hid_t dset = H5Dopen(fileID_, fullName->c_str(), H5P_DEFAULT);
  assert(dset >= 0);

  std::vector<hsize_t> partition;
  {
    // chunkDims attribute has always 1 dim
    std::array<hsize_t, 1> part_dims = {0};
    hid_t attr = H5Aopen_name(dset, "chunkDims");
    assert(attr >= 0);
    hid_t dspace = H5Aget_space(attr);
    assert(dspace >= 0);

    H5Sget_simple_extent_dims(dspace, part_dims.data(), nullptr);
    partition.resize(part_dims[0]);
    H5Aread(attr, H5T_NATIVE_LONG, partition.data());

    err = H5Sclose(dspace);
    assert(err >= 0);
    err = H5Aclose(attr);
    assert(err >= 0);
  }
  if (partition.size() != worldSize_) {
    LOG(ERROR) << "Wrong partition size";
    return false;
  }

  std::array<hsize_t, RANK> count = {1};
  std::array<hsize_t, RANK> stride = {1};
  std::array<hsize_t, RANK> block = {partition[ownRank_]};
  std::array<hsize_t, RANK> offset = {
      std::accumulate(partition.begin(), partition.begin() + ownRank_, 0ul)};

  hid_t filespace = H5Dget_space(dset);
  assert(filespace >= 0);
  const int ndims = H5Sget_simple_extent_ndims(filespace);
  if (ndims != 1) {
    LOG(ERROR) << "Dims not 1";
    err = H5Sclose(filespace);
    assert(err >= 0);
    err = H5Dclose(dset);
    assert(err >= 0);
    return false;
  }

  std::array<hsize_t, 1> dims = {0};
  err = H5Sget_simple_extent_dims(filespace, dims.data(), nullptr);
  assert(err >= 0);
  err = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(),
                            stride.data(), count.data(), block.data());
  assert(err >= 0);
  hid_t memspace = H5Screate_simple(RANK, block.data(), nullptr);
  assert(memspace >= 0);

  out.resize(block[0]);
  err = H5Dread(dset, memTypeId, memspace, filespace, H5P_DEFAULT, out.data());
  assert(err >= 0);

  err = H5Sclose(memspace);
  assert(err >= 0);
  err = H5Sclose(filespace);
  assert(err >= 0);
  err = H5Dclose(dset);
  assert(err >= 0);

  return true;
}
} // namespace HDF5
} // namespace InputReader
