#include "input_reader/hdf5/full_dataset.h"

namespace InputReader {
namespace HDF5 {
template <typename T>
bool FullDataset::readDataset(const char *name, const std::string &groupName,
                              hid_t memTypeId, std::vector<T> &out) const {
  const std::string *fullName = getFullDatasetName(name, groupName);
  if (!fullName) {
    return false;
  }
  herr_t err;
  hid_t dset = H5Dopen(fileID_, fullName->c_str(), H5P_DEFAULT);
  assert(dset >= 0);

  hid_t dspace = H5Dget_space(dset);
  assert(dspace >= 0);
  const int ndims = H5Sget_simple_extent_ndims(dspace);
  if (ndims != 1) {
    err = H5Sclose(dspace);
    assert(err >= 0);
    err = H5Dclose(dset);
    assert(err >= 0);
    return false;
  }
  std::array<hsize_t, 1> dims = {0};
  err = H5Sget_simple_extent_dims(dspace, dims.data(), nullptr);
  assert(err >= 0);
  err = H5Sclose(dspace);
  assert(err >= 0);

  out.resize(dims[0]);
  err = H5Dread(dset, memTypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
  assert(err >= 0);
  err = H5Dclose(dset);
  assert(err >= 0);

  return true;
}
} // namespace HDF5
} // namespace InputReader
