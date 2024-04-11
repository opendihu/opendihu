#include "output_writer/hdf5/hdf5.h"

#include "easylogging++.h"

#include "control/diagnostic_tool/performance_measurement.h"

namespace OutputWriter {

template <typename T>
void HDF5::writeCombinedValuesVector(hid_t fileID, const std::vector<T> &values,
                                     const char *dsname,
                                     bool writeFloatsAsInt) {
  VLOG(1) << "HDF5::writeCombinedValuesVector, " << values.size()
          << " values: " << values;
  VLOG(1) << "rankSubset: " << *this->rankSubset_;

  herr_t err;
  std::array<hsize_t, 1> dims = {values.size()};
  hid_t dspace = H5Screate_simple(1, dims.data(), nullptr);
  assert(dspace >= 0);
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  assert(plist >= 0);
  err = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  assert(err >= 0);

  // convert float values to int32 and insert into valuesVector
  if (std::is_same<T, double>::value || std::is_same<T, float>::value) {
    if (writeFloatsAsInt) {
      std::vector<int32_t> valuesVector(values.size());

      for (const T &v : values) {
        valuesVector.push_back((int32_t)(round(v)));
      }
    } else if (std::is_same<T, double>::value) {
      hid_t dset = H5Dcreate(fileID, dsname, H5T_IEEE_F64BE, dspace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      assert(dset >= 0);
      err = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist,
                     values.data());
      assert(err >= 0);
      err = H5Dclose(dset);
      assert(err >= 0);
    } else if (std::is_same<T, float>::value) {
      hid_t dset = H5Dcreate(fileID, dsname, H5T_IEEE_F32BE, dspace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      assert(dset >= 0);
      err = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, plist,
                     values.data());
      assert(err >= 0);
      err = H5Dclose(dset);
      assert(err >= 0);
    }
  } else {
    hid_t dset = H5Dcreate(fileID, dsname, H5T_STD_I32BE, dspace, H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
    assert(dset >= 0);
    err =
        H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist, values.data());
    assert(err >= 0);
    err = H5Dclose(dset);
    assert(err >= 0);
  }
  err = H5Pclose(plist);
  assert(err >= 0);
  err = H5Sclose(dspace);
  assert(err >= 0);
}
} // namespace OutputWriter
