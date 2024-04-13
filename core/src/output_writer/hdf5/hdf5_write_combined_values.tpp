#include "output_writer/hdf5/hdf5.h"

#include "easylogging++.h"

#include "control/diagnostic_tool/performance_measurement.h"

namespace OutputWriter {

template <typename T>
herr_t
HDF5::writeCombinedValuesVector(hid_t fileID, const std::vector<T> &values,
                                const char *dsname, bool writeFloatsAsInt) {
  VLOG(1) << "HDF5::writeCombinedValuesVector, " << values.size()
          << " values: " << values;
  VLOG(1) << "rankSubset: " << *this->rankSubset_;
  LOG(DEBUG) << "HDF5::writeCombinedValuesVector: fileID: " << fileID
             << " dsname: " << dsname;

  if (std::is_same<T, double>::value || std::is_same<T, float>::value) {
    if (writeFloatsAsInt) {
      std::vector<int32_t> valuesVector(values.size());

      for (const T &v : values) {
        valuesVector.push_back((int32_t)(round(v)));
      }
      return HDF5Utils::writeSimpleVec<int32_t>(fileID, valuesVector, dsname);
    } else if (std::is_same<T, double>::value) {
      return HDF5Utils::writeSimpleVec<T>(fileID, values, dsname);
    } else if (std::is_same<T, float>::value) {
      return HDF5Utils::writeSimpleVec<T>(fileID, values, dsname);
    }
  } else {
    return HDF5Utils::writeSimpleVec<T>(fileID, values, dsname);
  }
}
} // namespace OutputWriter
