#include "output_writer/hdf5/hdf5.h"

namespace OutputWriter {
HDF5::HDF5(DihuContext context, PythonConfig settings,
           std::shared_ptr<Partition::RankSubset> rankSubset)
    : Generic(context, settings, rankSubset) {
  combineFiles_ = settings.getOptionBool("combineFiles", false);
  writeMeta_ = settings.getOptionBool("writeMeta", true);
}

//! constructor, initialize nPoints and nCells to 0
HDF5::Piece::Piece() {
  properties.nPointsLocal = 0;
  properties.nCellsLocal = 0;
  properties.nPointsGlobal = 0;
  properties.nCellsGlobal = 0;
  properties.dimensionality = 0;
}

void HDF5::setCombineFiles(bool v) { combineFiles_ = v; }
void HDF5::setWriteMeta(bool v) { writeMeta_ = v; }

//! assign the correct values to firstScalarName and firstVectorName, only if
//! properties has been set
void HDF5::Piece::setVTKValues() {
  // set values for firstScalarName and firstVectorName from the values in
  // pointDataArrays
  for (auto pointDataArray : properties.pointDataArrays) {
    if (firstScalarName == "" && pointDataArray.nComponents != 3)
      firstScalarName = pointDataArray.name;
    if (firstVectorName == "" && pointDataArray.nComponents == 3)
      firstVectorName = pointDataArray.name;
  }
}

herr_t HDF5::writeCombinedTypesVector(hid_t fileID, uint64_t nValues,
                                      bool output3DMeshes, const char *dsname) {
  if (output3DMeshes) {
    std::vector<int32_t> values(nValues, 12);
    return HDF5Utils::writeSimpleVec<int32_t>(fileID, values, dsname);
  } else {
    std::vector<int32_t> values(nValues, 9);
    return HDF5Utils::writeSimpleVec<int32_t>(fileID, values, dsname);
  }
}

hid_t HDF5::openHDF5File(const char *filename, bool mpiio) {
  if (mpiio) {
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    herr_t err;
    assert(plist >= 0);
    err = H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(err >= 0);
    err = H5Pset_all_coll_metadata_ops(plist, true);
    assert(err >= 0);
    err = H5Pset_coll_metadata_write(plist, true);
    assert(err >= 0);

    hid_t fileID = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    assert(fileID >= 0);
    err = H5Pclose(plist);
    assert(err >= 0);

    return fileID;
  } else {
    hid_t fileID = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(fileID >= 0);
    return fileID;
  }
}

namespace HDF5Utils {
std::pair<int, int> __getMemSpace(hsize_t len, int32_t ownRank,
                                  int32_t worldSize) {
  if (worldSize > len) {
    return std::make_pair(len, len);
  }

  if (len % worldSize == 0) {
    int d = len / worldSize;
    return std::make_pair(d, d);
  } else {
    int32_t split = (len / worldSize);
    if (ownRank == (worldSize - 1)) {
      return std::make_pair((len - (split * worldSize)) + split, split);
    } else {
      return std::make_pair(split, split);
    }
  }
}
} // namespace HDF5Utils
} // namespace OutputWriter
