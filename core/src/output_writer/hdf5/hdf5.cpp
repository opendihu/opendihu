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

herr_t HDF5::writeCombinedTypesVector(HDF5Utils::Group &group, uint64_t nValues,
                                      bool output3DMeshes, const char *dsname) {
  if (output3DMeshes) {
    std::vector<int32_t> values(nValues, 12);
    return group.writeSimpleVec<int32_t>(values, dsname);
  } else {
    std::vector<int32_t> values(nValues, 9);
    return group.writeSimpleVec<int32_t>(values, dsname);
  }
}

namespace HDF5Utils {
File::File(const char *filename, bool mpiio)
    : filename_(filename), mpiio_(mpiio), fileID_(-1) {
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize_);
  MPI_Comm_rank(MPI_COMM_WORLD, &ownRank_);

  if (mpiio_) {
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    herr_t err;
    assert(plist >= 0);
    err = H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(err >= 0);
    err = H5Pset_all_coll_metadata_ops(plist, true);
    assert(err >= 0);
    err = H5Pset_coll_metadata_write(plist, true);
    assert(err >= 0);

    fileID_ = H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    assert(fileID_ >= 0);
    err = H5Pclose(plist);
    assert(err >= 0);
  } else {
    fileID_ =
        H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(fileID_ >= 0);
  }
}

File::~File() {
  if (fileID_ >= 0) {
    H5Fclose(fileID_);
  }
}

hid_t File::getFileID() const { return fileID_; }
int32_t File::getOwnRank() const { return ownRank_; }
int32_t File::getWorldSize() const { return worldSize_; }
bool File::isMPIIO() const { return mpiio_; }

Group File::newGroup(const char *name) const { return Group(this, name); }

herr_t File::writeAttrInt(const char *key, int32_t value) const {
  hsize_t dims = 1;
  std::array<int32_t, 1> data = {value};
  hid_t dspace = H5Screate_simple(1, &dims, nullptr);
  if (dspace < 0) {
    return dspace;
  }
  hid_t attr =
      H5Acreate(fileID_, key, H5T_STD_I32LE, dspace, H5P_DEFAULT, H5P_DEFAULT);
  if (attr < 0) {
    // ignore error here, everything is lost anyway at this point in time
    H5Sclose(dspace);
    return attr;
  }
  herr_t err = H5Awrite(attr, H5T_NATIVE_INT, data.data());
  if (err < 0) {
    return err;
  }
  err = H5Aclose(attr);
  if (err < 0) {
    return err;
  }
  err = H5Sclose(dspace);
  if (err < 0) {
    return err;
  }

  return 0;
}

herr_t File::writeAttrDouble(const char *key, double value) const {
  hsize_t dims = 1;
  std::array<double, 1> data = {value};
  hid_t dspace = H5Screate_simple(1, &dims, nullptr);
  if (dspace < 0) {
    return dspace;
  }
  hid_t attr =
      H5Acreate(fileID_, key, H5T_IEEE_F64LE, dspace, H5P_DEFAULT, H5P_DEFAULT);
  if (attr < 0) {
    // ignore error here, everything is lost anyway at this point in time
    H5Sclose(dspace);
    return attr;
  }
  herr_t err = H5Awrite(attr, H5T_NATIVE_DOUBLE, data.data());
  if (err < 0) {
    return err;
  }
  err = H5Aclose(attr);
  if (err < 0) {
    return err;
  }
  err = H5Sclose(dspace);
  if (err < 0) {
    return err;
  }

  return 0;
}

herr_t File::writeAttrStr(const char *key, const std::string &value) const {
  hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
  herr_t err = H5Tset_size(filetype, value.length());
  if (err < 0) {
    return err;
  }

  hid_t memtype = H5Tcopy(H5T_C_S1); // Datatype ID
  err = H5Tset_size(memtype, value.length() + 1);
  if (err < 0) {
    return err;
  }

  hsize_t dims = 1;
  hid_t dspace = H5Screate_simple(1, &dims, nullptr);
  if (dspace < 0) {
    return dspace;
  }
  hid_t attr =
      H5Acreate(fileID_, key, filetype, dspace, H5P_DEFAULT, H5P_DEFAULT);
  if (attr < 0) {
    // ignore error here, everything is lost anyway at this point in time
    H5Sclose(dspace);
    return attr;
  }
  err = H5Awrite(attr, memtype, value.c_str());
  if (err < 0) {
    H5Sclose(dspace);
    H5Aclose(attr);
    return err;
  }
  err = H5Aclose(attr);
  if (err < 0) {
    return err;
  }
  err = H5Sclose(dspace);
  if (err < 0) {
    return err;
  }

  return 0;
}

Group::Group(const File *f, const char *name) : file_(f), groupID_(-1) {
  groupID_ = H5Gcreate(file_->getFileID(), name, H5P_DEFAULT, H5P_DEFAULT,
                       H5P_DEFAULT);
  assert(groupID_ >= 0);
}

Group::Group(const File *f, hid_t id, const char *name)
    : file_(f), groupID_(-1) {
  groupID_ = H5Gcreate(id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(groupID_ >= 0);
}

Group::~Group() {
  if (groupID_ >= 0) {
    H5Gclose(groupID_);
  }
}

Group Group::newGroup(const char *name) const {
  return Group(file_, groupID_, name);
}

herr_t Group::writeVectorMPIIO(const void *data, const std::string &dsname,
                               const size_t dataSize, hid_t typeId,
                               hid_t memTypeId, size_t dsize) {
  const herr_t RANK = 1;

  std::vector<int64_t> sizes;
  sizes.resize(file_->getWorldSize());
  MPI_Allgather(&dataSize, 1, MPI_LONG, sizes.data(), 1, MPI_LONG,
                MPI_COMM_WORLD);

  herr_t err;
  std::array<hsize_t, RANK> dims = {
      std::accumulate(sizes.begin(), sizes.end(), 0ul)};
  std::array<hsize_t, RANK> chunkDims = {dataSize};

  hid_t filespace = H5Screate_simple(RANK, dims.data(), nullptr);
  if (filespace < 0) {
    return filespace;
  }
  std::string name = dsname;
  std::replace(name.begin(), name.end(), '/', '|');
  hid_t dset = H5Dcreate(groupID_, name.c_str(), typeId, filespace, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) {
    H5Sclose(filespace);
    return dset;
  }
  err = H5Sclose(filespace);
  if (err < 0) {
    H5Dclose(dset);
    return err;
  }

  filespace = H5Dget_space(dset);
  if (filespace < 0) {
    H5Dclose(dset);
    return filespace;
  }
  std::array<hsize_t, RANK> count = {1};
  std::array<hsize_t, RANK> stride = {1};
  std::array<hsize_t, RANK> block = chunkDims;
  std::array<hsize_t, RANK> offset = {
      std::accumulate(sizes.begin(), sizes.begin() + file_->getOwnRank(), 0ul)};

  err = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(),
                            stride.data(), count.data(), block.data());
  if (err < 0) {
    H5Dclose(dset);
    H5Sclose(filespace);
    return err;
  }

  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  if (plist < 0) {
    H5Dclose(dset);
    H5Sclose(filespace);
    return plist;
  }
  err = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  if (err < 0) {
    H5Dclose(dset);
    H5Sclose(filespace);
    H5Pclose(plist);
    return err;
  }

  hid_t memspace = H5Screate_simple(RANK, chunkDims.data(), nullptr);
  if (memspace < 0) {
    H5Dclose(dset);
    H5Sclose(filespace);
    H5Pclose(plist);
    return memspace;
  }
  err = H5Dwrite(dset, memTypeId, memspace, filespace, plist, data);
  if (err < 0) {
    H5Dclose(dset);
    H5Sclose(filespace);
    H5Pclose(plist);
    H5Sclose(memspace);
    return err;
  }

  err = writeAttrArray(dset, "chunkDims", sizes);
  if (err < 0) {
    H5Dclose(dset);
    H5Sclose(filespace);
    H5Pclose(plist);
    H5Sclose(memspace);
    return err;
  }

  err = H5Dclose(dset);
  if (err < 0) {
    return err;
  }
  err = H5Sclose(filespace);
  if (err < 0) {
    return err;
  }
  err = H5Pclose(plist);
  if (err < 0) {
    return err;
  }
  err = H5Sclose(memspace);
  if (err < 0) {
    return err;
  }

  return 0;
}

herr_t Group::writeVector(const void *data, const std::string &dsname,
                          const size_t dataSize, hid_t typeId, hid_t memTypeId,
                          size_t dsize) {
  const herr_t RANK = 1;
  std::array<hsize_t, 1> dims = {dataSize};

  herr_t err;
  hid_t filespace = H5Screate_simple(RANK, dims.data(), nullptr);
  if (filespace < 0) {
    return filespace;
  }
  std::string name = dsname;
  std::replace(name.begin(), name.end(), '/', '|');
  hid_t dset = H5Dcreate(groupID_, name.c_str(), typeId, filespace, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) {
    H5Sclose(filespace);
    return dset;
  }

  err = H5Dwrite(dset, memTypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (err < 0) {
    H5Dclose(dset);
    H5Sclose(filespace);
    return err;
  }

  err = H5Dclose(dset);
  if (err < 0) {
    return err;
  }
  err = H5Sclose(filespace);
  if (err < 0) {
    return err;
  }

  return 0;
}

herr_t Group::writeAttrArray(hid_t dsetID, const char *key,
                             const std::vector<int64_t> &data) const {
  hsize_t dims = data.size();
  hid_t dspace = H5Screate_simple(1, &dims, nullptr);
  if (dspace < 0) {
    return dspace;
  }
  hid_t attr =
      H5Acreate(dsetID, key, H5T_STD_I64LE, dspace, H5P_DEFAULT, H5P_DEFAULT);
  if (attr < 0) {
    // ignore error here, everything is lost anyway at this point in time
    H5Sclose(dspace);
    return attr;
  }
  herr_t err = H5Awrite(attr, H5T_NATIVE_LONG, data.data());
  if (err < 0) {
    return err;
  }
  err = H5Aclose(attr);
  if (err < 0) {
    return err;
  }
  err = H5Sclose(dspace);
  if (err < 0) {
    return err;
  }

  return 0;
}
} // namespace HDF5Utils
} // namespace OutputWriter
