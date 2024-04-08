#include "output_writer/hdf5/hdf5.h"

namespace OutputWriter {
HDF5::HDF5(DihuContext context, PythonConfig settings,
           std::shared_ptr<Partition::RankSubset> rankSubset)
    : Generic(context, settings, rankSubset) {
  combineFiles_ = settings.getOptionBool("combineFiles", false);
}

//! constructor, initialize nPoints and nCells to 0
HDF5::VTKPiece::VTKPiece() {
  properties.nPointsLocal = 0;
  properties.nCellsLocal = 0;
  properties.nPointsGlobal = 0;
  properties.nCellsGlobal = 0;
  properties.dimensionality = 0;
}

//! assign the correct values to firstScalarName and firstVectorName, only if
//! properties has been set
void HDF5::VTKPiece::setVTKValues() {
  // set values for firstScalarName and firstVectorName from the values in
  // pointDataArrays
  for (auto pointDataArray : properties.pointDataArrays) {
    if (firstScalarName == "" && pointDataArray.nComponents != 3)
      firstScalarName = pointDataArray.name;
    if (firstVectorName == "" && pointDataArray.nComponents == 3)
      firstVectorName = pointDataArray.name;
  }
}

void HDF5::writeCombinedTypesVector(hid_t fileID, int nValues,
                                    bool output3DMeshes, const char *dsname) {
  if (output3DMeshes) {
    std::vector<int> values(nValues, 12);
    writeCombinedValuesVector(fileID, values, dsname);
  } else {
    std::vector<int> values(nValues, 9);
    writeCombinedValuesVector(fileID, values, dsname);
  }
}

hid_t HDF5::openHDF5File(const char *filename, bool mpiio) {
  if (mpiio) {
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    herr_t err;
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

herr_t HDF5::writeAttrInt(hid_t fileID, const char *key, int32_t value) {
  std::array<hsize_t, 1> dims = {1};
  std::array<int32_t, 1> data = {value};
  hid_t dspace = H5Screate_simple(1, dims.data(), nullptr);
  if (dspace < 0) {
    return dspace;
  }
  hid_t attr =
      H5Acreate(fileID, key, H5T_STD_I32BE, dspace, H5P_DEFAULT, H5P_DEFAULT);
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

herr_t HDF5::writeAttrDouble(hid_t fileID, const char *key, double value) {
  std::array<hsize_t, 1> dims = {1};
  std::array<double, 1> data = {value};
  hid_t dspace = H5Screate_simple(1, dims.data(), nullptr);
  if (dspace < 0) {
    return dspace;
  }
  hid_t attr =
      H5Acreate(fileID, key, H5T_IEEE_F64BE, dspace, H5P_DEFAULT, H5P_DEFAULT);
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

herr_t HDF5::writeAttrString(hid_t fileID, const char *key,
                             const std::string &value) {
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

  std::array<hsize_t, 1> dims = {1};
  hid_t dspace = H5Screate_simple(1, dims.data(), nullptr);
  if (dspace < 0) {
    return dspace;
  }
  hid_t attr =
      H5Acreate(fileID, key, filetype, dspace, H5P_DEFAULT, H5P_DEFAULT);
  if (attr < 0) {
    // ignore error here, everything is lost anyway at this point in time
    H5Sclose(dspace);
    return attr;
  }
  err = H5Awrite(attr, memtype, value.c_str());
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

} // namespace OutputWriter
