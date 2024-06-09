#include "output_writer/hdf5/hdf5.h"

#include "output_writer/loop_count_n_field_variables_of_mesh.h"
#include "output_writer/hdf5/loop_output.h"

namespace OutputWriter {
template <typename DataType>
void HDF5::write(DataType &data, const char *filename, int timeStepNo,
                 double currentTime, int callCountIncrement) {
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime,
                             callCountIncrement)) {
    return;
  }
  this->innerWrite(data.getFieldVariablesForOutputWriter(), filename,
                   timeStepNo, currentTime, callCountIncrement);
}

template <typename FieldVariablesForOutputWriterType>
void HDF5::innerWrite(const FieldVariablesForOutputWriterType &variable,
                      const char *filename, int timeStepNo, double currentTime,
                      int callCountIncrement) {
  Control::PerformanceMeasurement::start("durationHDF5Output");

  std::set<std::string> combined1DMeshes;
  std::set<std::string> combined2DMeshes;
  std::set<std::string> combined3DMeshes;

  if (combineFiles_) {
    std::unique_ptr<HDF5Utils::File> file;
    if (!filename) {
      // determine filename, broadcast from rank 0
      std::stringstream filename;
      filename << this->filenameBaseWithNo_ << "_c.h5";
      int filenameLength = filename.str().length();

      // broadcast length of filename
      MPIUtility::handleReturnValue(
          MPI_Bcast(&filenameLength, 1, MPI_INT, 0,
                    this->rankSubset_->mpiCommunicator()),
          "MPI_Bcast");

      std::vector<char> receiveBuffer(filenameLength + 1, char(0));
      strcpy(receiveBuffer.data(), filename.str().c_str());
      MPIUtility::handleReturnValue(
          MPI_Bcast(receiveBuffer.data(), filenameLength, MPI_CHAR, 0,
                    this->rankSubset_->mpiCommunicator()),
          "MPI_Bcast");
      file = std::make_unique<HDF5Utils::File>(receiveBuffer.data(), true);
    } else {
      file = std::make_unique<HDF5Utils::File>(filename, true);
    }

    herr_t err;
    if (writeMeta_) {
      err = file->writeAttrStr("version", DihuContext::versionText());
      assert(err >= 0);
      err = file->writeAttrStr("meta", DihuContext::metaText());
      assert(err >= 0);
    }
    err = file->writeAttrDouble("currentTime", this->currentTime_);
    assert(err >= 0);
    err = file->writeAttrInt("timeStepNo", this->timeStepNo_);
    assert(err >= 0);

    Control::PerformanceMeasurement::start("durationHDF51D");

    LOG(DEBUG) << "FieldVariablesForOutputWriterType: "
               << StringUtility::demangle(
                      typeid(FieldVariablesForOutputWriterType).name());

    // create a PolyData file that combines all 1D meshes into one file
    {
      HDF5Utils::Group group = file->newGroup("1D");
      writePolyDataFile<FieldVariablesForOutputWriterType>(group, variable,
                                                           combined1DMeshes);
    }

    Control::PerformanceMeasurement::stop("durationHDF51D");
    Control::PerformanceMeasurement::start("durationHDF53D");

    // create an UnstructuredMesh file that combines all 3D meshes into one file
    {
      HDF5Utils::Group group = file->newGroup("3D");
      writeCombinedUnstructuredGridFile<FieldVariablesForOutputWriterType>(
          group, variable, combined3DMeshes, true);
    }

    Control::PerformanceMeasurement::stop("durationHDF53D");
    Control::PerformanceMeasurement::start("durationHDF52D");

    // create an UnstructuredMesh file that combines all 2D meshes into one file
    {
      HDF5Utils::Group group = file->newGroup("2D");
      writeCombinedUnstructuredGridFile<FieldVariablesForOutputWriterType>(
          group, variable, combined2DMeshes, false);
    }

    Control::PerformanceMeasurement::stop("durationHDF52D");
  }

  // output normal files, parallel or if combineFiles_, only the 2D and 3D
  // meshes, combined

  // collect all available meshes
  std::set<std::string> meshNames;
  LoopOverTuple::loopCollectMeshNames<FieldVariablesForOutputWriterType>(
      variable, meshNames);

  // remove 1D meshes that were already output by writePolyDataFile
  std::set<std::string> meshesWithout1D;
  std::set_difference(meshNames.begin(), meshNames.end(),
                      combined1DMeshes.begin(), combined1DMeshes.end(),
                      std::inserter(meshesWithout1D, meshesWithout1D.end()));

  // remove 3D meshes that were already output by
  // writeCombinedUnstructuredGridFile
  std::set<std::string> meshesWithout1D3D;
  std::set_difference(
      meshesWithout1D.begin(), meshesWithout1D.end(), combined3DMeshes.begin(),
      combined3DMeshes.end(),
      std::inserter(meshesWithout1D3D, meshesWithout1D3D.end()));

  // remove 2D meshes that were already output by
  // writeCombinedUnstructuredGridFile
  std::set<std::string> meshesToOutput;
  std::set_difference(meshesWithout1D3D.begin(), meshesWithout1D3D.end(),
                      combined2DMeshes.begin(), combined2DMeshes.end(),
                      std::inserter(meshesToOutput, meshesToOutput.end()));

  if (meshesToOutput.size() > 0) {
    std::stringstream s;
    if (filename) {
      s << filename;
      if (combineFiles_) {
        s << "p";
      }
    } else {
      s << this->filename_ << "_p.h5";
    }

    HDF5Utils::File file = HDF5Utils::File(s.str().c_str(), false);
    herr_t err;
    if (writeMeta_) {
      err = file.writeAttrStr("version", DihuContext::versionText());
      assert(err >= 0);
      err = file.writeAttrStr("meta", DihuContext::metaText());
      assert(err >= 0);
    }
    err = file.writeAttrDouble("currentTime", this->currentTime_);
    assert(err >= 0);
    err = file.writeAttrInt("timeStepNo", this->timeStepNo_);
    assert(err >= 0);
    for (const std::string &meshName : meshesToOutput) {
      HDF5Utils::Group group = file.newGroup(meshName.c_str());
      // loop over all field variables and output those that are associated with
      // the mesh given by meshName
      HDF5LoopOverTuple::loopOutput(group, variable, variable, meshName,
                                    specificSettings_, currentTime);
    }
  }

  Control::PerformanceMeasurement::stop("durationHDF5Output");
}

namespace HDF5Utils {
template <typename T, size_t RANK>
herr_t Group::writeSimpleVec(const std::vector<T> &data,
                             const std::string &dsname,
                             const std::array<hsize_t, RANK> &dims) {
  if (file_->isMPIIO()) {
    if (std::is_same<T, int32_t>::value) {
      return writeVectorMPIIO(data.data(), dsname, dims, H5T_STD_I32LE,
                              H5T_NATIVE_INT, sizeof(int32_t));
    } else if (std::is_same<T, double>::value) {
      return writeVectorMPIIO(data.data(), dsname, dims, H5T_IEEE_F64LE,
                              H5T_NATIVE_DOUBLE, sizeof(double));
    }
  } else {
    if (std::is_same<T, int32_t>::value) {
      return writeVector(data.data(), dsname, dims, H5T_STD_I32LE,
                         H5T_NATIVE_INT, sizeof(int32_t));
    } else if (std::is_same<T, double>::value) {
      return writeVector(data.data(), dsname, dims, H5T_IEEE_F64LE,
                         H5T_NATIVE_DOUBLE, sizeof(double));
    }
  }
  return 0;
}

template <typename T>
herr_t Group::writeSimpleVec(const std::vector<T> &data,
                             const std::string &dsname) {
  std::array<hsize_t, 1> dims = {data.size()};
  return writeSimpleVec(data, dsname, dims);
}

template <size_t RANK>
herr_t Group::writeVectorMPIIO(const void *data, const std::string &dsname,
                               const std::array<hsize_t, RANK> &dims,
                               hid_t typeId, hid_t memTypeId, size_t dsize) {
  herr_t err;
  std::array<hsize_t, RANK> chunkDims = dims;

  std::pair<int, int> chunkAndOffsetSpace = file_->getMemSpace(chunkDims[0]);
  // data into set_chunk need to be the same size for everything
  chunkDims[0] = chunkAndOffsetSpace.second;
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
  chunkDims[0] = chunkAndOffsetSpace.first;
  std::array<hsize_t, RANK> count = {1};
  std::array<hsize_t, RANK> stride = {1};
  std::array<hsize_t, RANK> block = chunkDims;
  std::array<hsize_t, RANK> offset = {0};
  offset[0] = chunkAndOffsetSpace.second * file_->getOwnRank();

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
  void *dataOffset = ((char *)data) +
                     (chunkAndOffsetSpace.second * file_->getOwnRank() * dsize);
  err = H5Dwrite(dset, memTypeId, memspace, filespace, plist, dataOffset);
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

template <size_t RANK>
herr_t Group::writeVector(const void *data, const std::string &dsname,
                          const std::array<hsize_t, RANK> &dims, hid_t typeId,
                          hid_t memTypeId, size_t dsize) {
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

template <typename FieldVariableType>
herr_t writeFieldVariable(Group &group, FieldVariableType &fieldVariable) {
  LOG(DEBUG) << "HDF5 write field variable " << fieldVariable.name();
  VLOG(1) << fieldVariable;

  int nComponentsH5 = fieldVariable.nComponents();

  const int nComponents = FieldVariableType::nComponents();
  std::string stringData;

  std::vector<double> values;
  std::array<std::vector<double>, nComponents> componentValues;

  // ensure that ghost values are in place
  Partition::values_representation_t old_representation =
      fieldVariable.currentRepresentation();
  fieldVariable.zeroGhostBuffer();
  fieldVariable.setRepresentationGlobal();
  fieldVariable.startGhostManipulation();

  // get all local values including ghosts for the components
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    std::vector<double> retrievedLocalValues;
    fieldVariable.getValues(componentNo,
                            fieldVariable.functionSpace()
                                ->meshPartition()
                                ->dofNosLocalNaturalOrdering(),
                            retrievedLocalValues);

    const int nDofsPerNode = FieldVariableType::FunctionSpace::nDofsPerNode();
    const node_no_t nNodesLocal =
        fieldVariable.functionSpace()->meshPartition()->nNodesLocalWithGhosts();

    // for Hermite only extract the non-derivative values
    componentValues[componentNo].resize(nNodesLocal);

    int index = 0;
    for (int i = 0; i < nNodesLocal; i++) {
      componentValues[componentNo][i] = retrievedLocalValues[index];
      index += nDofsPerNode;
    }
  }

  // reset variable to old representation, so external code is not suprised
  // eg. time stepping code usally uses representationContiguous and will be
  // suprised if this changed we did not write to the values
  fieldVariable.setRepresentation(old_representation,
                                  values_modified_t::values_unchanged);

  // copy values in consecutive order (x y z x y z) to output
  values.reserve(componentValues[0].size() * nComponentsH5);
  for (int i = 0; i < componentValues[0].size(); i++) {
    for (int componentNo = 0; componentNo < nComponentsH5; componentNo++) {
      if (nComponents == 2 && nComponentsH5 == 3 && componentNo == 2) {
        values.push_back(0.0);
      } else {
        values.push_back(componentValues[componentNo][i]);
      }
    }
  }

  std::array<hsize_t, 2> dims = {nComponents, componentValues[0].size()};
  return group.writeSimpleVec(values, fieldVariable.name(), dims);
}

template <typename FieldVariableType>
herr_t writePartitionFieldVariable(Group &group,
                                   FieldVariableType &geometryField) {
  const node_no_t nNodesLocal =
      geometryField.functionSpace()->meshPartition()->nNodesLocalWithGhosts();

  std::vector<int32_t> values(nNodesLocal,
                              (int32_t)DihuContext::ownRankNoCommWorld());
  return group.writeSimpleVec(values, "partitioning");
}
} // namespace HDF5Utils
} // namespace OutputWriter
