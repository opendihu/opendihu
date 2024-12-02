#include "output_writer/json/json.h"

#include "output_writer/loop_count_n_field_variables_of_mesh.h"
#include "output_writer/json/loop_output.h"

namespace OutputWriter {
template <typename DataType>
void Json::write(DataType &data, const char *filename, int timeStepNo,
                 double currentTime, int callCountIncrement) {
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime,
                             callCountIncrement)) {
    return;
  }
  if (useCheckpointData_) {
    this->innerWrite(data.getFieldVariablesForCheckpointing(), filename,
                     timeStepNo, currentTime, callCountIncrement);
  } else {
    this->innerWrite(data.getFieldVariablesForOutputWriter(), filename,
                     timeStepNo, currentTime, callCountIncrement);
  }
}

template <typename FieldVariablesForOutputWriterType>
void Json::innerWrite(const FieldVariablesForOutputWriterType &variables,
                      const char *filename, int timeStepNo, double currentTime,
                      int callCountIncrement) {
  Control::PerformanceMeasurement::start("durationJsonOutput");

  std::set<std::string> combined1DMeshes;
  std::set<std::string> combined2DMeshes;
  std::set<std::string> combined3DMeshes;

  std::unique_ptr<JsonUtils::File> file;
  {
    if (!filename) {
      std::string finalFileName;
      if (combineFiles_) {
        // determine filename, broadcast from rank 0
        std::stringstream filename;
        filename << this->filenameBaseWithNo_ << "_c.json";
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

        finalFileName = std::string(receiveBuffer.data());
      } else {
        finalFileName = std::string(this->filename_ + "_p.h5");
      }
      file = std::make_unique<JsonUtils::File>(finalFileName.c_str(),
                                               combineFiles_);
    } else {
      file = std::make_unique<JsonUtils::File>(filename, combineFiles_);
    }

    if (writeMeta_) {
      file->writeAttr("rawVersion", DihuContext::version());
      file->writeAttr("version", DihuContext::versionText());
      file->writeAttr("meta", DihuContext::metaText());
      file->writeAttr("worldSize", DihuContext::nRanksCommWorld());
    }
    file->writeAttr("currentTime", this->currentTime_);
    file->template writeAttr<int32_t>("timeStepNo", this->timeStepNo_);

    Control::PerformanceMeasurement::start("durationJson1D");

    if (useCheckpointData_) {
      LOG(DEBUG) << "FieldVariablesForCheckpointing: "
                 << StringUtility::demangle(
                        typeid(FieldVariablesForOutputWriterType).name());
    } else {
      LOG(DEBUG) << "FieldVariablesForOutputWriter: "
                 << StringUtility::demangle(
                        typeid(FieldVariablesForOutputWriterType).name());
    }

    // create a PolyData file that combines all 1D meshes into one file
    {
      JsonUtils::Group group = file->newGroup("1D");
      writePolyDataFile<FieldVariablesForOutputWriterType>(
          group, variables, combined1DMeshes, useCheckpointData_);
    }

    Control::PerformanceMeasurement::stop("durationJson1D");
    Control::PerformanceMeasurement::start("durationJson3D");

    // create an UnstructuredMesh file that combines all 3D meshes into one file
    {
      JsonUtils::Group group = file->newGroup("3D");
      writeCombinedUnstructuredGridFile<FieldVariablesForOutputWriterType>(
          group, variables, combined3DMeshes, true, useCheckpointData_);
    }

    Control::PerformanceMeasurement::stop("durationJson3D");
    Control::PerformanceMeasurement::start("durationJson2D");

    // create an UnstructuredMesh file that combines all 2D meshes into one file
    {
      JsonUtils::Group group = file->newGroup("2D");
      writeCombinedUnstructuredGridFile<FieldVariablesForOutputWriterType>(
          group, variables, combined2DMeshes, false, useCheckpointData_);
    }

    Control::PerformanceMeasurement::stop("durationJson2D");
  }

  // if we have a combined file, we are done with it
  if (combineFiles_) {
    file.reset();
  }

  // output normal files, parallel or if combineFiles_, only the 2D and 3D
  // meshes, combined

  // collect all available meshes
  std::set<std::string> meshNames;
  LoopOverTuple::loopCollectMeshNames<FieldVariablesForOutputWriterType>(
      variables, meshNames);

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
      s << this->filename_ << "_p.json";
    }

    // if we have a combine file we create a new file, otherwise we can reuse
    // the original file
    if (combineFiles_) {
      file = std::make_unique<JsonUtils::File>(s.str().c_str(), false);
      if (writeMeta_) {
        file->writeAttr("rawVersion", DihuContext::version());
        file->writeAttr("version", DihuContext::versionText());
        file->writeAttr("meta", DihuContext::metaText());
        file->writeAttr("worldSize", DihuContext::nRanksCommWorld());
      }
      file->writeAttr("currentTime", this->currentTime_);
      file->template writeAttr<int32_t>("timeStepNo", this->timeStepNo_);
    }

    for (const std::string &meshName : meshesToOutput) {
      JsonUtils::Group group = file->newGroup(meshName.c_str());
      // loop over all field variables and output those that are associated with
      // the mesh given by meshName
      JsonLoopOverTuple::loopOutput(group, variables, variables, meshName,
                                    specificSettings_, currentTime);
    }
  }

  Control::PerformanceMeasurement::stop("durationJsonOutput");
}

namespace JsonUtils {
template <typename T, std::enable_if_t<std::is_same<T, int32_t>::value, bool>>
void File::writeAttr(const char *key, const T &v) {
  content_["__attributes"][key] = v;
}

template <typename T, std::enable_if_t<std::is_same<T, double>::value, bool>>
void File::writeAttr(const char *key, const T &v) {
  content_["__attributes"][key] = v;
}

template <typename T,
          std::enable_if_t<std::is_same<T, std::string>::value, bool>>
void File::writeAttr(const char *key, const T &v) {
  content_["__attributes"][key] = v;
}

template <typename T>
void Group::writeSimpleVec(const std::vector<T> &data,
                           const std::string &dsname) {
  if (file_->isMPIIO()) {
    return writeVectorMPIIO(data, dsname);
  } else {
    return writeVector(data, dsname);
  }
}

template <typename T>
void Group::writeVectorMPIIO(const std::vector<T> &data,
                             const std::string &dsname) {
  MPI_Datatype sendReceiveType = MPI_DATATYPE_NULL;
  if (std::is_same<T, int>::value || std::is_same<T, int32_t>::value) {
    sendReceiveType = MPI_INT;
  } else if (std::is_same<T, double>::value) {
    sendReceiveType = MPI_DOUBLE;
  }

  std::string name = dsname;
  std::replace(name.begin(), name.end(), '/', '|');

  if (sendReceiveType == MPI_DATATYPE_NULL) {
    LOG(ERROR) << "writeVectorMPIIO was called with unsupported T: currently "
                  "only supports int32_t and double";
    return;
  }

  std::vector<int32_t> sizes, displacements;
  sizes.resize(file_->getWorldSize());
  displacements.resize(file_->getWorldSize());

  const int32_t sendSize = (int32_t)data.size();
  MPIUtility::handleReturnValue(MPI_Allgather(&sendSize, 1, MPI_INT,
                                              sizes.data(), 1, MPI_INT,
                                              MPI_COMM_WORLD),
                                "MPI_Allgather");
  displacements[0] = 0;
  for (size_t i = 1; i < file_->getWorldSize(); i++) {
    displacements[i] = displacements[i - 1] + sizes[i - 1];
  }

  std::vector<T> recv;
  recv.resize(std::accumulate(sizes.begin(), sizes.end(), 0ul));
  MPIUtility::handleReturnValue(
      MPI_Allgatherv(data.data(), data.size(), sendReceiveType, recv.data(),
                     sizes.data(), displacements.data(), sendReceiveType,
                     MPI_COMM_WORLD),
      "MPI_Allgatherv");

  // write combined data with metadata
  groupContent_[name]["__data"] = recv;
  groupContent_[name]["__attributes"]["chunkDims"] = sizes;
}

template <typename T>
void Group::writeVector(const std::vector<T> &data, const std::string &dsname) {
  std::string name = dsname;
  std::replace(name.begin(), name.end(), '/', '|');

  groupContent_[name]["__data"] = data;
}

template <typename FieldVariableType>
void writeFieldVariable(Group &group, FieldVariableType &fieldVariable) {
  LOG(DEBUG) << "Json write field variable " << fieldVariable.name();
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

  return group.writeSimpleVec(values, fieldVariable.uniqueName());
}

template <typename FieldVariableType>
void writePartitionFieldVariable(Group &group,
                                 FieldVariableType &geometryField) {
  const node_no_t nNodesLocal =
      geometryField.functionSpace()->meshPartition()->nNodesLocalWithGhosts();

  std::vector<int32_t> values(nNodesLocal,
                              (int32_t)DihuContext::ownRankNoCommWorld());
  group.writeSimpleVec(values, "partitioning");
}
} // namespace JsonUtils
} // namespace OutputWriter
