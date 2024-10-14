#include "output_writer/hdf5/hdf5.h"

#include "easylogging++.h"

#include "output_writer/loop_collect_mesh_properties.h"
#include "output_writer/loop_get_nodal_values.h"
#include "output_writer/loop_get_geometry_field_nodal_values.h"
#include "control/diagnostic_tool/performance_measurement.h"

namespace OutputWriter {

template <typename FieldVariablesForOutputWriterType>
void HDF5::writePolyDataFile(
    HDF5Utils::Group &group,
    const FieldVariablesForOutputWriterType &fieldVariables,
    std::set<std::string> &meshNames) {
  bool meshPropertiesInitialized = !meshPropertiesPolyDataFile_.empty();
  std::vector<std::string> meshNamesVector;

  if (!meshPropertiesInitialized) {
    Control::PerformanceMeasurement::start("durationHDF51DInit");

    // collect the size data that is needed to compute offsets for parallel file
    // output
    LoopOverTuple::loopCollectMeshProperties<FieldVariablesForOutputWriterType>(
        fieldVariables, meshPropertiesPolyDataFile_, meshNamesVector);

    Control::PerformanceMeasurement::stop("durationHDF51DInit");
  }

  VLOG(1) << "writePolyDataFile on rankSubset_: " << *this->rankSubset_;
  assert(this->rankSubset_);

  VLOG(1) << "meshPropertiesPolyDataFile_: " << meshPropertiesPolyDataFile_;

  if (!meshPropertiesInitialized) {
    Control::PerformanceMeasurement::start("durationHDF51DInit");

    // parse the collected properties of the meshes that will be output to the
    // file
    for (const std::pair<const std::string, PolyDataPropertiesForMesh> &props :
         meshPropertiesPolyDataFile_) {
      std::string meshName = props.first;

      // do not combine meshes other than 1D meshes
      if (props.second.dimensionality != 1) {
        continue;
      }

      // check if this mesh should be combined with other meshes
      bool combineMesh = true;

      // check if mesh can be merged into previous meshes
      if (!piece1D_.properties.pointDataArrays
               .empty()) // if properties are already assigned by an earlier
                         // mesh
      {
        if (piece1D_.properties.pointDataArrays.size() !=
            props.second.pointDataArrays.size()) {
          LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with "
                     << piece1D_.meshNamesCombinedMeshes
                     << ". Number of field variables mismatches for "
                     << meshName << " (is "
                     << props.second.pointDataArrays.size() << " instead of "
                     << piece1D_.properties.pointDataArrays.size() << ")";
          combineMesh = false;
        } else {
          for (int j = 0; j < props.second.pointDataArrays.size(); j++) {
            if (piece1D_.properties.pointDataArrays[j].name !=
                props.second.pointDataArrays[j]
                    .name) // if the name of the jth field variable is different
            {
              LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with "
                         << piece1D_.meshNamesCombinedMeshes
                         << ". Field variable names mismatch for " << meshName
                         << " (there is \""
                         << piece1D_.properties.pointDataArrays[j].name
                         << "\" instead of \""
                         << props.second.pointDataArrays[j].name << "\")";
              combineMesh = false;
            }
          }
        }

        if (combineMesh) {
          VLOG(1) << "Combine mesh " << meshName << " with "
                  << piece1D_.meshNamesCombinedMeshes << ", add "
                  << props.second.nPointsLocal << " points, "
                  << props.second.nCellsLocal << " elements to "
                  << piece1D_.properties.nPointsLocal << " points, "
                  << piece1D_.properties.nCellsLocal << " elements";

          piece1D_.properties.nPointsLocal += props.second.nPointsLocal;
          piece1D_.properties.nCellsLocal += props.second.nCellsLocal;

          piece1D_.properties.nPointsGlobal += props.second.nPointsGlobal;
          piece1D_.properties.nCellsGlobal += props.second.nCellsGlobal;
          piece1D_.setVTKValues();
        }
      } else {
        VLOG(1) << "this is the first 1D mesh";

        // properties are not yet assigned
        piece1D_.properties = props.second; // store properties
        piece1D_.setVTKValues();
      }

      VLOG(1) << "combineMesh: " << combineMesh;
      if (combineMesh) {
        piece1D_.meshNamesCombinedMeshes.insert(meshName);
      }
    }

    LOG(DEBUG) << "piece1D_: meshNamesCombinedMeshes: "
               << piece1D_.meshNamesCombinedMeshes
               << ", properties: " << piece1D_.properties
               << ", firstScalarName: " << piece1D_.firstScalarName
               << ", firstVectorName: " << piece1D_.firstVectorName;

    Control::PerformanceMeasurement::stop("durationHDF51DInit");
  }

  meshNames = piece1D_.meshNamesCombinedMeshes;

  // if there are no 1D meshes, return
  if (meshNames.empty()) {
    return;
  }

  if (!meshPropertiesInitialized) {
    // add field variable "partitioning" with 1 component
    PolyDataPropertiesForMesh::DataArrayName dataArrayName;
    dataArrayName.name = "partitioning";
    dataArrayName.nComponents = 1;
    dataArrayName.componentNames = std::vector<std::string>(1, "rankNo");

    piece1D_.properties.pointDataArrays.push_back(dataArrayName);
  }

  if (!meshPropertiesInitialized) {
    // exchange information about offset in terms of nCells and nPoints
    nCellsPreviousRanks1D_ = 0;
    nPointsPreviousRanks1D_ = 0;

    Control::PerformanceMeasurement::start("durationHDF51DInit");
    Control::PerformanceMeasurement::start("durationHDF51DReduction");
    MPIUtility::handleReturnValue(
        MPI_Exscan(&piece1D_.properties.nCellsLocal, &nCellsPreviousRanks1D_, 1,
                   MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()),
        "MPI_Exscan");
    MPIUtility::handleReturnValue(
        MPI_Exscan(&piece1D_.properties.nPointsLocal, &nPointsPreviousRanks1D_,
                   1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()),
        "MPI_Exscan");
    Control::PerformanceMeasurement::stop("durationHDF51DReduction");
    Control::PerformanceMeasurement::stop("durationHDF51DInit");
  }

  // get local data values
  // setup connectivity array
  std::vector<int> connectivityValues;
  connectivityValues.reserve(2 * piece1D_.properties.nCellsLocal);
  for (int i = 0; i < piece1D_.properties.nCellsLocal; i++) {
    connectivityValues[2 * i + 0] = nPointsPreviousRanks1D_ + i;
    connectivityValues[2 * i + 1] = nPointsPreviousRanks1D_ + i + 1;
  }

  // setup offset array
  std::vector<int> offsetValues;
  offsetValues.reserve(piece1D_.properties.nCellsLocal);
  for (int i = 0; i < piece1D_.properties.nCellsLocal; i++) {
    offsetValues[i] = 2 * nCellsPreviousRanks1D_ + 2 * i + 1;
  }

  // collect all data for the field variables, organized by field variable names
  std::map<std::string, std::vector<double>> fieldVariableValues;
  LoopOverTuple::loopGetNodalValues<FieldVariablesForOutputWriterType>(
      fieldVariables, piece1D_.meshNamesCombinedMeshes, fieldVariableValues);

  assert(!fieldVariableValues.empty());
  fieldVariableValues["partitioning"].resize(
      piece1D_.properties.nPointsLocal, (double)this->rankSubset_->ownRankNo());

  // if next assertion will fail, output why for debugging
  if (fieldVariableValues.size() !=
      piece1D_.properties.pointDataArrays.size()) {
    LOG(DEBUG) << "n field variable values: " << fieldVariableValues.size()
               << ", n point data arrays: "
               << piece1D_.properties.pointDataArrays.size();
    LOG(DEBUG) << "piece1D_.meshNamesCombinedMeshes: "
               << piece1D_.meshNamesCombinedMeshes;
    std::stringstream pointDataArraysNames;
    for (const auto &pointDataArray : piece1D_.properties.pointDataArrays) {
      pointDataArraysNames << pointDataArray.name << " ";
    }
    LOG(DEBUG) << "pointDataArraysNames: " << pointDataArraysNames.str();
  }

  assert(fieldVariableValues.size() ==
         piece1D_.properties.pointDataArrays.size());

#ifndef NDEBUG
  LOG(DEBUG) << "fieldVariableValues: ";
  for (const std::pair<const std::string, std::vector<double>> &kv :
       fieldVariableValues) {
    LOG(DEBUG) << kv.first;
  }
#endif

  // check if field variable names have changed since last initialization
  for (const PolyDataPropertiesForMesh::DataArrayName &pointDataArray :
       piece1D_.properties.pointDataArrays) {
    LOG(DEBUG) << "  field variable \"" << pointDataArray.name << "\".";

    // if there is a field variable with a name that was not present when
    // piece1D_ was created
    if (fieldVariableValues.find(pointDataArray.name) ==
        fieldVariableValues.end()) {
      LOG(DEBUG) << "Field variable names have changed, reinitialize HDF5 "
                    "output writer.";

      // reset now old variables
      meshPropertiesInitialized = false;
      meshPropertiesPolyDataFile_.clear();
      piece1D_ = Piece();

      // recursively call this method
      writePolyDataFile(group, fieldVariables, meshNames);

      return;
    }
  }

  // collect all data for the geometry field variable
  std::vector<double> geometryFieldValues;
  LoopOverTuple::loopGetGeometryFieldNodalValues<
      FieldVariablesForOutputWriterType>(
      fieldVariables, piece1D_.meshNamesCombinedMeshes, geometryFieldValues);

  // only continue if there is data to reduce
  if (piece1D_.meshNamesCombinedMeshes.empty()) {
    LOG(ERROR) << "There are no 1D meshes that could be combined, but HDF5 "
                  "output with combineFiles=True was specified. \n(This only "
                  "works for 1D meshes.)";
  }

  LOG(DEBUG) << "Combined mesh from " << piece1D_.meshNamesCombinedMeshes;

  Control::PerformanceMeasurement::start("durationHDF51DWrite");

  herr_t err;
  // write field variables
  for (PolyDataPropertiesForMesh::DataArrayName &pointDataArray :
       piece1D_.properties.pointDataArrays) {
    assert(fieldVariableValues.find(pointDataArray.name) !=
           fieldVariableValues.end());

    const std::vector<double> &dataToWrite =
        fieldVariableValues[pointDataArray.name];
    // write values
    // for partitioning, convert float values to integer values for output
    if (pointDataArray.name == "partitioning") {
      std::vector<int32_t> newVals(dataToWrite.size());

      for (const auto &v : dataToWrite) {
        newVals.push_back((int32_t)(round(v)));
      }
      err = group.writeSimpleVec<int32_t>(newVals, pointDataArray.name.c_str());
    } else {
      err = group.writeSimpleVec<double>(dataToWrite,
                                         pointDataArray.name.c_str());
    }
    assert(err >= 0);
  }

  err = group.writeSimpleVec<double>(geometryFieldValues, "geometry");
  assert(err >= 0);
  err = group.writeSimpleVec<int32_t>(connectivityValues, "connectivity");
  assert(err >= 0);
  err = group.writeSimpleVec<int32_t>(offsetValues, "offsets");
  assert(err >= 0);
  Control::PerformanceMeasurement::stop("durationHDF51DWrite");
}
} // namespace OutputWriter
