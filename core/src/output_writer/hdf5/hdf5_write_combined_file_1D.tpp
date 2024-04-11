#include "output_writer/hdf5/hdf5.h"

#include "easylogging++.h"

#include "output_writer/paraview/loop_collect_mesh_properties.h"
#include "control/diagnostic_tool/performance_measurement.h"

namespace OutputWriter {

template <typename FieldVariablesForOutputWriterType>
void HDF5::writePolyDataFile(
    hid_t fileID, const FieldVariablesForOutputWriterType &fieldVariables,
    std::set<std::string> &meshNames) {
  bool meshPropertiesInitialized = !meshPropertiesPolyDataFile_.empty();
  std::vector<std::string> meshNamesVector;

  if (!meshPropertiesInitialized) {
    Control::PerformanceMeasurement::start("durationParaview1DInit");

    // collect the size data that is needed to compute offsets for parallel file
    // output
    ParaviewLoopOverTuple::loopCollectMeshProperties<
        FieldVariablesForOutputWriterType>(
        fieldVariables, meshPropertiesPolyDataFile_, meshNamesVector);

    Control::PerformanceMeasurement::stop("durationParaview1DInit");
  }

  VLOG(1) << "writePolyDataFile on rankSubset_: " << *this->rankSubset_;
  assert(this->rankSubset_);

  VLOG(1) << "meshPropertiesPolyDataFile_: " << meshPropertiesPolyDataFile_;

  if (!meshPropertiesInitialized) {
    Control::PerformanceMeasurement::start("durationParaview1DInit");

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
      if (!vtkPiece1D_.properties.pointDataArrays
               .empty()) // if properties are already assigned by an earlier
                         // mesh
      {
        if (vtkPiece1D_.properties.pointDataArrays.size() !=
            props.second.pointDataArrays.size()) {
          LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with "
                     << vtkPiece1D_.meshNamesCombinedMeshes
                     << ". Number of field variables mismatches for "
                     << meshName << " (is "
                     << props.second.pointDataArrays.size() << " instead of "
                     << vtkPiece1D_.properties.pointDataArrays.size() << ")";
          combineMesh = false;
        } else {
          for (int j = 0; j < props.second.pointDataArrays.size(); j++) {
            if (vtkPiece1D_.properties.pointDataArrays[j].name !=
                props.second.pointDataArrays[j]
                    .name) // if the name of the jth field variable is different
            {
              LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with "
                         << vtkPiece1D_.meshNamesCombinedMeshes
                         << ". Field variable names mismatch for " << meshName
                         << " (there is \""
                         << vtkPiece1D_.properties.pointDataArrays[j].name
                         << "\" instead of \""
                         << props.second.pointDataArrays[j].name << "\")";
              combineMesh = false;
            }
          }
        }

        if (combineMesh) {
          VLOG(1) << "Combine mesh " << meshName << " with "
                  << vtkPiece1D_.meshNamesCombinedMeshes << ", add "
                  << props.second.nPointsLocal << " points, "
                  << props.second.nCellsLocal << " elements to "
                  << vtkPiece1D_.properties.nPointsLocal << " points, "
                  << vtkPiece1D_.properties.nCellsLocal << " elements";

          vtkPiece1D_.properties.nPointsLocal += props.second.nPointsLocal;
          vtkPiece1D_.properties.nCellsLocal += props.second.nCellsLocal;

          vtkPiece1D_.properties.nPointsGlobal += props.second.nPointsGlobal;
          vtkPiece1D_.properties.nCellsGlobal += props.second.nCellsGlobal;
          vtkPiece1D_.setVTKValues();
        }
      } else {
        VLOG(1) << "this is the first 1D mesh";

        // properties are not yet assigned
        vtkPiece1D_.properties = props.second; // store properties
        vtkPiece1D_.setVTKValues();
      }

      VLOG(1) << "combineMesh: " << combineMesh;
      if (combineMesh) {
        vtkPiece1D_.meshNamesCombinedMeshes.insert(meshName);
      }
    }

    LOG(DEBUG) << "vtkPiece1D_: meshNamesCombinedMeshes: "
               << vtkPiece1D_.meshNamesCombinedMeshes
               << ", properties: " << vtkPiece1D_.properties
               << ", firstScalarName: " << vtkPiece1D_.firstScalarName
               << ", firstVectorName: " << vtkPiece1D_.firstVectorName;

    Control::PerformanceMeasurement::stop("durationParaview1DInit");
  }

  meshNames = vtkPiece1D_.meshNamesCombinedMeshes;

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

    vtkPiece1D_.properties.pointDataArrays.push_back(dataArrayName);
  }

  if (!meshPropertiesInitialized) {
    // exchange information about offset in terms of nCells and nPoints
    nCellsPreviousRanks1D_ = 0;
    nPointsPreviousRanks1D_ = 0;

    Control::PerformanceMeasurement::start("durationParaview1DInit");
    Control::PerformanceMeasurement::start("durationParaview1DReduction");
    MPIUtility::handleReturnValue(
        MPI_Exscan(&vtkPiece1D_.properties.nCellsLocal, &nCellsPreviousRanks1D_,
                   1, MPI_INT, MPI_SUM, this->rankSubset_->mpiCommunicator()),
        "MPI_Exscan");
    MPIUtility::handleReturnValue(
        MPI_Exscan(&vtkPiece1D_.properties.nPointsLocal,
                   &nPointsPreviousRanks1D_, 1, MPI_INT, MPI_SUM,
                   this->rankSubset_->mpiCommunicator()),
        "MPI_Exscan");
    Control::PerformanceMeasurement::stop("durationParaview1DReduction");
    Control::PerformanceMeasurement::stop("durationParaview1DInit");
  }

  // get local data values
  // setup connectivity array
  std::vector<int> connectivityValues(2 * vtkPiece1D_.properties.nCellsLocal);
  for (int i = 0; i < vtkPiece1D_.properties.nCellsLocal; i++) {
    connectivityValues[2 * i + 0] = nPointsPreviousRanks1D_ + i;
    connectivityValues[2 * i + 1] = nPointsPreviousRanks1D_ + i + 1;
  }

  // setup offset array
  std::vector<int> offsetValues(vtkPiece1D_.properties.nCellsLocal);
  for (int i = 0; i < vtkPiece1D_.properties.nCellsLocal; i++) {
    offsetValues[i] = 2 * nCellsPreviousRanks1D_ + 2 * i + 1;
  }

  // collect all data for the field variables, organized by field variable names
  std::map<std::string, std::vector<double>> fieldVariableValues;
  ParaviewLoopOverTuple::loopGetNodalValues<FieldVariablesForOutputWriterType>(
      fieldVariables, vtkPiece1D_.meshNamesCombinedMeshes, fieldVariableValues);

  assert(!fieldVariableValues.empty());
  fieldVariableValues["partitioning"].resize(
      vtkPiece1D_.properties.nPointsLocal,
      (double)this->rankSubset_->ownRankNo());

  // if next assertion will fail, output why for debugging
  if (fieldVariableValues.size() !=
      vtkPiece1D_.properties.pointDataArrays.size()) {
    LOG(DEBUG) << "n field variable values: " << fieldVariableValues.size()
               << ", n point data arrays: "
               << vtkPiece1D_.properties.pointDataArrays.size();
    LOG(DEBUG) << "vtkPiece1D_.meshNamesCombinedMeshes: "
               << vtkPiece1D_.meshNamesCombinedMeshes;
    std::stringstream pointDataArraysNames;
    for (const auto &pointDataArray : vtkPiece1D_.properties.pointDataArrays) {
      pointDataArraysNames << pointDataArray.name << " ";
    }
    LOG(DEBUG) << "pointDataArraysNames: " << pointDataArraysNames.str();
  }

  assert(fieldVariableValues.size() ==
         vtkPiece1D_.properties.pointDataArrays.size());

#ifndef NDEBUG
  LOG(DEBUG) << "fieldVariableValues: ";
  for (const std::pair<const std::string, std::vector<double>> &kv :
       fieldVariableValues) {
    LOG(DEBUG) << kv.first;
  }
#endif

  // check if field variable names have changed since last initialization
  for (const PolyDataPropertiesForMesh::DataArrayName &pointDataArray :
       vtkPiece1D_.properties.pointDataArrays) {
    LOG(DEBUG) << "  field variable \"" << pointDataArray.name << "\".";

    // if there is a field variable with a name that was not present when
    // vtkPiece1D_ was created
    if (fieldVariableValues.find(pointDataArray.name) ==
        fieldVariableValues.end()) {
      LOG(DEBUG) << "Field variable names have changed, reinitialize Paraview "
                    "output writer.";

      // reset now old variables
      meshPropertiesInitialized = false;
      meshPropertiesPolyDataFile_.clear();
      vtkPiece1D_ = VTKPiece();

      // recursively call this method
      writePolyDataFile(fileID, fieldVariables, meshNames);

      return;
    }
  }

  // collect all data for the geometry field variable
  std::vector<double> geometryFieldValues;
  ParaviewLoopOverTuple::loopGetGeometryFieldNodalValues<
      FieldVariablesForOutputWriterType>(
      fieldVariables, vtkPiece1D_.meshNamesCombinedMeshes, geometryFieldValues);

  // only continue if there is data to reduce
  if (vtkPiece1D_.meshNamesCombinedMeshes.empty()) {
    LOG(ERROR) << "There are no 1D meshes that could be combined, but Paraview "
                  "output with combineFiles=True was specified. \n(This only "
                  "works for 1D meshes.)";
  }

  LOG(DEBUG) << "Combined mesh from " << vtkPiece1D_.meshNamesCombinedMeshes;

  Control::PerformanceMeasurement::start("durationHDF51DWrite");

  hid_t groupID =
      H5Gcreate(fileID, "1D", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(groupID >= 0);

  // write field variables
  for (PolyDataPropertiesForMesh::DataArrayName &pointDataArray :
       vtkPiece1D_.properties.pointDataArrays) {
    assert(fieldVariableValues.find(pointDataArray.name) !=
           fieldVariableValues.end());

    // write values
    // for partitioning, convert float values to integer values for output
    bool writeFloatsAsInt = pointDataArray.name == "partitioning";
    writeCombinedValuesVector(groupID, fieldVariableValues[pointDataArray.name],
                              pointDataArray.name.c_str(), writeFloatsAsInt);
  }

  writeCombinedValuesVector(groupID, geometryFieldValues, "geometry");
  writeCombinedValuesVector(groupID, connectivityValues, "connectivity");
  writeCombinedValuesVector(groupID, offsetValues, "offsets");

  herr_t err = H5Gclose(groupID);
  assert(err >= 0);
  (void)err; // err is unsed in release mode, silence warning with this
  Control::PerformanceMeasurement::stop("durationHDF51DWrite");
}
} // namespace OutputWriter
