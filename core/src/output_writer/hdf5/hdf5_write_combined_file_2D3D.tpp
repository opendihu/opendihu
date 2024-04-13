#include "output_writer/hdf5/hdf5.h"

#include "easylogging++.h"

#include "output_writer/loop_collect_mesh_properties.h"
#include "output_writer/loop_get_nodal_values.h"
#include "output_writer/loop_get_geometry_field_nodal_values.h"
#include "control/diagnostic_tool/performance_measurement.h"

namespace OutputWriter {

template <typename FieldVariablesForOutputWriterType>
void HDF5::writeCombinedUnstructuredGridFile(
    hid_t fileID, const FieldVariablesForOutputWriterType &fieldVariables,
    std::set<std::string> &meshNames, bool output3DMeshes) {
  std::map<std::string, PolyDataPropertiesForMesh>
      &meshPropertiesUnstructuredGridFile_ =
          (output3DMeshes ? meshPropertiesUnstructuredGridFile3D_
                          : meshPropertiesUnstructuredGridFile2D_);

  bool meshPropertiesInitialized = !meshPropertiesUnstructuredGridFile_.empty();
  std::vector<std::string> meshNamesVector;

  if (!meshPropertiesInitialized) {
    Control::PerformanceMeasurement::start("durationHDF53DInit");

    // collect the size data that is needed to compute offsets for parallel file
    // output
    LoopOverTuple::loopCollectMeshProperties<FieldVariablesForOutputWriterType>(
        fieldVariables, meshPropertiesUnstructuredGridFile_, meshNamesVector);

    Control::PerformanceMeasurement::stop("durationHDF53DInit");
  }

  VLOG(1) << "writeCombinedUnstructuredGridFile on rankSubset_: "
          << *this->rankSubset_;
  assert(this->rankSubset_);

  VLOG(1) << "meshPropertiesUnstructuredGridFile_: "
          << meshPropertiesUnstructuredGridFile_;

  const char *targetDimStr = "2D";
  int targetDimensionality = 2;
  if (output3DMeshes) {
    targetDimStr = "3D";
    targetDimensionality = 3;
  }

  VLOG(1) << "targetDimensionality: " << targetDimensionality;

  if (!meshPropertiesInitialized) {
    Control::PerformanceMeasurement::start("durationHDF52D3DInit");

    // parse the collected properties of the meshes that will be output to the
    // file
    for (const std::string &meshName : meshNamesVector) {
      PolyDataPropertiesForMesh &polyDataPropertiesForMesh =
          meshPropertiesPolyDataFile_[meshName];

      // do not combine meshes other than 1D meshes
      if (polyDataPropertiesForMesh.dimensionality != targetDimensionality) {
        continue;
      }

      // check if this mesh should be combined with other meshes
      bool combineMesh = true;

      // check if mesh can be merged into previous meshes
      if (!piece3D_.properties.pointDataArrays
               .empty()) // if properties are already assigned by an earlier
                         // mesh
      {
        if (piece3D_.properties.pointDataArrays.size() !=
            polyDataPropertiesForMesh.pointDataArrays.size()) {
          LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with "
                     << piece3D_.meshNamesCombinedMeshes
                     << ". Number of field variables mismatches for "
                     << meshName << " (is "
                     << polyDataPropertiesForMesh.pointDataArrays.size()
                     << " instead of "
                     << piece3D_.properties.pointDataArrays.size() << ")";
          combineMesh = false;
        } else {
          for (int j = 0; j < polyDataPropertiesForMesh.pointDataArrays.size();
               j++) {
            if (piece3D_.properties.pointDataArrays[j].name !=
                polyDataPropertiesForMesh.pointDataArrays[j]
                    .name) // if the name of the jth field variable is different
            {
              LOG(DEBUG) << "Mesh " << meshName << " cannot be combined with "
                         << piece3D_.meshNamesCombinedMeshes
                         << ". Field variable names mismatch for " << meshName
                         << " (there is \""
                         << piece3D_.properties.pointDataArrays[j].name
                         << "\" instead of \""
                         << polyDataPropertiesForMesh.pointDataArrays[j].name
                         << "\")";
              combineMesh = false;
            }
          }
        }

        if (combineMesh) {
          VLOG(1) << "Combine mesh " << meshName << " with "
                  << piece3D_.meshNamesCombinedMeshes << ", add "
                  << polyDataPropertiesForMesh.nPointsLocal << " points, "
                  << polyDataPropertiesForMesh.nCellsLocal << " elements to "
                  << piece3D_.properties.nPointsLocal << " points, "
                  << piece3D_.properties.nCellsLocal << " elements";

          piece3D_.properties.nPointsLocal +=
              polyDataPropertiesForMesh.nPointsLocal;
          piece3D_.properties.nCellsLocal +=
              polyDataPropertiesForMesh.nCellsLocal;

          piece3D_.properties.nPointsGlobal +=
              polyDataPropertiesForMesh.nPointsGlobal;
          piece3D_.properties.nCellsGlobal +=
              polyDataPropertiesForMesh.nCellsGlobal;
          piece3D_.setVTKValues();
        }
      } else {
        VLOG(1) << "this is the first " << targetDimensionality << "D mesh";

        // properties are not yet assigned
        piece3D_.properties = polyDataPropertiesForMesh; // store properties
        piece3D_.setVTKValues();
      }

      VLOG(1) << "combineMesh: " << combineMesh;
      if (combineMesh) {
        // if the mesh is not yet present in meshNamesCombinedMeshes, add it
        if (piece3D_.meshNamesCombinedMeshes.find(meshName) ==
            piece3D_.meshNamesCombinedMeshes.end()) {
          piece3D_.meshNamesCombinedMeshes.insert(meshName);
          piece3D_.meshNamesCombinedMeshesVector.push_back(meshName);
        }
      }
    }

    LOG(DEBUG) << "piece3D_: meshNamesCombinedMeshes: "
               << piece3D_.meshNamesCombinedMeshes
               << ", properties: " << piece3D_.properties
               << ", firstScalarName: " << piece3D_.firstScalarName
               << ", firstVectorName: " << piece3D_.firstVectorName;

    Control::PerformanceMeasurement::stop("durationHDF52D3DInit");
  }

  // write combined mesh
  if (!piece3D_.meshNamesCombinedMeshes.empty()) {
    const PolyDataPropertiesForMesh &polyDataPropertiesForMesh =
        piece3D_.properties;

    VLOG(1) << "polyDataPropertiesForMesh combined: "
            << polyDataPropertiesForMesh;

    // if the meshes have the right dimensionality (2D or 3D)
    if (piece3D_.properties.dimensionality == targetDimensionality) {
      meshNames = piece3D_.meshNamesCombinedMeshes;

      hid_t groupID = H5Gcreate(fileID, targetDimStr, H5P_DEFAULT, H5P_DEFAULT,
                                H5P_DEFAULT);
      assert(groupID >= 0);

      // write actual file
      writeCombinedUnstructuredGridFile<FieldVariablesForOutputWriterType>(
          groupID, fieldVariables, piece3D_.properties,
          meshPropertiesUnstructuredGridFile_,
          piece3D_.meshNamesCombinedMeshesVector, meshPropertiesInitialized);
      herr_t err = H5Gclose(groupID);
      assert(err >= 0);
      (void)err; // err is unsed in release mode, silence warning with this
    } else {
      VLOG(1) << "skip meshes " << piece3D_.meshNamesCombinedMeshes
              << " because " << piece3D_.properties.dimensionality
              << " != " << targetDimensionality;
    }
  }

  // loop over 3D or 2D meshes
  for (std::pair<const std::string, PolyDataPropertiesForMesh> &props :
       meshPropertiesUnstructuredGridFile_) {
    PolyDataPropertiesForMesh &polyDataPropertiesForMesh = props.second;

    // if the current mesh does not have the correct dimensionality, skip this
    // mesh
    if (polyDataPropertiesForMesh.dimensionality != targetDimensionality) {
      continue;
    }

    // if the current mesh was already output within the combined file, skip
    // this mesh
    if (meshNames.find(props.first) != meshNames.end()) {
      continue;
    }

    LOG(DEBUG) << "meshNames: " << meshNames
               << ", next mesh to write: " << props.first;
    hid_t groupID = H5Gcreate(fileID, props.first.c_str(), H5P_DEFAULT,
                              H5P_DEFAULT, H5P_DEFAULT);
    assert(groupID >= 0);

    // write actual file
    std::vector<std::string> currentMesh;
    currentMesh.push_back(props.first);
    meshNames.insert(props.first);

    writeCombinedUnstructuredGridFile<FieldVariablesForOutputWriterType>(
        groupID, fieldVariables, polyDataPropertiesForMesh,
        meshPropertiesUnstructuredGridFile_, currentMesh,
        meshPropertiesInitialized);
    herr_t err = H5Gclose(groupID);
    assert(err >= 0);
    (void)err; // err is unsed in release mode, silence warning with this

    // for the next meshes, reinitialize, the initialized information can only
    // be shared if there is exactly 1 mesh to write in this method
    meshPropertiesInitialized = false;
  }
}

template <typename FieldVariablesForOutputWriterType>
void HDF5::writeCombinedUnstructuredGridFile(
    hid_t fileID, const FieldVariablesForOutputWriterType &fieldVariables,
    PolyDataPropertiesForMesh &polyDataPropertiesForMesh,
    const std::map<std::string, PolyDataPropertiesForMesh>
        &meshPropertiesUnstructuredGridFile,
    std::vector<std::string> meshNames, bool meshPropertiesInitialized) {
  int targetDimensionality = polyDataPropertiesForMesh.dimensionality;
  bool output3DMeshes = targetDimensionality == 3;

  std::set<std::string> meshNamesSet(meshNames.begin(), meshNames.end());

  VLOG(1) << "writeCombinedUnstructuredGridFile, fileID=" << fileID
          << ", meshNames: " << meshNames
          << ", meshPropertiesInitialized=" << meshPropertiesInitialized
          << ", targetDimensionality: " << targetDimensionality;

  //! find out name of first scalar and vector field variables
  std::string firstScalarName;
  std::string firstVectorName;

  for (const PolyDataPropertiesForMesh::DataArrayName &pointDataArray :
       polyDataPropertiesForMesh.pointDataArrays) {
    if (pointDataArray.nComponents == 3 && firstVectorName.empty()) {
      firstVectorName = pointDataArray.name;
    } else if (firstScalarName.empty()) {
      firstScalarName = pointDataArray.name;
    }

    if (!firstScalarName.empty() && firstVectorName.empty())
      break;
  }

  // exchange information about offset in terms of nCells and nPoints
  if (!meshPropertiesInitialized) {
    VLOG(1) << "set nCellsPreviousRanks3D_, nPointsPreviousRanks3D_";

    Control::PerformanceMeasurement::start("durationHDF53DInit");
    Control::PerformanceMeasurement::start("durationHDF53DReduction");
    for (const std::string &meshName : meshNames) {
      nCellsPreviousRanks3D_[meshName] = 0;
      nPointsPreviousRanks3D_[meshName] = 0;

      MPIUtility::handleReturnValue(
          MPI_Exscan(
              &meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal,
              &nCellsPreviousRanks3D_[meshName], 1, MPI_INT, MPI_SUM,
              this->rankSubset_->mpiCommunicator()),
          "MPI_Exscan");
      MPIUtility::handleReturnValue(
          MPI_Exscan(
              &meshPropertiesUnstructuredGridFile.at(meshName).nPointsLocal,
              &nPointsPreviousRanks3D_[meshName], 1, MPI_INT, MPI_SUM,
              this->rankSubset_->mpiCommunicator()),
          "MPI_Exscan");

      VLOG(1) << "meshName \"" << meshName << "\"";
      VLOG(1) << "nCellsLocal local: "
              << meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal
              << ", prefix sum: " << nCellsPreviousRanks3D_[meshName];
      VLOG(1) << "nPointsLocal local: "
              << meshPropertiesUnstructuredGridFile.at(meshName).nPointsLocal
              << ", prefix sum: " << nPointsPreviousRanks3D_[meshName];
    }

    Control::PerformanceMeasurement::stop("durationHDF53DReduction");
    Control::PerformanceMeasurement::stop("durationHDF53DInit");
  } else {
    VLOG(1) << "nCellsPreviousRanks3D_, nPointsPreviousRanks3D_ already set to "
            << nCellsPreviousRanks3D_ << ", " << nPointsPreviousRanks3D_;
  }

#ifndef NDEBUG
  const std::vector<node_no_t> &nNodesLocalWithGhosts =
      polyDataPropertiesForMesh.nNodesLocalWithGhosts;
  assert(nNodesLocalWithGhosts.size() == targetDimensionality);
#endif

  int nNodesPerCell = 4;
  if (output3DMeshes) {
    nNodesPerCell = 8;
  }

  // get local data values
  // setup connectivity array, which gives the node numbers for every
  // element/cell
  int nConnectivityValues = 0;

  // loop over the meshes that will be combined into the current file
  // meshPropertiesUnstructuredGridFile[...] contains information for all
  // meshes, polyDataPropertiesForMesh contains already combined information of
  // the set of meshes
  for (const std::string &meshName : meshNames) {
    nConnectivityValues +=
        meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal *
        nNodesPerCell;
    VLOG(1) << "  mesh \"" << meshName << "\" has "
            << meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal
            << " local cells, " << nNodesPerCell << ", "
            << meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal *
                   nNodesPerCell
            << " connectivity values for this mesh";
  }
  VLOG(1) << "nConnectivityValues: " << nConnectivityValues;

  std::vector<int> connectivityValues(nConnectivityValues);

  VLOG(1)
      << "n connectivity values from unstructured: "
      << polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.size();
  VLOG(1) << "nCellsLocal: " << polyDataPropertiesForMesh.nCellsLocal;

  // if connectivity values are already explicitly given, this is the case if we
  // have an unstructured mesh to output
  if (!polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.empty()) {
    if (meshNames.size() > 1) {
      LOG(FATAL) << "Cannot combine file for unstructured meshes.";
    }
    assert(
        polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.size() ==
        connectivityValues.size());

    VLOG(1) << "connectivityValues is initialized to "
            << connectivityValues.size() << ", values: " << connectivityValues;
    VLOG(1)
        << "now copy "
        << polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.size()
        << ", values: "
        << polyDataPropertiesForMesh.unstructuredMeshConnectivityValues;
    std::copy(
        polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.begin(),
        polyDataPropertiesForMesh.unstructuredMeshConnectivityValues.end(),
        connectivityValues.begin());
  } else {
    int connectivityValuesIndexOffset = 0;
    int connectivityValuesOffset = 0;

    for (const std::string &meshName : meshNames) {
      connectivityValuesOffset += nPointsPreviousRanks3D_[meshName];
    }

    VLOG(1) << "connectivityValuesOffset start: " << connectivityValuesOffset;

    // loop over meshes
    for (const std::string &meshName : meshNames) {
      int nPointsPreviousRanks3D = connectivityValuesOffset;

      const std::vector<node_no_t> &nNodesLocalWithGhosts =
          meshPropertiesUnstructuredGridFile.at(meshName).nNodesLocalWithGhosts;

      VLOG(1) << "connectivity for mesh \"" << meshName
              << "\",  nPointsPreviousRanks3D: " << nPointsPreviousRanks3D
              << ", nNodesLocalWithGhosts: " << nNodesLocalWithGhosts;
      VLOG(1) << "connectivityValuesIndexOffset: "
              << connectivityValuesIndexOffset;
      VLOG(1) << "nNodesLocalWithGhosts: " << nNodesLocalWithGhosts;

      // for structured meshes create connectivity values now
      if (output3DMeshes) {
        // fill connectivityValues for 3D meshes
        element_no_t elementIndex = 0;
        for (int indexZ = 0; indexZ < nNodesLocalWithGhosts[2] - 1; indexZ++) {
          for (int indexY = 0; indexY < nNodesLocalWithGhosts[1] - 1;
               indexY++) {
            for (int indexX = 0; indexX < nNodesLocalWithGhosts[0] - 1;
                 indexX++, elementIndex++) {
              if (connectivityValuesIndexOffset + elementIndex * 8 + 7 >=
                  connectivityValues.size()) {
                LOG(FATAL) << connectivityValuesIndexOffset + elementIndex * 8 +
                                  7
                           << ">= " << connectivityValues.size()
                           << ", connectivityValues are not large enough: "
                           << connectivityValues.size() << ", but "
                           << connectivityValuesIndexOffset << "+"
                           << nNodesLocalWithGhosts[0] - 1 << "x"
                           << nNodesLocalWithGhosts[1] - 1 << "x"
                           << nNodesLocalWithGhosts[2] - 1 << " = "
                           << connectivityValuesIndexOffset +
                                  (nNodesLocalWithGhosts[0] - 1) *
                                      (nNodesLocalWithGhosts[1] - 1) *
                                      (nNodesLocalWithGhosts[2] - 1)
                           << " elements";
              }

              connectivityValues[connectivityValuesIndexOffset +
                                 elementIndex * 8 + 0] =
                  nPointsPreviousRanks3D +
                  indexZ * nNodesLocalWithGhosts[0] * nNodesLocalWithGhosts[1] +
                  indexY * nNodesLocalWithGhosts[0] + indexX;
              connectivityValues[connectivityValuesIndexOffset +
                                 elementIndex * 8 + 1] =
                  nPointsPreviousRanks3D +
                  indexZ * nNodesLocalWithGhosts[0] * nNodesLocalWithGhosts[1] +
                  indexY * nNodesLocalWithGhosts[0] + indexX + 1;
              connectivityValues[connectivityValuesIndexOffset +
                                 elementIndex * 8 + 2] =
                  nPointsPreviousRanks3D +
                  indexZ * nNodesLocalWithGhosts[0] * nNodesLocalWithGhosts[1] +
                  (indexY + 1) * nNodesLocalWithGhosts[0] + indexX + 1;
              connectivityValues[connectivityValuesIndexOffset +
                                 elementIndex * 8 + 3] =
                  nPointsPreviousRanks3D +
                  indexZ * nNodesLocalWithGhosts[0] * nNodesLocalWithGhosts[1] +
                  (indexY + 1) * nNodesLocalWithGhosts[0] + indexX;
              connectivityValues[connectivityValuesIndexOffset +
                                 elementIndex * 8 + 4] =
                  nPointsPreviousRanks3D +
                  (indexZ + 1) * nNodesLocalWithGhosts[0] *
                      nNodesLocalWithGhosts[1] +
                  indexY * nNodesLocalWithGhosts[0] + indexX;
              connectivityValues[connectivityValuesIndexOffset +
                                 elementIndex * 8 + 5] =
                  nPointsPreviousRanks3D +
                  (indexZ + 1) * nNodesLocalWithGhosts[0] *
                      nNodesLocalWithGhosts[1] +
                  indexY * nNodesLocalWithGhosts[0] + indexX + 1;
              connectivityValues[connectivityValuesIndexOffset +
                                 elementIndex * 8 + 6] =
                  nPointsPreviousRanks3D +
                  (indexZ + 1) * nNodesLocalWithGhosts[0] *
                      nNodesLocalWithGhosts[1] +
                  (indexY + 1) * nNodesLocalWithGhosts[0] + indexX + 1;
              connectivityValues[connectivityValuesIndexOffset +
                                 elementIndex * 8 + 7] =
                  nPointsPreviousRanks3D +
                  (indexZ + 1) * nNodesLocalWithGhosts[0] *
                      nNodesLocalWithGhosts[1] +
                  (indexY + 1) * nNodesLocalWithGhosts[0] + indexX;
            }
          }
        }
      } else {
        // fill connectivityValues for 2D meshes
        element_no_t elementIndex = 0;
        for (int indexY = 0; indexY < nNodesLocalWithGhosts[1] - 1; indexY++) {
          for (int indexX = 0; indexX < nNodesLocalWithGhosts[0] - 1;
               indexX++, elementIndex++) {
            if (connectivityValuesIndexOffset + elementIndex * 4 + 3 >=
                connectivityValues.size()) {
              LOG(FATAL) << connectivityValuesIndexOffset + elementIndex * 4 + 3
                         << ">= " << connectivityValues.size()
                         << ", connectivityValues are not large enough: "
                         << connectivityValues.size() << ", but "
                         << connectivityValuesIndexOffset << " + "
                         << nNodesLocalWithGhosts[0] - 1 << "x"
                         << nNodesLocalWithGhosts[1] - 1 << " = "
                         << connectivityValuesIndexOffset +
                                (nNodesLocalWithGhosts[0] - 1) *
                                    (nNodesLocalWithGhosts[1] - 1)
                         << " elements";
            }

            connectivityValues[connectivityValuesIndexOffset +
                               elementIndex * 4 + 0] =
                nPointsPreviousRanks3D + indexY * nNodesLocalWithGhosts[0] +
                indexX;
            connectivityValues[connectivityValuesIndexOffset +
                               elementIndex * 4 + 1] =
                nPointsPreviousRanks3D + indexY * nNodesLocalWithGhosts[0] +
                indexX + 1;
            connectivityValues[connectivityValuesIndexOffset +
                               elementIndex * 4 + 2] =
                nPointsPreviousRanks3D +
                (indexY + 1) * nNodesLocalWithGhosts[0] + indexX + 1;
            connectivityValues[connectivityValuesIndexOffset +
                               elementIndex * 4 + 3] =
                nPointsPreviousRanks3D +
                (indexY + 1) * nNodesLocalWithGhosts[0] + indexX;
          }
        }
      }

      connectivityValuesIndexOffset +=
          meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal *
          nNodesPerCell;
      connectivityValuesOffset +=
          meshPropertiesUnstructuredGridFile.at(meshName).nPointsLocal;
    }
  }

  VLOG(1) << "nPointsPreviousRanks3D_: " << nPointsPreviousRanks3D_
          << ", nCellsPreviousRanks3D_: " << nCellsPreviousRanks3D_;
  VLOG(1) << "connectivity: " << connectivityValues;

  // setup offset array
  std::vector<int> offsetValues(polyDataPropertiesForMesh.nCellsLocal);

  size_t offsetValuesIndexOffset = 0;
  int offsetValuesOffset = 0;

  for (const std::string &meshName : meshNames) {
    offsetValuesOffset += nCellsPreviousRanks3D_[meshName] * nNodesPerCell;
  }

  // loop over meshes
  for (const std::string &meshName : meshNames) {
    for (size_t i = 0;
         i < meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal; i++) {
      assert(offsetValuesIndexOffset + i < offsetValues.size());
      // specifies the end, i.e. one after the last, of the last of nodes for
      // each element
      offsetValues[offsetValuesIndexOffset + i] =
          offsetValuesOffset + (i + 1) * nNodesPerCell;
    }

    offsetValuesIndexOffset +=
        meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal;
    offsetValuesOffset +=
        meshPropertiesUnstructuredGridFile.at(meshName).nCellsLocal *
        nNodesPerCell;
  }

  VLOG(1) << "offsetValues: " << offsetValues;

  // collect all data for the field variables, organized by field variable names
  std::map<std::string, std::vector<double>> fieldVariableValues;
  LoopOverTuple::loopGetNodalValues<FieldVariablesForOutputWriterType>(
      fieldVariables, meshNamesSet, fieldVariableValues);

  if (!meshPropertiesInitialized) {
    // if next assertion fails, output why for debugging
    if (fieldVariableValues.size() !=
        polyDataPropertiesForMesh.pointDataArrays.size()) {
      LOG(DEBUG) << "n field variable values: " << fieldVariableValues.size()
                 << ", n point data arrays: "
                 << polyDataPropertiesForMesh.pointDataArrays.size();
      LOG(DEBUG) << "mesh names: " << meshNames;
      std::stringstream pointDataArraysNames;
      for (const PolyDataPropertiesForMesh::DataArrayName &pointDataArray :
           polyDataPropertiesForMesh.pointDataArrays) {
        pointDataArraysNames << pointDataArray.name << " ";
      }
      LOG(DEBUG) << "pointDataArraysNames: " << pointDataArraysNames.str();
      LOG(DEBUG) << "FieldVariablesForOutputWriterType: "
                 << StringUtility::demangle(
                        typeid(FieldVariablesForOutputWriterType).name());
    }

    assert(fieldVariableValues.size() ==
           polyDataPropertiesForMesh.pointDataArrays.size());
  }

  // output 3D or 2D mesh
  if (!meshPropertiesInitialized) {
    // add field variable "partitioning" with 1 component
    PolyDataPropertiesForMesh::DataArrayName dataArrayName;
    dataArrayName.name = "partitioning";
    dataArrayName.nComponents = 1;
    dataArrayName.componentNames = std::vector<std::string>(1, "rankNo");

    polyDataPropertiesForMesh.pointDataArrays.push_back(dataArrayName);
  }

  // set data for partitioning field variable
  assert(!fieldVariableValues.empty());
  fieldVariableValues["partitioning"].resize(
      polyDataPropertiesForMesh.nPointsLocal,
      (double)this->rankSubset_->ownRankNo());

  // for 2D field variables, add 0 every 2nd entry to generate 3D values,
  // because paraview cannot handle 2D vectors properly
  if (targetDimensionality == 2) {
    std::vector<double> buffer;

    // loop over all field variables
    for (const PolyDataPropertiesForMesh::DataArrayName &pointDataArray :
         polyDataPropertiesForMesh.pointDataArrays) {
      // if it is a 2D vector field
      if (pointDataArray.nComponents == 2) {
        std::vector<double> &values = fieldVariableValues[pointDataArray.name];

        // copy all values to a buffer
        buffer.assign(values.begin(), values.end());
        int nValues = buffer.size();

        // resize the values vector to 2/3 the size
        values.resize(nValues / 2 * 3);

        // loop over the new values vector and set the entries
        for (int i = 0; i < nValues / 2; i++) {
          values[3 * i + 0] = buffer[2 * i + 0];
          values[3 * i + 1] = buffer[2 * i + 1];
          values[3 * i + 2] = 0.0;
        }
      }
    }
  }

  // collect all data for the geometry field variable
  std::vector<double> geometryFieldValues;
  LoopOverTuple::loopGetGeometryFieldNodalValues<
      FieldVariablesForOutputWriterType>(fieldVariables, meshNamesSet,
                                         geometryFieldValues);

  VLOG(1) << "meshNames: " << meshNames << ", rank "
          << this->rankSubset_->ownRankNo()
          << ", n geometryFieldValues: " << geometryFieldValues.size();
  if (geometryFieldValues.size() == 0) {
    LOG(FATAL)
        << "There is no geometry field. You have to provide a geomteryField in "
           "the field variables returned by getFieldVariablesForOutputWriter!";
  }

  Control::PerformanceMeasurement::start("durationHDF53DWrite");

  herr_t err;
  // write field variables
  for (PolyDataPropertiesForMesh::DataArrayName &pointDataArray :
       polyDataPropertiesForMesh.pointDataArrays) {
    assert(fieldVariableValues.find(pointDataArray.name) !=
           fieldVariableValues.end());

    // write values
    // for partitioning, convert float values to integer values for output
    bool writeFloatsAsInt = pointDataArray.name == "partitioning";
    err = writeCombinedValuesVector(
        fileID, fieldVariableValues[pointDataArray.name],
        pointDataArray.name.c_str(), writeFloatsAsInt);
    assert(err >= 0);
  }

  err = writeCombinedValuesVector(fileID, geometryFieldValues, "geometry");
  assert(err >= 0);
  err = writeCombinedValuesVector(fileID, connectivityValues, "connectivity");
  assert(err >= 0);
  err = writeCombinedValuesVector(fileID, offsetValues, "offsets");
  assert(err >= 0);
  err = writeCombinedTypesVector(fileID, polyDataPropertiesForMesh.nCellsGlobal,
                                 output3DMeshes, "types");
  assert(err >= 0);

  Control::PerformanceMeasurement::stop("durationHDF53DWrite");
}
} // namespace OutputWriter
