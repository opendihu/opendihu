#include "control/precice/nested_solver.h"

namespace Control {
template <typename T1>
std::shared_ptr<typename PreciceAdapterNestedSolver<
    FastMonodomainSolver<T1>>::FunctionSpace>
PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::functionSpace(
    NestedSolverType &nestedSolver) {
  return nestedSolver.data().functionSpace();
}

template <typename T1>
template <typename PreciceVolumeData, typename PreciceVolumeMesh>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::
    initializePreciceVolumeData(
        PythonConfig &specificSettings, NestedSolverType &nestedSolver,
        std::shared_ptr<precice::Participant> &preciceParticipant,
        std::vector<PreciceVolumeData> &preciceVolumeData,
        std::vector<std::shared_ptr<PreciceVolumeMesh>> &preciceVolumeMeshes) {
  // parse settings for coupling participants
  // loop over items of the key "preciceData"
  std::string settingsKey("preciceVolumeData");
  PyObject *listPy = specificSettings.getOptionPyObject(settingsKey);
  std::vector<PyObject *> list =
      PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);
  PythonConfig preciceDataConfig(specificSettings, settingsKey);

  // get the slot connector data from the nested solver
  int slotNo = 0;
  using SlotConnectorDataType =
      typename NestedSolverType::SlotConnectorDataType;
  std::shared_ptr<SlotConnectorDataType> slotConnectorData =
      nestedSolver.getSlotConnectorData();

  // get all present slot names
  std::vector<std::string> slotNames;
  SlotConnectorDataHelper<SlotConnectorDataType>::getSlotNames(
      slotConnectorData, slotNames);

  // loop over items of the list under "preciceData", the order of the items
  // is also the order of the data connector slots
  for (int listIndex = 0; listIndex < list.size(); listIndex++) {
    // extract the list item as PythonConfig {"mode": "write", "variableName":
    // ...}
    PythonConfig currentPreciceData(preciceDataConfig, listIndex);

    // create a new preciceData instance that will be added to
    // preciceVolumeData_ at the end of this loop
    PreciceVolumeData preciceData;

    // parse options
    preciceData.isGeometryField =
        currentPreciceData.getOptionBool("isGeometryField", false);
    preciceData.opendihuMeshName =
        currentPreciceData.getOptionString("opendihuMeshName", "");

    // parse slot name
    preciceData.slotName = currentPreciceData.getOptionString("slotName", "");
    preciceData.slotNo = slotNo;

    LOG(DEBUG) << "isGeometryField: " << preciceData.isGeometryField
               << ", slotName: \"" << preciceData.slotName << "\".";

    if (preciceData.slotName == "") {
      preciceData.slotNo = slotNo;
      LOG(DEBUG) << "Using slotNo " << slotNo
                 << " because slotName was not given.";
    } else {
      LOG(DEBUG) << "The slotName \"" << preciceData.slotName << "\" is given.";
      std::vector<std::string>::iterator slotNameIter =
          std::find(slotNames.begin(), slotNames.end(), preciceData.slotName);

      // if the "from" and "to" slot names match to slot names in the current
      // slotConnectorData
      if (slotNameIter != slotNames.end()) {
        LOG(DEBUG) << "Slot found";
        // determine slot nos
        preciceData.slotNo = std::distance(slotNames.begin(), slotNameIter);
        LOG(DEBUG) << "Slot \"" << preciceData.slotName << "\" is slot no "
                   << preciceData.slotNo;
      } else {
        LOG(DEBUG) << "Slot not found";
        preciceData.slotNo = slotNo;
        LOG(WARNING)
            << "A slot with name \"" << preciceData.slotName
            << "\" does not exist. Available slots: " << slotNames
            << ". Use `None` or the empty string to specify the slot by its "
               "number (referenced slots will be in the order of "
            << "the items under \"preciceData\").\n Using slot No " << slotNo;
      }
    }

    // parse the mesh
    std::string currentMeshName =
        currentPreciceData.getOptionString("preciceMeshName", "");

    // find the mesh in the already parsed precice meshes
    typename std::vector<std::shared_ptr<PreciceVolumeMesh>>::iterator iter =
        std::find_if(
            preciceVolumeMeshes.begin(), preciceVolumeMeshes.end(),
            [&currentMeshName](std::shared_ptr<PreciceVolumeMesh> preciceMesh) {
              return preciceMesh->preciceMeshName == currentMeshName;
            });

    // if the mesh is not in preciceVolumeMeshes, create it and add it to
    // preciceVolumeMeshes
    if (iter == preciceVolumeMeshes.end()) {
      //   LOG(DEBUG) << "mesh was not yet initialized.";

      // create new precice mesh object
      std::shared_ptr<PreciceVolumeMesh> preciceMesh =
          std::make_shared<PreciceVolumeMesh>();
      preciceMesh->preciceMeshName = currentMeshName;

      std::vector<Vec3> geometryValues;
      int nArrayItems = 1;

      // get the function space that is used to initialize the mapping
      std::shared_ptr<typename NestedSolverType::FunctionSpace> functionSpace =
          nullptr;
      if (preciceData.opendihuMeshName != "") {
        functionSpace =
            DihuContext::meshManager()
                ->functionSpace<typename NestedSolverType::FunctionSpace>(
                    preciceData.opendihuMeshName);
        LOG(DEBUG) << "Using opendihu mesh with name \""
                   << preciceData.opendihuMeshName
                   << "\" to initialize precice mapping.";

        // get the node positions of the opendihu mesh to initialize the
        // mapping
        preciceMesh->nNodesLocal =
            functionSpace->geometryField().nDofsLocalWithoutGhosts();
        functionSpace->geometryField().getValuesWithoutGhosts(geometryValues);

        // store opendihu mesh name
        preciceMesh->opendihuMeshName = functionSpace->meshName();

        // the following cannot happen
        if (preciceMesh->opendihuMeshName != preciceData.opendihuMeshName)
          LOG(FATAL)
              << "Initializing precice mesh from settings opendihuMeshName=\""
              << preciceData.opendihuMeshName << "\" and the resulting mesh is "
              << preciceMesh->opendihuMeshName;
      } else {
        LOG(DEBUG) << "get mesh from mesh partition, slot No "
                   << preciceData.slotNo;
        LOG(DEBUG) << "slot connector data type: "
                   << StringUtility::demangle(
                          typeid(SlotConnectorDataType).name());

        // get the mesh partition
        std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
            SlotConnectorDataHelper<
                SlotConnectorDataType>::getMeshPartitionBase(slotConnectorData,
                                                             preciceData.slotNo,
                                                             0);

        if (!meshPartitionBase) {
          LOG(FATAL) << "Could not get mesh for slot No " << preciceData.slotNo;
        } else
          LOG(DEBUG) << "got mesh partition for slot No " << preciceData.slotNo;

        // get opendihu mesh name
        preciceMesh->opendihuMeshName =
            SlotConnectorDataHelper<SlotConnectorDataType>::getMeshName(
                slotConnectorData, preciceData.slotNo);

        int nDofsLocalWithoutGhosts =
            meshPartitionBase->nDofsLocalWithoutGhosts();
        nArrayItems =
            SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
                slotConnectorData,
                preciceData.slotNo); // number of fibers if there are fibers

        preciceMesh->nNodesLocal = nDofsLocalWithoutGhosts * nArrayItems;

        // get the vector of values [0,1,...,nDofsLocalWithGhosts]
        const std::vector<PetscInt> &dofNosLocalWithGhosts =
            meshPartitionBase->dofNosLocal();
        std::vector<PetscInt> dofNosLocalWithoutGhosts(
            dofNosLocalWithGhosts.begin(),
            dofNosLocalWithGhosts.begin() + nDofsLocalWithoutGhosts);

        // loop over fibers if there are any
        for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
          static std::vector<Vec3> nodePositionsFiber;
          nodePositionsFiber.clear();
          SlotConnectorDataHelper<SlotConnectorDataType>::slotGetGeometryValues(
              slotConnectorData, preciceData.slotNo, arrayIndex,
              dofNosLocalWithoutGhosts, nodePositionsFiber);
          geometryValues.insert(geometryValues.end(),
                                nodePositionsFiber.begin(),
                                nodePositionsFiber.end());
        }

        LOG(DEBUG) << "collected " << geometryValues.size()
                   << " node positions from the " << nArrayItems << " fibers";
      }

      // transform to contiguous memory layout for precice
      std::vector<double> geometryValuesContiguous(3 * geometryValues.size());
      for (int entryNo = 0; entryNo < geometryValues.size(); entryNo++)
        for (int componentNo = 0; componentNo < 3; componentNo++)
          geometryValuesContiguous[3 * entryNo + componentNo] =
              geometryValues[entryNo][componentNo];

      if (geometryValuesContiguous.size() != preciceMesh->nNodesLocal * 3)
        LOG(FATAL) << "size mismatch: " << geometryValuesContiguous.size()
                   << "!=" << preciceMesh->nNodesLocal * 3;

      // resize buffer for vertex ids
      preciceMesh->preciceVertexIds.resize(preciceMesh->nNodesLocal);

      // give the node positions to precice and get the vertex ids
      preciceParticipant->setMeshVertices(preciceMesh->preciceMeshName,
                                          geometryValuesContiguous,
                                          preciceMesh->preciceVertexIds);

      // store newly created mesh in preciceData
      preciceData.preciceMesh = preciceMesh;

      // store the precice mesh to the vector of meshes
      preciceVolumeMeshes.push_back(preciceMesh);

      // output message about mesh names
      std::stringstream message;
      message << "Initialized precice mesh \"" << preciceMesh->preciceMeshName
              << "\" from opendihu mesh \"" << preciceMesh->opendihuMeshName
              << "\" with " << preciceMesh->nNodesLocal << " local nodes";
      if (nArrayItems > 1)
        message << " (and other corresponding meshes, in total " << nArrayItems
                << " fibers/compartments)";
      LOG(INFO) << message.str() << ".";
    } else {
      LOG(DEBUG) << "Use existing precice mesh " << currentMeshName;
      preciceData.preciceMesh = *iter;
    }

    // parse mode
    std::string mode = currentPreciceData.getOptionString("mode", "");
    if (mode == "read") {
      preciceData.ioType =
          PreciceAdapterInitialize<NestedSolverType>::PreciceVolumeData::ioRead;
    } else if (mode == "write") {
      preciceData.ioType = PreciceAdapterInitialize<
          NestedSolverType>::PreciceVolumeData::ioWrite;
    } else {
      LOG(FATAL) << currentPreciceData << "[\"mode\"] is \"" << mode << "\", "
                 << "possible values are: \"read\", \"write\".";
    }

    // parse variable name
    preciceData.preciceDataName =
        currentPreciceData.getOptionString("preciceDataName", "variable");

    // increment slotNo, this value is used for slots that are not identified
    // by slotName in the config
    if (!preciceData.isGeometryField)
      slotNo++;

    // check if mesh size matches the specified slot
    // get the mesh partition
    std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
        SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase(
            slotConnectorData, preciceData.slotNo, 0);

    int nDofsLocalWithoutGhosts = meshPartitionBase->nDofsLocalWithoutGhosts();
    int nArrayItems =
        SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
            slotConnectorData,
            preciceData.slotNo); // number of fibers if there are fibers

    if (nDofsLocalWithoutGhosts * nArrayItems !=
        preciceData.preciceMesh->nNodesLocal) {
      LOG(DEBUG) << ", all available slots: "
                 << SlotConnectorDataHelper<SlotConnectorDataType>::getString(
                        slotConnectorData);
      LOG(FATAL)
          << currentPreciceData
          << ": Mesh does not match slot in "
             "PreciceAdapterVolumeCoupling.\n\n "
          << "The " << (listIndex + 1)
          << (listIndex == 0
                  ? "st"
                  : (listIndex == 1 ? "nd" : (listIndex == 2 ? "rd" : "th")))
          << " list item under \"preciceData\" uses slot " << preciceData.slotNo
          << " (slotName \"" << preciceData.slotName << "\")"
          << " and preciceMesh \"" << preciceData.preciceMesh->preciceMeshName
          << "\".\nThe slot has " << nDofsLocalWithoutGhosts << " dofs * "
          << nArrayItems
          << " array items (fibers) = " << nDofsLocalWithoutGhosts * nArrayItems
          << " dofs in total. The mesh was initialized from opendihu mesh \""
          << preciceData.preciceMesh->opendihuMeshName << "\" with "
          << preciceData.preciceMesh->nNodesLocal << " nodes. \n"
          << "(isGeometryField: " << std::boolalpha
          << preciceData.isGeometryField << ")\n"
          << "You can do the following:\n"
          << "- Rename the preciceMesh to a unique name and update the "
             "precice "
             "config xml file accordingly.\n"
          << "- Check the solver structure file to find out the meshes that "
             "are associated with the slots.\n"
          << "  Check that the correct mesh is specified and that the "
             "precice "
             "mesh was not initialized earlier with a different opendihu "
             "mesh.\n"
          << "  Every precice mesh is initialized the first time in appears "
             "under \"preciceData\", sometimes reordering the entries can "
             "help.\n"
          << "  If the desired mesh is not available at the current solver, "
             "maybe insert a MapDofs class.";
    }

    LOG(INFO) << "Precice data \"" << preciceData.preciceDataName
              << "\" maps to "
              << (preciceData.isGeometryField ? "the geometry field of " : "")
              << "slot " << preciceData.slotNo << " (\"" << preciceData.slotName
              << "\") and uses precice mesh \""
              << preciceData.preciceMesh->preciceMeshName
              << "\", which is opendihu mesh \""
              << preciceData.preciceMesh->opendihuMeshName << "\" with "
              << preciceData.preciceMesh->nNodesLocal << " local nodes.";

    // store preciceData to vector
    preciceVolumeData.push_back(preciceData);
  }
}

template <typename T1>
template <typename SurfaceDataVector, typename VolumeDataVector>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::preciceReadData(
    NestedSolverType &nestedSolver,
    std::shared_ptr<precice::Participant> &preciceParticipant,
    SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData,
    DihuContext &context) {
  LOG(DEBUG) << "read volume data from precice";
  double preciceDt = preciceParticipant->getMaxTimeStepSize();

  using SlotConnectorDataType =
      typename NestedSolverType::SlotConnectorDataType;
  std::shared_ptr<SlotConnectorDataType> slotConnectorData =
      nestedSolver.getSlotConnectorData();

  // loop over data
  for (auto &preciceData : preciceVolumeData) {
    if (preciceData.ioType ==
        PreciceAdapterInitialize<NestedSolverType>::PreciceVolumeData::ioRead) {
      int nEntries = preciceData.preciceMesh->nNodesLocal;

      if (preciceData.isGeometryField) {
        nEntries = preciceData.preciceMesh->nNodesLocal * 3;
      }

      // allocate temporary memory
      scalarValues_.resize(nEntries);

      // get all data at once
      preciceParticipant->readData(
          preciceData.preciceMesh->preciceMeshName, preciceData.preciceDataName,
          preciceData.preciceMesh->preciceVertexIds, preciceDt, scalarValues_);

      // get the mesh partition
      std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
          SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase(
              slotConnectorData, preciceData.slotNo, 0);

      int nDofsLocalWithoutGhosts =
          meshPartitionBase->nDofsLocalWithoutGhosts();

      // get the vector of values [0,1,...,nDofsLocalWithGhosts]
      const std::vector<PetscInt> &dofNosLocalWithGhosts =
          meshPartitionBase->dofNosLocal();
      std::vector<PetscInt> dofNosLocalWithoutGhosts(
          dofNosLocalWithGhosts.begin(),
          dofNosLocalWithGhosts.begin() + nDofsLocalWithoutGhosts);

      int nArrayItems =
          SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
              slotConnectorData,
              preciceData.slotNo); // number of fibers if there are fibers

      // store received data in field variable
      if (preciceData.isGeometryField) {
        // loop over fibers if there are any
        for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
          // fill the vector geometryValues_ with the geometry values of the
          // current fiber or mesh
          geometryValues_.resize(nDofsLocalWithoutGhosts);
          for (int dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts;
               dofNoLocal++) {
            for (int componentNo = 0; componentNo < 3; componentNo++) {
              geometryValues_[dofNoLocal][componentNo] =
                  scalarValues_[3 * (arrayIndex * nDofsLocalWithoutGhosts +
                                     dofNoLocal) +
                                componentNo];
              //
            }
          }

          SlotConnectorDataHelper<SlotConnectorDataType>::slotSetGeometryValues(
              slotConnectorData, preciceData.slotNo, arrayIndex,
              dofNosLocalWithoutGhosts, geometryValues_);
        }
      } else {
        // loop over fibers if there are any
        for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
          // fill the vector geometryValues_ with the geometry values of the
          // current fiber or mesh
          scalarValuesOfMesh_.assign(
              scalarValues_.begin() + arrayIndex * nDofsLocalWithoutGhosts,
              scalarValues_.begin() +
                  (arrayIndex + 1) * nDofsLocalWithoutGhosts);

          SlotConnectorDataHelper<SlotConnectorDataType>::slotSetValues(
              slotConnectorData, preciceData.slotNo, arrayIndex,
              dofNosLocalWithoutGhosts, scalarValuesOfMesh_);
        }
      }
    }
  }
}

template <typename T1>
template <typename SurfaceDataVector, typename VolumeDataVector>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::preciceWriteData(
    NestedSolverType &nestedSolver,
    std::shared_ptr<precice::Participant> &preciceParticipant,
    SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData,
    double scalingFactor) {
  // write data to precice
  LOG(DEBUG) << "write volume data to precice";

  using SlotConnectorDataType =
      typename NestedSolverType::SlotConnectorDataType;
  std::shared_ptr<SlotConnectorDataType> slotConnectorData =
      nestedSolver.getSlotConnectorData();

  // loop over data
  for (auto &preciceData : preciceVolumeData) {
    if (preciceData.ioType ==
        PreciceAdapterInitialize<
            NestedSolverType>::PreciceVolumeData::ioWrite) {
      scalarValues_.clear();

      // get the mesh partition
      std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
          SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase(
              slotConnectorData, preciceData.slotNo, 0);

      int nDofsLocalWithoutGhosts =
          meshPartitionBase->nDofsLocalWithoutGhosts();

      // get the vector of values [0,1,...,nDofsLocalWithGhosts]
      const std::vector<PetscInt> &dofNosLocalWithGhosts =
          meshPartitionBase->dofNosLocal();
      std::vector<PetscInt> dofNosLocalWithoutGhosts(
          dofNosLocalWithGhosts.begin(),
          dofNosLocalWithGhosts.begin() + nDofsLocalWithoutGhosts);

      int nArrayItems =
          SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
              slotConnectorData,
              preciceData.slotNo); // number of fibers if there are fibers

      // if it is a geometry field, get the node positions of a mesh
      if (preciceData.isGeometryField) {
        geometryValues_.clear();

        // loop over fibers if there are any
        for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
          static std::vector<Vec3> geometryValuesFiber;
          geometryValuesFiber.clear();
          SlotConnectorDataHelper<SlotConnectorDataType>::slotGetGeometryValues(
              slotConnectorData, preciceData.slotNo, arrayIndex,
              dofNosLocalWithoutGhosts, geometryValuesFiber);

          geometryValues_.insert(geometryValues_.end(),
                                 geometryValuesFiber.begin(),
                                 geometryValuesFiber.end());
        }
        LOG(DEBUG) << "Using geometry field of opendihu meshes of "
                   << nArrayItems << " array items (e.g., fibers).";

        // transform to contiguous memory layout for precice
        scalarValues_.resize(3 * geometryValues_.size());

        for (int entryNo = 0; entryNo < geometryValues_.size(); entryNo++)
          for (int componentNo = 0; componentNo < 3; componentNo++)
            scalarValues_[3 * entryNo + componentNo] =
                geometryValues_[entryNo][componentNo];

      } else {
        // the data is a normal slot, no geometry field

        // loop over fibers if there are any
        for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
          static std::vector<double> values;
          values.clear();
          // static void slotGetValues(std::shared_ptr<SlotConnectorDataType>
          // slotConnectorData,
          //   int slotNo, int arrayIndex, const std::vector<dof_no_t>
          //   &dofNosLocal, std::vector<double> &values);
          SlotConnectorDataHelper<SlotConnectorDataType>::slotGetValues(
              slotConnectorData, preciceData.slotNo, arrayIndex,
              dofNosLocalWithoutGhosts, values);
          scalarValues_.insert(scalarValues_.end(), values.begin(),
                               values.end());

          LOG(DEBUG) << "arrayIndex " << arrayIndex << ", add " << values.size()
                     << " values, now number: " << scalarValues_.size();
        }
      }

      // scale the values by a factor given in the config
      for (double &value : scalarValues_)
        value *= scalingFactor;

      // check dim for scalarValues_
      if (preciceData.isGeometryField) {
        if (scalarValues_.size() != 3 * preciceData.preciceMesh->nNodesLocal)
          LOG(FATAL) << "Wrong number of scalar values (isGeometryField): "
                     << scalarValues_.size()
                     << " != " << preciceData.preciceMesh->nNodesLocal
                     << ", nArrayItems: " << nArrayItems;
      } else {
        if (scalarValues_.size() != preciceData.preciceMesh->nNodesLocal)
          LOG(FATAL) << "Wrong number of scalar values: "
                     << scalarValues_.size()
                     << " != " << preciceData.preciceMesh->nNodesLocal
                     << ", nArrayItems: " << nArrayItems;
      }

      // write values in precice
      preciceParticipant->writeData(
          preciceData.preciceMesh->preciceMeshName, preciceData.preciceDataName,
          preciceData.preciceMesh->preciceVertexIds, scalarValues_);
    }
  }
  LOG(DEBUG) << "write volume data to precice complete";
}

template <typename T1>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::
    addDirichletBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::vector<
            typename SpatialDiscretization::DirichletBoundaryConditionsBase<
                FunctionSpace, 6>::ElementWithNodes>
            &dirichletBoundaryConditionElements) {}

//! update existing boundary conditions with new values
template <typename T1>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::
    updateDirichletBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::vector<std::pair<global_no_t, std::array<double, 6>>>
            newDirichletBoundaryConditionValues) {}

template <typename T1>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::
    updateNeumannBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>>
            neumannBoundaryConditions) {}

//! get the displacement and velocity vectors of the given local dof nos
template <typename T1>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::
    getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                  const std::vector<dof_no_t> &dofNosLocal,
                                  std::vector<double> &displacementValues,
                                  std::vector<double> &velocityValues) {}

//! get the traction vectors of the given local dof nos
template <typename T1>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::getTractionValues(
    NestedSolverType &nestedSolver, const std::vector<dof_no_t> &dofNosLocal,
    std::vector<double> &tractionValues) {}

template <typename T1>
Vec PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::currentState(
    NestedSolverType &nestedSolver) {}

template <typename T1>
std::shared_ptr<
    FieldVariable::FieldVariable<typename PreciceAdapterNestedSolver<
                                     FastMonodomainSolver<T1>>::FunctionSpace,
                                 9>>
PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::deformationGradientField(
    NestedSolverType &nestedSolver) {}

template <typename T1>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::reset(
    NestedSolverType &nestedSolver) {
  // nestedSolver.timeStepping1().reset();
  nestedSolver.reset();
}

//! save fibers checkpoint
template <typename T1>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::saveFiberData(
    NestedSolverType &nestedSolver) {}

//! load fibers checkpoint
template <typename T1>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::loadFiberData(
    NestedSolverType &nestedSolver) {}

} // namespace Control
