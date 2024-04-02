#include "control/precice/initialize.h"

#include <sstream>

namespace Control {

template <typename NestedSolver>
PreciceAdapterInitialize<NestedSolver>::PreciceAdapterInitialize(
    DihuContext context)
    : context_(context["PreciceAdapter"]), nestedSolver_(this->context_),
      maximumPreciceTimestepSize_(0), timeStepOutputInterval_(1),
      initialized_(false) {
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();
}

template <typename NestedSolver>
void PreciceAdapterInitialize<NestedSolver>::initialize() {
#ifdef HAVE_PRECICE

  LOG(DEBUG) << "initialize precice adapter, initialized_=" << initialized_;

  // make sure that we initialize only once, in the next call, initialized_ is
  // true
  if (initialized_)
    return;

  // initialize() will be called before the simulation starts.

  // add this solver to the solvers diagram, which is an ASCII art
  // representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver(
      "Control::PreciceAdapter",
      true); // hasInternalConnectionToFirstNestedSolver=true (the last
             // argument) means slot connector data is shared with the first
             // subsolver

  // indicate in solverStructureVisualizer that now a child solver will be
  // initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested solver
  nestedSolver_.initialize();

  // initialize function space
  functionSpace_ = this->functionSpace(nestedSolver_);

  // indicate in solverStructureVisualizer that the child solver initialization
  // is done
  DihuContext::solverStructureVisualizer()->endChild();

  // set the slotConnectorData for the solverStructureVisualizer to appear in
  // the solver diagram
  // DihuContext::solverStructureVisualizer()->setSlotConnectorData(this->getSlotConnectorData());

  // check if the coupling is enabled
  couplingEnabled_ =
      this->specificSettings_.getOptionBool("couplingEnabled", true);

  timeStepWidth_ = this->specificSettings_.getOptionDouble(
      "timestepWidth", 0.01, PythonUtility::Positive);

  // if not enabled, abort initialization
  if (!couplingEnabled_) {
    endTimeIfCouplingDisabled_ = this->specificSettings_.getOptionDouble(
        "endTimeIfCouplingDisabled", 1, PythonUtility::Positive);

    LOG(WARNING) << "Coupling in PreciceAdapterVolumeCoupling is disabled "
                    "(option \"couplingEnabled\": False), "
                 << "using end time \"endTimeIfCouplingDisabled\": "
                 << endTimeIfCouplingDisabled_ << ".";

    initialized_ = true;
    return;
  }

  // initialize precice
  preciceParticipantName_ = this->specificSettings_.getOptionString(
      "preciceParticipantName", "MuscleSolver");
  const std::string configFileName = this->specificSettings_.getOptionString(
      "preciceConfigFilename", "../precice-config.xml");
  outputOnlyConvergedTimeSteps_ = this->specificSettings_.getOptionBool(
      "outputOnlyConvergedTimeSteps", true);

  int rankNo = functionSpace_->meshPartition()->rankSubset()->ownRankNo();
  int nRanks = functionSpace_->meshPartition()->rankSubset()->size();

  // initialize interface to precice for the bottom surface mesh
  preciceParticipant_ = std::make_shared<precice::Participant>(
      preciceParticipantName_, configFileName, rankNo, nRanks);

  // parse the options in "preciceSurfaceMeshes" and initialize all meshes in
  // precice, store in variable preciceSurfaceMeshes_
  initializePreciceSurfaceMeshes();

  // parse the options in "preciceSurfaceData" and initialize all variables in
  // precice, store in variable preciceSurfaceData_
  initializePreciceSurfaceData();

  // initialize Dirichlet boundary conditions at all dofs that will get some
  // prescribed values during coupling
  initializeDirichletBoundaryConditions();

  this->initializePreciceVolumeData(this->specificSettings_, nestedSolver_,
                                    preciceParticipant_, preciceVolumeData_,
                                    preciceVolumeMeshes_);
  // initializePreciceVolumeData();
  //  parse scalingFactor from settings
  scalingFactor_ = this->specificSettings_.getOptionDouble("scalingFactor", 1);

  // preciceParticipant_->initialize();
  // // determine maximum timestep size
  // maximumPreciceTimestepSize_ = preciceParticipant_->getMaxTimeStepSize();

  // LOG(DEBUG) << "precice initialization done, dt: "
  //            << maximumPreciceTimestepSize_ << "," << timeStepWidth_;

  // initialized_ = true;

#else
  LOG(FATAL) << "Failed to initialize PreciceAdapter (surface coupling) "
                "because opendihu is not compiled with preCICE.";
#endif
}

#ifdef HAVE_PRECICE
template <typename NestedSolver>
void PreciceAdapterInitialize<NestedSolver>::initializePreciceSurfaceMeshes() {
  const int nNodesX = functionSpace_->nNodesLocalWithoutGhosts(0);
  const int nNodesY = functionSpace_->nNodesLocalWithoutGhosts(1);
  const int nNodesZ = functionSpace_->nNodesLocalWithoutGhosts(2);

  std::vector<Vec3> geometryValues;
  functionSpace_->geometryField().getValuesWithoutGhosts(geometryValues);

  // parse settings of meshes
  std::string settingsKey("preciceSurfaceMeshes");
  PyObject *listPy = this->specificSettings_.getOptionPyObject(settingsKey);
  std::vector<PyObject *> list =
      PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);
  PythonConfig preciceMeshConfig(this->specificSettings_, settingsKey);

  // loop over items of the list under "preciceMeshes"
  for (int i = 0; i < list.size(); i++) {
    PythonConfig currentMeshConfig(preciceMeshConfig, i);

    std::shared_ptr<PreciceSurfaceMesh> preciceMesh =
        std::make_shared<PreciceSurfaceMesh>();

    // parse name of mesh
    preciceMesh->preciceMeshName =
        currentMeshConfig.getOptionString("preciceMeshName", "");

    // parse face
    std::string face = currentMeshConfig.getOptionString("face", "2-");
    if (face == "2-") {
      preciceMesh->face = PreciceSurfaceMesh::face2Minus;
    } else if (face == "2+") {
      preciceMesh->face = PreciceSurfaceMesh::face2Plus;
    } else {
      LOG(FATAL) << currentMeshConfig << "[\"face\"] is \"" << face
                 << "\", valid values are: \"2-\", \"2+\".";
    }

    // check if there are any local nodes of the surface on the local partition
    bool localDomainHasPartOfSurface = true;
    if (preciceMesh->face == PreciceSurfaceMesh::face2Minus &&
        functionSpace_->meshPartition()->ownRankPartitioningIndex(2) > 0) {
      localDomainHasPartOfSurface = false;
    } else if (preciceMesh->face == PreciceSurfaceMesh::face2Plus &&
               functionSpace_->meshPartition()->ownRankPartitioningIndex(2) <
                   functionSpace_->meshPartition()->nRanks(2) - 1) {
      localDomainHasPartOfSurface = false;
    }

    if (localDomainHasPartOfSurface) {
      // store number of nodes
      preciceMesh->nNodesLocal = nNodesX * nNodesY;

      // collect node positions for all surface nodes of the coupling surface
      std::vector<double> geometryValuesSurface(3 * preciceMesh->nNodesLocal);

      // resize buffer for the local dof nos in the 3D mesh of the surface mesh
      preciceMesh->dofNosLocal.resize(preciceMesh->nNodesLocal);

      int nodeIndexZ = 0;
      if (face == "2+") {
        nodeIndexZ = (nNodesZ - 1);
      }

      // loop over nodes
      for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++) {
        for (int nodeIndexX = 0; nodeIndexX < nNodesX; nodeIndexX++) {
          node_no_t surfaceDofNo = nodeIndexY * nNodesX + nodeIndexX;
          dof_no_t dofNoLocal = nodeIndexZ * nNodesX * nNodesY +
                                nodeIndexY * nNodesX + nodeIndexX;

          preciceMesh->dofNosLocal[surfaceDofNo] = dofNoLocal;

          for (int i = 0; i < 3; i++) {
            geometryValuesSurface[3 * surfaceDofNo + i] =
                geometryValues[dofNoLocal][i];
          }
        }
      }

      preciceMesh->preciceVertexIds.resize(preciceMesh->nNodesLocal);
      // give the node positions to precice and get the vertex ids
      preciceParticipant_->setMeshVertices(preciceMesh->preciceMeshName,
                                           geometryValuesSurface,
                                           preciceMesh->preciceVertexIds);
    } else {
      // there are no local nodes of this surface
      preciceMesh->nNodesLocal = 0;
    }

    // store the precice mesh to the vector of meshes
    preciceSurfaceMeshes_.push_back(preciceMesh);
  }
}

template <typename NestedSolver>
void PreciceAdapterInitialize<NestedSolver>::initializePreciceSurfaceData() {
  // parse settings for coupling participants / tendons
  // loop over items of the key "preciceData"
  std::string settingsKey("preciceSurfaceData");
  PyObject *listPy = this->specificSettings_.getOptionPyObject(settingsKey);
  std::vector<PyObject *> list =
      PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);
  PythonConfig preciceDataConfig(this->specificSettings_, settingsKey);

  // loop over items of the list under "preciceData"
  for (int i = 0; i < list.size(); i++) {
    PythonConfig currentPreciceData(preciceDataConfig, i);

    PreciceSurfaceData preciceData;

    // parse the mesh
    std::string currentMeshName =
        currentPreciceData.getOptionString("preciceMeshName", "");

    // find the mesh in the already parsed precice meshes
    typename std::vector<std::shared_ptr<PreciceSurfaceMesh>>::iterator iter =
        std::find_if(preciceSurfaceMeshes_.begin(), preciceSurfaceMeshes_.end(),
                     [&currentMeshName](
                         std::shared_ptr<PreciceSurfaceMesh> preciceMesh) {
                       return preciceMesh->preciceMeshName == currentMeshName;
                     });

    if (iter == preciceSurfaceMeshes_.end()) {
      std::stringstream s;
      for (std::shared_ptr<PreciceSurfaceMesh> preciceMesh :
           preciceSurfaceMeshes_)
        s << " \"" << preciceMesh->preciceMeshName << "\"";
      LOG(FATAL) << currentPreciceData << "[\"preciceMeshName\"] = \""
                 << currentMeshName
                 << "\" could not be found, available precice meshes: "
                 << s.str();
    }
    preciceData.preciceMesh = *iter;

    // parse mode and variable names
    std::string mode = currentPreciceData.getOptionString("mode", "");
    if (mode == "read-displacements-velocities") {
      preciceData.ioType = PreciceSurfaceData::ioRead;
      preciceData.boundaryConditionType = PreciceSurfaceData::bcTypeDirichlet;

      // get precice names of the variables
      preciceData.displacementsName = currentPreciceData.getOptionString(
          "displacementsName", "Displacement");
      preciceData.velocitiesName =
          currentPreciceData.getOptionString("velocitiesName", "Velocity");

      // get precice data ids
      preciceData.preciceDataDimDisplacements =
          preciceParticipant_->getDataDimensions(currentMeshName,
                                                 preciceData.displacementsName);

      preciceData.preciceDataDimVelocities =
          preciceParticipant_->getDataDimensions(currentMeshName,
                                                 preciceData.velocitiesName);
    } else if (mode == "read-traction") {
      preciceData.ioType = PreciceSurfaceData::ioRead;
      preciceData.boundaryConditionType = PreciceSurfaceData::bcTypeNeumann;

      // get precice names of the variables
      preciceData.tractionName =
          currentPreciceData.getOptionString("tractionName", "Traction");

      // get precice data ids
      preciceData.preciceDataDimTraction =
          preciceParticipant_->getDataDimensions(currentMeshName,
                                                 preciceData.tractionName);
    } else if (mode == "write-displacements-velocities") {
      preciceData.ioType = PreciceSurfaceData::ioWrite;

      // get precice names of the variables
      preciceData.displacementsName = currentPreciceData.getOptionString(
          "displacementsName", "Displacement");
      preciceData.velocitiesName =
          currentPreciceData.getOptionString("velocitiesName", "Velocity");

      // get precice data ids
      preciceData.preciceDataDimDisplacements =
          preciceParticipant_->getDataDimensions(currentMeshName,
                                                 preciceData.displacementsName);
      preciceData.preciceDataDimVelocities =
          preciceParticipant_->getDataDimensions(currentMeshName,
                                                 preciceData.velocitiesName);

    } else if (mode == "write-traction") {
      preciceData.ioType = PreciceSurfaceData::ioWrite;

      // get precice names of the variables
      preciceData.tractionName =
          currentPreciceData.getOptionString("tractionName", "Traction");

      // get precice data ids
      preciceData.preciceDataDimTraction =
          preciceParticipant_->getDataDimensions(currentMeshName,
                                                 preciceData.tractionName);
    } else if (mode == "write-averaged-traction") {
      preciceData.ioType = PreciceSurfaceData::ioWrite;
      preciceData.average = true;
      LOG(INFO) << "Set average = true";

      // get precice names of the variables
      preciceData.tractionName =
          currentPreciceData.getOptionString("tractionName", "Traction");

      // get precice data ids
      preciceData.preciceDataDimTraction =
          preciceParticipant_->getDataDimensions(currentMeshName,
                                                 preciceData.tractionName);
    } else {
      LOG(FATAL) << currentPreciceData << "[\"mode\"] is \"" << mode << "\", "
                 << "possible values are: \"read-displacements-velocities\", "
                    "\"read-traction\", \"write-displacements-velocities\", "
                    "\"write-traction\", \"write-averaged-traction\".";
    }

    // store preciceData to vector
    preciceSurfaceData_.push_back(preciceData);
  }
}

template <typename NestedSolver>
void PreciceAdapterInitialize<
    NestedSolver>::initializeDirichletBoundaryConditions() {
  using ElementWithNodesType =
      typename SpatialDiscretization::DirichletBoundaryConditionsBase<
          typename PreciceAdapterNestedSolver<NestedSolver>::FunctionSpace,
          6>::ElementWithNodes;

  // initialize dirichlet boundary conditions, set all Dirichlet boundary
  // condition values that will be needed later to vector (0,0,0)
  std::vector<ElementWithNodesType> dirichletBoundaryConditionElements;

  int nElementsX = functionSpace_->meshPartition()->nElementsLocal(0);
  int nElementsY = functionSpace_->meshPartition()->nElementsLocal(1);
  int nElementsZ = functionSpace_->meshPartition()->nElementsLocal(2);

  std::set<int> elementIndicesZ;

  // loop over precice field variables to be transferred, collect surface meshes
  // at bottom or top
  for (PreciceSurfaceData &preciceData : preciceSurfaceData_) {
    if (preciceData.ioType == PreciceSurfaceData::ioRead &&
        preciceData.boundaryConditionType ==
            PreciceSurfaceData::bcTypeDirichlet) {
      if (preciceData.preciceMesh->face == PreciceSurfaceMesh::face2Minus) {
        elementIndicesZ.insert(0);
      } else if (preciceData.preciceMesh->face ==
                 PreciceSurfaceMesh::face2Plus) {
        elementIndicesZ.insert(nElementsZ - 1);
      }
    }
  }

  // for either or both of top and bottom surface coupling mesh
  for (int elementIndexZ : elementIndicesZ) {
    int indexZ = 0;
    if (elementIndexZ > 0)
      indexZ = 2;

    // loop over elements
    for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++) {
      for (int elementIndexX = 0; elementIndexX < nElementsX; elementIndexX++) {
        ElementWithNodesType elementWithNodes;
        elementWithNodes.elementNoLocal =
            elementIndexZ * nElementsX * nElementsY +
            elementIndexY * nElementsX + elementIndexX;

        for (int indexY = 0; indexY < 3; indexY++) {
          for (int indexX = 0; indexX < 3; indexX++) {
            int elementalDofIndex = indexZ * 9 + indexY * 3 + indexX;
            elementWithNodes.elementalDofIndex.insert(std::pair<int, VecD<6>>(
                elementalDofIndex, VecD<6>{0, 0, 0, 0, 0, 0}));
          }
        }
        dirichletBoundaryConditionElements.push_back(elementWithNodes);
      }
    }
  }

  // add dirichlet bc values for all nodes that will get a prescribed value
  // during coupling
  this->addDirichletBoundaryConditions(nestedSolver_,
                                       dirichletBoundaryConditionElements);
}

// template <typename NestedSolver>
// void PreciceAdapterInitialize<NestedSolver>::initializePreciceVolumeData() {
//   // parse settings for coupling participants
//   // loop over items of the key "preciceData"
//   std::string settingsKey("preciceVolumeData");
//   PyObject *listPy = this->specificSettings_.getOptionPyObject(settingsKey);
//   std::vector<PyObject *> list =
//       PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);
//   PythonConfig preciceDataConfig(this->specificSettings_, settingsKey);

//   // get the slot connector data from the nested solver
//   int slotNo = 0;
//   using SlotConnectorDataType = typename NestedSolver::SlotConnectorDataType;
//   std::shared_ptr<SlotConnectorDataType> slotConnectorData =
//       nestedSolver_.getSlotConnectorData();

//   // get all present slot names
//   std::vector<std::string> slotNames;
//   SlotConnectorDataHelper<SlotConnectorDataType>::getSlotNames(
//       slotConnectorData, slotNames);

//   // loop over items of the list under "preciceData", the order of the items
//   is
//   // also the order of the data connector slots
//   for (int listIndex = 0; listIndex < list.size(); listIndex++) {
//     // extract the list item as PythonConfig {"mode": "write",
//     "variableName":
//     // ...}
//     PythonConfig currentPreciceData(preciceDataConfig, listIndex);

//     // create a new preciceData instance that will be added to
//     // preciceVolumeData_ at the end of this loop
//     PreciceVolumeData preciceData;

//     // parse options
//     preciceData.isGeometryField =
//         currentPreciceData.getOptionBool("isGeometryField", false);
//     preciceData.opendihuMeshName =
//         currentPreciceData.getOptionString("opendihuMeshName", "");

//     // parse slot name
//     preciceData.slotName = currentPreciceData.getOptionString("slotName",
//     ""); preciceData.slotNo = slotNo;

//     LOG(DEBUG) << "isGeometryField: " << preciceData.isGeometryField
//                << ", slotName: \"" << preciceData.slotName << "\".";

//     if (preciceData.slotName == "") {
//       preciceData.slotNo = slotNo;
//       LOG(DEBUG) << "Using slotNo " << slotNo
//                  << " because slotName was not given.";
//     } else {
//       LOG(DEBUG) << "The slotName \"" << preciceData.slotName << "\" is
//       given."; std::vector<std::string>::iterator slotNameIter =
//           std::find(slotNames.begin(), slotNames.end(),
//           preciceData.slotName);

//       // if the "from" and "to" slot names match to slot names in the current
//       // slotConnectorData
//       if (slotNameIter != slotNames.end()) {
//         LOG(DEBUG) << "Slot found";
//         // determine slot nos
//         preciceData.slotNo = std::distance(slotNames.begin(), slotNameIter);
//         LOG(DEBUG) << "Slot \"" << preciceData.slotName << "\" is slot no "
//                    << preciceData.slotNo;
//       } else {
//         LOG(DEBUG) << "Slot not found";
//         preciceData.slotNo = slotNo;
//         LOG(WARNING)
//             << "A slot with name \"" << preciceData.slotName
//             << "\" does not exist. Available slots: " << slotNames
//             << ". Use `None` or the empty string to specify the slot by its "
//                "number (referenced slots will be in the order of "
//             << "the items under \"preciceData\").\n Using slot No " <<
//             slotNo;
//       }
//     }

//     // parse the mesh
//     std::string currentMeshName =
//         currentPreciceData.getOptionString("preciceMeshName", "");

//     // find the mesh in the already parsed precice meshes
//     typename std::vector<std::shared_ptr<PreciceVolumeMesh>>::iterator iter =
//         std::find_if(
//             preciceVolumeMeshes_.begin(), preciceVolumeMeshes_.end(),
//             [&currentMeshName](std::shared_ptr<PreciceVolumeMesh>
//             preciceMesh) {
//               return preciceMesh->preciceMeshName == currentMeshName;
//             });

//     // if the mesh is not in preciceVolumeMeshes_, create it and add it to
//     // preciceVolumeMeshes_
//     if (iter == preciceVolumeMeshes_.end()) {
//       LOG(DEBUG) << "mesh was not yet initialized.";

//       // create new precice mesh object
//       std::shared_ptr<PreciceVolumeMesh> preciceMesh =
//           std::make_shared<PreciceVolumeMesh>();
//       preciceMesh->preciceMeshName = currentMeshName;

//       std::vector<Vec3> geometryValues;
//       int nArrayItems = 1;

//       // get the function space that is used to initialize the mapping
//       std::shared_ptr<typename NestedSolver::FunctionSpace> functionSpace =
//           nullptr;
//       if (preciceData.opendihuMeshName != "") {
//         functionSpace =
//             DihuContext::meshManager()
//                 ->functionSpace<typename NestedSolver::FunctionSpace>(
//                     preciceData.opendihuMeshName);
//         LOG(DEBUG) << "Using opendihu mesh with name \""
//                    << preciceData.opendihuMeshName
//                    << "\" to initialize precice mapping.";

//         // get the node positions of the opendihu mesh to initialize the
//         mapping preciceMesh->nNodesLocal =
//             functionSpace->geometryField().nDofsLocalWithoutGhosts();
//         functionSpace->geometryField().getValuesWithoutGhosts(geometryValues);

//         // store opendihu mesh name
//         preciceMesh->opendihuMeshName = functionSpace->meshName();

//         // the following cannot happen
//         if (preciceMesh->opendihuMeshName != preciceData.opendihuMeshName)
//           LOG(FATAL)
//               << "Initializing precice mesh from settings
//               opendihuMeshName=\""
//               << preciceData.opendihuMeshName << "\" and the resulting mesh
//               is "
//               << preciceMesh->opendihuMeshName;
//       } else {
//         LOG(DEBUG) << "get mesh from mesh partition, slot No "
//                    << preciceData.slotNo;
//         LOG(DEBUG) << "slot connector data type: "
//                    << StringUtility::demangle(
//                           typeid(SlotConnectorDataType).name());

//         // get the mesh partition
//         std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
//             SlotConnectorDataHelper<
//                 SlotConnectorDataType>::getMeshPartitionBase(slotConnectorData,
//                                                              preciceData.slotNo,
//                                                              0);

//         if (!meshPartitionBase) {
//           LOG(FATAL) << "Could not get mesh for slot No " <<
//           preciceData.slotNo;
//         } else
//           LOG(DEBUG) << "got mesh partition for slot No " <<
//           preciceData.slotNo;

//         // get opendihu mesh name
//         preciceMesh->opendihuMeshName =
//             SlotConnectorDataHelper<SlotConnectorDataType>::getMeshName(
//                 slotConnectorData, preciceData.slotNo);

//         int nDofsLocalWithoutGhosts =
//             meshPartitionBase->nDofsLocalWithoutGhosts();
//         nArrayItems =
//             SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
//                 slotConnectorData,
//                 preciceData.slotNo); // number of fibers if there are fibers

//         preciceMesh->nNodesLocal = nDofsLocalWithoutGhosts * nArrayItems;

//         // get the vector of values [0,1,...,nDofsLocalWithGhosts]
//         const std::vector<PetscInt> &dofNosLocalWithGhosts =
//             meshPartitionBase->dofNosLocal();
//         std::vector<PetscInt> dofNosLocalWithoutGhosts(
//             dofNosLocalWithGhosts.begin(),
//             dofNosLocalWithGhosts.begin() + nDofsLocalWithoutGhosts);

//         // loop over fibers if there are any
//         for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
//           static std::vector<Vec3> nodePositionsFiber;
//           nodePositionsFiber.clear();
//           SlotConnectorDataHelper<SlotConnectorDataType>::slotGetGeometryValues(
//               slotConnectorData, preciceData.slotNo, arrayIndex,
//               dofNosLocalWithoutGhosts, nodePositionsFiber);
//           geometryValues.insert(geometryValues.end(),
//                                 nodePositionsFiber.begin(),
//                                 nodePositionsFiber.end());
//         }

//         LOG(DEBUG) << "collected " << geometryValues.size()
//                    << " node positions from the " << nArrayItems << "
//                    fibers";
//       }

//       // transform to contiguous memory layout for precice
//       std::vector<double> geometryValuesContiguous(3 *
//       geometryValues.size()); for (int entryNo = 0; entryNo <
//       geometryValues.size(); entryNo++)
//         for (int componentNo = 0; componentNo < 3; componentNo++)
//           geometryValuesContiguous[3 * entryNo + componentNo] =
//               geometryValues[entryNo][componentNo];

//       if (geometryValuesContiguous.size() != preciceMesh->nNodesLocal * 3)
//         LOG(FATAL) << "size mismatch: " << geometryValuesContiguous.size()
//                    << "!=" << preciceMesh->nNodesLocal * 3;

//       // resize buffer for vertex ids
//       preciceMesh->preciceVertexIds.resize(preciceMesh->nNodesLocal);

//       // give the node positions to precice and get the vertex ids
//       preciceParticipant_->setMeshVertices(preciceMesh->preciceMeshName,
//                                            geometryValuesContiguous,
//                                            preciceMesh->preciceVertexIds);

//       // store newly created mesh in preciceData
//       preciceData.preciceMesh = preciceMesh;

//       // store the precice mesh to the vector of meshes
//       preciceVolumeMeshes_.push_back(preciceMesh);

//       // output message about mesh names
//       std::stringstream message;
//       message << "Initialized precice mesh \"" <<
//       preciceMesh->preciceMeshName
//               << "\" from opendihu mesh \"" << preciceMesh->opendihuMeshName
//               << "\" with " << preciceMesh->nNodesLocal << " local nodes";
//       if (nArrayItems > 1)
//         message << " (and other corresponding meshes, in total " <<
//         nArrayItems
//                 << " fibers/compartments)";
//       LOG(INFO) << message.str() << ".";
//     } else {
//       LOG(DEBUG) << "Use existing precice mesh " << currentMeshName;
//       preciceData.preciceMesh = *iter;
//     }

//     // parse mode
//     std::string mode = currentPreciceData.getOptionString("mode", "");
//     if (mode == "read") {
//       preciceData.ioType = PreciceVolumeData::ioRead;
//     } else if (mode == "write") {
//       preciceData.ioType = PreciceVolumeData::ioWrite;
//     } else {
//       LOG(FATAL) << currentPreciceData << "[\"mode\"] is \"" << mode << "\",
//       "
//                  << "possible values are: \"read\", \"write\".";
//     }

//     // parse variable name
//     preciceData.preciceDataName =
//         currentPreciceData.getOptionString("preciceDataName", "variable");

//     // increment slotNo, this value is used for slots that are not identified
//     by
//     // slotName in the config
//     if (!preciceData.isGeometryField)
//       slotNo++;

//     // check if mesh size matches the specified slot
//     // get the mesh partition
//     std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
//         SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase(
//             slotConnectorData, preciceData.slotNo, 0);

//     int nDofsLocalWithoutGhosts =
//     meshPartitionBase->nDofsLocalWithoutGhosts(); int nArrayItems =
//         SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
//             slotConnectorData,
//             preciceData.slotNo); // number of fibers if there are fibers

//     if (nDofsLocalWithoutGhosts * nArrayItems !=
//         preciceData.preciceMesh->nNodesLocal) {
//       LOG(DEBUG) << ", all available slots: "
//                  <<
//                  SlotConnectorDataHelper<SlotConnectorDataType>::getString(
//                         slotConnectorData);
//       LOG(FATAL)
//           << currentPreciceData
//           << ": Mesh does not match slot in PreciceAdapterVolumeCoupling.\n\n
//           "
//           << "The " << (listIndex + 1)
//           << (listIndex == 0
//                   ? "st"
//                   : (listIndex == 1 ? "nd" : (listIndex == 2 ? "rd" : "th")))
//           << " list item under \"preciceData\" uses slot " <<
//           preciceData.slotNo
//           << " (slotName \"" << preciceData.slotName << "\")"
//           << " and preciceMesh \"" <<
//           preciceData.preciceMesh->preciceMeshName
//           << "\".\nThe slot has " << nDofsLocalWithoutGhosts << " dofs * "
//           << nArrayItems
//           << " array items (fibers) = " << nDofsLocalWithoutGhosts *
//           nArrayItems
//           << " dofs in total. The mesh was initialized from opendihu mesh \""
//           << preciceData.preciceMesh->opendihuMeshName << "\" with "
//           << preciceData.preciceMesh->nNodesLocal << " nodes. \n"
//           << "(isGeometryField: " << std::boolalpha
//           << preciceData.isGeometryField << ")\n"
//           << "You can do the following:\n"
//           << "- Rename the preciceMesh to a unique name and update the
//           precice "
//              "config xml file accordingly.\n"
//           << "- Check the solver structure file to find out the meshes that "
//              "are associated with the slots.\n"
//           << "  Check that the correct mesh is specified and that the precice
//           "
//              "mesh was not initialized earlier with a different opendihu "
//              "mesh.\n"
//           << "  Every precice mesh is initialized the first time in appears "
//              "under \"preciceData\", sometimes reordering the entries can "
//              "help.\n"
//           << "  If the desired mesh is not available at the current solver, "
//              "maybe insert a MapDofs class.";
//     }

//     LOG(INFO) << "Precice data \"" << preciceData.preciceDataName
//               << "\" maps to "
//               << (preciceData.isGeometryField ? "the geometry field of " :
//               "")
//               << "slot " << preciceData.slotNo << " (\"" <<
//               preciceData.slotName
//               << "\") and uses precice mesh \""
//               << preciceData.preciceMesh->preciceMeshName
//               << "\", which is opendihu mesh \""
//               << preciceData.preciceMesh->opendihuMeshName << "\" with "
//               << preciceData.preciceMesh->nNodesLocal << " local nodes.";

//     // store preciceData to vector
//     preciceVolumeData_.push_back(preciceData);
//   }
// }
#endif

} // namespace Control
