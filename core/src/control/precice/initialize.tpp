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
bool PreciceAdapterInitialize<NestedSolver>::inMuscleMeshTopA(const double x,
                                                              const double y,
                                                              const double z) {
  double planeX = 10.019566831165198, planeY = 17.85257919747023,
         planeZ = -37.500429089326744; // Point on the plane
  double normalX = 0.9989900836636769, normalY = -0.04480513952143465,
         normalZ = -0.0033633635106413177; // Normal vector to the plane

  // Compute the vector from point on plane to the point in question
  double vx = x - planeX;
  double vy = y - planeY;
  double vz = z - planeZ;

  // Calculate dot product of this vector with the normal vector of the plane
  double d = vx * normalX + vy * normalY + vz * normalZ;

  // If d > 0, the point is on the right side of the plane (where 'right' is in
  // the direction of the normal vector)
  return d < 0;
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

  // int rankNo = functionSpace_->meshPartition()->rankSubset()->ownRankNo();
  // int nRanks = functionSpace_->meshPartition()->rankSubset()->size();

  // get the union of all MPI ranks that occur for any fiber
  std::shared_ptr<Partition::RankSubset> rankSubset =
      this->context_.partitionManager()->rankSubsetForCollectiveOperations();
  int rankNo = rankSubset->ownRankNo();
  int nRanks = rankSubset->size();
  LOG(INFO) << "rankNo = " << rankNo;
  LOG(INFO) << "nRanks = " << nRanks;

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
  if (functionSpace_->dimension() < 3) {
    return;
  }

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

          if (preciceMesh->preciceMeshName == "MuscleMeshTopA") {
            if (inMuscleMeshTopA(geometryValuesSurface[3 * surfaceDofNo],
                                 geometryValuesSurface[3 * surfaceDofNo + 1],
                                 geometryValuesSurface[3 * surfaceDofNo + 2])) {
              LOG(DEBUG) << "Mesh MuscleMeshTopA has surfaceDofNo = "
                         << surfaceDofNo << ", with coordinates x = "
                         << geometryValuesSurface[3 * surfaceDofNo] << ", y = "
                         << geometryValuesSurface[3 * surfaceDofNo + 1]
                         << ",  z = "
                         << geometryValuesSurface[3 * surfaceDofNo + 2];
              preciceMesh->selectedDofNosLocal.push_back(surfaceDofNo);
            }

          } else if (preciceMesh->preciceMeshName == "MuscleMeshTopB") {
            if (!inMuscleMeshTopA(
                    geometryValuesSurface[3 * surfaceDofNo],
                    geometryValuesSurface[3 * surfaceDofNo + 1],
                    geometryValuesSurface[3 * surfaceDofNo + 2])) {
              LOG(DEBUG) << "Mesh MuscleMeshTopB has surfaceDofNo = "
                         << surfaceDofNo << ", and dofNoLocal = " << dofNoLocal
                         << ", with coordinates x = "
                         << geometryValuesSurface[3 * surfaceDofNo] << ", y = "
                         << geometryValuesSurface[3 * surfaceDofNo + 1]
                         << ",  z = "
                         << geometryValuesSurface[3 * surfaceDofNo + 2];
              preciceMesh->selectedDofNosLocal.push_back(surfaceDofNo);
            }
          } else {
            preciceMesh->selectedDofNosLocal.push_back(surfaceDofNo);
            LOG(DEBUG) << "MeshBottom has surfaceDofNo = " << surfaceDofNo
                       << ", and dofNoLocal = " << dofNoLocal
                       << ", with coordinates x = "
                       << geometryValuesSurface[3 * surfaceDofNo] << ", y = "
                       << geometryValuesSurface[3 * surfaceDofNo + 1]
                       << ",  z = "
                       << geometryValuesSurface[3 * surfaceDofNo + 2];
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
  if (preciceSurfaceData_.empty()) {
    return;
  }
  if (functionSpace_->dimension() < 3) {
    return;
  }

  // initialize dirichlet boundary conditions, set all Dirichlet boundary
  // condition values that will be needed later to vector (0,0,0)
  std::vector<ElementWithNodesType> dirichletBoundaryConditionElements;

  int nElementsX = functionSpace_->meshPartition()->nElementsLocal(0);
  int nElementsY = functionSpace_->meshPartition()->nElementsLocal(1);
  int nElementsZ = functionSpace_->meshPartition()->nElementsLocal(2);

  bool hasFace2Minus = false;
  bool hasFace2Plus = false;
  // loop over precice field variables to be transferred, collect surface
  // meshes at bottom or top
  for (const PreciceSurfaceData &preciceData : preciceSurfaceData_) {
    if (preciceData.ioType == PreciceSurfaceData::ioRead &&
        preciceData.boundaryConditionType ==
            PreciceSurfaceData::bcTypeDirichlet) {
      if (preciceData.preciceMesh->face == PreciceSurfaceMesh::face2Minus) {
        hasFace2Minus = true;
      } else if (preciceData.preciceMesh->face ==
                 PreciceSurfaceMesh::face2Plus) {
        hasFace2Plus = true;
      }
    }
  }

  const int nRanks = functionSpace_->meshPartition()->nRanks(2);
  const int ownRankPartitioningIndexZ =
      functionSpace_->meshPartition()->ownRankPartitioningIndex(2);
  LOG(DEBUG) << "overall position index: " << ownRankPartitioningIndexZ
             << " | with nRanks: " << nRanks;
  if (ownRankPartitioningIndexZ == 0 && hasFace2Minus) {
    const int elementIndexZ = 0;
    const int indexZ = 0;
    for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++) {
      for (int elementIndexX = 0; elementIndexX < nElementsX; elementIndexX++) {
        ElementWithNodesType face2MinusElementWithNodes;
        face2MinusElementWithNodes.elementNoLocal =
            elementIndexZ * nElementsX * nElementsY +
            elementIndexY * nElementsX + elementIndexX;
        for (int indexY = 0; indexY < 3; indexY++) {
          for (int indexX = 0; indexX < 3; indexX++) {
            int elementalDofIndex = indexZ * 9 + indexY * 3 + indexX;
            face2MinusElementWithNodes.elementalDofIndex.insert(
                std::pair<int, VecD<6>>(elementalDofIndex,
                                        VecD<6>{0, 0, 0, 0, 0, 0}));
          }
        }
        LOG(DEBUG) << "dirichletBoundaryConditionElements: [" << elementIndexX
                   << "," << elementIndexY << "," << indexZ
                   << "]: " << face2MinusElementWithNodes.elementNoLocal
                   << " | " << face2MinusElementWithNodes.elementalDofIndex;
        dirichletBoundaryConditionElements.push_back(
            face2MinusElementWithNodes);
      }
    }
  }

  if (ownRankPartitioningIndexZ == (nRanks - 1) && hasFace2Plus) {
    const int elementIndexZ = nElementsZ - 1;
    const int indexZ = 2;
    for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++) {
      for (int elementIndexX = 0; elementIndexX < nElementsX; elementIndexX++) {
        ElementWithNodesType face2PlusElementWithNodes;
        face2PlusElementWithNodes.elementNoLocal =
            elementIndexZ * nElementsX * nElementsY +
            elementIndexY * nElementsX + elementIndexX;
        for (int indexY = 0; indexY < 3; indexY++) {
          for (int indexX = 0; indexX < 3; indexX++) {
            int elementalDofIndex = indexZ * 9 + indexY * 3 + indexX;
            face2PlusElementWithNodes.elementalDofIndex.insert(
                std::pair<int, VecD<6>>(elementalDofIndex,
                                        VecD<6>{0, 0, 0, 0, 0, 0}));
          }
        }
        LOG(DEBUG) << "dirichletBoundaryConditionElements: [" << elementIndexX
                   << "," << elementIndexY << "," << indexZ
                   << "]: " << face2PlusElementWithNodes.elementNoLocal << " | "
                   << face2PlusElementWithNodes.elementalDofIndex;
        dirichletBoundaryConditionElements.push_back(face2PlusElementWithNodes);
      }
    }
  }

  // add dirichlet bc values for all nodes that will get a prescribed value
  // during coupling
  this->addDirichletBoundaryConditions(nestedSolver_,
                                       dirichletBoundaryConditionElements);
}

#endif

} // namespace Control
