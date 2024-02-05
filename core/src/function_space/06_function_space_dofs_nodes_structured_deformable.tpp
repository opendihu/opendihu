#include "function_space/06_function_space_dofs_nodes.h"

#include <Python.h> // has to be the first included header
#include <cmath>
#include <array>

#include "easylogging++.h"
#include "field_variable/field_variable.h"
#include "utility/petsc_utility.h"
//#include "function_space/00_function_space_base_dim.h"
#include "mesh/face_t.h"

namespace FunctionSpace {

// constructor
template <int D, typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::StructuredDeformableOfDimension<D>,
                       BasisFunctionType>::
    FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager,
                           PythonConfig specificSettings, bool noGeometryField)
    : FunctionSpaceDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,
                                       BasisFunctionType>(partitionManager,
                                                          specificSettings) {
  LOG(DEBUG) << "constructor FunctionSpaceDofsNodes StructuredDeformable, "
                "noGeometryField_="
             << this->noGeometryField_;

  this->noGeometryField_ = noGeometryField;
}

// constructor
template <int D, typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::StructuredDeformableOfDimension<D>,
                       BasisFunctionType>::
    FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager,
                           std::vector<double> &localNodePositions,
                           PythonConfig specificSettings, bool noGeometryField)
    : FunctionSpaceDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,
                                       BasisFunctionType>(partitionManager,
                                                          specificSettings) {
  LOG(DEBUG) << "constructor FunctionSpaceDofsNodes StructuredDeformable, "
                "noGeometryField_="
             << this->noGeometryField_;

  // local node positions are without ghost nodes
  localNodePositions_ = localNodePositions;
  LOG(DEBUG) << "store " << localNodePositions_.size() << " node positions";

  this->noGeometryField_ = noGeometryField;
}

template <int D, typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::StructuredDeformableOfDimension<D>,
                       BasisFunctionType>::
    FunctionSpaceDofsNodes(
        std::shared_ptr<Partition::Manager> partitionManager,
        const std::vector<Vec3> &localNodePositions,
        const std::array<element_no_t, D> nElementsPerCoordinateDirectionLocal,
        const std::array<int, D> nRanksPerCoordinateDirection)
    : FunctionSpaceDofsNodesStructured<Mesh::StructuredDeformableOfDimension<D>,
                                       BasisFunctionType>(partitionManager,
                                                          NULL) {
  LOG(DEBUG) << "constructor FunctionSpaceDofsNodes StructuredDeformable, from "
             << localNodePositions.size() << " localNodePositions";

  this->noGeometryField_ = false;
  this->nElementsPerCoordinateDirectionLocal_ =
      nElementsPerCoordinateDirectionLocal;
  this->nRanks_ = nRanksPerCoordinateDirection;
  this->forcePartitioningCreationFromLocalNumberOfElements_ =
      true; // this is defined in 03_function_space_partition.h

  // forcePartitioningCreationFromLocalNumberOfElements_ is set to true, this
  // means that the partitioning is created considering
  // this->nElementsPerCoordinateDirectionLocal_ and not depending on values of
  // inputMeshIsGlobal

  LOG(DEBUG) << "set local number of elements per coordinate direction: "
             << this->nElementsPerCoordinateDirectionLocal_
             << ", nRanks: " << this->nRanks_;
  LOG(DEBUG)
      << "set forcePartitioningCreationFromLocalNumberOfElements_ to true";

  // local node positions are without ghost nodes
  localNodePositions_.reserve(localNodePositions.size() * D);

  for (const Vec3 &vector : localNodePositions) {
    for (int i = 0; i < 3; i++)
      localNodePositions_.push_back(vector[i]);
  }
}

template <int D, typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::StructuredDeformableOfDimension<D>,
                       BasisFunctionType>::
    FunctionSpaceDofsNodes(
        std::shared_ptr<Partition::Manager> partitionManager,
        const std::vector<double> &nodePositionsFromBinaryFile,
        const std::vector<Vec3> &localNodePositions,
        const std::array<element_no_t, D> nElementsPerCoordinateDirectionLocal,
        const std::array<int, D> nRanksPerCoordinateDirection)
    : FunctionSpaceDofsNodes<Mesh::StructuredDeformableOfDimension<D>,
                             BasisFunctionType>::
          FunctionSpaceDofsNodes(partitionManager, localNodePositions,
                                 nElementsPerCoordinateDirectionLocal,
                                 nRanksPerCoordinateDirection) {}

template <int D, typename BasisFunctionType>
void FunctionSpaceDofsNodes<Mesh::StructuredDeformableOfDimension<D>,
                            BasisFunctionType>::initialize() {
  // create meshPartition and redistribute elements if necessary, this needs
  // information about mesh size
  FunctionSpacePartition<Mesh::StructuredDeformableOfDimension<D>,
                         BasisFunctionType>::initialize();

  // if this mesh does not have a geometry field, do nothing further
  if (this->noGeometryField_)
    return;

  // setup geometry field
  // if no node positions were given, e.g. by the constructor that takes node
  // positions
  if (localNodePositions_.empty()) {
    // parse node positions from python config
    this->parseNodePositionsFromSettings(this->specificSettings_);
  }

  // pass a shared "this" pointer to the geometryField
  // retrieve "this" pointer and convert to downwards pointer of most derived
  // class "FunctionSpace"
  std::shared_ptr<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,
                                BasisFunctionType>>
      thisFunctionSpace = std::static_pointer_cast<FunctionSpace<
          Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>>(
          this->shared_from_this());

  assert(thisFunctionSpace != nullptr);

  // create empty field variable for geometry field
  std::vector<std::string> componentNames{"x", "y", "z"};
  this->geometryField_ = std::make_shared<GeometryFieldType>(
      thisFunctionSpace, "geometry", componentNames, true);

  // assign values of geometry field
  this->setGeometryFieldValues();

  // set initalized_ to true which indicates that initialize has been called
  this->initialized_ = true;
}

// read in config nodes
template <int D, typename BasisFunctionType>
void FunctionSpaceDofsNodes<Mesh::StructuredDeformableOfDimension<D>,
                            BasisFunctionType>::
    parseNodePositionsFromSettings(PythonConfig specificSettings) {
  LOG(TRACE) << "FunctionSpaceDofsNodes<structuredDeformable> "
                "parseNodePositionsFromSettings";

  // compute number of nodes
  node_no_t nNodesLocal = this->nNodesLocalWithoutGhosts();
  global_no_t nNodesGlobal = this->nNodesGlobal();
  global_no_t nNodes;

  // if the given information about the mesh in config is for the global mesh
  bool inputMeshIsGlobal =
      specificSettings.getOptionBool("inputMeshIsGlobal", true);
  LOG(DEBUG) << "inputMeshIsGlobal: " << std::boolalpha << inputMeshIsGlobal;

  global_no_t vectorSize;

  if (inputMeshIsGlobal) {
    vectorSize = nNodesGlobal * 3;
    nNodes = nNodesGlobal;
  } else {
    vectorSize = nNodesLocal * 3;
    nNodes = nNodesLocal;
  }

  VLOG(1) << "nNodesLocal: " << nNodesLocal
          << ", nNodesGlobal: " << nNodesGlobal
          << ", vectorSize: " << vectorSize;

  localNodePositions_.resize(
      vectorSize); // resize vector and value-initialize to 0
  // The vector localNodePositions_ is filled with all available node positions
  // in this method, for inputMeshIsGlobal these are the global values of the
  // global domain. At the end of this method, non-local entries are removed if
  // inputMeshIsGlobal.

  // fill initial position from settings
  if (specificSettings.hasKey("nodePositions")) {
    bool nodesStoredAsLists = false;
    VLOG(1) << "specificSettings has \"nodePositions\"";

    // check if the node positions are stored as list, e.g. [[x,y,z],
    // [x,y,z],...]
    PyObject *nodePositionsListPy =
        specificSettings.getOptionPyObject("nodePositions");
    if (PyList_Check(nodePositionsListPy)) {
      if (PyList_Size(nodePositionsListPy) > 0) {
        PyObject *firstItemPy =
            PyList_GetItem(nodePositionsListPy, (Py_ssize_t)0);
        if (PyList_Check(firstItemPy)) {
          nodesStoredAsLists = true;
        }
      }
    }

    VLOG(1) << "nodePositions: " << nodePositionsListPy
            << ", nodesStoredAsLists=" << nodesStoredAsLists;

    if (nodesStoredAsLists) {
      node_no_t nNodesInList = PyList_Size(nodePositionsListPy);

      if (nNodesInList != nNodes) {
        LOG(ERROR) << specificSettings.getStringPath()
                   << "[\"nodePositions\"]: Number of nodes in list ("
                   << nNodesInList << ") "
                   << "does not match expected number of the mesh (" << nNodes
                   << "). inputMeshIsGlobal is " << std::boolalpha
                   << inputMeshIsGlobal;
      }

      node_no_t nodeNo = 0;
      for (; nodeNo < nNodesInList; nodeNo++) {
        // extract single node position, e.g. [x,y]
        PyObject *itemNodePositionPy =
            PyList_GetItem(nodePositionsListPy, (Py_ssize_t)nodeNo);
        if (PyList_Check(itemNodePositionPy)) {
          // parse components of node, e.g. x
          int i = 0;
          for (; i < std::min(3, (int)PyList_Size(itemNodePositionPy)); i++) {
            PyObject *pointComponentPy =
                PyList_GetItem(itemNodePositionPy, (Py_ssize_t)i);
            localNodePositions_[3 * nodeNo + i] =
                PythonUtility::convertFromPython<double>::get(pointComponentPy,
                                                              0.0);
          }

          // set the rest of the values that were not specified to 0.0, e.g.
          // z=0.0
          for (; i < 3; i++) {
            localNodePositions_[3 * nodeNo + i] = 0.0;
          }

          VLOG(2) << "(1) set node " << nodeNo << "["
                  << localNodePositions_[3 * nodeNo + 0] << ","
                  << localNodePositions_[3 * nodeNo + 1] << ","
                  << localNodePositions_[3 * nodeNo + 2] << "]";
        } else {
          // if the entry is not a list like [x,y,z] but a single value, assume
          // it is the x value
          double value = PythonUtility::convertFromPython<double>::get(
              itemNodePositionPy, 0.0);
          localNodePositions_[3 * nodeNo + 0] = value;
          localNodePositions_[3 * nodeNo + 1] = 0.0;
          localNodePositions_[3 * nodeNo + 2] = 0.0;

          VLOG(2) << "(2) set node " << nodeNo << "["
                  << localNodePositions_[3 * nodeNo + 0] << ","
                  << localNodePositions_[3 * nodeNo + 1] << ","
                  << localNodePositions_[3 * nodeNo + 2] << "]";
        }
      }

      if (nodeNo < nNodes) {
        LOG(WARNING) << "Expected " << nNodes
                     << " nodes, localNodePositions_ contains only " << nodeNo
                     << " nodes (" << nNodes - nodeNo << " missing).";
      }

      // fill rest of values with 0,0,0
      for (; nodeNo < nNodes; nodeNo++) {
        localNodePositions_[3 * nodeNo + 0] = 0.0;
        localNodePositions_[3 * nodeNo + 1] = 0.0;
        localNodePositions_[3 * nodeNo + 2] = 0.0;

        VLOG(2) << "(3) set node " << nodeNo << "["
                << localNodePositions_[3 * nodeNo + 0] << ","
                << localNodePositions_[3 * nodeNo + 1] << ","
                << localNodePositions_[3 * nodeNo + 2] << "]";
      }
    } else {
      // nodes are stored as contiguous array, e.g. [x,y,z,x,y,z] or
      // [x,y,x,y,x,y,...]

      int nodeDimension = specificSettings.getOptionInt(
          "nodeDimension", 3, PythonUtility::ValidityCriterion::Between1And3);

      int inputVectorSize = nNodesLocal * nodeDimension;
      specificSettings.getOptionVector("nodePositions", inputVectorSize,
                                       localNodePositions_);

      LOG(DEBUG) << "nodeDimension: " << nodeDimension
                 << ", expect input vector to have " << nNodesLocal << "*"
                 << nodeDimension << "=" << inputVectorSize << " entries.";

      // transform vector from (x,y) or (x) entries to (x,y,z)
      if (nodeDimension < 3) {
        localNodePositions_.resize(
            vectorSize); // resize vector and value-initialize to 0
        for (int i = nNodesLocal - 1; i >= 0; i--) {

          if (nodeDimension == 2) {
            localNodePositions_[i * 3 + 1] =
                localNodePositions_[i * nodeDimension + 1];
          } else {
            localNodePositions_[i * 3 + 1] = 0;
          }
          localNodePositions_[i * 3 + 0] =
              localNodePositions_[i * nodeDimension + 0];
          localNodePositions_[i * 3 + 2] = 0;
        }
      }
    }
  } else // there was no "nodePositions" given in config, use physicalExtent
         // instead
  {
    // if node positions are not given in settings but physicalExtent, generate
    // node positions such that physicalExtent is reached
    std::array<double, D> physicalExtent, meshWidth, physicalOffset;
    physicalExtent = specificSettings.getOptionArray<double, D>(
        "physicalExtent", 1.0, PythonUtility::NonNegative);
    physicalOffset =
        specificSettings.getOptionArray<double, D>("physicalOffset", 0.0);

    for (unsigned int dimNo = 0; dimNo < D; dimNo++) {
      // depending on inputMeshIsGlobal interpret physicalExtent as global or
      // local physical extent
      if (inputMeshIsGlobal) {
        meshWidth[dimNo] =
            physicalExtent[dimNo] /
            (this->nElementsPerCoordinateDirectionGlobal(dimNo) *
             (FunctionSpaceBaseDim<1, BasisFunctionType>::nNodesPerElement() -
              1));
      } else {
        meshWidth[dimNo] =
            physicalExtent[dimNo] /
            (this->nElementsPerCoordinateDirectionLocal(dimNo) *
             (FunctionSpaceBaseDim<1, BasisFunctionType>::nNodesPerElement() -
              1));
      }
      LOG(DEBUG) << "meshWidth[" << dimNo << "] = " << meshWidth[dimNo];
    }

    VLOG(1) << "specificSettings has no \"nodePositions\", use physicalExtent: "
            << physicalExtent << ", meshWidth: " << meshWidth;

    std::array<double, 3> position{0., 0., 0.};

    // compute absolute node positions
    global_no_t nNodesInXDirection = 0;
    global_no_t nNodesInYDirection = 0;

    global_no_t offsetX = 0;
    global_no_t offsetY = 0;
    global_no_t offsetZ = 0;

    // initialize variables
    if (inputMeshIsGlobal) {
      nNodesInXDirection = this->nNodesGlobal(0);
      if (D >= 2) {
        nNodesInYDirection = this->nNodesGlobal(1);
      }
    } else {
      nNodesInXDirection = this->nNodesLocalWithoutGhosts(0);
      if (D >= 2) {
        nNodesInYDirection = this->nNodesLocalWithoutGhosts(1);
      }

      offsetX = this->meshPartition_->beginNodeGlobalNatural(0);
      if (D >= 2)
        offsetY = this->meshPartition_->beginNodeGlobalNatural(1);
      if (D >= 3)
        offsetZ = this->meshPartition_->beginNodeGlobalNatural(2);
    }

    global_no_t nodeX, nodeY, nodeZ; // helper variables

    LOG(DEBUG) << "nNodes: " << nNodes << ", global: " << nNodesGlobal
               << ", local: " << nNodesLocal << ", vectorSize: " << vectorSize;
    VLOG(1) << "offsetX: " << offsetX << ", offsetY: " << offsetY
            << ", offsetZ: " << offsetZ;

    for (global_no_t nodeNo = 0; nodeNo < nNodes; nodeNo++) {
      switch (D) {
      case 3:
        nodeZ = global_no_t(nodeNo / (nNodesInXDirection * nNodesInYDirection));
        position[2] = physicalOffset[2] + meshWidth[2] * (offsetZ + nodeZ);
      case 2:
        nodeY = global_no_t(nodeNo / nNodesInXDirection) % nNodesInYDirection;
        position[1] = physicalOffset[1] + meshWidth[1] * (offsetY + nodeY);
      case 1:
        nodeX = nodeNo % nNodesInXDirection;
        position[0] = physicalOffset[0] + meshWidth[0] * (offsetX + nodeX);
      }

      VLOG(1) << "position: " << position;

      // store the position values in nodePositions
      for (int i = 0; i < 3; i++) {
        localNodePositions_[nodeNo * 3 + i] = position[i];
      }
    } // nodeNo
  }

  // if parsed node positions in vector localNodePositions_ actually contains
  // global node positions, extract local positions
  if (inputMeshIsGlobal) {
    this->meshPartition_->extractLocalNodesWithoutGhosts(localNodePositions_,
                                                         3);
  }
}

// create geometry field from config nodes
template <int D, typename BasisFunctionType>
void FunctionSpaceDofsNodes<Mesh::StructuredDeformableOfDimension<D>,
                            BasisFunctionType>::setGeometryFieldValues() {
  LOG(DEBUG) << " Mesh StructuredDeformable, setGeometryField, size of "
                "nodePositions vector: "
             << localNodePositions_.size()
             << " (=" << localNodePositions_.size() / 3 << " node positions)";

  // compute number of (local) dofs
  dof_no_t nDofsLocal = this->nDofsLocalWithoutGhosts();

  // fill geometry vector from nodePositions, initialize non-node position
  // entries to 0 (for Hermite)
  std::vector<Vec3> geometryValues(nDofsLocal, Vec3{0.0});

  if (this->nNodesLocalWithoutGhosts() * 3 > localNodePositions_.size()) {
    std::stringstream dimensionsString;
    if (D >= 1)
      dimensionsString << this->meshPartition_->nNodesLocalWithoutGhosts(0);
    if (D >= 2)
      dimensionsString << "x"
                       << this->meshPartition_->nNodesLocalWithoutGhosts(1);
    if (D == 3)
      dimensionsString << "x"
                       << this->meshPartition_->nNodesLocalWithoutGhosts(2);

    LOG(FATAL) << "Given number of node positions ("
               << localNodePositions_.size() << ", i.e. "
               << localNodePositions_.size() / 3
               << " points) is smaller than required number ("
               << this->nNodesLocalWithGhosts() * 3 << ", i.e. "
               << this->nNodesLocalWithoutGhosts() << " points), for mesh \""
               << this->meshName()
               << "\", dimension: " << dimensionsString.str();
  }

  int geometryValuesIndex = 0;
  int nodePositionsIndex = 0;
  // loop over nodes
  for (node_no_t nodeNo = 0; nodeNo < this->nNodesLocalWithoutGhosts();
       nodeNo++) {
    // assign node position as first dof of the node
    geometryValues[geometryValuesIndex] =
        Vec3{localNodePositions_[nodePositionsIndex + 0],
             localNodePositions_[nodePositionsIndex + 1],
             localNodePositions_[nodePositionsIndex + 2]};
    geometryValuesIndex++;
    nodePositionsIndex += 3;

    // set entries to 0 for rest of dofs at this node (derivatives). These
    // values are adjusted later for Hermite in setHermiteDerivatives
    for (int nodalDofIndex = 1; nodalDofIndex < this->nDofsPerNode();
         nodalDofIndex++) {
      geometryValuesIndex++;
    }
  }
  // set values for node positions as geometry field
  this->geometryField_->setValuesWithoutGhosts(geometryValues);

  // initialize Hermite derivative dofs such that geometry fields becomes "even"
  bool setHermiteDerivatives = false;
  if (std::is_same<BasisFunctionType, BasisFunction::Hermite>::value) {
    setHermiteDerivatives =
        this->specificSettings_.getOptionBool("setHermiteDerivatives", true);
  }
  if (setHermiteDerivatives) {
    this->setHermiteDerivatives();
  }

  // this->geometryField_->finishGhostManipulation();      // reduce ghost
  // values, not necessary
  this->geometryField_->setRepresentationGlobal();
  this->geometryField_->startGhostManipulation(); // distribute ghost values

  // output ghost values for debugging
#ifndef NDEBUG
  // get ghost values
  std::stringstream stream;
  std::vector<Vec3> geometryFieldValuesWithGhosts;
  this->geometryField().getValuesWithGhosts(geometryFieldValuesWithGhosts);
  for (int i = this->meshPartition_->nDofsLocalWithoutGhosts();
       i < geometryFieldValuesWithGhosts.size(); i++) {
    stream << " " << geometryFieldValuesWithGhosts[i];
  }
  LOG(DEBUG) << "in setGeometryFieldValues, geometry field ghost values: "
             << stream.str();
#endif

  VLOG(1) << "setGeometryField, geometryValues: " << geometryValues;
}

template <int D, typename BasisFunctionType>
void FunctionSpaceDofsNodes<
    Mesh::StructuredDeformableOfDimension<D>,
    BasisFunctionType>::refineMesh(std::array<int, D> refinementFactors) {
  LOG(DEBUG) << "refineMesh with refinementFactors: " << refinementFactors;

  std::vector<Vec3> oldGeometryValues;
  this->geometryField_->getValuesWithGhosts(oldGeometryValues);

  int averageNNodesPerElement1D =
      FunctionSpaceBaseDim<1, BasisFunctionType>::averageNNodesPerElement();

  // determine number of old nodes
  std::array<int, 3> nElementsLocalOld({1, 1, 1});
  std::array<int, 3> nElementsLocalNew({1, 1, 1});
  std::array<int, 3> nNodesWithGhostsOld({1, 1, 1});
  std::array<int, 3> nNodesNew({1, 1, 1});

  for (int i = 0; i < D; i++) {
    nElementsLocalOld[i] = this->meshPartition_->nElementsLocal(i);
    nElementsLocalNew[i] = nElementsLocalOld[i] * refinementFactors[i];
    nNodesNew[i] = nElementsLocalNew[i] * averageNNodesPerElement1D +
                   (this->meshPartition_->hasFullNumberOfNodes(i) ? 1 : 0);
    nNodesWithGhostsOld[i] = this->meshPartition_->nNodesLocalWithGhosts(i);
  }

  int nNodesLocalWithoutGhostsNew = nNodesNew[0] * nNodesNew[1] * nNodesNew[2];
  localNodePositions_.resize(nNodesLocalWithoutGhostsNew * 3);

  LOG(DEBUG) << "nElementsLocal: " << nElementsLocalOld << " -> "
             << nElementsLocalNew;
  LOG(DEBUG) << "oldGeometryValues.size: " << oldGeometryValues.size();
  LOG(DEBUG) << "nNodesNew: " << nNodesNew
             << ", nNodesWithGhostsOld: " << nNodesWithGhostsOld;

  // get local natural dof nos (ordering including ghost dofs)
  std::shared_ptr<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,
                                BasisFunctionType>>
      thisFunctionSpace = std::static_pointer_cast<FunctionSpace<
          Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>>(
          this->shared_from_this());

  const std::vector<dof_no_t> dofNosOldLocalNaturalOrdering(
      this->meshPartition_->dofNosLocalNaturalOrdering().begin(),
      this->meshPartition_->dofNosLocalNaturalOrdering().end());

  this->meshPartition_->refine(refinementFactors);
  for (int i = 0; i < D; i++) {
    this->nElementsPerCoordinateDirectionLocal_[i] =
        this->meshPartition_->nElementsLocal(i);
    this->nElementsPerCoordinateDirectionGlobal_[i] =
        this->meshPartition_->nElementsGlobal(i);
  }

  LOG(DEBUG) << "nNodesLocalWithoutGhosts: "
             << this->nNodesLocalWithoutGhosts();
  LOG(DEBUG) << "nNodesLocalWithoutGhosts: "
             << this->meshPartition_->nNodesLocalWithoutGhosts(0) << ","
             << this->meshPartition_->nNodesLocalWithoutGhosts(1) << ","
             << this->meshPartition_->nNodesLocalWithoutGhosts(2) << ": "
             << this->meshPartition_->nNodesLocalWithoutGhosts();

  std::array<std::array<std::array<Vec3, 2>, 2>, 2> neighbouringPoints({});
  LOG(DEBUG) << "set new node positions";
  LOG(DEBUG) << "dofNosOldLocalNaturalOrdering:"
             << dofNosOldLocalNaturalOrdering;

  // loop over new node positions
  for (int zIndexNew = 0; zIndexNew < nNodesNew[2]; zIndexNew++) {
    for (int yIndexNew = 0; yIndexNew < nNodesNew[1]; yIndexNew++) {
      for (int xIndexNew = 0; xIndexNew < nNodesNew[0]; xIndexNew++) {

        std::array<double, 3> alpha({1.0, 1.0, 1.0});

        alpha[0] =
            double(xIndexNew % refinementFactors[0]) / refinementFactors[0];
        alpha[1] =
            double(yIndexNew % refinementFactors[1]) / refinementFactors[1];
        alpha[2] =
            double(zIndexNew % refinementFactors[2]) / refinementFactors[2];

        // get neighbouring node positions
        int xIndexOld = xIndexNew / refinementFactors[0];
        int yIndexOld = yIndexNew / refinementFactors[1];
        int zIndexOld = zIndexNew / refinementFactors[2];

        VLOG(1) << xIndexNew << "," << yIndexNew << "," << zIndexNew << " / "
                << nNodesNew << ", alpha: " << alpha
                << ", old index: " << xIndexOld << "," << yIndexOld << ","
                << zIndexOld;

        // z- y- x-
        int indexOld =
            zIndexOld * nNodesWithGhostsOld[1] * nNodesWithGhostsOld[0] +
            yIndexOld * nNodesWithGhostsOld[0] + xIndexOld;
        if (indexOld >= dofNosOldLocalNaturalOrdering.size()) {
          LOG(DEBUG) << "indexOld: " << indexOld
                     << ", dofNosOldLocalNaturalOrdering.size: "
                     << dofNosOldLocalNaturalOrdering.size();
          LOG(DEBUG) << "dofNosOldLocalNaturalOrdering: "
                     << dofNosOldLocalNaturalOrdering;
        }
        assert(indexOld < dofNosOldLocalNaturalOrdering.size());
        int dofNoOld = dofNosOldLocalNaturalOrdering[indexOld];
        if (dofNoOld >= oldGeometryValues.size()) {
          LOG(DEBUG) << "(x,y,z)=(" << xIndexNew << "," << yIndexNew << ","
                     << zIndexNew << "), old:(" << xIndexOld << "," << yIndexOld
                     << "," << zIndexOld << ") dofNoOld: " << dofNoOld
                     << " (indexOld: " << indexOld
                     << "), oldGeometryValues.size(): "
                     << oldGeometryValues.size();
        }
        assert(dofNoOld < oldGeometryValues.size());
        neighbouringPoints[0][0][0] = oldGeometryValues[dofNoOld];
        Vec3 resultingPoint = (1.0 - alpha[2]) * (1.0 - alpha[1]) *
                              (1.0 - alpha[0]) * neighbouringPoints[0][0][0];

        // z- y- x+
        indexOld = zIndexOld * nNodesWithGhostsOld[1] * nNodesWithGhostsOld[0] +
                   yIndexOld * nNodesWithGhostsOld[0] + (xIndexOld + 1);
        if (indexOld < dofNosOldLocalNaturalOrdering.size()) {
          dofNoOld = dofNosOldLocalNaturalOrdering[indexOld];
          if (dofNoOld < oldGeometryValues.size()) {
            neighbouringPoints[0][0][1] = oldGeometryValues[dofNoOld];
            resultingPoint += (1.0 - alpha[2]) * (1.0 - alpha[1]) * alpha[0] *
                              neighbouringPoints[0][0][1];
          }
        }

        // z- y+ x-
        indexOld = zIndexOld * nNodesWithGhostsOld[1] * nNodesWithGhostsOld[0] +
                   (yIndexOld + 1) * nNodesWithGhostsOld[0] + xIndexOld;
        if (indexOld < dofNosOldLocalNaturalOrdering.size()) {
          dofNoOld = dofNosOldLocalNaturalOrdering[indexOld];
          if (dofNoOld < oldGeometryValues.size()) {
            neighbouringPoints[0][1][0] = oldGeometryValues[dofNoOld];
            resultingPoint += (1.0 - alpha[2]) * alpha[1] * (1.0 - alpha[0]) *
                              neighbouringPoints[0][1][0];
          }
        }

        // z- y+ x+
        indexOld = zIndexOld * nNodesWithGhostsOld[1] * nNodesWithGhostsOld[0] +
                   (yIndexOld + 1) * nNodesWithGhostsOld[0] + (xIndexOld + 1);
        if (indexOld < dofNosOldLocalNaturalOrdering.size()) {
          dofNoOld = dofNosOldLocalNaturalOrdering[indexOld];
          if (dofNoOld < oldGeometryValues.size()) {
            neighbouringPoints[0][1][1] = oldGeometryValues[dofNoOld];
            resultingPoint += (1.0 - alpha[2]) * alpha[1] * alpha[0] *
                              neighbouringPoints[0][1][1];
          }
        }

        // z+ y- x-
        indexOld =
            (zIndexOld + 1) * nNodesWithGhostsOld[1] * nNodesWithGhostsOld[0] +
            yIndexOld * nNodesWithGhostsOld[0] + xIndexOld;
        if (indexOld < dofNosOldLocalNaturalOrdering.size()) {
          dofNoOld = dofNosOldLocalNaturalOrdering[indexOld];
          if (dofNoOld < oldGeometryValues.size()) {
            neighbouringPoints[1][0][0] = oldGeometryValues[dofNoOld];
            resultingPoint += alpha[2] * (1.0 - alpha[1]) * (1.0 - alpha[0]) *
                              neighbouringPoints[1][0][0];
          }
        }

        // z+ y- x+
        indexOld =
            (zIndexOld + 1) * nNodesWithGhostsOld[1] * nNodesWithGhostsOld[0] +
            yIndexOld * nNodesWithGhostsOld[0] + (xIndexOld + 1);
        if (indexOld < dofNosOldLocalNaturalOrdering.size()) {
          dofNoOld = dofNosOldLocalNaturalOrdering[indexOld];
          if (dofNoOld < oldGeometryValues.size()) {
            neighbouringPoints[1][0][1] = oldGeometryValues[dofNoOld];
            resultingPoint += alpha[2] * (1.0 - alpha[1]) * alpha[0] *
                              neighbouringPoints[1][0][1];
          }
        }

        // z+ y+ x-
        indexOld =
            (zIndexOld + 1) * nNodesWithGhostsOld[1] * nNodesWithGhostsOld[0] +
            (yIndexOld + 1) * nNodesWithGhostsOld[0] + xIndexOld;
        if (indexOld < dofNosOldLocalNaturalOrdering.size()) {
          dofNoOld = dofNosOldLocalNaturalOrdering[indexOld];
          if (dofNoOld < oldGeometryValues.size()) {
            neighbouringPoints[1][1][0] = oldGeometryValues[dofNoOld];
            resultingPoint += alpha[2] * alpha[1] * (1.0 - alpha[0]) *
                              neighbouringPoints[1][1][0];
          }
        }

        // z+ y+ x+
        indexOld =
            (zIndexOld + 1) * nNodesWithGhostsOld[1] * nNodesWithGhostsOld[0] +
            (yIndexOld + 1) * nNodesWithGhostsOld[0] + (xIndexOld + 1);
        if (indexOld < dofNosOldLocalNaturalOrdering.size()) {
          dofNoOld = dofNosOldLocalNaturalOrdering[indexOld];
          if (dofNoOld < oldGeometryValues.size()) {
            neighbouringPoints[1][1][1] = oldGeometryValues[dofNoOld];
            resultingPoint +=
                alpha[2] * alpha[1] * alpha[0] * neighbouringPoints[1][1][1];
          }
        }

        int indexNew = zIndexNew * nNodesNew[1] * nNodesNew[0] +
                       yIndexNew * nNodesNew[0] + xIndexNew;
        assert(3 * indexNew + 2 < localNodePositions_.size());

        localNodePositions_[3 * indexNew + 0] = resultingPoint[0];
        localNodePositions_[3 * indexNew + 1] = resultingPoint[1];
        localNodePositions_[3 * indexNew + 2] = resultingPoint[2];

        VLOG(1) << "    neighbouringPoints: " << neighbouringPoints
                << ", resultingPoint: " << resultingPoint;
      }
    }
  }

  LOG(DEBUG) << "create new geometry field";

  // create new geometry field
  std::vector<std::string> componentNames{"x", "y", "z"};
  this->geometryField_ = std::make_shared<GeometryFieldType>(
      thisFunctionSpace, "geometry", componentNames, true);

  LOG(DEBUG) << "setGeometryFieldValues";

  // assign new values of geometry field
  this->setGeometryFieldValues();
}

} // namespace FunctionSpace
