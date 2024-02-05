#include "spatial_discretization/dirichlet_boundary_conditions/00_dirichlet_boundary_conditions_base.h"

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "utility/vector_operators.h"
#include "control/types.h"

namespace SpatialDiscretization {

template <typename FunctionSpaceType, int nComponents>
DirichletBoundaryConditionsBase<FunctionSpaceType, nComponents>::
    DirichletBoundaryConditionsBase(DihuContext context)
    : specificSettings_(NULL) {
  LOG(DEBUG) << " create new dirichlet boundary conditions object.";
}

template <typename FunctionSpaceType, int nComponents>
void DirichletBoundaryConditionsBase<FunctionSpaceType, nComponents>::
    initialize(PythonConfig specificSettings,
               std::shared_ptr<FunctionSpaceType> functionSpace,
               std::string boundaryConditionsConfigKey) {
  functionSpace_ = functionSpace;

  LOG(DEBUG)
      << "DirichletBoundaryConditionsBase::initialize, use functionSpace_: "
      << functionSpace_->meshName();

  specificSettings_ = specificSettings;
  printDebuggingInfo();

  // clear previous boundary condition values
  boundaryConditionElements_.clear();
  boundaryConditionNonGhostDofLocalNos_.clear();
  boundaryConditionValues_.clear();

  // parse new boundary condition dofs and values
  parseBoundaryConditionsForElements(boundaryConditionsConfigKey);

  // fill auxiliary ghost element data structures, this is only needed for
  // Dirichlet boundary conditions on scalar fields
  this->initializeGhostElements();

  // write the constraint dofs and the prescribed values as points in a vtk file
  writeOutput();
}

template <typename FunctionSpaceType, int nComponents>
void DirichletBoundaryConditionsBase<FunctionSpaceType, nComponents>::
    initialize(std::shared_ptr<FunctionSpaceType> functionSpace,
               std::vector<DirichletBoundaryConditionsBase<
                   FunctionSpaceType, nComponents>::ElementWithNodes>
                   &boundaryConditionElements,
               std::vector<dof_no_t> &boundaryConditionNonGhostDofLocalNos,
               std::vector<ValueType> &boundaryConditionValues) {
  functionSpace_ = functionSpace;

  LOG(DEBUG)
      << "DirichletBoundaryConditionsBase::initialize, use functionSpace_: "
      << functionSpace_->meshName();

  boundaryConditionElements_ = boundaryConditionElements;
  boundaryConditionNonGhostDofLocalNos_ = boundaryConditionNonGhostDofLocalNos;
  boundaryConditionValues_ = boundaryConditionValues;

  // create the boundaryConditionsByComponent_ data structure from
  // boundaryConditionNonGhostDofLocalNos_ and boundaryConditionValues_
  generateBoundaryConditionsByComponent();
  printDebuggingInfo();

  // fill auxiliary ghost element data structures, this is only needed for
  // Dirichlet boundary conditions on scalar fields
  this->initializeGhostElements();
}

template <typename FunctionSpaceType, int nComponents>
void DirichletBoundaryConditionsBase<FunctionSpaceType,
                                     nComponents>::printDebuggingInfo() {
  if (VLOG_IS_ON(1)) {
    VLOG(1) << "parsed boundary conditions:";
    VLOG(1) << "  dofsLocal of BC: " << boundaryConditionNonGhostDofLocalNos_
            << " with prescribed values: " << boundaryConditionValues_;

    for (auto boundaryConditionElement : boundaryConditionElements_) {
      VLOG(1) << "  elementNo: " << boundaryConditionElement.elementNoLocal
              << " has (dof,value): "
              << boundaryConditionElement.elementalDofIndex;
    }
  }
}

template <typename FunctionSpaceType, int nComponents>
void DirichletBoundaryConditionsBase<FunctionSpaceType, nComponents>::
    parseBoundaryConditions(
        PythonConfig settings, std::shared_ptr<FunctionSpaceType> functionSpace,
        std::string boundaryConditionsConfigKey,
        std::vector<std::pair<int, std::array<double, nComponents>>>
            &boundaryConditions) {
  LOG(TRACE) << "parseBoundaryConditions, key \"" << boundaryConditionsConfigKey
             << "\"";
  assert(functionSpace);

  if (VLOG_IS_ON(1)) {
    PythonUtility::printDict(settings.pyObject());
  }

  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // determine if the BC indices in the config are given for global or local dof
  // nos
  bool inputMeshIsGlobal = settings.getOptionBool("inputMeshIsGlobal", true);

  int nDofs = 0;
  if (inputMeshIsGlobal) {
    // For composite meshes, the maximum number of dofs is the sum of the
    // numbers of dofs of each sub mesh, i.e. counting the shared nodes multiple
    // times. For all other meshes it equals nDofsGlobal.
    nDofs = functionSpace->meshPartition()->nDofsGlobalForBoundaryConditions();
  } else {
    nDofs = functionSpace->nDofsLocalWithoutGhosts();
  }

  if (settings.hasKey("DirichletBoundaryCondition")) {
    LOG(ERROR) << "Option \"DirichletBoundaryCondition\" was renamed to "
                  "\"dirichletBoundaryConditions\".";
  }

  VLOG(1) << "inputMeshIsGlobal: " << inputMeshIsGlobal << ", nDofs: " << nDofs;

  // parse all boundary conditions that are given in config
  std::pair<int, std::array<double, nComponents>> boundaryCondition =
      settings.getOptionDictBegin<int, std::array<double, nComponents>>(
          boundaryConditionsConfigKey);
  for (; !settings.getOptionDictEnd(boundaryConditionsConfigKey);
       settings.getOptionDictNext<int, std::array<double, nComponents>>(
           boundaryConditionsConfigKey, boundaryCondition)) {
    // for negative indices add number of dofs such that -1 is the last dof, -2
    // is the second-last etc.
    if (boundaryCondition.first < 0) {
      boundaryCondition.first += nDofs;
    } else if (boundaryCondition.first > nDofs) {
      node_no_t nodeNoLocal = boundaryCondition.first / nDofsPerNode;
      node_no_t nNodesLocal = functionSpace->nNodesLocalWithoutGhosts();
      LOG(ERROR) << "Dirichlet boundary condition specified for index "
                 << boundaryCondition.first << " (node " << nodeNoLocal << "), "
                 << "but there are only "
                 << functionSpace->nDofsLocalWithoutGhosts() << " local dofs ("
                 << nNodesLocal << " nodes).\nFunctionSpaceType: "
                 << StringUtility::demangle(typeid(FunctionSpaceType).name());
    }

    boundaryConditions.push_back(boundaryCondition);
  }
  if (boundaryConditions.empty()) {
    LOG(DEBUG) << "no boundary conditions specified";
  }
}

template <typename FunctionSpaceType, int nComponents>
void DirichletBoundaryConditionsBase<FunctionSpaceType, nComponents>::
    parseBoundaryConditionsForElements(
        std::string boundaryConditionsConfigKey) {
  LOG(TRACE) << "parseBoundaryConditionsForElements, functionSpace: "
             << functionSpace_->meshName()
             << ", nElementsLocal: " << functionSpace_->nElementsLocal();

  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // determine if the BC indices in the config are given for global or local dof
  // nos
  bool inputMeshIsGlobal =
      this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);

  LOG(DEBUG) << "boundaryConditionNonGhostDofLocalNos: "
             << boundaryConditionNonGhostDofLocalNos_
             << ", boundaryConditionValues: " << boundaryConditionValues_;

  // Boundary conditions are specified for dof numbers, not nodes, such that for
  // Hermite it is possible to prescribe derivatives. However the ordering of
  // the dofs is not known in the config for unstructured meshes. Therefore the
  // ordering is special. For every node there are as many values as dofs, in
  // contiguous order. Example for 2D Hermite, unstructured grid, 2x2 elements:
  //
  // node numbering:
  //  6_7_8
  // 3|_4_|5
  // 0|_1_|2
  //
  // dof numbering:
  //  6_7_8
  // 2|_3_|5
  // 0|_1_|4
  //
  // To specify du/dn = 0 at the left boundary you would set:
  // bc[0*2+1] = 0, bc[3*2+1] = 0, bc[6*2+1] = 0
  //
  // To specifiy u=0 on the bottom, you would set:
  // bc[0] = 0, bc[2] = 0, bc[4] = 0

  // ValueType is the array with number of components:
  // std::array<double,nComponents>
  std::vector<std::pair<int, ValueType>> boundaryConditions; // (index, value)
  parseBoundaryConditions(this->specificSettings_, functionSpace_,
                          boundaryConditionsConfigKey, boundaryConditions);

  LOG(DEBUG) << "read in Dirichlet boundary conditions from config (None means "
                "there is no BC for this dof, "
             << "the internal representation for this is max double (\""
             << std::numeric_limits<double>::max() << "\")). ";
  LOG(DEBUG) << boundaryConditions;

  // sort all parsed boundary conditions for their index no
  auto compareFunction = [](const std::pair<int, ValueType> &item1,
                            const std::pair<int, ValueType> &item2) {
    return item1.first < item2.first;
  };
  std::sort(boundaryConditions.begin(), boundaryConditions.end(),
            compareFunction);

  LOG(DEBUG) << "sorted Dirichlet boundary conditions from config: "
             << boundaryConditions;

  // determine elements with nodes that have prescribed boundary conditions,
  // store them in the vector boundaryConditionElements_, which is organized by
  // local elements
  element_no_t lastBoundaryConditionElement = -1;
  std::set<dof_no_t>
      boundaryConditionNonGhostDofLocalNosSet; //< same data as in
                                               // boundaryConditionNonGhostDofLocalNos_,
                                               // but as set

  std::vector<std::pair<int, ValueType>>
      boundaryConditionsNonGhost; //< boundary condition dof nos and values for
                                  // non-ghost dofs, used to collect values and
                                  // sort them afterwards

  // loop over all local elements
  for (element_no_t elementNoLocal = 0;
       elementNoLocal < functionSpace_->nElementsLocal(); elementNoLocal++) {
    VLOG(2) << "elementNoLocal: " << elementNoLocal;

    // loop over the nodes of this element
    int elementalDofIndex = 0;
    for (int nodeIndex = 0; nodeIndex < FunctionSpaceType::nNodesPerElement();
         nodeIndex++) {
      // get global or local nodeNo of current node (depending on
      // inputMeshIsGlobal)
      global_no_t nodeNo = 0;
      if (inputMeshIsGlobal) {
        nodeNo = functionSpace_->getNodeNoGlobalNaturalFromElementNoLocal(
            elementNoLocal, nodeIndex);
        VLOG(2) << "elementNoLocal: " << elementNoLocal
                << ", nodeIndex: " << nodeIndex
                << " -> global nodeNo: " << nodeNo;
      } else {
        // we assume that ghost dofs were also specified in the config
        nodeNo = functionSpace_->getNodeNo(elementNoLocal, nodeIndex);
      }

      // loop over dofs of node and thus over the elemental dofs
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode;
           nodalDofIndex++, elementalDofIndex++) {
        global_no_t indexForBoundaryCondition =
            nodeNo * nDofsPerNode + nodalDofIndex;

        // check if a boundary condition for the current node is given
        bool boundaryConditionForDofFound = false;

        // loop over all stored boundary conditions
        ValueType boundaryConditionValue;
        for (typename std::vector<std::pair<int, ValueType>>::const_iterator
                 boundaryConditionIter = boundaryConditions.begin();
             boundaryConditionIter != boundaryConditions.end();
             boundaryConditionIter++) {
          if (boundaryConditionIter->first == indexForBoundaryCondition) {
            boundaryConditionValue = boundaryConditionIter->second;
            boundaryConditionForDofFound = true;
            break;
          }
          if (boundaryConditionIter->first >= indexForBoundaryCondition)
            break;
        }

        VLOG(2) << "boundaryConditionForDofFound: "
                << boundaryConditionForDofFound;

        // here boundaryConditionIter->first is equal to
        // indexForBoundaryCondition (then the current element/node/dof matches
        // the boundary condition) or boundaryConditionIter->first >
        // indexForBoundaryCondition then the current dofIndex does not have any
        // boundary condition or boundaryConditionIter ==
        // boundaryConditions.end() then, too

        // if the currently considered boundaryCondition entry from config
        // matches the current nodeNo and nodalDofIndex
        if (boundaryConditionForDofFound) {
          VLOG(2) << "elementNoLocal: " << elementNoLocal
                  << ", lastBoundaryConditionElement: "
                  << lastBoundaryConditionElement;

          // if there is not yet an entry in boundaryConditionElements with the
          // current element, create one
          if (elementNoLocal != lastBoundaryConditionElement) {
            // add current element
            boundaryConditionElements_.emplace_back();
            boundaryConditionElements_.back().elementNoLocal = elementNoLocal;

            lastBoundaryConditionElement = elementNoLocal;

            VLOG(2) << "add empty entry for elementNoLocal " << elementNoLocal;
          }

          // add current node and boundary condition value to list of boundary
          // conditions for current element
          ElementWithNodes &boundaryConditionElement =
              boundaryConditionElements_.back();

          VLOG(2) << " add (el-dof, value)"
                  << std::pair<int, ValueType>(elementalDofIndex,
                                               boundaryConditionValue)
                  << ", to boundaryConditionElement of element "
                  << boundaryConditionElement.elementNoLocal;
          boundaryConditionElement.elementalDofIndex.insert(
              std::pair<int, ValueType>(elementalDofIndex,
                                        boundaryConditionValue));

          // also store localDofNo
          dof_no_t dofLocalNo =
              functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);

          // if the dof is non-ghost
          if (dofLocalNo < functionSpace_->nDofsLocalWithoutGhosts()) {
            // if dofLocalNo is not already contained in
            // boundaryConditionsNonGhost
            if (boundaryConditionNonGhostDofLocalNosSet.find(dofLocalNo) ==
                boundaryConditionNonGhostDofLocalNosSet.end()) {
              boundaryConditionNonGhostDofLocalNosSet.insert(dofLocalNo);
              boundaryConditionsNonGhost.push_back(std::pair<int, ValueType>(
                  dofLocalNo, boundaryConditionValue));
            }
          }
        }
      }
    }
  }

  // sort list according to dof no.s
  std::sort(boundaryConditionsNonGhost.begin(),
            boundaryConditionsNonGhost.end(),
            [&](std::pair<int, ValueType> a, std::pair<int, ValueType> b) {
              return a.first < b.first;
            });

  // transfer data to other data structure
  boundaryConditionNonGhostDofLocalNos_.reserve(
      boundaryConditionsNonGhost.size());
  boundaryConditionValues_.reserve(boundaryConditionsNonGhost.size());

  for (typename std::vector<std::pair<int, ValueType>>::iterator iter =
           boundaryConditionsNonGhost.begin();
       iter != boundaryConditionsNonGhost.end(); iter++) {
    boundaryConditionNonGhostDofLocalNos_.push_back(iter->first);
    boundaryConditionValues_.push_back(iter->second);
  }

  // create the boundaryConditionsByComponent_ data structure from
  // boundaryConditionNonGhostDofLocalNos_ and boundaryConditionValues_
  generateBoundaryConditionsByComponent();

  LOG(DEBUG) << "boundaryConditionsNonGhost: " << boundaryConditionsNonGhost
             << ", boundaryConditionNonGhostDofLocalNos: "
             << boundaryConditionNonGhostDofLocalNos_
             << ", boundaryConditionValues: " << boundaryConditionValues_
             << ", boundaryConditionsByComponent_: ";
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    LOG(DEBUG) << "  component " << componentNo << ", dofs "
               << boundaryConditionsByComponent_[componentNo].dofNosLocal
               << ", values "
               << boundaryConditionsByComponent_[componentNo].values;
  }
}

template <typename FunctionSpaceType, int nComponents>
void DirichletBoundaryConditionsBase<
    FunctionSpaceType, nComponents>::generateBoundaryConditionsByComponent() {
  // fill the data structure
  //  struct BoundaryConditionsForComponent
  //  {
  //    std::vector<dof_no_t> dofsLocal;
  //    std::vector<double> values;
  //  };
  //  std::array<BoundaryConditionsForComponent, nComponents>
  //  boundaryConditionsByComponent_;   //< the boundary condition data
  //  organized by component

  // from the two data structures
  //   std::vector<dof_no_t> boundaryConditionNonGhostDofLocalNos_;        //<
  //   vector of all local (non-ghost) boundary condition dofs
  //   std::vector<ValueType> boundaryConditionValues_;               //< vector
  //   of the local prescribed values, related to
  //   boundaryConditionNonGhostDofLocalNos_

  // The boundaryConditionValues_ vector can contain NaN values, where in the
  // config None was set. This indicates that these dofs should not be set as
  // Dirichlet boundary conditions.

  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    // preallocate helper data structure, as if all values were valid boundary
    // conditions
    std::vector<std::pair<int, double>> dofNosValues;
    dofNosValues.reserve(boundaryConditionValues_.size());

    // iterate over all values
    for (int i = 0; i < boundaryConditionValues_.size(); i++) {
      // if value is not std::numeric_limits<double>::max() and therefore valid
      // if (std::isfinite(boundaryConditionValues_[i][componentNo]) && )
      if (boundaryConditionValues_[i][componentNo] !=
          std::numeric_limits<double>::max()) {
        dofNosValues.push_back(
            std::pair<int, double>(boundaryConditionNonGhostDofLocalNos_[i],
                                   boundaryConditionValues_[i][componentNo]));
      }
    }

    // sort list according to dof no.s
    std::sort(dofNosValues.begin(), dofNosValues.end(),
              [&](std::pair<int, double> a, std::pair<int, double> b) {
                return a.first < b.first;
              });

    // assign values in boundaryConditionsByComponent_
    boundaryConditionsByComponent_[componentNo].dofNosLocal.resize(
        dofNosValues.size());
    boundaryConditionsByComponent_[componentNo].values.resize(
        dofNosValues.size());

    for (int i = 0; i < dofNosValues.size(); i++) {
      boundaryConditionsByComponent_[componentNo].dofNosLocal[i] =
          dofNosValues[i].first;
      boundaryConditionsByComponent_[componentNo].values[i] =
          dofNosValues[i].second;
    }
  }
}

//! add boundary conditions to the currently present boundary conditions
template <typename FunctionSpaceType, int nComponents>
void DirichletBoundaryConditionsBase<FunctionSpaceType, nComponents>::
    addBoundaryConditions(
        std::vector<ElementWithNodes> &newBoundaryConditionElements,
        bool overwriteBcOnSameDof) {
  // merge newBoundaryConditionElements into boundaryConditionElements_
#ifndef NDEBUG
  LOG(DEBUG) << "addBoundaryConditions, overwriteBcOnSameDof="
             << overwriteBcOnSameDof << ", "
             << newBoundaryConditionElements.size() << " elements";
  for (ElementWithNodes newElementWithNodes : newBoundaryConditionElements) {
    std::stringstream s;
    for (std::pair<int, ValueType> dofIndex :
         newElementWithNodes.elementalDofIndex) {
      s << ", " << dofIndex.first << ": (";
      for (int i = 0; i < dofIndex.second.size(); i++) {
        if (i != 0)
          s << ",";
        s << dofIndex.second[i];
      }
      s << ")";
    }
    LOG(DEBUG) << "newElementWithNodes: elementNoLocal "
               << newElementWithNodes.elementNoLocal
               << ", dofVectors: " << s.str();
  }
#endif

  // loop over given boundaryConditionElements
  for (const ElementWithNodes &newElementWithNodes :
       newBoundaryConditionElements) {
    // find corresponding element in boundaryConditionElements_
    bool elementFound = false;
    for (ElementWithNodes &existingElementWithNodes :
         boundaryConditionElements_) {
      // if the element was found
      if (existingElementWithNodes.elementNoLocal ==
          newElementWithNodes.elementNoLocal) {
        LOG(DEBUG) << "element " << newElementWithNodes.elementNoLocal
                   << " already exists with elementalDofIndex: "
                   << newElementWithNodes.elementalDofIndex;

        // loop over elemental dof indices and add all entries
        // std::vector<std::pair<int,ValueType>> elementalDofIndex;
        for (const std::pair<const int, ValueType> &newElementalDofIndex :
             newElementWithNodes.elementalDofIndex) {
          bool elementalDofIndexFound = false;

          // if the dof already has a prescribed value
          if (existingElementWithNodes.elementalDofIndex.find(
                  newElementalDofIndex.first) !=
              existingElementWithNodes.elementalDofIndex.end()) {
            // if the bc values should be overwritten
            if (overwriteBcOnSameDof) {
              LOG(DEBUG) << "  dof " << newElementalDofIndex.first
                         << ", overwrite previous bc value "
                         << existingElementWithNodes
                                .elementalDofIndex[newElementalDofIndex.first]
                         << " with new value " << newElementalDofIndex.second;
              existingElementWithNodes
                  .elementalDofIndex[newElementalDofIndex.first] =
                  newElementalDofIndex.second;
            }

            elementalDofIndexFound = true;
          }

          if (!elementalDofIndexFound) {
            LOG(DEBUG) << "  add dof " << newElementalDofIndex;

            // if the newElementalDofIndex is not already present, add the whole
            // data structure
            existingElementWithNodes.elementalDofIndex.insert(
                newElementalDofIndex);
          }
        }

        elementFound = true;
        break;
      }
    }

    // if element was not existent yet, add it
    if (!elementFound) {
      LOG(DEBUG) << "add element " << newElementWithNodes.elementNoLocal;

      // add element
      boundaryConditionElements_.push_back(newElementWithNodes);
    }
  }

  // clear previous boundary condition values and create anew
  boundaryConditionNonGhostDofLocalNos_.clear();
  boundaryConditionValues_.clear();

  std::set<dof_no_t>
      boundaryConditionNonGhostDofLocalNosSet; //< same data as in
                                               // boundaryConditionNonGhostDofLocalNos_,
                                               // but as set
  std::vector<std::pair<int, ValueType>>
      boundaryConditionsNonGhost; //< boundary condition dof nos and values for
                                  // non-ghost dofs, used to collect values and
                                  // sort them afterwards

  // create boundaryConditionNonGhostDofLocalNos_ and boundaryConditionValues_
  // from boundaryConditionValues_ loop over existing boundaryConditionElements_
  for (const ElementWithNodes &elementWithNodes : boundaryConditionElements_) {
    int elementNoLocal = elementWithNodes.elementNoLocal;

    // loop over elemental dof indices
    for (const std::pair<const int, ValueType> &elementalDofIndexAndValues :
         elementWithNodes.elementalDofIndex) {
      int elementalDofIndex = elementalDofIndexAndValues.first;
      ValueType boundaryConditionValue = elementalDofIndexAndValues.second;

      // also store localDofNo
      dof_no_t dofLocalNo =
          functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);

      // if the dof is non-ghost
      if (dofLocalNo < functionSpace_->nDofsLocalWithoutGhosts()) {
        // if dofLocalNo is not already contained in
        // boundaryConditionGhostDofLocalNos_
        if (boundaryConditionNonGhostDofLocalNosSet.find(dofLocalNo) ==
            boundaryConditionNonGhostDofLocalNosSet.end()) {
          boundaryConditionNonGhostDofLocalNosSet.insert(dofLocalNo);
          boundaryConditionsNonGhost.push_back(
              std::pair<int, ValueType>(dofLocalNo, boundaryConditionValue));
        }
      }
    }
  }

  // sort list according to dof no.s
  std::sort(boundaryConditionsNonGhost.begin(),
            boundaryConditionsNonGhost.end(),
            [&](std::pair<int, ValueType> a, std::pair<int, ValueType> b) {
              return a.first < b.first;
            });

  // transfer data to other data structure
  boundaryConditionNonGhostDofLocalNos_.reserve(
      boundaryConditionsNonGhost.size());
  boundaryConditionValues_.reserve(boundaryConditionsNonGhost.size());

  for (typename std::vector<std::pair<int, ValueType>>::iterator iter =
           boundaryConditionsNonGhost.begin();
       iter != boundaryConditionsNonGhost.end(); iter++) {
    boundaryConditionNonGhostDofLocalNos_.push_back(iter->first);
    boundaryConditionValues_.push_back(iter->second);
  }

  // create the boundaryConditionsByComponent_ data structure from
  // boundaryConditionNonGhostDofLocalNos_ and boundaryConditionValues_
  generateBoundaryConditionsByComponent();

  // fill auxiliary ghost element data structures, this is only needed for
  // Dirichlet boundary conditions on scalar fields
  this->initializeGhostElements();

  LOG(DEBUG) << "boundaryConditionsNonGhost: " << boundaryConditionsNonGhost
             << ", boundaryConditionNonGhostDofLocalNos: "
             << boundaryConditionNonGhostDofLocalNos_
             << ", boundaryConditionValues: " << boundaryConditionValues_
             << ", boundaryConditionsByComponent_: ";
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    LOG(DEBUG) << "  component " << componentNo << ", dofs "
               << boundaryConditionsByComponent_[componentNo].dofNosLocal
               << ", values "
               << boundaryConditionsByComponent_[componentNo].values;
  }
}

//! write the constraint dofs and the prescribed values as points in a vtk file
template <typename FunctionSpaceType, int nComponents>
void DirichletBoundaryConditionsBase<FunctionSpaceType,
                                     nComponents>::writeOutput() {
  // get the filename
  filenameOutput_ =
      specificSettings_.getOptionString("dirichletOutputFilename", "");

  // if filename option was None or empty string
  if (filenameOutput_ != "") {
    // get geometry of boundary condition dofs
    std::vector<Vec3> geometryValues;
    this->functionSpace_->geometryField().getValues(
        boundaryConditionNonGhostDofLocalNos_, geometryValues);

    // output the dofs and their values
    // outputPointsVtp(std::string filename, double currentTime, const
    // std::vector<Vec3> &geometry,
    //           const std::vector<VecD<nComponents>> &vectorValues, std::string
    //           fieldVariableName, std::shared_ptr<Partition::RankSubset>
    //           rankSubset)
    outputPoints_.outputPointsVtp<nComponents>(
        filenameOutput_, 0, geometryValues, boundaryConditionValues_,
        "dirichletBoundaryConditions",
        functionSpace_->meshPartition()->rankSubset());
  }
}

//! get a reference to the vector of bc local dofs
template <typename FunctionSpaceType, int nComponents>
const std::vector<dof_no_t> &DirichletBoundaryConditionsBase<
    FunctionSpaceType, nComponents>::boundaryConditionNonGhostDofLocalNos()
    const {
  return boundaryConditionNonGhostDofLocalNos_;
}

//! get a reference to the vector of bc local dofs
template <typename FunctionSpaceType, int nComponents>
const std::vector<typename DirichletBoundaryConditionsBase<
    FunctionSpaceType, nComponents>::ValueType> &
DirichletBoundaryConditionsBase<FunctionSpaceType,
                                nComponents>::boundaryConditionValues() const {
  return boundaryConditionValues_;
}

//! get the boundary conditions data organized by component
template <typename FunctionSpaceType, int nComponents>
const std::array<
    typename DirichletBoundaryConditionsBase<
        FunctionSpaceType, nComponents>::BoundaryConditionsForComponent,
    nComponents> &
DirichletBoundaryConditionsBase<
    FunctionSpaceType, nComponents>::boundaryConditionsByComponent() const {
  return boundaryConditionsByComponent_;
}

template <typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(
    std::ostream &stream,
    const typename DirichletBoundaryConditionsBase<
        FunctionSpaceType, nComponents>::BoundaryConditionsForComponent rhs) {
  stream << "{dofNos: " << rhs.dofNosLocal << ", values: " << rhs.values << "}";
  return stream;
}

template <typename FunctionSpaceType, int nComponents>
std::ostream &
operator<<(std::ostream &stream,
           const typename DirichletBoundaryConditionsBase<
               FunctionSpaceType, nComponents>::ElementWithNodes rhs) {
  stream << "{el." << rhs.elementNoLocal
         << ", (dof,v):" << rhs.elementalDofIndex << "}";
  return stream;
}

} // namespace SpatialDiscretization
