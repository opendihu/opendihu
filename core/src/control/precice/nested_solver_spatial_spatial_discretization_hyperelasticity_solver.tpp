#include "control/precice/nested_solver.h"

namespace Control {
template <typename Material>
std::shared_ptr<typename SpatialDiscretization::HyperelasticitySolver<
    Material>::FunctionSpace>
PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::functionSpace(NestedSolverType &nestedSolver) {
  return nestedSolver.data().functionSpace();
}

template <typename Material>
template <typename PreciceVolumeData, typename PreciceVolumeMesh>
void PreciceAdapterNestedSolver<
    SpatialDiscretization::HyperelasticitySolver<Material>>::
    initializePreciceVolumeData(
        PythonConfig &specificSettings, NestedSolverType &nestedSolver,
        std::shared_ptr<precice::Participant> &preciceParticipant,
        std::vector<PreciceVolumeData> &preciceVolumeData,
        std::vector<std::shared_ptr<PreciceVolumeMesh>> &preciceVolumeMeshes) {}

template <typename Material>
template <typename SurfaceDataVector, typename VolumeDataVector>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::preciceReadData(NestedSolverType &nestedSolver,
                                std::shared_ptr<precice::Participant>
                                    &preciceParticipant,
                                SurfaceDataVector &preciceSurfaceData,
                                VolumeDataVector &preciceVolumeData,
                                DihuContext &context) {
  LOG(DEBUG) << "read surface data from precice";
  double preciceDt = preciceParticipant->getMaxTimeStepSize();
  // loop over data
  for (auto &preciceData : preciceSurfaceData) {
    if (preciceData.ioType ==
        PreciceAdapterInitialize<
            NestedSolverType>::PreciceSurfaceData::ioRead) {
      // allocate memory
      int nEntries = preciceData.preciceMesh->nNodesLocal * 3;

      // if the data is displacements and velocities
      if (!preciceData.displacementsName.empty()) {
        displacementValues_.resize(nEntries);
        velocityValues_.resize(nEntries);

        // get all data at once
        preciceParticipant->readData(preciceData.preciceMesh->preciceMeshName,
                                     preciceData.displacementsName,
                                     preciceData.preciceMesh->preciceVertexIds,
                                     preciceDt, displacementValues_);

        preciceParticipant->readData(preciceData.preciceMesh->preciceMeshName,
                                     preciceData.velocitiesName,
                                     preciceData.preciceMesh->preciceVertexIds,
                                     preciceDt, velocityValues_);

        setDirichletBoundaryConditions(preciceData, nestedSolver);
      }
      // if the data is traction
      else if (!preciceData.tractionName.empty()) {
        tractionValues_.resize(nEntries);
        preciceParticipant->readData(preciceData.preciceMesh->preciceMeshName,
                                     preciceData.tractionName,
                                     preciceData.preciceMesh->preciceVertexIds,
                                     preciceDt, tractionValues_);

        setNeumannBoundaryConditions(preciceData, nestedSolver, context);
      } else {
        LOG(FATAL) << "Unknown precice data (read), none of displacements, "
                      "velocities or traction is set.";
      }
    }
  };
}

template <typename Material>
template <typename SurfaceDataVector, typename VolumeDataVector>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::preciceWriteData(NestedSolverType &nestedSolver,
                                 std::shared_ptr<precice::Participant>
                                     &preciceParticipant,
                                 SurfaceDataVector &preciceSurfaceData,
                                 VolumeDataVector &preciceVolumeData,
                                 double scalingFactor) {
  // write data to precice
  LOG(DEBUG) << "write surface data to precice";

  // loop over data
  for (auto &preciceData : preciceSurfaceData) {
    if (preciceData.ioType ==
        PreciceAdapterInitialize<
            NestedSolverType>::PreciceSurfaceData::ioWrite) {
      // if the data is displacements and velocities
      if (!preciceData.displacementsName.empty()) {
        // convert geometry values to precice data layout
        displacementValues_.clear();
        velocityValues_.clear();

        this->getDisplacementVelocityValues(
            nestedSolver, preciceData.preciceMesh->dofNosLocal,
            displacementValues_, velocityValues_);

#ifndef NDEBUG
        LOG(DEBUG) << "write displacements data to precice: "
                   << displacementValues_;
#endif
        // scale displacement and velocity values
        for (double &value : displacementValues_)
          value *= scalingFactor;

        for (double &value : velocityValues_)
          value *= scalingFactor;

        // write displacement values in precice
        preciceParticipant->writeData(preciceData.preciceMesh->preciceMeshName,
                                      preciceData.displacementsName,
                                      preciceData.preciceMesh->preciceVertexIds,
                                      displacementValues_);

        // write velocity values in precice
        preciceParticipant->writeData(preciceData.preciceMesh->preciceMeshName,
                                      preciceData.velocitiesName,
                                      preciceData.preciceMesh->preciceVertexIds,
                                      velocityValues_);
      }
      // if the data is traction
      else if (!preciceData.tractionName.empty() &&
               preciceData.average == true) {
        LOG(INFO) << "Write averaged traction";
        // convert geometry values to precice data layout
        tractionValues_.clear();
        this->getTractionValues(nestedSolver,
                                preciceData.preciceMesh->dofNosLocal,
                                tractionValues_);
        // average z-values of traction
        double average_traction = 0.0;
        for (int i = 2; i < tractionValues_.size(); i += 3) {
          average_traction += tractionValues_[i];
        }
        average_traction /= (tractionValues_.size() / 3);
        for (int i = 2; i < tractionValues_.size(); i += 3) {
          tractionValues_[i] = average_traction;
        }
#ifndef NDEBUG
        LOG(DEBUG) << "Write averaged-traction data to precice: "
                   << tractionValues_[2];
#endif
        // scale traction values, they are always scaled by the factor of -1
        for (double &value : tractionValues_) {
          value *= scalingFactor;
          value *= -1;
        }

        preciceParticipant->writeData(
            preciceData.preciceMesh->preciceMeshName, preciceData.tractionName,
            preciceData.preciceMesh->preciceVertexIds, tractionValues_);
      } else if (!preciceData.tractionName.empty()) {

        LOG(INFO) << "Write non-averaged traction";

        // convert geometry values to precice data layout
        tractionValues_.clear();
        this->getTractionValues(nestedSolver,
                                preciceData.preciceMesh->dofNosLocal,
                                tractionValues_);

#ifndef NDEBUG
        LOG(DEBUG) << "write traction data to precice: " << tractionValues_;
        std::stringstream s;
        for (int i = 2; i < tractionValues_.size(); i += 3) {
          s << " " << tractionValues_[i];
        }
        LOG(DEBUG) << "z values of traction: " << s.str();
#endif
        // scale traction values, they are always scaled by the factor of -1
        for (double &value : tractionValues_) {
          value *= scalingFactor;
          value *= -1;
        }

        preciceParticipant->writeData(
            preciceData.preciceMesh->preciceMeshName, preciceData.tractionName,
            preciceData.preciceMesh->preciceVertexIds, tractionValues_);
      } else {
        LOG(FATAL) << "Unknown precice data (write), none of displacements, "
                      "velocities or traction is set.";
      }
    }
  }

  LOG(DEBUG) << "write traction data to precice complete";
}

template <typename Material>
template <typename SurfaceData>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::setDirichletBoundaryConditions(SurfaceData &preciceData,
                                               NestedSolverType &nestedSolver) {
  std::vector<std::pair<global_no_t, std::array<double, 6>>>
      newDirichletBCValues;

  auto functionSpace = this->functionSpace(nestedSolver);

  // if this rank has no data, do not set any boundary conditions
  if (preciceData.preciceMesh->nNodesLocal != 0) {
    // loop over nodes
    const int nNodesX = functionSpace->nNodesLocalWithoutGhosts(0);
    const int nNodesY = functionSpace->nNodesLocalWithoutGhosts(1);
    const int nNodesZ = functionSpace->nNodesLocalWithoutGhosts(2);

    // set node index in z direction for bottom surface
    int nodeIndexZ = 0;

    // for top surface
    if (preciceData.preciceMesh->face ==
        PreciceAdapterInitialize<
            NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
      nodeIndexZ = nNodesZ - 1;
    }

    // loop over nodes to set the received values
    newDirichletBCValues.reserve(nNodesX * nNodesY);
    int valueIndex = 0;

    // loop over nodes of surface mesh
    for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++) {
      for (int nodeIndexX = 0; nodeIndexX < nNodesX;
           nodeIndexX++, valueIndex++) {
        node_no_t nodeNoLocal =
            nodeIndexZ * nNodesX * nNodesY + nodeIndexY * nNodesX + nodeIndexX;

        dof_no_t dofNoLocal = nodeNoLocal;
        global_no_t dofNoGlobal =
            functionSpace->meshPartition()->getDofNoGlobalPetsc(dofNoLocal);

        // assign received values to dirichlet bc vector of size 6
        std::array<double, 6> newDirichletBCValue;

        for (int i = 0; i < 3; i++) {
          newDirichletBCValue[i] = displacementValues_[3 * valueIndex + i];
          newDirichletBCValue[3 + i] = velocityValues_[3 * valueIndex + i];
        }

        newDirichletBCValues.push_back(
            std::pair<global_no_t, std::array<double, 6>>(dofNoGlobal,
                                                          newDirichletBCValue));
      }
    }

    LOG(DEBUG) << "read data from precice complete, displacement values: "
               << displacementValues_
               << ", velocityValues: " << velocityValues_;
    LOG(DEBUG) << "read and set Dirichlet BC: " << newDirichletBCValues;
  }

  //! set new dirichlet boundary condition values
  this->updateDirichletBoundaryConditions(nestedSolver, newDirichletBCValues);
}

template <typename Material>
template <typename SurfaceData>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::setNeumannBoundaryConditions(SurfaceData &preciceData,
                                             NestedSolverType &nestedSolver,
                                             DihuContext &context) {
  // set traction values as neumann boundary conditions
  using FunctionSpace =
      typename PreciceAdapterNestedSolver<NestedSolverType>::FunctionSpace;
  using ElementWithFacesType =
      typename SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>::ElementWithFaces;

  std::vector<ElementWithFacesType> neumannBoundaryConditionElements;

  auto functionSpace = this->functionSpace(nestedSolver);

  // if this rank has no data, do not set any boundary conditions
  if (preciceData.preciceMesh->nNodesLocal != 0) {

    const int nNodesX = functionSpace->nNodesLocalWithoutGhosts(0);
    const int nNodesY = functionSpace->nNodesLocalWithoutGhosts(1);
    const int nNodesZ = functionSpace->nNodesLocalWithoutGhosts(2);

    int nElementsX = functionSpace->meshPartition()->nElementsLocal(0);
    int nElementsY = functionSpace->meshPartition()->nElementsLocal(1);
    int nElementsZ = functionSpace->meshPartition()->nElementsLocal(2);

    // set node and element indics in z direction for bottom surface
    int nodeIndexZ = 0;
    int elementalNodeIndexZ = 0;
    int elementIndexZ = 0;

    // for top surface
    if (preciceData.preciceMesh->face ==
        PreciceAdapterInitialize<
            NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
      nodeIndexZ = nNodesZ - 1;
      elementIndexZ = nElementsZ - 1;
      elementalNodeIndexZ = 2;
    }

    // loop over elements
    for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++) {
      for (int elementIndexX = 0; elementIndexX < nElementsX; elementIndexX++) {
        ElementWithFacesType elementWithFaces;
        element_no_t elementNoLocal = elementIndexZ * nElementsX * nElementsY +
                                      elementIndexY * nElementsX +
                                      elementIndexX;
        elementWithFaces.elementNoLocal = elementNoLocal;

        // set surface dofs
        Mesh::face_t face = Mesh::face_t::face2Minus;
        if (preciceData.preciceMesh->face ==
            PreciceAdapterInitialize<
                NestedSolverType>::PreciceSurfaceMesh::face2Plus)
          face = Mesh::face_t::face2Plus;

        elementWithFaces.face = face;

        // get dofs indices within the numbering of the volume element that
        // correspond to the selected face
        const int nDofsPerNode = FunctionSpace::nDofsPerNode();
        const int nSurfaceDofs =
            ::FunctionSpace::FunctionSpaceBaseDim<
                2, typename FunctionSpace::BasisFunction>::nNodesPerElement() *
            nDofsPerNode;
        std::array<dof_no_t, nSurfaceDofs> surfaceDofs;
        FunctionSpace::getFaceDofs(face, surfaceDofs);

        elementWithFaces.surfaceDofs.assign(surfaceDofs.begin(),
                                            surfaceDofs.end());

        // loop over the nodes of the element
        for (int elementalNodeIndexY = 0; elementalNodeIndexY < 3;
             elementalNodeIndexY++) {
          for (int elementalNodeIndexX = 0; elementalNodeIndexX < 3;
               elementalNodeIndexX++) {
            int elementalDofIndex = elementalNodeIndexZ * 9 +
                                    elementalNodeIndexY * 3 +
                                    elementalNodeIndexX;

            dof_no_t dofNoLocal =
                functionSpace->getDofNo(elementNoLocal, elementalDofIndex);

            // only iterate over local dofs, the ghost dofs are not considered
            // here (although it would be correct to communicate them with
            // precice)
            if (dofNoLocal >= functionSpace->nDofsLocalWithoutGhosts())
              continue;

            int valueIndex = dofNoLocal - nodeIndexZ * nNodesX * nNodesY;

            LOG(DEBUG) << "(x,y,z)=(" << elementalNodeIndexX << ","
                       << elementalNodeIndexY << "," << elementalNodeIndexZ
                       << ") dofNoLocal " << dofNoLocal
                       << ", valueIndex: " << valueIndex << "/"
                       << tractionValues_.size() / 3;

            assert(valueIndex >= 0);
            if (valueIndex >= preciceData.preciceMesh->nNodesLocal) {
              LOG(ERROR) << "valueIndex: " << valueIndex
                         << ", dofNoLocal: " << dofNoLocal
                         << ", nNodes: " << nNodesX << "," << nNodesY
                         << ", nodeIndexZ: " << nodeIndexZ << ", nNodesLocal: "
                         << preciceData.preciceMesh->nNodesLocal
                         << " elementNoLocal: " << elementNoLocal
                         << ", elementalDofIndex: " << elementalDofIndex
                         << ", elementalNodeIndexZ: " << elementalNodeIndexZ
                         << ", meshPartition: "
                         << *functionSpace->meshPartition();
            }
            assert(valueIndex < preciceData.preciceMesh->nNodesLocal);

            Vec3 traction;
            for (int i = 0; i < 3; i++) {
              traction[i] = tractionValues_[3 * valueIndex + i];
            }

            dof_no_t surfaceDof = elementalDofIndex;

            // for top surface
            if (preciceData.preciceMesh->face ==
                PreciceAdapterInitialize<
                    NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
              surfaceDof = 18 + elementalDofIndex;
            }
            elementWithFaces.dofVectors.push_back(
                std::pair<dof_no_t, Vec3>(surfaceDof, traction));

            // LOG(INFO) << "dofVectors: " << elementWithFaces.dofVectors <<
            // ", traction: " << traction;
          }
        }
        neumannBoundaryConditionElements.push_back(elementWithFaces);
      }
    }
  }
  // create new Neumann BC object
  using NeumannBoundaryConditionsType =
      SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,
                                                       Quadrature::Gauss<3>, 3>;
  std::shared_ptr<NeumannBoundaryConditionsType> neumannBoundaryConditions =
      std::make_shared<NeumannBoundaryConditionsType>(context);
  neumannBoundaryConditions->initialize(functionSpace,
                                        neumannBoundaryConditionElements);
  neumannBoundaryConditions->setDeformationGradientField(
      this->deformationGradientField(nestedSolver));

#ifndef NDEBUG
  std::stringstream s;
  for (int i = 0; i < neumannBoundaryConditionElements.size(); i++) {
    s << "{el. " << neumannBoundaryConditionElements[i].elementNoLocal << ", \""
      << Mesh::getString(neumannBoundaryConditionElements[i].face)
      << "\", dofVectors: [";
    for (int j = 0; j < neumannBoundaryConditionElements[i].dofVectors.size();
         j++) {
      s << "(" << neumannBoundaryConditionElements[i].dofVectors[j].first
        << ": (";
      for (int k = 0;
           k < neumannBoundaryConditionElements[i].dofVectors[j].second.size();
           k++) {
        if (k != 0)
          s << ",";
        s << neumannBoundaryConditionElements[i].dofVectors[j].second[k];
      }
      s << ")),";
    }
    s << "], surfaceDofs: ";
    for (int j = 0; j < neumannBoundaryConditionElements[i].surfaceDofs.size();
         j++)
      s << neumannBoundaryConditionElements[i].surfaceDofs[j] << ",";
    s << "} ";
  }
  LOG(DEBUG) << "read and set Neumann BC:\n" << s.str();
#endif

  // set Neumann BCs in the static hyperelasticity of the
  // TimeSteppingScheme::DynamicHyperelasticitySolver solver
  this->updateNeumannBoundaryConditions(nestedSolver,
                                        neumannBoundaryConditions);

  LOG(DEBUG) << "read data from precice complete, traction values: "
             << tractionValues_;
}

template <typename Material>
void PreciceAdapterNestedSolver<
    SpatialDiscretization::HyperelasticitySolver<Material>>::
    addDirichletBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::vector<
            typename SpatialDiscretization::DirichletBoundaryConditionsBase<
                FunctionSpace, 6>::ElementWithNodes>
            &dirichletBoundaryConditionElements) {
  // transform the ElementWithNodes variable from one with 6 components
  // (displacements + velocities) to one with only 3 components
  std::vector<typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 3>::ElementWithNodes>
      dirichletBoundaryConditionElements3;
  for (typename SpatialDiscretization::DirichletBoundaryConditionsBase<
           FunctionSpace, 6>::ElementWithNodes &element6 :
       dirichletBoundaryConditionElements) {
    typename SpatialDiscretization::DirichletBoundaryConditionsBase<
        FunctionSpace, 3>::ElementWithNodes element3;

    // copy elementNo from element6 to element3
    element3.elementNoLocal = element6.elementNoLocal;

    // copy all bc entries from element6 to element3, use the first 3 components
    // for every prescribed value
    for (std::map<int, std::array<double, 6>>::const_iterator iter =
             element6.elementalDofIndex.cbegin();
         iter != element6.elementalDofIndex.cend(); iter++) {
      int elementalDofIndex = iter->first;
      std::array<double, 3> value = {iter->second[0], iter->second[1],
                                     iter->second[2]};
      element3.elementalDofIndex.insert(
          std::pair<int, std::array<double, 3>>(elementalDofIndex, value));
    }

    dirichletBoundaryConditionElements3.push_back(element3);
  }

  // add the dirichlet bc values
  bool overwriteBcOnSameDof = true;
  nestedSolver.addDirichletBoundaryConditions(
      dirichletBoundaryConditionElements3, overwriteBcOnSameDof);
}

//! update existing boundary conditions with new values
template <typename Material>
void PreciceAdapterNestedSolver<
    SpatialDiscretization::HyperelasticitySolver<Material>>::
    updateDirichletBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::vector<std::pair<global_no_t, std::array<double, 6>>>
            newDirichletBoundaryConditionValues) {
  // transform the ElementWithNodes variable from one with 6 components
  // (displacements + velocities) to one with only 3 components
  std::vector<std::pair<global_no_t, std::array<double, 3>>>
      newDirichletBoundaryConditionValues3;
  for (std::pair<global_no_t, std::array<double, 6>> &element6 :
       newDirichletBoundaryConditionValues) {
    std::pair<global_no_t, std::array<double, 3>> element3;

    element3.first = element6.first;

    element3.second[0] = element6.second[0];
    element3.second[1] = element6.second[1];
    element3.second[2] = element6.second[2];

    newDirichletBoundaryConditionValues3.push_back(element3);
  }

  nestedSolver.updateDirichletBoundaryConditions(
      newDirichletBoundaryConditionValues3);
}

template <typename Material>
void PreciceAdapterNestedSolver<
    SpatialDiscretization::HyperelasticitySolver<Material>>::
    updateNeumannBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>>
            neumannBoundaryConditions) {
  nestedSolver.updateNeumannBoundaryConditions(neumannBoundaryConditions);
}

//! get the displacement and velocity vectors of the given local dof nos
template <typename Material>
void PreciceAdapterNestedSolver<
    SpatialDiscretization::HyperelasticitySolver<Material>>::
    getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                  const std::vector<dof_no_t> &dofNosLocal,
                                  std::vector<double> &displacementValues,
                                  std::vector<double> &velocityValues) {
  // get the displacement values
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.data().displacements()->getValues(dofNosLocal, values);

  // store displacement values in interleaved order (array of struct)
  int nVectors = values.size();
  displacementValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++) {
    displacementValues[3 * i + 0] = values[i][0];
    displacementValues[3 * i + 1] = values[i][1];
    displacementValues[3 * i + 2] = values[i][2];
  }

  // get the velocity values
  nVectors = values.size();
  velocityValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++) {
    velocityValues[3 * i + 0] = values[i][0];
    velocityValues[3 * i + 1] = values[i][1];
    velocityValues[3 * i + 2] = values[i][2];
  }
}

//! get the traction vectors of the given local dof nos
template <typename Material>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::getTractionValues(NestedSolverType &nestedSolver,
                                  const std::vector<dof_no_t> &dofNosLocal,
                                  std::vector<double> &tractionValues) {
  /*std::vector<Vec3> values0;
  nestedSolver.data().materialTraction()->getValuesWithoutGhosts(values0);

  std::stringstream s;
  for (int i = 0; i < values0.size(); i++)
  {
    s << " " << values0[i][2];
  }
  LOG(INFO) << values0.size() << " local traction values in total (), "
    << "\ndofNosLocal: " << dofNosLocal << "\nall values: " << values0 << "\nz
  values: " << s.str();
*/
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.data().materialTraction()->getValues(dofNosLocal, values);

  int nVectors = values.size();
  tractionValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++) {
    tractionValues[3 * i + 0] = values[i][0];
    tractionValues[3 * i + 1] = values[i][1];
    tractionValues[3 * i + 2] = values[i][2];
  }
}

template <typename Material>
Vec PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::currentState(NestedSolverType &nestedSolver) {
  return nestedSolver.currentState();
}

template <typename Material>
std::shared_ptr<FieldVariable::FieldVariable<
    typename PreciceAdapterNestedSolver<
        SpatialDiscretization::HyperelasticitySolver<Material>>::FunctionSpace,
    9>>
PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::deformationGradientField(NestedSolverType &nestedSolver) {
  return nestedSolver.data().deformationGradient();
}

template <typename Material>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::reset(NestedSolverType &nestedSolver) {
  nestedSolver.reset();
}

template <typename Material>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::saveFiberData(NestedSolverType &nestedSolver) {}

template <typename Material>
void PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::loadFiberData(NestedSolverType &nestedSolver) {}
} // namespace Control
