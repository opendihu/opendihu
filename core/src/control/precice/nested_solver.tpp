#include "control/precice/nested_solver.h"
#include "control/precice/initialize.h"

#include <sstream>

namespace Control {


template <typename NestedSolver>
void PreciceAdapterNestedSolver<NestedSolver>::preciceReadData_impl(){
    LOG(DEBUG) << "read data from precice";
    double preciceDt = this->preciceParticipant_->getMaxTimeStepSize();
    // loop over surface data
    for (auto &preciceData : this->preciceSurfaceData_) {
    // bool tmp = preciceData.ioType == ;
        // preciceData.ioType = "tmp"; 
        // if (preciceData.ioType == decltype(preciceData.ioType)::ioRead) {
        // using typename Control::PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData::ReadWrite;
        if (preciceData.ioType == PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData::ioRead) {

        // allocate memory
        int nEntries = preciceData.preciceMesh->nNodesLocal * 3;

        // if the data is displacements and velocities
        if (!preciceData.displacementsName.empty()) {
            displacementValues_.resize(nEntries);
            velocityValues_.resize(nEntries);

            // get all data at once
            this->preciceParticipant_->readData(
                preciceData.preciceMesh->preciceMeshName,
                preciceData.displacementsName,
                preciceData.preciceMesh->preciceVertexIds, preciceDt,
                displacementValues_);

            this->preciceParticipant_->readData(
                preciceData.preciceMesh->preciceMeshName,
                preciceData.velocitiesName,
                preciceData.preciceMesh->preciceVertexIds, preciceDt,
                velocityValues_);

            setDirichletBoundaryConditions(preciceData);
        }
        // if the data is traction
        else if (!preciceData.tractionName.empty()) {
            tractionValues_.resize(nEntries);
            this->preciceParticipant_->readData(
                preciceData.preciceMesh->preciceMeshName, preciceData.tractionName,
                preciceData.preciceMesh->preciceVertexIds, preciceDt,
                tractionValues_);

            setNeumannBoundaryConditions(preciceData);
        } else {
            LOG(FATAL) << "Unknown precice data (read), none of displacements, "
                        "velocities or traction is set.";
        }
        }
    }
};

template <typename NestedSolver>
void PreciceAdapterNestedSolver<NestedSolver>::preciceWriteData_impl(){
  // write data to precice
    LOG(DEBUG) << "write data to precice";

    // loop over data
    for (auto &preciceData : this->preciceSurfaceData_) {
        if (preciceData.ioType ==
             PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData::ioWrite) {
        // if the data is displacements and velocities
        if (!preciceData.displacementsName.empty()) {
            // convert geometry values to precice data layout
            displacementValues_.clear();
            velocityValues_.clear();

            this->getDisplacementVelocityValues(
                this->nestedSolver_, preciceData.preciceMesh->dofNosLocal,
                displacementValues_, velocityValues_);

    #ifndef NDEBUG
            LOG(DEBUG) << "write displacements data to precice: "
                    << displacementValues_;
    #endif
            // scale displacement and velocity values
            for (double &value : displacementValues_)
            value *= this->scalingFactor_;

            for (double &value : velocityValues_)
            value *= this->scalingFactor_;

            // write displacement values in precice
            this->preciceParticipant_->writeData(
                preciceData.preciceMesh->preciceMeshName,
                preciceData.displacementsName,
                preciceData.preciceMesh->preciceVertexIds, displacementValues_);

            // write velocity values in precice
            this->preciceParticipant_->writeData(
                preciceData.preciceMesh->preciceMeshName,
                preciceData.velocitiesName,
                preciceData.preciceMesh->preciceVertexIds, velocityValues_);
        }
        // if the data is traction
        else if (!preciceData.tractionName.empty() &&
                preciceData.average == true) {
            LOG(INFO) << "Write averaged traction";
            // convert geometry values to precice data layout
            tractionValues_.clear();
            this->getTractionValues(this->nestedSolver_,
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
            value *= this->scalingFactor_;
            value *= -1;
            }

            this->preciceParticipant_->writeData(
                preciceData.preciceMesh->preciceMeshName, preciceData.tractionName,
                preciceData.preciceMesh->preciceVertexIds, tractionValues_);
        } else if (!preciceData.tractionName.empty()) {

            LOG(INFO) << "Write non-averaged traction";

            // convert geometry values to precice data layout
            tractionValues_.clear();
            this->getTractionValues(this->nestedSolver_,
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
            value *= this->scalingFactor_;
            value *= -1;
            }

            this->preciceParticipant_->writeData(
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

template <typename NestedSolver>
void PreciceAdapterNestedSolver<NestedSolver>::setNeumannBoundaryConditions_impl(typename PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData &preciceData){
   using FunctionSpace =
                typename PreciceAdapterNestedSolver<NestedSolver>::FunctionSpace;
            using ElementWithFacesType =
                typename SpatialDiscretization::NeumannBoundaryConditions<
                    FunctionSpace, Quadrature::Gauss<3>, 3>::ElementWithFaces;
            /*
            struct ElementWithFaces
            {
                element_no_t elementNoLocal;                                     //< the
            local no of the element

                Mesh::face_t face;                                               //< face on
            which the Neumann BC is applied std::vector<std::pair<dof_no_t,
            VecD<nComponents>>> dofVectors;  //< <surface-local dof no, value>,
            nComponents == FunctionSpaceType::dim() for traction boundary condition or
            nComponents = 1 for flux BC std::vector<dof_no_t> surfaceDofs; //< dof nos of
            the volume element that correspond to the face / surface. These are different
            from the dofs in dofsVector which are numbered for the surface only,
            surfaceDofs are in the numbering of the volume element.
                // note, for flux BC, dofVectors[i].second is a VecD<1>
            };*/

            std::vector<ElementWithFacesType> neumannBoundaryConditionElements;

            // if this rank has no data, do not set any boundary conditions
            if (preciceData.preciceMesh->nNodesLocal != 0) {

                const int nNodesX = this->functionSpace_->nNodesLocalWithoutGhosts(0);
                const int nNodesY = this->functionSpace_->nNodesLocalWithoutGhosts(1);
                const int nNodesZ = this->functionSpace_->nNodesLocalWithoutGhosts(2);

                int nElementsX = this->functionSpace_->meshPartition()->nElementsLocal(0);
                int nElementsY = this->functionSpace_->meshPartition()->nElementsLocal(1);
                int nElementsZ = this->functionSpace_->meshPartition()->nElementsLocal(2);

                // set node and element indics in z direction for bottom surface
                int nodeIndexZ = 0;
                int elementalNodeIndexZ = 0;
                int elementIndexZ = 0;

                // for top surface
                if (preciceData.preciceMesh->face == PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData::face2Plus) {
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
                         PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData::face2Plus)
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

                        dof_no_t dofNoLocal = this->functionSpace_->getDofNo(
                            elementNoLocal, elementalDofIndex);

                        // only iterate over local dofs, the ghost dofs are not considered
                        // here (although it would be correct to communicate them with
                        // precice)
                        if (dofNoLocal >= this->functionSpace_->nDofsLocalWithoutGhosts())
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
                                    << *this->functionSpace_->meshPartition();
                        }
                        assert(valueIndex < preciceData.preciceMesh->nNodesLocal);

                        Vec3 traction;
                        for (int i = 0; i < 3; i++) {
                        traction[i] = tractionValues_[3 * valueIndex + i];
                        }

                        dof_no_t surfaceDof = elementalDofIndex;

                        // for top surface
                        if (preciceData.preciceMesh->face ==
                             PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData::face2Plus) {
                        surfaceDof = 18 + elementalDofIndex;
                        }
                        elementWithFaces.dofVectors.push_back(
                            std::pair<dof_no_t, Vec3>(surfaceDof, traction));

                        // LOG(INFO) << "dofVectors: " << elementWithFaces.dofVectors << ",
                        // traction: " << traction;
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
                std::make_shared<NeumannBoundaryConditionsType>(this->context_);
            neumannBoundaryConditions->initialize(this->functionSpace_,
                                                    neumannBoundaryConditionElements);
            neumannBoundaryConditions->setDeformationGradientField(
                this->deformationGradientField(this->nestedSolver_));

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
            this->updateNeumannBoundaryConditions(this->nestedSolver_,
                                                    neumannBoundaryConditions);

            LOG(DEBUG) << "read data from precice complete, traction values: "
                        << tractionValues_;

}

template <typename NestedSolver>
void PreciceAdapterNestedSolver<NestedSolver>::setDirichletBoundaryConditions_impl(typename PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData &preciceData){
     std::vector<std::pair<global_no_t, std::array<double, 6>>>
                newDirichletBCValues;

            // if this rank has no data, do not set any boundary conditions
            if (preciceData.preciceMesh->nNodesLocal != 0) {
                // loop over nodes
                const int nNodesX = this->functionSpace_->nNodesLocalWithoutGhosts(0);
                const int nNodesY = this->functionSpace_->nNodesLocalWithoutGhosts(1);
                const int nNodesZ = this->functionSpace_->nNodesLocalWithoutGhosts(2);

                // set node index in z direction for bottom surface
                int nodeIndexZ = 0;

                // for top surface
                if (preciceData.preciceMesh->face ==
                    decltype(preciceData.preciceMesh.ioType)::face2Mas) {
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
                        this->functionSpace_->meshPartition()->getDofNoGlobalPetsc(
                            dofNoLocal);

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
            this->updateDirichletBoundaryConditions(this->nestedSolver_,
                                                    newDirichletBCValues);
        

        std::vector<double> displacementValues_; 
        std::vector<double> velocityValues_; 
        std::vector<double> tractionValues_; 
        std::vector<Vec3> displacementVectors_; 
        std::vector<Vec3> velocityVectors_;
        std::vector<Vec3> tractionVectors_; 
        std::vector<double> scalarValues_; 
        std::vector<double> scalarValuesOfMesh_; 
        std::vector<Vec3> geometryValues_; 

}
// --------------------------------------------------
// FastMonodomainSolver<T1>
template <typename T1>
std::shared_ptr<typename PreciceAdapterNestedSolver<
    FastMonodomainSolver<T1>>::FunctionSpace>
PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::functionSpace(NestedSolverType
                                                             &nestedSolver) {
  return nestedSolver.data().functionSpace();
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
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::
    getTractionValues(NestedSolverType &nestedSolver,
                      const std::vector<dof_no_t> &dofNosLocal,
                      std::vector<double> &tractionValues) {}

template <typename T1>
Vec PreciceAdapterNestedSolver<
    FastMonodomainSolver<T1>>::currentState(NestedSolverType
                                                            &nestedSolver) {}

template <typename T1>
std::shared_ptr<FieldVariable::FieldVariable<
    typename PreciceAdapterNestedSolver<
        FastMonodomainSolver<T1>>::FunctionSpace,
    9>>
PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::
    deformationGradientField(NestedSolverType &nestedSolver) {}

template <typename T1>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::reset(NestedSolverType
                                                     &nestedSolver) {
  // nestedSolver.timeStepping1().reset();
  nestedSolver.reset();
}

//! save fibers checkpoint
template <typename T1>
void PreciceAdapterNestedSolver<
    FastMonodomainSolver<T1>>::saveFiberData(NestedSolverType
                                                             &nestedSolver) {}

//! load fibers checkpoint
template <typename T1>
void PreciceAdapterNestedSolver<FastMonodomainSolver<T1>>::loadFiberData(NestedSolverType
                                                             &nestedSolver) {}
                  
// --------------------------------------------------
// MuscleContractionSolver<T1,T2>
template <typename T1, typename T2>
std::shared_ptr<typename PreciceAdapterNestedSolver<
    MuscleContractionSolver<T1, T2>>::FunctionSpace>
PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>>::functionSpace(NestedSolverType
                                                             &nestedSolver) {
  return nestedSolver.data().functionSpace();
}

template <typename T1, typename T2>
void PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>>::
    addDirichletBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::vector<
            typename SpatialDiscretization::DirichletBoundaryConditionsBase<
                FunctionSpace, 6>::ElementWithNodes>
            &dirichletBoundaryConditionElements) {
  // add the dirichlet bc values
  bool overwriteBcOnSameDof = true;
  nestedSolver.dynamicHyperelasticitySolver()
      ->addDirichletBoundaryConditions(dirichletBoundaryConditionElements,
                                       overwriteBcOnSameDof);
}

//! update existing boundary conditions with new values
template <typename T1, typename T2>
void PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>>::
    updateDirichletBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::vector<std::pair<global_no_t, std::array<double, 6>>>
            newDirichletBoundaryConditionValues) {
  nestedSolver.dynamicHyperelasticitySolver()
      ->updateDirichletBoundaryConditions(newDirichletBoundaryConditionValues);
}

template <typename T1, typename T2>
void PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>>::
    updateNeumannBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>>
            neumannBoundaryConditions) {
  nestedSolver.dynamicHyperelasticitySolver()
      ->hyperelasticitySolver()
      .updateNeumannBoundaryConditions(neumannBoundaryConditions);
}

//! get the displacement and velocity vectors of the given local dof nos
template <typename T1, typename T2>
void PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>>::
    getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                  const std::vector<dof_no_t> &dofNosLocal,
                                  std::vector<double> &displacementValues,
                                  std::vector<double> &velocityValues) {
  // get the displacement values
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.data().displacements()->getValues(dofNosLocal,
                                                                 values);

  // store displacement values in interleaved order (array of struct)
  int nVectors = values.size();
  displacementValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++) {
    displacementValues[3 * i + 0] = values[i][0];
    displacementValues[3 * i + 1] = values[i][1];
    displacementValues[3 * i + 2] = values[i][2];
  }

  // get the velocity values
  values.clear();
  nestedSolver.data().velocities()->getValues(dofNosLocal,
                                                              values);

  // store velocity values in interleaved order (array of struct)
  nVectors = values.size();
  velocityValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++) {
    velocityValues[3 * i + 0] = values[i][0];
    velocityValues[3 * i + 1] = values[i][1];
    velocityValues[3 * i + 2] = values[i][2];
  }
}

//! get the traction vectors of the given local dof nos
template <typename T1, typename T2>
void PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>>::
    getTractionValues(NestedSolverType &nestedSolver,
                      const std::vector<dof_no_t> &dofNosLocal,
                      std::vector<double> &tractionValues) {
  /*std::vector<Vec3> values0;
  nestedSolver.timeStepping2().data().materialTraction()->getValuesWithoutGhosts(values0);

  std::stringstream s;
  for (int i = 0; i < values0.size(); i++)
  {
    s << " " << values0[i][2];
  }
  LOG(INFO) << values0.size() << " local traction values in total
  (MuscleContractionSolver), "
    << "\ndofNosLocal: " << dofNosLocal << "\nall values: " << values0 << "\nz
  values: " << s.str();
*/
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.data().materialTraction()->getValues(dofNosLocal,
                                                                    values);

  int nVectors = values.size();
  tractionValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++) {
    tractionValues[3 * i + 0] = values[i][0];
    tractionValues[3 * i + 1] = values[i][1];
    tractionValues[3 * i + 2] = values[i][2];
  }
}

template <typename T1, typename T2>
Vec PreciceAdapterNestedSolver<
    MuscleContractionSolver<T1, T2>>::currentState(NestedSolverType
                                                            &nestedSolver) {
  return nestedSolver.dynamicHyperelasticitySolver()
      ->currentState();
}

template <typename T1, typename T2>
std::shared_ptr<FieldVariable::FieldVariable<
    typename PreciceAdapterNestedSolver<
        MuscleContractionSolver<T1, T2>>::FunctionSpace,
    9>>
PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>>::
    deformationGradientField(NestedSolverType &nestedSolver) {
  return nestedSolver.dynamicHyperelasticitySolver()
      ->hyperelasticitySolver()
      .data()
      .deformationGradient();
}

template <typename T1, typename T2>
void PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>>::reset(NestedSolverType
                                                     &nestedSolver) {
  // nestedSolver.timeStepping1().reset();
  nestedSolver.reset();
}

//! save fibers checkpoint
template <typename T1, typename T2>
void PreciceAdapterNestedSolver<
    MuscleContractionSolver<T1, T2>>::saveFiberData(NestedSolverType
                                                             &nestedSolver) {}

//! load fibers checkpoint
template <typename T1, typename T2>
void PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>>::loadFiberData(NestedSolverType
                                                             &nestedSolver) {}

// --------------------------------------------------
// Coupling<T1,MuscleContractionSolver<T2,T3>>
template <typename T1, typename T2, typename T3>
std::shared_ptr<typename PreciceAdapterNestedSolver<
    Control::Coupling<T1, MuscleContractionSolver<T2, T3>>>::FunctionSpace>
PreciceAdapterNestedSolver<Control::Coupling<
    T1, MuscleContractionSolver<T2, T3>>>::functionSpace(NestedSolverType
                                                             &nestedSolver) {
  return nestedSolver.timeStepping2().data().functionSpace();
}

template <typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<
    Control::Coupling<T1, MuscleContractionSolver<T2, T3>>>::
    addDirichletBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::vector<
            typename SpatialDiscretization::DirichletBoundaryConditionsBase<
                FunctionSpace, 6>::ElementWithNodes>
            &dirichletBoundaryConditionElements) {
  // add the dirichlet bc values
  bool overwriteBcOnSameDof = true;
  nestedSolver.timeStepping2()
      .dynamicHyperelasticitySolver()
      ->addDirichletBoundaryConditions(dirichletBoundaryConditionElements,
                                       overwriteBcOnSameDof);
}

//! update existing boundary conditions with new values
template <typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<
    Control::Coupling<T1, MuscleContractionSolver<T2, T3>>>::
    updateDirichletBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::vector<std::pair<global_no_t, std::array<double, 6>>>
            newDirichletBoundaryConditionValues) {
  nestedSolver.timeStepping2()
      .dynamicHyperelasticitySolver()
      ->updateDirichletBoundaryConditions(newDirichletBoundaryConditionValues);
}

template <typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<
    Control::Coupling<T1, MuscleContractionSolver<T2, T3>>>::
    updateNeumannBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>>
            neumannBoundaryConditions) {
  nestedSolver.timeStepping2()
      .dynamicHyperelasticitySolver()
      ->hyperelasticitySolver()
      .updateNeumannBoundaryConditions(neumannBoundaryConditions);
}

//! get the displacement and velocity vectors of the given local dof nos
template <typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<
    Control::Coupling<T1, MuscleContractionSolver<T2, T3>>>::
    getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                  const std::vector<dof_no_t> &dofNosLocal,
                                  std::vector<double> &displacementValues,
                                  std::vector<double> &velocityValues) {
  // get the displacement values
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.timeStepping2().data().displacements()->getValues(dofNosLocal,
                                                                 values);

  // store displacement values in interleaved order (array of struct)
  int nVectors = values.size();
  displacementValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++) {
    displacementValues[3 * i + 0] = values[i][0];
    displacementValues[3 * i + 1] = values[i][1];
    displacementValues[3 * i + 2] = values[i][2];
  }

  // get the velocity values
  values.clear();
  nestedSolver.timeStepping2().data().velocities()->getValues(dofNosLocal,
                                                              values);

  // store velocity values in interleaved order (array of struct)
  nVectors = values.size();
  velocityValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++) {
    velocityValues[3 * i + 0] = values[i][0];
    velocityValues[3 * i + 1] = values[i][1];
    velocityValues[3 * i + 2] = values[i][2];
  }
}

//! get the traction vectors of the given local dof nos
template <typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<
    Control::Coupling<T1, MuscleContractionSolver<T2, T3>>>::
    getTractionValues(NestedSolverType &nestedSolver,
                      const std::vector<dof_no_t> &dofNosLocal,
                      std::vector<double> &tractionValues) {
  /*std::vector<Vec3> values0;
  nestedSolver.timeStepping2().data().materialTraction()->getValuesWithoutGhosts(values0);

  std::stringstream s;
  for (int i = 0; i < values0.size(); i++)
  {
    s << " " << values0[i][2];
  }
  LOG(INFO) << values0.size() << " local traction values in total
  (MuscleContractionSolver), "
    << "\ndofNosLocal: " << dofNosLocal << "\nall values: " << values0 << "\nz
  values: " << s.str();
*/
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.timeStepping2().data().materialTraction()->getValues(dofNosLocal,
                                                                    values);

  int nVectors = values.size();
  tractionValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++) {
    tractionValues[3 * i + 0] = values[i][0];
    tractionValues[3 * i + 1] = values[i][1];
    tractionValues[3 * i + 2] = values[i][2];
  }
}

template <typename T1, typename T2, typename T3>
Vec PreciceAdapterNestedSolver<Control::Coupling<
    T1, MuscleContractionSolver<T2, T3>>>::currentState(NestedSolverType
                                                            &nestedSolver) {
  return nestedSolver.timeStepping2()
      .dynamicHyperelasticitySolver()
      ->currentState();
}

template <typename T1, typename T2, typename T3>
std::shared_ptr<FieldVariable::FieldVariable<
    typename PreciceAdapterNestedSolver<
        Control::Coupling<T1, MuscleContractionSolver<T2, T3>>>::FunctionSpace,
    9>>
PreciceAdapterNestedSolver<
    Control::Coupling<T1, MuscleContractionSolver<T2, T3>>>::
    deformationGradientField(NestedSolverType &nestedSolver) {
  return nestedSolver.timeStepping2()
      .dynamicHyperelasticitySolver()
      ->hyperelasticitySolver()
      .data()
      .deformationGradient();
}

template <typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<
    T1, MuscleContractionSolver<T2, T3>>>::reset(NestedSolverType
                                                     &nestedSolver) {
  // nestedSolver.timeStepping1().reset();
  nestedSolver.timeStepping2().reset();
}

//! save fibers checkpoint
template <typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<
    T1, MuscleContractionSolver<T2, T3>>>::saveFiberData(NestedSolverType
                                                             &nestedSolver) {
  // nestedSolver.timeStepping1().reset();
  LOG(INFO) << "saveFiberDataCheckpoint";
  nestedSolver.timeStepping1().saveFiberDataCheckpoint();
}

//! load fibers checkpoint
template <typename T1, typename T2, typename T3>
void PreciceAdapterNestedSolver<Control::Coupling<
    T1, MuscleContractionSolver<T2, T3>>>::loadFiberData(NestedSolverType
                                                             &nestedSolver) {
  // nestedSolver.timeStepping1().reset();
  LOG(INFO) << "loadFiberDataCheckpoint";
  nestedSolver.timeStepping1().restoreFiberDataCheckpoint();
}

// --------------------------------------------------
// DynamicHyperelasticitySolver

template <typename Material>
std::shared_ptr<typename TimeSteppingScheme::DynamicHyperelasticitySolver<
    Material>::FunctionSpace>
PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<
    Material>>::functionSpace(NestedSolverType &nestedSolver) {
  return nestedSolver.data().functionSpace();
}

template <typename Material>
void PreciceAdapterNestedSolver<
    TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
    addDirichletBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::vector<
            typename SpatialDiscretization::DirichletBoundaryConditionsBase<
                FunctionSpace, 6>::ElementWithNodes>
            &dirichletBoundaryConditionElements) {
  // add the dirichlet bc values
  LOG(INFO) << "add dirichlet BC \n";

  bool overwriteBcOnSameDof = true;
  nestedSolver.addDirichletBoundaryConditions(
      dirichletBoundaryConditionElements, overwriteBcOnSameDof);
}

//! update existing boundary conditions with new values
template <typename Material>
void PreciceAdapterNestedSolver<
    TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
    updateDirichletBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::vector<std::pair<global_no_t, std::array<double, 6>>>
            newDirichletBoundaryConditionValues) {
  LOG(INFO) << "update dirichlet BC \n";
  nestedSolver.updateDirichletBoundaryConditions(
      newDirichletBoundaryConditionValues);
}

template <typename Material>
void PreciceAdapterNestedSolver<
    TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
    updateNeumannBoundaryConditions(
        NestedSolverType &nestedSolver,
        std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>>
            neumannBoundaryConditions) {
  LOG(INFO) << "update neumann BC \n";
  nestedSolver.hyperelasticitySolver().updateNeumannBoundaryConditions(
      neumannBoundaryConditions);
}

//! get the displacement and velocity vectors of the given local dof nos
template <typename Material>
void PreciceAdapterNestedSolver<
    TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
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
  values.clear();
  nestedSolver.data().velocities()->getValues(dofNosLocal, values);

  // store velocity values in interleaved order (array of struct)
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
void PreciceAdapterNestedSolver<
    TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
    getTractionValues(NestedSolverType &nestedSolver,
                      const std::vector<dof_no_t> &dofNosLocal,
                      std::vector<double> &tractionValues) {
  /*std::vector<Vec3> values0;
  nestedSolver.hyperelasticitySolver().data().materialTraction()->getValuesWithoutGhosts(values0);

  std::stringstream s;
  for (int i = 0; i < values0.size(); i++)
  {
    s << " " << values0[i][2];
  }
  LOG(INFO) << values0.size() << " local traction values in total
  (DynamicHyperelasticitySolver), "
    << "\ndofNosLocal: " << dofNosLocal << "\nall values: " << values0 << "\nz
  values: " << s.str();
*/
  static std::vector<Vec3> values;
  values.clear();
  nestedSolver.hyperelasticitySolver().data().materialTraction()->getValues(
      dofNosLocal, values);

  int nVectors = values.size();
  tractionValues.resize(nVectors * 3);

  for (int i = 0; i < nVectors; i++) {
    tractionValues[3 * i + 0] = values[i][0];
    tractionValues[3 * i + 1] = values[i][1];
    tractionValues[3 * i + 2] = values[i][2];
  }
}

template <typename Material>
Vec PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<
    Material>>::currentState(NestedSolverType &nestedSolver) {
  return nestedSolver.currentState();
}

template <typename Material>
std::shared_ptr<FieldVariable::FieldVariable<
    typename PreciceAdapterNestedSolver<
        TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
        FunctionSpace,
    9>>
PreciceAdapterNestedSolver<TimeSteppingScheme::DynamicHyperelasticitySolver<
    Material>>::deformationGradientField(NestedSolverType &nestedSolver) {
  return nestedSolver.hyperelasticitySolver().data().deformationGradient();
}

template <typename Material>
void PreciceAdapterNestedSolver<
    TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
    reset(NestedSolverType &nestedSolver) {
  // nestedSolver.reset();
}

template <typename Material>
void PreciceAdapterNestedSolver<
    TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
    saveFiberData(NestedSolverType &nestedSolver) {}

template <typename Material>
void PreciceAdapterNestedSolver<
    TimeSteppingScheme::DynamicHyperelasticitySolver<Material>>::
    loadFiberData(NestedSolverType &nestedSolver) {}

// --------------------------------------------------
// HyperelasticitySolver

template <typename Material>
std::shared_ptr<typename SpatialDiscretization::HyperelasticitySolver<
    Material>::FunctionSpace>
PreciceAdapterNestedSolver<SpatialDiscretization::HyperelasticitySolver<
    Material>>::functionSpace(NestedSolverType &nestedSolver) {
  return nestedSolver.data().functionSpace();
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