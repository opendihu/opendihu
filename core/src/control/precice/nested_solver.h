#pragma once

#include <Python.h> // has to be the first included header

#include "specialized_solver/muscle_contraction_solver.h"

#include "time_stepping_scheme/crank_nicolson.h"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver.h"

#include "precice/precice.hpp"
#include "control/precice/initialize.h"

namespace Control {

// Forward declaration, required for the anonymous enums
template <typename NestedSolver>
class PreciceAdapterInitialize;

// template <typename NestedSolver>
// struct PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData;
/** This is a base class of the precice adapter that contains functionality that
 * depends on the type of the nested solver. All solvers that should be able to
 * use precice surface coupling have to implement this interface.
 */

template <typename NestedSolver>
class PreciceAdapterNestedSolver : public Runnable {};

/** Partial specialization for tendon or pure mechanics solver in a coupling
 * scheme, muscle contraction solver (nonlinear elasticity with active stress)
 */
template <typename T1>
class PreciceAdapterNestedSolver<FastMonodomainSolver<T1>> {
public:

  std::vector<double> displacementValues_;
  std::vector<double> velocityValues_;
  std::vector<double> tractionValues_;
  std::vector<Vec3> displacementVectors_;
  std::vector<Vec3> velocityVectors_;
  std::vector<Vec3> tractionVectors_;
  std::vector<double> scalarValues_;
  std::vector<double> scalarValuesOfMesh_;
  std::vector<Vec3> geometryValues_;
  //! define the type of the nested solver
  typedef FastMonodomainSolver<T1> NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolverType::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<FunctionSpace> functionSpace(NestedSolverType &nestedSolver);

  template<typename SurfaceDataVector, typename VolumeDataVector>
  void preciceReadData(NestedSolverType &nestedSolver, std::shared_ptr<precice::Participant> &preciceParticipant, SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData) {
    LOG(DEBUG) << "read data from precice";
    double preciceDt = preciceParticipant->getMaxTimeStepSize();

    using SlotConnectorDataType = typename
    NestedSolverType::SlotConnectorDataType;
    std::shared_ptr<SlotConnectorDataType> slotConnectorData =
        nestedSolver.getSlotConnectorData();

    // loop over data
    for (auto &preciceData : preciceVolumeData) {
      if (preciceData.ioType == PreciceAdapterInitialize<NestedSolverType>::PreciceVolumeData::ioRead)
                                    {
        int nEntries = preciceData.preciceMesh->nNodesLocal;

        if (preciceData.isGeometryField) {
          nEntries = preciceData.preciceMesh->nNodesLocal * 3;
        }

        // allocate temporary memory
        scalarValues_.resize(nEntries);

        // get all data at once
        preciceParticipant->readData(
            preciceData.preciceMesh->preciceMeshName,
            preciceData.preciceDataName,
            preciceData.preciceMesh->preciceVertexIds, preciceDt,
            scalarValues_);

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
                std::cout<<geometryValues_[dofNoLocal][componentNo]<<std::endl;
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

  template<typename SurfaceDataVector, typename VolumeDataVector>
  void preciceWriteData(NestedSolverType &nestedSolver, std::shared_ptr<precice::Participant> &preciceParticipant, SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData) {
    // ReadWriteDataBase ReadWriteDataBase_;
    // ReadWriteDataBase_.preciceWriteVolumeData(preciceParticipant);
  }

  //! initialize dirichlet boundary conditions by adding new dofs and prescribed
  //! values for all bottom or top nodes
  void addDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<
          typename SpatialDiscretization::DirichletBoundaryConditionsBase<
              FunctionSpace, 6>::ElementWithNodes>
          &dirichletBoundaryConditionElements);

  //! update existing boundary conditions with new values
  void updateDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<std::pair<global_no_t, std::array<double, 6>>>
          newDirichletBoundaryConditionValues);

  //! update the neumann boundary conditions by replacing the complete object
  void updateNeumannBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>>
          neumannBoundaryConditions);

  //! get the displacement and velocity vectors of the given local dof nos
  void getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                     const std::vector<dof_no_t> &dofNosLocal,
                                     std::vector<double> &displacementValues,
                                     std::vector<double> &velocityValues);

  //! get the traction vectors of the given local dof nos
  void getTractionValues(NestedSolverType &nestedSolver,
                         const std::vector<dof_no_t> &dofNosLocal,
                         std::vector<double> &tractionValues);

  //! get at Petsc Vec that stores all values of the current state, to be used
  //! to store and restore checkpoints
  Vec currentState(NestedSolverType &nestedSolver);

  //! get the field variable of the deformation gradient
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 9>>
  deformationGradientField(NestedSolverType &nestedSolver);

  void reset(NestedSolverType &nestedSolver);

protected:
  void saveFiberData(NestedSolverType &nestedSolver);

  void loadFiberData(NestedSolverType &nestedSolver);
};

/** Partial specialization for tendon or pure mechanics solver in a coupling
 * scheme, muscle contraction solver (nonlinear elasticity with active stress)
 */
template <typename T1, typename T2>
class PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>> {
public:
  //! define the type of the nested solver
  typedef MuscleContractionSolver<T1, T2> NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolverType::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<FunctionSpace> functionSpace(NestedSolverType &nestedSolver);

  template<typename SurfaceDataVector, typename VolumeDataVector>
  void preciceReadData(NestedSolverType &nestedSolver, std::shared_ptr<precice::Participant> &preciceParticipant, SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData) {
    // ReadWriteDataBase ReadWriteDataBase_;
    // ReadWriteDataBase_.preciceReadVolumeData(preciceParticipant);
  }

  template<typename SurfaceDataVector, typename VolumeDataVector>
  void preciceWriteData(NestedSolverType &nestedSolver, std::shared_ptr<precice::Participant> &preciceParticipant, SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData) {
    // ReadWriteDataBase ReadWriteDataBase_;
    // ReadWriteDataBase_.preciceWriteVolumeData(preciceParticipant);
  }

  //! initialize dirichlet boundary conditions by adding new dofs and prescribed
  //! values for all bottom or top nodes
  void addDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<
          typename SpatialDiscretization::DirichletBoundaryConditionsBase<
              FunctionSpace, 6>::ElementWithNodes>
          &dirichletBoundaryConditionElements);

  //! update existing boundary conditions with new values
  void updateDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<std::pair<global_no_t, std::array<double, 6>>>
          newDirichletBoundaryConditionValues);

  //! update the neumann boundary conditions by replacing the complete object
  void updateNeumannBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>>
          neumannBoundaryConditions);

  //! get the displacement and velocity vectors of the given local dof nos
  void getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                     const std::vector<dof_no_t> &dofNosLocal,
                                     std::vector<double> &displacementValues,
                                     std::vector<double> &velocityValues);

  //! get the traction vectors of the given local dof nos
  void getTractionValues(NestedSolverType &nestedSolver,
                         const std::vector<dof_no_t> &dofNosLocal,
                         std::vector<double> &tractionValues);

  //! get at Petsc Vec that stores all values of the current state, to be used
  //! to store and restore checkpoints
  Vec currentState(NestedSolverType &nestedSolver);

  //! get the field variable of the deformation gradient
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 9>>
  deformationGradientField(NestedSolverType &nestedSolver);

  void reset(NestedSolverType &nestedSolver);

  void saveFiberData(NestedSolverType &nestedSolver);

  void loadFiberData(NestedSolverType &nestedSolver);
};

/** Partial specialization for tendon or pure mechanics solver in a coupling
 * scheme, muscle contraction solver (nonlinear elasticity with active stress)
 */
template <typename T1, typename T2, typename T3>
class PreciceAdapterNestedSolver<
    Control::Coupling<T1, MuscleContractionSolver<T2, T3>>> {
public:
  //! define the type of the nested solver
  typedef Control::Coupling<T1, MuscleContractionSolver<T2, T3>>
      NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef
      typename NestedSolverType::TimeStepping2Type::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<FunctionSpace> functionSpace(NestedSolverType &nestedSolver);
  
  template<typename SurfaceDataVector, typename VolumeDataVector>
  void preciceReadData(NestedSolverType &nestedSolver, std::shared_ptr<precice::Participant> &preciceParticipant, SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData) {
    // ReadWriteDataBase ReadWriteDataBase_;
    // ReadWriteDataBase_.preciceReadVolumeData(preciceParticipant);
  }

  template<typename SurfaceDataVector, typename VolumeDataVector>
  void preciceWriteData(NestedSolverType &nestedSolver, std::shared_ptr<precice::Participant> &preciceParticipant, SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData) {
    // ReadWriteDataBase ReadWriteDataBase_;
    // ReadWriteDataBase_.preciceWriteVolumeData(preciceParticipant);
  }

  //! initialize dirichlet boundary conditions by adding new dofs and prescribed
  //! values for all bottom or top nodes
  void addDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<
          typename SpatialDiscretization::DirichletBoundaryConditionsBase<
              FunctionSpace, 6>::ElementWithNodes>
          &dirichletBoundaryConditionElements);

  //! update existing boundary conditions with new values
  void updateDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<std::pair<global_no_t, std::array<double, 6>>>
          newDirichletBoundaryConditionValues);

  //! update the neumann boundary conditions by replacing the complete object
  void updateNeumannBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>>
          neumannBoundaryConditions);

  //! get the displacement and velocity vectors of the given local dof nos
  void getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                     const std::vector<dof_no_t> &dofNosLocal,
                                     std::vector<double> &displacementValues,
                                     std::vector<double> &velocityValues);

  //! get the traction vectors of the given local dof nos
  void getTractionValues(NestedSolverType &nestedSolver,
                         const std::vector<dof_no_t> &dofNosLocal,
                         std::vector<double> &tractionValues);

  //! get at Petsc Vec that stores all values of the current state, to be used
  //! to store and restore checkpoints
  Vec currentState(NestedSolverType &nestedSolver);

  //! get the field variable of the deformation gradient
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 9>>
  deformationGradientField(NestedSolverType &nestedSolver);

  void reset(NestedSolverType &nestedSolver);

  void saveFiberData(NestedSolverType &nestedSolver);

  void loadFiberData(NestedSolverType &nestedSolver);
};

/** Partial specialization for tendon or pure mechanics solver, dynamic
 * nonlinear elasticity
 */
template <typename Material>
class PreciceAdapterNestedSolver<
    TimeSteppingScheme::DynamicHyperelasticitySolver<Material>> {
public:
  //! define the type of the nested solver
  typedef TimeSteppingScheme::DynamicHyperelasticitySolver<Material>
      NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolverType::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<typename TimeSteppingScheme::DynamicHyperelasticitySolver<
      Material>::FunctionSpace>
  functionSpace(NestedSolverType &nestedSolver);
  
  template<typename SurfaceDataVector, typename VolumeDataVector>
  void preciceReadData(NestedSolverType &nestedSolver, std::shared_ptr<precice::Participant> &preciceParticipant, SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData) {
    // ReadWriteDataBase ReadWriteDataBase_;
    // ReadWriteDataBase_.preciceReadVolumeData(preciceParticipant);
  }

  template<typename SurfaceDataVector, typename VolumeDataVector>
  void preciceWriteData(NestedSolverType &nestedSolver, std::shared_ptr<precice::Participant> &preciceParticipant, SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData) {
    // ReadWriteDataBase ReadWriteDataBase_;
    // ReadWriteDataBase_.preciceWriteVolumeData(preciceParticipant);
  }
  //! initialize dirichlet boundary conditions by adding prescribed values for
  //! all bottom or top nodes
  void addDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<
          typename SpatialDiscretization::DirichletBoundaryConditionsBase<
              FunctionSpace, 6>::ElementWithNodes>
          &dirichletBoundaryConditionElements);

  //! update existing boundary conditions with new values
  void updateDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<std::pair<global_no_t, std::array<double, 6>>>
          newDirichletBoundaryConditionValues);

  //! update the neumann boundary conditions by replacing the complete object
  void updateNeumannBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>>
          neumannBoundaryConditions);

  //! get the displacement and velocity vectors of the given local dof nos
  void getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                     const std::vector<dof_no_t> &dofNosLocal,
                                     std::vector<double> &displacementValues,
                                     std::vector<double> &velocityValues);

  //! get the traction vectors of the given local dof nos
  void getTractionValues(NestedSolverType &nestedSolver,
                         const std::vector<dof_no_t> &dofNosLocal,
                         std::vector<double> &tractionValues);

  //! get at Petsc Vec that stores all values of the current state, to be used
  //! to store and restore checkpoints
  Vec currentState(NestedSolverType &nestedSolver);

  //! get the field variable of the deformation gradient
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 9>>
  deformationGradientField(NestedSolverType &nestedSolver);

  void reset(NestedSolverType &nestedSolver);

  void saveFiberData(NestedSolverType &nestedSolver);

  void loadFiberData(NestedSolverType &nestedSolver);
};

/** Partial specialization for tendon or pure mechanics solver, static nonlinear
 * elasticity
 */
template <typename Material>
class PreciceAdapterNestedSolver<
    SpatialDiscretization::HyperelasticitySolver<Material>> {
public:
  //! define the type of the nested solver
  typedef SpatialDiscretization::HyperelasticitySolver<Material>
      NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolverType::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<typename SpatialDiscretization::HyperelasticitySolver<
      Material>::FunctionSpace>
  functionSpace(NestedSolverType &nestedSolver);

  template<typename SurfaceDataVector, typename VolumeDataVector>
  void preciceReadData(NestedSolverType &nestedSolver, std::shared_ptr<precice::Participant> &preciceParticipant, SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData) {
    // ReadWriteDataBase ReadWriteDataBase_;
    // ReadWriteDataBase_.preciceReadVolumeData(preciceParticipant);
  }

  template<typename SurfaceDataVector, typename VolumeDataVector>
  void preciceWriteData(NestedSolverType &nestedSolver, std::shared_ptr<precice::Participant> &preciceParticipant, SurfaceDataVector &preciceSurfaceData, VolumeDataVector &preciceVolumeData) {
    // ReadWriteDataBase ReadWriteDataBase_;
    // ReadWriteDataBase_.preciceWriteVolumeData(preciceParticipant);
  }
  //! initialize dirichlet boundary conditions by adding prescribed values for
  //! all bottom or top nodes
  void addDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<
          typename SpatialDiscretization::DirichletBoundaryConditionsBase<
              FunctionSpace, 6>::ElementWithNodes>
          &dirichletBoundaryConditionElements);

  //! update existing boundary conditions with new values
  void updateDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<std::pair<global_no_t, std::array<double, 6>>>
          newDirichletBoundaryConditionValues);

  //! update the neumann boundary conditions by replacing the complete object
  void updateNeumannBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>>
          neumannBoundaryConditions);

  //! get the displacement and velocity vectors of the given local dof nos
  void getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                     const std::vector<dof_no_t> &dofNosLocal,
                                     std::vector<double> &displacementValues,
                                     std::vector<double> &velocityValues);

  //! get the traction vectors of the given local dof nos
  void getTractionValues(NestedSolverType &nestedSolver,
                         const std::vector<dof_no_t> &dofNosLocal,
                         std::vector<double> &tractionValues);

  //! get at Petsc Vec that stores all values of the current state, to be used
  //! to store and restore checkpoints
  Vec currentState(NestedSolverType &nestedSolver);

  //! get the field variable of the deformation gradient
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 9>>
  deformationGradientField(NestedSolverType &nestedSolver);

  void reset(NestedSolverType &nestedSolver);

  void saveFiberData(NestedSolverType &nestedSolver);

  void loadFiberData(NestedSolverType &nestedSolver);
};

} // namespace Control

#include "control/precice/nested_solver.tpp"