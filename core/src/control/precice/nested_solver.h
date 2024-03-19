#pragma once

#include <Python.h> // has to be the first included header

#include "control/precice/initialize.h"

#include "specialized_solver/muscle_contraction_solver.h"

#include "time_stepping_scheme/crank_nicolson.h"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver.h"

namespace Control {

/** This is a base class of the precice adapter that contains functionality that
 * depends on the type of the nested solver. All solvers that should be able to
 * use precice surface coupling have to implement this interface.
 */
template <typename NestedSolver>
class PreciceAdapterNestedSolver : public Runnable {
    public:

        virtual void preciceReadData()=0;
        virtual void preciceWriteData()=0;

        // void setNeumannBoundaryConditions(
        // typename PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData
        //   &preciceData)=0;
        // void setDirichletBoundaryConditions(
        // typename PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData
        //   &preciceData)=0;
    
    protected:

        void preciceReadData_impl(){};
        void preciceWriteData_impl(){};

        // void setNeumannBoundaryConditions_impl(
        // typename PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData
        //   &preciceData);
        // void setDirichletBoundaryConditions_impl(
        // typename PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData
        //   &preciceData);

        std::vector<double> displacementValues_; 
        std::vector<double> velocityValues_; 
        std::vector<double> tractionValues_; 
        std::vector<Vec3> displacementVectors_; 
        std::vector<Vec3> velocityVectors_;
        std::vector<Vec3> tractionVectors_; 
        std::vector<double> scalarValues_; 
        std::vector<double> scalarValuesOfMesh_; 
        std::vector<Vec3> geometryValues_; 

};

/** Partial specialization for tendon or pure mechanics solver in a coupling
 * scheme, muscle contraction solver (nonlinear elasticity with active stress)
 */
template <typename T1>
class PreciceAdapterNestedSolver<FastMonodomainSolver<T1>> {
public:
  //! define the type of the nested solver
  typedef FastMonodomainSolver<T1> NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef
      typename NestedSolverType::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<FunctionSpace> functionSpace(NestedSolverType &nestedSolver);

  void preciceReadData() override{
    this->preciceReadData_impl();}

  void preciceWriteData() override{
    this->preciceWriteData_impl();}

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
template <typename T1, typename T2>
class PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>> {
public:
  //! define the type of the nested solver
  typedef MuscleContractionSolver<T1, T2> NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef
      typename NestedSolverType::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<FunctionSpace> functionSpace(NestedSolverType &nestedSolver);

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