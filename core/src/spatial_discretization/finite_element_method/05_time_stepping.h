#pragma once

#include "spatial_discretization/finite_element_method/04_rhs.h"

#include "mesh/mesh.h"
#include "interfaces/discretizable_in_time.h"
#include "interfaces/multipliable.h"

namespace SpatialDiscretization {

/** class used for timestepping as for diffusion equation
 */
template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
class FiniteElementMethodTimeStepping
    : public AssembleRightHandSide<FunctionSpaceType, QuadratureType,
                                   nComponents_, Term>,
      public DiscretizableInTime<FunctionSpaceType, nComponents_>,
      public Splittable {
public:
  typedef Data::SlotConnectorData<FunctionSpaceType, nComponents_>
      SlotConnectorDataType; // type of return value of getSlotConnectorData

  //! constructor, if function space is not given, create new one according to
  //! settings if the function space is given as parameter, is has to be already
  //! initialize()d
  FiniteElementMethodTimeStepping(
      DihuContext context,
      std::shared_ptr<FunctionSpaceType> functionSpace = nullptr);

  using AssembleRightHandSide<FunctionSpaceType, QuadratureType, nComponents_,
                              Term>::initialize;

  //! return the compile-time constant number of variable components of the
  //! solution field variable
  static constexpr int nComponents();

  //! proceed time stepping by computing output = stiffnessMatrix*input, output
  //! back in strong form
  void evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output,
                                                 int timeStepNo,
                                                 double currentTime);

  //! initialize for use with timestepping
  void initialize();

  //! initialize for use with timestepping, this sets mass matrix and inverse
  //! lumped mass matrix
  void initializeForImplicitTimeStepping();

  //! reset the object to uninitialized state
  void reset();

  //! hook to set initial values for a time stepping from this FiniteElement
  //! context, return true if it has set the values or don't do anything and
  //! return false
  bool setInitialValues(
      std::shared_ptr<
          FieldVariable::FieldVariable<FunctionSpaceType, nComponents_>>
          initialValues);

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);

  //! enable or disable boundary condition handling on initialization, set to
  //! false to not care for boundary conditions
  void
  setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled);

  //! set the solution field variable in the data object, that actual data is
  //! stored in the timestepping scheme object
  void setSolutionVariable(
      std::shared_ptr<
          FieldVariable::FieldVariable<FunctionSpaceType, nComponents_>>
          solution);

  //! return the mesh that is stored in the data class
  std::shared_ptr<FunctionSpaceType> functionSpace();

  //! pass on the slot connector data object from the timestepping scheme object
  //! to be modified, this is needed for other DiscretizableInTime objects
  void setSlotConnectorData(
      std::shared_ptr<Data::SlotConnectorData<FunctionSpaceType, nComponents_>>
          slotConnectorDataTimeStepping);

  //! get the data that will be transferred in the operator splitting to the
  //! other term of the splitting the transfer is done by the
  //! slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  typedef FunctionSpaceType
      FunctionSpace; //< the FunctionSpace type needed for time stepping scheme

  friend class StiffnessMatrixTester; //< a class used for testing
protected:
  //! do nothing, needed for initialize of base class that is overridden anyway
  void setRightHandSide(){};

  //! Compute from the rhs in weak formulation the rhs vector in strong
  //! formulation
  void computeInverseMassMatrixTimesRightHandSide(Vec &result);

  //! initialize the linear solve that is needed for the solution of the
  //! implicit timestepping system
  void initializeLinearSolver();

  //! solves the linear system of equations resulting from the Implicit Euler
  //! method time discretization
  void solveLinearSystem(Vec &input, Vec &output);

  std::shared_ptr<Solver::Linear>
      linearSolver_; //< the linear solver used for inverting the mass matrix
  std::shared_ptr<KSP> ksp_; //< the linear solver context
};

} // namespace SpatialDiscretization

#include "spatial_discretization/finite_element_method/05_time_stepping.tpp"
#include "spatial_discretization/finite_element_method/05_time_stepping_explicit.tpp"
