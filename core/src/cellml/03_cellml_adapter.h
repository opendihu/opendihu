#pragma once

#include <Python.h> // has to be the first included header

#include <vector>

#include "interfaces/splittable.h"
#include "cellml/02_callback_handler.h"

/** This is a class that contains cellml equations and can be used with a time
 * stepping scheme. The nStates template parameter specifies the number of state
 * variables that should be used with the integrator. It is necessary that this
 * value is fixed at compile time because the timestepping scheme needs to know
 * which field variable types it has to construct. This class can also be
 * computed easily in multiple instances along the nodes of a mesh. The number
 * of instances is deduced from the mesh.
 *
 *  The states values are not stored inside the class but in the time stepping
 * scheme that is used to integrate the cellml problems.
 *
 *  Naming:
 *   Algebraic (opendihu) = KNOWN (OpenCMISS) = Algebraic (OpenCOR)
 *   Parameter (opendihu, OpenCMISS) = KNOWN (OpenCMISS), in OpenCOR also
 * algebraic Constant - these are constants that are only present in the source
 * files State: state variable Rate: the time derivative of the state variable,
 * i.e. the increment value in an explicit Euler stepping
 */
template <int nStates_, int nAlgebraics_ = 9,
          typename FunctionSpaceType = FunctionSpace::Generic>
class CellmlAdapter
    : public CallbackHandler<nStates_, nAlgebraics_, FunctionSpaceType>,
      public Splittable {
public:
  //! this class needs to define a function space in which its solution
  //! variables live. This does not matter at all for a CellML problem,
  //! therefore Generic is sufficient. But when using in an operator splitting
  //! with FEM as second operator part, it has to be compatible to that and thus
  //! needs to be set correctly.
  typedef FunctionSpaceType FunctionSpace; //< FunctionSpace type

  //! constructor from context object
  CellmlAdapter(DihuContext context);

  //! constructor from other CellmlAdapter with new functionSpace (and therefore
  //! different number of instances), preserves everything else, initialize does
  //! not need to be called afterwards this is needed by the LoadBalancing class
  CellmlAdapter(const CellmlAdapter &rhs,
                std::shared_ptr<FunctionSpace> functionSpace);

  //! return nStates_
  static constexpr int nStates();

  //! initialize callback functions and rhs
  void initialize();

  //! initialize timestepping
  void initializeForImplicitTimeStepping();

  //! reset the object to uninitialized state
  void reset();

  //! evaluate rhs
  void evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output,
                                                 int timeStepNo,
                                                 double currentTime);

  //! return the mesh
  std::shared_ptr<FunctionSpaceType> functionSpace();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);

  //! set initial values and return true or don't do anything and return false
  bool setInitialValues(
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, nStates_>>
          initialValues);

  //! get a vector with the names of the states
  void getComponentNames(std::vector<std::string> &stateNames) override;

  //! if the class should handle Dirichlet boundary conditions, this does not
  //! apply here
  void
  setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled){};

  //! after this call, getSlotConnectorData() will be called, transfer algebraic
  //! field variable to global representation
  void prepareForGetSlotConnectorData();

  //! the FastMonodomainSolver accesses the internals of CellmlAdapter
  template <int a, int b, typename c> friend class FastMonodomainSolverBase;

protected:
  //! check if the callback function "setSpecificParameters" needs to be called
  //! and if so, execute the call
  void checkCallbackParameters(double currentTime);

  //! check if the callback function "setSpecificStates" needs to be called and
  //! if so, execute the call
  void checkCallbackStates(double currentTime, double *statesLocal);

  //! check if the callback function "handleResult" needs to be called and if
  //! so, execute the call
  void checkCallbackAlgebraics(double currentTime, double *statesLocal,
                               double *algebraicsLocal);

  //! compute equilibrium of states for option "initializeStatesToEquilibrium"
  void initializeToEquilibriumValues(
      std::array<double, nStates_> &statesInitialValues);
};

#include "cellml/03_cellml_adapter.tpp"
