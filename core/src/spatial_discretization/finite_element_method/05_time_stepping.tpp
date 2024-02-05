#include "spatial_discretization/finite_element_method/05_time_stepping.h"

#include <Python.h>
#include <iostream>
#include <sstream>
#include <petscksp.h>
#include <vector>
#include <numeric>
#include <omp.h>

#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "solver/linear.h"

namespace SpatialDiscretization {

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
FiniteElementMethodTimeStepping<
    FunctionSpaceType, QuadratureType, nComponents_,
    Term>::FiniteElementMethodTimeStepping(DihuContext context,
                                           std::shared_ptr<FunctionSpaceType>
                                               functionSpace)
    : AssembleRightHandSide<FunctionSpaceType, QuadratureType, nComponents_,
                            Term>(context, functionSpace),
      Splittable(), linearSolver_(nullptr), ksp_(nullptr) {}

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType,
                                     nComponents_, Term>::
    setBoundaryConditionHandlingEnabled(bool boundaryConditionHandlingEnabled) {
  BoundaryConditions<FunctionSpaceType, QuadratureType, nComponents_, Term,
                     Term>::
      setBoundaryConditionHandlingEnabled(boundaryConditionHandlingEnabled);
}

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType,
                                     nComponents_, Term>::initialize() {
  LOG(DEBUG) << "FiniteElementMethodTimeStepping::initialize";

  // call initialize of the parent class
  FiniteElementMethodBase<FunctionSpaceType, QuadratureType, nComponents_,
                          Term>::initialize();

  // initialize the linear solver
  this->initializeLinearSolver();

  // print a warning if this finite element class has output writers, because we
  // do not have solution data to write
  if (this->outputWriterManager_.hasOutputWriters()) {
    LOG(WARNING)
        << "You have specified output writers for a FiniteElementMethod which "
           "is used for a time stepping problem. "
           "The output will not contain any solution data, only geometry. "
           "Probably you want to get output from the time stepping scheme, "
           "then define the output writers there.";
  }
}

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
void FiniteElementMethodTimeStepping<
    FunctionSpaceType, QuadratureType, nComponents_,
    Term>::initializeForImplicitTimeStepping() {
  LOG(DEBUG)
      << "FiniteElementMethodTimeStepping::initializeForImplicitTimeStepping()";

  // initialize everything needed for implicit time stepping
  // currently this is executed regardless of explicit or implicit time stepping
  // scheme

  // initialize matrices
  this->data_.initializeMassMatrix();
  this->data_.initializeInverseLumpedMassMatrix();

  // compute the mass matrix
  this->setMassMatrix();

  // compute inverse lumped mass matrix
  this->setInverseLumpedMassMatrix();
}

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType,
                                     nComponents_, Term>::reset() {
  LOG(DEBUG) << " FiniteElementMethodTimeStepping::reset";
  FiniteElementMethodBase<FunctionSpaceType, QuadratureType, nComponents_,
                          Term>::reset();
  linearSolver_ = nullptr;
  ksp_ = nullptr;
}

//! hook to set initial values for a time stepping from this FiniteElement
//! context, return true if it has set the values or don't do anything and
//! return false
template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
bool FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType,
                                     nComponents_, Term>::
    setInitialValues(
        std::shared_ptr<
            FieldVariable::FieldVariable<FunctionSpaceType, nComponents_>>
            initialValues) {
  // Do not set initial values from within the "FiniteElements" section of the
  // config. (therefore return false) The initial values are set by the time
  // stepping scheme under its section.
  return false;
}

//! set the solution field variable in the data object, that actual data is
//! stored in the timestepping scheme object
template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType,
                                     nComponents_, Term>::
    setSolutionVariable(
        std::shared_ptr<
            FieldVariable::FieldVariable<FunctionSpaceType, nComponents_>>
            solution) {
  // this will be called by the time stepping scheme after initialize()
  this->data_.setSolutionVariable(solution);
}

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
void FiniteElementMethodTimeStepping<
    FunctionSpaceType, QuadratureType, nComponents_,
    Term>::setRankSubset(Partition::RankSubset rankSubset) {
  FiniteElementMethodBase<FunctionSpaceType, QuadratureType, nComponents_,
                          Term>::setRankSubset(rankSubset);
}

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType,
                                     nComponents_,
                                     Term>::initializeLinearSolver() {
  if (linearSolver_ == nullptr) {
    std::stringstream s;
    s << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
    LOG(DEBUG) << s.str()
               << ": FiniteElementMethodTimeStepping: initialize linearSolver";

    // retrieve linear solver
    // The mpiCommunicator is needed such that the solver knowns which ranks to
    // use (it could be a subset of all ranks).
    linearSolver_ =
        this->context_.solverManager()->template solver<Solver::Linear>(
            this->specificSettings_,
            this->data_.functionSpace()->meshPartition()->mpiCommunicator());
    ksp_ = linearSolver_->ksp();
  } else {
    std::stringstream s;
    s << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
    VLOG(2) << s.str() << ": linearSolver_ already set";
  }
}

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
constexpr int
FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, nComponents_,
                                Term>::nComponents() {
  return nComponents_;
}

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType,
                                     nComponents_, Term>::
    setSlotConnectorData(
        std::shared_ptr<
            Data::SlotConnectorData<FunctionSpaceType, nComponents_>>
            slotConnectorDataTimeStepping) {
  // This method is called once in initialize() of the timestepping scheme. It
  // provides the slotConnectorData of the timestepping.

  // Set the first slot name if given
  std::vector<std::string> &slotNames =
      slotConnectorDataTimeStepping->slotNames;
  std::string ownSolutionVariableSlotName =
      getSlotConnectorData()->slotNames.front();

  LOG(DEBUG) << "slotNames of timestepping initially: " << slotNames
             << " own slot names: " << getSlotConnectorData()->slotNames;

  // remove the first slot name which is always the one for the solution
  // variable
  slotNames.erase(slotNames.begin());

  // add all own slot names
  slotNames.insert(slotNames.begin(), ownSolutionVariableSlotName);

  LOG(DEBUG) << "in setSlotConnectorData of FiniteElementMethod, insert slot "
                "name of own solution variable: \""
             << ownSolutionVariableSlotName
             << "\", now slotNames of timestepping scheme: " << slotNames;
}

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
std::shared_ptr<typename FiniteElementMethodTimeStepping<
    FunctionSpaceType, QuadratureType, nComponents_,
    Term>::SlotConnectorDataType>
FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, nComponents_,
                                Term>::getSlotConnectorData() {
  // check for nans or infs
  return this->data_.getSlotConnectorData();
}

template <typename FunctionSpaceType, typename QuadratureType, int nComponents_,
          typename Term>
std::shared_ptr<FunctionSpaceType>
FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, nComponents_,
                                Term>::functionSpace() {
  return FiniteElementMethodBase<FunctionSpaceType, QuadratureType,
                                 nComponents_, Term>::functionSpace();
}

} // namespace SpatialDiscretization
