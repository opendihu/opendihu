#include "time_stepping_scheme/03_time_stepping_implicit.h"

#include <Python.h> // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include <petscksp.h>
#include "solver/solver_manager.h"
#include "solver/linear.h"
#include "data_management/time_stepping/time_stepping_implicit.h"

namespace TimeSteppingScheme {

template <typename DiscretizableInTimeType>
TimeSteppingImplicit<DiscretizableInTimeType>::TimeSteppingImplicit(
    DihuContext context, std::string name)
    : TimeSteppingSchemeOde<DiscretizableInTimeType>(context, name) {}

template <typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::initialize() {
  LOG_SCOPE_FUNCTION;

  if (this->initialized_)
    return;

  // initialize data objects that are needed for
  // TimeSteppingSchemeOde<DiscretizableInTimeType>::initialize();
  this->data_ = std::make_shared<Data::TimeSteppingImplicit<
      typename DiscretizableInTimeType::FunctionSpace,
      DiscretizableInTimeType::nComponents()>>(
      this->context_); // create data object for implicit euler
  this->dataImplicit_ = std::static_pointer_cast<Data::TimeSteppingImplicit<
      typename DiscretizableInTimeType::FunctionSpace,
      DiscretizableInTimeType::nComponents()>>(this->data_);

  TimeSteppingSchemeOde<DiscretizableInTimeType>::initialize();
  LOG(TRACE) << "TimeSteppingImplicit::initialize";

  timeStepWidthRelativeTolerance_ = this->specificSettings().getOptionDouble(
      "timeStepWidthRelativeTolerance", 1e-10, PythonUtility::NonNegative);
  if (this->specificSettings().hasKey("timeStepWidthRelativeToleranceAsKey")) {
    std::string relTolKey = this->specificSettings().getOptionString(
        "timeStepWidthRelativeToleranceAsKey",
        "timeStepWidthRelativeTolerance");
    Control::PerformanceMeasurement::setParameter(
        relTolKey, timeStepWidthRelativeTolerance_);
  }

  if (this->specificSettings().hasKey("durationInitTimeStepLogKey")) {
    this->durationInitTimeStepLogKey_ =
        this->specificSettings().getOptionString("durationInitTimeStepLogKey",
                                                 "");
  }

  this->initialized_ = true;
}

template <typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::initializeWithTimeStepWidth(
    double timeStepWidth) {
  LOG(TRACE) << "TimeSteppingImplicit::initializeWithTimeStepWidth("
             << timeStepWidth << ")";

  // check if the time step changed and a new initialization is neccessary
  if (this->initializedTimeStepWidth_ < 0.0 ||
      !this->dataImplicit_->systemMatrix()) {
    // first initialization
    LOG(DEBUG) << "initializeWithTimeStepWidth(" << timeStepWidth << ")";
  } else {
    // check if the time step size changed
    const double eps = this->timeStepWidthRelativeTolerance_;

    const double relDiff = (this->initializedTimeStepWidth_ - timeStepWidth) /
                           this->initializedTimeStepWidth_;
    if (-eps <= relDiff && relDiff <= eps) {
      LOG(DEBUG) << "do not re-initializeWithTimeStepWidth as relative "
                    "difference of time steps is small (tolerance: "
                 << eps << "): " << relDiff
                 << ". Old: " << this->initializedTimeStepWidth_
                 << ", new: " << timeStepWidth;
      return;
    }
    LOG(DEBUG) << "re-initializeWithTimeStepWidth as relative difference of "
                  "time steps is too large (tolerance: "
               << eps << "): " << relDiff
               << ". Old: " << this->initializedTimeStepWidth_
               << ", new: " << timeStepWidth;
  }

  if (this->durationInitTimeStepLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationInitTimeStepLogKey_);

  // perform actual (re-)initialization
  this->initializeWithTimeStepWidth_impl(timeStepWidth);

  if (this->durationInitTimeStepLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationInitTimeStepLogKey_);

  this->initializedTimeStepWidth_ = timeStepWidth;
}

template <typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::
    initializeWithTimeStepWidth_impl(double timeStepWidth) {
  // compute the system matrix
  this->setSystemMatrix(timeStepWidth);

  LOG(DEBUG) << "time_stepping_implicit applyInSystemMatrix, from "
                "TimeSteppingImplicit::initialize";
  // set the boundary conditions to system matrix, i.e. zero rows and columns of
  // Dirichlet BC dofs and set diagonal to 1
  this->dirichletBoundaryConditions_->applyInSystemMatrix(
      this->dataImplicit_->systemMatrix(), this->dataImplicit_->systemMatrix(),
      this->dataImplicit_->boundaryConditionsRightHandSideSummand());

  // initialize the linear solver that is used for solving the implicit system
  initializeLinearSolver();

  // set matrix used for linear system and preconditioner to ksp context
  Mat &systemMatrix = this->dataImplicit_->systemMatrix()->valuesGlobal();
  assert(this->ksp_);
  PetscErrorCode ierr;
  ierr = KSPSetOperators(*ksp_, systemMatrix, systemMatrix);
  CHKERRV(ierr);
}

template <typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::reset() {
  TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::reset();

  LOG(DEBUG) << "set linearSolver_ to nullptr";
  if (linearSolver_) {
    LOG(DEBUG) << "delete linear solver";
    this->context_.solverManager()->deleteSolver(linearSolver_->name());
  }

  linearSolver_ = nullptr;
}

template <typename DiscretizableInTimeType>
Data::TimeSteppingImplicit<typename DiscretizableInTimeType::FunctionSpace,
                           DiscretizableInTimeType::nComponents()> &
TimeSteppingImplicit<DiscretizableInTimeType>::dataImplicit() {
  return *dataImplicit_;
}

template <typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::solveLinearSystem(
    Vec &input, Vec &output) {
  if (!dataImplicit_)
    LOG(FATAL) << this->name_
               << ", solveLinearSystem, implicit data is not initialized, "
                  "initialized_="
               << this->initialized_;

  if (!this->dataImplicit_->systemMatrix())
    LOG(FATAL) << this->name_
               << ", solveLinearSystem, system matrix is not set, initialized_="
               << this->initialized_;

  // solve systemMatrix*output = input for output
  Mat &systemMatrix = this->dataImplicit_->systemMatrix()->valuesGlobal();

  PetscUtility::checkDimensionsMatrixVector(systemMatrix, input);

  if (VLOG_IS_ON(1)) {
    linearSolver_->solve(input, output,
                         "Linear system of implicit time stepping solved");
  } else {
    linearSolver_->solve(input, output);
  }
}

template <typename DiscretizableInTimeType>
void TimeSteppingImplicit<DiscretizableInTimeType>::initializeLinearSolver() {
  LOG(DEBUG) << "initializeLinearSolver, linearSolver_ == nullptr: "
             << (linearSolver_ == nullptr);
  if (linearSolver_ == nullptr) {
    LOG(DEBUG) << "Implicit time stepping: initialize linearSolver";

    // retrieve linear solver
    linearSolver_ =
        this->context_.solverManager()->template solver<Solver::Linear>(
            this->specificSettings_,
            this->data_->functionSpace()->meshPartition()->mpiCommunicator());
    ksp_ = linearSolver_->ksp();
  } else {
    VLOG(2) << ": linearSolver_ already set";
  }
}

/*
//! output the given data for debugging
template<typename DiscretizableInTimeType>
std::string TimeSteppingImplicit<DiscretizableInTimeType>::
getString(typename
TimeSteppingSchemeOde<DiscretizableInTimeType>::SlotConnectorDataType &data)
{
  return dataImplicit_->getString(data);
}*/

} // namespace TimeSteppingScheme
