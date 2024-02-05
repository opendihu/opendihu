#include "specialized_solver/static_bidomain_solver.h"

#include <Python.h> // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/diagnostic_tool/performance_measurement.h"

namespace TimeSteppingScheme {

template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
StaticBidomainSolver<
    FiniteElementMethodPotentialFlow,
    FiniteElementMethodDiffusion>::StaticBidomainSolver(DihuContext context)
    : context_(context["StaticBidomainSolver"]), data_(this->context_),
      finiteElementMethodPotentialFlow_(this->context_["PotentialFlow"]),
      finiteElementMethodDiffusionTransmembrane_(this->context_["Activation"]),
      finiteElementMethodDiffusionExtracellular_(this->context_["Activation"]),
      rankSubset_(this->context_.rankSubset()), initialized_(false) {
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey")) {
    this->durationLogKey_ =
        specificSettings_.getOptionString("durationLogKey", "");
  }

  this->initialGuessNonzero_ =
      specificSettings_.getOptionBool("initialGuessNonzero", true);

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_,
                                        this->specificSettings_);
}

template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
void StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                          FiniteElementMethodDiffusion>::
    advanceTimeSpan(bool withOutputWritersEnabled) {
  LOG_SCOPE_FUNCTION;

  // start duration measurement, the name of the output variable can be set by
  // "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  /*
  Bidomain equation:
       K(sigma_i) Vm + K(sigma_i+sigma_e) phi_e = 0
    => K(sigma_i+sigma_e) phi_e = -K(sigma_i) Vm
   */

  // update right hand side: transmembraneFlow = -K(sigma_i) Vm
  PetscErrorCode ierr;
  ierr =
      MatMult(data_.rhsMatrix(), data_.transmembranePotential()->valuesGlobal(),
              data_.transmembraneFlow()->valuesGlobal());
  CHKERRV(ierr);

  // solve K(sigma_i+sigma_e) phi_e = transmembraneFlow for phi_e
  this->solveLinearSystem();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  if (withOutputWritersEnabled)
    this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
}

template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
void StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                          FiniteElementMethodDiffusion>::run() {
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
void StaticBidomainSolver<
    FiniteElementMethodPotentialFlow,
    FiniteElementMethodDiffusion>::setTimeSpan(double startTime,
                                               double endTime) {
  endTime_ = endTime;
}

//! start time of time interval to be simulated
template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
double StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                            FiniteElementMethodDiffusion>::startTime() {
  return endTime_; // there is no start time as the static bidomain solver
                   // solves a static equation
}

//! end time of simulation
template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
double StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                            FiniteElementMethodDiffusion>::endTime() {
  return endTime_;
}

template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
void StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                          FiniteElementMethodDiffusion>::initialize() {
  LOG_SCOPE_FUNCTION;

  if (this->initialized_)
    return;

  LOG(DEBUG) << "initialize static_bidomain_solver";
  assert(this->specificSettings_.pyObject());

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("StaticBidomainSolver");

  // indicate in solverStructureVisualizer that now a child solver will be
  // initialized
  DihuContext::solverStructureVisualizer()->beginChild("PotentialFlow");

  // initialize the potential flow finite element method, this also creates the
  // function space
  finiteElementMethodPotentialFlow_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization
  // is done
  DihuContext::solverStructureVisualizer()->endChild();

  // initialize the data object
  data_.setFunctionSpace(finiteElementMethodPotentialFlow_.functionSpace());
  data_.initialize();

  LOG(INFO) << "Run potential flow simulation for fiber directions.";

  // avoid that solver structure file is created, this should only be done after
  // the whole simulation has finished
  DihuContext::solverStructureVisualizer()->disable();

  // solve potential flow Laplace problem
  finiteElementMethodPotentialFlow_.run();

  // enable again
  DihuContext::solverStructureVisualizer()->enable();
  DihuContext::solverStructureVisualizer()->beginChild(
      "Activation Transmembrane");

  LOG(DEBUG) << "fem flow solution: "
             << *finiteElementMethodPotentialFlow_.data().solution();
  LOG(DEBUG) << "compute gradient field";

  // compute a gradient field from the solution of the potential flow
  data_.flowPotential()->setValues(
      *finiteElementMethodPotentialFlow_.data().solution());
  LOG(DEBUG) << "flow potential: " << *data_.flowPotential();

  // enable computation of the condition number of the jacobian in every element
  std::shared_ptr<FieldVariableType> jacobianConditionNumber = nullptr;
  if (this->specificSettings_.hasKey("enableJacobianConditionNumber"))
    if (this->specificSettings_.getOptionBool("enableJacobianConditionNumber",
                                              true))
      jacobianConditionNumber = data_.jacobianConditionNumber();

  data_.flowPotential()->computeGradientField(data_.fiberDirection(),
                                              jacobianConditionNumber);
  // note, fiberDirection is not normalized, it gets normalized in the
  // initialize method of the finite element data class
  // (diffusion_tensor_directional.tpp)

  VLOG(1) << "flow potential: " << *data_.flowPotential();
  VLOG(1) << "fiber direction: " << *data_.fiberDirection();

  // initialize the finite element class, from which only the stiffness matrix
  // is needed diffusion object without prefactor, for normal diffusion (1st
  // multidomain eq.)

  // initialize(direction, spatiallyVaryingPrefactor,
  // useAdditionalDiffusionTensor)
  finiteElementMethodDiffusionTransmembrane_.initialize(data_.fiberDirection(),
                                                        nullptr);

  // indicate in solverStructureVisualizer that the child solver initialization
  // is done
  DihuContext::solverStructureVisualizer()->endChild();
  DihuContext::solverStructureVisualizer()->beginChild(
      "Activation Extracellular");

  // direction, spatiallyVaryingPrefactor, useAdditionalDiffusionTensor=true
  finiteElementMethodDiffusionExtracellular_.initialize(data_.fiberDirection(),
                                                        nullptr, true);

  // indicate in solverStructureVisualizer that the child solver initialization
  // is done
  DihuContext::solverStructureVisualizer()->endChild();

  // initialize the matrix to be used for computing the rhs
  data_.rhsMatrix() = finiteElementMethodDiffusionTransmembrane_.data()
                          .stiffnessMatrix()
                          ->valuesGlobal();
  PetscErrorCode ierr;
  ierr = MatScale(data_.rhsMatrix(), -1);
  CHKERRV(ierr);

  LOG(DEBUG) << "initialize linear solver";

  // initialize linear solver
  if (linearSolver_ == nullptr) {
    // retrieve linear solver
    linearSolver_ =
        this->context_.solverManager()->template solver<Solver::Linear>(
            this->specificSettings_, this->rankSubset_->mpiCommunicator());
  }

  LOG(DEBUG) << "set system matrix to linear solver";

  // set matrix used for linear solver and preconditioner to ksp context
  Mat systemMatrix = finiteElementMethodDiffusionExtracellular_.data()
                         .stiffnessMatrix()
                         ->valuesGlobal();
  assert(this->linearSolver_->ksp());
  ierr =
      KSPSetOperators(*this->linearSolver_->ksp(), systemMatrix, systemMatrix);
  CHKERRV(ierr);

  // set the nullspace of the matrix
  // as we have Neumann boundary conditions, constant functions are in the
  // nullspace of the matrix
  MatNullSpace constantFunctions;
  MatNullSpaceCreate(rankSubset_->mpiCommunicator(), PETSC_TRUE, 0, nullptr,
                     &constantFunctions);
  MatSetNullSpace(systemMatrix, constantFunctions);
  MatSetNearNullSpace(systemMatrix, constantFunctions); // for multigrid methods
  MatNullSpaceDestroy(&constantFunctions);

  // set the slotConnectorData for the solverStructureVisualizer to appear in
  // the solver diagram
  DihuContext::solverStructureVisualizer()->setSlotConnectorData(
      getSlotConnectorData());

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
void StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                          FiniteElementMethodDiffusion>::reset() {
  this->initialized_ = false;
}

//! call the output writer on the data object, output files will contain
//! currentTime, with callCountIncrement !=1 output timesteps can be skipped
template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
void StaticBidomainSolver<
    FiniteElementMethodPotentialFlow,
    FiniteElementMethodDiffusion>::callOutputWriter(int timeStepNo,
                                                    double currentTime,
                                                    int callCountIncrement) {
  this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
}

template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
void StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                          FiniteElementMethodDiffusion>::solveLinearSystem() {

  VLOG(1) << "in solveLinearSystem";

  // configure that the initial value for the iterative solver is the value in
  // solution, not zero
  if (initialGuessNonzero_) {
    PetscErrorCode ierr;
    ierr = KSPSetInitialGuessNonzero(*this->linearSolver_->ksp(), PETSC_TRUE);
    CHKERRV(ierr);
  }

  // rename the involved vectors
  Vec rightHandSide = data_.transmembraneFlow()->valuesGlobal();
  Vec solution = data_.extraCellularPotential()->valuesGlobal();

  // check if there are nans
  data_.transmembraneFlow()->checkNansInfs();

  // dump vectors to be able to later check values
  // debugDumpData();

  // solve the system, KSPSolve(ksp,b,x)
#ifndef NDEBUG
  this->linearSolver_->solve(rightHandSide, solution,
                             "Linear system of bidomain problem solved");
#else
  this->linearSolver_->solve(rightHandSide, solution);
#endif
}

//! return whether the underlying discretizableInTime object has a specified
//! mesh type and is not independent of the mesh type
template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
void StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                          FiniteElementMethodDiffusion>::debugDumpData() {
  static int counter = 0;

  // compute matrix norm
  double norm1, normFrobenius, normInfinity;
  MatNorm(finiteElementMethodDiffusionExtracellular_.data()
              .stiffnessMatrix()
              ->valuesGlobal(),
          NORM_1, &norm1);
  MatNorm(finiteElementMethodDiffusionExtracellular_.data()
              .stiffnessMatrix()
              ->valuesGlobal(),
          NORM_FROBENIUS, &normFrobenius);
  MatNorm(finiteElementMethodDiffusionExtracellular_.data()
              .stiffnessMatrix()
              ->valuesGlobal(),
          NORM_INFINITY, &normInfinity);

  double directionNorm;
  VecNorm(data_.fiberDirection()->valuesGlobal(), NORM_2, &directionNorm);

  std::stringstream filename;
  filename << "norm_" << counter << ".txt";
  std::ofstream file0;
  file0.open(filename.str().c_str(), std::ios::out | std::ios::app);

  file0 << Control::PerformanceMeasurement::getParameter("scenarioName") << ";"
        << DihuContext::ownRankNoCommWorld() << "/"
        << DihuContext::nRanksCommWorld() << ";" << norm1 << ";"
        << normFrobenius << ";" << normInfinity << ";directionNorm;"
        << directionNorm << std::endl;
  file0.close();

  filename.str("");
  filename << "rhs_"
           << Control::PerformanceMeasurement::getParameter("scenarioName")
           << "_" << counter << "." << DihuContext::ownRankNoCommWorld()
           << ".bin";
  std::ofstream file;
  file.open(filename.str().c_str(),
            std::ios::out | std::ios::trunc | std::ios::binary);

  if (file.is_open()) {
    std::vector<double> values;
    data_.transmembraneFlow()->getValuesWithoutGhosts(values);

    // loop over rhs vector
    for (int i = 0; i < values.size(); i++) {
      file.write((char *)(&values[i]), 8);
    }
    file.close();
  } else {
    LOG(INFO) << "Could not open file";
  }

  filename.str("");
  filename << "initial_value_"
           << Control::PerformanceMeasurement::getParameter("scenarioName")
           << "_" << counter << "." << DihuContext::ownRankNoCommWorld()
           << ".bin";
  file.open(filename.str().c_str(),
            std::ios::out | std::ios::trunc | std::ios::binary);

  if (file.is_open()) {
    std::vector<double> values;
    data_.extraCellularPotential()->getValuesWithoutGhosts(values);

    // loop over rhs vector
    for (int i = 0; i < values.size(); i++) {
      file.write((char *)(&values[i]), 8);
    }
    file.close();
  } else {
    LOG(INFO) << "Could not open file";
  }
  counter++;
}

template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
typename StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                              FiniteElementMethodDiffusion>::Data &
StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                     FiniteElementMethodDiffusion>::data() {
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other
//! term of the splitting the transfer is done by the
//! slot_connector_data_transfer class
template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
std::shared_ptr<typename StaticBidomainSolver<
    FiniteElementMethodPotentialFlow,
    FiniteElementMethodDiffusion>::SlotConnectorDataType>
StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                     FiniteElementMethodDiffusion>::getSlotConnectorData() {
  return this->data_.getSlotConnectorData(); // transmembranePotential
}

//! output the given data for debugging
template <typename FiniteElementMethodPotentialFlow,
          typename FiniteElementMethodDiffusion>
std::string StaticBidomainSolver<FiniteElementMethodPotentialFlow,
                                 FiniteElementMethodDiffusion>::
    getString(std::shared_ptr<typename StaticBidomainSolver<
                  FiniteElementMethodPotentialFlow,
                  FiniteElementMethodDiffusion>::SlotConnectorDataType>
                  data) {
  std::stringstream s;
  s << "<StaticBidomain:" << data << ">";
  return s.str();
}

} // namespace TimeSteppingScheme
