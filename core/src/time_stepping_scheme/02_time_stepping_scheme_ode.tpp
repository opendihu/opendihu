#include "time_stepping_scheme/02_time_stepping_scheme_ode.h"

#include <Python.h> // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme {

template <typename DiscretizableInTimeType>
TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::
    TimeSteppingSchemeOdeBaseDiscretizable(DihuContext context,
                                           std::string name)
    : TimeSteppingSchemeOdeBase<typename DiscretizableInTimeType::FunctionSpace,
                                DiscretizableInTimeType::nComponents()>::
          TimeSteppingSchemeOdeBase(context, name),
      discretizableInTime_(this->context_), initialized_(false) {
  // initialize output writers
  this->outputWriterManager_.initialize(this->context_,
                                        this->specificSettings_);

  // create dirichlet Boundary conditions object
  this->dirichletBoundaryConditions_ =
      std::make_shared<SpatialDiscretization::DirichletBoundaryConditions<
          typename DiscretizableInTimeType::FunctionSpace,
          DiscretizableInTimeType::nComponents()>>(this->context_);
}

template <typename DiscretizableInTimeType>
DiscretizableInTimeType &TimeSteppingSchemeOdeBaseDiscretizable<
    DiscretizableInTimeType>::discretizableInTime() {
  return this->discretizableInTime_;
}

template <typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBaseDiscretizable<
    DiscretizableInTimeType>::setRankSubset(Partition::RankSubset rankSubset) {
  TimeSteppingSchemeOdeBase<
      typename DiscretizableInTimeType::Functionspace,
      DiscretizableInTimeType::nComponents()>::setRankSubset(rankSubset);
  discretizableInTime_.setRankSubset(rankSubset);
}

template <typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::reset() {
  TimeSteppingSchemeOdeBase<typename DiscretizableInTimeType::FunctionSpace,
                            DiscretizableInTimeType::nComponents()>::reset();
  discretizableInTime_.reset();
  this->data_.reset();

  initialized_ = false;
}

template <typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBaseDiscretizable<
    DiscretizableInTimeType>::initialize() {
  if (initialized_)
    return;

  TimeSteppingScheme::initialize();
  LOG(TRACE) << "TimeSteppingSchemeOdeBase::initialize";

  // disable boundary condition handling in finite element method, because
  // Dirichlet BC have to be handled in the system matrix here
  discretizableInTime_.setBoundaryConditionHandlingEnabled(false);

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver(this->name_);

  // parse description for solverStructureVisualizer, if there was any
  std::string description;
  if (this->specificSettings_.hasKey("description"))
    description = this->specificSettings_.getOptionString("description", "");
  DihuContext::solverStructureVisualizer()->setSolverDescription(description);

  // indicate in solverStructureVisualizer that now a child solver will be
  // initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // initialize underlying DiscretizableInTime object, also with time step width
  discretizableInTime_.initialize();
  discretizableInTime_
      .initializeForImplicitTimeStepping(); // this performs extra
                                            // initialization for implicit
                                            // timestepping methods, i.e. it
                                            // sets the inverse lumped mass
                                            // matrix

  // indicate in solverStructureVisualizer that the child solver initialization
  // is done
  DihuContext::solverStructureVisualizer()->endChild();

  // retrieve the function space from the discretizable in time object, this is
  // used for the data object
  std::shared_ptr<typename DiscretizableInTimeType::FunctionSpace>
      functionSpace = discretizableInTime_.functionSpace();

  assert(functionSpace->meshPartition()); // assert that the function space was
                                          // retrieved correctly
  assert(this->data_);
  this->data_->setFunctionSpace(functionSpace);

  // set component names in data
  std::vector<std::string> componentNames;
  discretizableInTime_.getComponentNames(componentNames);
  this->data_->setComponentNames(componentNames);

  // create the vectors in the data object
  this->data_->initialize();

  // pass the solution field variable on to the discretizableInTime object, by
  // this e.g. CellMLAdapter gets the solution field variable
  discretizableInTime_.setSolutionVariable(this->data_->solution());

  // pass the full slot connector data to the discretizableInTime object, by
  // this e.g. CellMLAdapter can set additional algebraic variables that will be
  // used further
  discretizableInTime_.setSlotConnectorData(
      this->data_->getSlotConnectorData());

  // parse boundary conditions, needs functionSpace set
  // initialize dirichlet boundary conditions object which parses dirichlet
  // boundary condition dofs and values from config
  this->dirichletBoundaryConditions_->initialize(this->specificSettings_,
                                                 this->data_->functionSpace(),
                                                 "dirichletBoundaryConditions");

  // set initial values from settings

  // load initial values as specified in config under the "CellML" section
  if (!discretizableInTime_.setInitialValues(this->data_->solution())) {
    LOG(DEBUG) << "initial values were not set by DiscretizableInTime, set now";

    // if it did not initialize it,
    // load initial values from config under the timestepping section
    this->setInitialValues();
  } else {
    // the "CellML" section already specified initial values. Keep them.
    LOG(DEBUG) << "initial values were set by DiscretizableInTime";
  }
  VLOG(1) << "initial solution vector: " << *this->data_->solution();

  // output initial values but don't increment counter
  this->outputWriterManager_.writeOutput(*this->data_, 0, 0.0, 0);

  this->data_->print();

  // set the slot connector data to the solverStructureVisualizer which can then
  // include it in the structured solvers diagram
  DihuContext::solverStructureVisualizer()->setSlotConnectorData(
      this->getSlotConnectorData());

  initialized_ = true;
}

template <typename DiscretizableInTimeType>
void TimeSteppingSchemeOdeBaseDiscretizable<
    DiscretizableInTimeType>::prepareForGetSlotConnectorData() {
  discretizableInTime_.prepareForGetSlotConnectorData();
}

template <typename DiscretizableInTimeType>
std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<
    typename DiscretizableInTimeType::FunctionSpace,
    DiscretizableInTimeType::nComponents()>>
TimeSteppingSchemeOdeBaseDiscretizable<
    DiscretizableInTimeType>::dirichletBoundaryConditions() {
  return dirichletBoundaryConditions_;
}

} // namespace TimeSteppingScheme
