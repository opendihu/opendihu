#include "specialized_solver/my_new_solver/my_new_static_solver.h"

#include <omp.h>
#include <sstream>

#ifdef HAVE_XBRAID
#include "braid.h"
#endif

template <typename NestedSolver>
MyNewStaticSolver<NestedSolver>::MyNewStaticSolver(DihuContext context)
    : Runnable(), context_(context["MyNewStaticSolver"]),
      nestedSolver_(this->context_), data_(this->context_),
      initialized_(false) {
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_,
                                        this->specificSettings_);

  // parse options
  int myOption = this->specificSettings_.getOptionInt("myOption", 1,
                                                      PythonUtility::Positive);

  LOG(DEBUG) << "myOption: " << myOption;
}

template <typename NestedSolver>
void MyNewStaticSolver<NestedSolver>::initialize() {
  LOG_SCOPE_FUNCTION;

  // make sure that we initialize only once, in the next call, initialized_ is
  // true
  if (initialized_)
    return;

  // initialize() will be called before the simulation starts.

  // add this solver to the solvers diagram, which is an ASCII art
  // representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver(
      "MyNewStaticSolver",
      true); // hasInternalConnectionToFirstNestedSolver=true
             // (the last argument) means slot connector
             // data is shared with the first subsolver
  // if you have your own slot connector data rather than the one of the
  // subsolver, call "addSolver" with false as second argument

  // indicate in solverStructureVisualizer that now a child solver will be
  // initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested solver
  nestedSolver_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization
  // is done
  DihuContext::solverStructureVisualizer()->endChild();

  // In order to initialize the data object and actuall create all variables, we
  // first need to assign a function space to the data object. A function space
  // object of type FunctionSpace<MeshType,BasisFunctionType> (see
  // "function_space/function_space.h") is an object that stores the mesh (e.g.,
  // all nodes and elements) as well as the basis function (e.g. linear Lagrange
  // basis functions). The NestedSolver solver already created a function space
  // that we should use. We already have a typedef "FunctionSpace" that is the
  // class of NestedSolver's function space type.
  std::shared_ptr<FunctionSpace> functionSpace =
      nestedSolver_.data().functionSpace();

  // Pass the function space to the data object. data_ stores field variables.
  // It needs to know the number of nodes and degrees of freedom (dof) from the
  // function space in order to create the vectors with the right size.
  data_.setFunctionSpace(functionSpace);

  // now call initialize, data will then create all variables (Petsc Vec's)
  data_.initialize();

  // it is also possible to pass some field variables from the data of the
  // NestedSolver to own data object
  data_.setSolutionVariable(nestedSolver_.data().solution());

  // set the slotConnectorData for the solverStructureVisualizer to appear in
  // the solver diagram
  DihuContext::solverStructureVisualizer()->setSlotConnectorData(
      getSlotConnectorData());

  // here is the space to initialize anything else that is needed for your
  // solver

#if 0
  // for example, initialize Braid here (not like this)
  braid_Core    core;
  MPI_Comm comm = MPI_COMM_WORLD;
  braid_Init(MPI_COMM_WORLD, comm, 0, 0, 0, nullptr,
        NULL, NULL, NULL, NULL, NULL, NULL,
        NULL, NULL, NULL, NULL, &core);
#endif
  initialized_ = true;
}

template <typename NestedSolver> void MyNewStaticSolver<NestedSolver>::run() {
  // initialize everything
  initialize();

  // perform the computation of this solver

  // call the nested solver
  nestedSolver_.run();

  // do something else
  executeMyHelperMethod();

  // write current output values using the output writers
  this->outputWriterManager_.writeOutput(this->data_);
}

template <typename NestedSolver> void MyNewStaticSolver<NestedSolver>::reset() {
  nestedSolver_.reset();

  initialized_ = false;
  // "uninitialize" everything
}

template <typename NestedSolver>
void MyNewStaticSolver<NestedSolver>::executeMyHelperMethod() {
  // this is the template for an own private method

  // Get the solution field variable from the data object. The Petsc object can
  // be retrieved with "valuesGlobal()".
  Vec solution = data_.solution()->valuesGlobal();

  // As an example we invert all the solution values.
  // Because "Vec"'s are actually pointers, this effects the actual data, no
  // copy-back is needed. Note the error handling with Petsc functions. Always
  // use "ierr = PetscFunction(); CHKERRV(ierr);"
  PetscErrorCode ierr;
  ierr = VecScale(solution, -1);
  CHKERRV(ierr);

  // get a field variable from data object:
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 1>>
      fieldVariableB = data_.fieldVariableB();

  // copy the solution from the nestedSolver_ to fieldVariableB.
  ierr = VecCopy(solution, fieldVariableB->valuesGlobal());
  CHKERRV(ierr);

  // add 1.0 to every entry in fieldVariableA, this already updates
  // fieldVariableA in data because it is a pointer. There is no need to copy
  // the values back.
  ierr = VecShift(fieldVariableB->valuesGlobal(), 1.0);
  CHKERRV(ierr);

  // For example you can also get the stiffness matrix of the nestedSolver_
  std::shared_ptr<PartitionedPetscMat<FunctionSpace>> stiffnessMatrix =
      nestedSolver_.data().stiffnessMatrix();
  Mat m = stiffnessMatrix->valuesGlobal();

  // e.g. add identity to m
  ierr = MatShift(m, 1.0);
  CHKERRV(ierr);
}

template <typename NestedSolver>
typename MyNewStaticSolver<NestedSolver>::Data &
MyNewStaticSolver<NestedSolver>::data() {
  // get a reference to the data object
  return data_;

  // The nestedSolver_ object also has a data object, we could also directly use
  // this and avoid having an own data object:
  //  return nestedSolver_.data();
}

//! get the data that will be transferred in the operator splitting to the other
//! term of the splitting the transfer is done by the
//! slot_connector_data_transfer class
template <typename NestedSolver>
std::shared_ptr<typename MyNewStaticSolver<NestedSolver>::SlotConnectorDataType>
MyNewStaticSolver<NestedSolver>::getSlotConnectorData() {
  //! This is relevant only, if this solver is part of a splitting or coupling
  //! scheme. Then this method returns the values/variables that will be
  // transferred to the other solvers. We can just reuse the values of the
  // nestedSolver_.
  return nestedSolver_.getSlotConnectorData();
}
