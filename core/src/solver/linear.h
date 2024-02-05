#pragma once

#include <Python.h> // has to be the first included header
#include "solver/solver.h"

#include <petscksp.h>
#include <memory>

namespace Solver {

/**
 * A linear solver configuration that contains all parameters for PETSc and the
 * KSP object. Usually you do not create an object of this class, but use
 * SolverManager for this task
 */
class Linear : public Solver {
public:
  //! construct solver from python settings
  Linear(PythonConfig specificSettings, MPI_Comm mpiCommunicator,
         std::string name);

  //! destructor
  virtual ~Linear();

  //! initialize the KSP (Krylov space solver in Petsc) and PC (preconditioner
  //! in Petsc) objects, has to be called before solve
  virtual void initialize();

  //! return the KSP object that is used for solving
  std::shared_ptr<KSP> ksp();

  //! perform the solve, @return: if the solution succeeded
  bool solve(Vec rightHandSide, Vec solution, std::string message = "");

  //! return the number of iterations of the last solve
  int lastNumberOfIterations();

  //! dump files containing rhs, solution and system matrix
  void dumpMatrixRightHandSideSolution(Vec rightHandSide, Vec solution);

  //! parse the solver and preconditioner type from strings that come from the
  //! settings
  //! @param solverType [in]
  //! @param preconditionerType [in]
  //! @param kspType [out]
  //! @param pcType [out]
  static void parseSolverTypes(std::string solverType,
                               std::string preconditionerType, KSPType &kspType,
                               PCType &pcType);

protected:
  //! parse options from settings
  void parseOptions();

  //! set options for KSP object
  void setupKsp(KSP ksp);

  std::shared_ptr<KSP> ksp_; //< the PETSc KSP (Krylov subspace) object
  double
      relativeTolerance_; //< relative solver tolerance of the residuum norm
                          //relative to the initial value of the residual norm
  double absoluteTolerance_; //< absolute solver tolerance of the residuum norm
  long int maxIterations_;   //< maximum number of iterations
  int lastNumberOfIterations_; //< the number of iterations of the previous
                               //solve

  std::string dumpFormat_;   //< format to use for dumping matrices and vectors
  std::string dumpFilename_; //< filename used for dumping matrices and vectors,
                             //empty for no dump

  KSPType kspType_;          //< the solver type
  PCType pcType_;            //< the preconditioner type
  MPI_Comm mpiCommunicator_; //< the MPI communcator to use

  std::shared_ptr<Vec>
      temporaryVectorLeft_; //< temporary vector for computation of residual for
                            //direct solvers
  std::shared_ptr<Vec>
      temporaryVectorRight_; //< temporary vector for computation of residual
                             //for direct solvers
  std::shared_ptr<Vec> residual_; //< residual vector for direct solvers

  std::string nIterationsLogKey_;  //< the keyword for the log with which the
                                   //number of iterations will be stored
  std::string residualNormLogKey_; //< the keyword for the log with which the
                                   //residual norm gets stored
  std::string
      nIterationsTotalLogKey_; //< the keyword for the log with which the total
                               //number of iterations gets stored

  std::string solverType_; //< the type of the solver as given in the settings
  std::string preconditionerType_; //< the type of the preconditioner, as given
                                   //in the settings
};

} // namespace Solver
