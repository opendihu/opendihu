#include "data_management/parallel_fiber_estimation.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data {

template <typename FunctionSpaceType>
ParallelFiberEstimation<FunctionSpaceType>::ParallelFiberEstimation(
    DihuContext context)
    : Data<FunctionSpaceType>(context) {}

template <typename FunctionSpaceType>
ParallelFiberEstimation<FunctionSpaceType>::~ParallelFiberEstimation() {
  // free PETSc objects
  if (this->initialized_) {
    // PetscErrorCode ierr;
    // ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template <typename FunctionSpaceType>
void ParallelFiberEstimation<FunctionSpaceType>::setProblem(
    std::shared_ptr<FiniteElementMethodType> problem) {
  problem_ = problem;
}

template <typename FunctionSpaceType>
bool ParallelFiberEstimation<FunctionSpaceType>::restoreState(
    const InputReader::Generic &r) {
  std::vector<double> gradient;
  std::vector<double> dirichletValues;
  std::vector<double> jacobianConditionNumber;
  if (!r.readDoubleVector(this->gradient_->uniqueName().c_str(), gradient)) {
    return false;
  }
  if (!r.readDoubleVector(this->dirichletValues_->uniqueName().c_str(),
                          dirichletValues)) {
    return false;
  }
  if (!r.readDoubleVector(this->jacobianConditionNumber_->uniqueName().c_str(),
                          jacobianConditionNumber)) {
    return false;
  }

  this->gradient_->setValues(gradient);
  this->dirichletValues_->setValues(dirichletValues);
  this->jacobianConditionNumber_->setValues(jacobianConditionNumber);

  return this->problem_->data().restoreState(r);
}

template <typename FunctionSpaceType>
void ParallelFiberEstimation<FunctionSpaceType>::createPetscObjects() {
  LOG(DEBUG)
      << "ParallelFiberEstimation<FunctionSpaceType>::createPetscObjects()"
      << std::endl;
  assert(this->functionSpace_);

  // create field variables on local partition
  this->gradient_ =
      this->functionSpace_->template createFieldVariable<3>("gradient");
  this->gradient_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "parallel_fiber_estimation_") +
      "gradient");
  this->dirichletValues_ =
      this->functionSpace_->template createFieldVariable<1>("dirichletValues");
  this->dirichletValues_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "parallel_fiber_estimation_") +
      "dirichletValues");
  this->jacobianConditionNumber_ =
      this->functionSpace_->template createFieldVariable<1>(
          "jacobianConditionNumber");
  this->jacobianConditionNumber_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "parallel_fiber_estimation_") +
      "jacobianConditionNumber");
}

template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
ParallelFiberEstimation<FunctionSpaceType>::gradient() {
  return this->gradient_;
}

template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 1>>
ParallelFiberEstimation<FunctionSpaceType>::dirichletValues() {
  return this->dirichletValues_;
}

template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 1>>
ParallelFiberEstimation<FunctionSpaceType>::jacobianConditionNumber() {
  return this->jacobianConditionNumber_;
}

template <typename FunctionSpaceType>
void ParallelFiberEstimation<FunctionSpaceType>::print() {
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->gradient_;
  VLOG(4) << "======================";
}

template <typename FunctionSpaceType>
typename ParallelFiberEstimation<
    FunctionSpaceType>::FieldVariablesForOutputWriter
ParallelFiberEstimation<FunctionSpaceType>::getFieldVariablesForOutputWriter() {
  assert(problem_);
  return std::tuple_cat(
      problem_->data().getFieldVariablesForOutputWriter(),
      std::tuple<
          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>>(
          this->gradient_),
      std::tuple<
          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 1>>>(
          this->dirichletValues_),
      std::tuple<
          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 1>>>(
          this->jacobianConditionNumber_));
}

template <typename FunctionSpaceType>
typename ParallelFiberEstimation<
    FunctionSpaceType>::FieldVariablesForCheckpointing
ParallelFiberEstimation<
    FunctionSpaceType>::getFieldVariablesForCheckpointing() {
  return this->getFieldVariablesForOutputWriter();
}
} // namespace Data
