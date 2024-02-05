#include "data_management/time_stepping/time_stepping.h"

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
#include "partition/mesh_partition/01_mesh_partition.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

namespace Data {

template <typename FunctionSpaceType, int nComponents>
TimeStepping<FunctionSpaceType, nComponents>::TimeStepping(DihuContext context)
    : Data<FunctionSpaceType>(context) {
  this->debuggingName_ = "timestepping";
}

template <typename FunctionSpaceType, int nComponents>
TimeStepping<FunctionSpaceType, nComponents>::~TimeStepping() {
  // free PETSc objects
  if (this->initialized_) {
    // PetscErrorCode ierr;
    // ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template <typename FunctionSpaceType, int nComponents>
void TimeStepping<FunctionSpaceType, nComponents>::createPetscObjects() {
  LOG(DEBUG)
      << "TimeStepping<FunctionSpaceType,nComponents>::createPetscObjects("
      << nComponents << ")";
  assert(this->functionSpace_);
  assert(this->functionSpace_->meshPartition());

  if (componentNames_.empty()) {
    this->solution_ =
        this->functionSpace_->template createFieldVariable<nComponents>(
            "solution");
    this->increment_ =
        this->functionSpace_->template createFieldVariable<nComponents>(
            "increment");
  } else {
    // if there are component names stored, use them for construction of the
    // field variables
    this->solution_ =
        this->functionSpace_->template createFieldVariable<nComponents>(
            "solution", componentNames_);
    this->increment_ =
        this->functionSpace_->template createFieldVariable<nComponents>(
            "increment", componentNames_);
  }

  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();
  slotConnectorData_->addFieldVariable(this->solution_);

  // create additional field variables that appear as connector slots and can be
  // connected to discretizableInTime_ and enclosing solvers
  int nAdditionalFieldVariables = this->context_.getPythonConfig().getOptionInt(
      "nAdditionalFieldVariables", 0, PythonUtility::NonNegative);
  additionalFieldVariables_.resize(nAdditionalFieldVariables);

  for (int i = 0; i < nAdditionalFieldVariables; i++) {
    std::stringstream name;
    name << "additionalFieldVariable" << i;
    additionalFieldVariables_[i] =
        this->functionSpace_->template createFieldVariable<1>(name.str());

    slotConnectorData_->addFieldVariable2(additionalFieldVariables_[i]);
    LOG(DEBUG) << "  add field variable " << name.str();
  }

  LOG(DEBUG) << debuggingName_
             << ": initial slot names: " << slotConnectorData_->slotNames;

  // parse slot names of the additional field variables
  this->context_.getPythonConfig().getOptionVector(
      "additionalSlotNames", slotConnectorData_->slotNames);
  slotConnectorData_->slotNames.insert(
      slotConnectorData_->slotNames.begin(),
      std::string("")); // add a dummy slot name, it will be replaced by the
                        // nested solver in their setSlotConnectorData method

  // make sure that there are as many slot names as slots
  slotConnectorData_->slotNames.resize(slotConnectorData_->nSlots());

  LOG(DEBUG) << debuggingName_
             << ": final slot names: " << slotConnectorData_->slotNames;
}

template <typename FunctionSpaceType, int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>
TimeStepping<FunctionSpaceType, nComponents>::solution() {
  return this->solution_;
}

template <typename FunctionSpaceType, int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>
TimeStepping<FunctionSpaceType, nComponents>::increment() {
  return this->increment_;
}

template <typename FunctionSpaceType, int nComponents>
void TimeStepping<FunctionSpaceType, nComponents>::print() {
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->increment_;
  VLOG(4) << *this->solution_;
  VLOG(4) << "======================";
}

template <typename FunctionSpaceType, int nComponents>
void TimeStepping<FunctionSpaceType, nComponents>::setComponentNames(
    std::vector<std::string> componentNames) {
  componentNames_ = componentNames;
}

template <typename FunctionSpaceType, int nComponents>
std::shared_ptr<typename TimeStepping<FunctionSpaceType,
                                      nComponents>::SlotConnectorDataType>
TimeStepping<FunctionSpaceType, nComponents>::getSlotConnectorData() {
  return slotConnectorData_;
}

template <typename FunctionSpaceType, int nComponents>
typename TimeStepping<FunctionSpaceType,
                      nComponents>::FieldVariablesForOutputWriter
TimeStepping<FunctionSpaceType,
             nComponents>::getFieldVariablesForOutputWriter() {
  // recover additional field variables from slotConnectorData_, they may have
  // been changed by transfer
  assert(slotConnectorData_->variable2.size() >=
         additionalFieldVariables_.size());
  for (int i = 0; i < additionalFieldVariables_.size(); i++) {
    LOG(DEBUG) << " Data::TimeStepping::getFieldVariablesForOutputWriter(), "
               << " get field variable "
               << slotConnectorData_->variable2[i].values << ", \""
               << slotConnectorData_->variable2[i].values->name()
               << "\" for additionalFieldVariables_[" << i << "]";
    additionalFieldVariables_[i] = slotConnectorData_->variable2[i].values;
  }

  // these field variables will be written to output files
  return FieldVariablesForOutputWriter(
      std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType, 3>>(
          this->functionSpace_->geometryField()),
      solution_, additionalFieldVariables_);
}

//! output the given data for debugging
template <typename FunctionSpaceType, int nComponents>
std::string TimeStepping<FunctionSpaceType, nComponents>::getString(
    std::shared_ptr<typename TimeStepping<FunctionSpaceType,
                                          nComponents>::SlotConnectorDataType>
        data) {
  std::stringstream s;
  s << "<" << debuggingName_ << ":" << data << ">";

  return s.str();
}

} // namespace Data
