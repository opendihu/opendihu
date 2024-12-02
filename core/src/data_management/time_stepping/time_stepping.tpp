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
bool TimeStepping<FunctionSpaceType, nComponents>::restoreState(
    const InputReader::Generic &r) {
  std::vector<double> solution, increment;
  if (!r.readDoubleVector(this->solution_->uniqueName().c_str(), solution)) {
    return false;
  }
  if (!r.readDoubleVector(this->increment_->uniqueName().c_str(), increment)) {
    return false;
  }
  std::vector<std::vector<double>> additionalValues;
  for (int i = 0; i < additionalFieldVariables_.size(); i++) {
    std::vector<double> additional;
    if (!r.readDoubleVector(
            this->additionalFieldVariables_[i]->uniqueName().c_str(),
            additional)) {
      return false;
    }
    additionalValues.push_back(additional);
  }
  std::array<std::vector<double>, 3> geometryValues;
  if (!r.template readDoubleVecD<3>(
          this->functionSpace_->geometryField().name().c_str(), geometryValues,
          "3D/")) {
    return false;
  }

  this->solution_->setValues(solution);
  this->increment_->setValues(increment);
  for (int i = 0; i < additionalFieldVariables_.size(); i++) {
    slotConnectorData_->variable2[i].values->setValues(additionalValues[i]);
  }

  // for (size_t i = 0; i < 3; i++) {
  //   this->functionSpace_->geometryField().setValuesWithGhosts(
  //       i, geometryValues[i], INSERT_VALUES);
  // }

  return true;
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
    this->solution_->setUniqueName(
        StringUtility::getFirstNE(this->uniquePrefix_, "time_stepping_") +
        "solution");
    this->increment_ =
        this->functionSpace_->template createFieldVariable<nComponents>(
            "increment");
    this->increment_->setUniqueName(
        StringUtility::getFirstNE(this->uniquePrefix_, "time_stepping_") +
        "increment");
  } else {
    // if there are component names stored, use them for construction of the
    // field variables
    this->solution_ =
        this->functionSpace_->template createFieldVariable<nComponents>(
            "solution", componentNames_);
    this->solution_->setUniqueName(
        StringUtility::getFirstNE(this->uniquePrefix_, "time_stepping_") +
        "solution");
    this->increment_ =
        this->functionSpace_->template createFieldVariable<nComponents>(
            "increment", componentNames_);
    this->increment_->setUniqueName(
        StringUtility::getFirstNE(this->uniquePrefix_, "time_stepping_") +
        "increment");
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
    additionalFieldVariables_[i]->setUniqueName(
        StringUtility::getFirstNE(this->uniquePrefix_, "time_stepping_") +
        name.str());

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

template <typename FunctionSpaceType, int nComponents>
typename TimeStepping<FunctionSpaceType,
                      nComponents>::FieldVariablesForCheckpointing
TimeStepping<FunctionSpaceType,
             nComponents>::getFieldVariablesForCheckpointing() {
  // recover additional field variables from slotConnectorData_, they may have
  // been changed by transfer
  assert(slotConnectorData_->variable2.size() >=
         additionalFieldVariables_.size());
  for (int i = 0; i < additionalFieldVariables_.size(); i++) {
    LOG(DEBUG) << " Data::TimeStepping::getFieldVariablesForCheckpointing(), "
               << " get field variable "
               << slotConnectorData_->variable2[i].values << ", \""
               << slotConnectorData_->variable2[i].values->name()
               << "\" for additionalFieldVariables_[" << i << "]";
    additionalFieldVariables_[i] = slotConnectorData_->variable2[i].values;
  }
  auto geometryField =
      std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType, 3>>(
          this->functionSpace_->geometryField());

  // these field variables will be written to output files
  return FieldVariablesForCheckpointing(geometryField, solution_, increment_,
                                        additionalFieldVariables_);
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
