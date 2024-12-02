#include "data_management/control/multiple_instances.h"

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

template <typename FunctionSpaceType, typename NestedSolverType>
MapDofs<FunctionSpaceType, NestedSolverType>::MapDofs(DihuContext context)
    : Data<FunctionSpaceType>(context) {}

//! create the additionalFieldVariables_
template <typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType, NestedSolverType>::initialize(
    int nAdditionalFieldVariables, NestedSolverType &nestedSolver) {
  nAdditionalFieldVariables_ = nAdditionalFieldVariables;

  // call initialize of base class
  Data<FunctionSpaceType>::initialize();

  // store the nested solver
  nestedSolver_ = std::make_shared<NestedSolverType>(nestedSolver);

  // initialize slot connector data
  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();

  std::get<0>(*slotConnectorData_) = nestedSolver_->getSlotConnectorData();
  std::get<1>(*slotConnectorData_) =
      std::make_shared<SlotConnectorData<FunctionSpaceType, 1>>();

  // add all additional field variables
  for (int additionalFieldVariableNo = 0;
       additionalFieldVariableNo < nAdditionalFieldVariables_;
       additionalFieldVariableNo++) {
    std::get<1>(*slotConnectorData_)
        ->addFieldVariable(
            additionalFieldVariables_[additionalFieldVariableNo]);
  }

  // parse slot names of the additional field variables
  this->context_.getPythonConfig().getOptionVector(
      "additionalSlotNames", std::get<1>(*slotConnectorData_)->slotNames);
}

template <typename FunctionSpaceType, typename NestedSolverType>
bool MapDofs<FunctionSpaceType, NestedSolverType>::restoreState(
    const InputReader::Generic &r) {
  std::vector<std::vector<double>> stateData;
  stateData.resize(additionalFieldVariables_.size());
  for (size_t i = 0; i < additionalFieldVariables_.size(); i++) {
    if (!r.readDoubleVector(additionalFieldVariables_[i]->uniqueName().c_str(),
                            stateData[i])) {
      return false;
    }
  }
  for (size_t i = 0; i < additionalFieldVariables_.size(); i++) {
    additionalFieldVariables_[i]->setValues(stateData[i]);
  }
  return true;
}

template <typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType, NestedSolverType>::createPetscObjects() {
  for (int i = 0; i < nAdditionalFieldVariables_; i++) {
    std::stringstream name;
    name << "additionalFieldVariable" << i;
    std::shared_ptr<FieldVariableType> additionalFieldVariable =
        this->functionSpace_->template createFieldVariable<1>(name.str());
    additionalFieldVariable->setUniqueName(
        StringUtility::getFirstNE(this->uniquePrefix_, "control_map_dofs_") +
        name.str());

    additionalFieldVariables_.push_back(additionalFieldVariable);
  }
}

//! get a reference to the additional field variables
template <typename FunctionSpaceType, typename NestedSolverType>
std::vector<std::shared_ptr<
    typename MapDofs<FunctionSpaceType, NestedSolverType>::FieldVariableType>> &
MapDofs<FunctionSpaceType, NestedSolverType>::additionalFieldVariables() {
  return additionalFieldVariables_;
}

//! Get the data that will be transferred in the operator splitting or coupling
//! to the other term of the splitting/coupling. the transfer is done by the
//! slot_connector_data_transfer class
template <typename FunctionSpaceType, typename NestedSolverType>
std::shared_ptr<typename MapDofs<FunctionSpaceType,
                                 NestedSolverType>::SlotConnectorDataType>
MapDofs<FunctionSpaceType, NestedSolverType>::getSlotConnectorData() {
  // call getSlotConnectorData of the nested solver such that it can prepare the
  // connector slots
  nestedSolver_->getSlotConnectorData();

  return slotConnectorData_;
}

template <typename FunctionSpaceType, typename NestedSolverType>
typename MapDofs<FunctionSpaceType,
                 NestedSolverType>::FieldVariablesForOutputWriter
MapDofs<FunctionSpaceType,
        NestedSolverType>::getFieldVariablesForOutputWriter() {
  return std::make_tuple(additionalFieldVariables_);
}

template <typename FunctionSpaceType, typename NestedSolverType>
typename MapDofs<FunctionSpaceType,
                 NestedSolverType>::FieldVariablesForCheckpointing
MapDofs<FunctionSpaceType,
        NestedSolverType>::getFieldVariablesForCheckpointing() {
  return this->getFieldVariablesForOutputWriter();
}
} // namespace Data
