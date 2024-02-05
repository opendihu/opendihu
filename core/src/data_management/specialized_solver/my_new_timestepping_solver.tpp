#include "data_management/specialized_solver/my_new_timestepping_solver.h"

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
MyNewTimesteppingSolver<FunctionSpaceType>::MyNewTimesteppingSolver(
    DihuContext context)
    : Data<FunctionSpaceType>::Data(context) {}

template <typename FunctionSpaceType>
void MyNewTimesteppingSolver<FunctionSpaceType>::initialize() {
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();

  // create th slot connector data object
  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();

  // add all needed field variables to be transferred

  // add component 0 of fieldvariableA_
  slotConnectorData_->addFieldVariable(this->fieldVariableA_, 0);

  // There is addFieldVariable(...) and addFieldVariable2(...) for the two
  // different field variable types, Refer to
  // "slot_connection/slot_connector_data.h" for details.

  // you can also access settings from the python config here:
  std::string option1 = this->context_.getPythonConfig().getOptionString(
      "option1", "default string");
  LOG(DEBUG) << "In data object, parsed option1: [" << option1 << "].";

  // parse slot names for all slot connector data slots. here we only expect one
  // value because we have set one slot (fieldVariableA)
  this->context_.getPythonConfig().getOptionVector(
      "slotNames", slotConnectorData_->slotNames);

  // make sure that there are as many slot names as slots
  slotConnectorData_->slotNames.resize(slotConnectorData_->nSlots());
}

template <typename FunctionSpaceType>
void MyNewTimesteppingSolver<FunctionSpaceType>::createPetscObjects() {
  assert(this->functionSpace_);

  // Here, the actual field variables will be created.
  // Make sure, the number of components matches. The string is the name of the
  // field variable. It will also be used in the VTK output files.
  this->fieldVariableA_ =
      this->functionSpace_->template createFieldVariable<1>("a");
  this->fieldVariableB_ =
      this->functionSpace_->template createFieldVariable<3>("b");
}

// ... add a "getter" method for each fieldvariable with the same name as the
// field variable (but without underscore)
template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 1>>
MyNewTimesteppingSolver<FunctionSpaceType>::fieldVariableA() {
  return this->fieldVariableA_;
}

template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
MyNewTimesteppingSolver<FunctionSpaceType>::fieldVariableB() {
  return this->fieldVariableB_;
}

template <typename FunctionSpaceType>
std::shared_ptr<
    typename MyNewTimesteppingSolver<FunctionSpaceType>::SlotConnectorDataType>
MyNewTimesteppingSolver<FunctionSpaceType>::getSlotConnectorData() {
  // return the slot connector data object
  return this->slotConnectorData_;
}

template <typename FunctionSpaceType>
typename MyNewTimesteppingSolver<
    FunctionSpaceType>::FieldVariablesForOutputWriter
MyNewTimesteppingSolver<FunctionSpaceType>::getFieldVariablesForOutputWriter() {
  // these field variables will be written to output files by the output writer

  // get the geometry field, which is always needed, from the function space
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
      geometryField =
          std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType, 3>>(
              this->functionSpace_->geometryField());

  return std::make_tuple(
      geometryField, this->fieldVariableA_,
      this->fieldVariableB_ // add all field variables that should appear in the
                            // output file. Of course, this list has to match
                            // the type in the header file.
  );
}

} // namespace Data
