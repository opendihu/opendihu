#include "data_management/specialized_solver/muscle_contraction_solver.h"

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
MuscleContractionSolver<FunctionSpaceType>::MuscleContractionSolver(
    DihuContext context)
    : Data<FunctionSpaceType>::Data(context) {}

template <typename FunctionSpaceType>
void MuscleContractionSolver<FunctionSpaceType>::initialize() {
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();

  // create the slot connector data object
  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();

  // add all needed field variables to be transferred

  // add λ, λdot and γ as slot connectors
  slotConnectorData_->addFieldVariable(this->lambda_);
  slotConnectorData_->addFieldVariable(this->lambdaDot_);
  slotConnectorData_->addFieldVariable(this->gamma_);

  // There is addFieldVariable(...) and addFieldVariable2(...) for the two
  // different field variable types, Refer to
  // "slot_connection/slot_connector_data.h" for details.

  // parse slot names of the field variables
  this->context_.getPythonConfig().getOptionVector(
      "slotNames", slotConnectorData_->slotNames);

  // make sure that there are as many slot names as slots
  slotConnectorData_->slotNames.resize(slotConnectorData_->nSlots());
}

template <typename FunctionSpaceType>
bool MuscleContractionSolver<FunctionSpaceType>::restoreState(
    const InputReader::Generic &r) {
  std::vector<double> lambda, lambdaDot, gamma, displacements, velocities,
      activePK2Stress, pK2Stress, fiberDirection, materialTraction;
  if (!r.readDoubleVector(this->gamma_->uniqueName().c_str(), gamma)) {
    return false;
  }
  if (!r.readDoubleVector(this->lambda_->uniqueName().c_str(), lambda)) {
    return false;
  }
  if (!r.readDoubleVector(this->lambdaDot_->uniqueName().c_str(), lambdaDot)) {
    return false;
  }
  if (!r.readDoubleVector(this->displacements_->uniqueName().c_str(),
                          displacements)) {
    return false;
  }
  if (!r.readDoubleVector(this->velocities_->uniqueName().c_str(),
                          velocities)) {
    return false;
  }
  if (!r.readDoubleVector(this->activePK2Stress_->uniqueName().c_str(),
                          activePK2Stress)) {
    return false;
  }
  if (!r.readDoubleVector(this->pK2Stress_->uniqueName().c_str(), pK2Stress)) {
    return false;
  }
  if (!r.readDoubleVector(this->fiberDirection_->uniqueName().c_str(),
                          fiberDirection)) {
    return false;
  }
  if (!r.readDoubleVector(this->materialTraction_->uniqueName().c_str(),
                          materialTraction)) {
    return false;
  }

  std::array<std::vector<double>, 3> geometryValues;
  if (!r.template readDoubleVecD<3>(
          this->functionSpace_->geometryField().name().c_str(), geometryValues,
          "3D/")) {
    return false;
  }

  this->lambda_->setValues(lambda);
  this->lambdaDot_->setValues(lambdaDot);
  this->gamma_->setValues(gamma);
  LOG(INFO) << "=> setting displacements in muscle contraction";
  this->displacements_->setValues(displacements);
  LOG(INFO) << "=> setting displacements in muscle contraction done";
  this->velocities_->setValues(velocities);
  this->activePK2Stress_->setValues(activePK2Stress);
  this->pK2Stress_->setValues(pK2Stress);
  this->fiberDirection_->setValues(fiberDirection);
  this->materialTraction_->setValues(materialTraction);

  // for (size_t i = 0; i < 3; i++) {
  //   this->functionSpace_->geometryField().setValuesWithGhosts(
  //       i, geometryValues[i], INSERT_VALUES);
  // }

  return true;
}

template <typename FunctionSpaceType>
void MuscleContractionSolver<FunctionSpaceType>::createPetscObjects() {
  assert(this->functionSpace_);

  // Here, the actual field variables will be created.
  // The string is the name of the field variable. It will also be used in the
  // VTK output files.
  this->gamma_ = this->functionSpace_->template createFieldVariable<1>("γ");
  this->gamma_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "muscle_contraction_solver_") +
      "γ");
  this->lambda_ = this->functionSpace_->template createFieldVariable<1>("λ");
  this->lambda_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "muscle_contraction_solver_") +
      "λ");
  this->lambdaDot_ =
      this->functionSpace_->template createFieldVariable<1>("λdot");
  this->lambdaDot_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "muscle_contraction_solver_") +
      "λdot");
}

template <typename FunctionSpaceType>
void MuscleContractionSolver<FunctionSpaceType>::setFieldVariables(
    std::shared_ptr<
        MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType>
        displacements,
    std::shared_ptr<
        MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType>
        velocities,
    std::shared_ptr<
        MuscleContractionSolver<FunctionSpaceType>::StressFieldVariableType>
        activePK2Stress,
    std::shared_ptr<
        MuscleContractionSolver<FunctionSpaceType>::StressFieldVariableType>
        pK2Stress,
    std::shared_ptr<
        MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType>
        fiberDirection,
    std::shared_ptr<
        MuscleContractionSolver<FunctionSpaceType>::VectorFieldVariableType>
        materialTraction,
    bool setGeometryFieldForTransfer) {
  displacements_ = displacements;
  velocities_ = velocities;
  activePK2Stress_ = activePK2Stress;
  pK2Stress_ = pK2Stress;
  fiberDirection_ = fiberDirection;
  materialTraction_ = materialTraction;

  if (setGeometryFieldForTransfer) {
    slotConnectorData_->addGeometryField(
        std::make_shared<typename FunctionSpaceType::GeometryFieldType>(
            this->displacements_->functionSpace()->geometryField()));
  }

  // add material traction field variable (is stored in hyperelasticity solver)
  slotConnectorData_->addFieldVariable2(this->materialTraction_);

  // add displacements in x,y and z directions
  slotConnectorData_->addFieldVariable2(this->displacements_, 0);
  slotConnectorData_->addFieldVariable2(this->displacements_, 1);
  slotConnectorData_->addFieldVariable2(this->displacements_, 2);

  // There is addFieldVariable(...) and addFieldVariable2(...) for the two
  // different field variable types, Refer to
  // "slot_connection/slot_connector_data.h" for details.

  // parse slot names of the field variables
  this->context_.getPythonConfig().getOptionVector(
      "slotNames", slotConnectorData_->slotNames);

  // make sure that there are as many slot names as slots
  slotConnectorData_->slotNames.resize(slotConnectorData_->nSlots());
}

template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 1>>
MuscleContractionSolver<FunctionSpaceType>::lambda() {
  return this->lambda_;
}

template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 1>>
MuscleContractionSolver<FunctionSpaceType>::lambdaDot() {
  return this->lambdaDot_;
}

template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 1>>
MuscleContractionSolver<FunctionSpaceType>::gamma() {
  return this->gamma_;
}

template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
MuscleContractionSolver<FunctionSpaceType>::materialTraction() {
  return this->materialTraction_;
}

template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
MuscleContractionSolver<FunctionSpaceType>::displacements() {
  return this->displacements_;
}

template <typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
MuscleContractionSolver<FunctionSpaceType>::velocities() {
  return this->velocities_;
}

template <typename FunctionSpaceType>
std::shared_ptr<
    typename MuscleContractionSolver<FunctionSpaceType>::SlotConnectorDataType>
MuscleContractionSolver<FunctionSpaceType>::getSlotConnectorData() {
  // return the slot connector data object
  return this->slotConnectorData_;
}

template <typename FunctionSpaceType>
typename MuscleContractionSolver<
    FunctionSpaceType>::FieldVariablesForOutputWriter
MuscleContractionSolver<FunctionSpaceType>::getFieldVariablesForOutputWriter() {
  // these field variables will be written to output files by the output writer

  // get the geometry field, which is always needed, from the function space
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
      geometryField =
          std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType, 3>>(
              this->functionSpace_->geometryField());

  return std::make_tuple(
      geometryField,
      this->lambda_,          //< relative fiber stretch
      this->lambdaDot_,       //< contraction velocity
      this->gamma_,           //< gamma, the homogenized stress
      this->displacements_,   //< u, the displacements
      this->velocities_,      //< v, the velocities
      this->activePK2Stress_, //< the symmetric PK2 stress tensor of the active
                              // contribution in Voigt notation
      this->pK2Stress_, //< the symmetric PK2 stress tensor in Voigt notation
      this->fiberDirection_,  //< direction of fibers at current point
      this->materialTraction_ //< material traction
  );
}

template <typename FunctionSpaceType>
typename MuscleContractionSolver<
    FunctionSpaceType>::FieldVariablesForCheckpointing
MuscleContractionSolver<
    FunctionSpaceType>::getFieldVariablesForCheckpointing() {
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
      geometryField =
          std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType, 3>>(
              this->functionSpace_->geometryField());
  geometryField->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "muscle_contraction_solver_") +
      geometryField->name());

  return std::make_tuple(
      geometryField,
      this->lambda_,          //< relative fiber stretch
      this->lambdaDot_,       //< contraction velocity
      this->gamma_,           //< gamma, the homogenized stress
      this->displacements_,   //< u, the displacements
      this->velocities_,      //< v, the velocities
      this->activePK2Stress_, //< the symmetric PK2 stress tensor of the active
                              // contribution in Voigt notation
      this->pK2Stress_, //< the symmetric PK2 stress tensor in Voigt notation
      this->fiberDirection_,  //< direction of fibers at current point
      this->materialTraction_ //< material traction
  );
}
} // namespace Data
