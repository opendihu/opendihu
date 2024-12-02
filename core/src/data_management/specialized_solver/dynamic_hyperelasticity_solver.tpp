#include "data_management/specialized_solver/dynamic_hyperelasticity_solver.h"

namespace Data {

template <typename FunctionSpaceType>
DynamicHyperelasticitySolver<FunctionSpaceType>::DynamicHyperelasticitySolver(
    DihuContext context)
    : Data<FunctionSpaceType>::Data(context) {}

template <typename FunctionSpaceType>
bool DynamicHyperelasticitySolver<FunctionSpaceType>::restoreState(
    const InputReader::Generic &r) {
  std::vector<double> displacements, velocities, internalVirtualWork,
      accelerationTerm, externalVirtualWorkDead;
  if (!r.readDoubleVector(this->displacements_->uniqueName().c_str(),
                          displacements)) {
    return false;
  }
  if (!r.readDoubleVector(this->velocities_->uniqueName().c_str(),
                          velocities)) {
    return false;
  }
  if (!r.readDoubleVector(this->internalVirtualWork_->uniqueName().c_str(),
                          internalVirtualWork)) {
    return false;
  }
  if (!r.readDoubleVector(this->accelerationTerm_->uniqueName().c_str(),
                          accelerationTerm)) {
    return false;
  }
  if (!r.readDoubleVector(this->externalVirtualWorkDead_->uniqueName().c_str(),
                          externalVirtualWorkDead)) {
    return false;
  }

  std::array<std::vector<double>, 3> geometryValues;
  if (!r.template readDoubleVecD<3>(
          this->functionSpace_->geometryField().name().c_str(), geometryValues,
          "3D/")) {
    return false;
  }

  LOG(INFO) << "=> setting displacements in dynamic";
  displacements_->setValues(displacements);
  LOG(INFO) << "=> setting displacements in dynamic done";
  velocities_->setValues(velocities);
  internalVirtualWork_->setValues(internalVirtualWork);
  accelerationTerm_->setValues(accelerationTerm);
  externalVirtualWorkDead_->setValues(externalVirtualWorkDead);

  // for (size_t i = 0; i < 3; i++) {
  //   this->functionSpace_->geometryField().setValuesWithGhosts(
  //       i, geometryValues[i], INSERT_VALUES);
  // }

  return true;
}

template <typename FunctionSpaceType>
void DynamicHyperelasticitySolver<FunctionSpaceType>::createPetscObjects() {
  LOG(DEBUG) << "DynamicHyperelasticitySolver::createPetscObjects";

  assert(this->functionSpace_);

  std::vector<std::string> displacementsComponentNames({"x", "y", "z"});
  displacements_ = this->functionSpace_->template createFieldVariable<3>(
      "u", displacementsComponentNames);
  displacements_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "dynamic_hyperelasticity_solver_") +
      "u");
  velocities_ = this->functionSpace_->template createFieldVariable<3>(
      "v", displacementsComponentNames);
  velocities_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "dynamic_hyperelasticity_solver_") +
      "v");
  internalVirtualWork_ = this->functionSpace_->template createFieldVariable<3>(
      "δWint_displacement", displacementsComponentNames);
  internalVirtualWork_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "dynamic_hyperelasticity_solver_") +
      "δWint_displacement");
  accelerationTerm_ = this->functionSpace_->template createFieldVariable<3>(
      "δWint_acceleration", displacementsComponentNames);
  accelerationTerm_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "dynamic_hyperelasticity_solver_") +
      "δWint_acceleration");
  externalVirtualWorkDead_ =
      this->functionSpace_->template createFieldVariable<3>(
          "δWext", displacementsComponentNames);
  externalVirtualWorkDead_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "dynamic_hyperelasticity_solver_") +
      "δWext");
}

//! field variable of u
template <typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<
    FunctionSpaceType>::DisplacementsFieldVariableType>
DynamicHyperelasticitySolver<FunctionSpaceType>::displacements() {
  return this->displacements_;
}

//! field variable of v
template <typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<
    FunctionSpaceType>::DisplacementsFieldVariableType>
DynamicHyperelasticitySolver<FunctionSpaceType>::velocities() {
  return this->velocities_;
}

//! field variable of u_compressible
template <typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<
    FunctionSpaceType>::DisplacementsFieldVariableType>
DynamicHyperelasticitySolver<FunctionSpaceType>::externalVirtualWorkDead() {
  return this->externalVirtualWorkDead_;
}

//! field variable of ∂W_int_compressible
template <typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<
    FunctionSpaceType>::DisplacementsFieldVariableType>
DynamicHyperelasticitySolver<FunctionSpaceType>::internalVirtualWork() {
  return this->internalVirtualWork_;
}

//! field variable of
template <typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<
    FunctionSpaceType>::DisplacementsFieldVariableType>
DynamicHyperelasticitySolver<FunctionSpaceType>::accelerationTerm() {
  return this->accelerationTerm_;
}

template <typename FunctionSpaceType>
typename DynamicHyperelasticitySolver<
    FunctionSpaceType>::FieldVariablesForOutputWriter
DynamicHyperelasticitySolver<
    FunctionSpaceType>::getFieldVariablesForOutputWriter() {
  return std::make_tuple(
      std::shared_ptr<DisplacementsFieldVariableType>(
          std::make_shared<typename FunctionSpaceType::GeometryFieldType>(
              this->functionSpace_->geometryField())), // geometry
      std::shared_ptr<DisplacementsFieldVariableType>(
          this->displacements_), // displacements_
      std::shared_ptr<DisplacementsFieldVariableType>(this->velocities_),
      std::shared_ptr<DisplacementsFieldVariableType>(
          this->internalVirtualWork_),
      std::shared_ptr<DisplacementsFieldVariableType>(this->accelerationTerm_),
      std::shared_ptr<DisplacementsFieldVariableType>(
          this->externalVirtualWorkDead_));
}

template <typename FunctionSpaceType>
typename DynamicHyperelasticitySolver<
    FunctionSpaceType>::FieldVariablesForCheckpointing
DynamicHyperelasticitySolver<
    FunctionSpaceType>::getFieldVariablesForCheckpointing() {
  auto geometryField = std::shared_ptr<DisplacementsFieldVariableType>(
      std::make_shared<typename FunctionSpaceType::GeometryFieldType>(
          this->functionSpace_->geometryField()));
  geometryField->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "dynamic_hyperelasticity_solver_") +
      geometryField->name());

  return std::make_tuple(
      geometryField,
      std::shared_ptr<DisplacementsFieldVariableType>(
          this->displacements_), // displacements_
      std::shared_ptr<DisplacementsFieldVariableType>(this->velocities_),
      std::shared_ptr<DisplacementsFieldVariableType>(
          this->internalVirtualWork_),
      std::shared_ptr<DisplacementsFieldVariableType>(this->accelerationTerm_),
      std::shared_ptr<DisplacementsFieldVariableType>(
          this->externalVirtualWorkDead_));
}
} // namespace Data
