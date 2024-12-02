#include "data_management/specialized_solver/quasi_static_nonlinear_elasticity_febio.h"

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data {

QuasiStaticNonlinearElasticityFebio::QuasiStaticNonlinearElasticityFebio(
    DihuContext context)
    : Data<FunctionSpace>::Data(context) {}

void QuasiStaticNonlinearElasticityFebio::initialize() {
  // call initialize of base class
  Data<FunctionSpace>::initialize();

  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();
  slotConnectorData_->addFieldVariable(activation());

  // parse slot names for all slot connector data slots, only one slot here
  this->context_.getPythonConfig().getOptionVector(
      "slotNames", slotConnectorData_->slotNames);

  // make sure that there are as many slot names as slots
  slotConnectorData_->slotNames.resize(slotConnectorData_->nSlots());
}

bool QuasiStaticNonlinearElasticityFebio::restoreState(
    const InputReader::Generic &r) {
  std::vector<double> activation, displacements, reactionForce, cauchyStress,
      pk2Stress, greenLagrangeStrain, relativeVolume;
  if (!r.readDoubleVector(this->activation_->uniqueName().c_str(),
                          activation)) {
    return false;
  }
  if (!r.readDoubleVector(this->displacements_->uniqueName().c_str(),
                          displacements)) {
    return false;
  }
  if (!r.readDoubleVector(this->reactionForce_->uniqueName().c_str(),
                          reactionForce)) {
    return false;
  }
  if (!r.readDoubleVector(this->cauchyStress_->uniqueName().c_str(),
                          cauchyStress)) {
    return false;
  }
  if (!r.readDoubleVector(this->pk2Stress_->uniqueName().c_str(), pk2Stress)) {
    return false;
  }
  if (!r.readDoubleVector(this->greenLagrangeStrain_->uniqueName().c_str(),
                          greenLagrangeStrain)) {
    return false;
  }
  if (!r.readDoubleVector(this->relativeVolume_->uniqueName().c_str(),
                          relativeVolume)) {
    return false;
  }

  std::array<std::vector<double>, 3> geometryValues;
  if (!r.template readDoubleVecD<3>(
          this->functionSpace_->geometryField().name().c_str(), geometryValues,
          "3D/")) {
    return false;
  }

  this->activation_->setValues(activation);
  this->displacements_->setValues(displacements);
  this->reactionForce_->setValues(reactionForce);
  this->cauchyStress_->setValues(cauchyStress);
  this->pk2Stress_->setValues(pk2Stress);
  this->greenLagrangeStrain_->setValues(greenLagrangeStrain);
  this->relativeVolume_->setValues(relativeVolume);

  // for (size_t i = 0; i < 3; i++) {
  //   this->functionSpace_->geometryField().setValuesWithGhosts(
  //       i, geometryValues[i], INSERT_VALUES);
  // }

  // Note we do not need to hande referenceGeometry here
  return true;
}

void QuasiStaticNonlinearElasticityFebio::createPetscObjects() {
  LOG(DEBUG) << "QuasiStaticNonlinearElasticityFebio::createPetscObjects";

  assert(this->functionSpace_);

  activation_ =
      this->functionSpace_->template createFieldVariable<1>("activation");
  activation_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "quasi_static_nonlinear_elasticity_febio_") +
      "activation");
  displacements_ = this->functionSpace_->template createFieldVariable<3>("u");
  displacements_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "quasi_static_nonlinear_elasticity_febio_") +
      "u");
  reactionForce_ =
      this->functionSpace_->template createFieldVariable<3>("reactionForce");
  reactionForce_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "quasi_static_nonlinear_elasticity_febio_") +
      "reactionForce");
  cauchyStress_ = this->functionSpace_->template createFieldVariable<6>(
      "sigma (Cauchy stress)");
  cauchyStress_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "quasi_static_nonlinear_elasticity_febio_") +
      "sigma (Cauchy stress)");
  pk2Stress_ =
      this->functionSpace_->template createFieldVariable<6>("S (Pk2 stress)");
  pk2Stress_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "quasi_static_nonlinear_elasticity_febio_") +
      "S (Pk2 stress)");
  greenLagrangeStrain_ = this->functionSpace_->template createFieldVariable<6>(
      "E (Green-Lagrange strain)");
  greenLagrangeStrain_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "quasi_static_nonlinear_elasticity_febio_") +
      "E (Green-Lagrange strain)");
  relativeVolume_ = this->functionSpace_->template createFieldVariable<1>("J");
  relativeVolume_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "quasi_static_nonlinear_elasticity_febio_") +
      "J");

  // copy initial geometry to referenceGeometry
  referenceGeometry_ = std::make_shared<FieldVariableTypeVector>(
      this->functionSpace_->geometryField(), "referenceGeometry");
  referenceGeometry_->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "quasi_static_nonlinear_elasticity_febio_") +
      "referenceGeometry");

  LOG(DEBUG) << "pointer referenceGeometry: "
             << referenceGeometry_->partitionedPetscVec();
  LOG(DEBUG) << "pointer geometryField: "
             << this->functionSpace_->geometryField().partitionedPetscVec();
}

std::shared_ptr<typename QuasiStaticNonlinearElasticityFebio::FieldVariableType>
QuasiStaticNonlinearElasticityFebio::activation() {
  return this->activation_;
}

//! return the field variable
std::shared_ptr<
    typename QuasiStaticNonlinearElasticityFebio::FieldVariableTypeVector>
QuasiStaticNonlinearElasticityFebio::referenceGeometry() {
  return this->referenceGeometry_;
}

//! return the field variable
std::shared_ptr<
    typename QuasiStaticNonlinearElasticityFebio::FieldVariableTypeVector>
QuasiStaticNonlinearElasticityFebio::displacements() {
  return this->displacements_;
}

//! return the field variable
std::shared_ptr<
    typename QuasiStaticNonlinearElasticityFebio::FieldVariableTypeVector>
QuasiStaticNonlinearElasticityFebio::reactionForce() {
  return this->reactionForce_;
}

//! return the field variable
std::shared_ptr<
    typename QuasiStaticNonlinearElasticityFebio::FieldVariableTypeTensor>
QuasiStaticNonlinearElasticityFebio::cauchyStress() {
  return this->cauchyStress_;
}

//! return the field variable
std::shared_ptr<
    typename QuasiStaticNonlinearElasticityFebio::FieldVariableTypeTensor>
QuasiStaticNonlinearElasticityFebio::pk2Stress() {
  return this->pk2Stress_;
}

//! return the field variable
std::shared_ptr<
    typename QuasiStaticNonlinearElasticityFebio::FieldVariableTypeTensor>
QuasiStaticNonlinearElasticityFebio::greenLagrangeStrain() {
  return this->greenLagrangeStrain_;
}

//! return the field variable
std::shared_ptr<typename QuasiStaticNonlinearElasticityFebio::FieldVariableType>
QuasiStaticNonlinearElasticityFebio::relativeVolume() {
  return this->relativeVolume_;
}

void QuasiStaticNonlinearElasticityFebio::print() {}

std::shared_ptr<
    typename QuasiStaticNonlinearElasticityFebio::SlotConnectorDataType>
QuasiStaticNonlinearElasticityFebio::getSlotConnectorData() {
  return slotConnectorData_;
}

typename QuasiStaticNonlinearElasticityFebio::FieldVariablesForOutputWriter
QuasiStaticNonlinearElasticityFebio::getFieldVariablesForOutputWriter() {
  // these field variables will be written to output files

  std::shared_ptr<FieldVariableTypeVector> geometryField =
      std::make_shared<FieldVariableTypeVector>(
          this->functionSpace_->geometryField());

  return std::tuple_cat(
      std::tuple<std::shared_ptr<FieldVariableTypeVector>>(geometryField),
      std::tuple<std::shared_ptr<FieldVariableType>>(this->activation_),
      std::tuple<std::shared_ptr<FieldVariableTypeVector>>(
          this->displacements_),
      std::tuple<std::shared_ptr<FieldVariableTypeVector>>(
          this->reactionForce_),
      std::tuple<std::shared_ptr<FieldVariableTypeTensor>>(this->cauchyStress_),
      std::tuple<std::shared_ptr<FieldVariableTypeTensor>>(this->pk2Stress_),
      std::tuple<std::shared_ptr<FieldVariableTypeTensor>>(
          this->greenLagrangeStrain_),
      std::tuple<std::shared_ptr<FieldVariableType>>(this->relativeVolume_));
}

typename QuasiStaticNonlinearElasticityFebio::FieldVariablesForCheckpointing
QuasiStaticNonlinearElasticityFebio::getFieldVariablesForCheckpointing() {
  std::shared_ptr<FieldVariableTypeVector> geometryField =
      std::make_shared<FieldVariableTypeVector>(
          this->functionSpace_->geometryField());
  geometryField->setUniqueName(
      StringUtility::getFirstNE(this->uniquePrefix_,
                                "quasi_static_nonlinear_elasticity_febio_") +
      geometryField->name());

  return std::tuple_cat(
      std::tuple<std::shared_ptr<FieldVariableTypeVector>>(geometryField),
      std::tuple<std::shared_ptr<FieldVariableType>>(this->activation_),
      std::tuple<std::shared_ptr<FieldVariableTypeVector>>(
          this->displacements_),
      std::tuple<std::shared_ptr<FieldVariableTypeVector>>(
          this->reactionForce_),
      std::tuple<std::shared_ptr<FieldVariableTypeTensor>>(this->cauchyStress_),
      std::tuple<std::shared_ptr<FieldVariableTypeTensor>>(this->pk2Stress_),
      std::tuple<std::shared_ptr<FieldVariableTypeTensor>>(
          this->greenLagrangeStrain_),
      std::tuple<std::shared_ptr<FieldVariableType>>(this->relativeVolume_));
}
} // namespace Data
