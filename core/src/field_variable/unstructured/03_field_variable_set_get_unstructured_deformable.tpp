#include "field_variable/unstructured/03_field_variable_set_get_unstructured_deformable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <map>
#include <fstream>
#include <iomanip>
#include <cassert>

namespace FieldVariable {

using namespace StringUtility;

//! for a specific component, get all values
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    getValuesWithGhosts(int componentNo, std::vector<double> &values,
                        bool onlyNodalValues) const {
  assert(componentNo >= 0 && componentNo < nComponents);

  this->component_[componentNo].getValues(values, onlyNodalValues);
}

//! for a specific component, get all values
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    getValuesWithoutGhosts(int componentNo, std::vector<double> &values,
                           bool onlyNodalValues) const {
  this->getValuesWithGhosts(componentNo, values, onlyNodalValues);
}

//! for a specific component, get all values
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    getValuesWithGhosts(std::vector<std::array<double, nComponents>> &values,
                        bool onlyNodalValues) const {
  std::vector<double> buffer;
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    // get values into buffer
    this->component_[componentNo].getValues(buffer, onlyNodalValues);

    values.resize(buffer.size());

    // copy values from buffer to output vector
    for (int valueIndex = 0; valueIndex < buffer.size(); valueIndex++) {
      values[valueIndex][componentNo] = buffer[valueIndex];
    }
  }
}

//! for a specific component, get all values
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    getValuesWithGhosts(std::array<std::vector<double>, nComponents> &values,
                        bool onlyNodalValues) const {
  std::vector<double> buffer;
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    // get values into buffer
    this->component_[componentNo].getValues(buffer, onlyNodalValues);

    values.resize(buffer.size());

    // copy values from buffer to output vector
    for (int valueIndex = 0; valueIndex < buffer.size(); valueIndex++) {
      values[componentNo][valueIndex] = buffer[valueIndex];
    }
  }
}

//! for a specific component, get all values
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    getValuesWithoutGhosts(std::vector<std::array<double, nComponents>> &values,
                           bool onlyNodalValues) const {
  this->getValuesWithGhosts(values, onlyNodalValues);
}

//! for a specific component, get all values
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    getValuesWithoutGhosts(std::array<std::vector<double>, nComponents> &values,
                           bool onlyNodalValues) const {
  this->getValuesWithGhosts(values, onlyNodalValues);
}

//! get values from their local dof no.s for all components, this eventually
//! does not get all values if there are multiple versions
template <typename FunctionSpaceType, int nComponents>
template <int N>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::getValues(
    std::array<dof_no_t, N> dofLocalNo,
    std::array<std::array<double, nComponents>, N> &values) const {
  assert(this->values_);

  std::array<double, N * nComponents> resultVector;
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    // get values and assign them to result values vector
    this->values_->getValues(componentIndex, N, dofLocalNo.data(),
                             resultVector.data() + componentIndex * N);
  }

  // transform local dof no.s to vector indices of first component
  for (int dofIndex = 0; dofIndex < N; dofIndex++) {
    // copy retrieved values to output array
    for (int componentIndex = 0; componentIndex < nComponents;
         componentIndex++) {
      values[dofIndex][componentIndex] =
          resultVector[componentIndex * N + dofIndex];
    }
  }
}

//! for a specific component, get values from their local dof no.s
template <typename FunctionSpaceType, int nComponents>
template <int N>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::getValues(
    int componentNo, std::array<dof_no_t, N> dofLocalNo,
    std::array<double, N> &values) const {
  assert(componentNo >= 0 && componentNo < nComponents);

  this->component_[componentNo].template getValues<N>(dofLocalNo, values);
}

template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::getValues(
    int componentNo, int nValues, const dof_no_t *dofLocalNo,
    std::vector<double> &values) const {
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);

  int valuesPreviousSize = values.size();
  values.resize(valuesPreviousSize + nValues);
  this->values_->getValues(componentNo, nValues, (PetscInt *)dofLocalNo,
                           values.data() + valuesPreviousSize);
}

//! for a specific component, get values from their local dof no.s
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::getValues(
    int componentNo, const std::vector<dof_no_t> &dofLocalNo,
    std::vector<double> &values) const {
  assert(componentNo >= 0 && componentNo < nComponents);

  this->component_[componentNo].getValues(dofLocalNo, values);
}

template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::getValues(
    const std::vector<dof_no_t> &dofLocalNo,
    std::vector<double> &values) const {
  int nValues = dofLocalNo.size();
  values.resize(nValues * nComponents);

  // get values into values vector
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    std::vector<double> componentValues;
    this->component_[componentIndex].getValues(dofLocalNo, componentValues);
    std::copy(componentValues.begin(), componentValues.end(),
              values.begin() + componentIndex * nValues);
  }
}

//! get values from their local dof no.s for all components
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::getValues(
    std::vector<dof_no_t> dofLocalNo,
    std::vector<std::array<double, nComponents>> &values) const {
  assert(this->values_);
  const int nValues = dofLocalNo.size();
  std::vector<double> result(nValues * nComponents); // temporary result buffer

  int initialSize = values.size();
  values.resize(initialSize + nValues);

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    std::vector<double> componentValues;
    this->component_[componentIndex].getValues(dofLocalNo, componentValues);
    std::copy(componentValues.begin(), componentValues.end(),
              result.begin() + componentIndex * nValues);
  }

  // copy result to output values
  for (int dofIndex = 0; dofIndex < nValues; dofIndex++) {
    for (int componentIndex = 0; componentIndex < nComponents;
         componentIndex++) {
      values[initialSize + dofIndex][componentIndex] =
          result[componentIndex * nValues + dofIndex];
    }
  }
}

//! for a specific component, get a single value from local dof no.
template <typename FunctionSpaceType, int nComponents>
double
FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::getValue(
    int componentNo, node_no_t dofLocalNo) const {
  assert(componentNo >= 0 && componentNo < nComponents);
  return this->component_[componentNo].getValue(dofLocalNo);
}

//! get a single value from local dof no. for all components
template <typename FunctionSpaceType, int nComponents>
std::array<double, nComponents>
FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::getValue(
    node_no_t dofLocalNo) const {
  std::array<double, nComponents> result;
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    result[componentNo] = this->component_[componentNo].getValue(dofLocalNo);
  }
  return result;
}

//! for a specific component, get the values corresponding to all element-local
//! dofs
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    getElementValues(int componentNo, element_no_t elementNo,
                     std::array<double, FunctionSpaceType::nDofsPerElement()>
                         &values) const {
  assert(elementNo >= 0 && elementNo < this->functionSpace_->nElementsLocal());
  assert(componentNo >= 0 && componentNo < nComponents);

  this->component_[componentNo].getElementValues(elementNo, values);
}

//! get the values corresponding to all element-local dofs for all components
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    getElementValues(
        element_no_t elementNo,
        std::array<std::array<double, nComponents>,
                   FunctionSpaceType::nDofsPerElement()> &values) const {
  assert(elementNo >= 0 && elementNo < this->functionSpace_->nElementsLocal());
  assert(this->values_);

  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  const std::vector<dof_no_t> &dofLocalNo =
      this->elementToDofMapping_->getElementDofs(elementNo);

  std::array<double, nDofsPerElement * nComponents> resultVector;
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    // get values and assign them to result values vector
    this->values_->getValues(componentIndex, nDofsPerElement, dofLocalNo.data(),
                             resultVector.data() +
                                 componentIndex * nDofsPerElement);
  }

  // transform local dof no.s to vector indices of first component
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++) {
    // copy retrieved values to output array
    for (int componentIndex = 0; componentIndex < nComponents;
         componentIndex++) {
      values[dofIndex][componentIndex] =
          resultVector[componentIndex * nDofsPerElement + dofIndex];
    }
  }
}

//! vectorized version of getElementValues
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    getElementValues(
        Vc::int_v elementNoLocal,
        std::array<std::array<Vc::double_v, nComponents>,
                   FunctionSpaceType::nDofsPerElement()> &values) const {
  assert(this->functionSpace_);
  assert(this->values_);

  const int nVcComponents = Vc::double_v::size();
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  std::array<PetscInt, nDofsPerElement * nVcComponents> indices;
  std::array<Vc::double_v, nDofsPerElement * nComponents> result;

  VLOG(2) << "getElementValues (vectorized) element " << elementNoLocal
          << ", nComponents=" << nComponents
          << ", nDofsPerElement=" << nDofsPerElement;

  // prepare lookup indices for PETSc vector values_
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++) {
      for (int vcComponent = 0; vcComponent < nVcComponents; vcComponent++) {
        if (elementNoLocal[vcComponent] == -1)
          indices[dofIndex * nVcComponents + vcComponent] =
              0; // set index to 0, then here the dof 0 is retrieved and also
                 // further used in computation but it is discarded later in
                 // setValue
        else
          indices[dofIndex * nVcComponents + vcComponent] =
              this->functionSpace_->getDofNo(elementNoLocal[vcComponent],
                                             dofIndex);
      }
    }

    // get the values for the current component
    this->values_->getValues(
        componentIndex, nDofsPerElement * nVcComponents, indices.data(),
        (double *)&result[componentIndex * nDofsPerElement]);
  }

  // copy result to output values
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++) {
    for (int componentIndex = 0; componentIndex < nComponents;
         componentIndex++) {
      values[dofIndex][componentIndex] =
          result[componentIndex * nDofsPerElement + dofIndex];
    }
  }
}

template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValues(
    int componentNo, Vec petscVector) {
  const int nDofs = this->nDofs();
  std::vector<double> values(nDofs);
  std::vector<PetscInt> indices(nDofs);
  std::iota(indices.begin(), indices.end(), 0);

  PetscErrorCode ierr;
  ierr = VecGetValues(petscVector, nDofs, indices.data(), values.data());
  CHKERRV(ierr);

  this->setValues(componentNo, values);
}

template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValues(
    int componentNo, const std::vector<dof_no_t> &dofNosLocal,
    const std::vector<double> &values, InsertMode petscInsertMode) {
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);
  assert(dofNosLocal.size() == values.size());

  // set the values for the current component
  this->values_->setValues(componentNo, dofNosLocal.size(), dofNosLocal.data(),
                           values.data(), petscInsertMode);
}

template <typename FunctionSpaceType, int nComponents>
template <int N>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValues(
    int componentNo, const std::array<dof_no_t, N> &dofNosLocal,
    const std::array<double, N> &values, InsertMode petscInsertMode) {
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);

  // set the values for the current component
  this->values_->setValues(componentNo, N, dofNosLocal.data(), values.data(),
                           petscInsertMode);
}

//! set values for a given component for given dofs, using raw pointers
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValues(
    int componentNo, int nValues, const dof_no_t *dofNosLocal,
    const double *values, InsertMode petscInsertMode) {
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->values_);

  // set the values for the current component
  this->values_->setValues(componentNo, nValues, dofNosLocal, values,
                           petscInsertMode);
}

template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValues(
    int componentNo,
    std::shared_ptr<FieldVariable<FunctionSpaceType, 1>> fieldVariable) {
  std::vector<double> values;
  fieldVariable->getValuesWithoutGhosts(0, values, false);
  this->setValues(componentNo, values);
}

template <typename FunctionSpaceType, int nComponents>
template <int N>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValues(
    const std::array<dof_no_t, N> &dofNosLocal,
    const std::array<std::array<double, nComponents>, N> &values,
    InsertMode petscInsertMode) {
  assert(this->values_);

  std::array<double, N> valuesBuffer;

  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    // loop over dofs and prepare values of current component
    for (int dofIndex = 0; dofIndex < N; dofIndex++) {
      valuesBuffer[dofIndex] = values[dofIndex][componentIndex];
    }

    // set the values for the current component
    this->values_->setValues(componentIndex, N, dofNosLocal.data(),
                             valuesBuffer.data(), petscInsertMode);
  }
}

template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    extractComponentCopy(int componentNo,
                         std::shared_ptr<FieldVariable<FunctionSpaceType, 1>>
                             extractedFieldVariable) {
  std::vector<double> values;
  this->getValuesWithoutGhosts(componentNo, values, false);
  extractedFieldVariable->setValuesWithoutGhosts(values, INSERT_VALUES);
}

template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    extractComponentShared(int componentNo,
                           std::shared_ptr<FieldVariable<FunctionSpaceType, 1>>
                               extractedFieldVariable) {
  extractComponentCopy(componentNo, extractedFieldVariable);
}

//! copy the values from another field variable of the same type
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValues(
    FieldVariable<FunctionSpaceType, nComponents> &rhs) {
  this->values_->setValues(*rhs.partitionedPetscVec());
}

template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValues(
    const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values,
    InsertMode petscInsertMode) {
  assert(this->nComponents == 1 && nComponents == 1);
  const int nValues = values.size();
  assert(this->values_);

  this->values_->setValues(0, nValues, (const int *)dofNosLocal.data(),
                           values.data(), petscInsertMode);

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e.
  // finishGhostManipulation must be called
}

//! set value for all dofs
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValues(
    double value) {
  // get number of dofs
  assert(this->functionSpace_);
  const dof_no_t nDofs = this->functionSpace_->nDofsLocalWithoutGhosts();

  std::vector<double> valueBuffer(nDofs, value);

  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    this->setValuesWithoutGhosts(componentIndex, valueBuffer, INSERT_VALUES);
  }
}

//! set values for all components for dofs, after all calls to setValue(s),
//! finishGhostManipulation has to be called to apply the cached changes
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValues(
    const std::vector<dof_no_t> &dofNosLocal,
    const std::vector<std::array<double, nComponents>> &values,
    InsertMode petscInsertMode) {
  assert(dofNosLocal.size() == values.size());
  assert(this->values_);

  const int nValues = values.size();
  std::vector<double> valuesBuffer(nValues);

  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    // loop over dofs and prepare values of current component
    for (int dofIndex = 0; dofIndex < nValues; dofIndex++) {
      valuesBuffer[dofIndex] = values[dofIndex][componentIndex];
    }

    this->values_->setValues(componentIndex, nValues, dofNosLocal.data(),
                             valuesBuffer.data(), petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e.
  // finishGhostManipulation must be called
}

//! set a single dof (all components) , after all calls to setValue(s),
//! finishGhostManipulation has to be called to apply the cached changes
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValue(
    dof_no_t dofLocalNo, const std::array<double, nComponents> &value,
    InsertMode petscInsertMode) {
  assert(this->values_);

  // loop over components and set single value for each component
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    LOG(DEBUG) << "set value of \"" << this->name_ << "\" for component "
               << componentIndex << " to " << value[componentIndex];
    this->values_->setValues(componentIndex, 1, &dofLocalNo,
                             value.data() + componentIndex, petscInsertMode);
  }

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e.
  // finishGhostManipulation must be called
}

template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValue(
    Vc::int_v dofLocalNo, const std::array<Vc::double_v, nComponents> &value,
    InsertMode petscInsertMode) {
  // loop over components and set vectorized value for each component
  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    this->setValue(componentIndex, dofLocalNo, value[componentIndex],
                   petscInsertMode);
  }
}

template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValue(
    int componentNo, dof_no_t dofLocalNo, double value,
    InsertMode petscInsertMode) {
  assert(this->values_);

  this->values_->setValues(componentNo, 1, &dofLocalNo, &value,
                           petscInsertMode);
}

//! set a given component of Vc::double_v::size() dofs with the vectorized
//! value, after all calls to setValue(s), finishGhostManipulation has to be
//! called to apply the cached changes
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValue(
    int componentNo, Vc::int_v dofLocalNo, Vc::double_v value,
    InsertMode petscInsertMode) {
  assert(this->values_);

  // count number of non-negative indices in dofLocalNo, it is assumed that they
  // occur all before the negative indices
  int nEntries = Vc::double_v::size() - Vc::count(Vc::isnegative(dofLocalNo));

  this->values_->setValues(componentNo, nEntries, (PetscInt *)&dofLocalNo,
                           (double *)&value, petscInsertMode);
}

//! set a given component of Vc::double_v::size() dofs with the same value
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::setValue(
    int componentNo, Vc::int_v dofLocalNo, double value,
    InsertMode petscInsertMode) {
  assert(this->values_);
  /*  std::array<double,Vc::double_v::size()> data;
    data.fill(value);

    // store Vc vectors in order to get the raw memory
    std::array<int,Vc::double_v::size()> indices;
    dofLocalNo.store(indices.data());
  */
  // count number of non-negative indices in dofLocalNo, it is assumed that they
  // occur all before the negative indices
  int nEntries = Vc::double_v::size() - Vc::count(Vc::isnegative(dofLocalNo));

  this->values_->setValues(componentNo, nEntries, (PetscInt *)&dofLocalNo,
                           (double *)&value, petscInsertMode);
}

//! set values for the specified component for all local dofs, after all calls
//! to setValue(s), finishGhostManipulation has to be called to apply the cached
//! changes
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    setValuesWithGhosts(int componentNo, const std::vector<double> &values,
                        InsertMode petscInsertMode) {
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(values.size() ==
         this->functionSpace_->meshPartition()->nDofsLocalWithGhosts());
  assert(this->values_);

  this->values_->setValues(
      componentNo, values.size(),
      this->functionSpace_->meshPartition()->dofNosLocal().data(),
      values.data(), petscInsertMode);

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e.
  // finishGhostManipulation must be called
}

//! set values for the specified component for all local dofs, after all calls
//! to setValue(s), finishGhostManipulation has to be called to apply the cached
//! changes
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    setValuesWithoutGhosts(int componentNo, const std::vector<double> &values,
                           InsertMode petscInsertMode) {
  assert(componentNo >= 0 && componentNo < nComponents);
  assert(this->functionSpace_);
  assert(this->functionSpace_->meshPartition());
  assert(values.size() ==
         this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());
  assert(this->values_);

  // set the values, this is the same call as setValuesWithGhosts, but the
  // number of values is smaller and therefore the last dofs which are the
  // ghosts are not touched
  this->values_->setValues(
      componentNo, values.size(),
      this->functionSpace_->meshPartition()->dofNosLocal().data(),
      values.data(), petscInsertMode);
}

//! set values for all components for all local dofs, after all calls to
//! setValue(s), finishGhostManipulation has to be called to apply the cached
//! changes
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    setValuesWithGhosts(
        const std::vector<std::array<double, nComponents>> &values,
        InsertMode petscInsertMode) {
  assert(values.size() ==
         this->functionSpace_->meshPartition()->nDofsLocalWithGhosts());

  this->setValues(this->functionSpace_->meshPartition()->dofNosLocal(), values,
                  petscInsertMode);
}

//! set values for all components for all local dofs, after all calls to
//! setValue(s), finishGhostManipulation has to be called to apply the cached
//! changes
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType, nComponents>::
    setValuesWithoutGhosts(
        const std::vector<std::array<double, nComponents>> &values,
        InsertMode petscInsertMode) {
  assert(values.size() ==
         this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());

  this->setValues(this->functionSpace_->meshPartition()->dofNosLocal(), values,
                  petscInsertMode);
}

//! set value to zero for all dofs
template <typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetUnstructured<FunctionSpaceType,
                                     nComponents>::zeroEntries() {
  assert(this->values_);
  this->values_->zeroEntries();
}

} // namespace FieldVariable
