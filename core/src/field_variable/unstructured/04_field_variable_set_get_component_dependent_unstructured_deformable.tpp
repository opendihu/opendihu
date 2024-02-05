#include "field_variable/unstructured/04_field_variable_set_get_component_dependent_unstructured_deformable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <map>
#include <fstream>
#include <iomanip>
#include <cassert>

namespace FieldVariable {

using namespace StringUtility;

//! get a single value from local dof no. for all components
template <int D, typename BasisFunctionType, int nComponents>
std::array<double, nComponents> FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    nComponents>::getValue(node_no_t dofLocalNo) const {
  std::array<double, nComponents> resultVector;
  assert(this->values_);

  for (int componentIndex = 0; componentIndex < nComponents; componentIndex++) {
    // get values and assign them to result values vector
    this->values_->getValues(componentIndex, 1, &dofLocalNo,
                             resultVector.data() + componentIndex);
  }

  return resultVector;
}

//! get the values corresponding to all element-local dofs for all components
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::getElementValues(element_no_t elementNo,
                         std::array<double,
                                    FunctionSpaceType::nDofsPerElement()>
                             &values) const {
  assert(elementNo >= 0 && elementNo < this->functionSpace_->nElementsLocal());
  assert(this->values_);

  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  const std::vector<dof_no_t> &elementDofs =
      this->elementToDofMapping_->getElementDofs(elementNo);

  // VLOG(2) << "getElementValues element " << elementNo << ", nComponents=" <<
  // nComponents;

  this->values_->getValues(0, nDofsPerElement, (PetscInt *)elementDofs.data(),
                           values.data());
}

//! get the values corresponding to all element-local dofs for all components
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::getElementValues(Vc::int_v elementNoLocal,
                         std::array<Vc::double_v,
                                    FunctionSpaceType::nDofsPerElement()>
                             &values) const {
  assert(this->values_);

  const int nVcComponents = Vc::double_v::size();
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  std::array<PetscInt, nDofsPerElement * nVcComponents> indices;
  std::array<Vc::double_v, nDofsPerElement> result;

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
  this->values_->getValues(0, nDofsPerElement * nVcComponents, indices.data(),
                           (double *)&result[nDofsPerElement]);

  // copy result to output values
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++) {
    values[dofIndex] = result[dofIndex];
  }
}

//! get a single value from local dof no. for the single component
template <int D, typename BasisFunctionType>
double FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::getValue(node_no_t dofLocalNo) const {
  assert(this->values_);
  double result;
  this->values_->getValues(0, 1, (PetscInt *)&dofLocalNo, &result);
  return result;
}

//! get a single value from local dof no. for the single component
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::getValues(const std::vector<dof_no_t> &dofLocalNo,
                  std::vector<double> &values) const {
  this->getValues(0, dofLocalNo, values);
}

//! get all stored local values
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::getValuesWithGhosts(std::vector<double> &values,
                            bool onlyNodalValues) const {
  this->getValuesWithGhosts(0, values, onlyNodalValues);
}

//! get all stored local values
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::getValuesWithoutGhosts(std::vector<double> &values,
                               bool onlyNodalValues) const {
  this->getValuesWithoutGhosts(0, values, onlyNodalValues);
}

//! set a single dof (one components) , after all calls to setValue(s),
//! finishGhostManipulation has to be called to apply the cached changes
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::setValue(dof_no_t dofLocalNo, double value,
                 InsertMode petscInsertMode) {
  assert(this->values_);
  this->values_->setValues(0, 1, (PetscInt *)&dofLocalNo, &value,
                           petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e.
  // finishGhostManipulation must be called
}

//! set a single dof (all components), after all calls to setValue(s),
//! finishGhostManipulation has to be called to apply the cached changes
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::setValue(const Vc::int_v &dofLocalNo, const Vc::double_v &value,
                 InsertMode petscInsertMode) {
  assert(this->values_);

  // count number of non-negative indices in dofLocalNo, it is assumed that they
  // occur all before the negative indices
  int nEntries = Vc::double_v::size() - Vc::count(Vc::isnegative(dofLocalNo));
  this->values_->setValues(0, nEntries, (PetscInt *)&dofLocalNo,
                           (double *)&value, petscInsertMode);

  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e.
  // finishGhostManipulation must be called
}

//! set Vc::double_v::size() dofs with the same value
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::setValue(const Vc::int_v &dofLocalNo, double value,
                 InsertMode petscInsertMode) {
  assert(this->values_);
  std::array<double, Vc::double_v::size()> data;
  data.fill(value);

  // count number of non-negative indices in dofLocalNo, it is assumed that they
  // occur all before the negative indices
  int nEntries = Vc::double_v::size() - Vc::count(Vc::isnegative(dofLocalNo));

  this->values_->setValues(0, nEntries, (PetscInt *)&dofLocalNo, data.data(),
                           petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e.
  // finishGhostManipulation must be called
}

//! set the values from a petsc Vec
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::setValues(Vec petscVector) {
  this->setValues(0, petscVector);
}

//! set values for one components for dofs, after all calls to setValue(s),
//! finishGhostManipulation has to be called to apply the cached changes
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::setValues(const std::vector<dof_no_t> &dofNosLocal,
                  const std::vector<double> &values,
                  InsertMode petscInsertMode) {
  assert(this->values_);
  this->values_->setValues(0, dofNosLocal.size(),
                           (PetscInt *)dofNosLocal.data(), values.data(),
                           petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e.
  // finishGhostManipulation must be called
}

//! set values for one components for dofs, after all calls to setValue(s),
//! finishGhostManipulation has to be called to apply the cached changes
template <int D, typename BasisFunctionType>
template <int nValues>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::setValues(const std::array<dof_no_t, nValues> dofNosLocal,
                  const std::array<double, nValues> values,
                  InsertMode petscInsertMode) {
  assert(this->values_);
  this->values_->setValues(0, nValues, (PetscInt *)dofNosLocal.data(),
                           values.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e.
  // finishGhostManipulation must be called
}

//! set values for the single component for all local dofs, after all calls to
//! setValue(s), finishGhostManipulation has to be called to apply the cached
//! changes
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::setValuesWithGhosts(const std::vector<double> &values,
                            InsertMode petscInsertMode) {
  assert(this->values_);
  this->values_->setValues(0, values, petscInsertMode);
}

//! set values for the single component for all local dofs, after all calls to
//! setValue(s), finishGhostManipulation has to be called to apply the cached
//! changes
template <int D, typename BasisFunctionType>
void FieldVariableSetGet<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    1>::setValuesWithoutGhosts(const std::vector<double> &values,
                               InsertMode petscInsertMode) {
  assert(this->values_);
  this->values_->setValues(0, values, petscInsertMode);
}

} // namespace FieldVariable
