#include "output_writer/exfile/loop_check_if_new_exelem_header_necessary.h"

#include <cstdlib>
#include "field_variable/field_variable.h"

namespace OutputWriter {

namespace ExfileLoopOverTuple {

/** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template <typename FieldVariablesForOutputWriterType, int i>
    inline typename std::enable_if <
    i<std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
    loopCheckIfNewExelemHeaderNecessary(
        const FieldVariablesForOutputWriterType &fieldVariables,
        std::string meshName, element_no_t currentFieldVariableGlobalNo,
        bool &newHeaderNecessary) {
  // call what to do in the loop body
  if (checkIfNewExelemHeaderNecessary<typename std::tuple_element<
          i, FieldVariablesForOutputWriterType>::type>(
          std::get<i>(fieldVariables), meshName, currentFieldVariableGlobalNo,
          newHeaderNecessary))
    return;

  // advance iteration to next tuple element
  loopCheckIfNewExelemHeaderNecessary<FieldVariablesForOutputWriterType, i + 1>(
      fieldVariables, meshName, currentFieldVariableGlobalNo,
      newHeaderNecessary);
}

// current element is of pointer type (not vector)
template <typename CurrentFieldVariableType>
typename std::enable_if<
    !TypeUtility::isTuple<CurrentFieldVariableType>::value &&
        !TypeUtility::isVector<CurrentFieldVariableType>::value &&
        !Mesh::isComposite<CurrentFieldVariableType>::value,
    bool>::type
checkIfNewExelemHeaderNecessary(CurrentFieldVariableType currentFieldVariable,
                                std::string meshName,
                                element_no_t currentFieldVariableGlobalNo,
                                bool &newHeaderNecessary) {
  // if mesh name is not the specified meshName step over this field variable
  // but do not exit the loop over field variables
  if (currentFieldVariable->functionSpace()->meshName() != meshName) {
    return false; // do not break iteration
  }

  VLOG(2) << "check if field variable " << currentFieldVariable->name()
          << " has the same exfileRepr for elements "
          << currentFieldVariableGlobalNo - 1 << " and "
          << currentFieldVariableGlobalNo;

  if (!currentFieldVariable->haveSameExfileRepresentation(
          currentFieldVariableGlobalNo - 1, currentFieldVariableGlobalNo)) {
    newHeaderNecessary = true;

    // break iteration
    return true;
  }
  return false; // do not break iteration
}

// element i is of vector type
template <typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
checkIfNewExelemHeaderNecessary(VectorType currentFieldVariableGradient,
                                std::string meshName,
                                element_no_t currentFieldVariableGlobalNo,
                                bool &newHeaderNecessary) {
  for (auto &currentFieldVariable : currentFieldVariableGradient) {
    // call function on all vector entries
    if (checkIfNewExelemHeaderNecessary<typename VectorType::value_type>(
            currentFieldVariable, meshName, currentFieldVariableGlobalNo,
            newHeaderNecessary))
      return true;
  }

  return false; // do not break iteration
}

// element i is of vector type
template <typename TupleType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
checkIfNewExelemHeaderNecessary(TupleType currentFieldVariableTuple,
                                std::string meshName,
                                element_no_t currentFieldVariableGlobalNo,
                                bool &newHeaderNecessary) {
  // call for tuple element
  loopCheckIfNewExelemHeaderNecessary<TupleType>(
      currentFieldVariableTuple, meshName, currentFieldVariableGlobalNo,
      newHeaderNecessary);

  return false; // do not break iteration
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template <typename CurrentFieldVariableType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value,
                        bool>::type
checkIfNewExelemHeaderNecessary(CurrentFieldVariableType currentFieldVariable,
                                std::string meshName,
                                element_no_t currentFieldVariableGlobalNo,
                                bool &newHeaderNecessary) {
  const int D = CurrentFieldVariableType::element_type::FunctionSpace::dim();
  typedef typename CurrentFieldVariableType::element_type::FunctionSpace::
      BasisFunction BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,
                                       BasisFunctionType>
      SubFunctionSpaceType;
  const int nComponents = CurrentFieldVariableType::element_type::nComponents();

  typedef FieldVariable::FieldVariable<SubFunctionSpaceType, nComponents>
      SubFieldVariableType;

  std::vector<std::shared_ptr<SubFieldVariableType>> subFieldVariables;
  currentFieldVariable->getSubFieldVariables(subFieldVariables);

  for (auto &currentSubFieldVariable : subFieldVariables) {
    // call function on all vector entries
    if (checkIfNewExelemHeaderNecessary<std::shared_ptr<SubFieldVariableType>>(
            currentSubFieldVariable, meshName, currentFieldVariableGlobalNo,
            newHeaderNecessary))
      return true;
  }

  return false; // do not break iteration
}
} // namespace ExfileLoopOverTuple
} // namespace OutputWriter
