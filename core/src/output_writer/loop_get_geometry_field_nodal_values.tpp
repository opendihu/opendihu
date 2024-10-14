#include <cstdlib>
#include "field_variable/field_variable.h"

namespace OutputWriter {

namespace LoopOverTuple {

/** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template <typename FieldVariablesForOutputWriterType, int i>
    inline typename std::enable_if <
    i<std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
    loopGetGeometryFieldNodalValues(
        const FieldVariablesForOutputWriterType &fieldVariables,
        std::set<std::string> meshNames, std::vector<double> &values) {
  // call what to do in the loop body
  if (getGeometryFieldNodalValues<
          typename std::tuple_element<i,
                                      FieldVariablesForOutputWriterType>::type,
          FieldVariablesForOutputWriterType>(std::get<i>(fieldVariables),
                                             fieldVariables, meshNames, values))
    return;

  // advance iteration to next tuple element
  loopGetGeometryFieldNodalValues<FieldVariablesForOutputWriterType, i + 1>(
      fieldVariables, meshNames, values);
}

// current element is of pointer type (not vector)
template <typename CurrentFieldVariableType,
          typename FieldVariablesForOutputWriterType>
typename std::enable_if<
    !TypeUtility::isTuple<CurrentFieldVariableType>::value &&
        !TypeUtility::isVector<CurrentFieldVariableType>::value &&
        !Mesh::isComposite<CurrentFieldVariableType>::value,
    bool>::type
getGeometryFieldNodalValues(
    CurrentFieldVariableType currentFieldVariable,
    const FieldVariablesForOutputWriterType &fieldVariables,
    std::set<std::string> meshNames, std::vector<double> &values) {
  VLOG(1) << "getGeometryFieldNodalValues meshNames: " << meshNames
          << ", own: " << currentFieldVariable->functionSpace()->meshName()
          << ", fieldVariable name \"" << currentFieldVariable->name() << "\""
          << ", isGeometryField: " << currentFieldVariable->isGeometryField();

  // if mesh name is one of the specified meshNames and it is the geometry field
  if (meshNames.find(currentFieldVariable->functionSpace()->meshName()) !=
          meshNames.end() &&
      currentFieldVariable->isGeometryField()) {
    const int nComponents =
        CurrentFieldVariableType::element_type::nComponents();
    std::array<std::vector<double>, nComponents> componentValues;

    // ensure that ghost values are in place
    Partition::values_representation_t old_representation =
        currentFieldVariable->currentRepresentation();
    currentFieldVariable->zeroGhostBuffer();
    currentFieldVariable->setRepresentationGlobal();
    currentFieldVariable->startGhostManipulation();

    VLOG(1) << "nComponents: " << nComponents;

    // get all local values with ghosts for the components
    for (int componentNo = 0; componentNo < nComponents; componentNo++) {
      std::vector<double> retrievedLocalValues;
      currentFieldVariable->getValues(componentNo,
                                      currentFieldVariable->functionSpace()
                                          ->meshPartition()
                                          ->dofNosLocalNaturalOrdering(),
                                      retrievedLocalValues);

      const int nDofsPerNode =
          CurrentFieldVariableType::element_type::FunctionSpace::nDofsPerNode();
      const node_no_t nNodesLocal = currentFieldVariable->functionSpace()
                                        ->meshPartition()
                                        ->nNodesLocalWithGhosts();

      VLOG(1) << "nNodesLocal: " << nNodesLocal;

      // for Hermite only extract the non-derivative values
      componentValues[componentNo].resize(nNodesLocal);

      int index = 0;
      for (int i = 0; i < nNodesLocal; i++) {
        componentValues[componentNo][i] = retrievedLocalValues[index];
        index += nDofsPerNode;
      }
    }

    // reset variable to old representation, so external code is not suprised
    // eg. time stepping code usally uses representationContiguous and will be
    // suprised if this changed we did not write to the values
    currentFieldVariable->setRepresentation(
        old_representation, values_modified_t::values_unchanged);

    // copy values in consecutive order (x y z x y z) to output
    values.reserve(values.size() + componentValues[0].size() * nComponents);
    for (int i = 0; i < componentValues[0].size(); i++) {
      for (int componentNo = 0; componentNo < nComponents; componentNo++) {
        values.push_back(componentValues[componentNo][i]);
      }
    }
  }

  return false; // do not break iteration
}

// element i is of vector type
template <typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
getGeometryFieldNodalValues(
    VectorType currentFieldVariableGradient,
    const FieldVariablesForOutputWriterType &fieldVariables,
    std::set<std::string> meshNames, std::vector<double> &values) {
  for (auto &currentFieldVariable : currentFieldVariableGradient) {
    // call function on all vector entries
    if (getGeometryFieldNodalValues<typename VectorType::value_type,
                                    FieldVariablesForOutputWriterType>(
            currentFieldVariable, fieldVariables, meshNames, values))
      return true; // break iteration
  }
  return false; // do not break iteration
}

// element i is of tuple type
template <typename TupleType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
getGeometryFieldNodalValues(
    TupleType currentFieldVariableTuple,
    const FieldVariablesForOutputWriterType &fieldVariables,
    std::set<std::string> meshNames, std::vector<double> &values) {
  // call for tuple element
  loopGetGeometryFieldNodalValues<TupleType>(currentFieldVariableTuple,
                                             meshNames, values);

  return false; // do not break iteration
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template <typename CurrentFieldVariableType,
          typename FieldVariablesForOutputWriterType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value,
                        bool>::type
getGeometryFieldNodalValues(
    CurrentFieldVariableType currentFieldVariable,
    const FieldVariablesForOutputWriterType &fieldVariables,
    std::set<std::string> meshNames, std::vector<double> &values) {
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
    if (getGeometryFieldNodalValues<std::shared_ptr<SubFieldVariableType>,
                                    FieldVariablesForOutputWriterType>(
            currentSubFieldVariable, fieldVariables, meshNames, values))
      return true;
  }

  return false; // do not break iteration
}
} // namespace LoopOverTuple
} // namespace OutputWriter
