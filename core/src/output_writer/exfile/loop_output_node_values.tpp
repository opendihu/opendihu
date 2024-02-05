#include "output_writer/exfile/loop_output_node_values.h"

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
    loopOutputNodeValues(
        const FieldVariablesForOutputWriterType &fieldVariables,
        std::string meshName, std::ostream &stream, node_no_t nodeGlobalNo) {
  // call what to do in the loop body
  if (outputNodeValues<typename std::tuple_element<
          i, FieldVariablesForOutputWriterType>::type>(
          std::get<i>(fieldVariables), meshName, stream, nodeGlobalNo))
    return;

  // advance iteration to next tuple element
  loopOutputNodeValues<FieldVariablesForOutputWriterType, i + 1>(
      fieldVariables, meshName, stream, nodeGlobalNo);
}

// current element is of pointer type (not vector)
template <typename CurrentFieldVariableType>
typename std::enable_if<
    !TypeUtility::isTuple<CurrentFieldVariableType>::value &&
        !TypeUtility::isVector<CurrentFieldVariableType>::value &&
        !Mesh::isComposite<CurrentFieldVariableType>::value,
    bool>::type
outputNodeValues(CurrentFieldVariableType currentFieldVariable,
                 std::string meshName, std::ostream &stream,
                 node_no_t nodeGlobalNo) {
  VLOG(2) << "loop_output_node_values.tpp:34, outputNodeValues, field variable "
          << currentFieldVariable->name();

  // if mesh name is not the specified meshName step over this field variable
  // but do not exit the loop over field variables
  if (currentFieldVariable->functionSpace()->meshName() != meshName) {
    return false; // do not break iteration
  }

  // output the values of a node for the current field variable

  // get the class of the current field variable
  typedef CurrentFieldVariableType FieldVariableType;

  // get the current field variable
  auto fieldVariable =
      currentFieldVariable; // this type is std::shared_ptr<FieldVariable<..>>

  const int nDofsPerNode =
      FieldVariableType::element_type::FunctionSpace::nDofsPerNode();

  // get all dofs that are associated with the current node
  // std::array<dof_no_t,nDofsPerNode> dofGlobalNos;
  std::vector<dof_no_t> dofGlobalNos;
  dofGlobalNos.reserve(nDofsPerNode);
  fieldVariable->functionSpace()->getNodeDofs(nodeGlobalNo, dofGlobalNos);

  // loop over components of the field variable
  for (int componentNo = 0; componentNo < fieldVariable->nComponents();
       componentNo++) {
    // get the values of the dofs for the current component
    // std::array<double,nDofsPerNode> values;
    // fieldVariable->template getValues<nDofsPerNode>(componentNo,
    // dofGlobalNos, values);

    std::vector<double> values;
    values.reserve(nDofsPerNode);

    VLOG(2) << "get dofGlobalNos: " << dofGlobalNos;
    VLOG(2) << "field variable is geometry field: "
            << fieldVariable->isGeometryField();

    fieldVariable->getValues(componentNo, dofGlobalNos, values);

    // output values
    for (double value : values) {
      stream << "  " << std::scientific << std::setprecision(17) << value
             << std::endl;
    }
  }

  return false; // do not break iteration
}

// element i is of vector type
template <typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
outputNodeValues(VectorType currentFieldVariableGradient, std::string meshName,
                 std::ostream &stream, node_no_t nodeGlobalNo) {
  for (auto &currentFieldVariable : currentFieldVariableGradient) {
    // call function on all vector entries
    if (outputNodeValues<typename VectorType::value_type>(
            currentFieldVariable, meshName, stream, nodeGlobalNo))
      return true;
  }

  return false; // do not break iteration
}

// element i is of tuple type
template <typename TupleType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
outputNodeValues(TupleType currentFieldVariableTuple, std::string meshName,
                 std::ostream &stream, node_no_t nodeGlobalNo) {
  // call for tuple element
  loopOutputNodeValues<TupleType>(currentFieldVariableTuple, meshName, stream,
                                  nodeGlobalNo);

  return false; // do not break iteration
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template <typename CurrentFieldVariableType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value,
                        bool>::type
outputNodeValues(CurrentFieldVariableType currentFieldVariable,
                 std::string meshName, std::ostream &stream,
                 node_no_t nodeGlobalNo) {
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
    if (outputNodeValues<std::shared_ptr<SubFieldVariableType>>(
            currentSubFieldVariable, meshName, stream, nodeGlobalNo))
      return true;
  }

  return false; // do not break iteration
}
} // namespace ExfileLoopOverTuple
} // namespace OutputWriter
