#include "output_writer/hdf5/loop_output.h"

#include "output_writer/loop_count_n_field_variables_of_mesh.h"
#include "output_writer/hdf5/hdf5_writer.h"

#include <cstdlib>
#include "field_variable/field_variable.h"

namespace OutputWriter {

namespace HDF5LoopOverTuple {

/** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template <typename FieldVariablesForOutputWriterType,
          typename AllFieldVariablesForOutputWriterType, int i>
    inline typename std::enable_if <
    i<std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
    loopOutput(hid_t fileID,
               const FieldVariablesForOutputWriterType &fieldVariables,
               const AllFieldVariablesForOutputWriterType &allFieldVariables,
               const std::string &meshName,
               const PythonConfig &specificSettings, double currentTime) {
  // call what to do in the loop body
  if (output<typename std::tuple_element<
                 i, FieldVariablesForOutputWriterType>::type,
             AllFieldVariablesForOutputWriterType>(
          fileID, std::get<i>(fieldVariables), allFieldVariables, meshName,
          specificSettings, currentTime))
    return;

  // advance iteration to next tuple element
  loopOutput<FieldVariablesForOutputWriterType,
             AllFieldVariablesForOutputWriterType, i + 1>(
      fileID, fieldVariables, allFieldVariables, meshName, specificSettings,
      currentTime);
}

// current element is of pointer type (not vector)
template <typename CurrentFieldVariableType,
          typename FieldVariablesForOutputWriterType>
typename std::enable_if<
    !TypeUtility::isTuple<CurrentFieldVariableType>::value &&
        !TypeUtility::isVector<CurrentFieldVariableType>::value &&
        !Mesh::isComposite<CurrentFieldVariableType>::value,
    bool>::type
output(hid_t fileID, CurrentFieldVariableType currentFieldVariable,
       const FieldVariablesForOutputWriterType &fieldVariables,
       const std::string &meshName, const PythonConfig &specificSettings,
       double currentTime) {
  // if mesh name is the specified meshName
  if (currentFieldVariable->functionSpace()->meshName() == meshName) {
    // here we have the type of the mesh with meshName (which is typedef to
    // FunctionSpace)
    typedef typename CurrentFieldVariableType::element_type::FunctionSpace
        FunctionSpace;

    int nFieldVariablesInMesh = 0;
    LoopOverTuple::loopCountNFieldVariablesOfMesh(fieldVariables, meshName,
                                                  nFieldVariablesInMesh);

    // call exfile writer to output all field variables with the meshName
    HDF5Writer<FunctionSpace, FieldVariablesForOutputWriterType>::outputFile(
        fileID, fieldVariables, meshName, currentFieldVariable->functionSpace(),
        nFieldVariablesInMesh, specificSettings, currentTime);

    return true; // break iteration
  }

  return false; // do not break iteration
}

// element i is of vector type
template <typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
output(hid_t fileID, VectorType currentFieldVariableGradient,
       const FieldVariablesForOutputWriterType &fieldVariables,
       const std::string &meshName, const PythonConfig &specificSettings,
       double currentTime) {
  for (auto &currentFieldVariable : currentFieldVariableGradient) {
    // call function on all vector entries
    if (output<typename VectorType::value_type,
               FieldVariablesForOutputWriterType>(
            fileID, currentFieldVariable, fieldVariables, meshName,
            specificSettings, currentTime))
      return true; // break iteration
  }
  return false; // do not break iteration
}

// element i is of tuple type
template <typename TupleType, typename AllFieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
output(hid_t fileID, TupleType currentFieldVariableTuple,
       const AllFieldVariablesForOutputWriterType &fieldVariables,
       const std::string &meshName, const PythonConfig &specificSettings,
       double currentTime) {
  // call for tuple element
  loopOutput<TupleType, AllFieldVariablesForOutputWriterType>(
      fileID, currentFieldVariableTuple, fieldVariables, meshName,
      specificSettings, currentTime);

  return false; // do not break iteration
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template <typename CurrentFieldVariableType,
          typename AllFieldVariablesForOutputWriterType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value,
                        bool>::type
output(hid_t fileID, CurrentFieldVariableType currentFieldVariable,
       const AllFieldVariablesForOutputWriterType &fieldVariables,
       const std::string &meshName, const PythonConfig &specificSettings,
       double currentTime) {
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
    if (output<std::shared_ptr<SubFieldVariableType>,
               AllFieldVariablesForOutputWriterType>(
            fileID, currentSubFieldVariable, fieldVariables, meshName,
            specificSettings, currentTime))
      return true;
  }

  return false; // do not break iteration
}
} // namespace HDF5LoopOverTuple
} // namespace OutputWriter
