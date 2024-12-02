#pragma once

#include <iostream>
#include <vector>

#include "control/types.h"
#include "output_writer/generic.h"
#include "output_writer/hdf5/hdf5.h"
#include "mesh/structured_regular_fixed.h"

namespace OutputWriter {

/** Base class of HDF5Writer that writes vtk files of given field variables.
 *  FieldVariablesForOutputWriterType is a std::tuple<std::shared_ptr<>,
 * std::shared_ptr<>, ...> of field variables. Only field variables which are
 * defined on the specified mesh will be output. The FunctionSpaceType has to be
 * the type of the field variables given in meshName.
 */
template <typename FunctionSpaceType,
          typename FieldVariablesForOutputWriterType>
class HDF5Writer {
public:
  //! write paraview file to given fileID, only output fieldVariables that are
  //! on a mesh with the given meshName
  static void outputFile(HDF5Utils::Group &group,
                         FieldVariablesForOutputWriterType fieldVariables,
                         const std::string &meshName,
                         std::shared_ptr<FunctionSpaceType> mesh,
                         int nFieldVariablesOfMesh,
                         const PythonConfig &specificSettings) {}
};

/** Partial specialization for regular fixed mesh.
 *  Outputs a rectilinear grid.
 */
template <int D, typename BasisFunctionType,
          typename FieldVariablesForOutputWriterType>
class HDF5Writer<
    FunctionSpace::FunctionSpace<::Mesh::StructuredRegularFixedOfDimension<D>,
                                 BasisFunctionType>,
    FieldVariablesForOutputWriterType> {
public:
  //! write paraview file to given fileID, only output fieldVariables that are
  //! on a mesh with the given meshName
  static void outputFile(
      HDF5Utils::Group &group, FieldVariablesForOutputWriterType fieldVariables,
      const std::string &meshName,
      std::shared_ptr<FunctionSpace::FunctionSpace<
          ::Mesh::StructuredRegularFixedOfDimension<D>, BasisFunctionType>>
          mesh,
      int nFieldVariablesOfMesh, const PythonConfig &specificSettings,
      double currentTime);
};

/** Partial specialization for structured mesh.
 *  Outputs a structured grid.
 */
template <int D, typename BasisFunctionType,
          typename FieldVariablesForOutputWriterType>
class HDF5Writer<
    FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    FieldVariablesForOutputWriterType> {
public:
  //! write paraview file to given fileID, only output fieldVariables that are
  //! on a mesh with the given meshName
  static void
  outputFile(HDF5Utils::Group &group,
             FieldVariablesForOutputWriterType fieldVariables,
             const std::string &meshName,
             std::shared_ptr<FunctionSpace::FunctionSpace<
                 ::Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>>
                 mesh,
             int nFieldVariablesOfMesh, const PythonConfig &specificSettings,
             double currentTime);
};

/** Partial specialization for unstructured mesh.
 *  Outputs an unstructured grid.
 */
template <int D, typename BasisFunctionType,
          typename FieldVariablesForOutputWriterType>
class HDF5Writer<
    FunctionSpace::FunctionSpace<::Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    FieldVariablesForOutputWriterType> {
public:
  //! write paraview file to given fileID, only output fieldVariables that are
  //! on a mesh with the given meshName
  static void outputFile(
      HDF5Utils::Group &group, FieldVariablesForOutputWriterType fieldVariables,
      const std::string &meshName,
      std::shared_ptr<FunctionSpace::FunctionSpace<
          ::Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>
          mesh,
      int nFieldVariablesOfMesh, const PythonConfig &specificSettings,
      double currentTime);
};
} // namespace OutputWriter

#include "output_writer/hdf5/hdf5_writer.tpp"
