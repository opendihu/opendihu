#pragma once

#include <Python.h> // has to be the first included header

#include <array>
#include "mesh/face_t.h"
#include "mesh/type_traits.h"
#include "function_space/06_function_space_dofs_nodes.h"

namespace FunctionSpace {

/** Adds functionality to get dof numbers corresponding to the faces of the
 * element.
 */
template <typename MeshType, typename BasisFunctionType, typename = MeshType>
class FunctionSpaceFaces {};

/** Partial specialization for 1D meshes. In this case a face is a single point
 * (1 dof)
 */
template <typename MeshType, typename BasisFunctionType>
class FunctionSpaceFaces<MeshType, BasisFunctionType, Mesh::isDim<1, MeshType>>
    : public FunctionSpaceDofsNodes<MeshType, BasisFunctionType> {
public:
  //! inherit constructor
  using FunctionSpaceDofsNodes<MeshType,
                               BasisFunctionType>::FunctionSpaceDofsNodes;

  //! get all dof indices of a face, note: dimension in FunctionSpaceBaseDim is
  //! current-1 (=0), in this case the dofIndices array has exactly so many
  //! entries as there are dofs for a node
  static void getFaceDofs(
      Mesh::face_t face,
      std::array<dof_no_t,
                 FunctionSpaceBaseDim<1, BasisFunctionType>::nDofsPerNode()>
          &dofIndices);

  //! get the neighbouring elemental node index in given direction inside one
  //! element or -1 if there is no such node in the element in that direction
  static int getNeighbourNodeIndex(int nodeIndex, Mesh::face_t face);
};

/** Partial specialization for 2D meshes. In this case a face is a line
 */
template <typename MeshType, typename BasisFunctionType>
class FunctionSpaceFaces<MeshType, BasisFunctionType, Mesh::isDim<2, MeshType>>
    : public FunctionSpaceDofsNodes<MeshType, BasisFunctionType> {
public:
  //! inherit constructor
  using FunctionSpaceDofsNodes<MeshType,
                               BasisFunctionType>::FunctionSpaceDofsNodes;

  //! get all dof indices of a face, note: dimension in FunctionSpaceBaseDim is
  //! current-1 (=1)
  static void getFaceDofs(
      Mesh::face_t face,
      std::array<
          dof_no_t,
          FunctionSpaceBaseDim<1, BasisFunctionType>::nNodesPerElement() *
              FunctionSpaceBaseDim<MeshType::dim(),
                                   BasisFunctionType>::nDofsPerNode()>
          &dofIndices);

  //! get the neighbouring elemental node index in given direction inside one
  //! element or -1 if there is no such node in the element in that direction
  static int getNeighbourNodeIndex(int nodeIndex, Mesh::face_t face);
};

/** Partial specialization for 3D meshes. In this case a face is a real 2D face.
 */
template <typename MeshType, typename BasisFunctionType>
class FunctionSpaceFaces<MeshType, BasisFunctionType, Mesh::isDim<3, MeshType>>
    : public FunctionSpaceDofsNodes<MeshType, BasisFunctionType> {
public:
  //! inherit constructor
  using FunctionSpaceDofsNodes<MeshType,
                               BasisFunctionType>::FunctionSpaceDofsNodes;

  //! get all dof indices of a face, note: dimension in FunctionSpaceBaseDim is
  //! current-1 (=2)
  static void getFaceDofs(
      Mesh::face_t face,
      std::array<
          dof_no_t,
          FunctionSpaceBaseDim<2, BasisFunctionType>::nNodesPerElement() *
              FunctionSpaceBaseDim<MeshType::dim(),
                                   BasisFunctionType>::nDofsPerNode()>
          &dofIndices);

  //! get the neighbouring elemental node index in given direction inside one
  //! element or -1 if there is no such node in the element in that direction
  static int getNeighbourNodeIndex(int nodeIndex, Mesh::face_t face);
};

} // namespace FunctionSpace

#include "function_space/07_function_space_faces.tpp"
