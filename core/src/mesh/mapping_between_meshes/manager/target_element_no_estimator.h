#pragma once

#include <Python.h> // has to be the first included header

namespace MappingBetweenMeshes {

/** Helper class that gives a guess for the target function space element no
 * that contains a given source dof no. It helps if the actual element is a
 * neighbor to the predicted one (or is the predicted element itself), because
 * then the algorithm does not have to iterate over all elements to find the
 * correct element. For structured 3D elements, this a heuristic is implemented,
 * for other meshes, nothing is done.
 *
 *  This class is the version for general meshes and does nothing.
 */
template <typename SourceFunctionSpaceType, typename TargetFunctionSpaceType>
class TargetElementNoEstimator {
public:
  TargetElementNoEstimator(
      std::shared_ptr<SourceFunctionSpaceType> sourceFunctionSpace,
      std::shared_ptr<TargetFunctionSpaceType> targetFunctionSpace){};

  void estimateElementNo(dof_no_t sourceDofNoLocal,
                         element_no_t &targetElementNo){};

protected:
};

/** This is the partial specialization for two 3D structured meshes. It computes
 * the target element no.
 */
template <typename BasisFunctionType1, typename BasisFunctionType2>
class TargetElementNoEstimator<
    FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,
                                 BasisFunctionType1>,
    FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,
                                 BasisFunctionType2>> {
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,
                                       BasisFunctionType1>
      SourceFunctionSpaceType;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,
                                       BasisFunctionType2>
      TargetFunctionSpaceType;

  //! constructor
  TargetElementNoEstimator(
      std::shared_ptr<SourceFunctionSpaceType> sourceFunctionSpace,
      std::shared_ptr<TargetFunctionSpaceType> targetFunctionSpace);

  //! try to improve the estimation for targetElementNo, in which
  //! sourceDofNoLocal should be
  void estimateElementNo(dof_no_t sourceDofNoLocal,
                         element_no_t &targetElementNo);

protected:
  std::shared_ptr<SourceFunctionSpaceType>
      sourceFunctionSpace_; //< the source function space of the mapping
  std::shared_ptr<TargetFunctionSpaceType>
      targetFunctionSpace_; //< the target function space of the mapping

  element_no_t targetElementForYChange_; //< target element that will be used
                                         //when the source y coordinate changes
  element_no_t targetElementForZChange_; //< target element that will be used
                                         //when the source z coordinate changes
};

/** This is the partial specialization for two 3D structured meshes. It computes
 * the target element no.
 */
template <typename BasisFunctionType1, typename BasisFunctionType2>
class TargetElementNoEstimator<
    FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,
                                 BasisFunctionType1>,
    FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>,
                                 BasisFunctionType2>> {
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,
                                       BasisFunctionType1>
      SourceFunctionSpaceType;
  typedef FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>,
                                       BasisFunctionType2>
      TargetFunctionSpaceType;

  //! constructor
  TargetElementNoEstimator(
      std::shared_ptr<SourceFunctionSpaceType> sourceFunctionSpace,
      std::shared_ptr<TargetFunctionSpaceType> targetFunctionSpace);

  //! try to improve the estimation for targetElementNo, in which
  //! sourceDofNoLocal should be
  void estimateElementNo(dof_no_t sourceDofNoLocal,
                         element_no_t &targetElementNo);

protected:
  std::shared_ptr<SourceFunctionSpaceType>
      sourceFunctionSpace_; //< the source function space of the mapping
  std::shared_ptr<TargetFunctionSpaceType>
      targetFunctionSpace_; //< the target function space of the mapping

  element_no_t targetElementForYChange_; //< target element that will be used
                                         //when the source y coordinate changes
  element_no_t targetElementForZChange_; //< target element that will be used
                                         //when the source z coordinate changes
};

} // namespace MappingBetweenMeshes

#include "mesh/mapping_between_meshes/manager/target_element_no_estimator.tpp"
