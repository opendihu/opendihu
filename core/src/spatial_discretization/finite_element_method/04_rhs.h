#pragma once

#include "spatial_discretization/finite_element_method/03_assemble_rhs.h"

namespace SpatialDiscretization {
/** base class implementing right hand side, that can be set by user for poisson
 * equation
 */
template <typename FunctionSpaceType, typename QuadratureType, int nComponents,
          typename Term>
class FiniteElementMethodRhs
    : public AssembleRightHandSide<FunctionSpaceType, QuadratureType,
                                   nComponents, Term> {
public:
  //! use constructor of base class
  using AssembleRightHandSide<FunctionSpaceType, QuadratureType, nComponents,
                              Term>::AssembleRightHandSide;

  friend class StiffnessMatrixTester; //< a class used for testing

  //! read in rhs values from config and transfer to weak form
  void setRightHandSide();
};

/** rhs featuring active stress, from activation parameter field gamma
 */
template <typename FunctionSpaceType, typename QuadratureType, int nComponents>
class FiniteElementMethodRhs<FunctionSpaceType, QuadratureType, nComponents,
                             Equation::Static::LinearElasticityActiveStress>
    : public AssembleRightHandSide<
          FunctionSpaceType, QuadratureType, nComponents,
          Equation::Static::LinearElasticityActiveStress> {
public:
  //! use constructor of base class
  using AssembleRightHandSide<
      FunctionSpaceType, QuadratureType, nComponents,
      Equation::Static::LinearElasticityActiveStress>::AssembleRightHandSide;

  friend class StiffnessMatrixTester; //< a class used for testing

  //! assemble active stress values into rhs
  void setRightHandSide();
};

} // namespace SpatialDiscretization

#include "spatial_discretization/finite_element_method/04_rhs.tpp"
#include "spatial_discretization/finite_element_method/04_rhs_active_stress.tpp"
