#pragma once

#include "equation/type_traits.h"
#include "control/types.h"

namespace SpatialDiscretization
{

/** partial specialization for finite elasticity, dimension 1
 */
template<typename EvaluationsType,typename FunctionSpaceType>
class IntegrandStiffnessMatrix<1,EvaluationsType,FunctionSpaceType,1,Equation::Static::LinearElasticity>
{
public:
  static EvaluationsType evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,1,Equation::Static::LinearElasticity> &data,
                                           const std::array<Vec3,1> &jacobian, element_no_t elementNoLocal,
                                           const std::array<double,1> xi);
};

/** partial specialization for finite elasticity, dimension 2
 */
template<typename EvaluationsType,typename FunctionSpaceType>
class IntegrandStiffnessMatrix<2,EvaluationsType,FunctionSpaceType,2,Equation::Static::LinearElasticity>
{
public:
  static EvaluationsType evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,2,Equation::Static::LinearElasticity> &data,
                                           const std::array<Vec3,2> &jacobian, element_no_t elementNoLocal,
                                           const std::array<double,2> xi);
};


/** partial specialization for finite elasticity, dimension 3
 */
template<typename EvaluationsType,typename FunctionSpaceType>
class IntegrandStiffnessMatrix<3,EvaluationsType,FunctionSpaceType,3,Equation::Static::LinearElasticity>
{
public:
  static EvaluationsType evaluateIntegrand(const Data::FiniteElements<FunctionSpaceType,3,Equation::Static::LinearElasticity> &data,
                                           const std::array<Vec3,3> &jacobian, element_no_t elementNoLocal,
                                           const std::array<double,3> xi);
};

} // namespace

#include "spatial_discretization/finite_element_method/integrand/integrand_stiffness_matrix_linear_elasticity.tpp"