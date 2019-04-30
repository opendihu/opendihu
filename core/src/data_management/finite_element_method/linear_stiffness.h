#pragma once

#include <Python.h>  // has to be the first included header

namespace Data
{

/** class that holds the linear stiffness parameters
 */
template<typename FunctionSpaceType,int nComponents>
class LinearStiffness :
  public FiniteElementsBase<FunctionSpaceType,nComponents>
{
public:
  //! constructor
  using FiniteElementsBase<FunctionSpaceType,nComponents>::FiniteElementsBase;

  //! initialize stifness parameters, then call the initialize method of the base class
  virtual void initialize();

  //! get the value of the 2nd order stiffness tensor, C_abcd
  double linearStiffness(int a, int b, int c, int d) const;

  //! set values for active stress
  void setActiveStress(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,9>> activeStress);

protected:
  double bulkModulus_;   ///< material parameter bulk modulus, symbol K
  double shearModulus_;  ///< material parameter shear modulus, symbol μ

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,9>> activeStress_;  ///< active stress field variable, this is a 3x3 tensor
};

}  // namespace

#include "data_management/finite_element_method/linear_stiffness.tpp"