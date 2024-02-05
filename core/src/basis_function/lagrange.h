#pragma once

#include "basis_function/basis_function.h"
#include "control/types.h"

namespace BasisFunction {

/**  Defines Lagrange ansatz functions on the interval xi in [0,1] of given
 * order. Space Q_k in literature. Transformation to the actual physical space
 * is done using a jacobian, see finite element method.
 */
template <int order = 1> class LagrangeOfOrder : public BasisFunction {
public:
  typedef LagrangeOfOrder<order>
      BasisFunctionUsingOnlyNodalValues; //< the same basis function, but only
                                         //using nodal values, for Lagrange this
                                         //is the same.  This construct is used
                                         //for Neumann BC, where the integration
                                         //over the surface/line is done using
                                         //Lagrange ansatz functions.

  //! number of degrees of freedom of this basis
  static constexpr int nDofsPerBasis();

  //! number of degrees of freedom associated with a node in world space
  static constexpr int nDofsPerNode();

  //! evaluate the 1D basis function corresponding to element-local dof i at xi,
  //! interval for xi is [0,1]
  static double phi(int i, double xi);

  //! evaluate the first derivative of the 1D basis function corresponding to
  //! element-local dof i at xi, interval for xi is [0,1]
  static double dphi_dxi(int i, double xi);

  //! return the basis order value as used in python files and callbacks, e.g. 2
  static constexpr int getBasisOrder();

  //! return a basis function type string as used in python files and callbacks,
  //! e.g. "Lagrange"
  static std::string getBasisFunctionString();

  static constexpr bool isNodalBased =
      true; //< specify that this basis function is nodal based
private:
};

} // namespace BasisFunction

#include "basis_function/lagrange.tpp"
