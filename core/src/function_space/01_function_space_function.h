#pragma once

#include <Python.h> // has to be the first included header

#include <array>
#include "control/types.h"

#include "function_space/00_function_space_base_dim.h"
#include "mesh/mesh.h"

namespace FunctionSpace {
/** This class adds the basis function to the mesh.
 *  The basis function can be accessed by phi() and the derivative by
 * dphi_dxi().
 */
template <typename MeshType, typename BasisFunctionType>
class FunctionSpaceFunction
    : public FunctionSpaceBaseDim<MeshType::dim(), BasisFunctionType>,
      public MeshType {
public:
  using MeshType::MeshType;

  //! evaluate the basis function corresponding to element-local dof dofIndex at
  //! xi, xi lives in [0,1]^D
  static double phi(int dofIndex, std::array<double, MeshType::dim()> xi);

  //! evaluate the derivative of Phi(xi) w.r.t xi_i, where i is given by
  //! derivativeIdx, i.e. Phi_{dofIndex,derivativeIdx}(xi)
  static double dphi_dxi(int dofIndex, int derivativeIdx,
                         std::array<double, MeshType::dim()> xi);

  //! evaluate the first derivative of the basis function corresponding to
  //! element-local dof dofIndex at xi, interval for xi is [0,1]^D
  static std::array<double, MeshType::dim()>
  gradPhi(int dofIndex, std::array<double, MeshType::dim()> xi);

private:
  //! given an element-local dofIndex and dimension No (0 <= dimNo < D), return
  //! the basis function index in that direction
  static int getBasisFunctionIndex1D(int dofIndex, int dimNo);
};

} // namespace FunctionSpace

#include "function_space/01_function_space_function.tpp"
