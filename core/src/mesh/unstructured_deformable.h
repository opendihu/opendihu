#pragma once

#include <Python.h>
#include <vector>
#include <petscmat.h>

#include "control/types.h"
#include "mesh/mesh.h"
#include "mesh/deformable.h"
#include "function_space/00_function_space_base_dim.h"

namespace Mesh {

/**
 * An arbitrary mesh where each element can be adjacent to any other element or
 * at the boundary of the computational domain. There is no restriction that the
 * total domain must be quadratic or cubic.
 */
template <int D>
class UnstructuredDeformableOfDimension : public MeshOfDimension<D>,
                                          public Deformable {
public:
  //! constructor of base class
  using MeshOfDimension<D>::MeshOfDimension;

  //! get the total number of elements, this is implemented in
  //! function_space/04_function_space_data_unstructured.h
  // element_no_t nElementsLocal() const;

private:
};

} // namespace Mesh
