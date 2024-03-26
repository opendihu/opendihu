#pragma once

#include <Python.h> // has to be the first included header

#include "interfaces/runnable.h"
#include "control/precice/checkpoint.h"

namespace Control {

/** Generic Precice adapter for surface coupling, can be configured to either
 * prescribe Neumann or Dirichlet boundary conditions.
 */
template <typename NestedSolver>
class PreciceAdapter : public Runnable,
                       public PreciceAdapterCheckpoint<NestedSolver> {
public:
  //! define the type of the data object
  typedef typename NestedSolver::Data Data;

  //! constructor
  using PreciceAdapterCheckpoint<NestedSolver>::PreciceAdapterCheckpoint;

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! return the data object, with the call to this method the output writers
  //! get the data to create their output files
  Data &data();
};

} // namespace Control

#include "control/precice/precice_adapter.tpp"