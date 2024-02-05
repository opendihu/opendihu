#pragma once

#include <Python.h> // has to be the first included header

#include "control/precice/surface_coupling/02_read_write.h"

namespace Control {

/** Generic Precice adapter, can be configured to either prescribe Neumann or
 * Dirichlet boundary conditions.
 */
template <typename NestedSolver>
class PreciceAdapterCheckpoint : public PreciceAdapterReadWrite<NestedSolver> {
public:
  //! constructor
  using PreciceAdapterReadWrite<NestedSolver>::PreciceAdapterReadWrite;

  //! save the current state of the simulation to a checkpoint
  void saveCheckpoint(double currentTime);

  //! load all values from the last saved checkpoint
  //! @return current simulation time
  double loadCheckpoint();

protected:
  double savedCurrentTime_;     //< the saved value of the simulation time
  Vec savedState_ = PETSC_NULL; //< Petsc Vec that stores the whole data
};

} // namespace Control

#include "control/precice/surface_coupling/03_checkpoint.tpp"
