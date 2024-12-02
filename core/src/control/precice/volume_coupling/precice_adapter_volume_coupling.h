#pragma once

#include <Python.h> // has to be the first included header

#include "interfaces/runnable.h"
#include "control/precice/volume_coupling/01_read_write.h"

namespace Control {

/** Generic Precice adapter for volume coupling.
 */
template <typename NestedSolver>
class PreciceAdapterVolumeCoupling
    : public PreciceAdapterVolumeCouplingReadWrite<NestedSolver> {
public:
  //! define the type of the data object
  typedef typename NestedSolver::Data Data;

  //! constructor
  using PreciceAdapterVolumeCouplingReadWrite<
      NestedSolver>::PreciceAdapterVolumeCouplingReadWrite;

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary
  //! ("uninitialize")
  void reset();

  //! set unique data prefix
  void setUniqueDataPrefix(const std::string &prefix);

  //! return the data object, with the call to this method the output writers
  //! get the data to create their output files
  Data &data();

  //! return reference to the full data object that stores everything for a
  //! checkpoint
  Data &fullData();

private:
  std::string uniqueDataPrefix_;
};

} // namespace Control

#include "control/precice/volume_coupling/precice_adapter_volume_coupling.tpp"
