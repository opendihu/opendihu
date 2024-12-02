#pragma once

#include <Python.h> // has to be the first included header

#include "control/dihu_context.h"
#include "checkpointing/generic.h"
#include "output_writer/hdf5/hdf5.h"

namespace Checkpointing {
namespace HDF5 {
//! The concreate implementation for HDF5 combined checkpointing
class Combined : public Generic {
public:
  //! Constructs a new combined HDF5 checkpointing implementation
  Combined(DihuContext context,
           std::shared_ptr<Partition::RankSubset> rankSubset = nullptr,
           const std::string &prefix = ".");

  //! Create a new checkpoint with the problem data as well as the timestep and
  //! the current time.
  template <typename DataType>
  void createCheckpoint(DataType &problemData, int timeStepNo = -1,
                        double currentTime = 0.0) const;
  //! Restore a checkpoint into the output variables of the problemData,
  //! timeStepNo and currentTime. This will return true if it successfully
  //! restored a checkpoint and will return false if it did not.
  //! If it failed to restore a checkpoint then it guarantees that the
  //! problemData, timeStepNo and currentTime were not changed.
  template <typename DataType>
  bool restore(DataType &problemData, int &timeStepNo, double &currentTime,
               bool autoRestore,
               const std::string &checkpointingToRestore) const;

private:
  std::unique_ptr<OutputWriter::HDF5>
      writer_; //< Constructed output writer that is used to create new
               // checkpoints
};
} // namespace HDF5
} // namespace Checkpointing

#include "checkpointing/hdf5/combined.tpp"
