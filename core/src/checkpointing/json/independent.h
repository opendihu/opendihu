#pragma once

#include <Python.h> // has to be the first included header

#include "control/dihu_context.h"
#include "checkpointing/generic.h"
#include "output_writer/json/json.h"

namespace Checkpointing {
namespace Json {
//! The concreate implementation for Json independent checkpointing
class Independent : public Generic {
public:
  //! Constructs a new independent Json checkpointing implementation
  Independent(DihuContext context,
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
  std::unique_ptr<OutputWriter::Json>
      writer_; //< Constructed output writer that is used to create new
               // checkpoints
};
} // namespace Json
} // namespace Checkpointing

#include "checkpointing/json/independent.tpp"
