#pragma once

#include <Python.h> // has to be the first included header
#include <memory>

#include "control/types.h"
#include "checkpointing/generic.h"

namespace Checkpointing {
// forward declaration
class Handle;

//! This class handles reading and parsing configuration values used in
//! checkpointing. It also provides an interface to create a
//! `Checkpointing::Handle`, a construct that is needed because the
//! Checkpointing::Manager can not own a DiHuContext, even indirect, because the
//! Manager is created as a static variable. Holding an DiHuContext will cause
//! the application to crash as this would result in the cleanup of the
//! DiHuContext after we already cleaned up the static variables.
class Manager {
public:
  //! Constructor for initializing a new Manager with their respective settings.
  //! If a key is missing in the settings objects, we will provide a sensible
  //! default.
  Manager(PythonConfig specificSettings);

  //! This will actually initialize a new checkpointing object, called `Handle`,
  //! that can be used to interact with the underlying checkpointing
  //! implementation. This `Handle` will own the underlying implementation and
  //! indirect the DiHuContext.
  std::shared_ptr<Handle>
  initialize(DihuContext context,
             std::shared_ptr<Partition::RankSubset> rankSubset = nullptr);

  //! Returns the timestepping interval.
  int32_t getInterval() const;
  //! Returns the prefix directory in which we store timesteps.
  const char *getPrefix() const;
  //! Returns the name of the checkpoint we want to restore. Will be an empty
  //! string if we do not want to restore a specific checkpoint.
  const std::string &getCheckpointToRestore() const;

private:
  PythonConfig specificSettings_; //< config for this object

  int32_t interval_; //< Timestep interval in which we create a new checkpoint.
  std::string prefix_; //< Prefix of directory where we store all checkpoints.
  std::string type_;   //< Type of the checkpoint implementation
  bool autoRestore_;   //< If set it restores a checkpoint if we find a loadable
                       // checkpoint
  std::string checkpointToRestore_; //< Specific checkpoint we'd like to restore
};

//! Handle to the actual checkpointing implementation that wraps the calls for
//! creating a new checkpoint as well as restoring one. This handle is used by
//! the different timestepping implementation. It implicitly owns an
//! DiHuContext.
class Handle {
public:
  //! Construct a new Handle by passing in a generic checkpoint implementation
  //! as well if we want to automatically restore a checkpoint and if we want to
  //! restore a specific checkpoint. Additionally, we also pass in the prefix of
  //! the directory where we store all checkpoints.
  Handle(std::shared_ptr<Generic> checkpointing, bool autoRestore,
         const std::string &checkpointToRestore, const std::string &prefix);

  //! Create a new checkpoint with the problem data as well as the timestep and
  //! the current time.
  //! This is a wrapper for the underlying checkpointing implementation and
  //! calls `checkpointing_`.
  template <typename DataType>
  void createCheckpoint(DataType &problemData, int timeStepNo = -1,
                        double currentTime = 0.0) const;

  //! Restore a checkpoint into the output variables of the problemData,
  //! timeStepNo and currentTime. This will return true if it successfully
  //! restored a checkpoint and will return false if it did not.
  //! If it failed to restore a checkpoint then it guarantees that the
  //! problemData, timeStepNo and currentTime were not changed.
  //! This is a wrapper for the underlying checkpointing implementation and
  //! calls `checkpointing_`.
  template <typename DataType>
  bool restore(DataType &problemData, int &timeStepNo,
               double &currentTime) const;

  //! Will return if we need to create a checkpoint for the current timestep.
  //! This is a wrapper for the underlying checkpointing implementation and
  //! calls `checkpointing_`.
  bool needCheckpoint();
  //! Returns true if the application needs to be terminated, in that case the
  //! timestepping loop should be exited.
  //! This is a wrapper for the underlying checkpointing implementation and
  //! calls `checkpointing_`.
  bool shouldExit();

private:
  std::shared_ptr<Generic> checkpointing_; //< Handle to the current underlying
                                           // checkpointing implementation.
  std::string prefix_; //< Prefix of directory where we store all checkpoints.
  bool autoRestore_;   //< If set it restores a checkpoint if we find a loadable
                       // checkpoint
  std::string checkpointToRestore_; //< Specific checkpoint we'd like to restore
};
} // namespace Checkpointing

#include "checkpointing/manager.tpp"
