#pragma once

#include <Python.h> // has to be the first included header

#include "control/dihu_context.h"

namespace Checkpointing {

//! Base checkpointing class that all underlying checkpointing implementation
//! need to implement. It provides a default implementation for `needCheckpoint`
//! and `shouldExit` that use `scr` in the background, but this can be changed
//! if needed.
class Generic {
public:
  //! Constructor that takes in a prefix directory. Defaults to the current
  //! directory. Constructor also ensure that we don't end up with a empty
  //! prefix string.
  Generic(const std::string &prefix = ".");

  //! Default desctructor, does nothing
  virtual ~Generic() = default;

  //! Default implementation that wraps `SCR_Need_checkpoint`
  bool needCheckpoint();
  //! Default implementation that wraps `SCR_Should_exit`
  bool shouldExit();

protected:
  const int32_t maxAttempt =
      3; //< We can have multiple attempts to load a single checkpoint. The max
         // attempts are currently hardcoded in the base implementation.
  std::string prefix_; //< Prefix of directory where we store all checkpoints.
};
} // namespace Checkpointing
