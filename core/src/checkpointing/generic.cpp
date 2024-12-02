#include "checkpointing/generic.h"

#include <scr.h>

namespace Checkpointing {
Generic::Generic(const std::string &prefix) : prefix_(prefix) {
  if (prefix_ == "") {
    prefix_ = ".";
  }
}

bool Generic::needCheckpoint() {
  int32_t v;
  SCR_Need_checkpoint(&v);
  return v;
}

bool Generic::shouldExit() {
  int32_t v;
  SCR_Should_exit(&v);
  return v;
}
} // namespace Checkpointing
