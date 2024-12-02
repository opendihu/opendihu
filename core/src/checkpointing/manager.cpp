#include "checkpointing/manager.h"

#include "checkpointing/hdf5/combined.h"
#include "checkpointing/hdf5/independent.h"
#include "checkpointing/json/combined.h"
#include "checkpointing/json/independent.h"

namespace Checkpointing {
Manager::Manager(PythonConfig specificSettings)
    : specificSettings_(specificSettings),
      interval_(specificSettings.getOptionInt("interval", 1)),
      prefix_(specificSettings.getOptionString("directory", "state")),
      type_(specificSettings.getOptionString("type", "hdf5-combined")),
      autoRestore_(specificSettings.getOptionBool("autoRestore", true)),
      checkpointToRestore_(
          specificSettings.getOptionString("checkpointToRestore", "")) {}

std::shared_ptr<Handle>
Manager::initialize(DihuContext context,
                    std::shared_ptr<Partition::RankSubset> rankSubset) {
  std::shared_ptr<Generic> checkpointing = nullptr;
  if (type_ == "hdf5-combined") {
    checkpointing =
        std::make_shared<HDF5::Combined>(context, rankSubset, this->prefix_);
  } else if (type_ == "hdf5-independent") {
    checkpointing =
        std::make_shared<HDF5::Independent>(context, rankSubset, this->prefix_);
  } else if (type_ == "json-combined") {
    checkpointing =
        std::make_shared<Json::Combined>(context, rankSubset, this->prefix_);
  } else if (type_ == "json-independent") {
    checkpointing =
        std::make_shared<Json::Independent>(context, rankSubset, this->prefix_);
  } else {
    LOG(ERROR) << "checkpointing type: " << type_
               << " is not a valid type. Make sure to either configure it "
                  "with the type 'hdf5-combined', 'hdf5-independent', "
                  "'json-combined' or 'json-independent'";
  }

  if (checkpointing) {
    return std::make_shared<Handle>(checkpointing, autoRestore_,
                                    checkpointToRestore_, this->prefix_);
  } else {
    return nullptr;
  }
}

int32_t Manager::getInterval() const { return interval_; }
const char *Manager::getPrefix() const { return prefix_.c_str(); }
const std::string &Manager::getCheckpointToRestore() const {
  return checkpointToRestore_;
}

Handle::Handle(std::shared_ptr<Generic> checkpointing, bool autoRestore,
               const std::string &checkpointToRestore,
               const std::string &prefix)
    : checkpointing_(checkpointing), prefix_(prefix), autoRestore_(autoRestore),
      checkpointToRestore_(checkpointToRestore) {}

bool Handle::needCheckpoint() {
  if (checkpointing_) {
    return checkpointing_->needCheckpoint();
  }
  return false;
}

bool Handle::shouldExit() {
  if (checkpointing_) {
    return checkpointing_->shouldExit();
  }
  return false;
}

} // namespace Checkpointing
