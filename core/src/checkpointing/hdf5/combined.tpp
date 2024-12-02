#include "checkpointing/hdf5/combined.h"

#include <scr.h>
#include "input_reader/hdf5/partial.h"
#include "utility/path.h"

namespace Checkpointing {
namespace HDF5 {
template <typename DataType>
void Combined::createCheckpoint(DataType &data, int timeStepNo,
                                double currentTime) const {
  Control::PerformanceMeasurement::start("durationWriteCheckpoint");

  char ckpt_name[SCR_MAX_FILENAME];
  snprintf(ckpt_name, sizeof(ckpt_name), "timestep.%04d", timeStepNo);
  SCR_Start_output(ckpt_name, SCR_FLAG_CHECKPOINT);

  char checkpoint_file[256];
  sprintf(checkpoint_file, "%s/timestep.%04d", this->prefix_.c_str(),
          timeStepNo);

  char scr_file[SCR_MAX_FILENAME];
  SCR_Route_file(checkpoint_file, scr_file);

  writer_->write(data, scr_file, timeStepNo, currentTime);

  SCR_Complete_output(1);

  Control::PerformanceMeasurement::stop("durationWriteCheckpoint");
}

template <typename DataType>
bool Combined::restore(DataType &data, int &timeStepNo, double &currentTime,
                       bool autoRestore,
                       const std::string &checkpointToRestore) const {
  bool restarted = false;
  bool found = false;
  int32_t attempt = 1;
  while (!restarted) {
    if (!autoRestore) {
      break;
    }

    char ckpt_name[SCR_MAX_FILENAME];
    if (checkpointToRestore != "") {
      int32_t r = SCR_Current(checkpointToRestore.c_str());
      if (r != SCR_SUCCESS) {
        LOG(ERROR) << "Failed to specify specific checkpoint. Please check the "
                      "name for typos. Return value"
                   << r;
        return false;
      }
    }

    int32_t have_restart = 0;
    SCR_Have_restart(&have_restart, ckpt_name);
    if (!have_restart) {
      LOG(WARNING) << "We did not have a checkpoint to restart";
      break;
    }

    found = true;
    char checkpointing_dir[SCR_MAX_FILENAME];
    SCR_Start_restart(checkpointing_dir);

    char scr_file[SCR_MAX_FILENAME];
    int ret = 0;
    if (checkpointToRestore != "") {
      ret = SCR_Route_file(checkpointToRestore.c_str(), scr_file);
    } else {
      ret = SCR_Route_file(ckpt_name, scr_file);
    }

    if (ret == 1) {
      LOG(ERROR) << "File to restore not found.";
      break;
    }

    int valid = 1;
    InputReader::HDF5::Partial r(scr_file);

    std::string rawVersion;
    if (!r.hasAttribute("rawVersion") ||
        !r.readAttr("rawVersion", rawVersion)) {
      LOG(ERROR) << "rawVersion was not found in file, or we were not able to "
                    "load this attribute";
      break;
    }

    if (strcmp(rawVersion.c_str(), DihuContext::version().c_str()) != 0) {
      LOG(ERROR) << "Checkpoint was created with a different OpenDiHu version "
                    "and can not be restored! Current: "
                 << DihuContext::version() << " checkpoint: " << rawVersion;
      break;
    }
    int32_t worldSize;
    if (!r.hasAttribute("worldSize") || !r.readAttr("worldSize", worldSize)) {
      LOG(ERROR) << "worldSize was not found in file, or we were not able to "
                    "load this attribute";
      break;
    }
    if (worldSize != DihuContext::nRanksCommWorld()) {
      LOG(ERROR) << "Checkpoint was created with a different worldSize and can "
                    "not be restored! Current: "
                 << DihuContext::nRanksCommWorld()
                 << " checkpoint: " << worldSize;
      break;
    }

    int32_t step;
    double newTime;
    if (!r.hasAttribute("timeStepNo") || !r.readAttr("timeStepNo", step)) {
      LOG(ERROR) << "timeStepNo was not found in file, or we were not able to "
                    "load this attribute";
      break;
    }
    if (!r.hasAttribute("currentTime") || !r.readAttr("currentTime", newTime)) {
      LOG(ERROR) << "currentTime was not found in file, or we were not able to "
                    "load this attribute";
      break;
    }
    if (!data.restoreState(r)) {
      LOG(ERROR) << "Failed to restore state data";
      break;
    }

    LOG(INFO) << "Successfully restored checkpointing timeStepNo: " << step
              << " | currentTime: " << newTime;
    timeStepNo = step;
    currentTime = newTime;

    int rc = SCR_Complete_restart(valid);
    restarted = (rc == SCR_SUCCESS);

    if (!restarted) {
      if (attempt >= this->maxAttempt) {
        LOG(FATAL) << "We were not able to Successfully restore the checkpoint "
                      "on all ranks, only on some. Because of this we might a "
                      "have broken state on some ranks and can't continue "
                      "execution. Aborting application!";
      }

      LOG(ERROR) << "Fatal error, not all ranks were able to restore the "
                    "checkpoint. Trying again. This was attempt: "
                 << attempt << "/" << this->maxAttempt;
      attempt++;
    }
  }

  if (!restarted && found) {
    LOG(ERROR) << "Failed to reload previous state";
    SCR_Complete_restart(0);
  }

  return restarted;
}
} // namespace HDF5
} // namespace Checkpointing
