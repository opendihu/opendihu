#include "checkpointing/json/combined.h"

#include <scr.h>
#include "input_reader/json/partial.h"
#include "utility/path.h"

namespace Checkpointing {
namespace Json {
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
      break;
    }

    found = true;
    char checkpointing_dir[SCR_MAX_FILENAME];
    SCR_Start_restart(checkpointing_dir);

    char scr_file[SCR_MAX_FILENAME];
    if (checkpointToRestore != "") {
      SCR_Route_file(checkpointToRestore.c_str(), scr_file);
    } else {
      SCR_Route_file(ckpt_name, scr_file);
    }

    if (!Path::fileExists(scr_file)) {
      break;
    }

    int valid = 1;
    InputReader::Json::Partial r(scr_file);

    int32_t step;
    double newTime;
    if (!r.hasAttribute("timeStepNo") || !r.readAttr("timeStepNo", step)) {
      break;
    }
    if (!r.hasAttribute("currentTime") || !r.readAttr("currentTime", newTime)) {
      break;
    }
    if (!data.restoreState(r)) {
      break;
    }

    LOG(DEBUG) << "Successfully restored checkpointing timeStepNo: " << step
               << " | currentTime: " << newTime;
    timeStepNo = step;
    currentTime = newTime;

    int rc = SCR_Complete_restart(valid);
    restarted = (rc == SCR_SUCCESS);
  }

  if (!restarted && found) {
    LOG(ERROR) << "Failed to reload previous state";
    SCR_Complete_restart(0);
  }

  return restarted;
}
} // namespace Json
} // namespace Checkpointing
