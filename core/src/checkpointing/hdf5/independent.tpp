#include "checkpointing/hdf5/independent.h"

#include <scr.h>
#include "input_reader/hdf5/full_dataset.h"
#include "utility/path.h"

namespace Checkpointing {
namespace HDF5 {
template <typename DataType>
void Independent::createCheckpoint(DataType &data, int timeStepNo,
                                   double currentTime) const {
  Control::PerformanceMeasurement::start("durationWriteCheckpoint");

  int32_t ownRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &ownRank);

  std::stringstream ckptName;
  ckptName << "timestep." << std::setw(4) << std::setfill('0') << timeStepNo
           << "d";
  SCR_Start_output(ckptName.str().c_str(), SCR_FLAG_CHECKPOINT);

  std::stringstream checkpointDirStream;
  checkpointDirStream << prefix_ << "/timestep." << std::setw(4)
                      << std::setfill('0') << timeStepNo << "d";
  std::string checkpointDir = checkpointDirStream.str();

  if (ownRank == 0) {
    Path::mkpath(checkpointDir.c_str());
  }

  // hold all processes until directory is created
  MPI_Barrier(MPI_COMM_WORLD);

  std::stringstream checkpointFile;
  checkpointFile << checkpointDir << "/rank_" << ownRank << ".ckpt";

  char scr_file[SCR_MAX_FILENAME];
  SCR_Route_file(checkpointFile.str().c_str(), scr_file);

  writer_->write(data, scr_file, timeStepNo, currentTime);

  SCR_Complete_output(1);

  Control::PerformanceMeasurement::stop("durationWriteCheckpoint");
}

template <typename DataType>
bool Independent::restore(DataType &data, int &timeStepNo, double &currentTime,
                          bool autoRestore,
                          const std::string &checkpointToRestore) const {
  int32_t ownRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &ownRank);

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
    char checkpointDir[SCR_MAX_FILENAME];
    SCR_Start_restart(checkpointDir);

    std::stringstream checkpointFile;
    if (checkpointToRestore != "") {
      checkpointFile << checkpointToRestore << "/rank_" << ownRank << ".ckpt";
    } else {
      checkpointFile << checkpointDir << "/rank_" << ownRank << ".ckpt";
    }

    char scr_file[SCR_MAX_FILENAME];
    SCR_Route_file(checkpointFile.str().c_str(), scr_file);

    bool exists = Path::fileExists(scr_file);
    bool all_exists = false;
    MPIUtility::handleReturnValue(MPI_Allreduce(&exists, &all_exists, 1,
                                                MPI_C_BOOL, MPI_LAND,
                                                MPI_COMM_WORLD),
                                  "MPI_Allreduce");
    if (all_exists) {
      break;
    }

    int valid = 1;
    InputReader::HDF5::FullDataset r(scr_file);

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
} // namespace HDF5
} // namespace Checkpointing
