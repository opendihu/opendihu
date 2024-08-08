#include "control/precice/precice_adapter.h"

#include <sstream>

namespace Control {

template <typename NestedSolver> void PreciceAdapter<NestedSolver>::run() {
#ifdef HAVE_PRECICE

  // initialize everything
  this->initialize();

  double currentTime = 0;

  // if precice coupling is disabled in settings, run the timestep of the nested
  // solver until endTimeIfCouplingDisabled_ is reached
  if (!this->couplingEnabled_) {
    const int nTimeSteps =
        this->endTimeIfCouplingDisabled_ / this->timeStepWidth_;
    for (int timeStepNo = 0; timeStepNo < nTimeSteps; timeStepNo++) {
      if (timeStepNo % this->timeStepOutputInterval_ == 0 &&
          (this->timeStepOutputInterval_ <= 10 ||
           timeStepNo > 0)) // show first timestep only if
                            // timeStepOutputInterval is <= 10
      {
        LOG(INFO)
            << "PreCICE surface coupling (disabled for debugging), timestep "
            << timeStepNo << "/" << nTimeSteps << ", t=" << currentTime;
      }

      // set time span in nested solver
      this->nestedSolver_.setTimeSpan(currentTime,
                                      currentTime + this->timeStepWidth_);

      // call the nested solver to proceed with the simulation for the assigned
      // time span
      this->nestedSolver_.advanceTimeSpan();

      // increase current simulation time
      currentTime += this->timeStepWidth_;
    }
    return;
  }

  // initialize
  if (this->preciceParticipant_->requiresInitialData()) {
    this->preciceWriteData(this->nestedSolver_, this->preciceParticipant_,
                           this->preciceSurfaceData_, this->preciceVolumeData_,
                           this->scalingFactor_);
  }
  this->preciceParticipant_->initialize();
  this->maximumPreciceTimestepSize_ =
      this->preciceParticipant_->getMaxTimeStepSize();

  LOG(DEBUG) << "precice initialization done, dt: "
             << this->maximumPreciceTimestepSize_ << ","
             << this->timeStepWidth_;

  this->initialized_ = true;

  // assert that precice is properly initialized and the interface is available
  assert(this->preciceParticipant_);

  // main simulation loop of adapter
  for (int timeStepNo = 0; this->preciceParticipant_->isCouplingOngoing();
       timeStepNo++) {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 &&
        (this->timeStepOutputInterval_ <= 10 ||
         timeStepNo >
             0)) // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "preCICE coupling, timestep " << timeStepNo
                << ", t=" << currentTime;
    }

    // determine if checkpoint needs to be written
    if (this->preciceParticipant_->requiresWritingCheckpoint()) {
      // save checkpoint
      this->saveCheckpoint(currentTime);
      this->saveFiberData(this->nestedSolver_);
    }

    // read incoming values
    this->preciceReadData(this->nestedSolver_, this->preciceParticipant_,
                          this->preciceSurfaceData_, this->preciceVolumeData_,
                          this->context_);

    // compute the time step width such that it fits in the remaining time in
    // the current time window
    this->maximumPreciceTimestepSize_ =
        this->preciceParticipant_->getMaxTimeStepSize();
    double timeStepWidth =
        std::min(this->maximumPreciceTimestepSize_, this->timeStepWidth_);

    // set time span in nested solver
    this->nestedSolver_.setTimeSpan(currentTime, currentTime + timeStepWidth);

    // call the nested solver to proceed with the simulation for the assigned
    // time span the parameter specifies whether the output writers are enabled
    this->nestedSolver_.advanceTimeSpan(!this->outputOnlyConvergedTimeSteps_);

    // write outgoing data to precice
    this->preciceWriteData(this->nestedSolver_, this->preciceParticipant_,
                           this->preciceSurfaceData_, this->preciceVolumeData_,
                           this->scalingFactor_);

    // increase current simulation time
    currentTime += timeStepWidth;

    // advance timestepping in precice
    this->preciceParticipant_->advance(timeStepWidth);

    LOG(DEBUG) << "precice::advance(" << timeStepWidth
               << "), maximumPreciceTimestepSize_: "
               << this->maximumPreciceTimestepSize_;

    // if coupling did not converge, reset to previously stored checkpoint
    if (this->preciceParticipant_->requiresReadingCheckpoint()) {
      // set variables back to last checkpoint
      currentTime = this->loadCheckpoint();
      this->loadFiberData(this->nestedSolver_);
    }

    // if the current time step did converge and subcycling is complete
    if (this->preciceParticipant_->isTimeWindowComplete()) {
      if (this->outputOnlyConvergedTimeSteps_) {
        // output all data in the nested solvers
        this->nestedSolver_.callOutputWriter(timeStepNo, currentTime);
      }
    }

  } // loop over time steps

  // finalize precice interface
  this->preciceParticipant_->finalize();

#else
  LOG(FATAL) << "Not compiled with preCICE!";
#endif
}

template <typename NestedSolver>
typename PreciceAdapter<NestedSolver>::Data &
PreciceAdapter<NestedSolver>::data() {
  // get a reference to the data object
  return this->nestedSolver_.data();
}

} // namespace Control
