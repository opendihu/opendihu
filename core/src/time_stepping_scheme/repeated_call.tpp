#include "time_stepping_scheme/repeated_call.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename Solver>
RepeatedCall<Solver>::RepeatedCall(DihuContext context) :
  TimeSteppingScheme(context["RepeatedCall"]), solver_(context_)
{
  this->specificSettings_ = this->context_.getPythonConfig();
}

template<typename Solver>
void RepeatedCall<Solver>::
initialize()
{
  if (initialized_)
    return;
 
  TimeSteppingScheme::initialize();
  LOG(TRACE) << "RepeatedCall::initialize";

  // initialize underlying Solver object, also with time step width
  solver_.initialize();
}

template<typename Solver>
void RepeatedCall<Solver>::advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "RepeatedCall::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "RepeatedCall, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // set sub time step in solver
    this->solver_.setTimeSpan(currentTime, currentTime+this->timeStepWidth_);

    // advance solver
    this->solver_.advanceTimeSpan();

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename Solver>
void RepeatedCall<Solver>::
run()
{
  initialize();
  advanceTimeSpan();
}

} // namespace