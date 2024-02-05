#include "time_stepping_scheme/repeated_call_static.h"

#include <Python.h> // has to be the first included header
#include <vector>

#include "utility/python_utility.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

namespace TimeSteppingScheme {

template <typename Solver>
RepeatedCallStatic<Solver>::RepeatedCallStatic(DihuContext context)
    : TimeSteppingScheme(context["RepeatedCallStatic"]), solver_(context_) {
  this->specificSettings_ = this->context_.getPythonConfig();
}

template <typename Solver> void RepeatedCallStatic<Solver>::initialize() {
  if (initialized_)
    return;

  TimeSteppingScheme::initialize();
  LOG(TRACE) << "RepeatedCallStatic::initialize";

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver(
      "RepeatedCall", true); // hasInternalConnectionToFirstNestedSolver=true
                             // (the last argument) means slot connector data is
                             // shared with the first subsolver

  // indicate in solverStructureVisualizer that now a child solver will be
  // initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // initialize underlying Solver object, also with time step width
  solver_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization
  // is done
  DihuContext::solverStructureVisualizer()->endChild();
}

template <typename Solver>
void RepeatedCallStatic<Solver>::advanceTimeSpan(
    bool withOutputWritersEnabled) {
  // avoid that solver structure file is created, this should only be done after
  // the whole simulation has finished
  DihuContext::solverStructureVisualizer()->disable();

  // run the solver which does not have advanceTimeSpan
  solver_.run();

  // enable again
  DihuContext::solverStructureVisualizer()->enable();
}

template <typename Solver> void RepeatedCallStatic<Solver>::run() {
  initialize();
  advanceTimeSpan();
}

template <typename Solver>
void RepeatedCallStatic<Solver>::callOutputWriter(int timeStepNo,
                                                  double currentTime,
                                                  int callCountIncrement) {
  this->solver_.callOutputWriter(timeStepNo, currentTime, callCountIncrement);
}

template <typename Solver>
typename RepeatedCallStatic<Solver>::Data &RepeatedCallStatic<Solver>::data() {
  return solver_.data();
}

template <typename Solver>
std::shared_ptr<typename RepeatedCallStatic<Solver>::SlotConnectorDataType>
RepeatedCallStatic<Solver>::getSlotConnectorData() {
  return solver_.getSlotConnectorData();
}

} // namespace TimeSteppingScheme
