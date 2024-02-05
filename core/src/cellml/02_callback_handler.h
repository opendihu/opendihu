#pragma once

#include <vector>

#include "interfaces/runnable.h"
#include "control/dihu_context.h"
#include "output_writer/manager.h"
#include "function_space/function_space.h"
#include "basis_function/lagrange.h"
#include "cellml/01_rhs_routine_handler.h"

/** The is a class that contains cellml equations and can be used with a time
 * stepping scheme. The nStates template parameter specifies the number of state
 * variables that should be used with the integrator. It is necessary that this
 * value is fixed at compile time because the timestepping scheme needs to know
 * which field variable types is has to construct. This class can also be
 * computed easily in multiple instances along the nodes of a mesh.
 *
 *  Naming:
 *   Algebraic (opendihu) = KNOWN (OpenCMISS) = Algebraic (OpenCOR)
 *   Parameter (opendihu, OpenCMISS) = KNOWN (OpenCMISS), in OpenCOR also
 * algebraic Constant - these are constants that are only present in the source
 * files State: state variable Rate: the time derivative of the state variable,
 * i.e. the increment value in an explicit Euler stepping
 */
template <int nStates, int nAlgebraics_, typename FunctionSpaceType>
class CallbackHandler
    : public RhsRoutineHandler<nStates, nAlgebraics_, FunctionSpaceType> {
public:
  //! constructor
  CallbackHandler(DihuContext context);

  //! constructor
  CallbackHandler(DihuContext context,
                  const typename CellmlAdapterBase<
                      nStates, nAlgebraics_, FunctionSpaceType>::Data &rhsData);

  //! destructor
  virtual ~CallbackHandler();

  //! register a callback function setParameters that can set parameter values
  //! before each computation
  void registerSetParameters(void (*setParameters)(
      void *context, int nInstances, int timeStepNo, double currentTime,
      std::vector<double> &parameters));

  //! register a callback function setSpecificParameters that can set parameter
  //! values before each computation
  void registerSetSpecificParameters(void (*setSpecificParameters)(
      void *context, int nInstances, int timeStepNo, double currentTime,
      std::vector<double> &localParameters));

  //! register a callback function setSpecificStates that can set state values
  //! before each computation
  void registerSetSpecificStates(void (*setSpecificStatesParameters)(
      void *context, int nInstances, int timeStepNo, double currentTime,
      double *states));

  //! register a callbackfunction handleResult that gets called after each new
  //! values are available
  void registerHandleResult(void (*handleResult)(
      void *context, int nInstances, int timeStepNo, double currentTime,
      double *states, double algebraics[]));

  //! directly call the python callback if it exists
  void callPythonSetParametersFunction(int nInstances, int timeStepNo,
                                       double currentTime,
                                       double *parameterValues,
                                       int nParameters);

  //! directly call the python callback if it exists
  void callPythonSetSpecificParametersFunction(int nInstances, int timeStepNo,
                                               double currentTime,
                                               double *localParameterValues,
                                               int nLocalParameters);

  //! directly call the python callback if it exists
  void callPythonSetSpecificStatesFunction(int nInstances, int timeStepNo,
                                           double currentTime, double *states);

  //! directly call the python callback if it exists
  void callPythonHandleResultFunction(int nInstances, int timeStepNo,
                                      double currentTime, double *states,
                                      double *algebraics);

  //! get the values of this->lastCallSpecificStatesTime
  double lastCallSpecificStatesTime();

  //! set the value of this->lastCallSpecificStatesTime
  void setLastCallSpecificStatesTime(double lastCallSpecificStatesTime);

protected:
  //! construct the python call back functions from config
  virtual void initializeCallbackFunctions();

  //! call Py_CLEAR on all python objects
  void clearPyObjects();

  int setSpecificParametersCallInterval_; //< setSpecificParameters_ will be
                                          // called every callInterval_ time
                                          // steps
  int setSpecificStatesCallInterval_;     //< setSpecificStates_ will be called
                                          // every callInterval_ time steps
  int handleResultCallInterval_;          //< handleResult will be called every
                                          // callInterval_ time steps

  double
      setSpecificStatesCallFrequency_; //< frequency, after which the
                                       // setSpecificStates callback function
                                       // will be called, either this condition
                                       // or the condition with
                                       // setSpecificStatesCallInterval_ is used
  std::vector<double>
      setSpecificStatesFrequencyJitter_; //< relative jitter values: factors of
                                         // setSpecificStatesCallFrequency_,
                                         // random jitter to add or substract
                                         // from frequency
  double currentJitter_; //< the absolute value of the current jitter
  int jitterIndex_;      //< which of the stored jitter values in
                         // setSpecificStatesFrequencyJitter_ to use
  int fiberNoGlobal_;    //< the additionalArgument converted to an integer,
                         // interpreted as the global fiber no and used in the
                         // stimulation log

  double lastCallSpecificStatesTime_; //< last time the setSpecificStates_
                                      // method was called
  double
      setSpecificStatesRepeatAfterFirstCall_; //< duration of continuation of
                                              // calling the setSpecificStates
                                              // callback after it was triggered
  double
      setSpecificStatesCallEnableBegin_; //< first time when
                                         // setSpecificStates should be called

  PyObject
      *pythonSetSpecificParametersFunction_; //< Python function handle that is
                                             // called to set parameters to the
                                             // CellML problem from the python
                                             // config
  PyObject
      *pythonSetSpecificStatesFunction_; //< Python function handle that is
                                         // called to set states to the CellML
                                         // problem from the python config
  PyObject
      *pythonHandleResultFunction_; //< Python function handle that is
                                    // called to process results from CellML
                                    // problem from the python config

  PyObject
      *pySetFunctionAdditionalParameter_; //< an additional python object that
                                          // will be passed as last argument to
                                          // the setParameters,
                                          // setSpecificParameters and
                                          // setSpecificStates callback function
  PyObject
      *pyHandleResultFunctionAdditionalParameter_; //< an additional python
                                                   // object that will be passed
                                                   // as last argument to the
                                                   // handleResult callback
                                                   // function

  PyObject *pyGlobalNaturalDofsList_; //< python list of global dof nos
};

#include "cellml/02_callback_handler.tpp"
