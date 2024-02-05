#pragma once

#include <Python.h> // has to be the first included header
//#include "time_stepping_scheme/02_time_stepping_scheme_ode.h"
#include "time_stepping_scheme/03_time_stepping_implicit.h"
//#include "interfaces/runnable.h"
//#include "data_management/time_stepping/time_stepping.h"
#include "control/dihu_context.h"

namespace TimeSteppingScheme {

/** The implicit Euler integration scheme (backward Euler), u_{t+1} = u_{t} +
 * dt*f(t+1)
 */
template <typename DiscretizableInTimeType>
class ImplicitEuler : public TimeSteppingImplicit<DiscretizableInTimeType> {
public:
  //! constructor
  ImplicitEuler(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_] with
  //! given numberTimeSteps, data in solution is used, afterwards new data is in
  //! solution
  void advanceTimeSpan(bool withOutputWritersEnabled = true);

  //! precomputes the integration matrix A=I-dtM^(-1)K for the implicit euler
  //! scheme
  void setSystemMatrix(double timeStepWidth);

  //! call the output writer on the data object, output files will contain
  //! currentTime, with callCountIncrement !=1 output timesteps can be skipped
  virtual void callOutputWriter(int timeStepNo, double currentTime,
                                int callCountIncrement = 1) override;

protected:
};

} // namespace TimeSteppingScheme

#include "time_stepping_scheme/implicit_euler.tpp"
