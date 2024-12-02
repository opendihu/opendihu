#pragma once

#include <Python.h> // has to be the first included header
#include "time_stepping_scheme/02_time_stepping_scheme_ode.h"
#include "interfaces/runnable.h"
#include "data_management/time_stepping/time_stepping_implicit.h"
#include "control/dihu_context.h"
#include "solver/linear.h"

namespace TimeSteppingScheme {

template <typename DiscretizableInTimeType>
class FullImplicitDataForCheckpointing {
  typedef DiscretizableInTimeType DiscretizableInTime;
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;
  typedef typename DiscretizableInTimeType::FullData DiscretizableInTimeData;

  typedef typename Data::TimeSteppingImplicit<
      FunctionSpace, DiscretizableInTimeType::nComponents()>
      Data;

public:
  FullImplicitDataForCheckpointing(
      std::shared_ptr<Data> data,
      DiscretizableInTimeData &discretizableInTimeData);

  //! field variables that will be output by checkpointing
  typedef decltype(std::tuple_cat(
      std::declval<typename Data::FieldVariablesForCheckpointing>(),
      std::declval<
          typename DiscretizableInTimeData::FieldVariablesForCheckpointing>()))
      FieldVariablesForCheckpointing;

  //! get pointers to all field variables that can be written by checkpointing
  FieldVariablesForCheckpointing getFieldVariablesForCheckpointing();

  //! field variables that will be output by checkpointing
  typedef FieldVariablesForCheckpointing FieldVariablesForOutputWriter;

  //! Not needed for this implementation, shadowing checkpointing function
  FieldVariablesForCheckpointing getFieldVariablesForOutputWriter();

  bool restoreState(const InputReader::Generic &r);

  const std::shared_ptr<FunctionSpace> functionSpace() const;

private:
  std::shared_ptr<Data> data_; //< data object
  DiscretizableInTimeData
      &discretizableInTimeData_; //< the object to be discretized
};

/** The implicit time integration scheme
 */
template <typename DiscretizableInTimeType>
class TimeSteppingImplicit
    : public TimeSteppingSchemeOde<DiscretizableInTimeType> {
public:
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;
  typedef Data::TimeSteppingImplicit<
      typename DiscretizableInTimeType::FunctionSpace,
      DiscretizableInTimeType::nComponents()>
      DataImplicit;
  typedef typename DataImplicit::SlotConnectorDataType SlotConnectorDataType;
  typedef FullImplicitDataForCheckpointing<DiscretizableInTimeType> FullData;

  //! constructor
  TimeSteppingImplicit(DihuContext context, const std::string name);

  //! advance simulation by the given time span [startTime_, endTime_] with
  //! given numberTimeSteps, data in solution is used, afterwards new data is in
  //! solution
  virtual void advanceTimeSpan(
      bool withOutputWritersEnabled = true,
      std::shared_ptr<Checkpointing::Handle> checkpointing = nullptr) = 0;

  virtual void initialize();

  //! setup the system matrix and other things that depend on the step width
  void initializeWithTimeStepWidth(double timeStepWidth);

  //! reset the object's state, i.e. delete the linear solver
  virtual void reset();

  //! data for implicit timestepping
  DataImplicit &dataImplicit();

  //! output the given data for debugging
  // virtual std::string getString(typename
  // TimeSteppingSchemeOde<DiscretizableInTimeType>::SlotConnectorDataType
  // &data);

  //! there should be no getSlotConnectorData here!

  FullData fullData();

protected:
  //! actual implementation of initializeWithTimeStepWidth
  virtual void initializeWithTimeStepWidth_impl(double timeStepWidth);

  //! precomputes the integration matrix for example A = (I-dtM^(-1)K) for the
  //! implicit euler scheme
  virtual void setSystemMatrix(double timeStepWidth) = 0;

  //! initialize the linear solve that is needed for the solution of the
  //! implicit timestepping system
  void initializeLinearSolver();

  //! solves the linear system of equations resulting from the Implicit Euler
  //! method time discretization
  void solveLinearSystem(Vec &input, Vec &output);

  std::shared_ptr<Data::TimeSteppingImplicit<
      typename DiscretizableInTimeType::FunctionSpace,
      DiscretizableInTimeType::nComponents()>>
      dataImplicit_; //< a pointer to the data_ object but of type
                     // Data::TimeSteppingImplicit
  std::shared_ptr<Solver::Linear>
      linearSolver_;         //< the linear solver used for solving the system
  std::shared_ptr<KSP> ksp_; //< the ksp object of the linear solver

  double initializedTimeStepWidth_ =
      -1.0; //< the time step width that was used for the initialization, or
            // negative if the step width has not been initialized
  double timeStepWidthRelativeTolerance_;  //< tolerance for the time step width
                                           // to rebuild the system matrix and
                                           // integrationMatrixRHS
  std::string durationInitTimeStepLogKey_; //< log key for the duration of the
                                           //(re)initialization of the system
                                           // matrix and integrationMatrixRHS
};

} // namespace TimeSteppingScheme

#include "time_stepping_scheme/03_time_stepping_implicit.tpp"
