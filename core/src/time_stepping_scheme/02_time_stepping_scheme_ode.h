#pragma once

#include <Python.h> // has to be the first included header

#include "control/dihu_context.h"
#include "data_management/time_stepping/time_stepping.h"
#include "interfaces/discretizable_in_time.h"
#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "data_management/data.h"
#include "data_management/time_stepping/time_stepping.h"
#include "cellml/03_cellml_adapter.h"
#include "spatial_discretization/dirichlet_boundary_conditions/01_dirichlet_boundary_conditions.h"
#include "time_stepping_scheme/01_time_stepping_scheme_ode_base.h"

namespace TimeSteppingScheme {
template <typename DiscretizableInTimeType> class FullOdeDataForCheckpointing {
  typedef DiscretizableInTimeType DiscretizableInTime;
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;
  typedef typename DiscretizableInTimeType::FullData DiscretizableInTimeData;

  typedef typename Data::TimeStepping<FunctionSpace,
                                      DiscretizableInTimeType::nComponents()>
      Data; // type of Data object

public:
  FullOdeDataForCheckpointing(std::shared_ptr<Data> data,
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

/** This is the base class for all ode solvers.
 */
template <typename DiscretizableInTimeType>
class TimeSteppingSchemeOdeBaseDiscretizable
    : public TimeSteppingSchemeOdeBase<
          typename DiscretizableInTimeType::FunctionSpace,
          DiscretizableInTimeType::nComponents()> {
public:
  typedef DiscretizableInTimeType DiscretizableInTime;
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;
  typedef FullOdeDataForCheckpointing<DiscretizableInTimeType> FullData;

  //! constructor
  TimeSteppingSchemeOdeBaseDiscretizable(DihuContext context, std::string name);

  //! destructor
  virtual ~TimeSteppingSchemeOdeBaseDiscretizable() {}

  //! initialize discretizableInTime
  virtual void initialize();

  //! discretizable in time object
  DiscretizableInTimeType &discretizableInTime();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);

  //! reset state such that new initialization becomes necessary
  virtual void reset();

  //! object that stores Dirichlet boundary condition values
  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<
      FunctionSpace, DiscretizableInTimeType::nComponents()>>
  dirichletBoundaryConditions();

  FullData fullData();

protected:
  //! prepare the discretizableInTime object for the following call to
  //! getSlotConnectorData()
  virtual void prepareForGetSlotConnectorData() override;

  // int timeStepOutputInterval_;    //< time step number and time is output
  // every timeStepOutputInterval_ time steps
  DiscretizableInTimeType discretizableInTime_; //< the object to be discretized
  bool initialized_; //< if initialize() was already called

  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<
      FunctionSpace, DiscretizableInTimeType::nComponents()>>
      dirichletBoundaryConditions_; //< object that stores Dirichlet boundary
                                    // condition values
};

template <typename DiscretizableInTimeType>
class TimeSteppingSchemeOde
    : public TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType> {
public:
  //! use constructor of parent class
  using TimeSteppingSchemeOdeBaseDiscretizable<
      DiscretizableInTimeType>::TimeSteppingSchemeOdeBaseDiscretizable;

  // using
  // TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::initialize;
};
} // namespace TimeSteppingScheme

#include "time_stepping_scheme/02_time_stepping_scheme_ode.tpp"
