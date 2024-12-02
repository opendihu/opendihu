#pragma once

#include "time_stepping_scheme/03_time_stepping_explicit.h"
#include "interfaces/runnable.h"
#include "data_management/time_stepping/time_stepping_heun.h"
#include "control/dihu_context.h"

namespace TimeSteppingScheme {

template <typename DiscretizableInTimeType> class FullHeunDataForCheckpointing {
  typedef DiscretizableInTimeType DiscretizableInTime;
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;
  typedef typename DiscretizableInTimeType::FullData DiscretizableInTimeData;

  typedef
      typename Data::TimeSteppingHeun<FunctionSpace,
                                      DiscretizableInTimeType::nComponents()>
          Data;

public:
  FullHeunDataForCheckpointing(
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

} // namespace TimeSteppingScheme

#include "time_stepping_scheme/full_heun_data.tpp"
