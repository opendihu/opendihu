#include "time_stepping_scheme/full_heun_data.h"

namespace TimeSteppingScheme {
template <typename DiscretizableInTimeType>
FullHeunDataForCheckpointing<DiscretizableInTimeType>::
    FullHeunDataForCheckpointing(
        std::shared_ptr<Data> data,
        DiscretizableInTimeData &discretizableInTimeData)
    : data_(data), discretizableInTimeData_(discretizableInTimeData) {}

template <typename DiscretizableInTimeType>
typename FullHeunDataForCheckpointing<
    DiscretizableInTimeType>::FieldVariablesForCheckpointing
FullHeunDataForCheckpointing<
    DiscretizableInTimeType>::getFieldVariablesForCheckpointing() {
  return std::tuple_cat(
      data_->getFieldVariablesForCheckpointing(),
      discretizableInTimeData_.getFieldVariablesForCheckpointing());
}

template <typename DiscretizableInTimeType>
typename FullHeunDataForCheckpointing<
    DiscretizableInTimeType>::FieldVariablesForOutputWriter
FullHeunDataForCheckpointing<
    DiscretizableInTimeType>::getFieldVariablesForOutputWriter() {
  return getFieldVariablesForCheckpointing();
}

template <typename DiscretizableInTimeType>
bool FullHeunDataForCheckpointing<DiscretizableInTimeType>::restoreState(
    const InputReader::Generic &r) {
  return data_->restoreState(r) && discretizableInTimeData_.restoreState(r);
}

template <typename DiscretizableInTimeType>
const std::shared_ptr<typename FullHeunDataForCheckpointing<
    DiscretizableInTimeType>::FunctionSpace>
FullHeunDataForCheckpointing<DiscretizableInTimeType>::functionSpace() const {
  return data_->functionSpace();
}
} // namespace TimeSteppingScheme
