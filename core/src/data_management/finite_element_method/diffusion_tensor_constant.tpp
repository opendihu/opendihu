#include "data_management/finite_element_method/diffusion_tensor_constant.h"

namespace Data {

template <typename FunctionSpaceType>
DiffusionTensorConstant<FunctionSpaceType>::DiffusionTensorConstant(
    PythonConfig settings)
    : DiffusionTensorBase<FunctionSpaceType>::DiffusionTensorBase(settings) {
  LOG(DEBUG) << "construct Diffusion tensor";
  if (VLOG_IS_ON(1)) {
    PythonUtility::printDict(settings.pyObject());
  }
}

template <typename FunctionSpaceType>
void DiffusionTensorConstant<FunctionSpaceType>::initialize(
    std::shared_ptr<FunctionSpaceType> functionSpace) {
  LOG(DEBUG) << "DiffusionTensorConstant::initialize";
  this->dataFunctionSpace_ = functionSpace;

  // initialize diffusion tensor
  const int D = FunctionSpaceType::dim();

  // create identity matrix as default values
  MathUtility::Matrix<D, D> defaultValue({0});

  for (int i = 0; i < D; i++)
    defaultValue(i, i) = 1.0;

  // parse diffusion tensor
  this->diffusionTensor_.initialize(this->specificSettings_, "diffusionTensor",
                                    defaultValue, this->dataFunctionSpace_);
}

template <typename FunctionSpaceType>
template <typename double_v_t, typename element_no_v_t>
const MathUtility::Matrix<FunctionSpaceType::dim(), FunctionSpaceType::dim(),
                          double_v_t>
DiffusionTensorConstant<FunctionSpaceType>::diffusionTensor(
    element_no_v_t elementNoLocal,
    const std::array<double, FunctionSpaceType::dim()> xi) const {
  return this->diffusionTensor_.value(elementNoLocal);
}

} // namespace Data
