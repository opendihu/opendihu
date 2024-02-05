#include "data_management/time_stepping/time_stepping.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data {

template <typename FunctionSpaceType, typename BaseDataType>
StreamlineTracer<FunctionSpaceType, BaseDataType>::StreamlineTracer(
    DihuContext context)
    : Data<FunctionSpaceType>(context), fiberNo_(0) {}

template <typename FunctionSpaceType, typename BaseDataType>
StreamlineTracer<FunctionSpaceType, BaseDataType>::~StreamlineTracer() {
  // free PETSc objects
  if (this->initialized_) {
    // PetscErrorCode ierr;
    // ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template <typename FunctionSpaceType, typename BaseDataType>
void StreamlineTracer<FunctionSpaceType, BaseDataType>::setBaseData(
    std::shared_ptr<BaseDataType> baseData) {
  baseData_ = baseData;

  // set function space
  this->setFunctionSpace(baseData_->functionSpace());
}

template <typename FunctionSpaceType, typename BaseDataType>
void StreamlineTracer<FunctionSpaceType, BaseDataType>::createPetscObjects() {
  LOG(DEBUG) << "StreamlineTracer<FunctionSpaceType,BaseDataType>::"
                "createPetscObjects()"
             << std::endl;
  assert(this->functionSpace_);

  // create field variables on local partition
  this->gradient_ =
      this->functionSpace_->template createFieldVariable<3>("gradient");
}

template <typename FunctionSpaceType, typename BaseDataType>
void StreamlineTracer<FunctionSpaceType, BaseDataType>::createFiberMesh(
    const std::vector<Vec3> &nodePositions) {
  std::shared_ptr<FunctionSpaceFiber> meshPtr;

  // create name for fiber mesh
  std::stringstream name;
  name << "fiber" << std::setw(5) << std::setfill('0') << fiberNo_;
  fiberNo_++;

  // set number of elements. We have 1D linear Lagrange elements, i.e. 2 nodes
  // per element
  const int nElements = nodePositions.size() - 1;
  std::array<element_no_t, 1> nElementsPerCoordinateDirection{nElements};

  // create mesh by meshManager
  std::array<int, 1> nRanks({1}); // this is completely serial
  meshPtr = this->context_.meshManager()
                ->template createFunctionSpace<FunctionSpaceFiber>(
                    name.str(), nodePositions, nElementsPerCoordinateDirection,
                    nRanks);

  // get geometry field
  std::shared_ptr<FieldVariableFiberGeometry> geometryField =
      std::make_shared<FieldVariableFiberGeometry>(meshPtr->geometryField());

  // add geometry field
  this->fiberGeometry_.push_back(geometryField);
}

template <typename FunctionSpaceType, typename BaseDataType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
StreamlineTracer<FunctionSpaceType, BaseDataType>::gradient() {
  return this->gradient_;
}

template <typename FunctionSpaceType, typename BaseDataType>
dof_no_t
StreamlineTracer<FunctionSpaceType, BaseDataType>::nNodesLocalWithGhosts() {
  return this->functionSpace_->nNodesLocalWithGhosts();
}

template <typename FunctionSpaceType, typename BaseDataType>
dof_no_t
StreamlineTracer<FunctionSpaceType, BaseDataType>::nNodesLocalWithoutGhosts() {
  return this->functionSpace_->nNodesLocalWithoutGhosts();
}

template <typename FunctionSpaceType, typename BaseDataType>
void StreamlineTracer<FunctionSpaceType, BaseDataType>::print() {
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->gradient_;
  VLOG(4) << "======================";
}

template <typename FunctionSpaceType, typename BaseDataType>
typename StreamlineTracer<FunctionSpaceType,
                          BaseDataType>::FieldVariablesForOutputWriter
StreamlineTracer<FunctionSpaceType,
                 BaseDataType>::getFieldVariablesForOutputWriter() {
  return std::tuple_cat(
      baseData_->getFieldVariablesForOutputWriter(),
      std::tuple<
          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>>(
          gradient_),
      std::tuple<std::vector<std::shared_ptr<FieldVariableFiberGeometry>>>(
          fiberGeometry_));
}

} // namespace Data
