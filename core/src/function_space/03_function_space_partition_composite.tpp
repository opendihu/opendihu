#include "function_space/03_function_space_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace FunctionSpace
{
 
// forward declaration
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;

template<int D,typename BasisFunctionType>
void FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
initialize()
{
  PyObject *pyObject = this->specificSettings.pyObject();

  int nSubmeshes = 1;
  if (PyList_Check(pyObject))
  {
    nSubmeshes = PyList_Size(pyObject);
  }
  LOG(DEBUG) << "initialize " << nSubmeshes << " submeshes";
  subFunctionSpaces_.resize(nSubmeshes);

  // created sub meshes
  for (int i = 0; i < nSubmeshes; i++)
  {
    subFunctionSpaces_[i] = std::make_shared<Mesh::StructuredDeformableOfDimension<D>>(PythonConfig(this->specificSettings,i));
  }

  // Creation of the partitioning is only possible after the number of elements is known.
  // Because this may need file I/O (e.g. reading from exfiles)
 
  // This initialize() method is called from parseFromSettings(), 04_function_space_data_unstructured_exfile_io.tpp
  
  VLOG(1) << "FunctionSpacePartition<Composite>::initialize()";

  // initialize number of local and global elements
  this->nElementsLocal_ = 0;
  this->nElementsGlobal_ = 0;

  typedef FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> FunctionSpaceType;

  // iterate over submeshes
  for(std::shared_ptr<FunctionSpaceType> &subFunctionSpace : subFunctionSpaces_)
  {
    subFunctionSpace->initialize();

    // count number of elements
    this->nElementsLocal_ += subFunctionSpace->nElementsLocal();
    this->nElementsGlobal_ += subFunctionSpace->nElementsGlobal();
  }

  // create partitioning
  assert(this->partitionManager_ != nullptr);
  this->meshPartition_ = this->partitionManager_->template createPartitioningComposite<FunctionSpaceType>(subFunctionSpaces_);

  assert(this->meshPartition_);

  // set initalized_ to true which indicates that initialize has been called
  this->initialized_ = true;
}

} // namespace