#include "function_space/03_function_space_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace FunctionSpace
{
 
template<typename MeshType,typename BasisFunctionType>
void FunctionSpacePartition<MeshType,BasisFunctionType>::
initialize()
{
  // if meshPartition was already created earlier, do nothing
  if (this->meshPartition_) 
  {
    LOG(DEBUG) << "in FunctionSpacePartition<structured>: meshPartition already set";

    // Set number of global and local elements and number of ranks to be consistent with the values stored in meshPartition.
    // These local values of function space are the buffer to which the values from python config are parsed. This data duplication is needed because when the config is parsed there is no meshPartition yet.
    // The values are used e.g. in the python output writer.
    // At all other places, the values stored in mesh partition are used.
    for (int i = 0; i < MeshType::dim(); i++)
    {
      this->nElementsPerCoordinateDirectionGlobal_[i] = this->meshPartition_->nElementsGlobal(i);
      this->nElementsPerCoordinateDirectionLocal_[i] = this->meshPartition_->nElementsLocal(i);
      this->nRanks_[i] = this->meshPartition_->nRanks(i);
    }

    return;
  }
  
  // Creation of the partitioning is only possible after the number of elements is known.
  // Because this may need file I/O (e.g. reading from exfiles)
 
  LOG(DEBUG) << "FunctionSpacePartition<Structured>::initialize(), create meshPartition";
  
  // create partitioning
  assert(this->partitionManager_ != nullptr);
  
  // get whether the mesh information in config specifies local or global domain, when this function space is not generated by settings
  // but directly by a constructor with node positions (e.g. meshManager->createFunctionSpace) then local partitioning may be necessary
  bool inputMeshIsGlobal = true;
  if (this->forcePartitioningCreationFromLocalNumberOfElements_)
  {
    inputMeshIsGlobal = false;
    LOG(DEBUG) << "set inputMeshIsGlobal = false";
  }
  else
  {
    inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
    LOG(DEBUG) << "got inputMeshIsGlobal = " << std::boolalpha << inputMeshIsGlobal << " from config";
  }

  if (inputMeshIsGlobal)
  {
    this->meshPartition_ = this->partitionManager_->template createPartitioningStructuredGlobal<FunctionSpace<MeshType,BasisFunctionType>>(
      this->nElementsPerCoordinateDirectionGlobal_, this->nElementsPerCoordinateDirectionLocal_, this->nRanks_);
  }
  else 
  {
    this->meshPartition_ = this->partitionManager_->template createPartitioningStructuredLocal<FunctionSpace<MeshType,BasisFunctionType>>(
      this->nElementsPerCoordinateDirectionGlobal_, this->nElementsPerCoordinateDirectionLocal_, this->nRanks_);
  }
  assert(this->meshPartition_);
}

} // namespace
