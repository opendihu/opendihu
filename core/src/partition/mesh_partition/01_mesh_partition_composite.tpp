#include "partition/mesh_partition/01_mesh_partition_composite.h"

namespace Partition
{

//! number of elements in the current partition
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nElementsLocal() const
{
  return nElementsLocal_;
}

//! number of elements in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nElementsGlobal() const
{
  return nElementsGlobal_;
}

//! number of dofs in the local partition
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nDofsLocalWithGhosts() const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();
  return nNodesLocalWithGhosts_ * nDofsPerNode;
}

//! number of dofs in the local partition, without ghosts
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nDofsLocalWithoutGhosts() const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();
  return nNodesLocalWithoutGhosts_ * nDofsPerNode;
}

//! number of dofs in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nDofsGlobal() const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();
  return nNodesGlobal_ * nDofsPerNode;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesLocalWithGhosts() const
{
  return nNodesLocalWithGhosts_;
}

//! number of nodes in the local partition
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesLocalWithoutGhosts() const
{
  return nNodesLocalWithoutGhosts_;
}

//! number of nodes in total
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
nNodesGlobal() const
{
  return nNodesGlobal_;
}

//! get the number of nodes in the global Petsc ordering that are in partitions prior to the own rank
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
beginNodeGlobalPetsc() const
{
  return nonDuplicateNodeNoGlobalBegin_;
}

//! get the local to global mapping for the current partition, for the dof numbering
template<int D, typename BasisFunctionType>
ISLocalToGlobalMapping MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
localToGlobalMappingDofs()
{
  return localToGlobalPetscMappingDofs_;
}

//! get the global natural element no for a local element no
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getElementNoGlobalNatural(element_no_t elementNoLocal) const
{
  // from a local element no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
  int subMeshNo = 0;
  element_no_t elementOnMeshNoLocal = 0;
  getSubMeshNoAndElementNoLocal(elementNoLocal, subMeshNo, elementOnMeshNoLocal);

  // add the global number on previous sub meshes
  global_no_t elementNoGlobalNatural = 0;
  for (int subMeshIndex = 0; subMeshIndex < subMeshNo; subMeshIndex++)
  {
    elementNoGlobalNatural += subFunctionSpaces_[subMeshIndex]->nElementsGlobal();
  }

  // add the elementNo in the current submesh
  elementNoGlobalNatural += subFunctionSpaces_[subMeshNo]->meshPartition()->getElementNoGlobalNatural(elementOnMeshNoLocal);

  return elementNoGlobalNatural;
}

//! get the global natural element no for a local element no
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoGlobalNatural(element_no_t elementNoLocal, int nodeIndex) const
{
  // from a local element no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
  int subMeshNo = 0;
  element_no_t elementOnMeshNoLocal = 0;
  getSubMeshNoAndElementNoLocal(elementNoLocal, subMeshNo, elementOnMeshNoLocal);

  // get node no on subMesh
  node_no_t nodeNoLocal = subFunctionSpaces_[subMeshNo]->getNodeNo(elementOnMeshNoLocal, nodeIndex);

  // get global natural no on subMesh
  std::array<global_no_t,D> coordinatesGlobal = subFunctionSpaces_[subMeshNo]->meshPartition()->getCoordinatesGlobal(nodeNoLocal);
  global_no_t nodeNoGlobalNaturalSubMesh = subFunctionSpaces_[subMeshNo]->meshPartition()->getNodeNoGlobalNatural(coordinatesGlobal);

  // add the global number on previous sub meshes
  global_no_t nodeNoGlobalNatural = 0;
  for (int subMeshIndex = 0; subMeshIndex < subMeshNo; subMeshIndex++)
  {
    nodeNoGlobalNatural += subFunctionSpaces_[subMeshIndex]->nNodesGlobal();
  }

  // add the nodeNoGlobalNatural in the current submesh
  nodeNoGlobalNatural += nodeNoGlobalNaturalSubMesh;

  return nodeNoGlobalNatural;
}

//! get the node no in global petsc ordering from a local node no
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoGlobalPetsc(node_no_t nodeNoLocal) const
{
  std::pair<int,node_no_t> result = nodeNoNonDuplicateLocalToMeshAndDuplicateLocal_[nodeNoLocal];
  int subMeshNo = result.first;
  nodeNoLocal = result.second;

  return meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoLocal];
}

//! transfer the local nos in global dof nos, using the PETSc localToGlobal mapping for the dofs
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getDofNoGlobalPetsc(const std::vector<dof_no_t> &dofNosLocal, std::vector<PetscInt> &dofNosGlobalPetsc) const
{
  dofNosGlobalPetsc.resize(dofNosLocal.size());

  // transfer the local indices to global indices
  PetscErrorCode ierr;
  ierr = ISLocalToGlobalMappingApply(localToGlobalPetscMappingDofs_, dofNosLocal.size(), dofNosLocal.data(), dofNosGlobalPetsc.data()); CHKERRV(ierr);
}

//! get the global petsc dof no for the local no, using the PETSc localToGlobal mapping for the dofs
template<int D, typename BasisFunctionType>
global_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getDofNoGlobalPetsc(dof_no_t dofNoLocal) const
{
  PetscInt dofNoGlobal;
  PetscErrorCode ierr;
  ierr = ISLocalToGlobalMappingApply(localToGlobalPetscMappingDofs_, 1, &dofNoLocal, &dofNoGlobal); CHKERRQ(ierr);
  return (global_no_t)dofNoGlobal;
}

//! get the local element no. from the global no., set isOnLocalDomain to true if the node with global coordinates is in the local domain
template<int D, typename BasisFunctionType>
element_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getElementNoLocal(global_no_t elementNoGlobalPetsc, bool &isOnLocalDomain) const
{
  if (elementNoGlobalBegin_ <= elementNoGlobalPetsc && elementNoGlobalPetsc < elementNoGlobalBegin_+nElementsLocal_)
  {
    isOnLocalDomain = true;
    return elementNoGlobalPetsc - elementNoGlobalBegin_;
  }
  else
  {
    isOnLocalDomain = false;
    return -1;
  }
}

//! get the local node no for a global petsc node no, does not work for ghost nodes
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoLocal(global_no_t nodeNoGlobalPetsc, bool &isLocal) const
{
  node_no_t nodeNoLocal = (node_no_t)(nodeNoGlobalPetsc - nonDuplicateNodeNoGlobalBegin_);
  isLocal = (nodeNoLocal >= 0) && (nodeNoLocal < nDofsLocalWithoutGhosts());
  return nodeNoLocal;
}

//! get the local dof no for a global petsc dof no, does not work for ghost nodes
template<int D, typename BasisFunctionType>
dof_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getDofNoLocal(global_no_t dofNoGlobalPetsc, bool &isLocal) const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();
  global_no_t nodeNoGlobalPetsc = dofNoGlobalPetsc / nDofsPerNode;
  int nodalDofIndex = dofNoGlobalPetsc % nDofsPerNode;
  return getNodeNoLocal(nodeNoGlobalPetsc, isLocal) * nDofsPerNode + nodalDofIndex;
}

//! from a vector of values of global/natural node numbers remove all that are non-local, nComponents consecutive values for each dof are assumed
template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
extractLocalNodesWithoutGhosts(std::vector<T> &values, int nComponents) const
{
  std::vector<T> result(nNodesLocalWithoutGhosts() + nNodesSharedLocal_);
  global_no_t resultIndex = 0;

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    int nNodesLocal = subFunctionSpaces_[subMeshNo]->nNodesLocalWithoutGhosts();
    std::vector<T> vectorSubMesh(values.begin()+resultIndex, values.begin()+resultIndex+nNodesLocal);

    subFunctionSpaces_[subMeshNo]->extractLocalNodesWithoutGhosts(vectorSubMesh, nComponents);

    std::copy(vectorSubMesh.begin(), vectorSubMesh.end(), result.begin()+resultIndex);
    resultIndex += nNodesLocal;
  }

  values.clear();

  // only copy non-shared nodes
  resultIndex = 0;

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    // loop over duplicate nodes of submesh
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < subFunctionSpaces_[subMeshNo]->nNodesLocalWithoutGhosts(); nodeNoLocal++, resultIndex++)
    {
      if (meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo][nodeNoLocal] != -1)
      {
        values.push_back(result[resultIndex]);
      }
    }
  }
}

//! from a vector of values of global/natural dofs remove all that are non-local
template<int D, typename BasisFunctionType>
template <typename T>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
extractLocalDofsWithoutGhosts(std::vector<T> &values) const
{
  const int nDofsPerNode = FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>::nDofsPerNode();

  std::vector<T> result(nDofsLocalWithoutGhosts() + nNodesSharedLocal_*nDofsPerNode);
  global_no_t resultIndex = 0;

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    int nDofsLocal = subFunctionSpaces_[subMeshNo]->nDofsLocalWithoutGhosts();
    std::vector<T> vectorSubMesh(values.begin()+resultIndex, values.begin()+resultIndex+nDofsLocal);

    subFunctionSpaces_[subMeshNo]->meshPartition()->extractLocalDofsWithoutGhosts(vectorSubMesh);

    std::copy(vectorSubMesh.begin(), vectorSubMesh.end(), result.begin()+resultIndex);
    resultIndex += nDofsLocal;
  }

  values.clear();

  // only copy non-shared dofs
  resultIndex = 0;

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    // loop over duplicate nodes of submesh
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < subFunctionSpaces_[subMeshNo]->nDofsLocalWithoutGhosts(); nodeNoLocal++, resultIndex+=nDofsPerNode)
    {
      if (meshAndNodeNoLocalToNodeNoNonDuplicateLocal_[subMeshNo][nodeNoLocal] != -1)
      {
        for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
        {
          values.push_back(result[resultIndex + nodalDofIndex]);
        }
      }
    }
  }
}

template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
extractLocalDofsWithoutGhosts(std::vector<double> &vector) const
{
  this->template extractLocalDofsWithoutGhosts<double>(vector);
}

//! output to stream for debugging
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
output(std::ostream &stream)
{
  stream << "CompositeMesh ";
  for (int subMeshNo = 0; subMeshNo < nSubMeshes_; subMeshNo++)
  {
    stream << "subMesh " << subMeshNo << "/" << nSubMeshes_ << ": ";
    subFunctionSpaces_[subMeshNo]->meshPartition()->output(stream);
  }
}

//! get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the dofs without ghost dofs, the whole vector are the dofs with ghost dofs
//! @param onlyNodalValues: if for Hermite only get every second dof such that derivatives are not returned
template<int D, typename BasisFunctionType>
const std::vector<PetscInt> &MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
dofNosLocal(bool onlyNodalValues) const
{
  if (onlyNodalValues)
  {
    return onlyNodalDofLocalNos_;
  }
  else
  {
    return this->dofNosLocal_;
  }
}

// use getDofNoGlobalPetsc(dofNosLocal(), ...) to get dofNosGlobalPetsc

//! get the global dof nos of the ghost dofs in the local partition
template<int D, typename BasisFunctionType>
const std::vector<PetscInt> &MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
ghostDofNosGlobalPetsc() const
{
  return ghostDofNosGlobalPetsc_;
}

//! check if the given dof is owned by the own rank, then return true, if not, neighbourRankNo is set to the rank by which the dof is owned
template<int D, typename BasisFunctionType>
bool MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
isNonGhost(node_no_t nodeNoLocal, int &neighbourRankNo) const
{
  LOG(FATAL) << "\"isNonGhost\" is not implemented for composite mesh.";
  return true;
}

//! get information about neighbouring rank and boundary elements for specified face,
//! @param neighbourRankNo: the rank of the neighbouring process that shares the face, @param nElements: Size of one-layer mesh that contains boundary elements that touch the neighbouring process
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getBoundaryElements(Mesh::face_t face, int &neighbourRankNo, std::array<element_no_t,D> &nBoundaryElements, std::vector<dof_no_t> &dofNos)
{
  LOG(FATAL) << "\"getBoundaryElements\" is not implemented for composite mesh.";
}

//! get the rank no of the neighbour in direction face, -1 if there is no such neighbour
template<int D, typename BasisFunctionType>
int MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
neighbourRank(Mesh::face_t face)
{
  LOG(FATAL) << "\"getBoundaryElements\" is not implemented for composite mesh.";
  return 0;
}

//! from a local node no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getSubMeshNoAndNodeNoLocal(node_no_t nodeNoLocal, int &subMeshNo, node_no_t &nodeOnMeshNoLocal)
{
  // loop over submeshes
  node_no_t nNodesPreviousSubMeshes = 0;
  for (int subMeshIndex = 0; subMeshIndex < nSubMeshes_; subMeshIndex++)
  {
    node_no_t nNodesCurrentSubMesh = nNonDuplicateNodesWithoutGhosts_[subMeshIndex];

    if (nodeNoLocal < nNodesPreviousSubMeshes + nNodesCurrentSubMesh)
    {
      subMeshNo = subMeshIndex;
      nodeOnMeshNoLocal = nodeNoLocal - nNodesPreviousSubMeshes;
      break;
    }

    nNodesPreviousSubMeshes += nNodesCurrentSubMesh;
  }

  assert (nNodesPreviousSubMeshes < nNodesLocalWithoutGhosts_);
}

//! from a local element no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getSubMeshNoAndElementNoLocal(element_no_t elementNoLocal, int &subMeshNo, element_no_t &elementOnMeshNoLocal) const
{
  element_no_t nElementsPreviousSubMeshes = 0;

  // loop over all meshes
  for (int subMeshIndex = 0; subMeshIndex < nSubMeshes_; subMeshIndex++)
  {
    element_no_t nElementsLocal = subFunctionSpaces_[subMeshIndex]->nElementsLocal();

    if (elementNoLocal >= nElementsPreviousSubMeshes && elementNoLocal < nElementsPreviousSubMeshes + nElementsLocal)
    {
      subMeshNo = subMeshIndex;
      elementOnMeshNoLocal = elementNoLocal - nElementsPreviousSubMeshes;
      return;
    }
    nElementsPreviousSubMeshes += nElementsLocal;
  }
}

//! from a local node no in the composite numbering get the subMeshNo and the no in the submesh-based numbering
template<int D, typename BasisFunctionType>
void MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getSubMeshesWithNodes(node_no_t nodeNoLocal, std::vector<std::pair<int,node_no_t>> &subMeshesWithNodes)
{
  int subMeshNo = 0;
  node_no_t nodeOnMeshNoLocal = 0;
  getSubMeshNoAndNodeNoLocal(nodeNoLocal, subMeshNo, nodeOnMeshNoLocal);

  subMeshesWithNodes.push_back(std::make_pair(subMeshNo, nodeOnMeshNoLocal));

  // loop over all meshes and shared nodes
  for (int subMeshIndex = 0; subMeshIndex < nSubMeshes_; subMeshIndex++)
  {
    for (std::map<node_no_t,std::pair<int,node_no_t>>::iterator iter = removedSharedNodes_[subMeshIndex]; iter != removedSharedNodes_[subMeshIndex].end(); iter++)
    {
      int meshNo = iter->first;
      int nodeNo = iter->second;

      if (meshNo == subMeshNo && nodeNo == nodeOnMeshNoLocal)
      {
        subMeshesWithNodes.push_back(std::make_pair(meshNo, nodeNo));
      }
    }
  }
}


//! from the submesh no and the local node no in the submesh numbering get the local node no in the composite numbering
template<int D, typename BasisFunctionType>
node_no_t MeshPartition<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,Mesh::CompositeOfDimension<D>>::
getNodeNoLocalFromSubmesh(int subMeshNo, int nodeNoDuplicateOnSubmesh)
{
  assert(subMeshNo >= 0 && subMeshNo < meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_.size());
  assert(nodeNoDuplicateOnSubmesh >= 0 && nodeNoDuplicateOnSubmesh < meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo].size());

  return meshAndNodeNoLocalToNodeNoNonDuplicateGlobal_[subMeshNo][nodeNoDuplicateOnSubmesh];
}


}  // namespace

