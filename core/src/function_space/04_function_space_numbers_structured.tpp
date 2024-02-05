#include "function_space/04_function_space_numbers_structured.h"

#include <Python.h> // has to be the first included header
#include <cmath>
#include <array>

#include "easylogging++.h"

namespace FunctionSpace {

// These numberings run over the whole locally stored information, i.e. there is
// no distinguishing between ghost and interior nodes here.

// element-local dofIndex to local dofNo for 1D
template <typename MeshType, typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                              Mesh::isStructuredWithDim<1, MeshType>>::
    getDofNo(element_no_t elementNoLocal, int dofIndex) const {
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0,1 2,3
  // averageNDofsPerElement:
  // 1         2            2
  // VLOG(3) << "getDofNo<1D>(elementNoLocal=" << elementNoLocal << ",
  // dofIndex=" << dofIndex << ") = " <<
  // FunctionSpaceFunction<MeshType,BasisFunctionType>::averageNDofsPerElement()
  // * elementNoLocal + dofIndex;

  return FunctionSpaceFunction<MeshType,
                               BasisFunctionType>::averageNDofsPerElement() *
             elementNoLocal +
         dofIndex;
}

//! get all dofs of a specific node for 1D
template <typename MeshType, typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType, BasisFunctionType,
                          Mesh::isStructuredWithDim<1, MeshType>>::
    getNodeDofs(node_no_t nodeGlobalNo,
                std::vector<dof_no_t> &dofGlobalNos) const {
  dofGlobalNos.reserve(
      dofGlobalNos.size() +
      FunctionSpaceBaseDim<1, BasisFunctionType>::nDofsPerNode());
  for (int dofIndex = 0;
       dofIndex < FunctionSpaceBaseDim<1, BasisFunctionType>::nDofsPerNode();
       dofIndex++) {
    dofGlobalNos.push_back(getNodeDofNo(nodeGlobalNo, dofIndex));
  }
}

//! get all dofs of a specific node for 1D
template <typename MeshType, typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType, BasisFunctionType,
                          Mesh::isStructuredWithDim<1, MeshType>>::
    getNodeDofs(
        node_no_t nodeGlobalNo,
        std::array<dof_no_t,
                   FunctionSpaceBaseDim<1, BasisFunctionType>::nDofsPerNode()>
            &dofGlobalNos) const {
  for (int dofIndex = 0;
       dofIndex < FunctionSpaceBaseDim<1, BasisFunctionType>::nDofsPerNode();
       dofIndex++) {
    dofGlobalNos[dofIndex] = getNodeDofNo(nodeGlobalNo, dofIndex);
  }
}

//! get the dof no of the specified dof at the node for 1D
template <typename MeshType, typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                              Mesh::isStructuredWithDim<1, MeshType>>::
    getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const {
  return FunctionSpaceBaseDim<1, BasisFunctionType>::nDofsPerNode() *
             nodeGlobalNo +
         dofIndex;
}

// element-local dofIndex to local dofNo for 2D
template <typename MeshType, typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                              Mesh::isStructuredWithDim<2, MeshType>>::
    getDofNo(element_no_t elementNoLocal, int dofIndex) const {
  const int nDofsPerNode = this->nDofsPerNode();
  int nodeIndex = dofIndex / nDofsPerNode;
  int dofOnNodeIndex = dofIndex % nDofsPerNode;

  return getNodeNo(elementNoLocal, nodeIndex) * nDofsPerNode + dofOnNodeIndex;
}

//! get all dofs of a specific node for 2D
template <typename MeshType, typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType, BasisFunctionType,
                          Mesh::isStructuredWithDim<2, MeshType>>::
    getNodeDofs(node_no_t nodeGlobalNo,
                std::vector<dof_no_t> &dofGlobalNos) const {
  dofGlobalNos.reserve(
      dofGlobalNos.size() +
      FunctionSpaceBaseDim<2, BasisFunctionType>::nDofsPerNode());
  for (int dofIndex = 0;
       dofIndex < FunctionSpaceBaseDim<2, BasisFunctionType>::nDofsPerNode();
       dofIndex++) {
    dofGlobalNos.push_back(getNodeDofNo(nodeGlobalNo, dofIndex));
  }
}

//! get all dofs of a specific node for 2D
template <typename MeshType, typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType, BasisFunctionType,
                          Mesh::isStructuredWithDim<2, MeshType>>::
    getNodeDofs(
        node_no_t nodeGlobalNo,
        std::array<dof_no_t,
                   FunctionSpaceBaseDim<2, BasisFunctionType>::nDofsPerNode()>
            &dofGlobalNos) const {
  for (int dofIndex = 0;
       dofIndex < FunctionSpaceBaseDim<2, BasisFunctionType>::nDofsPerNode();
       dofIndex++) {
    dofGlobalNos[dofIndex] = getNodeDofNo(nodeGlobalNo, dofIndex);
  }
}

//! get the dof no of the specified dof at the node for 2D
template <typename MeshType, typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                              Mesh::isStructuredWithDim<2, MeshType>>::
    getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const {
  return FunctionSpaceBaseDim<2, BasisFunctionType>::nDofsPerNode() *
             nodeGlobalNo +
         dofIndex;
}

// element-local dofIndex to local dofNo for 3D
template <typename MeshType, typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                              Mesh::isStructuredWithDim<3, MeshType>>::
    getDofNo(element_no_t elementNoLocal, int dofIndex) const {
  const int nDofsPerNode = this->nDofsPerNode();
  int nodeIndex = dofIndex / nDofsPerNode;
  int dofOnNodeIndex = dofIndex % nDofsPerNode;

  return getNodeNo(elementNoLocal, nodeIndex) * nDofsPerNode + dofOnNodeIndex;
}

//! get all dofs of a specific node for 3D
template <typename MeshType, typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType, BasisFunctionType,
                          Mesh::isStructuredWithDim<3, MeshType>>::
    getNodeDofs(node_no_t nodeGlobalNo,
                std::vector<dof_no_t> &dofGlobalNos) const {
  dofGlobalNos.reserve(
      dofGlobalNos.size() +
      FunctionSpaceBaseDim<3, BasisFunctionType>::nDofsPerNode());
  for (int dofIndex = 0;
       dofIndex < FunctionSpaceBaseDim<3, BasisFunctionType>::nDofsPerNode();
       dofIndex++) {
    dofGlobalNos.push_back(getNodeDofNo(nodeGlobalNo, dofIndex));
  }
}

//! get all dofs of a specific node for 3D
template <typename MeshType, typename BasisFunctionType>
void FunctionSpaceNumbers<MeshType, BasisFunctionType,
                          Mesh::isStructuredWithDim<3, MeshType>>::
    getNodeDofs(
        node_no_t nodeGlobalNo,
        std::array<dof_no_t,
                   FunctionSpaceBaseDim<3, BasisFunctionType>::nDofsPerNode()>
            &dofGlobalNos) const {
  for (int dofIndex = 0;
       dofIndex < FunctionSpaceBaseDim<3, BasisFunctionType>::nDofsPerNode();
       dofIndex++) {
    dofGlobalNos[dofIndex] = getNodeDofNo(nodeGlobalNo, dofIndex);
  }
}

//! get the dof no of the specified dof at the node for 3D
template <typename MeshType, typename BasisFunctionType>
dof_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                              Mesh::isStructuredWithDim<3, MeshType>>::
    getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const {
  return FunctionSpaceBaseDim<3, BasisFunctionType>::nDofsPerNode() *
             nodeGlobalNo +
         dofIndex;
}

// element-local nodeIndex to local nodeNo for 1D
template <typename MeshType, typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                               Mesh::isStructuredWithDim<1, MeshType>>::
    getNodeNo(element_no_t elementNoLocal, int nodeIndex) const {
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0 1
  // averageNDofsPerElement:
  // 1         2            2
  // nDofsPerBasis:
  // 2         3            4
  // nDofsPerNode:
  // 1         1            2
  // nNodesPerElement:
  // 2         3            2

  /*
  VLOG(3) << "getNodeNo<1D>(elementNoLocal=" << elementNoLocal << ", nodeIndex="
  << nodeIndex << ") = "
    <<
  FunctionSpaceFunction<MeshType,BasisFunctionType>::averageNNodesPerElement() *
  elementNoLocal + nodeIndex;
  */

  return FunctionSpaceFunction<MeshType,
                               BasisFunctionType>::averageNNodesPerElement() *
             elementNoLocal +
         nodeIndex;
}

// element-local nodeIndex to local nodeNo for 2D
template <typename MeshType, typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                               Mesh::isStructuredWithDim<2, MeshType>>::
    getNodeNo(element_no_t elementNoLocal, int nodeIndex) const {
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      2 3
  // 0 1       0 1 2      0 1
  // nNodesPerElement:
  // 4         9          4

  // since this implementation is for structured meshes only, the number of
  // elements in each coordinate direction is given

  const std::array<element_no_t, MeshType::dim()> &nElements =
      this->nElementsPerCoordinateDirectionLocal_;
  int averageNNodesPerElement1D =
      FunctionSpaceBaseDim<1, BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D =
      FunctionSpaceBaseDim<1, BasisFunctionType>::nNodesPerElement();

  // the number of non-ghost nodes in different rows
  node_no_t nodesPerRow =
      averageNNodesPerElement1D * nElements[0] +
      (this->meshPartition_->hasFullNumberOfNodes(0) ? 1 : 0);
  element_no_t elementX = element_no_t(elementNoLocal % nElements[0]);
  element_no_t elementY = element_no_t(elementNoLocal / nElements[0]);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;
  dof_no_t localY = dof_no_t(nodeIndex / nNodesPerElement1D);

  /*
  VLOG(3) << "elementY=" << elementY << ", elementX=" << elementX;
  VLOG(3) << "localX=" << localX << ", localY=" << localY;
*/

  // check for ghost nodes which have nos after the normal contiguous numbering
  // scheme
  if (!this->meshPartition_->hasFullNumberOfNodes(1) &&
      elementY == nElements[1] - 1 && localY == nNodesPerElement1D - 1) {
    // node is a ghost node on the top boundary
    // if there are ghost nodes on the right boundary
    if (!this->meshPartition_->hasFullNumberOfNodes(0)) {
      // VLOG(3) << "getNodeNo<2D>b(elementNoLocal=" << elementNoLocal << ",
      // nodeIndex=" << nodeIndex << ") = "
      //   << this->meshPartition_->nNodesLocalWithoutGhosts() +
      //   averageNNodesPerElement1D * nElements[1] + averageNNodesPerElement1D
      //   * elementX + localX;

      return this->meshPartition_->nNodesLocalWithoutGhosts() +
             averageNNodesPerElement1D * nElements[1] +
             averageNNodesPerElement1D * elementX + localX;
    } else {
      // VLOG(3) << "getNodeNo<2D>c(elementNoLocal=" << elementNoLocal << ",
      // nodeIndex=" << nodeIndex << ") = "
      //   << this->meshPartition_->nNodesLocalWithoutGhosts() +
      //   averageNNodesPerElement1D * elementX + localX;

      return this->meshPartition_->nNodesLocalWithoutGhosts() +
             averageNNodesPerElement1D * elementX + localX;
    }
  }

  if (!this->meshPartition_->hasFullNumberOfNodes(0) &&
      elementX == nElements[0] - 1 && localX == nNodesPerElement1D - 1) {
    // VLOG(3) << "getNodeNo<2D>d(elementNoLocal=" << elementNoLocal << ",
    // nodeIndex=" << nodeIndex << ") = "
    //   << this->meshPartition_->nNodesLocalWithoutGhosts() +
    //   averageNNodesPerElement1D * elementY + localY;

    // node is a ghost node on the right boundary
    return this->meshPartition_->nNodesLocalWithoutGhosts() +
           averageNNodesPerElement1D * elementY + localY;
  }

  // VLOG(3) << "getNodeNo<2D>a(elementNoLocal=" << elementNoLocal << ",
  // nodeIndex=" << nodeIndex << ") = "
  //   << nodesPerRow * (elementY * averageNNodesPerElement1D + localY) +
  //   averageNNodesPerElement1D * elementX + localX;

  // compute local node no for non-ghost node
  return nodesPerRow * (elementY * averageNNodesPerElement1D + localY) +
         averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex to local nodeNo for 3D
template <typename MeshType, typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                               Mesh::isStructuredWithDim<3, MeshType>>::
    getNodeNo(element_no_t elementNoLocal, int nodeIndex) const {
  // since this implementation is for structured meshes only, the number of
  // elements in each coordinate direction is given

  const std::array<element_no_t, MeshType::dim()> &nElements =
      this->nElementsPerCoordinateDirectionLocal_;
  const int averageNNodesPerElement1D =
      FunctionSpaceBaseDim<1, BasisFunctionType>::averageNNodesPerElement();
  const int nNodesPerElement1D =
      FunctionSpaceBaseDim<1, BasisFunctionType>::nNodesPerElement();

  node_no_t nodesPerRow0 =
      (averageNNodesPerElement1D * nElements[0] +
       (this->meshPartition_->hasFullNumberOfNodes(0) ? 1 : 0));
  node_no_t nodesPerPlane =
      (averageNNodesPerElement1D * nElements[1] +
       (this->meshPartition_->hasFullNumberOfNodes(1) ? 1 : 0)) *
      nodesPerRow0;

  element_no_t elementZ =
      element_no_t(elementNoLocal / (nElements[0] * nElements[1]));
  element_no_t elementY = element_no_t(
      (elementNoLocal % (nElements[0] * nElements[1])) / nElements[0]);
  element_no_t elementX = elementNoLocal % nElements[0];
  dof_no_t localZ = dof_no_t(nodeIndex / MathUtility::sqr(nNodesPerElement1D));
  dof_no_t localY = dof_no_t(
      (nodeIndex % MathUtility::sqr(nNodesPerElement1D)) / nNodesPerElement1D);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;

  // check for ghost nodes which have nos after the normal contiguous numbering
  // scheme handle ghosts on z+ boundary
  if (!this->meshPartition_->hasFullNumberOfNodes(2) &&
      elementZ == nElements[2] - 1 && localZ == nNodesPerElement1D - 1) {
    // node is a ghost node on the z+ boundary, i.e. lies in a x-y plane on top
    // of the non-ghost dofs get number of first dof on z+ boundary
    node_no_t nodeNo = this->meshPartition_->nNodesLocalWithoutGhosts();
    if (this->meshPartition_->hasFullNumberOfNodes(
            0)) // if there are no ghosts on x+
    {
      if (!this->meshPartition_->hasFullNumberOfNodes(
              1)) // if there are ghosts on y+
      {
        nodeNo +=
            nodesPerRow0 * (elementZ * averageNNodesPerElement1D + localZ);
      }
    } else // if there are ghosts on x+
    {
      if (this->meshPartition_->hasFullNumberOfNodes(
              1)) // if there are no ghosts on y+
      {
        node_no_t nodesPerRow1 = averageNNodesPerElement1D * nElements[1] + 1;
        nodeNo +=
            nodesPerRow1 * (elementZ * averageNNodesPerElement1D + localZ);
      } else // if there are ghosts on y+
      {
        node_no_t nodesPerRow1 = averageNNodesPerElement1D * nElements[1];
        nodeNo += (nodesPerRow0 + nodesPerRow1 + 1) *
                  (elementZ * averageNNodesPerElement1D + localZ);
      }
    }

    // the z+ plane always has the full number of nodes in x and y direction (in
    // lower z levels these are either ghost or non-ghost nodes)
    nodeNo += (averageNNodesPerElement1D * nElements[0] + 1) *
              (elementY * averageNNodesPerElement1D + localY);
    nodeNo += averageNNodesPerElement1D * elementX + localX;

    // VLOG(3) << "getNodeNo<3D>(elementNoLocal=" << elementNoLocal << ",
    // nodeIndex=" << nodeIndex << ") = " << nodeNo;
    return nodeNo;
  }

  // handle ghosts on y+ boundary
  if (!this->meshPartition_->hasFullNumberOfNodes(1) &&
      elementY == nElements[1] - 1 && localY == nNodesPerElement1D - 1) {
    // node is a ghost node on the y+ boundary, i.e. lies in a x-z plane
    node_no_t nodeNo = this->meshPartition_->nNodesLocalWithoutGhosts();
    if (this->meshPartition_->hasFullNumberOfNodes(0)) {
      nodeNo += nodesPerRow0 * (elementZ * averageNNodesPerElement1D + localZ);
    } else {
      node_no_t nodesPerRow1 = averageNNodesPerElement1D * nElements[1];
      nodeNo += (nodesPerRow0 + nodesPerRow1 + 1) *
                (elementZ * averageNNodesPerElement1D + localZ);
      nodeNo += nodesPerRow1;
    }

    nodeNo += averageNNodesPerElement1D * elementX + localX;

    // VLOG(3) << "getNodeNo<3D>(elementNoLocal=" << elementNoLocal << ",
    // nodeIndex=" << nodeIndex << ") = " << nodeNo;
    return nodeNo;
  }

  // handle ghosts on x+ boundary
  if (!this->meshPartition_->hasFullNumberOfNodes(0) &&
      elementX == nElements[0] - 1 && localX == nNodesPerElement1D - 1) {
    // node is a ghost node on the x+ boundary, i.e. lies in a y-z plane
    node_no_t nodeNo = this->meshPartition_->nNodesLocalWithoutGhosts();
    if (this->meshPartition_->hasFullNumberOfNodes(1)) {
      node_no_t nodesPerRow1 = averageNNodesPerElement1D * nElements[1] + 1;
      nodeNo += nodesPerRow1 * (elementZ * averageNNodesPerElement1D + localZ);
    } else {
      node_no_t nodesPerRow1 = averageNNodesPerElement1D * nElements[1];
      nodeNo += (nodesPerRow0 + nodesPerRow1 + 1) *
                (elementZ * averageNNodesPerElement1D + localZ);
    }
    nodeNo += averageNNodesPerElement1D * elementY + localY;

    // VLOG(3) << "getNodeNo<3D>(elementNoLocal=" << elementNoLocal << ",
    // nodeIndex=" << nodeIndex << ") = " << nodeNo;
    return nodeNo;
  }

  // VLOG(3) << "getNodeNo<3D>(elementNoLocal=" << elementNoLocal << ",
  // nodeIndex=" << nodeIndex << ") = "
  //   << nodesPerPlane * (elementZ * averageNNodesPerElement1D + localZ) +
  //   nodesPerRow0 * (elementY * averageNNodesPerElement1D + localY) +
  //   averageNNodesPerElement1D * elementX + localX;

  // compute local node no for non-ghost node
  return nodesPerPlane * (elementZ * averageNNodesPerElement1D + localZ) +
         nodesPerRow0 * (elementY * averageNNodesPerElement1D + localY) +
         averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex of global element to global nodeNo for 1D
template <typename MeshType, typename BasisFunctionType>
global_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                                 Mesh::isStructuredWithDim<1, MeshType>>::
    getNodeNoGlobalNatural(global_no_t elementNoLocalGlobalNatural,
                           int nodeIndex) const {
  return FunctionSpaceFunction<MeshType,
                               BasisFunctionType>::averageNNodesPerElement() *
             elementNoLocalGlobalNatural +
         nodeIndex;
}

// element-local nodeIndex of global element to global nodeNo for 2D
template <typename MeshType, typename BasisFunctionType>
global_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                                 Mesh::isStructuredWithDim<2, MeshType>>::
    getNodeNoGlobalNatural(global_no_t elementNoLocalGlobalNatural,
                           int nodeIndex) const {
  const std::array<global_no_t, MeshType::dim()> &nElements =
      this->nElementsPerCoordinateDirectionGlobal_;
  int averageNNodesPerElement1D =
      FunctionSpaceBaseDim<1, BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D =
      FunctionSpaceBaseDim<1, BasisFunctionType>::nNodesPerElement();

  // the number of non-ghost nodes in different rows
  global_no_t nodesPerRow = averageNNodesPerElement1D * nElements[0] + 1;
  element_no_t elementX =
      element_no_t(elementNoLocalGlobalNatural % nElements[0]);
  element_no_t elementY =
      element_no_t(elementNoLocalGlobalNatural / nElements[0]);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;
  dof_no_t localY = dof_no_t(nodeIndex / nNodesPerElement1D);

  return nodesPerRow * (elementY * averageNNodesPerElement1D + localY) +
         averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex of global element to global nodeNo for 3D
template <typename MeshType, typename BasisFunctionType>
global_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                                 Mesh::isStructuredWithDim<3, MeshType>>::
    getNodeNoGlobalNatural(global_no_t elementNoLocalGlobalNatural,
                           int nodeIndex) const {
  const std::array<global_no_t, MeshType::dim()> &nElements =
      this->nElementsPerCoordinateDirectionGlobal_;
  int averageNNodesPerElement1D =
      FunctionSpaceBaseDim<1, BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D =
      FunctionSpaceBaseDim<1, BasisFunctionType>::nNodesPerElement();
  global_no_t nodesPerRow0 = averageNNodesPerElement1D * nElements[0] + 1;
  global_no_t nodesPerPlane =
      (averageNNodesPerElement1D * nElements[1] + 1) * nodesPerRow0;

  element_no_t elementZ =
      element_no_t(elementNoLocalGlobalNatural / (nElements[0] * nElements[1]));
  element_no_t elementY = element_no_t(
      (elementNoLocalGlobalNatural % (nElements[0] * nElements[1])) /
      nElements[0]);
  element_no_t elementX = elementNoLocalGlobalNatural % nElements[0];
  node_no_t localZ =
      node_no_t(nodeIndex / MathUtility::sqr(nNodesPerElement1D));
  node_no_t localY = node_no_t(
      (nodeIndex % MathUtility::sqr(nNodesPerElement1D)) / nNodesPerElement1D);
  node_no_t localX = nodeIndex % nNodesPerElement1D;

  // compute local node no for non-ghost node
  return nodesPerPlane * (elementZ * averageNNodesPerElement1D + localZ) +
         nodesPerRow0 * (elementY * averageNNodesPerElement1D + localY) +
         averageNNodesPerElement1D * elementX + localX;
}

// local node no of neighbour node, may be a ghost node, for 1D
template <typename MeshType, typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                               Mesh::isStructuredWithDim<1, MeshType>>::
    getNeighbourNodeNoLocal(node_no_t nodeNoLocal,
                            Mesh::face_t direction) const {
  assert(nodeNoLocal < this->nNodesLocalWithoutGhosts());

  if (direction == Mesh::face_t::face0Minus) // left node
  {
    // if there is no neighbouring node
    if (nodeNoLocal == 0) {
      return -1;
    }

    // get previous node no
    return nodeNoLocal - 1;
  } else if (direction == Mesh::face_t::face0Plus) // right node
  {
    // if there is no neighbouring node
    if (nodeNoLocal == this->meshPartition_->nNodesLocalWithGhosts(0) - 1) {
      return -1;
    }

    // get next node no, this is correct for ghost as well as non-ghost nodes
    return nodeNoLocal + 1;
  } else {
    assert(false);
  }

  return 0; // should not happen, but cray compiler does not recognize it
}

// local node no of neighbour node, may be a ghost node, for 2D
template <typename MeshType, typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                               Mesh::isStructuredWithDim<2, MeshType>>::
    getNeighbourNodeNoLocal(node_no_t nodeNoLocal,
                            Mesh::face_t direction) const {
  assert(nodeNoLocal < this->nNodesLocalWithoutGhosts());

  dof_no_t localX =
      nodeNoLocal % this->meshPartition_->nNodesLocalWithoutGhosts(0);
  dof_no_t localY =
      dof_no_t(nodeNoLocal / this->meshPartition_->nNodesLocalWithoutGhosts(0));

  if (direction == Mesh::face_t::face0Minus) // left node
  {
    // if there is no neighbouring node
    if (localX == 0) {
      return -1;
    }
    // get previous node no
    return nodeNoLocal - 1;
  } else if (direction == Mesh::face_t::face0Plus) // right node
  {
    // if the node is on the right row, there is no neighbouring node
    if (localX == this->meshPartition_->nNodesLocalWithGhosts(0) - 1) {
      return -1;
    }

    // check if right node is a ghost node
    if (localX == this->meshPartition_->nNodesLocalWithoutGhosts(0) - 1) {
      // there is a row of ghost nodes on the right
      return this->meshPartition_->nNodesLocalWithoutGhosts() + localY;
    }

    // return next node no
    return nodeNoLocal + 1;
  } else if (direction == Mesh::face_t::face1Minus) // bottom node
  {
    // if there is no neighbouring node
    if (localY == 0) {
      return -1;
    }

    // get node below
    return nodeNoLocal - this->meshPartition_->nNodesLocalWithoutGhosts(0);
  } else if (direction == Mesh::face_t::face1Plus) // top node
  {
    // if the node is at the top row, there is no neighbouring node
    if (localY == this->meshPartition_->nNodesLocalWithGhosts(1) - 1) {
      return -1;
    }

    if (localY == this->meshPartition_->nNodesLocalWithoutGhosts(1) - 1) {
      // there is one row of ghost nodes above the current
      if (this->meshPartition_->hasFullNumberOfNodes(0)) {
        // there are only top ghost nodes
        return this->meshPartition_->nNodesLocalWithoutGhosts() + localX;
      } else {
        // there are right and top ghost nodes
        return this->meshPartition_->nNodesLocalWithoutGhosts() +
               this->meshPartition_->nNodesLocalWithGhosts(1) - 1 + localX;
      }
    }

    // get node above
    return nodeNoLocal + this->meshPartition_->nNodesLocalWithoutGhosts(0);
  } else {
    assert(false);
  }
  return 0; // should not happen, but cray compiler does not recognize it
}

// local node no of neighbour node, may be a ghost node, for 3D
template <typename MeshType, typename BasisFunctionType>
node_no_t FunctionSpaceNumbers<MeshType, BasisFunctionType,
                               Mesh::isStructuredWithDim<3, MeshType>>::
    getNeighbourNodeNoLocal(node_no_t nodeNoLocal,
                            Mesh::face_t direction) const {
  assert(nodeNoLocal < this->nNodesLocalWithoutGhosts());

  node_no_t localZ = node_no_t(
      nodeNoLocal / (this->meshPartition_->nNodesLocalWithoutGhosts(0) *
                     this->meshPartition_->nNodesLocalWithoutGhosts(1)));
  node_no_t localY = node_no_t(
      (nodeNoLocal % (this->meshPartition_->nNodesLocalWithoutGhosts(0) *
                      this->meshPartition_->nNodesLocalWithoutGhosts(1))) /
      this->meshPartition_->nNodesLocalWithoutGhosts(0));
  node_no_t localX =
      nodeNoLocal % this->meshPartition_->nNodesLocalWithoutGhosts(0);

  if (direction == Mesh::face_t::face0Minus) // left node
  {
    // if there is no neighbouring node
    if (localX == 0) {
      return -1;
    }

    // get previous node no
    return nodeNoLocal - 1;
  } else if (direction == Mesh::face_t::face0Plus) // right node
  {
    // if the node is on the right row, there is no neighbouring node
    if (localX == this->meshPartition_->nNodesLocalWithGhosts(0) - 1) {
      return -1;
    }

    // if right node is a ghost node, this implies
    // !this->meshPartition_->hasFullNumberOfNodes(0)
    if (localX == this->meshPartition_->nNodesLocalWithoutGhosts(0) - 1) {
      // there is a y-z plane of ghost nodes on the right
      node_no_t neighbourNodeNo =
          this->meshPartition_->nNodesLocalWithoutGhosts();

      if (this->meshPartition_->hasFullNumberOfNodes(1)) {
        neighbourNodeNo +=
            (this->meshPartition_->nNodesLocalWithGhosts(1)) * localZ;
      } else {
        neighbourNodeNo +=
            (this->meshPartition_->nNodesLocalWithGhosts(0) +
             this->meshPartition_->nNodesLocalWithGhosts(1) - 1) *
            localZ;
      }
      return neighbourNodeNo + localY;
    }

    // return next node no
    return nodeNoLocal + 1;
  } else if (direction == Mesh::face_t::face1Minus) // y- node
  {
    // if there is no neighbouring node
    if (localY == 0) {
      return -1;
    }

    // get node below
    return nodeNoLocal - this->meshPartition_->nNodesLocalWithoutGhosts(0);
  } else if (direction == Mesh::face_t::face1Plus) // y+ node
  {
    // LOG(DEBUG) << "direction Y+, localX=" << localX << ", localY=" << localY
    // << ", localZ=" << localZ << ", pf" <<
    // this->meshPartition_->nNodesLocalWithoutGhosts(1)-1;

    // if the node is at the top row, there is no neighbouring node
    if (localY == this->meshPartition_->nNodesLocalWithGhosts(1) - 1) {
      return -1;
    }

    // if y+ node is a ghost node, this implies
    // !this->meshPartition_->hasFullNumberOfNodes(1)
    if (localY == this->meshPartition_->nNodesLocalWithoutGhosts(1) - 1) {
      // there is one row of ghost nodes behind the current
      node_no_t neighbourNodeNo =
          this->meshPartition_->nNodesLocalWithoutGhosts();

      if (this->meshPartition_->hasFullNumberOfNodes(0)) {
        // there are only ghosts at y+
        neighbourNodeNo +=
            (this->meshPartition_->nNodesLocalWithoutGhosts(0)) * localZ;
      } else {
        // there are ghosts at y+ and x+
        neighbourNodeNo +=
            (this->meshPartition_->nNodesLocalWithGhosts(0) +
             this->meshPartition_->nNodesLocalWithGhosts(1) - 1) *
            localZ;
        neighbourNodeNo += this->meshPartition_->nNodesLocalWithoutGhosts(1);
      }
      return neighbourNodeNo + localX;
    }

    // get node behind
    return nodeNoLocal + this->meshPartition_->nNodesLocalWithoutGhosts(0);
  } else if (direction == Mesh::face_t::face2Minus) // z- node
  {
    // if there is no node below
    if (localZ == 0) {
      return -1;
    }

    // get node below
    return nodeNoLocal - this->meshPartition_->nNodesLocalWithoutGhosts(0) *
                             this->meshPartition_->nNodesLocalWithoutGhosts(1);
  } else if (direction == Mesh::face_t::face2Plus) // z+ node
  {
    // if the node is at the top row, there is no neighbouring node
    if (localZ == this->meshPartition_->nNodesLocalWithGhosts(2) - 1) {
      return -1;
    }

    // if z+ node is a ghost node, this implies
    // !this->meshPartition_->hasFullNumberOfNodes(2)
    if (localZ == this->meshPartition_->nNodesLocalWithoutGhosts(2) - 1) {
      // there is one row of ghost nodes above the current
      node_no_t neighbourNodeNo =
          this->meshPartition_->nNodesLocalWithoutGhosts();

      // LOG(DEBUG) << "direction Z+, localX=" << localX << ", localY=" <<
      // localY << ", localZ=" << localZ << ", starting at " << neighbourNodeNo
      //   << ", nNodesLocalWithoutGhosts: " <<
      //   this->meshPartition_->nNodesLocalWithoutGhosts(0) << "," <<
      //   this->meshPartition_->nNodesLocalWithoutGhosts(1) << "," <<
      //   this->meshPartition_->nNodesLocalWithoutGhosts(2)
      //   << ", hasFullNumberOfNodes: " <<
      //   this->meshPartition_->hasFullNumberOfNodes(0) << "," <<
      //   this->meshPartition_->hasFullNumberOfNodes(1) << "," <<
      //   this->meshPartition_->hasFullNumberOfNodes(2);

      if (this->meshPartition_->hasFullNumberOfNodes(0)) {
        if (!this->meshPartition_->hasFullNumberOfNodes(1)) {
          // there are ghosts at y+ and z+
          neighbourNodeNo += this->meshPartition_->nNodesLocalWithoutGhosts(0) *
                             this->meshPartition_->nNodesLocalWithoutGhosts(2);
        }
      } else {
        if (this->meshPartition_->hasFullNumberOfNodes(1)) {
          // there are ghosts at x+ and z+
          neighbourNodeNo += this->meshPartition_->nNodesLocalWithoutGhosts(1) *
                             this->meshPartition_->nNodesLocalWithoutGhosts(2);
        } else {
          // there are ghosts at x+, y+ and z+
          neighbourNodeNo +=
              (this->meshPartition_->nNodesLocalWithGhosts(0) +
               this->meshPartition_->nNodesLocalWithGhosts(1) - 1) *
              this->meshPartition_->nNodesLocalWithoutGhosts(2);
        }
      }

      if (this->meshPartition_->hasFullNumberOfNodes(0)) {
        // there are no ghosts at x+
        return neighbourNodeNo +
               this->meshPartition_->nNodesLocalWithoutGhosts(0) * localY +
               localX;
      } else {
        // there are ghosts at x+
        return neighbourNodeNo +
               this->meshPartition_->nNodesLocalWithGhosts(0) * localY + localX;
      }
    }

    // get node above
    return nodeNoLocal + this->meshPartition_->nNodesLocalWithoutGhosts(0) *
                             this->meshPartition_->nNodesLocalWithoutGhosts(1);
  } else {
    assert(false);
  }
  return 0; // should not happen, but cray compiler does not recognize it
}

//! get the node no in the global natural ordering
template <typename MeshType, typename BasisFunctionType>
global_no_t FunctionSpaceNumbersCommon<MeshType, BasisFunctionType,
                                       Mesh::isStructured<MeshType>>::
    getNodeNoGlobalNaturalFromElementNoLocal(element_no_t elementNoLocal,
                                             int nodeIndex) const {
  global_no_t elementNoGlobalNatural =
      this->meshPartition_->getElementNoGlobalNatural(elementNoLocal);
  return this->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);
}

} // namespace FunctionSpace
