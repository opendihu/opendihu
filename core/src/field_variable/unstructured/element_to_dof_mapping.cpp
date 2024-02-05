#include "field_variable/unstructured/element_to_dof_mapping.h"

#include <cassert>
#include <set>

#include "utility/string_utility.h"
#include "utility/math_utility.h"
#include "control/types.h"
#include "easylogging++.h"

namespace FieldVariable {

void ElementToDofMapping::setNumberElements(element_no_t nElements) {
  elementDofs_.resize(nElements);
}

std::shared_ptr<NodeToDofMapping> ElementToDofMapping::setup(
    std::shared_ptr<ExfileRepresentation> exfileRepresentation,
    std::shared_ptr<ElementToNodeMapping> elementToNodeMapping,
    const int nDofsPerNode) {
  node_no_t dofGlobalNo = 0;

  // for setup to work we need the number of elements already set (by a previous
  // call to setNumberElements)
  assert(elementDofs_.size() != 0);
  element_no_t nElements = elementDofs_.size();

  // create node to dof mapping
  std::shared_ptr<NodeToDofMapping> nodeToDofMapping =
      std::make_shared<NodeToDofMapping>();

  VLOG(1) << "setup element to dof mapping and node to dof mapping from "
             "element to node mapping";

  // loop over elements
  for (element_no_t elementGlobalNo = 0; elementGlobalNo < nElements;
       elementGlobalNo++) {
    // get available information about element
    ElementToNodeMapping::Element &element =
        elementToNodeMapping->getElement(elementGlobalNo);
    std::shared_ptr<ExfileElementRepresentation> exfileElement =
        exfileRepresentation->getExfileElementRepresentation(elementGlobalNo);

    // element contains the following fields:
    //   std::vector<int> nodeGlobalNo;
    //   std::vector<double> scaleFactors;

    // resize dofs vector for element
    unsigned int nNodesInElement = element.nodeGlobalNo.size();
    elementDofs_[elementGlobalNo].resize(nNodesInElement * nDofsPerNode);

    VLOG(1) << "element " << elementGlobalNo << ", nNodes: " << nNodesInElement;

    int elementDofIndex = 0;

    // loop over nodes of element
    for (unsigned int nodeIndex = 0; nodeIndex < nNodesInElement; nodeIndex++) {
      node_no_t nodeGlobalNo = element.nodeGlobalNo[nodeIndex];

      VLOG(1) << "   node global " << nodeGlobalNo;

      ExfileElementRepresentation::Node exfileNode =
          exfileElement->getNode(nodeIndex);
      /*
        struct Node
        {
          std::vector<int> valueIndices;          //< the indices of the dof
        values of this node in the exnode file node values (sub-)block (for the
        particular field variable/component) of the node. If there are not
        multiple versions, this is simply 0,1,...,ndofs-1. If there are e.g. 2
        versions and 8 dofs per node, this can be 0,1,...,7 if the elements uses
        the 1st version, or 8,...,15 if the element uses the second version.
        Note, that the real index of the dofs inside the values block may be
        different when this is not the first component of the block.
          std::vector<int> scaleFactorIndices;   //< the indices of all scale
        factor entries for this node in the exelem element scale factors block.
        Thus this is kind of a node to element block mapping.
        };
       */

      // get the version no of this element at the node
      unsigned int versionNo = int(exfileNode.valueIndices[0] / nDofsPerNode);

      // if node was not visited yet
      if (!nodeToDofMapping->containsNode(nodeGlobalNo)) {
        // set nodeDofInformation at nodeGlobalNo
        NodeToDofMapping::NodeDofInformation &nodeDofInformation =
            nodeToDofMapping->getNodeDofInformation(nodeGlobalNo);

        /*
          struct NodeDofInformation
          {
            std::vector<int> dofs;               //< the dofs of the node
            std::vector<double> scaleFactors;    //< the scale factors for
          Hermite basis functions, not used yet, but set maybe for later use

            struct elementNoLocalde
            {
              element_no_t elementGlobalNo;      //< element global no
              unsigned int nodeIdx;              //< node index local to the
          element
            };
            std::vector<std::vector<elementNoLocalde>> elementsOfVersion;  //<
          for each version the global element no. and node index of the element
          that is adjacent to this node and use the specified version. The
          number of versions is thus the size of the outer vector
          };
        */

        nodeDofInformation.dofs.resize(exfileNode.valueIndices.back() + 1);
        nodeDofInformation.elementsOfVersion.resize(
            versionNo + 1); // create version with given no.
        // add adjacent element to version
        // NodeToDofMapping::NodeDofInformation::elementNoLocalde
        // elementNoLocalde;
        nodeDofInformation.elementsOfVersion[versionNo].push_back(
            {elementGlobalNo, nodeIndex});

        // copy exfileNode valueIndices to nodeDofInformation
        // std::copy(exfileNode.valueIndices.begin(),
        // exfileNode.valueIndices.end(),
        // nodeDofInformation.exfileValueIndices.begin());

        // assign global dof nos
        for (int dofIndex = 0; dofIndex < nDofsPerNode; dofIndex++) {
          elementDofs_[elementGlobalNo][elementDofIndex++] = dofGlobalNo;
          nodeDofInformation.dofs[exfileNode.valueIndices[dofIndex]] =
              dofGlobalNo;
          VLOG(1) << "     dofIndex " << dofIndex << " valueIndex "
                  << exfileNode.valueIndices[dofIndex] << " dof global "
                  << dofGlobalNo;
          dofGlobalNo++;
        }

        VLOG(1) << "      not yet visited, min."
                << nodeDofInformation.dofs.size()
                << " value indices: " << exfileNode.valueIndices
                << ", dofs at node: " << nodeDofInformation.dofs
                << ", nDofsPerNode=" << nDofsPerNode
                << ", versionNo: " << versionNo;
      } else {
        // if node was already visited, i.e. it is adjacent to an element with
        // an earlier element no.
        NodeToDofMapping::NodeDofInformation &nodeDofInformation =
            nodeToDofMapping->getNodeDofInformation(nodeGlobalNo);

        // check if dofs of version are already assigned
        bool dofsAlreadyAssigned = false;
        if (nodeDofInformation.elementsOfVersion.size() > versionNo)
          dofsAlreadyAssigned =
              !nodeDofInformation.elementsOfVersion[versionNo].empty();

        // node was visited earlier and has assigned dofs for this version
        if (dofsAlreadyAssigned) {
          // add element to element list of version
          nodeDofInformation.elementsOfVersion[versionNo].push_back(
              {elementGlobalNo, nodeIndex});

          // get dofs
          for (int dofIndex = 0; dofIndex < nDofsPerNode; dofIndex++) {
            dof_no_t dofNo =
                nodeDofInformation.dofs[exfileNode.valueIndices[dofIndex]];
            elementDofs_[elementGlobalNo][elementDofIndex++] = dofNo;
          }

          VLOG(1) << "      already visited on this version, elements of this "
                     "version: "
                  << nodeDofInformation.elementsOfVersion[versionNo];
        } else {
          // node was visited earlier but has no dofs for this version
          if ((int)nodeDofInformation.dofs.size() <
              exfileNode.valueIndices.back() + 1)
            nodeDofInformation.dofs.resize(exfileNode.valueIndices.back() + 1);

          if (nodeDofInformation.elementsOfVersion.size() < versionNo + 1)
            nodeDofInformation.elementsOfVersion.resize(
                versionNo + 1); // create version with given no.

          // add adjacent element to version
          nodeDofInformation.elementsOfVersion[versionNo].push_back(
              {elementGlobalNo, nodeIndex});

          // assign global dof nos
          for (int dofIndex = 0; dofIndex < nDofsPerNode; dofIndex++) {
            elementDofs_[elementGlobalNo][elementDofIndex++] = dofGlobalNo;
            nodeDofInformation.dofs[exfileNode.valueIndices[dofIndex]] =
                dofGlobalNo;
            dofGlobalNo++;
          }

          VLOG(1) << "      already visited, not on this version, elements of "
                     "this version: "
                  << nodeDofInformation.elementsOfVersion[versionNo];
        }
      }

      // debugging output
      NodeToDofMapping::NodeDofInformation &nodeDofInformation =
          nodeToDofMapping->getNodeDofInformation(nodeGlobalNo);

      VLOG(1) << "     dofs: " << nodeDofInformation.dofs;
      for (unsigned int versionIdx = 0;
           versionIdx < nodeDofInformation.elementsOfVersion.size();
           versionIdx++) {
        std::stringstream s;
        s << "(";
        for (int i = 0;
             i < nodeDofInformation.elementsOfVersion[versionIdx].size(); i++) {
          s << nodeDofInformation.elementsOfVersion[versionIdx][i]
                   .elementGlobalNo
            << "."
            << nodeDofInformation.elementsOfVersion[versionIdx][i].nodeIdx
            << " ";
        }
        s << ")";
        VLOG(1) << "     version " << versionIdx << " elements: " << s.str();
      }

#ifndef NDEBUG
      // check if the element has the same nodes multiple times
      std::set<node_no_t> nodesOfElement;
      for (unsigned int nodeIndex = 0; nodeIndex < nNodesInElement;
           nodeIndex++) {
        node_no_t nodeGlobalNo = element.nodeGlobalNo[nodeIndex];

        if (nodesOfElement.find(nodeGlobalNo) != nodesOfElement.end()) {
          LOG(ERROR) << "Element " << elementGlobalNo << " contains node "
                     << nodeGlobalNo << " multiple times!";
        }
        nodesOfElement.insert(nodeGlobalNo);
      }
#endif

    } // nodeIdx
  }   // elementGlobalNo

  nDofs_ = dofGlobalNo;

  return nodeToDofMapping;
}

dof_no_t ElementToDofMapping::nDofsLocal() const { return nDofs_; }

element_no_t ElementToDofMapping::nElementsLocal() const {
  return elementDofs_.size();
}

const std::vector<dof_no_t> &
ElementToDofMapping::getElementDofs(element_no_t elementGlobalNo) const {
  return elementDofs_[elementGlobalNo];
}

bool ElementToDofMapping::operator==(const ElementToDofMapping &rhs) {
  if (nDofs_ != rhs.nDofs_)
    return false;

  if (elementDofs_.size() != rhs.elementDofs_.size())
    return false;

  for (unsigned int i = 0; i < elementDofs_.size(); i++) {
    if (elementDofs_[i].size() != rhs.elementDofs_[i].size())
      return false;
    for (unsigned int j = 0; j < elementDofs_[i].size(); j++)
      if (elementDofs_[i][j] != rhs.elementDofs_[i][j])
        return false;
  }
  return true;
}

} // namespace FieldVariable
