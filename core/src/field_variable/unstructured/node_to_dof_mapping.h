#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <iostream>
#include <memory>
#include <map>

#include "control/types.h"
#include "field_variable/unstructured/exfile_representation.h"
#include "field_variable/unstructured/element_to_node_mapping.h"

namespace FieldVariable
{

/** For every node the dofs, the scale factors and the adjacent elements and the versions and the adjacent element-version mapping.
 *  This data structure is generated by the setup method of ElementToDofMapping.
 */
class NodeToDofMapping
{
public:
 
  /** all the dofs of a node, also the indices at which position of the node in the exfiles the value occurs
   */
  struct NodeDofInformation
  {
    std::vector<dof_no_t> dofs;               ///< the dofs of the node
    std::vector<double> scaleFactors;    ///< the scale factors for Hermite basis functions, not used yet, but set maybe for later use
    
    struct ElementLocalNode 
    {
      element_no_t elementGlobalNo;    ///< element global no
      unsigned int nodeIdx;       ///< node index local to the element
    };
    std::vector<std::vector<ElementLocalNode>> elementsOfVersion;  ///< for each version the global element no. and node index of the element that is adjacent to this node and use the specified version. The number of versions is thus the size of the outer vector
  };
  
  //! get the dof information of a node
  NodeDofInformation &getNodeDofInformation(node_no_t nodeGlobalNo);
  
  //! get the dofs of a node
  std::vector<dof_no_t> &getNodeDofs(node_no_t nodeGlobalNo);
  
  //! get the dofs of a node
  std::vector<double> &getNodeScaleFactors(node_no_t nodeGlobalNo);
  
  //! return the number of nodes
  node_no_t nNodes() const;
  
  //! check if a node is already contained in the internal node dof information map
  bool containsNode(node_no_t nodeGlobalNo) const;
  
  //! return the number of versions of a particular node. A version is an OpenCMISS iron construct that allows multiple dofs of the same type on a single node, useful for modeling discontinuities
  int nVersions(node_no_t nodeGlobalNo);
  
  //! output string representation
  void output(std::ostream &stream) const;
  
private:
  std::map<node_no_t, NodeDofInformation> nodeDofInformation_;   ///< for global node number the associated dofs, node-to-dof mapping
};

std::ostream &operator<<(std::ostream &stream, const NodeToDofMapping &rhs);
std::ostream &operator<<(std::ostream &stream, const NodeToDofMapping::NodeDofInformation::ElementLocalNode &rhs);

};  // namespace
