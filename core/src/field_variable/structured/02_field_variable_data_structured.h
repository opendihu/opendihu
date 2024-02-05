#pragma once

#include <Python.h> // has to be the first included header
#include <iostream>
#include <array>
#include <map>
#include <petscvec.h>
#include "easylogging++.h"

#include "field_variable/01_field_variable_components.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"

namespace FieldVariable {

// forward declaration
template <typename FunctionSpaceType, int nComponents> class FieldVariable;

/** Field variable for a structured mesh, i.e. dof and node information are
 * purely implicit. This is used for RegularFixed and StructuredDeformable
 * meshes.
 */
template <typename FunctionSpaceType, int nComponents_>
class FieldVariableDataStructured
    : public FieldVariableComponents<FunctionSpaceType, nComponents_> {
public:
  //! normal constructor without arguments
  FieldVariableDataStructured();

  //! contructor as data copy (reuseData=false) or reusing the Petsc Vec's
  //! (reuseData=true), with a different name (component names are the same).
  //! Note, it is not possible to make rhs const, because VecCopy needs
  //! globalValues() and this may change rhs
  FieldVariableDataStructured(
      FieldVariable<FunctionSpaceType, nComponents_> &rhs, std::string name,
      bool reuseData = false);

  //! contructor as data copy (reuseData=false) or reusing the Petsc Vec's
  //! (reuseData=true), with a different name and different components,
  //! @param rhsComponentNoBegin the first component of the rhs to reuse for
  //! this field variable
  template <int nComponents2>
  FieldVariableDataStructured(
      FieldVariable<FunctionSpaceType, nComponents2> &rhs, std::string name,
      std::vector<std::string> componentNames, bool reuseData = false,
      int rhsComponentNoBegin = 0);

  //! constructor with functionSpace, name and components and if it is a
  //! geometry field. This constructs a complete field variable
  FieldVariableDataStructured(std::shared_ptr<FunctionSpaceType> functionSpace,
                              std::string name,
                              std::vector<std::string> componentNames,
                              bool isGeometryField = false);

  //! destructor
  virtual ~FieldVariableDataStructured();

  //! write a exelem file header to a stream, for a particular element,
  //! fieldVariableNo is the field index x) in the exelem file header. For
  //! parallel program execution this writes headers for the local exelem files
  //! on every rank.
  void outputHeaderExelem(std::ostream &file,
                          element_no_t currentElementGlobalNo,
                          int fieldVariableNo = -1);

  //! write a exelem file header to a stream, for a particular node, for
  //! parallel program execution this writes headers for the local exnodes files
  //! on every rank.
  void outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo,
                          int &valueIndex, int fieldVariableNo = -1);

  //! tell if 2 elements have the same exfile representation, i.e. same number
  //! of versions
  bool haveSameExfileRepresentation(element_no_t element1,
                                    element_no_t element2);

  //! get the internal PETSc vector values, the local vector for the specified
  //! component
  Vec &valuesLocal(int componentNo = 0);

  //! get the internal PETSc vector values, the global vector for the specified
  //! component
  Vec &valuesGlobal(int componentNo);

  //! if the vector has multiple components, return a nested Vec of the global
  //! vector, else return the global vector
  Vec &valuesGlobal();

  //! fill a contiguous vector with all components after each other, "struct of
  //! array"-type data layout. after manipulation of the vector has finished one
  //! has to call restoreValuesContiguous
  Vec &getValuesContiguous();

  //! copy the values back from a contiguous representation where all components
  //! are in one vector to the standard internal format of PartitionedPetscVec
  //! where there is one local vector with ghosts for each component. this has
  //! to be called
  void restoreValuesContiguous();

  //! output string representation to stream for debugging
  void output(std::ostream &stream) const;

  //! not implemented interface methods

  //! parse current component's exfile representation from file contents
  virtual void parseHeaderFromExelemFile(std::string content) {}

  //! parse single element from exelem file
  virtual void parseElementFromExelemFile(std::string content) {}

  //! read in values from exnode file
  virtual void parseFromExnodeFile(std::string content) {}

  //! resize internal representation variable to number of elements
  virtual void setNumberElements(element_no_t nElements) {}

  //! reduce memory consumption by removing duplicates in ExfileRepresentations
  virtual void
  unifyMappings(std::shared_ptr<ElementToNodeMapping> elementToNodeMapping,
                const int nDofsPerNode) {}

  //! eliminate duplicate elementToDof and exfileRepresentation objects in
  //! components of two field variables (this and one other)
  virtual void unifyMappings(
      std::shared_ptr<FieldVariableBaseFunctionSpace<FunctionSpaceType>>
          fieldVariable2) {}

  //! initialize PETSc vector with size of total number of dofs for all
  //! components of this field variable
  virtual void initializeValuesVector() {}

  //! return the component by index
  virtual std::shared_ptr<Component<FunctionSpaceType, nComponents_>>
  component(int componentNo) {
    return nullptr;
  } // return empty Component

  //! get the element to dof mapping object
  virtual std::shared_ptr<ElementToDofMapping> elementToDofMapping() const {
    return nullptr;
  }

  //! get the node to dof mapping object
  virtual std::shared_ptr<NodeToDofMapping> nodeToDofMapping() const {
    return nullptr;
  }

  //! get the number of scale factors
  virtual int getNumberScaleFactors(element_no_t elementGlobalNo) const {
    return 0;
  }

  //! return the internal partitioned petsc vec
  std::shared_ptr<PartitionedPetscVec<FunctionSpaceType, nComponents_>>
  partitionedPetscVec();

protected:
  std::shared_ptr<PartitionedPetscVec<FunctionSpaceType, nComponents_>>
      values_ = nullptr; //< Petsc vector containing the values, the values for
                         // the components are stored as struct of array, e.g.
                         // (comp1val1, comp1val2, comp1val3, ..., comp2val1,
                         // comp2val2, comp2val3,
                         //...). Dof ordering proceeds fastest over dofs of a
                         // node, then over nodes, node numbering is along whole
                         // domain, fastes in x, then in y,z direction.
};

} // namespace FieldVariable

#include "field_variable/structured/02_field_variable_data_structured.tpp"
