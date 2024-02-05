#pragma once

#include <Python.h> // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "field_variable/field_variable.h"

namespace Data {

/**  The datastructures used for streamline tracer.
 *   BaseDataType is a Data class that provides the solution field variable for
 * the streamline tracer to operate on.
 */
template <typename FunctionSpaceType, typename BaseDataType>
class StreamlineTracer : public Data<FunctionSpaceType> {
public:
  //! constructor
  StreamlineTracer(DihuContext context);

  //! destructur
  ~StreamlineTracer();

  //! return a reference to the solution vector, the PETSc Vec can be obtained
  //! via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
  gradient();

  //! print all stored data to stdout
  virtual void print();

  //! set the data object that holds the field frorm which stream lines are
  //! generated
  void setBaseData(std::shared_ptr<BaseDataType> baseData);

  //! return the total number of degrees of freedom, this can be a multiple of
  //! the number of nodes of the mesh
  virtual dof_no_t nNodesLocalWithGhosts();

  //! return the total number of degrees of freedom, this can be a multiple of
  //! the number of nodes of the mesh
  virtual dof_no_t nNodesLocalWithoutGhosts();

  //! create a fiber mesh from the given node positions, store it in mesh
  //! manager and store a pointer to the geometry field in fiberGeometry
  void createFiberMesh(const std::vector<Vec3> &nodePositions);

  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>,
                                       BasisFunction::LagrangeOfOrder<1>>
      FunctionSpaceFiber;
  typedef FieldVariable::FieldVariable<FunctionSpaceFiber, 3>
      FieldVariableFiberGeometry;

  //! field variables that will be output by outputWriters
  typename BaseDataType::FieldVariablesForOutputWriter dummy;
  typedef decltype(std::tuple_cat(
      dummy,
      std::tuple<std::shared_ptr<FieldVariable::FieldVariable<
                     FunctionSpaceType, 3>>, // gradient field
                 std::vector<std::shared_ptr<
                     FieldVariableFiberGeometry>> // geometry fields of meshes
                 >())) FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

protected:
  //! initializes the vectors with size
  virtual void createPetscObjects();

  std::shared_ptr<BaseDataType>
      baseData_; //< the data object that holds the field frorm which stream
                 //lines are generated
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, 3>>
      gradient_; //< the gradient field of the solution field variable

  std::vector<std::shared_ptr<FieldVariableFiberGeometry>>
      fiberGeometry_; //< geometry fields of fibers

  int fiberNo_ = 0; //< counter for names of generated fibers
};

} // namespace Data

#include "data_management/streamline_tracer.tpp"
