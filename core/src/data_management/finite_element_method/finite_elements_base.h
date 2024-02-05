#pragma once

#include <Python.h> // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <tuple>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "function_space/mixed_function_space.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"
#include "slot_connection/slot_connector_data.h"

namespace Data {

/**  Base class storing data for all finite element computations, this mainly
 * includes the rhs and solution vectors and some matrices. nComponents is 1 for
 * normal scalar FE calculations and > 1 for structural mechanics
 */
template <typename FunctionSpaceType, int nComponents>
class FiniteElementsBase : public Data<FunctionSpaceType> {
public:
  //! type of data that will be transferred to nested solvers
  typedef SlotConnectorData<FunctionSpaceType, nComponents>
      SlotConnectorDataType;

  //! constructor
  FiniteElementsBase(DihuContext context);

  //! destructor
  virtual ~FiniteElementsBase();

  //! initialize the object, create all stored data
  virtual void initialize() override;

  //! reset the object and deallocate matrices
  virtual void reset() override;

  //! return reference to a right hand side vector, the PETSc Vec can be
  //! obtained via fieldVariable.valuesGlobal()
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>
  rightHandSide();

  //! return reference to solution of the system, the PETSc Vec can be obtained
  //! via fieldVariable.valuesGlobal()
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>
  solution();

  //! return reference to rhsNeumannBoundaryConditions, the PETSc Vec can be
  //! obtained via fieldVariable.valuesGlobal()
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>
  negativeRightHandSideNeumannBoundaryConditions();

  //! set the field variable rightHandSideNeumannBoundaryConditions
  void setNegativeRightHandSideNeumannBoundaryConditions(
      std::shared_ptr<
          FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>
          rightHandSideNeumannBoundaryConditions);

  //! print all stored data to stdout
  void print();

  //! create PETSc matrix
  void initializeMassMatrix();

  //! create the inverse of the lumped mass matrix
  void initializeInverseLumpedMassMatrix();

  //! set the solution variable if it is initialized externally, such as in a
  //! timestepping scheme
  void setSolutionVariable(
      std::shared_ptr<
          FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>
          solution);

  //! return reference to the stiffness matrix
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix();

  //! return reference to the stiffness matrix, this matrix stays the initial
  //! stifness matrix and will not be changed for boundary conditions
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>>
  stiffnessMatrixWithoutBc();

  //! get the mass matrix
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> massMatrix();

  //! get the inversed lumped mass matrix
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>>
  inverseLumpedMassMatrix();

  //! get maximum number of expected non-zeros in stiffness matrix
  static void getPetscMemoryParameters(int &nNonZerosDiagonal,
                                       int &nNonZerosOffdiagonal);

  //! get the data that will be transferred in the operator splitting to the
  //! other term of the splitting the transfer is done by the
  //! slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
      std::shared_ptr<
          FieldVariable::FieldVariable<FunctionSpaceType, 3>>, // geometry
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,
                                                   nComponents>>, // solution
      std::shared_ptr<
          FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>, // rhs
      std::shared_ptr<FieldVariable::FieldVariable<
          FunctionSpaceType, nComponents>> // neumann BC rhs
      >
      FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:
  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();

  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>>
      stiffnessMatrixWithoutBc_; //< the standard stiffness matrix of the finite
                                 // element formulation, without Dirichlet
                                 // boundary conditions applied
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>>
      stiffnessMatrix_; //< the standard stiffness matrix of the finite element
                        // formulation, dofs with Dirichlet BCs will get the
                        // columns and rows cleared and diagonal set to 1.
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>>
      massMatrix_; //< the standard mass matrix, which is a matrix that, applied
                   // to a rhs vector f, gives the rhs vector in weak
                   // formulation
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>>
      inverseLumpedMassMatrix_; //< the inverse lumped mass matrix that has only
                                // entries on the diagonal, they are the
                                // reciprocal of the row sums of the mass matrix

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>
      rhs_; //< the rhs vector in weak formulation
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>
      negativeRhsNeumannBoundaryConditions_; //< the rhs vector in weak
                                             // formulation, only contribution
                                             // from neumann boundary conditions
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType, nComponents>>
      solution_; //< the vector of the quantity of interest, e.g. displacement

  std::shared_ptr<SlotConnectorDataType>
      slotConnectorData_; //< the object that holds all slot connector
                          // components of field variables
};

} // namespace Data

#include "data_management/finite_element_method/finite_elements_base.tpp"
