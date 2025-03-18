#pragma once

#include <Python.h> // has to be the first included header

#include "control/precice/nested_solver.h"

#include "precice/precice.hpp"

namespace Control {

/** Generic Precice adapter, can be configured to either prescribe Neumann or
 * Dirichlet boundary conditions.
 */
template <typename NestedSolver>
class PreciceAdapterInitialize
    : public PreciceAdapterNestedSolver<NestedSolver> {
public:
  //! constructor, gets the DihuContext object which contains all python
  //! settings
  PreciceAdapterInitialize(DihuContext context);
  // PreciceAdapterInitialize(){};

  //! initialize the object
  void initialize();

  using PreciceAdapterNestedSolver<NestedSolver>::FunctionSpace;
  // typedef typename NestedSolver::FunctionSpace VolumeFunctionSpace;

  /** a coupling mesh of precice, this does not have to coincide with an
   * opendihu mesh
   */
  struct PreciceSurfaceMesh {
    std::string preciceMeshName; //< name of the precice mesh as used in the
                                 // precice config XML file
    std::vector<int>
        preciceVertexIds; //< the vertex ids in precice of the geometry values
    std::vector<dof_no_t> dofNosLocal; //< the local dof nos in the 3D mesh of
                                       // the surface mesh nodes
    std::vector<dof_no_t> selectedDofNosLocal;
    std::vector<double>
        geometryValuesSurface; //< the geometry values, i.e., node positions

    int nNodesLocal; //< local number of nodes on the surface where the tendon
                     // is coupled
    enum {
      face2Minus,
      face2Plus
    } face; //< the face of the 3D mesh where the 2D coupling mesh is, 2- =
            // bottom, 2+ = top
  };

  struct PreciceVolumeMesh {
    std::string preciceMeshName; //< name of the precice mesh as used in the
                                 // precice config XML file
    std::vector<int>
        preciceVertexIds; //< the vertex ids in precice of the geometry values
    std::vector<dof_no_t> dofNosLocal; //< the local dof nos in the 3D mesh of
                                       // the surface mesh nodes
    std::vector<double>
        geometryValues;           //< the geometry values, i.e., node positions
    int nNodesLocal;              //< local number of nodes
    std::string opendihuMeshName; //< opendihu mesh name corresponding to the
                                  // precice mesh
  };

  /** a precice coupling participant to which the current solver is coupled to
   */
  struct PreciceSurfaceData {

    std::string displacementsName; //< precice name of the displacements
                                   // variable, if any
    std::string
        velocitiesName; //< precice name of the velocities variable, if any
    std::string tractionName; //< precice name of the traction variable, if any

    int preciceDataDimDisplacements; //< precice dim of the displacements
                                     // variable
    int preciceDataDimVelocities;    //< precice dim of the velocities variable
    int preciceDataDimTraction;      //< precice dim of the traction variable

    bool average;

    enum {
      ioRead,
      ioWrite
    } ioType; //< if this variable is to be written or read to other
              // participants over precice

    enum {
      bcTypeDirichlet,
      bcTypeNeumann
    } boundaryConditionType; //< the type of the boundary condition to set, if
                             // ioType == ioRead

    std::shared_ptr<PreciceSurfaceMesh>
        preciceMesh; //< the coupling mesh, this is derived from the option
                     // preciceMeshName
  };

  struct PreciceVolumeData {
    std::string preciceDataName; //< precice name of the variable, if any
    std::string slotName; //< slot name as given in config, this is used to
                          // determine slotNo
    int slotNo;           //< slot no that corresponds to this field variable

    std::string opendihuMeshName; //< opendihu mesh name that is used for the
                                  // geometry initialization
    bool isGeometryField; //< if the corresponding field variable is a geometry
                          // field

    enum {
      ioRead,
      ioWrite
    } ioType; //< if this variable is to be written or read to other
              // participants over precice

    std::shared_ptr<PreciceVolumeMesh>
        preciceMesh; //< the coupling mesh, this is derived from the option
                     // preciceMeshName
  };

  //! parse the options in "preciceMeshes" and store in variable
  //! preciceSurfaceMeshes_
  void initializePreciceSurfaceMeshes();

  //! initialize all meshes in precice from the variable preciceSurfaceMeshes_
  void setMeshesInPrecice();

  //! parse the options in "preciceSurfaceData" and initialize all variables in
  //! precice, store in variable preciceSurfaceData_
  void initializePreciceSurfaceData();

  //! initialize Dirichlet boundary conditions at all dofs that will get some
  //! prescribed values during coupling
  void initializeDirichletBoundaryConditions();

  bool inMuscleMeshTopA(const double x, const double y, const double z);

  DihuContext context_; //< object that contains the python config for the
                        // current context and the global singletons meshManager
                        // and solverManager
  PythonConfig specificSettings_; //< python object containing the value of the
                                  // python config dict with corresponding key

  NestedSolver
      nestedSolver_; //< the nested solver that is controlled by this class

  std::string preciceParticipantName_; //< name of the participant as given in
                                       // the precice config
  std::shared_ptr<precice::Participant>
      preciceParticipant_; //< the precice solver interface that makes all
                           // preCICE functionality accessible
  std::vector<std::shared_ptr<PreciceSurfaceMesh>>
      preciceSurfaceMeshes_; //< all surface coupling meshes
  std::vector<std::shared_ptr<PreciceVolumeMesh>>
      preciceVolumeMeshes_; //< all coupling meshes
  std::vector<PreciceSurfaceData>
      preciceSurfaceData_; //< all surface precice variables "data"
  std::vector<PreciceVolumeData>
      preciceVolumeData_; //< all precice variables "data"

  std::shared_ptr<
      typename PreciceAdapterNestedSolver<NestedSolver>::FunctionSpace>
      functionSpace_; //< the function space of the nested solver

  bool ownRankIsInvolved_; //< if the own rank has part of a coupling surface
                           // and is involved in the coupling
  bool couplingEnabled_;   //< if the coupling is enabled, if not it can be used
                           // for debugging, without precice
  double maximumPreciceTimestepSize_; //< maximum timestep size that precice
                                      // will allow for the current time step
  double timeStepWidth_;              //< timestep width of the solver
  int timeStepOutputInterval_; //< interval in which to output current time
  double scalingFactor_;       //< a factor to scale the exchanged data prior to
                               // communication
  bool outputOnlyConvergedTimeSteps_; //< option if the output should be written
                                      // only for converged timesteps
  double endTimeIfCouplingDisabled_;  //< the end time that is used if the
                                      // coupling is not enabled

  bool initialized_; //< if initialize() was already called
};

} // namespace Control

#include "control/precice/initialize.tpp"
