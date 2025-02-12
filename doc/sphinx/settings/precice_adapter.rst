PreciceAdapter
=================

`preCICE <https://precice.org/>_` is a coupling library for partitioned multi-physics simulations. 
Using preCICE allows us to couple OpenDiHu simulations to other OpenDiHu simulations, or also to simulations that use other softwares, e.g., `deal.II <https://www.dealii.org/>_` or FEBio. 
The adapter can be used to couple both surface data and volume data between participants (simulations). 

- An example of surface coupling would be a muscle-tendon system, where the muscle and the tendon, share a contact surface and exchange data via preCICE.
- An example of volume coupling would be using preCICE to couple a muscle's mechanics simulation with a muscle's electrophysiology simulation. 
Coupling muscle fibers with the mechanics is possible without preCICE, but using preCICE might be helpful on some cases, i.e., if we want to avoid restrictions on the partitioning.

preCICE configuration
----------------------

To run a simulation using preCICE you will need a `precice-config.xml <https://precice.org/configuration-overview.html>_` file. This file is used at run time by OpenDiHu and it is where the user can configure what data is coupled and how it is coupled. 
Here is an example for how we configure volume coupling between a participant named *PartitionedFibers* and a participant named *MuscleContraction*. In this case, the fibers send *Gamma* to the mechanics participant, which sends back *Geometry*. The coupling takes place using a *serial-explicit* scheme and the mapping between the meshes follows a kernel approach. 
Currently, *serial-explicit* is the only option for OpenDiHu participants that do not include a solid mechanics solver. 

.. code-block:: xml

  <?xml version="1.0"?>

  <precice-configuration>
      
    <data:vector name="Geometry"/>
    <data:scalar name="Gamma"/>

    <mesh name="PartitionedFibersMesh" dimensions="3">
      <use-data name="Geometry"/>
      <use-data name="Gamma"/>
    </mesh>

    <mesh name="MuscleContractionMesh" dimensions="3">
      <use-data name="Geometry"/>
      <use-data name="Gamma"/>
    </mesh>
    
    <participant name="PartitionedFibers">
      <provide-mesh name="PartitionedFibersMesh"/>
      <receive-mesh name="MuscleContractionMesh" from="MuscleContraction"/>
      
      <write-data name="Gamma" mesh="PartitionedFibersMesh"/>
      <read-data  name="Geometry" mesh="PartitionedFibersMesh"/>
      
      <mapping:rbf-pum-direct direction="read" from="MuscleContractionMesh" to="PartitionedFibersMesh" constraint="consistent" relative-overlap="0.15" vertices-per-cluster="20" project-to-input="false" polynomial="separate">
        <basis-function:compact-polynomial-c6 support-radius="1" />
      </mapping:rbf-pum-direct>
    </participant>
    
    <participant name="MuscleContraction">
      <provide-mesh name="MuscleContractionMesh"/>
      <receive-mesh name="PartitionedFibersMesh" from="PartitionedFibers"/>
      <write-data name="Geometry" mesh="MuscleContractionMesh"/>
      <read-data  name="Gamma"    mesh="MuscleContractionMesh"/>
      
      <mapping:nearest-neighbor constraint="consistent" direction="read" from="MuscleFibersMesh" to="MuscleMechanicsMesh"/>
      
    </participant>

    <m2n:sockets acceptor="PartitionedFibers" connector="MuscleContraction" network="lo" />

    <coupling-scheme:serial-explicit>
      <participants first="PartitionedFibers" second="MuscleContraction"/>
      <max-time value="1000.0"/>           <!-- end time of the whole simulation -->
      <time-window-size value="1e-1"/>   <!-- timestep width dt_3D -->
      
      <exchange data="Gamma"    mesh="PartitionedFibersMesh" from="PartitionedFibers" to="MuscleContraction"/>
      <exchange data="Geometry" mesh="MuscleContractionMesh" from="MuscleContraction" to="PartitionedFibers" initialize="yes"/> 
    </coupling-scheme:serial-explicit>
    
  </precice-configuration>


OpenDiHu configuration
----------------------

To use preCICE together with OpenDiHu, you will need to modify your solver's source file and settings file. 

In the source file, you simply have to add a wrapper to your nested solvers. You can find which nested solvers are supported by looking at the template specializations of the precice nested solver.

C++ code:

.. code-block:: c

  // surface coupling adapter
  Control::PreciceAdapter<
    /*nested solver*/
  >
  

In the python settings, you will also have to add the *"PreciceAdapter":* wrapper on top of your nested solver settings. Here we show how this would look like for a muscle that exchanges data with two tendons. The tendons are located at the face 2+ and 2- of the muscle.  The python settings if you want to do surface coupling are as follows:

.. code-block:: python

  "PreciceAdapter": {        # precice adapter for muscle
    "timeStepOutputInterval":   100,                        # interval in which to display current timestep and time in console
    "timestepWidth":            1,                          # coupling time step width, must match the value in the precice config
    "couplingEnabled":          variables.enable_coupling,  # if the precice coupling is enabled, if not, it simply calls the nested solver, for debugging
    "preciceConfigFilename":    "precice-config.xml",    # the preCICE configuration file
    "preciceParticipantName":   "MuscleSolver",             # name of the own precice participant, has to match the name given in the precice xml config file
    "scalingFactor":            1,                          # a factor to scale the exchanged data, prior to communication
    "outputOnlyConvergedTimeSteps": True,                   # if the output writers should be called only after a time window of precice is complete, this means the timestep has converged
    "preciceSurfaceMeshes": [                                      # the precice meshes get created as the top or bottom surface of the main geometry mesh of the nested solver
      {
        "preciceMeshName":      "MuscleMeshBottom",         # precice name of the 2D coupling mesh
        "face":                 "2-",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
      },
      {
        "preciceMeshName":      "MuscleMeshTop",           # precice name of the 2D coupling mesh
        "face":                 "2+",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
      }
    ],
    "preciceSurfaceData": [
      {
        "mode":                 "read-displacements-velocities",    # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "MuscleMeshBottom",                 # name of the precice coupling surface mesh, as given in the precice xml settings file
        "displacementsName":    "Displacement",                     # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
        "velocitiesName":       "Velocity",                         # name of the velocity "data", i.e. field variable, as given in the precice xml settings file
      },
      {
        "mode":                 "read-displacements-velocities",    # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "MuscleMeshTop",                   # name of the precice coupling surface mesh, as given in the precice xml settings file
        "displacementsName":    "Displacement",                     # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
        "velocitiesName":       "Velocity",                         # name of the velocity "data", i.e. field variable, as given in the precice xml settings file
      },
      {
        "mode":                 "write-traction",                   # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "MuscleMeshBottom",                 # name of the precice coupling surface mesh, as given in the precice xml settings 
        "tractionName":         "Traction",                         # name of the traction "data", i.e. field variable, as given in the precice xml settings file
      },
      {
        "mode":                 "write-traction",                   # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "MuscleMeshTop",                   # name of the precice coupling surface mesh, as given in the precice xml settings 
        "tractionName":         "Traction",                         # name of the traction "data", i.e. field variable, as given in the precice xml settings file
      }
    ],
    # options of the nested solver
    

If you want to do volume coupling instead, then the settings file would like that:

.. code-block:: python

  "PreciceAdapter": {
    "timeStepOutputInterval":   100,                        # interval in which to display current timestep and time in console
    "timestepWidth":            1,                          # coupling time step width, must match the value in the precice config
    "couplingEnabled":          True,                       # if the precice coupling is enabled, if not, it simply calls the nested solver, for debugging
    "endTimeIfCouplingDisabled": variables.end_time,        # if "couplingEnabled" is set to False, use this end time for the simulation
    "preciceConfigFilename":    "../precice_config.xml",    # the preCICE configuration file
    "preciceParticipantName":   "PartitionedFibers",        # name of the own precice participant, has to match the name given in the precice xml config file
    "scalingFactor":            1,                          # a factor to scale the exchanged data, prior to communication
    "outputOnlyConvergedTimeSteps": True,                   # if the output writers should be called only after a time window of precice is complete, this means the timestep has converged
    
    "preciceVolumeData": [
      {
        "mode":                 "read",                     # mode is one of "read" or "write"
        "preciceDataName":      "Geometry",                 # name of the vector or scalar to transfer, as given in the precice xml settings file
        "preciceMeshName":      "PartitionedFibersMesh",    # name of the precice coupling mesh, as given in the precice xml settings file
        "opendihuMeshName":     None,                       # extra specification of the opendihu mesh that is used for the initialization of the precice mapping. If None or "", the mesh of the field variable is used.
        "slotName":             None,                       # name of the existing slot of the opendihu data connector to which this variable is associated to (only relevant if not isGeometryField)
        "isGeometryField":      True,                       # if this is the geometry field of the mesh
      },
      {
        "mode":                 "write",                    # mode is one of "read" or "write"
        "preciceDataName":      "Gamma",                    # name of the vector or scalar to transfer, as given in the precice xml settings file
        "preciceMeshName":      "PartitionedFibersMesh",    # name of the precice coupling mesh, as given in the precice xml settings file
        "opendihuMeshName":     None,                       # extra specification of the opendihu mesh that is used for the initialization of the precice mapping. If None or "", the mesh of the field variable is used.
        "slotName":             "gamma",                    # name of the existing slot of the opendihu data connector to which this variable is associated to (only relevant if not isGeometryField)
        "isGeometryField":      False,                      # if this is the geometry field of the mesh
      },
    ],
    
    # options of the nested solver
  }


Note that if we do surface coupling, only a fraction of the opendihu mesh is exchanged, whereas when we do volume coupling, the whole mesh is exchanged. It is possible to perform surface and volume coupling simultaneously. 