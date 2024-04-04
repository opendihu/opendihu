import sys
import os
import importlib

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# add folders to python path
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables

# if first argument contains "*.py", it is a custom variable definition file, load these values
if ".py" in sys.argv[0]:
  variables_path_and_filename = sys.argv[0]
  variables_path,variables_filename = os.path.split(variables_path_and_filename)  # get path and filename 
  sys.path.insert(0, os.path.join(script_path,variables_path))                    # add the directory of the variables file to python path
  variables_module,_ = os.path.splitext(variables_filename)                       # remove the ".py" extension to get the name of the module
  
  if rank_no == 0:
    print("Loading variables from \"{}\".".format(variables_path_and_filename))
    
  custom_variables = importlib.import_module(variables_module, package=variables_filename)    # import variables module
  variables.__dict__.update(custom_variables.__dict__)
  sys.argv = sys.argv[1:]     # remove first argument, which now has already been parsed
else:
  if rank_no == 0:
    print("Warning: There is no variables file, e.g:\n ./muscle ../settings_muscle_mechanics.py muscle.py\n")
  exit(0)

# define config
config = {
  "scenarioName":                   variables.scenario_name,

  "logFormat":                      "csv",
  "mappingsBetweenMeshesLogFile":   "out/" + variables.scenario_name + "/muscle_mappings_between_meshes_log.txt",
  "solverStructureDiagramFile":     "out/" + variables.scenario_name + "/muscle_solver_structure_diagram.txt",

  "Meshes":                         variables.meshes,
  "MappingsBetweenMeshes":          {},

  "Solvers": {
    "mechanicsSolver": {
      "solverType":                     "preonly",
      "preconditionerType":             "lu",
      "relativeTolerance":              1e-10,
      "absoluteTolerance":              1e-10,
      "maxIterations":                  1e4,
      "snesLineSearchType":             "l2",
      "snesRelativeTolerance":          1e-5,
      "snesAbsoluteTolerance":          1e-5,
      "snesMaxIterations":              10,
      "snesMaxFunctionEvaluations":     1e8,
      "snesRebuildJacobianFrequency":   5,
      "dumpFilename":                   "",
      "dumpFormat":                     "matlab"
    }
  },

  "PreciceAdapter": {
    "preciceConfigFilename":        "../precice_config.xml",
    "preciceParticipantName":       "Muscle",
    "couplingEnabled":              True,
    "timeStepOutputInterval":       100,
    "timestepWidth":                variables.dt_3D,
    "scalingFactor":                1,
    "outputOnlyConvergedTimeSteps": True,

    "preciceVolumeData": [
      {
        "mode":             "write",
        "preciceDataName":  "Geometry",
        "preciceMeshName":  "VolumeMuscleMesh",
        "opendihuMeshName": None,
        "slotName":         None,
        "isGeometryField":  True
      },
      {
        "mode":             "read",
        "preciceDataName":  "Gamma",
        "preciceMeshName":  "VolumeMuscleMesh",
        "opendihuMeshName": None,
        "slotName":         "gamma",
        "isGeometryField":  False
       }
    ],
    "preciceSurfaceMeshes": [                                      # the precice meshes get created as the top or bottom surface of the main geometry mesh of the nested solver
      {
        "preciceMeshName":      "SurfaceMuscleMesh",            # precice name of the 2D coupling mesh
        "face":                 "2-",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
      }   
    ],
    "preciceSurfaceData": [  
      {
        "mode":                 "read-displacements-velocities",   # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "SurfaceMuscleMesh",                    # name of the precice coupling surface mesh, as given in the precice xml settings file
        "displacementsName":    "Displacement",                     # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
        "velocitiesName":       "Velocity",                     # name of the velocities "data", i.e. field variable, as given in the precice xml settings file

      },
      {
        "mode":                 "write-traction",                    # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "SurfaceMuscleMesh",                    # name of the precice coupling surface mesh, as given in the precice xml settings 
        "tractionName":         "Traction",                         # name of the traction "data", i.e. field variable, as given in the precice xml settings file
      }
    ],

    "MuscleContractionSolver": {
      "Pmax":                         variables.pmax,
      "slotNames":                    ["lambda", "ldot", "gamma", "T"],
      "dynamic":                      True,

      "numberTimeSteps":              1,
      "timeStepOutputInterval":       100,
      "lambdaDotScalingFactor":       1,
      "enableForceLengthRelation":    True,
      "mapGeometryToMeshes":          [],

      "OutputWriter": [
        {
          "format":             "Paraview",
          "outputInterval":     int(1.0 / variables.dt_3D * variables.output_interval),
          "filename":           "out/" + variables.scenario_name + "/muscle_mechanics",
          "fileNumbering":      "incremental", 
          "binary":             True,
          "fixedFormat":        False,
          "onlyNodalValues":    True,
          "combineFiles":       True
        }
      ],

      "DynamicHyperelasticitySolver": { 
        "durationLogKey":         "duration_3D",
        "logTimeStepWidthAsKey":  "dt_3D",
        "numberTimeSteps":        1,
        "materialParameters":     variables.material_parameters,
        "density":                variables.rho,
        "timeStepOutputInterval": 1,

        "meshName":                 "muscleMesh3D",
        "fiberDirectionInElement":  variables.fiber_direction,
        "inputMeshIsGlobal":        True,
        "fiberMeshNames":           [],
        "fiberDirection":           None,

        "solverName":                 "mechanicsSolver",
        "displacementsScalingFactor": 1.0,
        "useAnalyticJacobian":        True,
        "useNumericJacobian":         False,
        "dumpDenseMatlabVariables":   False,
        "loadFactorGiveUpThreshold":  4e-2,
        "loadFactors":                [],
        "scaleInitialGuess":          False,
        "extrapolateInitialGuess":    True,
        "nNonlinearSolveCalls":       1,

        "dirichletBoundaryConditions":                            variables.dirichlet_bc,
        "neumannBoundaryConditions":                              variables.neumann_bc,
        "updateDirichletBoundaryConditionsFunction":              None,
        "updateDirichletBoundaryConditionsFunctionCallInterval":  1,
        "divideNeumannBoundaryConditionValuesByTotalArea":        True,

        "initialValuesDisplacements": [[0, 0, 0] for _ in range(variables.bs_x * variables.bs_y * variables.bs_z)],
        "initialValuesVelocities":    [[0, 0, 0] for _ in range(variables.bs_x * variables.bs_y * variables.bs_z)],
        "constantBodyForce":          (0, 0, 0),

        "dirichletOutputFilename":    "out/" + variables.scenario_name + "/muscle_dirichlet_output",
        "residualNormLogFilename":    "out/" + variables.scenario_name + "/muscle_residual_norm_log.txt",
        "totalForceLogFilename":      "out/" + variables.scenario_name + "/muscle_total_force_log.txt",

        "OutputWriter": [],
        "pressure":       { "OutputWriter": [] },
        "dynamic":        { "OutputWriter": [] },
        "LoadIncrements": { "OutputWriter": [] }
      }
    }
  }
}