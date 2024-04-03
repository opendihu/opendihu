import sys
import os
import importlib


# parse rank arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

#add variables subfolder to python path where the variables script is located
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
    print("Warning: There is no variables file, e.g:\n ./tendon ../settings_tendon.py tendon.py\n")
  exit(0)


# define config
config = {
  "scenarioName":                   variables.scenario_name,      # scenario name to identify the simulation runs in the log file
  "logFormat":                      "csv",                        # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",       # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes_log.txt",    # log file for mappings 
  "Meshes":                         variables.meshes,
  
  "PreciceAdapter": {        
      "timeStepOutputInterval":   1,                        # interval in which to display current timestep and time in console
      "timestepWidth":            variables.dt_3D,                          # coupling time step width, must match the value in the precice config
      "couplingEnabled":          True,                       # if the precice coupling is enabled, if not, it simply calls the nested solver, for debugging
      "preciceConfigFilename":    "../precice_config.xml",    # the preCICE configuration file
      "preciceParticipantName":   "Tendon",             # name of the own precice participant, has to match the name given in the precice xml config file
      "preciceSurfaceMeshes": [                                      # the precice meshes get created as the top or bottom surface of the main geometry mesh of the nested solver
        {
          "preciceMeshName":      "SurfaceTendonMesh",            # precice name of the 2D coupling mesh
          "face":                 "2+",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
        }   
      ],
      "preciceSurfaceData": [  
        {
          "mode":                 "write-displacements-velocities",   # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
          "preciceMeshName":      "SurfaceTendonMesh",                    # name of the precice coupling surface mesh, as given in the precice xml settings file
          "displacementsName":    "Displacement",                     # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
          "velocitiesName":       "Velocity",                     # name of the velocities "data", i.e. field variable, as given in the precice xml settings file

        },
        {
          "mode":                 "read-traction",                    # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
          "preciceMeshName":      "SurfaceTendonMesh",                    # name of the precice coupling surface mesh, as given in the precice xml settings 
          "tractionName":         "Traction",                         # name of the traction "data", i.e. field variable, as given in the precice xml settings file
        }
      ],

  "DynamicHyperelasticitySolver": {
      "timeStepWidth":              variables.dt_3D,      #variables.dt_3D,      # time step width 
      "endTime":                    variables.end_time,           # end time of the simulation time span    
      "durationLogKey":             "duration_mechanics",         # key to find duration of this solver in the log file
      "timeStepOutputInterval":     1,                            # how often the current time step should be printed to console
      
      "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
      "density":                    variables.rho,  
      "displacementsScalingFactor": 1.0,                          # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
      "residualNormLogFilename":    "log_residual_norm.txt",      # log file where residual norm values of the nonlinear solver will be written
      "useAnalyticJacobian":        True,                         # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
      "useNumericJacobian":         False,                        # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
      
      "dumpDenseMatlabVariables":   False,                        # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
      # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
      
      # mesh
      "meshName":                   "tendonMesh3D",           # mesh with quadratic Lagrange ansatz functions
      "inputMeshIsGlobal":          True,                         # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
      
      "fiberMeshNames":             [],                           # fiber meshes that will be used to determine the fiber direction
      "fiberDirection":             [0,0,1],                      # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
      
      # nonlinear solver
      "relativeTolerance":          variables.linear_relative_tolerance,                         # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":          variables.linear_absolute_tolerance,                        # 1e-10 absolute tolerance of the residual of the linear solver       
      "solverType":                 variables.elasticity_solver_type,                    # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
      "preconditionerType":         variables.elasticity_preconditioner_type,                         # type of the preconditioner
      "maxIterations":              variables.linear_max_iterations,                          # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": variables.snes_max_function_evaluations,                          # maximum number of function iterations
      "snesMaxIterations":          variables.snes_max_iterations,                           # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance":      variables.snes_relative_tolerance,                         # relative tolerance of the nonlinear solver
      "snesLineSearchType":         variables.snes_linear_search,                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesAbsoluteTolerance":      variables.snes_absolute_tolerance,                         # absolute tolerance of the nonlinear solver
      "snesRebuildJacobianFrequency":   variables.snes_rebuild_jacobian_frequency,      "relativeTolerance":          variables.linear_relative_tolerance,                         # 1e-10 relative tolerance of the linear solver
      
      #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),          # dump system matrix and right hand side after every solve
      "dumpFilename":               "",                           # dump disabled
      "dumpFormat":                 "matlab",                     # default, ascii, matlab
      
      #"loadFactors":                [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
      #"loadFactors":                [0.5, 1.0],                   # load factors for every timestep
      "loadFactors":                [],                           # no load factors, solve problem directly
      "loadFactorGiveUpThreshold":  1e-3, 
      "nNonlinearSolveCalls":       1,                            # how often the nonlinear solve should be called
      
      # boundary and initial conditions
      "dirichletBoundaryConditions": variables.dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
      "neumannBoundaryConditions": variables.neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
      "divideNeumannBoundaryConditionValuesByTotalArea": False,    # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
      "updateDirichletBoundaryConditionsFunction": None, #update_dirichlet_bc,   # function that updates the dirichlet BCs while the simulation is running
      "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # stide every which step the update function should be called, 1 means every time step
      "updateNeumannBoundaryConditionsFunction": None,       # a callback function to periodically update the Neumann boundary conditions
      "updateNeumannBoundaryConditionsFunctionCallInterval": 1,           # every which step the update function should be called, 1 means every time step 
      
        "initialValuesDisplacements": [[0, 0, 0] for _ in range(variables.bs_x * variables.bs_y * variables.bs_z)],
        "initialValuesVelocities":    [[0, 0, 0] for _ in range(variables.bs_x * variables.bs_y * variables.bs_z)],
        "constantBodyForce":          (0, 0, 0),
    
        "extrapolateInitialGuess":     True, 

      "dirichletOutputFilename":     "out/" + variables.scenario_name + "_dirichlet_boundary_conditions_tendon",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable

      "OutputWriter" : [
          {"format": "Paraview", "outputInterval": variables.output_interval, "filename": "out/" + variables.scenario_name + "/tendon", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
      ],
      "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
          "OutputWriter" : [
              # {"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ]
      },
      "dynamic": {    # output of the dynamic solver, has additional virtual work values 
          "OutputWriter" : [   # output files for displacements function space (quadratic elements)
          
          ],
      },
      "LoadIncrements": {   
          "OutputWriter" : [
              #{"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
          ]
      },
    }
  }
}
