# sensory organs, interneurons, motoneurons, multidomain mockup (no electrophyisiology), mechanics
#
# [muscle spindle model] -> [motor neuron model] <- mapping by callback <- [interneuron model] <- mapping by callback <- [golgi tendon organ model]
#                            +-> electro-mechanics model

import numpy as np
import scipy.stats
import pickle
import sys,os 
import timeit
import argparse
import importlib
import struct
sys.path.insert(0, "..")

# own MPI rank no and number of MPI ranks
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# add variables subfolder to python path where the variables script is located
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_path)
sys.path.insert(0, os.path.join(script_path,'variables'))

import variables              # file variables.py, defined default values for all parameters, you can set the parameters there  
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain

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
    print("Warning: There is no variables file, e.g:\n ./neurons_with_contraction ../settings_neurons_with_contraction.py spindles.py\n")
  exit(0)
  
# -------------------------------------------------------- begin parameters ---------------------------------------------------------
# -------------------------------------------------------- end parameters ------------------------------------------------------------
parser = argparse.ArgumentParser(description='neurons_with_contraction')
parser.add_argument('--scenario_name',                       help='The name to identify this run in the log.',   default=variables.scenario_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--potential_flow_solver_type',          help='The solver for the potential flow (non-spd matrix).', default=variables.potential_flow_solver_type, choices=["gmres","cg","lu","gamg","richardson","chebyshev","cholesky","jacobi","sor","preonly"])
parser.add_argument('--potential_flow_preconditioner_type',  help='The preconditioner for the potential flow.',  default=variables.potential_flow_preconditioner_type, choices=["jacobi","sor","lu","ilu","gamg","none"])
parser.add_argument('--paraview_output',                     help='Enable the paraview output writer.',          default=variables.paraview_output, action='store_true')
parser.add_argument('--adios_output',                        help='Enable the MegaMol/ADIOS output writer.',          default=variables.adios_output, action='store_true')
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('--firing_times_file',                   help='The filename of the file that contains the cellml model.', default=variables.firing_times_file)
parser.add_argument('--end_time', '--tend', '-t',            help='The end simulation time.',                    type=float, default=variables.end_time)
parser.add_argument('--output_timestep_multidomain',         help='The timestep for writing outputs.',           type=float, default=variables.output_timestep_multidomain)
parser.add_argument('--multidomain_solver_type',             help='Solver for multidomain system.',              default=variables.multidomain_solver_type)
parser.add_argument('--multidomain_preconditioner_type',     help='Preconditioner for multidomain system.',      default=variables.multidomain_preconditioner_type)
parser.add_argument('--multidomain_relative_tolerance',      help='relative residual tol. for multidomain solver.',           type=float, default=variables.multidomain_relative_tolerance)
parser.add_argument('--multidomain_absolute_tolerance',      help='absolute residual tol. for multidomain solver.',           type=float, default=variables.multidomain_absolute_tolerance)
parser.add_argument('--theta',                               help='parameter of Crank-Nicholson in multidomain system.',      type=float, default=variables.theta)
parser.add_argument('--use_lumped_mass_matrix',              help='If the formulation with lumbed mass matrix should be used.',           default=variables.use_lumped_mass_matrix, action='store_true')
parser.add_argument('--use_symmetric_preconditioner_matrix', help='If the preconditioner matrix should be symmetric (only diagnoal part of system matrix).',  default=variables.use_symmetric_preconditioner_matrix, action='store_true')
parser.add_argument('--initial_guess_nonzero',               help='If the initial guess to the linear solver should be the last solution.',  default=variables.initial_guess_nonzero, action='store_true')

# parse command line arguments and assign values to variables module
args, other_args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)
if len(other_args) != 0 and rank_no == 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)

# initialize some dependend variables
if variables.n_subdomains is not None:
  variables.n_subdomains_x = variables.n_subdomains[0]
  variables.n_subdomains_y = variables.n_subdomains[1]
  variables.n_subdomains_z = variables.n_subdomains[2]
  
variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z

# automatically initialize partitioning if it has not been set
if n_ranks != variables.n_subdomains:
  
  # create all possible partitionings to the given number of ranks
  optimal_value = n_ranks**(1/3)
  possible_partitionings = []
  for i in range(1,n_ranks+1):
    for j in [1]:
      if i*j <= n_ranks and n_ranks % (i*j) == 0:
        k = (int)(n_ranks / (i*j))
        performance = (k-optimal_value)**2 + (j-optimal_value)**2 + 1.1*(i-optimal_value)**2
        possible_partitionings.append([i,j,k,performance])
        
  # if no possible partitioning was found
  if len(possible_partitionings) == 0:
    if rank_no == 0:
      print("\n\n\033[0;31mError! Number of ranks {} does not match given partitioning {} x {} x {} = {} and no automatic partitioning could be done.\n\n\033[0m".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
    quit()
    
  # select the partitioning with the lowest value of performance which is the best
  lowest_performance = possible_partitionings[0][3]+1
  for i in range(len(possible_partitionings)):
    if possible_partitionings[i][3] < lowest_performance:
      lowest_performance = possible_partitionings[i][3]
      variables.n_subdomains_x = possible_partitionings[i][0]
      variables.n_subdomains_y = possible_partitionings[i][1]
      variables.n_subdomains_z = possible_partitionings[i][2]

# output information of run
if rank_no == 0:
  print("scenario_name: {},  n_subdomains: {} {} {},  n_ranks: {},  end_time: {}".format(variables.scenario_name, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, n_ranks, variables.end_time))
  print("dt_0D:           {:0.0e}    multidomain solver:         {}, lumped mass matrix: {}".format(variables.dt_0D, variables.multidomain_solver_type, variables.use_lumped_mass_matrix))
  print("dt_multidomain:  {:0.0e}    multidomain preconditioner: {}, symmetric precond.: {}".format(variables.dt_multidomain, variables.multidomain_preconditioner_type, variables.use_symmetric_preconditioner_matrix))
  print("dt_splitting:    {:0.0e}    theta: {}, solver tolerances, abs: {}, rel: {}".format(variables.dt_splitting, variables.theta, variables.multidomain_absolute_tolerance, variables.multidomain_relative_tolerance))
  print("dt_elasticity:   {:0.0e}    elasticity solver: {}, preconditioner: {}".format(variables.dt_elasticity, variables.elasticity_solver_type, variables.elasticity_preconditioner_type))
  print("fiber_file:              {}".format(variables.fiber_file))
  print("fat_mesh_file:           {}".format(variables.fat_mesh_file))
  print("cellml_file:             {}".format(variables.cellml_file))
  print("firing_times_file:       {}".format(variables.firing_times_file))
  print("********************************************************************************")
  
  # start timer to measure duration of parsing of this script  
  t_start_script = timeit.default_timer()
    
# initialize all helper variables
from helper import *


# callback function that receives the whole result values and produces plots while the simulation is running
def handle_result(n_instances, time_step_no, current_time, states, algebraics, name_information, additional_argument):
    
  print("handle result of muscle spindle")
  #print(name_information)
    
  # asign some states to variables
  Vm = states[name_information["stateNames"].index("membrane/V")]
  print("Vm: {}".format(Vm))
    
# add neuron meshes
variables.meshes.update(
{
  "motoneuronMesh": {
    "nElements" :         variables.n_motoneurons-1 if n_ranks == 1 else variables.n_motoneurons*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "motoneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "interneuronMesh": {
    "nElements" :         variables.n_interneurons-1 if n_ranks == 1 else variables.n_interneurons*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "interneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "muscleSpindleMesh": {
    "nElements" :         variables.n_muscle_spindles-1 if n_ranks == 1 else variables.n_muscle_spindles*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "muscle_spindle",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "muscleSpindleAndInterneuronMesh": {
    "nElements" :         variables.n_muscle_spindles+variables.n_interneurons-1 if n_ranks == 1 else (variables.n_muscle_spindles+variables.n_interneurons)*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "muscle_spindle_interneuron",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
  "golgiTendonOrganMesh": {
    "nElements" :         variables.n_golgi_tendon_organs-1 if n_ranks == 1 else variables.n_golgi_tendon_organs*n_ranks,
    "physicalExtent":     0,
    "physicalOffset":     0,
    "logKey":             "golgi_tendon_organ",
    "inputMeshIsGlobal":  True,
    "nRanks":             n_ranks
  },
})
  
config = {
  "scenarioName":          variables.scenario_name,
  "logFormat":                      "csv",                                # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",               # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes_log.txt",    # log file for mappings 
  "meta": {                                                               # additional fields that will appear in the log
    "partitioning":         [variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z]
  },
  "Meshes":                variables.meshes,
  "MappingsBetweenMeshes": {
    "3Dmesh": [
       {"name": "3Dmesh_elasticity_quadratic",                                 "xiTolerance": 1.5, "enableWarnings": True, "compositeUseOnlyInitializedMappings": True, "fixUnmappedDofs": True},
    ]
  },
  "Solvers": {

    "mechanicsSolver": {   # solver for the dynamic mechanics problem
      "relativeTolerance":   variables.snes_relative_tolerance,           # 1e-10 relative tolerance of the linear solver
      "absoluteTolerance":   variables.snes_absolute_tolerance,           # 1e-10 absolute tolerance of the residual of the linear solver
      "solverType":          variables.elasticity_solver_type,            # type of the linear solver
      "preconditionerType":  variables.elasticity_preconditioner_type,    # type of the preconditioner
      "maxIterations":       1e4,                                         # maximum number of iterations in the linear solver
      "snesMaxFunctionEvaluations": 1e8,                                  # maximum number of function iterations
      "snesMaxIterations":   variables.snes_max_iterations,               # maximum number of iterations in the nonlinear solver
      "snesRelativeTolerance": variables.snes_relative_tolerance,         # relative tolerance of the nonlinear solver
      "snesAbsoluteTolerance": variables.snes_absolute_tolerance,         # absolute tolerance of the nonlinear solver
      "snesLineSearchType": "l2",                                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
      "snesRebuildJacobianFrequency": variables.snes_rebuild_jacobian_frequency,    # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
      "hypreOptions":        "",                                          # additional options for the hypre solvers could be given here
      "dumpFilename":        "",                                          # dump system matrix and right hand side after every solve
      "dumpFormat":          "matlab",                                    # default, ascii, matlab
    }
  },
  
  # total solver
  "Coupling": {
    "description":            "overall solver",    # description that will be shown in solver structure visualization
    "endTime":                variables.end_time,
    "timeStepWidth":          variables.dt_stimulation_check,
    "logTimeStepWidthAsKey":  "dt_neurons",
    "durationLogKey":         "duration_total",
    "timeStepOutputInterval": 100,
    "connectedSlotsTerm1To2": {},          # connect motoneuron input
    "connectedSlotsTerm2To1": {8:1},    # connect input for muscle spindles and golgi tendon organs

    # muscle spindles, golgi tendon organs and interneurons
    "Term1": {
      # muscle spindles, golgi tendon organs and interneurons = motoneuron input solver

                   
      # muscle spindle model solver
      "Heun" : {
        "description":                  "muscle spindle",
        "timeStepWidth":                variables.dt_muscle_spindles,
        "logTimeStepWidthAsKey":        "dt_muscle_spindles",
        "durationLogKey":               "duration_muscle_spindles",
        "initialValues":                [],
        "timeStepOutputInterval":       1e4,
        "inputMeshIsGlobal":            True,
        "dirichletBoundaryConditions":  {},
        "dirichletOutputFilename":      None,                                 # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
        "checkForNanInf":               True,             # check if the solution vector contains nan or +/-inf values, if yes, an error is printed. This is a time-consuming check.
        "nAdditionalFieldVariables":    0,
        "additionalSlotNames":          [],
            
        # cellml model of muscle spindle
        "CellML" : {
          "modelFilename":                          variables.muscle_spindle_cellml_file,           # input C++ source file or cellml XML file
          "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
          "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
          
          # optimization parameters
          "optimizationType":                       "vc",                                           # "vc", "simd", "openmp" type of generated optimizated source file
          "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
          "compilerFlags":                          "-fPIC -O3 -march=native -shared ",             # compiler flags used to compile the optimized model code
          "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
          
          # stimulation callbacks, motor neuron is not stimulated by a callback function, but has a constant stimulation current
          "setSpecificStatesFunction":              None,                                           # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
          #"setSpecificStatesCallInterval":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
          "setSpecificStatesCallInterval":          0,                                              # 0 means disabled
          "setSpecificStatesCallFrequency":         0,                                              # set_specific_states should be called variables.stimulation_frequency times per ms
          "setSpecificStatesFrequencyJitter":       0,                                              # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
          "setSpecificStatesRepeatAfterFirstCall":  0.01,                                           # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
          "setSpecificStatesCallEnableBegin":       0,                                              # [ms] first time when to call setSpecificStates
          "additionalArgument":                     None,
                                
          "handleResultFunction":                   None, #handle_result,              # callback function that gets all current values and can do something with them
          "handleResultCallInterval":               1,                          # interval in which handle_result will be called
          "handleResultFunctionAdditionalParameter": None,                      # additional last argument for handle_result
          
          "mappings":                               variables.muscle_spindle_mappings,              # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
          "parametersInitialValues":                variables.muscle_spindle_parameters_initial_values,  # initial values for the parameters, either once for all instances, or for all instances in array of struct ordering with nParameters_ parameters per instance: [inst0p0, inst0p1, ... inst0pn, inst1p0, inst1p1, ...]
          
          "meshName":                               "muscleSpindleMesh",                                       
          "stimulationLogFilename":                 "out/stimulation.log",

          # output writer for states, algebraics and parameters                
          "OutputWriter" : [
            {"format": "PythonFile", "outputInterval": (int)(2./variables.dt_muscle_spindles*variables.output_timestep_neurons), "filename": "out/muscle_spindles", "binary": True, "fixedFormat": False, "combineFiles": True, "onlyNodalValues": True, "fileNumbering": "incremental"}
          ]
        }
      }
          
        },
        
    
    # motoneuron with electro-mechanics
    "Term2": {
      
      # motoneuron with electro-mechanics

      # map from λ in the 3D mesh to muscle spindles input
      "MapDofs": {
        "description":                "muscle_spindles_input",        # description that will be shown in solver structure visualization
        "nAdditionalFieldVariables":  1,                              # number of additional field variables that are defined by this object. They have 1 component, use the templated function space and mesh given by meshName.
        "additionalSlotNames":        [],
        "meshName":                   "muscleSpindleMesh",            # the mesh on which the additional field variables will be defined
        "beforeComputation": [                                        # transfer/mapping of dofs that will be performed before the computation of the nested solver
          {                                                 
            "fromConnectorSlot":              2,
            "toConnectorSlots":               8,
            "fromSlotConnectorArrayIndex":    0,                    # which fiber/compartment, this does not matter here because all compartment meshes have the same displacements
            "toSlotConnectorArrayIndex":      0,
            "mode":                             "callback",           # "copyLocal", "copyLocalIfPositive", "localSetIfAboveThreshold" or "communicate"
            "fromDofNosNumbering":              "global",
            "toDofNosNumbering":                "local",
            "dofsMapping":                      None,
            "inputDofs":                        muscle_spindle_node_nos,          # nodes are set at bottom in helper.py
            "outputDofs":                       [list(range(variables.n_muscle_spindles))],   # [0,1,...,n_muscle_spindles]
            "callback":                         variables.callback_muscle_spindles_input,
            #"thresholdValue":                   20,                  # if mode is "localSetIfAboveThreshold", this is the threshold, if the value is above it, set the value `valueToSet`
            #"valueToSet":                       20,                  # if mode is "localSetIfAboveThreshold", this is the value to set the target dof to, if the source dof is above thresholdValue.
          }
        ],
        "afterComputation":  None,  
            
                  
      # electro-mechanics solver
        "Coupling": {
          "description":            "electro-mechanics solver",   # description that will be shown in solver structure visualization
          "timeStepWidth":          variables.dt_elasticity,  # 1e-1
          "logTimeStepWidthAsKey":  "dt_elasticity",
          "durationLogKey":         "duration_total",
          "timeStepOutputInterval": 1,
          "endTime":                variables.end_time,
          "connectedSlotsTerm1To2": {1:2},          # transfer gamma to MuscleContractionSolver, the receiving slots are λ, λdot, γ
          "connectedSlotsTerm2To1":  None,       # transfer nothing back
          "Term1": {        # multidomain
            
              # electrophysiology solver
            "MultipleInstances": {
                "logKey":                     "duration_subdomains_xy",
                "ranksAllComputedInstances":  list(range(n_ranks)),
                "nInstances":                 variables.n_subdomains_xy,
                "instances": 
                [{
                  "ranks": list(range(subdomain_coordinate_y*variables.n_subdomains_x + subdomain_coordinate_x, n_ranks, variables.n_subdomains_x*variables.n_subdomains_y)),
                  
                  # this is for the actual model with fibers
                  "StrangSplitting": {
                    #"numberTimeSteps": 1,
                    "timeStepWidth":          variables.dt_splitting,  # 1e-1
                    "logTimeStepWidthAsKey":  "dt_splitting",
                    "durationLogKey":         "duration_monodomain",
                    "timeStepOutputInterval": 100,
                    "endTime":                variables.dt_splitting,
                    "connectedSlotsTerm1To2": None, #[0,1,2],   # transfer slot 0 = state Vm from Term1 (CellML) to Term2 (Diffusion)
                    "connectedSlotsTerm2To1": None, #[0,None,2],   # transfer the same back, this avoids data copy

                    "Term1": {      # CellML, i.e. reaction term of Monodomain equation
                      "MultipleInstances": {
                        "logKey":             "duration_subdomains_z",
                        "nInstances":         n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                        "instances": 
                        [{
                          "ranks":                          list(range(variables.n_subdomains_z)),   # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                          "Heun" : {
                            "timeStepWidth":                variables.dt_0D,                         # timestep width of 0D problem
                            "logTimeStepWidthAsKey":        "dt_0D",                                 # key under which the time step width will be written to the log file
                            "durationLogKey":               "duration_0D",                           # log key of duration for this solver
                            "timeStepOutputInterval":       1e4,                                     # how often to print the current timestep
                            "initialValues":                [],                                      # no initial values are specified
                            "dirichletBoundaryConditions":  {},                                      # no Dirichlet boundary conditions are specified
                            "dirichletOutputFilename":      None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                            
                            "inputMeshIsGlobal":            True,                                    # the boundary conditions and initial values would be given as global numbers
                            "checkForNanInf":               True,                                    # abort execution if the solution contains nan or inf values
                            "nAdditionalFieldVariables":    0,                                       # number of additional field variables
                            "additionalSlotNames":          [],                                      # names for the additional slots
                              
                            "CellML" : {
                              "modelFilename":                          variables.cellml_file,                          # input C++ source file or cellml XML file
                              #"statesInitialValues":                   [],                                             # if given, the initial values for the the states of one instance
                              "statesInitialValues":                    variables.states_initial_values,                # initial values for new_slow_TK
                              "initializeStatesToEquilibrium":          False,                                          # if the equilibrium values of the states should be computed before the simulation starts
                              "initializeStatesToEquilibriumTimestepWidth": 1e-4,                                       # if initializeStatesToEquilibrium is enable, the timestep width to use to solve the equilibrium equation
                              
                              # optimization parameters
                              "optimizationType":                       "vc" if variables.use_vc else "simd",           # "vc", "simd", "openmp" type of generated optimizated source file
                              "approximateExponentialFunction":         True,                                           # if optimizationType is "vc", whether the exponential function exp(x) should be approximate by (1+x/n)^n with n=1024
                              "compilerFlags":                          "-fPIC -O3 -march=native -Wno-deprecated-declarations -shared ",             # compiler flags used to compile the optimized model code
                              "maximumNumberOfThreads":                 0,                                              # if optimizationType is "openmp", the maximum number of threads to use. Default value 0 means no restriction.
                              
                              # stimulation callbacks
                              #"libraryFilename":                       "cellml_simd_lib.so",                           # compiled library
                              #"setSpecificParametersFunction":         set_specific_parameters,                        # callback function that sets parameters like stimulation current
                              #"setSpecificParametersCallInterval":     int(1./variables.stimulation_frequency/variables.dt_0D),         # set_specific_parameters should be called every 0.1, 5e-5 * 1e3 = 5e-2 = 0.05
                              "setSpecificStatesFunction":              set_specific_states,                                             # callback function that sets states like Vm, activation can be implemented by using this method and directly setting Vm values, or by using setParameters/setSpecificParameters
                              #"setSpecificStatesCallInterval":         2*int(1./variables.stimulation_frequency/variables.dt_0D),       # set_specific_states should be called variables.stimulation_frequency times per ms, the factor 2 is needed because every Heun step includes two calls to rhs
                              "setSpecificStatesCallInterval":          0,                                                               # 0 means disabled
                              "setSpecificStatesCallFrequency":         variables.get_specific_states_call_frequency(fiber_no, motor_unit_no),   # set_specific_states should be called variables.stimulation_frequency times per ms
                              "setSpecificStatesFrequencyJitter":       variables.get_specific_states_frequency_jitter(fiber_no, motor_unit_no), # random value to add or substract to setSpecificStatesCallFrequency every stimulation, this is to add random jitter to the frequency
                              "setSpecificStatesRepeatAfterFirstCall":  0.01,                                                            # [ms] simulation time span for which the setSpecificStates callback will be called after a call was triggered
                              "setSpecificStatesCallEnableBegin":       variables.get_specific_states_call_enable_begin(fiber_no, motor_unit_no),# [ms] first time when to call setSpecificStates
                              "additionalArgument":                     fiber_no,                                       # last argument that will be passed to the callback functions set_specific_states, set_specific_parameters, etc.
                              
                              # parameters to the cellml model
                              "mappings":                               variables.mappings,                             # mappings between parameters and algebraics/constants and between outputConnectorSlots and states, algebraics or parameters, they are defined in helper.py
                              "parametersInitialValues":                variables.parameters_initial_values,            #[0.0, 1.0],      # initial values for the parameters: I_Stim, l_hs
                              
                              "meshName":                               "MeshFiber_{}".format(fiber_no),                # reference to the fiber mesh
                              "stimulationLogFilename":                 "out/stimulation.log",                          # a file that will contain the times of stimulations
                            },      
                            "OutputWriter" : []
                            #[
                            #  {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.scenario_name + "/0D_states({},{})".format(fiber_in_subdomain_coordinate_x,fiber_in_subdomain_coordinate_y), "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"}
                            #] if variables.states_output else []
                            
                          },
                        } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                            for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                              for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                                for motor_unit_no in [get_motor_unit_no(fiber_no)]],
                      }
                    },
                    "Term2": {     # Diffusion
                      "MultipleInstances": {
                        "nInstances": n_fibers_in_subdomain_x(subdomain_coordinate_x)*n_fibers_in_subdomain_y(subdomain_coordinate_y),
                        "instances": 
                        [{
                          "ranks":                         list(range(variables.n_subdomains_z)),   # these rank nos are local nos to the outer instance of MultipleInstances, i.e. from 0 to number of ranks in z direction
                          "CrankNicolson" : {
                            "initialValues":               [],                                      # no initial values are given
                            #"numberTimeSteps":            1,
                            "timeStepWidth":               variables.dt_1D,                         # timestep width for the diffusion problem
                            "timeStepWidthRelativeTolerance": 1e-10,
                            "logTimeStepWidthAsKey":       "dt_1D",                                 # key under which the time step width will be written to the log file
                            "durationLogKey":              "duration_1D",                           # log key of duration for this solver
                            "timeStepOutputInterval":      1e4,                                     # how often to print the current timestep
                            "dirichletBoundaryConditions": {},                                      # old Dirichlet BC that are not used in FastMonodomainSolver: {0: -75.0036, -1: -75.0036},
                            "dirichletOutputFilename":     None,                                    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
                            "inputMeshIsGlobal":           True,                                    # initial values would be given as global numbers
                            "solverName":                  "diffusionTermSolver",                   # reference to the linear solver
                            "nAdditionalFieldVariables":   4,                                       # number of additional field variables that will be written to the output file, here for stress
                            "additionalSlotNames":         ["stress", "alpha", "lambda", "ldot"],   # names for the additional slots, maximum length is 10 characters per slot name
                            "checkForNanInf":              True,                                    # abort execution if the solution contains nan or inf values
                            
                            "FiniteElementMethod" : {
                              "inputMeshIsGlobal":         True,
                              "meshName":                  "MeshFiber_{}".format(fiber_no),
                              "solverName":                "diffusionTermSolver",
                              "prefactor":                 get_diffusion_prefactor(fiber_no, motor_unit_no),  # resolves to Conductivity / (Am * Cm)
                              "slotName":                  "vm",
                            },
                            "OutputWriter" : [
                              #{"format": "Paraview", "outputInterval": int(1./variables.dt_1D*variables.output_timestep), "filename": "out/fiber_"+str(fiber_no), "binary": True, "fixedFormat": False, "combineFiles": True},
                              #{"format": "Paraview", "outputInterval": 1./variables.dt_1D*variables.output_timestep, "filename": "out/fiber_"+str(i)+"_txt", "binary": False, "fixedFormat": False},
                              #{"format": "ExFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "sphereSize": "0.02*0.02*0.02"},
                              #{"format": "PythonFile", "filename": "out/fiber_"+str(i), "outputInterval": 1./variables.dt_1D*variables.output_timestep, "binary":True, "onlyNodalValues":True},
                            ]
                          },
                        } for fiber_in_subdomain_coordinate_y in range(n_fibers_in_subdomain_y(subdomain_coordinate_y)) \
                            for fiber_in_subdomain_coordinate_x in range(n_fibers_in_subdomain_x(subdomain_coordinate_x)) \
                              for fiber_no in [get_fiber_no(subdomain_coordinate_x, subdomain_coordinate_y, fiber_in_subdomain_coordinate_x, fiber_in_subdomain_coordinate_y)] \
                                for motor_unit_no in [get_motor_unit_no(fiber_no)]],
                        "OutputWriter" : variables.output_writer_fibers,
                      },
                    },
                  }
                  } if (subdomain_coordinate_x,subdomain_coordinate_y) == (variables.own_subdomain_coordinate_x,variables.own_subdomain_coordinate_y) else None
                  for subdomain_coordinate_y in range(variables.n_subdomains_y)
                    for subdomain_coordinate_x in range(variables.n_subdomains_x)]
          },
            "fiberDistributionFile":    variables.fiber_distribution_file,   # for FastMonodomainSolver, e.g. MU_fibre_distribution_3780.txt
            "firingTimesFile":          variables.firing_times_file,         # for FastMonodomainSolver, e.g. MU_firing_times_real.txt
            "onlyComputeIfHasBeenStimulated": variables.fast_monodomain_solver_optimizations,                          # only compute fibers after they have been stimulated for the first time
            "disableComputationWhenStatesAreCloseToEquilibrium": variables.fast_monodomain_solver_optimizations,       # optimization where states that are close to their equilibrium will not be computed again      
            "valueForStimulatedPoint":  variables.vm_value_stimulated,       # to which value of Vm the stimulated node should be set      
            "neuromuscularJunctionRelativeSize": 0.1,                        # range where the neuromuscular junction is located around the center, relative to fiber length. The actual position is draws randomly from the interval [0.5-s/2, 0.5+s/2) with s being this option. 0 means sharply at the center, 0.1 means located approximately at the center, but it can vary 10% in total between all fibers.
            "generateGPUSource":        True,                                # (set to True) only effective if optimizationType=="gpu", whether the source code for the GPU should be generated. If False, an existing source code file (which has to have the correct name) is used and compiled, i.e. the code generator is bypassed. This is useful for debugging, such that you can adjust the source code yourself. (You can also add "-g -save-temps " to compilerFlags under CellMLAdapter)
            "useSinglePrecision":       False,   
          },
          "Term2": {        # solid mechanics
          "MuscleContractionSolver": {
            "numberTimeSteps":              1,                         # only use 1 timestep per interval
            "timeStepOutputInterval":       100,                       # do not output time steps
            "Pmax":                         variables.pmax,            # maximum PK2 active stress
            "slotNames":                    [],
            "OutputWriter" : [
              {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/mechanics_3D", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles": True, "fileNumbering": "incremental"},
            ],
            "mapGeometryToMeshes":          ["3Dmesh"],    # the mesh names of the meshes that will get the geometry transferred
            "dynamic":                      True,                      # if the dynamic solid mechanics solver should be used, else it computes the quasi-static problem
            
            # the actual solid mechanics solver, this is either "DynamicHyperelasticitySolver" or "HyperelasticitySolver", depending on the value of "dynamic"
            "DynamicHyperelasticitySolver": {
              "timeStepWidth":              variables.dt_elasticity,           # time step width 
              "durationLogKey":             "duration_mechanics",               # key to find duration of this solver in the log file
              "timeStepOutputInterval":     1,                         # how often the current time step should be printed to console
              
              "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
              "density":                    variables.rho,             # density of the material
              "displacementsScalingFactor": 1.0,                       # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
              "residualNormLogFilename":    "log_residual_norm.txt",   # log file where residual norm values of the nonlinear solver will be written
              "useAnalyticJacobian":        True,                      # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
              "useNumericJacobian":         False,                     # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
                
              "dumpDenseMatlabVariables":   False,                     # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
              # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
              
              # mesh
              "inputMeshIsGlobal":          True,                     # the mesh is given locally
              "meshName":                   ["3Dmesh_elasticity_quadratic"],       # name of the 3D mesh, it is defined under "Meshes" at the beginning of this config
              "fiberMeshNames":             [],                       # fiber meshes that will be used to determine the fiber direction, there are no fibers in multidomain, so this is empty
              "fiberDirection":             [0,0,1],                  # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
        
              # solving
              "solverName":                 "mechanicsSolver",         # name of the nonlinear solver configuration, it is defined under "Solvers" at the beginning of this config
              #"loadFactors":                [0.5, 1.0],                # load factors for every timestep
              "loadFactors":                [],                        # no load factors, solve problem directly
              "loadFactorGiveUpThreshold":   1,                        # when to abort the solve
              "nNonlinearSolveCalls":       1,                         # how often the nonlinear solve should be repeated
              
              # boundary and initial conditions
              "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
              "neumannBoundaryConditions":   variables.elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
              "divideNeumannBoundaryConditionValuesByTotalArea": True,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
              "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
              "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
              
              "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
              "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
              "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
              "constantBodyForce":           variables.constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
              
              "dirichletOutputFilename":     "out/dirichlet_boundary_conditions",             # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
              
              # define which file formats should be written
              # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
              "OutputWriter" : [
                
                # Paraview files
                {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/mechanics_uv_stress", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                
                # Python callback function "postprocess"
                #{"format": "PythonCallback", "outputInterval": 1, "callback": postprocess, "onlyNodalValues":True, "filename": ""},
              ],
              # 2. additional output writer that writes also the hydrostatic pressure
              "pressure": {   # output files for pressure function space (linear elements), contains pressure values, as well as displacements and velocities
                "OutputWriter" : [
                  #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/p", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                ]
              },
              # 3. additional output writer that writes virtual work terms
              "dynamic": {    # output of the dynamic solver, has additional virtual work values 
                "OutputWriter" : [   # output files for displacements function space (quadratic elements)
                  {"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/mechanics_virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                  #{"format": "Paraview", "outputInterval": 1, "filename": "out/"+variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                ],
              },
              # 4. output writer for debugging, outputs files after each load increment, the geometry is not changed but u and v are written
              "LoadIncrements": {   
                "OutputWriter" : [
                  #{"format": "Paraview", "outputInterval": 1, "filename": "out/load_increments", "binary": False, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
                ]
              },
            }
                    }
                  }
                
              
            
          
        
      }
    
        }
      
    }
  }
}