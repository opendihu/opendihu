# Transversely-isotropic Mooney Rivlin on a tendon geometry
# Note, this is not possible to be run in parallel because the fibers cannot be initialized without MultipleInstances class.
import sys, os
import numpy as np
import pickle
import argparse
import sys
sys.path.insert(0, '.')
import variables              # file variables.py, defines default values for all parameters, you can set the parameters there
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings with helper functions about own subdomain

# material parameters
# --------------------
# quantities in mechanics unit system
variables.rho = 10          # [1e-4 kg/cm^3] density of the muscle (density of water)

# material parameters for tendon material, see Carniel 2017 "A transversely isotropic coupled hyperelastic model for the mechanical behavior of tendons"
c = 9.98e-1                 # [N/cm^2 = 10kPa]    parameter for Fung model
ca = 14.92e-1               # [N/cm^2 = 10kPa]    normal stress axial
ct = 14.7e-1                # [N/cm^2 = 10kPa]    normal stress transversal
cat = 9.64e-1               # [N/cm^2 = 10kPa]    shear stress axial-transversal
ctt = 11.24e-1              # [N/cm^2 = 10kPa]    shear stress transversal-transversal
mu = 3.76e-1                # [N/cm^2 = 10kPa]    parameter for Neo-Hookean model
#k1 = 42.217e2               # [N/cm^2 = 1e-2 MPa]   shape parameter for fiber stress, human achilles (Couppe 2015)
#k2 = 411.360e2              # [N/cm^2 = 1e-2 MPa]   shape parameter for fiber stress, human achilles (Couppe 2015)
#k1 = 0.010e2                # [N/cm^2 = 1e-2 MPa]   shape parameter for fiber stress, human achilles (Csapo 2010) does not converge well
#k2 = 197.34e2               # [N/cm^2 = 1e-2 MPa]   shape parameter for fiber stress, human achilles (Csapo 2010)
#k1 = 2.893e2                # [N/cm^2 = 1e-2 MPa]   shape parameter for fiber stress, Equine(Horse) Digital Flexor (Thorpe 2012)
#k2 = 357.23e2               # [N/cm^2 = 1e-2 MPa]   shape parameter for fiber stress, Equine Digital Flexor (Thorpe 2012)
k1 = 92.779e2               # [N/cm^2 = 1e-2 MPa]   shape parameter for fiber stress, Equine Digital Flexor (Vergari 2011)
k2 = 305.87e2               # [N/cm^2 = 1e-2 MPa]   shape parameter for fiber stress, Equine Digital Flexor (Vergari 2011)
variables.material_parameters = [c, ca, ct, cat, ctt, mu, k1, k2]

variables.constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
variables.force = 1000.0       # [N]

print("expected stress: {} N/cm^2 = {} MPa".format(variables.force, variables.force * 1e-2))

variables.dt_elasticity = 0.01      # [ms] time step width for elasticity
variables.end_time      = 10     # [ms] simulation time
variables.scenario_name = "tendon_bottom"
variables.is_bottom_tendon = True        # whether the tendon is at the bottom (negative z-direction), this is important for the boundary conditions

# input mesh file
variables.fiber_file = "../../../../input/left_biceps_brachii_tendon1.bin"        # bottom tendon
#variables.fiber_file = "../../../../input/left_biceps_brachii_tendon2a.bin"
#variables.fiber_file = "../../../../input/left_biceps_brachii_tendon2b.bin"
#variables.fiber_file = "../../../../input/left_biceps_brachii_7x7fibers.bin"
#variables.fiber_file = "../../../../input/left_biceps_brachii_7x7fibers.bin"

load_fiber_data = False             # If the fiber geometry data should be loaded completely in the python script. If True, this reads the binary file and assigns the node positions in the config. If False, the C++ code will read the binary file and only extract the local node positions. This is more performant for highly parallel runs.

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

# define command line arguments
parser = argparse.ArgumentParser(description='tendon')
parser.add_argument('--precice_config_file',                       help='The name to identify this run in the log.',   default=variables.precice_config_file)
parser.add_argument('--case_name',                       help='The name to identify this run in the log.',   default=variables.case_name)
parser.add_argument('--n_subdomains', nargs=3,               help='Number of subdomains in x,y,z direction.',    type=int)
parser.add_argument('--n_subdomains_x', '-x',                help='Number of subdomains in x direction.',        type=int, default=variables.n_subdomains_x)
parser.add_argument('--n_subdomains_y', '-y',                help='Number of subdomains in y direction.',        type=int, default=variables.n_subdomains_y)
parser.add_argument('--n_subdomains_z', '-z',                help='Number of subdomains in z direction.',        type=int, default=variables.n_subdomains_z)
parser.add_argument('--fiber_file',                          help='The filename of the file that contains the fiber data.', default=variables.fiber_file)
parser.add_argument('-vmodule', help='ignore')

# parse command line arguments and assign values to variables module
args = parser.parse_known_args(args=sys.argv[:-2], namespace=variables)

# partitioning
# ------------
# this has to match the total number of processes
if variables.n_subdomains is not None:
  variables.n_subdomains_x = variables.n_subdomains[0]
  variables.n_subdomains_y = variables.n_subdomains[1]
  variables.n_subdomains_z = variables.n_subdomains[2]

# compute partitioning
if rank_no == 0:
  if n_ranks != variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z:
    print("\n\nError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
    sys.exit(-1)
    
# stride for sampling the 3D elements from the fiber data
# here any number is possible
sampling_stride_x = 1
sampling_stride_y = 1
sampling_stride_z = 2

# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, 
    variables.fiber_file, load_fiber_data,
    sampling_stride_x, sampling_stride_y, sampling_stride_z, True, True)

[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result

n_points_3D_mesh_linear_global_x = sum([n_sampled_points_in_subdomain_x(subdomain_coordinate_x) for subdomain_coordinate_x in range(variables.n_subdomains_x)])
n_points_3D_mesh_linear_global_y = sum([n_sampled_points_in_subdomain_y(subdomain_coordinate_y) for subdomain_coordinate_y in range(variables.n_subdomains_y)])
n_points_3D_mesh_linear_global_z = sum([n_sampled_points_in_subdomain_z(subdomain_coordinate_z) for subdomain_coordinate_z in range(variables.n_subdomains_z)])
n_points_3D_mesh_linear_global = n_points_3D_mesh_linear_global_x*n_points_3D_mesh_linear_global_y*n_points_3D_mesh_linear_global_z
nx = n_points_3D_mesh_linear_global_x-1
ny = n_points_3D_mesh_linear_global_y-1
nz = n_points_3D_mesh_linear_global_z-1

node_positions = variables.meshes["3Dmesh_quadratic"]["nodePositions"]

# boundary conditions (for quadratic elements)
# --------------------------------------------
[mx, my, mz] = variables.meshes["3Dmesh_quadratic"]["nPointsGlobal"]
[nx, ny, nz] = variables.meshes["3Dmesh_quadratic"]["nElements"]

# set Dirichlet BC, fix one end
variables.elasticity_dirichlet_bc = {}
k = mz-1

# fix z value on the whole x-y-plane
for j in range(my):
  for i in range(mx):
    variables.elasticity_dirichlet_bc[k*mx*my + j*mx + i] = [None,None,0.0,None,None,None]

# fix left edge 
for j in range(my):
  variables.elasticity_dirichlet_bc[k*mx*my + j*mx + 0][0] = 0.0
  
# fix front edge 
for i in range(mx):
  variables.elasticity_dirichlet_bc[k*mx*my + 0*mx + i][1] = 0.0
       
# set Neumann BC, set traction at the end of the tendon that is attached to the muscle
k = 0
traction_vector = [0, 0, -variables.force]     # the traction force in specified in the reference configuration
#traction_vector = [0, 0.1*variables.force, -0.2*variables.force]     # the traction force in specified in the reference configuration
face = "2-"
variables.elasticity_neumann_bc = [{"element": k*nx*ny + j*nx + i, "constantVector": traction_vector, "face": face} for j in range(ny) for i in range(nx)]
variables.elasticity_dirichlet_bc = {}

print("nRanks: ",variables.meshes["3Dmesh_quadratic"]["nRanks"])


config_hyperelasticity = {    # for both "HyperelasticitySolver" and "DynamicHyperelasticitySolver"
  "timeStepWidth":              variables.dt_elasticity,      # time step width 
  "endTime":                    variables.end_time,           # end time of the simulation time span    
  "durationLogKey":             "duration_mechanics",         # key to find duration of this solver in the log file
  "timeStepOutputInterval":     1,                            # how often the current time step should be printed to console
  
  "materialParameters":         variables.material_parameters,  # material parameters of the Mooney-Rivlin material
  "density":                    variables.rho,                # density of the material
  "displacementsScalingFactor": 1.0,                          # scaling factor for displacements, only set to sth. other than 1 only to increase visual appearance for very small displacements
  "residualNormLogFilename":    "log_residual_norm.txt",      # log file where residual norm values of the nonlinear solver will be written
  "useAnalyticJacobian":        True,                         # whether to use the analytically computed jacobian matrix in the nonlinear solver (fast)
  "useNumericJacobian":         False,                        # whether to use the numerically computed jacobian matrix in the nonlinear solver (slow), only works with non-nested matrices, if both numeric and analytic are enable, it uses the analytic for the preconditioner and the numeric as normal jacobian
    
  "dumpDenseMatlabVariables":   False,                        # whether to have extra output of matlab vectors, x,r, jacobian matrix (very slow)
  # if useAnalyticJacobian,useNumericJacobian and dumpDenseMatlabVariables all all three true, the analytic and numeric jacobian matrices will get compared to see if there are programming errors for the analytic jacobian
  
  # mesh
  "meshName":                   "3Dmesh_quadratic",           # mesh with quadratic Lagrange ansatz functions
  "inputMeshIsGlobal":          True,                         # boundary conditions are specified in global numberings, whereas the mesh is given in local numberings
  
  "fiberMeshNames":             [],                           # fiber meshes that will be used to determine the fiber direction
  "fiberDirection":             [0,0,1],                      # if fiberMeshNames is empty, directly set the constant fiber direction, in element coordinate system
  
  # nonlinear solver
  "relativeTolerance":          1e-10,                         # 1e-10 relative tolerance of the linear solver
  "absoluteTolerance":          1e-10,                        # 1e-10 absolute tolerance of the residual of the linear solver       
  "solverType":                 "preonly",                    # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
  "preconditionerType":         "lu",                         # type of the preconditioner
  "maxIterations":              1e4,                          # maximum number of iterations in the linear solver
  "snesMaxFunctionEvaluations": 1e8,                          # maximum number of function iterations
  "snesMaxIterations":          240,                           # maximum number of iterations in the nonlinear solver
  "snesRelativeTolerance":      1e-2,                         # relative tolerance of the nonlinear solver
  "snesLineSearchType":         "l2",                         # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"
  "snesAbsoluteTolerance":      1e-5,                         # absolute tolerance of the nonlinear solver
  "snesRebuildJacobianFrequency": 1,                          # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
  
  #"dumpFilename": "out/r{}/m".format(sys.argv[-1]),          # dump system matrix and right hand side after every solve
  "dumpFilename":               "",                           # dump disabled
  "dumpFormat":                 "matlab",                     # default, ascii, matlab
  
  #"loadFactors":                [0.1, 0.2, 0.35, 0.5, 1.0],   # load factors for every timestep
  #"loadFactors":                [0.5, 1.0],                   # load factors for every timestep
  "loadFactors":                [],                           # no load factors, solve problem directly
  "loadFactorGiveUpThreshold":  1e-3,                         # a threshold for the load factor, when to abort the solve of the current time step. The load factors are adjusted automatically if the nonlinear solver diverged. If the load factors get too small, it aborts the solve.
  "nNonlinearSolveCalls":       1,                            # how often the nonlinear solve should be called
  
  # boundary and initial conditions
  "dirichletBoundaryConditions": variables.elasticity_dirichlet_bc,   # the initial Dirichlet boundary conditions that define values for displacements u and velocity v
  "neumannBoundaryConditions":   variables.elasticity_neumann_bc,     # Neumann boundary conditions that define traction forces on surfaces of elements
  "divideNeumannBoundaryConditionValuesByTotalArea": False,            # if the given Neumann boundary condition values under "neumannBoundaryConditions" are total forces instead of surface loads and therefore should be scaled by the surface area of all elements where Neumann BC are applied
  "updateDirichletBoundaryConditionsFunction": None,                  # function that updates the dirichlet BCs while the simulation is running
  "updateDirichletBoundaryConditionsFunctionCallInterval": 1,         # every which step the update function should be called, 1 means every time step
  
  "initialValuesDisplacements":  [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the displacements, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  "initialValuesVelocities":     [[0.0,0.0,0.0] for _ in range(mx*my*mz)],     # the initial values for the velocities, vector of values for every node [[node1-x,y,z], [node2-x,y,z], ...]
  "extrapolateInitialGuess":     True,                                # if the initial values for the dynamic nonlinear problem should be computed by extrapolating the previous displacements and velocities
  "constantBodyForce":           variables.constant_body_force,       # a constant force that acts on the whole body, e.g. for gravity
  
  "dirichletOutputFilename":     "out/" + variables.case_name + '/' +variables.scenario_name+"/dirichlet_boundary_conditions_tendon",    # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
  
  # define which file formats should be written
  # 1. main output writer that writes output files using the quadratic elements function space. Writes displacements, velocities and PK2 stresses.
  "OutputWriter" : [
    
    # Paraview files
    {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.case_name + '/' +variables.scenario_name+"/tendon", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
    
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
      {"format": "Paraview", "outputInterval": 1, "filename": "out/" + variables.case_name + '/' +variables.scenario_name+"/virtual_work", "binary": True, "fixedFormat": False, "onlyNodalValues":True, "combineFiles":True, "fileNumbering": "incremental"},
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

config = {
  "scenarioName":                   variables.scenario_name,      # scenario name to identify the simulation runs in the log file
  "logFormat":                      "csv",                        # "csv" or "json", format of the lines in the log file, csv gives smaller files
  "solverStructureDiagramFile":     "solver_structure.txt",       # output file of a diagram that shows data connection between solvers
  "mappingsBetweenMeshesLogFile":   "mappings_between_meshes_log.txt",    # log file for mappings 
  "Meshes":                         variables.meshes,
  
  "PreciceAdapter": {        # precice adapter for bottom tendon
    "timeStepOutputInterval":   100,                        # interval in which to display current timestep and time in console
    "timestepWidth":            1,                          # coupling time step width, must match the value in the precice config
    "couplingEnabled":          True,                       # if the precice coupling is enabled, if not, it simply calls the nested solver, for debugging
    "preciceConfigFilename":    variables.precice_config_file,    # the preCICE configuration file
    "preciceParticipantName":   "TendonSolver",             # name of the own precice participant, has to match the name given in the precice xml config file
    "preciceMeshes": [                                      # the precice meshes get created as the top or bottom surface of the main geometry mesh of the nested solver
      {
        "preciceMeshName":      "TendonMeshTop",            # precice name of the 2D coupling mesh
        "face":                 "2+",                       # face of the 3D mesh where the 2D mesh is located, "2-" = bottom, "2+" = top
      }
    ],
    "preciceData": [
      {
        "mode":                 "write-displacements-velocities",   # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "TendonMeshTop",                    # name of the precice coupling surface mesh, as given in the precice xml settings file
        "displacementsName":    "Displacement",                     # name of the displacements "data", i.e. field variable, as given in the precice xml settings file
        "velocitiesName":       "Velocity",                         # name of the velocity "data", i.e. field variable, as given in the precice xml settings file
      },
      {
        "mode":                 "read-traction",                    # mode is one of "read-displacements-velocities", "read-traction", "write-displacements-velocities", "write-traction"
        "preciceMeshName":      "TendonMeshTop",                    # name of the precice coupling surface mesh, as given in the precice xml settings 
        "tractionName":         "Traction",                         # name of the traction "data", i.e. field variable, as given in the precice xml settings file
      }
    ],
    "HyperelasticitySolver": config_hyperelasticity,
    "DynamicHyperelasticitySolver": config_hyperelasticity,
  }
}
