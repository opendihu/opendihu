scenario_name = ""

# timing parameters [ms]
# ----------------------
dt_elasticity = 0.05
end_time = 10.0

# pulling force
# -------------
force = 2e4 
force_time_span = 10 # [ms] time at which full force is reached

output_timestep_elasticity = 1 # 1 means to output every timestep     

# material parameters
# -------------------
tendon_material = "nonLinear"
#tendon_material = "SaintVenantKirchoff"         
rho = 10   ## [1e-4 kg/cm^3] density of the water
constant_body_force = [0,0,0]

# single tendon case
# ------------------
tendon_extent = [2.0, 2.0, 6.0] # [cm, cm, cm] 
n_elements_tendon = [4,4,12] 
tendon_offset = [0.0, 0.0, 0.0]

# two tendon case
# ---------------
precice_file = "../precice_config.xml"

top_tendon_extent = [2.0, 2.0, 3.0] # [cm, cm, cm] 
top_n_elements_tendon = [4,4,6] 
top_tendon_offset = [0.0, 0.0, 3.0]

bottom_tendon_extent = [2.0, 2.0, 3.0] # [cm, cm, cm] 
bottom_n_elements_tendon = [4,4,6] 
bottom_tendon_offset = [0.0, 0.0, 0.0]

# boundary conditions
# -------------------
divideNeumannBoundaryConditionValuesByTotalArea = False
elasticity_dirichlet_bc = {}
elasticity_neumann_bc = []
meshes = {}

# mechanics solvers settings
# --------------------------
elasticity_solver_type = "preonly"  # type of the linear solver: cg groppcg pipecg pipecgrr cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres tsirm cgls
elasticity_preconditioner_type = "lu"

snes_max_iterations = 10                  # maximum number of iterations in the nonlinear solver
snes_rebuild_jacobian_frequency = 10       # how often the jacobian should be recomputed, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
snes_relative_tolerance = 1e-5            # relative tolerance of the nonlinear solver
snes_absolute_tolerance = 1e-5            # absolute tolerance of the nonlinear solver
snes_max_function_evaluations = 1000       # maximum number of function iterations
snes_linear_search = "l2"                 # type of linesearch, possible values: "bt" "nleqerr" "basic" "l2" "cp" "ncglinear"


linear_relative_tolerance = 1e-5           # relative tolerance of the residual of the linear solver
linear_absolute_tolerance = 1e-10          # absolute tolerance of the residual of the linear solver
linear_max_iterations = 1000                # maximum number of iterations in the linear solver

# partitioning
# ------------
# this has to match the total number of processes
n_subdomains_x = 1
n_subdomains_y = 1
n_subdomains_z = 1

