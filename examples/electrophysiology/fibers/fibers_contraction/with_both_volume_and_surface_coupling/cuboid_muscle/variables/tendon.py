import sys

rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

scenario_name = "tendon"

# timing parameters [ms]
# ----------------------
dt_3D = 0.05
end_time = 10.0

output_interval = dt_3D # time interval between outputs

# material parameters
# -------------------
c = 9.98                    # [N/cm^2=kPa]
ca = 14.92                  # [-]
ct = 14.7                   # [-]
cat = 9.64                  # [-]
ctt = 11.24                 # [-]
mu = 3.76                   # [N/cm^2=kPa]
k1 = 42.217e3               # [N/cm^2=kPa]
k2 = 411.360e3              # [N/cm^2=kPa]

material_parameters = [c, ca, ct, cat, ctt, mu, k1, k2]
rho = 10   ## [1e-4 kg/cm^3] density of the water
constant_body_force = [0,0,0]

# Meshes
ex_x, ex_y, ex_z = 3.0, 3.0, 4.0               # extent of muscle
el_x, el_y, el_z = 3, 3, 4                     # number of elements
bs_x, bs_y, bs_z = 2*el_x+1, 2*el_y+1, 2*el_z+1 # quadratic basis functions

meshes = { # create 3D mechanics mesh
    "tendonMesh3D": {
        "nElements":            [el_x, el_y, el_z],
        "physicalExtent":       [ex_x, ex_y, ex_z],
        "physicalOffset":       [0, 0, -4.0],
        "logKey":               "mesh3D",
        "inputMeshIsGlobal":    True,
        "nRanks":               n_ranks
    }
}
# boundary conditions
# -------------------
# Boundary conditions
dirichlet_bc = {} # fix z=0 with dirichlet boundary conditions
for x in range(bs_x):
    for y in range(bs_y):
        dirichlet_bc[x + y*bs_x] = [0.0, 0.0, 0.0, None, None, None]

neumann_bc = []

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

