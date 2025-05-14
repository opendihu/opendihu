import sys
import itertools
import numpy as np
import json

rank_no = int(sys.argv[-2])
n_ranks = int(sys.argv[-1])

# PreCICE
precice_config = "../../variables/precice-config.xml"

# Time stepping
dt_3D = 1e-1            # time step of 3D mechanics
dt_splitting = 2e-3     # time step of strang splitting
dt_1D = 2e-3            # time step of 1D fiber diffusion
dt_0D = 1e-3            # time step of 0D cellml problem
end_time = 20         # end time of the simulation 
output_interval = dt_3D # time interval between outputs

# Material parameters
pmax = 7.3                                                  # maximum active stress
rho = 10                                                    # density of the muscle
material_parameters = [3.176e-10, 1.813, 1.075e-2, 1.0]     # [c1, c2, b, d]
diffusion_prefactor = 3.828 / (500.0 * 0.58)                # Conductivity / (Am * Cm)

# Meshes
ex_x, ex_y, ex_z = 6.0, 6.0, 10.0               # extent of muscle
el_x, el_y, el_z = 2, 2, 4                     # number of elements
# bs_x, bs_y, bs_z = 2*el_x+1, 2*el_y+1, 2*el_z+1 # quadratic basis functions

fiber_direction = [0, 0, 1] # direction of fiber in element

mesh3D_nodes = []

ox= 0.6
oxx = 0.15
oxxx = 0.05

for z in np.linspace(0, 12, 9):
    # for y in np.linspace(0, 4, 5):
    #     for x in np.linspace(0, 4, 5):
    #         mesh3D_nodes.append([x,y,z])
    mesh3D_nodes.append([0.0 + ox, 0.0+ox, z])
    mesh3D_nodes.append([1.0+oxx, 0.0+oxx, z])
    mesh3D_nodes.append([2.0, 0.0, z])
    mesh3D_nodes.append([3.0-oxx, 0.0+oxx, z])
    mesh3D_nodes.append([4.0-ox, 0.0+ox, z])

    mesh3D_nodes.append([0.0 + oxx, 1.0 + oxx, z])
    mesh3D_nodes.append([1.0+oxxx, 1.0+oxxx, z])
    mesh3D_nodes.append([2.0, 1.0, z])
    mesh3D_nodes.append([3.0-oxxx, 1.0+oxxx, z])
    mesh3D_nodes.append([4.0-oxx, 1.0+oxx, z])

    mesh3D_nodes.append([0.0, 2.0, z])
    mesh3D_nodes.append([1.0, 2.0, z])
    mesh3D_nodes.append([2.0, 2.0, z])
    mesh3D_nodes.append([3.0, 2.0, z])
    mesh3D_nodes.append([4.0, 2.0, z])

    mesh3D_nodes.append([0.0 + oxx, 3.0-oxx, z])
    mesh3D_nodes.append([1.0+oxx, 3.0-oxxx, z])
    mesh3D_nodes.append([2.0, 3.0, z])
    mesh3D_nodes.append([3.0-oxxx, 3.0-oxxx, z])
    mesh3D_nodes.append([4.0-oxx, 3.0-oxx, z])

    mesh3D_nodes.append([0.0+ox, 4.0-ox, z])
    mesh3D_nodes.append([1.0+oxx, 4.0-oxx, z])
    mesh3D_nodes.append([2.0, 4.0, z])
    mesh3D_nodes.append([3.0-oxx, 4.0-oxx, z])
    mesh3D_nodes.append([4.0-ox, 4.0-ox, z])
# for z in [0,1.25,2.5,3.75,5,6.25,7.5,8.75,10]:
    # for y in [1.5,3.0,4.5,6.0,7.5]:
    #      for x in [1.5,3.0,4.5,6.0,7.5]:
    #           mesh3D_nodes.append([x,y,z])
    # mesh3D_nodes.append([2.3,2.3,z])
    # mesh3D_nodes.append([4.5,1.5,z])
    # mesh3D_nodes.append([6.7,2.3,z])
    # mesh3D_nodes.append([1.5,4.5,z])
    # mesh3D_nodes.append([4.5,4.5,z])
    # mesh3D_nodes.append([7.5,4.5,z])
    # mesh3D_nodes.append([2.3,6.7,z])
    # mesh3D_nodes.append([4.5,7.5,z])
    # mesh3D_nodes.append([6.7,6.7,z])
    # ox = 0.8
    # oxx = 0.2
    # oxxx = 0.15

    # mesh3D_nodes.append([1.5 + ox,1.5 + ox,z])
    # mesh3D_nodes.append([3.0 + oxx,1.5+oxx,z])
    # mesh3D_nodes.append([4.5,1.5,z])
    # mesh3D_nodes.append([6.0-oxx,1.5+oxx,z])
    # mesh3D_nodes.append([7.5-ox,1.5+ox,z])

    # mesh3D_nodes.append([1.5+oxx,3.0+oxx,z])
    # mesh3D_nodes.append([3.0+oxxx,3.0+oxxx,z])
    # mesh3D_nodes.append([4.5,3.0,z])
    # mesh3D_nodes.append([6.0-oxxx,3.0+oxxx,z])
    # mesh3D_nodes.append([7.5-oxx,3.0+oxx,z])

    # mesh3D_nodes.append([1.5,4.5,z])
    # mesh3D_nodes.append([3.0,4.5,z])
    # mesh3D_nodes.append([4.5,4.5,z])
    # mesh3D_nodes.append([6.0,4.5,z])
    # mesh3D_nodes.append([7.5,4.5,z])

    # mesh3D_nodes.append([1.5+oxx,6.0-oxx,z])
    # mesh3D_nodes.append([3.0+oxxx,6.0-oxxx,z])
    # mesh3D_nodes.append([4.5,6.0,z])
    # mesh3D_nodes.append([6.0-oxxx,6.0-oxxx,z])
    # mesh3D_nodes.append([7.5-oxx,6.0-oxx,z])

    # mesh3D_nodes.append([1.5+ox,7.5-ox,z])
    # mesh3D_nodes.append([3.0+oxx,7.5-oxx,z])
    # mesh3D_nodes.append([4.5,7.5,z])
    # mesh3D_nodes.append([6.0-oxx,7.5-oxx,z])
    # mesh3D_nodes.append([7.5-ox,7.5-ox,z])

print(mesh3D_nodes)

meshes = { # create 3D mechanics mesh
    "mesh3D": {
        "nElements":            [el_x, el_y, el_z],
        # "nodeDimension": 3,
        "nodePositions":        mesh3D_nodes,
        "logKey":               "mesh3D",
        "inputMeshIsGlobal":    True,
        "nRanks":               n_ranks
    }
}



# Boundary conditions
dirichlet_bc = {} # fix z=0 with dirichlet boundary conditions
for x in range(el_x+1):
    for y in range(el_y+1):
        dirichlet_bc[x + y*(el_x+1)] = [0.0, 0.0, 0.0, None, None, None]

neumann_bc = [] # add pulling force to z=el_z with neumann boundary conditions
neumann_force = 0
for x in range(el_x):
    for y in range(el_y):
        neumann_bc += [{
            "element": x + y*el_x + (el_z-1)*el_y*el_x, 
            "constantVector": [0, 0, neumann_force], 
            "face": "2+"
        }]

# Define directory for cellml files
import os
input_dir = os.path.join(os.environ.get('OPENDIHU_HOME', '../../../../../../../'), "examples/electrophysiology/input/")

# Fiber activation
fiber_distribution_file = input_dir + "MU_fibre_distribution_3780.txt"
firing_times_file = input_dir + "MU_firing_times_always.txt"
specific_states_call_enable_begin = 1.0                     # time of first fiber activation
specific_states_call_frequency = 1e-3                       # frequency of fiber activation