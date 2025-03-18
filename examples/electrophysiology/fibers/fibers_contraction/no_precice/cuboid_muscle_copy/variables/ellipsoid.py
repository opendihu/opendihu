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
end_time = 0.0         # end time of the simulation 
output_interval = dt_3D # time interval between outputs

# Material parameters
pmax = 7.3                                                  # maximum active stress
rho = 10                                                    # density of the muscle
material_parameters = [3.176e-10, 1.813, 1.075e-2, 1.0]     # [c1, c2, b, d]
diffusion_prefactor = 3.828 / (500.0 * 0.58)                # Conductivity / (Am * Cm)

# Meshes
el_x, el_y, el_z = 1, 1, 4                     # number of elements

fiber_direction = [0, 0, 1] # direction of fiber in element

mesh3D_nodes = []

a = 4
c = 12
xc = 2
yc = 2

for z in np.linspace(-5.5,7,9):
    # analytical solution
    r=np.sqrt(((1-z**2/c**2)/2*a**2))
    f = 1/np.sqrt(2)
    mesh3D_nodes.append([-r*f, -r*f, z])
    mesh3D_nodes.append([0.0, -r, z])
    mesh3D_nodes.append([r*f, -r*f, z])

    mesh3D_nodes.append([-r, 0.0, z])
    mesh3D_nodes.append([0.0, 0.0, z])
    mesh3D_nodes.append([r, 0.0, z])

    mesh3D_nodes.append([-r*f, r*f, z])
    mesh3D_nodes.append([0.0, r, z])
    mesh3D_nodes.append([r*f, r*f, z])


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

with open("../ellipsoid3d.json","r") as f:
	fdata = json.load(f)

fb_points = 100           # number of points per fiber

fiber_idx = 0
for fiber in fdata:
	fdict = fdata[fiber]
	npos = [[fdict[ii]['x'],fdict[ii]['y'],fdict[ii]['z']] for ii in range(len(fdict)-1) ]
	meshName = "fiber{}".format(fiber_idx)
	meshes[meshName] = {
			"nElements":		    [fb_points-1],
			"nodePositions":	    npos,
			"inputMeshIsGlobal":	True,
			"nRanks":				n_ranks
	}
	fiber_idx += 1
     
n_fibers = fiber_idx
n_fibers_x =5
n_fibers_y =6

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