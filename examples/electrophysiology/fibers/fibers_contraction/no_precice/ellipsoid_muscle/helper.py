# Multiple 1D fibers (monodomain) with 3D contraction, biceps geometry
# This is a helper script that sets a lot of the internal variables which are all defined in variables.py
#
# if variables.fiber_file=cuboid.bin, it uses a small cuboid test example

import numpy as np
import pickle
import sys,os
import struct
import argparse
sys.path.insert(0, '..')
import variables    # file variables.py
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])


# create the partitioning using the script in create_partitioned_meshes_for_settings.py
result = create_partitioned_meshes_for_settings(
    2,2,4, 
    variables.fiber_file, variables.load_fiber_data,
    variables.sampling_stride_x, variables.sampling_stride_y, variables.sampling_stride_z, variables.generate_linear_3d_mesh, variables.generate_quadratic_3d_mesh)
[variables.meshes, variables.own_subdomain_coordinate_x, variables.own_subdomain_coordinate_y, variables.own_subdomain_coordinate_z, variables.n_fibers_x, variables.n_fibers_y, variables.n_points_whole_fiber] = result
  

# create mappings between meshes
#variables.mappings_between_meshes = {"MeshFiber_{}".format(i) : "3Dmesh" for i in range(variables.n_fibers_total)}
variables.mappings_between_meshes = {"fiber_{}".format(i) : {"name": "3Dmesh", "xiTolerance": 3e-1, "defaultValue": 0} for i in range(variables.n_fibers)}

# a higher tolerance includes more fiber dofs that may be almost out of the 3D mesh
variables.mappings_between_meshes = {
  "fiber_{}".format(i) : {
    "name": "3Dmesh_quadratic",
    "xiTolerance": variables.mapping_tolerance,
    "enableWarnings": False, 
    "compositeUseOnlyInitializedMappings": False,
    "fixUnmappedDofs": True,
    "defaultValue": 0,
  } for i in range(variables.n_fibers)
}

