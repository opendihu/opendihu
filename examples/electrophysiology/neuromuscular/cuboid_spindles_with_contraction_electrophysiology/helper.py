# Multiple 1D fibers (monodomain) with 3D contraction, biceps geometry
# This is a helper script that sets a lot of the internal variables which are all defined in variables.py

import numpy as np
import scipy
import pickle
import sys,os
import struct
import argparse
import random
import time
sys.path.insert(0, '..')
import variables    # file variables.py
from create_partitioned_meshes_for_settings import *   # file create_partitioned_meshes_for_settings

# parse arguments
rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

variables.n_subdomains = variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z

if variables.n_subdomains != n_ranks:
  print("\n\n\033[0;31mError! Number of ranks {} does not match given partitioning {} x {} x {} = {}.\033[0m\n\n".format(n_ranks, variables.n_subdomains_x, variables.n_subdomains_y, variables.n_subdomains_z, variables.n_subdomains_x*variables.n_subdomains_y*variables.n_subdomains_z))
  quit()

#variables.relative_factors_file = "compartments_relative_factors.{}.{}_mus_partitioning_{}x{}x{}".\
#  format(os.path.basename(variables.fiber_file),len(variables.motor_units),variables.n_subdomains_x,variables.n_subdomains_y,variables.n_subdomains_z)
#
include_global_node_positions = False
#if not os.path.exists(variables.relative_factors_file) and rank_no == 0:
#  include_global_node_positions = True

  
#############################

# set output writer    
variables.output_writer_fibers = []
variables.output_writer_elasticity = []
variables.output_writer_emg = []
variables.output_writer_0D_states = []

if True:  
  subfolder = ""
  if variables.paraview_output:
    if variables.adios_output:
      subfolder = "paraview/"
    variables.output_writer_elasticity.append({"format": "Paraview", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "binary": True, "fixedFormat": False, "combineFiles": True})
    variables.output_writer_fibers.append({"format": "Paraview", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers/100), "filename": "out/fibers", "binary": True, "fixedFormat": False, "combineFiles": True, "fileNumbering": "incremental"})
    if variables.states_output:
      variables.output_writer_0D_states.append({"format": "Paraview", "outputInterval": 1, "filename": "out/" + subfolder + variables.scenario_name + "/0D_states", "binary": True, "fixedFormat": False, "combineFiles": True})

  if variables.adios_output:
    if variables.paraview_output:
      subfolder = "adios/"
    variables.output_writer_elasticity.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "useFrontBackBuffer": False})
    variables.output_writer_fibers.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "combineNInstances": variables.n_subdomains_xy, "useFrontBackBuffer": False})
    #variables.output_writer_fibers.append({"format": "MegaMol", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + variables.scenario_name + "/fibers", "combineNInstances": 1, "useFrontBackBuffer": False}

  if variables.python_output:
    if variables.adios_output:
      subfolder = "python/"
    variables.output_writer_elasticity.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity", "binary": True})
    variables.output_writer_fibers.append({"format": "PythonFile", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers", "binary": True})

  if variables.exfile_output:
    if variables.adios_output:
      subfolder = "exfile/"
    variables.output_writer_elasticity.append({"format": "Exfile", "outputInterval": int(1./variables.dt_elasticity*variables.output_timestep_elasticity), "filename": "out/" + subfolder + variables.scenario_name + "/elasticity"})
    variables.output_writer_fibers.append({"format": "Exfile", "outputInterval": int(1./variables.dt_splitting*variables.output_timestep_fibers), "filename": "out/" + subfolder + variables.scenario_name + "/fibers"})

# set variable mappings for cellml model
if "hodgkin_huxley" in variables.cellml_file and "hodgkin_huxley-razumova" not in variables.cellml_file:
  # parameters: I_stim
  variables.mappings = {
    ("parameter", 0):           ("constant", "membrane/i_Stim"),      # parameter 0 is constant 2 = I_stim
    ("connectorSlot", 0): ("state", "membrane/V"),              # expose state 0 = Vm to the operator splitting
  }
  variables.parameters_initial_values = [0.0]                         # initial value for stimulation current
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 20.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)

elif "shorten" in variables.cellml_file:
  # parameters: stimulation current I_stim, fiber stretch λ
  variables.mappings = {
    ("parameter", 0):           ("algebraic", "wal_environment/I_HH"), # parameter is algebraic 32
    ("parameter", 1):           ("constant", "razumova/L_x"),             # parameter is constant 65, fiber stretch λ, this indicates how much the fiber has stretched, 1 means no extension
    ("connectorSlot", 0): ("state", "wal_environment/vS"),          # expose state 0 = Vm to the operator splitting
  }
  variables.parameters_initial_values = [0.0, 1.0]                        # stimulation current I_stim, fiber stretch λ
  variables.nodal_stimulation_current = 1200.                             # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "slow_TK_2014" in variables.cellml_file:   # this is (3a, "MultiPhysStrain", old tomo mechanics) in OpenCMISS
  # parameters: I_stim, fiber stretch λ
  variables.mappings = {
    ("parameter", 0):           ("constant", "wal_environment/I_HH"), # parameter 0 is constant 54 = I_stim
    ("parameter", 1):           ("constant", "razumova/L_S"),         # parameter 1 is constant 67 = fiber stretch λ
    ("connectorSlot", 0): ("state", "wal_environment/vS"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): ("algebraic", "razumova/stress"),  # expose algebraic 12 = γ to the operator splitting
  }
  variables.parameters_initial_values = [0.0, 1.0]                    # wal_environment/I_HH = I_stim, razumova/L_S = λ
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "Aliev_Panfilov_Razumova_2016_08_22" in variables.cellml_file :   # this is (3, "MultiPhysStrain", numerically more stable) in OpenCMISS, this only computes A1,A2,x1,x2 not the stress
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.mappings = {
    ("parameter", 0):           ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
    ("parameter", 1):           ("constant", "Razumova/l_hs"),        # parameter 1 is constant 8 = fiber stretch λ
    ("parameter", 2):           ("constant", "Razumova/velo"),        # parameter 2 is constant 9 = fiber contraction velocity \dot{λ}
    ("connectorSlot", 0): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): ("algebraic", "Razumova/sigma"),   # expose algebraic 0 = γ to the operator splitting
  }
  variables.parameters_initial_values = [0, 1, 0]                     # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/velo = \dot{λ}
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "Aliev_Panfilov_Razumova_Titin" in variables.cellml_file:   # this is (4, "Titin") in OpenCMISS
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.mappings = {
    ("parameter", 0):           ("constant", "Aliev_Panfilov/I_HH"),  # parameter 0 is constant 0 = I_stim
    ("parameter", 1):           ("constant", "Razumova/l_hs"),        # parameter 1 is constant 11 = fiber stretch λ
    ("parameter", 2):           ("constant", "Razumova/rel_velo"),    # parameter 2 is constant 12 = fiber contraction velocity \dot{λ}
    ("connectorSlot", 0): ("state", "Aliev_Panfilov/V_m"),      # expose state 0 = Vm to the operator splitting
    ("connectorSlot", 1): ("algebraic", "Razumova/ActiveStress"),   # expose algebraic 4 = γ to the operator splitting
    ("connectorSlot", 2): ("algebraic", "Razumova/Activation"),     # expose algebraic 5 = α to the operator splitting
  }
  variables.parameters_initial_values = [0, 1, 0]                     # Aliev_Panfilov/I_HH = I_stim, Razumova/l_hs = λ, Razumova/rel_velo = \dot{λ}
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 40.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)
  
elif "hodgkin_huxley-razumova" in variables.fiber_cellml_file:   # this is (4, "Titin") in OpenCMISS
  # parameters: I_stim, fiber stretch λ, fiber contraction velocity \dot{λ}
  variables.fiber_mappings = {
    ("parameter", 0):           "membrane/i_Stim",          # parameter 0 is I_stim
    ("parameter", 1):           "Razumova/l_hs",            # parameter 1 is fiber stretch λ
    ("connectorSlot", 0): "membrane/V",               # expose Vm to the operator splitting
    ("connectorSlot", 2): "Razumova/activation",      # expose activation .
    ("connectorSlot", 1): "Razumova/activestress",
  }
  variables.parameters_initial_values = [0, 1]
  variables.nodal_stimulation_current = 40.                           # not used
  variables.vm_value_stimulated = 20.                                 # to which value of Vm the stimulated node should be set (option "valueForStimulatedPoint" of FastMonodomainSolver)

else:
  pass
  #print("\033[0;31mCellML file {} has no mappings implemented in helper.py\033[0m".format(variables.cellml_file))
  #quit()

# load MU distribution and firing times
variables.firing_times = np.genfromtxt(variables.firing_times_file)

# ---------------------
# callback functions
def get_motor_unit_no(compartment_no):
  return compartment_no

def get_diffusion_prefactor(fiber_no, mu_no):
  diffusion_prefactor = variables.get_conductivity(fiber_no, mu_no) / (variables.get_am(fiber_no, mu_no) * variables.get_cm(fiber_no, mu_no))
  #print("diffusion_prefactor: {}/({}*{}) = {}".format(variables.get_conductivity(fiber_no, mu_no), variables.get_am(fiber_no, mu_no), variables.get_cm(fiber_no, mu_no), diffusion_prefactor))
  return diffusion_prefactor

def compartment_gets_stimulated(compartment_no, frequency, current_time):
  """
  determine if compartment compartment_no gets stimulated at simulation time current_time
  """

  # determine motor unit
  alpha = 1.0   # 0.8
  mu_no = (int)(get_motor_unit_no(compartment_no)*alpha)
  
  # determine if fiber fires now
  index = int(np.round(current_time * frequency))
  n_firing_times = np.size(variables.firing_times,0)
  
  #if variables.firing_times[index % n_firing_times, mu_no] == 1:
    #print("{}: fiber {} is mu {}, t = {}, row: {}, stimulated: {} {}".format(rank_no, fiber_no, mu_no, current_time, (index % n_firing_times), variables.firing_times[index % n_firing_times, mu_no], "true" if variables.firing_times[index % n_firing_times, mu_no] == 1 else "false"))
  
  return variables.firing_times[index % n_firing_times, mu_no] == 1
  
# callback function that can set states, i.e. prescribed values for stimulation
def set_specific_states(n_nodes_global, time_step_no, current_time, states, compartment_no):
  
  # determine if fiber gets stimulated at the current time
  is_compartment_gets_stimulated = compartment_gets_stimulated(compartment_no, variables.stimulation_frequency, current_time)

  if is_compartment_gets_stimulated:  
    
    n_nodes_x = variables.n_points_3D_mesh_global_x
    n_nodes_y = variables.n_points_3D_mesh_global_y
    n_nodes_z = variables.n_points_3D_mesh_global_z
    z_index_center = (int)(n_nodes_z/2)
    y_index_center = (int)(n_nodes_y/2)
    x_index_center = (int)(n_nodes_x/2)
    
    for k in range(n_nodes_z):
      if z_index_center-1 <= k <= z_index_center+1:  # use only 3 nodes in z direction from center
        for j in range(n_nodes_y):
          #if y_index_center-1 <= j <= y_index_center+1:  # use only 3 nodes in y direction from center
          if True:                                        # use all nodes
            for i in range(n_nodes_x):
              #if x_index_center-1 <= i <= x_index_center+1:  # use only 3 nodes in x direction from center
              if True:                                        # use all nodes
                
                key = ((i,j,k),0,0)        # key: ((x,y,z),nodal_dof_index,state_no)
                states[key] = variables.vm_value_stimulated
                #print("set states at ({},{},{}) to 40".format(i,j,k))

    #print("states: {}".format(states))
    #print("n_nodes: ({},{},{})".format(n_nodes_x, n_nodes_y, n_nodes_z))
    #print("n_nodes_global: {}, time_step_no: {}, current_time: {}, compartment_no: {}".format(n_nodes_global, time_step_no, current_time, compartment_no))
    #wait = input("Press any key to continue...")

  
# for debugging output show when the first 20 fibers will fire
if rank_no == 0 and not variables.disable_firing_output:
  print("\nDebugging output about compartment firing: Taking input from file \"{}\"".format(variables.firing_times_file))
  import timeit
  t_start = timeit.default_timer()
  
  first_stimulation_info = []
  
  n_firing_times = np.size(variables.firing_times,0)
  for compartment_no_index in range(variables.n_compartments):
    if compartment_no_index % 100 == 0:
      t_algebraic = timeit.default_timer()
      if t_algebraic - t_start > 100:
        print("Note: break after {}/{} compartments ({:.0f}%) because it already took {:.3f}s".format(compartment_no_index,variables.n_compartments,100.0*compartment_no_index/(variables.n_compartments-1.),t_algebraic - t_start))
        break
    
    first_stimulation = None
    for current_time in np.linspace(0,1./variables.stimulation_frequency*n_firing_times,n_firing_times):
      if compartment_gets_stimulated(compartment_no_index, variables.stimulation_frequency, current_time):
        first_stimulation = current_time
        break
    mu_no = get_motor_unit_no(compartment_no_index)
    first_stimulation_info.append([compartment_no_index,mu_no,first_stimulation])
  
  first_stimulation_info.sort(key=lambda x: 1e6+1e-6*x[1]+1e-12*x[0] if x[2] is None else x[2]+1e-6*x[1]+1e-12*x[0])
  
  print("First stimulation times")
  print("    Time  MU compartments")
  n_stimulated_mus = 0
  n_not_stimulated_mus = 0
  stimulated_fibers = []
  last_time = 0
  last_mu_no = first_stimulation_info[0][1]
  for stimulation_info in first_stimulation_info:
    mu_no = stimulation_info[1]
    fiber_no = stimulation_info[0]
    if mu_no == last_mu_no:
      stimulated_fibers.append(fiber_no)
    else:
      if last_time is not None:
        if len(stimulated_fibers) > 10:
          print("{:8.2f} {:3} {} (only showing first 10, {} total)".format(last_time,last_mu_no,str(stimulated_fibers[0:10]),len(stimulated_fibers)))
        else:
          print("{:8.2f} {:3} {}".format(last_time,last_mu_no,str(stimulated_fibers)))
        n_stimulated_mus += 1
      else:
        if len(stimulated_fibers) > 10:
          print("  never stimulated: MU {:3}, fibers {} (only showing first 10, {} total)".format(last_mu_no,str(stimulated_fibers[0:10]),len(stimulated_fibers)))
        else:
          print("  never stimulated: MU {:3}, fibers {}".format(last_mu_no,str(stimulated_fibers)))
        n_not_stimulated_mus += 1
      stimulated_fibers = [fiber_no]

    last_time = stimulation_info[2]
    last_mu_no = mu_no
    
  print("stimulated MUs: {}, not stimulated MUs: {}".format(n_stimulated_mus,n_not_stimulated_mus))

  t_end = timeit.default_timer()
  print("duration of assembling this list: {:.3f} s\n".format(t_end-t_start))  
  
####################################

# compute relative factors fr for compartments
def compute_compartment_relative_factors(mesh_node_positions, n_mesh_points_xy, n_mesh_points_z, fiber_data, motor_units):
  """
  Compute the relative factors, f_r, that are needed in the multidomain formulation as a weighting for compartments.
  Result is relative_factors[motor_unit_no][node_no] for the 3D mesh.
  :param mesh_node_positions:  list of (x,y,z) values, global node positions of the 3D mesh
  :param fiber_data: list of fibers, each fiber is a list of points, i.e. point = fiber_data[xy_index][z_index]
  :param motor_units: a list of dicts, settings for the motor units, [{"fiber_no": 0, "standard_deviation": 0.5, "maximum": 1}]
  """
  
  # list of fibers, fiber = list of points, point = list with 3 coordinate entries
  n_compartments = len(motor_units)
  n_points_fiber = len(fiber_data[0])

  # create relative factors for compartments
  #if rank_no == 0:
  #  print("determine relative factors for {} motor units:\n{}".format(n_compartments, motor_units))

  # determine approximate diameter of muscle at every point is z direction
  diameters = []
    
  # loop over points in z direction
  for z_index_mesh in range(n_mesh_points_z):
    
    z_index_fiber = int(z_index_mesh / (float)(n_mesh_points_z) * n_points_fiber)
    
    # get point on first and last fiber
    point0 = np.array(fiber_data[0][z_index_fiber])
    point4 = np.array(fiber_data[(variables.n_fibers_x-1)//2][z_index_fiber])
    point1 = np.array(fiber_data[variables.n_fibers_x-1][z_index_fiber])
    point2 = np.array(fiber_data[-variables.n_fibers_x][z_index_fiber])
    point5 = np.array(fiber_data[(-variables.n_fibers_x)//2][z_index_fiber])
    point3 = np.array(fiber_data[-1][z_index_fiber])
    
    # their distance is an approximation for the diameter
    distance01 = np.linalg.norm(point0 - point1)
    distance02 = np.linalg.norm(point0 - point2)
    distance03 = np.linalg.norm(point0 - point3)
    distance04 = np.linalg.norm(point0 - point4)
    distance05 = np.linalg.norm(point0 - point5)
    distance12 = np.linalg.norm(point1 - point2)
    distance13 = np.linalg.norm(point1 - point3)
    distance14 = np.linalg.norm(point1 - point4)
    distance15 = np.linalg.norm(point1 - point5)
    distance23 = np.linalg.norm(point2 - point3)
    distance24 = np.linalg.norm(point2 - point4)
    distance25 = np.linalg.norm(point2 - point5)
    distance34 = np.linalg.norm(point3 - point4)
    distance35 = np.linalg.norm(point3 - point5)
    distance45 = np.linalg.norm(point4 - point5)
    distance = max(distance01,distance02,distance03,distance04,distance05,distance12,distance13,distance14,distance15,distance23,distance24,distance25,distance34,distance35,distance45)
    diameters.append(distance)

  #print("diameters: {}".format(diameters))

  # create data structure with 0
  relative_factors = np.zeros((n_compartments, len(mesh_node_positions)))   # each row is one compartment

  # loop over nodes of mesh
  for node_no,node_position in enumerate(mesh_node_positions):
    node_position = np.array(node_position)
    
    z_index_mesh = int((float)(node_no) / n_mesh_points_xy)
    z_index_fiber = int(z_index_mesh / (float)(n_mesh_points_z) * n_points_fiber)
    
    # loop over motor units
    for motor_unit_no,motor_unit in enumerate(motor_units):
      
      # find point on fiber that is closest to current node
      fiber_no = motor_unit["fiber_no"]
      if fiber_no >= len(fiber_data):
        new_fiber_no = fiber_no % len(fiber_data)
        if node_no == 0:
          print("\033[0;31mError with motor unit {} around fiber {}, only {} fibers available, now using fiber {} % {} = {} instead.\033[0m".format(motor_unit_no, fiber_no, len(fiber_data), fiber_no, len(fiber_data), new_fiber_no))
        fiber_no = new_fiber_no
      
      min_distance = None
      search_range = int(1 / (float)(n_mesh_points_z) * n_points_fiber)
      search_range = max(10,search_range)
      z_start = max(0,z_index_fiber - search_range)
      z_end = min(n_points_fiber, z_index_fiber + search_range)
      
      #print("node_position: {}, z_index_fiber: {}, fiber at z index: {}, fiber: {}".format(node_position, z_index_fiber, fiber_data[fiber_no][z_index_fiber], fiber_data[fiber_no][z_start:z_end]))
      #print("search_range: {}".format(search_range))
      
      for k,fiber_point in enumerate(fiber_data[fiber_no][z_start:z_end]):
        d = np.array(fiber_point) - node_position
        distance = np.inner(d,d)
        if min_distance is None or distance < min_distance:
          min_distance = distance
          #print("node_position {}, fiber_point {}, d={}, |d|={}".format(node_position, fiber_point, d, np.sqrt(distance)))
      
      distance = np.sqrt(min_distance)
      
      # compute value as gaussian with given standard_deviation and maximum
      standard_deviation = motor_unit["standard_deviation"]*diameters[z_index_mesh]
      gaussian = scipy.stats.norm(loc = 0., scale = standard_deviation)
      value = gaussian.pdf(distance)*standard_deviation*np.sqrt(2*np.pi)*motor_unit["maximum"]
      relative_factors[motor_unit_no][node_no] += value
      
      #print("motor unit {}, fiber {}, distance {}, value {}".format(motor_unit_no, fiber_no, distance, value))

  return relative_factors

####################################
# load relative factors for motor units

if False:
  # determine relative factor fields fr(x) for compartments
  if not os.path.exists(variables.relative_factors_file):

    # the file does not yet exist, create it on rank 0
    if rank_no == 0: 
      
      mesh_node_positions = variables.meshes["3Dmesh"]["globalNodePositions"]
      n_points_global = variables.meshes["3Dmesh"]["nPointsGlobal"]
      n_mesh_points_xy = n_points_global[0]*n_points_global[1]
      n_mesh_points_z = n_points_global[2]
      
      print("Computing the relative MU factors, f_r, for {} motor units and {} mesh nodes, {} fibers. This may take a while ...".format(len(variables.motor_units), len(mesh_node_positions), len(variables.fibers)))
      variables.relative_factors = compute_compartment_relative_factors(mesh_node_positions, n_mesh_points_xy, n_mesh_points_z, variables.fibers, variables.motor_units)
      if rank_no == 0:
        print("Save relative factors to file \"{}\".".format(variables.relative_factors_file))
        with open(variables.relative_factors_file, "wb") as f:
          pickle.dump(variables.relative_factors, f)
    else:
      # wait until file is created on rank 0
      while not os.path.exists(variables.relative_factors_file):
        time.sleep(1)

  if os.path.exists(variables.relative_factors_file):
    with open(variables.relative_factors_file, "rb") as f:
      if rank_no == 0:
        print("Load relative factors, f_r, from file \"{}\"".format(variables.relative_factors_file))
      variables.relative_factors = pickle.load(f, encoding='latin1')
  else:
    print("\033[0;31mError: Could not load relative factors file \"{}\"\033[0m".format(variables.relative_factors_file))
    quit()

  # debugging output
  if rank_no == 0 and not variables.disable_firing_output:
    for i,factors_list in enumerate(variables.relative_factors.tolist()):
      print("MU {}, maximum fr: {}".format(i,max(factors_list)))


muscle_spindle_node_nos =  [4, 22, 58, 75] #[12,112,212,312,412]
# f = open("spindles_nofs.csv", "a")

# for muscle_spindle_no in range(variables.n_muscle_spindles):
#   i = random.randrange(0,variables.el_x+1)
#   j = random.randrange(0,variables.el_y+1)
#   k = random.randrange(0,variables.el_z+1)
  
#   dof_no_global = k*(variables.el_x+1)*(variables.el_y+1) + j*(variables.el_x+1) + i
#   muscle_spindle_node_nos.append(dof_no_global)
#   f.write(str(dof_no_global) + " out of " + str((variables.el_x+1)*(variables.el_y+1)*(variables.el_z+1)) + "\n")
# f.close()
# determine positions of golgi tendon organs
golgi_tendon_organ_node_nos = []
for golgi_tendon_organ_no in range(variables.n_golgi_tendon_organs):
  i = random.randrange(0,variables.bs_x)
  j = random.randrange(0,variables.bs_y)
  if golgi_tendon_organ_no % 2 == 0:
    k = int(0.1*variables.bs_z)
  else:
    k = int(0.9*variables.bs_z)
  
  dof_no_global = k*variables.bs_x*variables.bs_y + j*variables.bs_x + i
  golgi_tendon_organ_node_nos.append(dof_no_global)

