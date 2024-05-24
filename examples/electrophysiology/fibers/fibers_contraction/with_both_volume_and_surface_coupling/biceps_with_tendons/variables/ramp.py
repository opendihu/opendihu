
# scenario name for log file and output directory under "out"
case_name = "ramp100N"
scenario_name = "ramp"


# time parameters
# --------------------
end_time = 100                    # [ms] end time of the simulation
dt_0D = 2e-4                      # [ms] timestep width of ODEs (2e-3)
dt_1D = 4e-4                      # [ms] timestep width of diffusion (4e-3)
dt_splitting = 4e-4               # [ms] overall timestep width of strang splitting (4e-3)
dt_3D = 1e-2                      # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG, this has to be the same value as in the precice_config.xml
output_timestep = 50              # [ms] timestep for output files, 5.0
output_timestep_fibers = int(dt_3D/dt_splitting)*output_timestep   # [ms] timestep for fiber output
output_timestep_big = 1.0            # [ms] timestep for output big files of 3D EMG, 100

# material parameters
# --------------------
# quantities in CellML unit system
sigma_f = 8.93              # [mS/cm] conductivity in fiber direction (f)
Conductivity = sigma_f      # [mS/cm] sigma, conductivity
Am = 500.0                  # [cm^-1] surface area to volume ratio
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)

rho = 10                    # [1e-4 kg/cm^3] density of the muscle (density of water)

c1 = 3.176e-10              # [N/cm^2]
c2 = 1.813                  # [N/cm^2]
b  = 1.075e-2               # [N/cm^2] anisotropy parameter
d  = 9.1733                 # [-] anisotropy parameter
muscle_material_parameters = [c1, c2, b, d]   # material parameters
pmax = 7.3                  # [N/cm^2] maximum isometric active stress


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
tendon_material_parameters = [c, ca, ct, cat, ctt, mu, k1, k2]
# for debugging, b = 0 leads to normal Mooney-Rivlin
#b = 0
force = 1000.0
constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
bottom_traction = [0.0,0.0,0.0]        # [1 N]

# timing and activation parameters
# -----------------
import random
random.seed(0)  # ensure that random numbers are the same on every rank
motor_units = [
  {"radius": 40.00, "activation_start_time": 0.0, "stimulation_frequency": 23.92, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #1    # low number of fibers
  {"radius": 42.35, "activation_start_time": 0.2, "stimulation_frequency": 23.36, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #2
  {"radius": 45.00, "activation_start_time": 0.4, "stimulation_frequency": 23.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #3
  {"radius": 48.00, "activation_start_time": 0.6, "stimulation_frequency": 22.46, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #4
  {"radius": 51.42, "activation_start_time": 0.8, "stimulation_frequency": 20.28, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #5
  {"radius": 55.38, "activation_start_time": 1.0, "stimulation_frequency": 16.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #6 
  {"radius": 60.00, "activation_start_time": 1.2, "stimulation_frequency": 12.05, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #7
  {"radius": 65.45, "activation_start_time": 1.4, "stimulation_frequency": 10.03, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #8
  {"radius": 72.00, "activation_start_time": 1.6, "stimulation_frequency": 8.32,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #9
  {"radius": 80.00, "activation_start_time": 1.8, "stimulation_frequency": 7.66,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]}, #10   # high number of fibers
]
# note: negative start time is the same as zero, it is just there for debugging. Delete the minus signs to get a ramp


# The values of dt_3D and end_time have to be also defined in "precice-config.xml" with the same value (the value is only significant in the precice-config.xml, the value here is used for output writer time intervals)
# <max-time value="100.0"/>           <!-- end time of the whole simulation -->
# <time-window-size value="1e0"/>   <!-- timestep width dt_3D -->

# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements

import opendihu

# parameters for the contraction program
if "contraction" in opendihu.program_name:
  sampling_stride_x = 2
  sampling_stride_y = 2
  sampling_stride_z = 74
  # good values: divisors of 1480: 1480 = 1*1480 = 2*740 = 4*370 = 5*296 = 8*185 = 10*148 = 20*74 = 37*40 

else:
  # parameters for the fibers_with_3d program
  
  sampling_stride_x = 1
  sampling_stride_y = 1
  sampling_stride_z = 20
  # good values: divisors of 1480: 1480 = 1*1480 = 2*740 = 4*370 = 5*296 = 8*185 = 10*148 = 20*74 = 37*40 

distribute_nodes_equally = False     # (default: False)
# True: set high priority to make subdomains have approximately equal number of fibers but creates tiny remainder elements inside the subdomains
# False: make elements more equally sized, this can lead to a slight imbalance in the number of fibers per subdomain

# input files
import os
input_directory   = os.path.join(os.environ["OPENDIHU_HOME"], "examples/electrophysiology/input")

#fiber_file        = input_directory + "/left_biceps_brachii_7x7fibers.bin"
fiber_file        = input_directory + "/left_biceps_brachii_9x9fibers.bin"
firing_times_file = input_directory + "/MU_firing_times_always.txt"
fiber_distribution_file = input_directory + "/MU_fibre_distribution_10MUs.txt"
#fiber_distribution_file = input_directory + "/MU_fibre_distribution_10MUs_13x13.txt"
cellml_file       = input_directory + "/new_slow_TK_2014_12_08.cellml"

# EMG solver parameters
emg_solver_type = "cg"              # solver and preconditioner for the 3D static Bidomain equation that solves the intra-muscular EMG signal
emg_preconditioner_type = "none"    # preconditioner
emg_initial_guess_nonzero = True    # If the initial guess for the emg linear system should be set to the previous solution
emg_solver_maxit = 1e4              # maximum number of iterations for the static bidomain solver
emg_solver_reltol = 1e-5            # relative tolerance for solver
emg_solver_abstol = 1e-5            # absolute tolerance for solver

# other options
states_output = False           # if the 0D states should be written (this produces large output files)
paraview_output = True         # produce output files for paraview
adios_output = False           # produce adios output files
exfile_output = False          # produce exfiles output files
python_output = False          # produce python output files


# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(fiber_no, mu_no):
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  r = motor_units[mu_no]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  return 2./r

def get_cm(fiber_no, mu_no):
  return Cm
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity

def get_specific_states_call_frequency(fiber_no, mu_no):
  stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
  return stimulation_frequency*1e-3

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  return motor_units[mu_no % len(motor_units)]["jitter"]

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  return motor_units[mu_no % len(motor_units)]["activation_start_time"]
