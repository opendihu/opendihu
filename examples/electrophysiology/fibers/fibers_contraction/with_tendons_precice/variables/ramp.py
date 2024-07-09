
case_name = "default"
scenario_name = "muscle"
precice_config_file ="precice_config.xml"

# material parameters
# --------------------
# quantities in mechanics unit system
rho = 10                    # [1e-4 kg/cm^3] density of the muscle (density of water)

# Mooney-Rivlin parameters [c1,c2,b,d] of c1*(Ibar1 - 3) + c2*(Ibar2 - 3) + b/d (λ - 1) - b*ln(λ)
# Heidlauf13: [6.352e-10 kPa, 3.627 kPa, 2.756e-5 kPa, 43.373] = [6.352e-11 N/cm^2, 3.627e-1 N/cm^2, 2.756e-6 N/cm^2, 43.373], pmax = 73 kPa = 7.3 N/cm^2
# Heidlauf16: [3.176e-10 N/cm^2, 1.813 N/cm^2, 1.075e-2 N/cm^2, 9.1733], pmax = 7.3 N/cm^2

c1 = 3.176e-10              # [N/cm^2]
c2 = 1.813                  # [N/cm^2]
b  = 1.075e-2               # [N/cm^2] anisotropy parameter
d  = 9.1733                 # [-] anisotropy parameter

# for debugging, b = 0 leads to normal Mooney-Rivlin
b = 0

material_parameters = [c1, c2, b, d]   # material parameters
pmax = 7.3                  # [N/cm^2] maximum isometric active stress (30-40)

# load
constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
bottom_traction = [0.0,0.0,0.0]        # [N]

# Monodomain parameters
# --------------------
# quantities in CellML unit system
Conductivity = 3.828      # [mS/cm] sigma, conductivity
Am = 500.0                  # [cm^-1] surface area to volume ratio (this is not used, instead values of motor_units are used)
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)
# diffusion prefactor = Conductivity/(Am*Cm)

# timing and activation parameters
# -----------------
# motor units from paper Klotz2019 "Modelling the electrical activity of skeletal muscle tissue using a multi‐domain approach"
import random
random.seed(0)  # ensure that random numbers are the same on every rank
# radius: [μm], stimulation frequency [Hz], jitter [-]
motor_units = [
  {"radius": 40.00, "activation_start_time": 0.0, "stimulation_frequency": 23.92, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},    # low number of fibers
  {"radius": 42.35, "activation_start_time": 0.2, "stimulation_frequency": 23.36, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 45.00, "activation_start_time": 0.4, "stimulation_frequency": 23.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 48.00, "activation_start_time": 0.6, "stimulation_frequency": 22.46, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 51.42, "activation_start_time": 0.8, "stimulation_frequency": 20.28, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 55.38, "activation_start_time": 1.0, "stimulation_frequency": 16.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 60.00, "activation_start_time": 1.2, "stimulation_frequency": 12.05, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 65.45, "activation_start_time": 1.4, "stimulation_frequency": 10.03, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 72.00, "activation_start_time": 1.6, "stimulation_frequency": 8.32,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 80.00, "activation_start_time": 1.8, "stimulation_frequency": 7.66,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},    # high number of fibers
]

# timing parameters
# -----------------
end_time = 100                    # [ms] end time of the simulation
dt_0D = 2e-4                      # [ms] timestep width of ODEs (2e-3)
dt_1D = 4e-4                      # [ms] timestep width of diffusion (4e-3)
dt_splitting = 4e-4               # [ms] overall timestep width of strang splitting (4e-3)
dt_3D = 1e-2                      # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG, this has to be the same value as in the precice_config.xml
output_timestep = 50              # [ms] timestep for output files, 5.0
output_timestep_fibers = int(dt_3D/dt_splitting)*output_timestep   # [ms] timestep for fiber output
output_timestep_big = 1.0            # [ms] timestep for output big files of 3D EMG, 100

# input files
fiber_file = "../../../../input/left_biceps_brachii_9x9fibers.bin"
#fiber_file = "../../../../input/left_biceps_brachii_13x13fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"
firing_times_file = "../../../../input/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file = "../../../../input/MU_fibre_distribution_10MUs.txt"
cellml_file             = "../../../../input/new_slow_TK_2014_12_08.c"

# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
sampling_stride_x = 2
sampling_stride_y = 2
sampling_stride_z = 74        # good values: divisors of 1480: 1480 = 1*1480 = 2*740 = 4*370 = 5*296 = 8*185 = 10*148 = 20*74 = 37*40 

# Tolerance value in the element coordinate system of the 3D elements, [0,1]^3
# when a fiber point is still considered part of the element.
# Try to increase this such that all mappings have all points.
mapping_tolerance = 0.5

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
disable_firing_output = False

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(fiber_no, mu_no):
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  r = motor_units[mu_no]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  return 2./r
  #return Am

def get_cm(fiber_no, mu_no):
  return Cm
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity

def get_specific_states_call_frequency(fiber_no, mu_no):
  stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
  return stimulation_frequency*1e-3

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  #return 0
  return motor_units[mu_no % len(motor_units)]["jitter"]

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  return motor_units[mu_no % len(motor_units)]["activation_start_time"]
