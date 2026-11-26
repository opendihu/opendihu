# Configuration for scons build system
#
# For each package the following variables are available:
# <PACKAGE>_DIR         Location of the package, must contain subfolders "include" and "lib" or "lib64" with header and library files.
# <PACKAGE>_INC_DIR     Location of (*.h) header files
# <PACKAGE>_LIB_DIR     Location of (*.a) libraries
# <PACKAGE>_LIBS        List of libraries, optional since the standard names are already hardcoded.
# <PACKAGE>_DOWNLOAD    Download, build and use a local copy of the package.
# <PACKAGE>_REDOWNLOAD  Force update of previously downloaded copy. For that <PACKAGE>_DOWNLOAD has to be also true.
# <PACKAGE>_REBUILD     Force a new build of the package without redownloading it if already has been downloaded earlier.
#
# You do one of the following:
# 1. Not specify any of the variables. Then standard locations in dependencies as well as /usr, /usr/local are searched.
# 2. Specify <PACKAGE>_DIR to directly give the base directory to the package's location. Do this to e.g. use system provided libraries.
# 3. Specify <PACKAGE>_INC_DIR and <PACKAGE>_LIB_DIR to point to the header and library directories. They are usually named "include" and "lib".
# 4. Set <PACKAGE>_DOWNLOAD=True or additionally <PACKAGE>_REDOWNLOAD=True to let the build system download and install everything on their own.

# set compiler to use
cc = "gcc"         # C compiler
CC = "g++"         # C++ compiler
mpiCC = "mpic++"   # MPI C++ wrapper
cmake="cmake"      # cmake command

# PETSc, True: download and build PETSc into dependencies/petsc, this also downloads and installs MUMPS (direct solver package)
# and its dependencies PT-Scotch, SCAlapack, ParMETIS, METIS (takes a while).
# False: use an existing PETSc given by PETSC_DIR. On Debian/Ubuntu install the system PETSc with `sudo apt install petsc-dev`,
# it includes MUMPS and HYPRE. Without root rights, set PETSC_DOWNLOAD = True.
PETSC_DOWNLOAD = False
PETSC_DIR = "/usr/lib/petsc"

HDF5_DOWNLOAD = True
NLOHMANNJSON_DOWNLOAD = False

# Python 3.12, True: download and build python into dependencies/python (takes a while).
# False: use the existing system python3 (or PYTHON_DIR/bin/python3 if PYTHON_DIR is set). This needs the
# python C-API headers and shared library, e.g. `sudo apt install python3-dev` on Ubuntu.
PYTHON_DOWNLOAD = False

# Python packages - they are now all combined with the option PYTHONPACKAGES_DOWNLOAD
# True: pip-install numpy, scipy, matplotlib, ... into the downloaded python (only works with PYTHON_DOWNLOAD = True).
# False: use the packages of the existing python, install them yourself (at least numpy is required for the check).
PYTHONPACKAGES_DOWNLOAD = False

# Base64, encoding library for binary vtk (paraview) output files
BASE64_DOWNLOAD = True

# Google Test, testing framework, not needed on Hazelhen
GOOGLETEST_DOWNLOAD = True

# SEMT, library for symbolic differentiation
SEMT_DOWNLOAD = True

# EasyLoggingPP, provides logging facilities
EASYLOGGINGPP_DOWNLOAD = True

# MegaMol, visualization framework of VISUS, optional, needs ADIOS2
ADIOS_DOWNLOAD = False # ADIOS2, adaptable I/O library, needed for interfacing MegaMol
MEGAMOL_DOWNLOAD = False    # install MegaMol from official git repo, but needed is the private repo, ask Tobias Rau for access to use MegaMol with opendihu

# Vc, vectorization types and C++ utility to produce vectorized code (but does not support AVX-512)
# std::experimental::simd supports AVX-512, but requires C++17. Therefore the package std_simd includes a compatibility script that falls back to Vc, if C++17 is not available.
VC_DOWNLOAD = True
STD_SIMD_DOWNLOAD = True

# xbraid, used for parallel-in time methods (currently only on branch `xbraid`)
XBRAID_DOWNLOAD = True

# OpenCOR, utility to view CellML models and to convert them from xml format to c code
OPENCOR_DOWNLOAD = True

# preCICE coupling library
# libxml2 is a dependency required by preCICE which can be installed using `sudo apt install libxml2-dev`. For historical reasons we provide the option to download from a tarball by setting LIBXML2_DOWNLOAD = True.  
# preCICE is an optional dependency to run partitioned simulations, i.e., multiple OpenDiHu solvers or OpenDiHu solvers coupled to external software. 
# You can install preCICE by setting PRECICE_DOWNLOAD = True. If you already have built preCICE in your system you can point to it to avoid having multiple installations of preCICE. In that case keep PRECICE_DOWNLOAD = False and set PRECICE_INC_DIR PRECICE_LIB_DIR
LIBXML2_DOWNLOAD = False
PRECICE_DOWNLOAD = False
#PRECICE_INC_DIR = "/path/to/precice/include"
#PRECICE_LIB_DIR = "/path/to/precice/build"

# MPI
# MPI is normally detected by running the mpicc command. If this is not available, you can provide the MPI_DIR manually.
#MPI_DIR = "/usr/lib/openmpi"    # standard path for openmpi on ubuntu 16.04
MPI_DIR = "/usr/lib/x86_64-linux-gnu/openmpi"    # standard path for openmpi on ubuntu >= 18.04

# Vectorized code for matrix assembly
# Set to True for fastest code, set to False for faster compilation, note this works only using the Vc code, not std::simd code (only for GCC < 9)
USE_VECTORIZED_FE_MATRIX_ASSEMBLY = False
if USE_VECTORIZED_FE_MATRIX_ASSEMBLY:
  print("Note, USE_VECTORIZED_FE_MATRIX_ASSEMBLY is True in user-variables.scons.py, this means faster programs but longer compilation times.\n")

# Use the implementation of std::simd instead of Vc to support AVX-512. This automatically sets the C++ standard from C++14 to C++17
USE_STDSIMD = False
if USE_STDSIMD:
  print("Note, USE_STDSIMD is True and, thus, c++17 will be used.");

# -------------------------------------------------------------------------
# automatically set MPI_DIR for other systems, like ubuntu 16.04 and Debian
try:
  import lsb_release
  lsb_info = lsb_release.get_lsb_information()   # get information about ubuntu version, if available
  if "RELEASE" in lsb_info:
    if lsb_info["RELEASE"] == "16.04":
      MPI_DIR="/usr/lib/openmpi"   # this is the standard path on ubuntu 16.04
except:
  pass

try:
  import platform
  if 'debian' in platform.dist():
    MPI_DIR="/usr/lib/x86_64-linux-gnu/openmpi"    # path for debian (on Aaron's workstation)
except:
  pass

try:
  # use value of environment variable 'MPI_HOME' if it is set
  import os
  if os.environ.get("MPI_HOME") is not None:
    MPI_DIR = os.environ.get("MPI_HOME")
    
  # for Travis CI, build MPI ourselves
  if os.environ.get("TRAVIS") is not None:
    print("Travis CI detected, del MPI_DIR")
    del MPI_DIR
    MPI_DOWNLOAD=True
  
  import socket

  # special settings on cluster "lead"
  if "lead" in socket.gethostname():
    MPI_DIR = os.environ["MPI_HOME"]
 
  # special settings on supercomputer Hawk 
  elif "hawk" in os.environ["SITE_PLATFORM_NAME"]:
    if "MPT_ROOT" in os.environ:
      MPI_DIR = os.environ["MPT_ROOT"]
    else:
      MPI_DIR = os.environ["MPI_ROOT"]
    
    MPI_IGNORE_MPICC = True
    LAPACK_DOWNLOAD = False
    PETSC_DOWNLOAD = False
    PETSC_DIR = os.environ["PETSC_ROOT"]
    PYTHONPACKAGES_DOWNLOAD = False
    GOOGLETEST_DOWNLOAD = True 
    XBRAID_DOWNLOAD = True
    ADIOS_DOWNLOAD = False
    ADIOS_DIR = os.environ["ADIOS2_ROOT"]
except:
  pass

# download and build debugging MPI version
if False:
  del MPI_DIR
  MPI_DOWNLOAD = True
  MPI_IGNORE_MPICC = True    # this downloads and builds openmpi

#PETSC_DEBUG = True            # this enables debugging flags such that valgrind memcheck can track MPI errors

# specialized settings for supercomputer HazelHen, this is left here in case there will be another Cray supercomputer
import os
if os.environ.get("PE_ENV") is not None:  # if on hazelhen
  cc = "cc"   # C compiler wrapper
  CC = "CC"   # C++ compiler wrapper
  mpiCC = "CC"  # mpi C++ compiler wrapper
  cmake = "/lustre/cray/ws8/ws/icbbnmai-opendihu1/cmake/cmake-3.13.2-Linux-x86_64/bin/cmake"

  # use cray-pat for profiling
  USE_CRAY_PAT = False

  # use -hpl option with cray compiler to create an optimization program library
  USE_HPL = False

  # do not use googletest
  GOOGLETEST_DOWNLOAD = False  

  # do not use buggy python packages
  PYTHONPACKAGES_DOWNLOAD = False
  
  # steps for getting started on HazelHen:
  #   module swap PrgEnv-cray/6.0.4 PrgEnv-gnu  # to switch to GNU programming environment, however also Intel and Cray environments work
  #   module load cray-libsci
  #   module load cray-petsc  (or cray-petsc-64 for big data)

# Optionally, you can create a "user-variables.local.scons.py" file in this directory and
# use this for settings that are specific to your own machine, e.g. paths to own installations of preCICE.
# That file is not tracked by git and overrides the values above.
import os
if os.path.isfile("user-variables.local.scons.py"):
  exec(open("user-variables.local.scons.py").read())


