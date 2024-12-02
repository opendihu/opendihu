# Modified diffusion 1D
# Compute

import sys
import argparse

rank_no = (int)(sys.argv[-2])
n_ranks = (int)(sys.argv[-1])

parser = argparse.ArgumentParser(description='settings_my_dynamic_solver')

parser.add_argument("--enable_checkpointing", help="Enable checkpointing", default=False, action="store_true")
parser.add_argument("--checkpointing_dir", help="Set a Checkpointing directory, only used if Checkpointing is enabled", default="states")
parser.add_argument("--checkpointing_interval", help="Set a specific Checkpointing Timestep Interval, only used if Checkpointing is enabled", default=10)
parser.add_argument("--checkpointing_system", help="Set a specific Checkpointing System, only used if Checkpointing is enabled", default="hdf5-combined")
parser.add_argument("--checkpointing_autorestore", help="Enable automatic restore if a Checkpoint is found, only used if Checkpointing is enabled", default=False, action="store_true")
parser.add_argument("--checkpointing_restore_checkpoint", help="Restore a specific Checkpoint, only used if Checkpointing is enabled", default=None)

args, other_args = parser.parse_known_args(args=sys.argv[:-2])
if len(other_args) != 0 and rank_no == 0:
    print("Warning: These arguments were not parsed by the settings python file\n  " + "\n  ".join(other_args), file=sys.stderr)

n = 5   # number of elements

config = {
  "logFormat": "csv",
  "Solvers": {
    "linearSolver": {
      "solverType": "gmres",          # the solver type, refer to PETSc documentation about implemented solvers and preconditioners (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html)
      "preconditionerType": "none",   # the preconditioner type (https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType)
      "relativeTolerance": 1e-15,     # the relative tolerance of the linear solver
      "absoluteTolerance": 1e-10,     # 1e-10 absolute tolerance of the residual
      "maxIterations": 1e4,           # maximum number of iterations of the linear solver
      "dumpFilename": "",             # a filename to dump the system matrix and rhs in every timestep, "" to disable
      "dumpFormat": "matlab",         # the file format of the dump file, "matlab", "ascii" or "default"
    }
  },
  "MyNewTimesteppingSolver": {        # this is the name of the solver, as given in the constructor to the timestepping object
    "myOption": 42,                   # example option that is parsed in the constructor
    "option1": "blabla",              # another example option that is parsed in the data object

    # option for the timestepping of MyNewTimesteppingSolver
    "endTime": 60.0,                  # end time of the simulation
    "timeStepWidth": 2.0,             # time step width
    "timeStepOutputInterval": 1,      # how often to print the current timestep to console

    # settings for the nested solver
    "ImplicitEuler" : {
      "numberTimeSteps": 10,
      "endTime": 0.1,
      "initialValues": [2,5,10,4,-2,2],    # the initial values
      "dirichletBoundaryConditions": {}, # Dirichlet boundary conditions as dict
      "dirichletOutputFilename":  None,  # filename for a vtp file that contains the Dirichlet boundary condition nodes and their values, set to None to disable
      "inputMeshIsGlobal": True,         # initial values and BC's are given for all dofs, even if executed in parallel
      "timeStepOutputInterval": 1,       # how often to print the current timestep to console
      "nAdditionalFieldVariables": 0,    # for more complex nested solvers, the number of additional field variables that will be transferred without being touched
      "solverName": "linearSolver",      # the solver to use, referes to what was defined under "Solvers"

      "FiniteElementMethod" : {
        # mesh
        "nElements": n,                 # number of elements
        "physicalExtent": 4.0,          # the physical size of the domain

        "solverName": "linearSolver",   # the solver to use, referes to what was defined under "Solvers"
        "prefactor": 5.0,               # the prefactor 'c' of 'du/dt = c du^2/dx^2'
        "inputMeshIsGlobal": True,      # boundary conditions are given as global indices
        "nAdditionalFieldVariables": 0, # for more complex nested solvers, the number of additional field variables that will be transferred without being touched
      },
    },

    "OutputWriter" : [
      {"format": "Paraview", "outputInterval": 1, "filename": "out/paraview", "binary": True, "fixedFormat": False, "onlyNodalValues": True, "combineFiles": False},
      {"format": "PythonFile", "filename": "out/python", "outputInterval": 1, "binary":False, "onlyNodalValues": True}
    ]
  }
}

if args.enable_checkpointing:
    config["checkpointing"] = {
        "directory": args.checkpointing_dir,
        "interval": args.checkpointing_interval,
        "type": args.checkpointing_system,
        "autoRestore": args.checkpointing_autorestore,
        "checkpointToRestore": args.checkpointing_restore_checkpoint,
    }
