#pragma once

#include <Python.h>  // has to be the first included header
#include <petscsys.h>
#include <list>
#include <memory>
#include <map>
#include "partition/partition_manager.h"

// forward declaration
namespace Mesh { class Manager; }
namespace Solver { class Manager; }

/** This class contains global variables (mesh manager and solver manager) and a python config dictionary
 *  which can be a sub-dict of the global config. Different objects of this class are used for the methods
 *  throughout this program which have a different context in terms of their config.
 */
class DihuContext
{
public:
  ///! constructor, initialize context, parse command line parameters and input file
  DihuContext(int argc, char *argv[], bool doNotFinalizeMpi=false, bool settingsFromFile=true);

  ///! constructor for test cases
  DihuContext(int argc, char *argv[], std::string pythonSettings, bool doNotFinalizeMpi=true);

  ///! copy-constructor
  DihuContext(const DihuContext &rhs);

  ///! return a context object with config originated at child node with given key
  DihuContext operator[](std::string keyString);

  ///! return the top-level python config object
  PyObject *getPythonConfig() const;
  
  ///! create a context object, like with the operator[] but with given settings
  DihuContext createSubContext(PyObject *settings);

  ///! return the mesh manager object that contains all meshes
  std::shared_ptr<Mesh::Manager> meshManager() const;

  ///! return the solver manager object that contains all solvers
  std::shared_ptr<Solver::Manager> solverManager() const;

  ///! return the partition manager object that creates partitionings
  std::shared_ptr<Partition::Manager> partitionManager() const;

  ///! get the own MPI rank no in the world communicator
  int ownRankNo() const;

  ///! destructor
  ~DihuContext();

private:
  ///! read in file and execute python script and store global variables
  void loadPythonScriptFromFile(std::string filename);

  ///! execute python script and store global variables
  void loadPythonScript(std::string text);

  ///! initiaize the easylogging++ framework
  void initializeLogging(int argc, char *argv[]);

  PyObject *pythonConfig_;    ///< the top level python config dictionary of the current context (i.e. may be a sub-dict of the global config)

  static std::shared_ptr<Mesh::Manager> meshManager_;   ///< object that saves all meshes that are used
//  static std::shared_ptr<Solver::Manager> solverManager_; ///< object that saves all solver configurations that are used
  static std::map<int, std::shared_ptr<Solver::Manager>> solverManagerForThread_;  ///< object that saves all solver configurations that are used, different for each thread

  static std::shared_ptr<Partition::Manager> partitionManager_;  ///< partition manager object that creates and manages partitionings
  
  static bool initialized_;  ///< if MPI, Petsc and easyloggingPP is already initialized. This needs to be done only once in the program.
  static int nObjects_;   ///< number of objects of DihuContext, if the last object gets destroyed, call MPI_Finalize or MPI_Barrier, depending on doNotFinalizeMpi
  bool doNotFinalizeMpi_;  ///< when the last object gets destroyed, either MPI_Finalize() is called (should be used) or MPI_Barrier (only needed in testcases where MPI context needs to be used for the next test cases)
};
