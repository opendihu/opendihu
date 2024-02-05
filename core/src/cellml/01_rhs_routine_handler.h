#pragma once

#include <Python.h> // has to be the first included header

#include <Python.h>
#include <vector>
#include <list>

#include "control/dihu_context.h"
#include "output_writer/manager.h"
#include "cellml/00_cellml_adapter_base.h"

/** This is a class that handles the right hand side routine of CellML
 * It needs access on the system size and some C-code of the right hand side.
 * It provides a method to call the right hand side.
 *
 *  Naming:
 *   Algebraic (opendihu) = KNOWN (OpenCMISS) = Algebraic (OpenCOR)
 *   Parameter (opendihu, OpenCMISS) = KNOWN (OpenCMISS), in OpenCOR also
 * algebraic Constant - these are constants that are only present in the source
 * files State: state variable Rate: the time derivative of the state variable,
 * i.e. the increment value in an explicit Euler stepping
 */
template <int nStates, int nAlgebraics_, typename FunctionSpaceType>
class RhsRoutineHandler
    : public CellmlAdapterBase<nStates, nAlgebraics_, FunctionSpaceType> {
public:
  //! constructor
  using CellmlAdapterBase<nStates, nAlgebraics_,
                          FunctionSpaceType>::CellmlAdapterBase;

  //! if when using "vc" as optimizationType_, the exp() function should be
  //! approximated, this is faster
  bool approximateExponentialFunction();

  //! get the optimization type as it was specified in the settings
  std::string optimizationType();

  //! load a given shared object library (<file>.so) and return the handle
  static void *loadRhsLibraryGetHandle(std::string libraryFilename);

protected:
  //! given a normal cellml source file for rhs routine create a second file for
  //! multiple instances. @return: if successful
  bool createSimdSourceFile(std::string &simdSourceFilename);

  //! given a normal cellml source file for rhs routine create a third file for
  //! gpu acceloration. @return: if successful
  bool createGPUSourceFile(std::string &gpuSourceFilename);

  //! initialize the rhs routine, either directly from a library or compile it
  void initializeRhsRoutine();

  //! load the library (<file>.so) that was created earlier, store
  bool loadRhsLibrary(std::string libraryFilename);

  //! create the source filename using the CellmlSourceCodeGenerator, then
  //! compile to library
  void createLibraryOnOneRank(std::string libraryFilename,
                              const std::vector<int> &nInstancesRanks);

  std::string
      sourceToCompileFilename_;  //< filename of the processed source file that
                                 //will be used to compile the library
  std::string optimizationType_; //< type of generated file, e.g. "simd", "gpu",
                                 //"openmp"
  bool approximateExponentialFunction_; //< when using "vc" as
                                        //optimizationType_, the exp() function
                                        //should be approximated, this is faster
  int maximumNumberOfThreads_; //< when using "openmp" as optimizationType_, the
                               //maximum number of threads to use, 0 means no
                               //restriction
  bool useAoVSMemoryLayout_;   //< which memory layout to use for the vc
                               //optimization type,
                               //true=Array-of-Vectorized-Struct,
                               //false=Struct-of-Vectorized-Array

  void (*rhsRoutine_)(
      void *context, double t, double *states, double *rates,
      double *algebraics,
      double *
          parameters); //< function pointer to the rhs routine that can compute
                       //several instances of the problem in parallel. Data is
                       //assumed to contain values for a state contiguously,
                       //e.g. (state[1], state[1], state[1], state[2], state[2],
                       //state[2], ...). The first parameter is a this pointer.

  //! helper rhs routines
  void (*rhsRoutineGPU_)(
      void *context, double t, double *states, double *rates,
      double *algebraics,
      double *parameters); //< function pointer to a gpu processed rhs function
                           //that is passed as dynamic library. Data is assumed
                           //to contain values for a state contiguously, e.g.
                           //(state[1], state[1], state[1], state[2], state[2],
                           //state[2], ...).
  void (*rhsRoutineSingleInstance_)(
      void *context, double t, double *states, double *rates,
      double *algebraics,
      double *
          parameters); //< function pointer to the rhs routine that can compute
                       //several instances of the problem in parallel. Data is
                       //assumed to contain values for a state contiguously,
                       //e.g. (state[1], state[1], state[1], state[2], state[2],
                       //state[2], ...). The first parameter is a this pointer.
  void (*initConstsOpenCOR_)(
      double *CONSTANTS, double *RATES,
      double *STATES); //< function pointer to a function that initializes
                       //everything, generated by OpenCOR
  void (*computeRatesOpenCOR_)(
      double VOI, double *CONSTANTS, double *RATES, double *STATES,
      double *ALGEBRAIC); //< function pointer to a function that computes the
                          //rates from the current states, generated by OpenCOR
  void (*computeVariablesOpenCOR_)(
      double VOI, double *CONSTANTS, double *RATES, double *STATES,
      double *ALGEBRAIC); //< function pointer to a function that computes the
                          //algebraic variables from the current states,
                          //generated by OpenCOR
};

#include "cellml/01_rhs_routine_handler.tpp"
