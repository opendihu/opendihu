#pragma once

#include <Python.h> // has to be the first included header

#include "interfaces/runnable.h"
#include "data_management/control/map_dofs.h"
#include "control/dihu_context.h"
#include "control/map_dofs/value_communicator.h"

namespace Control {

/** Object that defines an own number of field variable with given function
 * space. It can copy dofs between field variables, also reducing among
 * processes.
 *
 *  FunctionSpaceType is the function space of the additional field variables
 * that are created. Apart from this function space also field variables of
 * NestedSolverType that may have different function spaces are used for the
 * mapping.
 */
template <typename FunctionSpaceType, typename NestedSolverType>
class MapDofs : public Runnable {
public:
  //! make the FunctionSpace available
  typedef FunctionSpaceType FunctionSpace;

  //! define the type of the data object,
  typedef ::Data::MapDofs<FunctionSpaceType, NestedSolverType> Data;
  typedef typename Data::SlotConnectorDataType SlotConnectorDataType;

  //! constructor, gets the DihuContext object which contains all python
  //! settings
  MapDofs(DihuContext context);

  //! advance simulation by the given time span
  void advanceTimeSpan(bool withOutputWritersEnabled = true);

  //! initialize field variables and everything needed for the dofs mapping
  void initialize();

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary
  //! ("uninitialize")
  void reset();

  //! set the time span of the nested solver
  void setTimeSpan(double startTime, double endTime);

  //! start time of time interval to be simulated
  double startTime();

  //! end time of simulation
  double endTime();

  //! call the output writer on the data object, output files will contain
  //! currentTime, with callCountIncrement !=1 output timesteps can be skipped
  void callOutputWriter(int timeStepNo, double currentTime,
                        int callCountIncrement = 1);

  //! return the data object of the timestepping scheme
  Data &data();

  //! Get the data that will be transferred in the operator splitting or
  //! coupling to the other term of the splitting/coupling. the transfer is done
  //! by the slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

protected:
  /** all the data that is stored for a mapping, as in settings
   * "beforeComputation" or "afterComputation"
   */
  struct DofsMappingType {
    int connectorSlotNoFrom; //< slot no, i.e. a field variable from which to
                             // get the dofs
    std::vector<int> connectorSlotNosTo; //< slots nos, i.e. field variables to
                                         // which to write the dofs to
    bool dofNoIsGlobalFrom; //< if the keys in dofsMapping specify global nos
    bool dofNoIsGlobalTo;   //< if the values in dofsMapping specify global nos
    int slotConnectorArrayIndexFrom; //< array index if the connector slot
                                     // consists of a vector of multiple layers,
                                     // e.g. Multidomain with multiple
                                     // compartments or fibers with even two
                                     // nested MultipleInstances
    int slotConnectorArrayIndexTo;   //< array index if the connector slot
                                     // consists of a vector of multiple layers,
    // e.g. Multidomain with multiple compartments
    // or fibers with even two nested
    // MultipleInstances

    double thresholdValue; //< if mode is "localSetIfAboveThreshold", this is
                           // the threshold, if the value is above it, set the
                           // value `valueToSet`
    double valueToSet; //< if mode is "localSetIfAboveThreshold", this is the
                       // value to set the target dof to, if the source dof is
                       // above thresholdValue.

    enum {
      modeCopyLocal,
      modeCopyLocalIfPositive,
      modeLocalSetIfAboveThreshold,
      modeCallback,
      modeCommunicate
    } mode; // how to handle multiple dofs that map on a single dof

    std::map<int, std::vector<int>>
        dofsMapping; //< dofNoFrom : list of dofNosTo

    ValueCommunicator
        valueCommunicator; //< the object that performs the MPI communication
    std::map<int, std::vector<dof_no_t>>
        dofNosLocalOfValuesToSendToRanks; //< Only for modeCommunicate. For
                                          // every rank the local dof nos of the
                                          // values that will be sent to the
                                          // rank
    std::vector<dof_no_t>
        allDofNosToSetLocal; //< For all modes except modeCallback. For the
                             // received values the local dof nos where to store
                             // the values in the field variable
    std::vector<dof_no_t>
        dofNosToSetLocal; //< For all modes except modeCallback. Temporary
                          // variable, the local dofs where to set the values
                          // that were retrieved from the dofs.

    std::vector<dof_no_t>
        inputDofs; //< Only for modeCallback. The input dofs for the callback.
    std::vector<std::vector<dof_no_t>>
        outputDofs; //< Only for modeCallback. outputDofs[slotIndex], the output
                    // dofs for the callback.

    PyObject *callback; //< python callback to manipulate the input data
    PyObject *
        slotNosPy; //< list [fromSlotNo, toSlotNo, fromArrayIndex, toArrayIndex]
    PyObject *buffer;         //< user buffer of the callback function
    PyObject *outputValuesPy; //< python list of output values
  };

  //! execute the configured mapping of dofs, this gets called before and after
  //! advanceTimeSpan of the nested solver
  void performMappings(std::vector<DofsMappingType> &mappings,
                       double currentTime);

  //! parse the python settings, one of "beforeComputation", "afterComputation"
  //! (this is settingsKey) and store in mappings
  void parseMappingFromSettings(std::string settingsKey,
                                std::vector<DofsMappingType> &mappings);

  //! initialize the valueCommunicator and dof vectors in the given mappings
  void initializeCommunication(std::vector<DofsMappingType> &mappings);

  //! get the mesh partition of the function space of the field variable that
  //! corresponds to the given slot. Note, it is of type MeshPartitionBase,
  //! because the FunctionSpaceType of that slot is not fixed.
  std::shared_ptr<Partition::MeshPartitionBase>
  getMeshPartitionBase(int slotNo, int arrayIndex);

  //! set the values at given dofs at the field variable given by slotNo, return
  //! value: if values were set
  bool slotSetValues(int slotNo, int arrayIndex,
                     const std::vector<dof_no_t> &dofNosLocal,
                     const std::vector<double> &values,
                     InsertMode petscInsertMode = INSERT_VALUES);

  //! get the values at given dofs at the field variable given by slotNo
  void slotGetValues(int slotNo, int arrayIndex,
                     const std::vector<dof_no_t> &dofNosLocal,
                     std::vector<double> &values);

  //! call slotSetRepresentationGlobal on the field variable of the indicated
  //! slot
  void slotSetRepresentationGlobal(int slotNo, int arrayIndex);

  DihuContext context_; //< context object that gives access to global singleton
                        // and has python settings
  PythonConfig specificSettings_; //< the python settings for this object
  NestedSolverType nestedSolver_; //< the nested solver that will be called in
                                  // advanceTimeSpan

  Data data_; //< the data object that stores the additional field variables
  std::vector<std::string> slotNames_; //< names of all slots

  std::vector<DofsMappingType>
      mappingsBeforeComputation_; //< settings for all mappings that should be
                                  // done before the computation
  std::vector<DofsMappingType>
      mappingsAfterComputation_; //< settings for all mappings that should be
                                 // done after the computation
};

} // namespace Control

#include "control/map_dofs/map_dofs.tpp"
#include "control/map_dofs/map_dofs_initialize.tpp"
