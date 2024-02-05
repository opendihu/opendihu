#pragma once

#include "control/dihu_context.h"
#include "data_management/time_stepping/time_stepping.h"
#include "interfaces/discretizable_in_time.h"
#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "data_management/data.h"
#include "data_management/time_stepping/time_stepping.h"
#include "cellml/03_cellml_adapter.h"

namespace TimeSteppingScheme {

/**
 * Normal timestepping scheme to solve ODES, like Heun
 */
template <typename FunctionSpaceType, int nComponents>
class TimeSteppingSchemeOdeBase : public TimeSteppingScheme {
public:
  typedef FunctionSpaceType FunctionSpace;
  typedef typename Data::TimeStepping<FunctionSpaceType, nComponents>
      Data; // type of Data object
  typedef typename Data::SlotConnectorDataType SlotConnectorDataType;

  //! constructor
  TimeSteppingSchemeOdeBase(DihuContext context, std::string name);

  //! destructor
  virtual ~TimeSteppingSchemeOdeBase() {}

  //! run simulation
  virtual void run();

  //! return the data object
  Data &data();

  //! output the given data for debugging
  // virtual std::string getString(std::shared_ptr<SlotConnectorDataType> data);

  //! Get the data that will be transferred in the operator splitting to the
  //! other term of the splitting. The transfer is done by the
  //! slot_connector_data_transfer class. The data is passed on from the
  //! DiscretizableInTime object.
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! initialize discretizableInTime
  virtual void initialize();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);

  //! reset state such that new initialization becomes necessary
  virtual void reset();

  //! call the output writer on the data object, output files will contain
  //! currentTime, with callCountIncrement !=1 output timesteps can be skipped
  virtual void callOutputWriter(int timeStepNo, double currentTime,
                                int callCountIncrement = 1);

protected:
  //! read initial values from settings and set field accordingly
  void setInitialValues();

  //! check if the current solution contains nan or inf values, if so, output an
  //! error
  void checkForNanInf(int timeStepNo, double currentTime);

  //! prepare the discretizableInTime object for the following call to
  //! getSlotConnectorData()
  virtual void prepareForGetSlotConnectorData() = 0;

  std::shared_ptr<Data>
      data_; //< data object that holds all PETSc vectors and matrices

  bool initialized_; //< if initialize() was already called
  std::string name_; //< the name given to this time stepping scheme
  double prefactor_; //< a factor with which the result is multiplied when the
                     //data is used in a splitting scheme
  bool checkForNanInf_; //< if the solution should be checked for nan and inf
                        //values
};

} // namespace TimeSteppingScheme

#include "time_stepping_scheme/01_time_stepping_scheme_ode_base.tpp"
