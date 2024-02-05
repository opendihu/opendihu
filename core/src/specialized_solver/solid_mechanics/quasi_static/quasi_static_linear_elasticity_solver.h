#pragma once

#include <Python.h> // has to be the first included header
#include "time_stepping_scheme/02_time_stepping_scheme_ode.h"
#include "interfaces/runnable.h"
#include "data_management/finite_element_method/finite_elements.h"
#include "control/dihu_context.h"
#include "partition/rank_subset.h"
#include "equation/linear_elasticity.h"
#include "data_management/specialized_solver/quasi_static_linear_elasticity.h"
#include "spatial_discretization/finite_element_method/finite_element_method.h"
#include "slot_connection/slot_connector_data.h"

namespace TimeSteppingScheme {

/** A specialized solver for 3D linear elasticity, as quasi-static timestepping
 * scheme (a new static solution every timestep)
 */
template <typename FiniteElementMethod>
class QuasiStaticLinearElasticitySolver : public Runnable {
public:
  typedef typename FiniteElementMethod::FunctionSpace FunctionSpace;
  typedef typename Data::FiniteElements<
      FunctionSpace, 3, Equation::Static::LinearElasticityActiveStress>
      DataLinearElasticityType;
  typedef Data::QuasiStaticLinearElasticity<DataLinearElasticityType> Data;
  typedef FieldVariable::FieldVariable<FunctionSpace, 1> FieldVariableType;
  typedef typename Data::SlotConnectorDataType SlotConnectorDataType;

  typedef ::SpatialDiscretization::
      FiniteElementMethod< // FEM for initial potential flow, fiber directions
          typename FunctionSpace::Mesh, typename FunctionSpace::BasisFunction,
          Quadrature::Gauss<3>, Equation::Static::Laplace>
          FiniteElementMethodPotentialFlow;

  //! constructor
  QuasiStaticLinearElasticitySolver(DihuContext context);

  //! advance simulation by the given time span, data in solution is used,
  //! afterwards new data is in solution
  void advanceTimeSpan(bool withOutputWritersEnabled = true);

  //! initialize components of the simulation
  void initialize();

  //! dummy method, set endTime as current output time
  void setTimeSpan(double startTime, double endTime);

  //! run the simulation
  void run();

  //! reset state
  void reset();

  //! call the output writer on the data object, output files will contain
  //! currentTime, with callCountIncrement !=1 output timesteps can be skipped
  void callOutputWriter(int timeStepNo, double currentTime,
                        int callCountIncrement = 1);

  //! return the data object
  Data &data();

  //! get the data that will be transferred in the operator splitting to the
  //! other term of the splitting the transfer is done by the
  //! slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! output the given data for debugging
  std::string getString(std::shared_ptr<SlotConnectorDataType> data);

protected:
  // from the activation scalar field (symbol gamma), compute the active stress
  // tensor field in fiber direction
  void computeActiveStress();

  DihuContext context_; //< object that contains the python config for the
                        // current context and the global singletons meshManager
                        // and solverManager

  OutputWriter::Manager
      outputWriterManager_; //< manager object holding all output writer
  Data data_;               //< data object

  FiniteElementMethodPotentialFlow
      finiteElementMethodPotentialFlow_; //< the finite element object that is
                                         // used for the Laplace problem of the
                                         // potential flow, needed for the fiber
                                         // directions
  FiniteElementMethod
      finiteElementMethodLinearElasticity_; //< the finite element object that
                                            // solves the linear elasticity
                                            // equation

  std::string durationLogKey_; //< key with with the duration of the computation
                               // is written to the performance measurement log

  SpatialParameter<FunctionSpace, MathUtility::Matrix<3, 3, double>>
      anisotropyTensor_; //< the tensor that defines the anisotropy in the
                         // material model

  bool initialized_;              //< if this object was already initialized
  PythonConfig specificSettings_; //< python object containing the value of the
                                  // python config dict with corresponding key
  double endTime_;                //< end time of current time step
  double maximumActiveStress_; //< parameter value of the maximum active stress,
                               // this is the scaling factor of the activation
                               // value to get the active stress tensor
  double
      strainScalingCurveWidth_; //< width of a parabola that scales the stress
                                // dependend on the relative sarcomere length
  double scalingFactor_;        //< factor with which to scale the displacement
  Vec3 fiberDirection_; //< direction of fibers, needed for anisotropy of
                        // material
  bool usePotentialFlowForFiberDirection_; //< if a computation of a potential
                                           // flow should be used for
                                           // determining the fiber direction
};

} // namespace TimeSteppingScheme

#include "specialized_solver/solid_mechanics/quasi_static/quasi_static_linear_elasticity_solver.tpp"
