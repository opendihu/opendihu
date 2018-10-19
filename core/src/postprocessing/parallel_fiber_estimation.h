#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>

#include "discretizable_in_time/discretizable_in_time.h"
#include "control/runnable.h"
#include "data_management/parallel_fiber_estimation.h"
#include "quadrature/gauss.h"

namespace Postprocessing
{

/** A class that creates a parallel mesh and at the same time traces streamlines through a Δu=0 field.
 */
template<typename BasisFunctionType>
class ParallelFiberEstimation : public Runnable
{
public:
  //! constructor
  ParallelFiberEstimation(DihuContext context);

  //! initialize
  void initialize();

  //! run tracing of stream lines
  void run();

  //! function space to use, i.e. 3D structured deformable grid
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType> FunctionSpace;
  typedef SpatialDiscretization::FiniteElementMethod<
    Mesh::StructuredDeformableOfDimension<3>,
    BasisFunctionType,
    Quadrature::Gauss<3>,
    Equation::Static::Laplace
  > FiniteElementMethodType;

protected:

  //! perform the algorithm to recursively, collectively refine the mesh and trace fibers for the boundaries of the subdomains
  void generateParallelMesh();

  const DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  FiniteElementMethodType problem_;   ///< the DiscretizableInTime object that is managed by this class

  Data::ParallelFiberEstimation<FunctionSpace> data_;    ///< the data object that holds the gradient field variable

  std::shared_ptr<FunctionSpace> functionSpace_;   ///< current function space / mesh

  PyObject *specificSettings_;   ///< the specific python config for this module
  std::vector<Vec3> seedPositions_;  ///< the seed points from where the streamlines start

  std::string stlFilename_;   ///< the filename of the STL file
  int bottomZClip_;   ///< bottom z-value of the volume to consider
  int topZClip_;   ///< top z-value of the volume to consider
  int nElementsZPerSubdomain_;   ///< number of elements per subdomain in z direction

  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
};

};  // namespace

#include "postprocessing/parallel_fiber_estimation.tpp"
