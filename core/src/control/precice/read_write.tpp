#include "control/precice/read_write.h"

#include <sstream>

#include "spatial_discretization/neumann_boundary_conditions/00_neumann_boundary_conditions_base.h"
#include "utility/vector_operators.h"

namespace Control {

#ifdef HAVE_PRECICE
template <typename NestedSolver>
void PreciceAdapterReadWrite<NestedSolver>::preciceReadData() {}

template <typename NestedSolver>
void PreciceAdapterReadWrite<NestedSolver>::setDirichletBoundaryConditions(
    typename PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData &preciceData) {}

template <typename NestedSolver>
void PreciceAdapterReadWrite<NestedSolver>::setNeumannBoundaryConditions(
    typename PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData &preciceData) {}

template <typename NestedSolver>
void PreciceAdapterReadWrite<NestedSolver>::preciceWriteData() {}

#endif

} // namespace Control
