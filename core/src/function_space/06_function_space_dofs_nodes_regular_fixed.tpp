#include "function_space/06_function_space_dofs_nodes.h"

#include <Python.h> // has to be the first included header
#include "easylogging++.h"
#include <array>

#include "utility/python_utility.h"
#include "control/types.h"
#include "utility/vector_operators.h"
//#include "field_variable/field_variable.h"
//#include "field_variable/00_field_variable_base.h"

namespace FunctionSpace {

template <int D, typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,
                       BasisFunctionType>::
    FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager,
                           PythonConfig specificSettings)
    : FunctionSpaceDofsNodesStructured<
          Mesh::StructuredRegularFixedOfDimension<D>, BasisFunctionType>::
          FunctionSpaceDofsNodesStructured(partitionManager, specificSettings),
      physicalExtent_({-1.0}) {
  this->meshWidth_ = 0;
  // meshWidth will be initialized in initialize
}

template <int D, typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,
                       BasisFunctionType>::
    FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager,
                           std::vector<double> &null,
                           PythonConfig specificSettings)
    : FunctionSpaceDofsNodes<
          Mesh::StructuredRegularFixedOfDimension<D>,
          BasisFunctionType>::FunctionSpaceDofsNodes(partitionManager,
                                                     specificSettings) {}

template <int D, typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,
                       BasisFunctionType>::
    FunctionSpaceDofsNodes(
        std::shared_ptr<Partition::Manager> partitionManager,
        std::vector<double> &null, std::array<element_no_t, D> nElements,
        std::array<double, D> physicalExtent,
        const std::array<int, D> nRanksPerCoordinateDirection,
        bool inputMeshIsGlobal)
    : FunctionSpaceDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,
                             BasisFunctionType>::
          FunctionSpaceDofsNodes(partitionManager, nElements, physicalExtent,
                                 nRanksPerCoordinateDirection,
                                 inputMeshIsGlobal) {}

template <int D, typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,
                       BasisFunctionType>::
    FunctionSpaceDofsNodes(
        std::shared_ptr<Partition::Manager> partitionManager,
        std::array<element_no_t, D> nElements,
        std::array<double, D> physicalExtent,
        const std::array<int, D> nRanksPerCoordinateDirection,
        bool inputMeshIsGlobal)
    : FunctionSpaceDofsNodesStructured<
          Mesh::StructuredRegularFixedOfDimension<D>,
          BasisFunctionType>::FunctionSpaceDofsNodesStructured(partitionManager,
                                                               nullptr),
      physicalExtent_(physicalExtent) {
  // compute mesh width from physical extent and number of elements in the
  // coordinate directions note for quadratic elements the mesh width is the
  // distance between the nodes, not length of elements
  this->meshWidth_ = 0;
  this->nRanks_ = nRanksPerCoordinateDirection;

  if (inputMeshIsGlobal) {
    std::copy(nElements.begin(), nElements.end(),
              this->nElementsPerCoordinateDirectionGlobal_.begin());
  } else {
    std::copy(nElements.begin(), nElements.end(),
              this->nElementsPerCoordinateDirectionLocal_.begin());

    this->forcePartitioningCreationFromLocalNumberOfElements_ =
        true; // this is defined in 03_function_space_partition.h

    // forcePartitioningCreationFromLocalNumberOfElements_ is set to true, this
    // means that the partitioning is created considering
    // this->nElementsPerCoordinateDirectionLocal_ and not depending on values
    // of inputMeshIsGlobal
  }
}

template <int D, typename BasisFunctionType>
void FunctionSpaceDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,
                            BasisFunctionType>::computeMeshWidth() {
  // if physicalExtent_ is not yet set
  if (physicalExtent_[0] < 0.0) {
    // only get physicalExtent if it is not a 1-node mesh with 0 elements
    if (D > 1 || this->nElementsPerCoordinateDirectionLocal_[0] != 0 ||
        this->nElementsPerCoordinateDirectionGlobal_[0] != 0) {
      if (this->specificSettings_.hasKey("physicalExtend")) {
        LOG(ERROR)
            << "You misspelled \"physicalExtent\" as \"physicalExtend\"!";
      }

      physicalExtent_ =
          this->specificSettings_.template getOptionArray<double, D>(
              "physicalExtent", 1.0, PythonUtility::NonNegative);
    } else {
      physicalExtent_[0] = 0.0; // set to 0 as a sign that this function space
                                // has no real spatial mesh meaning
    }
    VLOG(1) << "get physicalExtent: " << physicalExtent_;
  }

  std::array<element_no_t, D> nElements;
  bool inputMeshIsGlobal =
      this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
  if (inputMeshIsGlobal) {
    for (int coordinateIndex = 0; coordinateIndex < D; coordinateIndex++) {
      nElements[coordinateIndex] =
          this->meshPartition_->nElementsGlobal(coordinateIndex);
    }
  } else {
    for (int coordinateIndex = 0; coordinateIndex < D; coordinateIndex++) {
      nElements[coordinateIndex] =
          this->meshPartition_->nElementsLocal(coordinateIndex);
    }
  }

  // compute mesh width from physical extent and number of elements in the
  // coordinate directions note for quadratic elements the mesh width is the
  // distance between the nodes, not length of elements
  if (D > 1 || nElements[0] != 0) {
    for (int coordinateIndex = 0; coordinateIndex != D; coordinateIndex++) {
      double meshWidthCurrentDirection =
          physicalExtent_[coordinateIndex] /
          double(nElements[coordinateIndex] *
                 FunctionSpaceBaseDim<
                     1, BasisFunctionType>::averageNNodesPerElement());

      if (this->meshWidth_ == 0) {
        this->meshWidth_ = meshWidthCurrentDirection;
      } else if (fabs(this->meshWidth_ - meshWidthCurrentDirection) > 1e-14) {
        LOG(ERROR) << "Mesh has no uniform mesh width, use a "
                      "StructuredDeformableOfDimension<"
                   << D << "> mesh instead.";
        LOG(ERROR) << "mesh width: " << this->meshWidth_
                   << ", other: " << physicalExtent_[coordinateIndex] << "/"
                   << double(
                          nElements[coordinateIndex] *
                          FunctionSpaceBaseDim<
                              1, BasisFunctionType>::averageNNodesPerElement())
                   << "=" << meshWidthCurrentDirection;
      }
    }
  } else {
    // 1D 1-node mesh
    this->meshWidth_ = 1.0;
  }

  LOG(DEBUG) << "  FunctionSpaceDofsNodes Mesh::RegularFixed constructor, D="
             << D << ", nElements: " << nElements;
  LOG(DEBUG) << "  physicalExtent: " << physicalExtent_;
  LOG(DEBUG) << "  meshWidth: " << this->meshWidth_;
}

template <int D, typename BasisFunctionType>
void FunctionSpaceDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,
                            BasisFunctionType>::initialize() {
  // create meshPartition and redistribute elements if necessary, this needs
  // information about mesh size
  FunctionSpacePartition<Mesh::StructuredRegularFixedOfDimension<D>,
                         BasisFunctionType>::initialize();

  // initialize meshWidth from settings (physicalExtent) and number of elements
  computeMeshWidth();

  // pass a shared "this" pointer to the geometryField
  if (this->noGeometryField_)
    return;

  // retrieve "this" pointer and convert to downwards pointer of most derived
  // class "FunctionSpace"
  std::shared_ptr<FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,
                                BasisFunctionType>>
      thisMesh = std::static_pointer_cast<FunctionSpace<
          Mesh::StructuredRegularFixedOfDimension<D>, BasisFunctionType>>(
          this->shared_from_this());

  assert(thisMesh != nullptr);

  // create empty field variable for geometry field
  std::vector<std::string> componentNames{"x", "y", "z"};
  this->geometryField_ = std::make_shared<GeometryFieldType>(
      thisMesh, "geometry", componentNames, true);

  // no need to set values of the geometry field, because there is no data
  // explicitly stored the derivative values are appropriately returned by the
  // field variable

  // set initalized_ to true which indicates that initialize has been called
  this->initialized_ = true;
}

template <int D, typename BasisFunctionType>
double FunctionSpaceDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,
                              BasisFunctionType>::meshWidth() const {
  return this->meshWidth_;
}

} // namespace FunctionSpace
