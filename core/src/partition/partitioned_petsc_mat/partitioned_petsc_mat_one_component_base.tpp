#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_one_component_base.h"

template <typename RowsFunctionSpaceType, typename ColumnsFunctionSpaceType>
PartitionedPetscMatOneComponentBase<RowsFunctionSpaceType,
                                    ColumnsFunctionSpaceType>::
    PartitionedPetscMatOneComponentBase(
        std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>>
            meshPartitionRows,
        std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>>
            meshPartitionColumns,
        std::string name)
    : meshPartitionRows_(meshPartitionRows),
      meshPartitionColumns_(meshPartitionColumns), name_(name) {}

template <typename RowsFunctionSpaceType, typename ColumnsFunctionSpaceType>
std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>>
PartitionedPetscMatOneComponentBase<
    RowsFunctionSpaceType, ColumnsFunctionSpaceType>::meshPartitionRows() {
  return meshPartitionRows_;
}

template <typename RowsFunctionSpaceType, typename ColumnsFunctionSpaceType>
std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>>
PartitionedPetscMatOneComponentBase<
    RowsFunctionSpaceType, ColumnsFunctionSpaceType>::meshPartitionColumns() {
  return meshPartitionColumns_;
}
