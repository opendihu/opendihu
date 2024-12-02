#include "checkpointing/hdf5/independent.h"

namespace Checkpointing {
namespace HDF5 {
Independent::Independent(DihuContext context,
                         std::shared_ptr<Partition::RankSubset> rankSubset,
                         const std::string &prefix)
    : Generic(prefix), writer_(std::make_unique<OutputWriter::HDF5>(
                           context, PythonConfig(nullptr), rankSubset)) {
  writer_->setCombineFiles(false);
  writer_->setWriteMeta(true);
  writer_->setUseCheckpointData(true);
}
} // namespace HDF5
} // namespace Checkpointing
