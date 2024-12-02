#include "checkpointing/json/independent.h"

namespace Checkpointing {
namespace Json {
Independent::Independent(DihuContext context,
                         std::shared_ptr<Partition::RankSubset> rankSubset,
                         const std::string &prefix)
    : Generic(prefix), writer_(std::make_unique<OutputWriter::Json>(
                           context, PythonConfig(nullptr), rankSubset)) {
  writer_->setCombineFiles(false);
  writer_->setWriteMeta(true);
  writer_->setUseCheckpointData(true);
}
} // namespace Json
} // namespace Checkpointing
