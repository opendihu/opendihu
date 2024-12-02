#include "input_reader/json/partial.h"

namespace InputReader {
namespace Json {
Partial::Partial(const char *file) : File(file) {
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize_);
  MPI_Comm_rank(MPI_COMM_WORLD, &ownRank_);
}

bool Partial::readIntVector(const char *name, std::vector<int32_t> &out,
                            const std::string &groupName) const {
  const std::string *path = getPathToDataset(name, groupName);
  if (!path) {
    return false;
  }
  std::vector<int32_t> chunks =
      content_[nlohmann::json::json_pointer(*path)][name]["__attributes"]
              ["chunkDims"]
                  .template get<std::vector<int32_t>>();
  std::vector<int32_t> data =
      content_[nlohmann::json::json_pointer(*path)][name]["__data"]
          .template get<std::vector<int32_t>>();
  int32_t blockSize = chunks[ownRank_];
  int32_t offset =
      std::accumulate(chunks.begin(), chunks.begin() + ownRank_, 0);

  out.resize(blockSize);
  std::copy(data.begin() + offset, data.begin() + offset + blockSize,
            out.begin());

  return true;
}

bool Partial::readDoubleVector(const char *name, std::vector<double> &out,
                               const std::string &groupName) const {
  const std::string *path = getPathToDataset(name, groupName);
  if (!path) {
    return false;
  }
  std::vector<int32_t> chunks =
      content_[nlohmann::json::json_pointer(*path)][name]["__attributes"]
              ["chunkDims"]
                  .template get<std::vector<int32_t>>();
  std::vector<double> data =
      content_[nlohmann::json::json_pointer(*path)][name]["__data"]
          .template get<std::vector<double>>();
  int32_t blockSize = chunks[ownRank_];
  int32_t offset =
      std::accumulate(chunks.begin(), chunks.begin() + ownRank_, 0);

  out.resize(blockSize);
  std::copy(data.begin() + offset, data.begin() + offset + blockSize,
            out.begin());
  return true;
}
} // namespace Json
} // namespace InputReader
