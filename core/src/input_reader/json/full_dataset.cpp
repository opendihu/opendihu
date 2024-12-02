#include "input_reader/json/full_dataset.h"

namespace InputReader {
namespace Json {
FullDataset::FullDataset(const char *file) : File(file) {}

bool FullDataset::readIntVector(const char *name, std::vector<int32_t> &out,
                                const std::string &groupName) const {
  const std::string *path = getPathToDataset(name, groupName);
  if (!path) {
    return false;
  }
  content_[nlohmann::json::json_pointer(*path)][name]["__data"].get_to(out);
  return true;
}

bool FullDataset::readDoubleVector(const char *name, std::vector<double> &out,
                                   const std::string &groupName) const {
  const std::string *path = getPathToDataset(name, groupName);
  if (!path) {
    return false;
  }
  content_[nlohmann::json::json_pointer(*path)][name]["__data"].get_to(out);
  return true;
}
} // namespace Json
} // namespace InputReader
