#include "input_reader/json/file.h"

#include <sstream>

namespace InputReader {
namespace Json {
static void extractObjectsInto(const nlohmann::json &object,
                               std::vector<Object> &datasets,
                               const std::string &path) {
  for (const auto &el : object.items()) {
    if (el.key().rfind("__", 0)) {
      continue;
    }
    if (el.value().contains("__data")) {
      datasets.emplace_back(el.key(), path);
    } else {
      std::stringstream ss;
      ss << path << "/" << el.key();
      extractObjectsInto(el.value(), datasets, ss.str());
    }
  }
}

File::File(const char *file) : Generic() {
  std::ifstream ifs(file);
  content_ = nlohmann::json::parse(ifs);

  for (auto &el : content_["__attributes"].items()) {
    attributes_.emplace_back(el.key());
  }

  extractObjectsInto(content_, datasets_, "");
}

bool File::hasAttribute(const char *name) const {
  for (const auto &e : attributes_) {
    if (e.name == name) {
      return true;
    }
  }
  return false;
}

bool File::hasDataset(const char *name) const {
  for (const auto &e : datasets_) {
    if (e.name == name) {
      return true;
    }
  }
  return false;
}

bool File::readIntVector(const char *name, std::vector<int32_t> &out,
                         const std::string &groupName) const {
  return false;
}

bool File::readDoubleVector(const char *name, std::vector<double> &out,
                            const std::string &groupName) const {
  return false;
}

const std::string *File::getPathToDataset(const char *name,
                                          const std::string &groupName) const {
  std::string dsname = name;
  std::replace(dsname.begin(), dsname.end(), '/', '|');

  for (const auto &e : datasets_) {
    std::string sp = e.name.substr(e.name.find_last_of("/") + 1);
    if (sp == dsname) {
      if (groupName == "" || e.path.find(groupName) != std::string::npos) {
        // Return a pointer into the attributes set, this does not need to be
        // deleted
        return &e.path;
      }
    }
  }
  return nullptr;
}
} // namespace Json
} // namespace InputReader
