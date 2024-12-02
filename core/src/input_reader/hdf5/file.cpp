#include "input_reader/hdf5/file.h"

namespace InputReader {
namespace HDF5 {
File::File(const char *file) : Generic() {
  herr_t err;
  fileID_ = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(fileID_ >= 0);

  err = H5Aiterate(fileID_, H5_INDEX_NAME, H5_ITER_INC, nullptr, attr_iter,
                   (void *)&attributes_);
  assert(err >= 0);
  err = H5Ovisit(fileID_, H5_INDEX_NAME, H5_ITER_NATIVE, obj_iter,
                 (void *)&datasets_, H5O_INFO_BASIC);
  assert(err >= 0);

  for (const auto &e : attributes_) {
    LOG(DEBUG) << "InputReader - Attribute: " << e.name << " - " << e.size;
  }
  for (const auto &e : datasets_) {
    LOG(DEBUG) << "InputReader - Dataset: " << e.name << " - " << e.type;
  }
}

File::~File() {
  herr_t err = H5Fclose(fileID_);
  assert(err >= 0);
  (void)err;
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
    std::string sp = e.name.substr(e.name.find("/") + 1);
    if (sp == name) {
      return true;
    }
  }
  return false;
}

herr_t File::readAttribute(const char *name, void *out) const {
  hid_t attr = H5Aopen_name(fileID_, name);
  if (attr < 0) {
    return attr;
  }
  hid_t type = H5Aget_type(attr);
  if (type < 0) {
    H5Aclose(attr);
    return type;
  }
  herr_t err = H5Aread(attr, type, out);
  if (err < 0) {
    H5Aclose(attr);
    return err;
  }
  err = H5Aclose(attr);
  return err;
}

bool File::readIntVector(const char *name, std::vector<int32_t> &out,
                         const std::string &groupName) const {
  return false;
}

bool File::readDoubleVector(const char *name, std::vector<double> &out,
                            const std::string &groupName) const {
  return false;
}

const std::string *
File::getFullDatasetName(const char *name, const std::string &groupName) const {
  // same as output writer, sanatize name so it doesn't contain any `/`
  std::string dsname = name;
  std::replace(dsname.begin(), dsname.end(), '/', '|');

  for (const auto &e : datasets_) {
    std::string sp = e.name.substr(e.name.find_last_of("/") + 1);
    if (sp == dsname) {
      if (groupName == "" || e.name.find(groupName) != std::string::npos) {
        // Return a pointer into the attributes set, this does not need to be
        // deleted
        return &e.name;
      }
    }
  }
  return nullptr;
}
} // namespace HDF5
} // namespace InputReader
