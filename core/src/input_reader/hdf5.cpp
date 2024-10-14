#include "input_reader/hdf5.h"

static herr_t attr_iter(hid_t loc_id, const char *name, const H5A_info_t *info,
                        void *ctx) {
  std::vector<InputReader::Attribute> *attrs =
      (std::vector<InputReader::Attribute> *)ctx;
  attrs->emplace_back(name, info->cset, info->data_size);

  return 0;
}

static herr_t obj_iter(hid_t loc_id, const char *name, const H5O_info_t *info,
                       void *ctx) {
  std::vector<InputReader::Object> *datasets =
      (std::vector<InputReader::Object> *)ctx;

  switch (info->type) {
  case H5O_TYPE_GROUP:
    break;
  case H5O_TYPE_DATASET:
    datasets->emplace_back(name, info->type);
    break;
  case H5O_TYPE_NAMED_DATATYPE:
    LOG(WARNING) << "Unsed HDF5 type found: Datatype | " << name;
    break;
  default:
    LOG(WARNING) << "Unknown HDF5 type found: " << name;
    break;
  }

  return 0;
}

namespace InputReader {
Object::Object(const char *name, H5O_type_t type) : name(name), type(type) {}

Attribute::Attribute(const char *name, H5T_cset_t cset, hsize_t size)
    : name(name), cset(cset), size(size) {}

HDF5::HDF5(const char *file) {
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

HDF5::~HDF5() {
  herr_t err = H5Fclose(fileID_);
  assert(err >= 0);
  (void)err;
}

bool HDF5::hasAttribute(const char *name) const {
  for (const auto &e : attributes_) {
    if (e.name == name) {
      return true;
    }
  }
  return false;
}

bool HDF5::hasDataset(const char *name) const {
  for (const auto &e : datasets_) {
    std::string sp = e.name.substr(e.name.find("/") + 1);
    if (sp == name) {
      return true;
    }
  }
  return false;
}

int32_t HDF5::readInt32Attribute(const char *name) const {
  hid_t attr = H5Aopen_name(fileID_, name);
  assert(attr >= 0);

  int32_t out;
  herr_t err = H5Aread(attr, H5T_NATIVE_INT, &out);
  assert(err >= 0);
  err = H5Aclose(attr);
  assert(err >= 0);
  return out;
}

double HDF5::readDoubleAttribute(const char *name) const {
  hid_t attr = H5Aopen_name(fileID_, name);
  assert(attr >= 0);

  double out;
  herr_t err = H5Aread(attr, H5T_NATIVE_DOUBLE, &out);
  assert(err >= 0);
  err = H5Aclose(attr);
  assert(err >= 0);
  return out;
}

std::string HDF5::readStringAttribute(const char *name) const {
  hid_t attr = H5Aopen_name(fileID_, name);
  assert(attr >= 0);

  std::string out;
  hid_t atype = H5Tcopy(H5T_C_S1);
  for (const auto &e : attributes_) {
    if (e.name == name) {
      H5Tset_size(atype, e.size + 1);
      out.resize(e.size + 1);
      break;
    }
  }

  herr_t err = H5Aread(attr, atype, (void *)out.c_str());
  assert(err >= 0);
  err = H5Tclose(atype);
  assert(err >= 0);
  err = H5Aclose(attr);
  assert(err >= 0);
  return out;
}

bool HDF5::readIntVector(const char *name, std::vector<int32_t> &out) const {
  const std::string *fullName = getFullDatasetName(name);
  if (!fullName) {
    return false;
  }
  herr_t err;
  hid_t dset = H5Dopen(fileID_, fullName->c_str(), H5P_DEFAULT);
  assert(dset >= 0);

  hid_t dspace = H5Dget_space(dset);
  assert(dspace >= 0);
  const int ndims = H5Sget_simple_extent_ndims(dspace);
  if (ndims != 1) {
    err = H5Sclose(dspace);
    assert(err >= 0);
    err = H5Dclose(dset);
    assert(err >= 0);
    return false;
  }
  std::array<hsize_t, 1> dims = {0};
  err = H5Sget_simple_extent_dims(dspace, dims.data(), nullptr);
  assert(err >= 0);
  err = H5Sclose(dspace);
  assert(err >= 0);

  out.resize(dims[0]);
  err =
      H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
  assert(err >= 0);
  err = H5Dclose(dset);
  assert(err >= 0);

  return true;
}

bool HDF5::readDoubleVector(const char *name, std::vector<double> &out) const {
  const std::string *fullName = getFullDatasetName(name);
  if (!fullName) {
    return false;
  }
  herr_t err;
  hid_t dset = H5Dopen(fileID_, fullName->c_str(), H5P_DEFAULT);
  assert(dset >= 0);

  hid_t dspace = H5Dget_space(dset);
  assert(dspace >= 0);
  const int ndims = H5Sget_simple_extent_ndims(dspace);
  if (ndims != 1) {
    err = H5Sclose(dspace);
    assert(err >= 0);
    err = H5Dclose(dset);
    assert(err >= 0);
    return false;
  }
  std::array<hsize_t, 1> dims = {0};
  err = H5Sget_simple_extent_dims(dspace, dims.data(), nullptr);
  assert(err >= 0);
  err = H5Sclose(dspace);
  assert(err >= 0);

  out.resize(dims[0]);
  err = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                out.data());
  assert(err >= 0);
  err = H5Dclose(dset);
  assert(err >= 0);

  return true;
}

const std::string *HDF5::getFullDatasetName(const char *name) const {
  for (const auto &e : datasets_) {
    std::string sp = e.name.substr(e.name.find("/") + 1);
    if (sp == name) {
      return &e.name;
    }
  }
  return nullptr;
}
} // namespace InputReader
