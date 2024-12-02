#include "input_reader/hdf5/utility.h"

namespace InputReader {
namespace HDF5 {
Object::Object(const char *name, H5O_type_t type) : name(name), type(type) {}

Attribute::Attribute(const char *name, H5T_cset_t cset, hsize_t size)
    : name(name), cset(cset), size(size) {}

herr_t attr_iter(hid_t loc_id, const char *name, const H5A_info_t *info,
                 void *ctx) {
  std::vector<Attribute> *attrs = (std::vector<Attribute> *)ctx;
  attrs->emplace_back(name, info->cset, info->data_size);

  return 0;
}

herr_t obj_iter(hid_t loc_id, const char *name, const H5O_info_t *info,
                void *ctx) {
  std::vector<Object> *datasets = (std::vector<Object> *)ctx;

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
} // namespace HDF5
} // namespace InputReader
