#pragma once

#include <hdf5.h>

#include "control/types.h"
#include "control/dihu_context.h"

namespace InputReader {
namespace HDF5 {
//! Helper objects that are used to cache Datasets
struct Object {
  //! Construct a new Object with a given name and type
  Object(const char *name, H5O_type_t type);

  std::string name; //< name of the object
  H5O_type_t type;  //< type of the object
};

//! Helper objects that are used to cache Attributes
struct Attribute {
  //! Construct a new attribute with a given name, cset and size
  Attribute(const char *name, H5T_cset_t cset, hsize_t size);

  std::string name; //< name of the attribute
  H5T_cset_t cset;  //< cset of the attribute
  hsize_t size;     //< size of the attribute, useful if we have a string
};

//! Callback function that can be used to iterate over all attributes in a hdf5
//! file. `ctx` is expected to be a `std::vector<Attribute>` and the new
//! Attribute will be emplaced back into that vector
herr_t attr_iter(hid_t loc_id, const char *name, const H5A_info_t *info,
                 void *ctx);

//! Callback function that can be used to iterate over all objects in a hdf5
//! file. `ctx` is expected to be a `std::vector<Object>` and the new
//! Object will be emplaced back into that vector. This currently will only
//! append Datasets to the vector. Groups e.g. are ignored.
herr_t obj_iter(hid_t loc_id, const char *name, const H5O_info_t *info,
                void *ctx);

} // namespace HDF5
} // namespace InputReader
