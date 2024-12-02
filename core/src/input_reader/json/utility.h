#pragma once

#include <hdf5.h>

#include "control/types.h"
#include "control/dihu_context.h"

namespace InputReader {
namespace Json {
//! Helper objects that are used to cache Datasets
struct Object {
  //! Construct a new Object with a given name and type
  Object(const std::string &name, const std::string &path);

  std::string name; //< name of the object
  std::string path; //< json path as json pointer to reach this object
};

//! Helper objects that are used to cache Attributes
struct Attribute {
  //! Construct a new attribute with a given name, cset and size
  Attribute(const std::string &name);

  std::string name; //< name of the attribute
};
} // namespace Json
} // namespace InputReader
