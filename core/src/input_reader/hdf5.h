#pragma once

#include <hdf5.h>

#include "control/types.h"
#include "control/dihu_context.h"

namespace InputReader {
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

//! HDF5 File abstraction that automatically opens and closes a file, extracts
//! all written datasets and caches them as well as provide a simple interface
//! for accessing anything within the hdf5 file.
class HDF5 {
public:
  //! constructor that opens a HDF5 file
  HDF5(const char *file);
  //! closes the open file
  ~HDF5();

  //! Check if the file has a attribute on the root node with a given name
  bool hasAttribute(const char *name) const;
  //! Check if the file has a dataset with a given name
  bool hasDataset(const char *name) const;

  //! Read a i32 attribute on the root node and return it
  int32_t readInt32Attribute(const char *name) const;
  //! Read a double attribute on the root node and return it
  double readDoubleAttribute(const char *name) const;
  //! Read a string attribute on the root node and return it
  std::string readStringAttribute(const char *name) const;

  //! Read a dataset of integers in a flat vector passed in using an out
  //! variable
  bool readIntVector(const char *name, std::vector<int32_t> &out) const;
  //! Read a dataset of doubles in a flat vector passed in using an out variable
  bool readDoubleVector(const char *name, std::vector<double> &out) const;

  //! Read a dataset of doubles in a nested vector given a specific amount of
  //! components
  template <int nComponents>
  bool readNestedDoubleVector(
      const char *name,
      std::vector<std::array<double, nComponents>> &out) const;

private:
  //! construct the full dataset name based on a given suffix
  const std::string *getFullDatasetName(const char *name) const;

  hid_t fileID_; //< id of open file

  std::vector<Object> datasets_;      //< cached datasets
  std::vector<Attribute> attributes_; //< cached attributes
};
} // namespace InputReader

#include "input_reader/hdf5.tpp"
