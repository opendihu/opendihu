#pragma once

#include <hdf5.h>

#include "control/types.h"
#include "control/dihu_context.h"

#include "input_reader/generic.h"
#include "input_reader/hdf5/utility.h"

namespace InputReader {
namespace HDF5 {
//! HDF5 File abstraction that automatically opens and closes a file, extracts
//! all written datasets and caches them as well as provide a simple interface
//! for accessing anything within the hdf5 file.
class File : public Generic {
public:
  //! constructor that opens a HDF5 file
  File(const char *file);
  //! closes the open file
  virtual ~File();

  //! Check if the file has a attribute on the root node with a given name
  bool hasAttribute(const char *name) const;
  //! Check if the file has a dataset with a given name
  bool hasDataset(const char *name) const;

  //! Read a i32 attribute on the root node and return it
  template <typename T,
            std::enable_if_t<std::is_same<T, int32_t>::value, bool> = true>
  bool readAttr(const char *name, T &v) const;
  //! Read a double attribute on the root node and return it
  template <typename T,
            std::enable_if_t<std::is_same<T, double>::value, bool> = true>
  bool readAttr(const char *name, T &v) const;
  //! Read a string attribute on the root node and return it
  template <typename T,
            std::enable_if_t<std::is_same<T, std::string>::value, bool> = true>
  bool readAttr(const char *name, T &v) const;

  //! Read a dataset of integers in a flat vector passed in using an out
  //! variable
  virtual bool readIntVector(const char *name, std::vector<int32_t> &out,
                             const std::string &groupName = "") const;
  //! Read a dataset of doubles in a flat vector passed in using an out variable
  virtual bool readDoubleVector(const char *name, std::vector<double> &out,
                                const std::string &groupName = "") const;

protected:
  //! construct the full dataset name based on a given suffix
  //! It returns a pointer of a dataset name within the `datasets_` vector, so
  //! this does not need to be deleted.
  const std::string *getFullDatasetName(const char *name,
                                        const std::string &groupName) const;

  hid_t fileID_; //< id of open file

  std::vector<Object> datasets_;      //< cached datasets
  std::vector<Attribute> attributes_; //< cached attributes

private:
  //! Helper function that reads an attribute with a given name into a output
  //! variable. The type is determined by reading it from the attribute metainfo
  //! found in the file.
  herr_t readAttribute(const char *name, void *out) const;
};
} // namespace HDF5
} // namespace InputReader

#include "input_reader/hdf5/file.tpp"
