#pragma once

#include <hdf5.h>

#include "control/types.h"
#include "control/dihu_context.h"

#include "input_reader/hdf5/file.h"
#include "input_reader/hdf5/utility.h"

namespace InputReader {
namespace HDF5 {
//! HDF5 File abstraction that automatically opens and closes a file, extracts
//! all written datasets and caches them as well as provide a simple interface
//! for accessing anything within the hdf5 file.
//!
//! This interface will always read the full dataset into the output
//! variablesvariables.
class FullDataset : public File {
public:
  //! constructor that opens a HDF5 file
  FullDataset(const char *file);
  //! closes the open file
  virtual ~FullDataset() = default;

  //! Read a dataset of integers in a flat vector passed in using an out
  //! variable
  bool readIntVector(const char *name, std::vector<int32_t> &out,
                     const std::string &groupName = "") const;
  //! Read a dataset of doubles in a flat vector passed in using an out variable
  bool readDoubleVector(const char *name, std::vector<double> &out,
                        const std::string &groupName = "") const;

private:
  //! Helper function for reading any type into a flat vector witha a given
  //! name. Function is here to reduce code duplication.
  template <typename T>
  bool readDataset(const char *name, const std::string &groupName,
                   hid_t memTypeId, std::vector<T> &out) const;
};
} // namespace HDF5
} // namespace InputReader

#include "input_reader/hdf5/full_dataset.tpp"
