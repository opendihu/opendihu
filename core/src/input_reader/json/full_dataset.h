#pragma once

#include "control/types.h"
#include "control/dihu_context.h"

#include "input_reader/json/file.h"
#include "input_reader/json/utility.h"

namespace InputReader {
namespace Json {
//! Json File abstraction that automatically opens and closes a file, extracts
//! all written datasets and caches them as well as provide a simple interface
//! for accessing anything within the hdf5 file.
//!
//! This interface will always read the full dataset into the output
//! variablesvariables.
class FullDataset : public File {
public:
  //! constructor that opens a Json file
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
};
} // namespace Json
} // namespace InputReader
