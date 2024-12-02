#pragma once

#include "control/types.h"
#include "control/dihu_context.h"

#include "input_reader/generic.h"
#include "input_reader/json/utility.h"

#include <nlohmann/json.hpp>

namespace InputReader {
namespace Json {
//! Json file abstraction that automatically opens and closes a file, extracts
//! all written datasets and caches them as well as provide a simple interface
//! for accessing anything within the json file.
class File : public Generic {
public:
  //! constructor that opens a Json file
  File(const char *file);
  //! closes the open file
  virtual ~File() = default;

  //! Check if the file has a attribute on the root node with a given name
  bool hasAttribute(const char *name) const;
  //! Check if the file has a dataset with a given name
  bool hasDataset(const char *name) const;

  //! Read a attribute on the root node and return it
  template <typename T> bool readAttr(const char *name, T &out) const;

  //! Read a dataset of integers in a flat vector passed in using an out
  //! variable
  virtual bool readIntVector(const char *name, std::vector<int32_t> &out,
                             const std::string &groupName = "") const;
  //! Read a dataset of doubles in a flat vector passed in using an out variable
  virtual bool readDoubleVector(const char *name, std::vector<double> &out,
                                const std::string &groupName = "") const;

protected:
  //! Returns a jsonpointer path to the dataset
  const std::string *getPathToDataset(const char *name,
                                      const std::string &groupName) const;

  std::vector<Object> datasets_;      //< cached datasets
  std::vector<Attribute> attributes_; //< cached attributes

  nlohmann::json content_; //< Parsed content of the json file
};
} // namespace Json
} // namespace InputReader

#include "input_reader/json/file.tpp"
