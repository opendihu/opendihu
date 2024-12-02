#include "output_writer/json/json.h"

#include <nlohmann/json.hpp>

namespace OutputWriter {
Json::Json(DihuContext context, PythonConfig settings,
           std::shared_ptr<Partition::RankSubset> rankSubset)
    : Generic(context, settings, rankSubset) {
  combineFiles_ = settings.getOptionBool("combineFiles", false);
  writeMeta_ = settings.getOptionBool("writeMeta", true);
  useCheckpointData_ = false;
}

//! constructor, initialize nPoints and nCells to 0
Json::Piece::Piece() {
  properties.nPointsLocal = 0;
  properties.nCellsLocal = 0;
  properties.nPointsGlobal = 0;
  properties.nCellsGlobal = 0;
  properties.dimensionality = 0;
}

void Json::setCombineFiles(bool v) { combineFiles_ = v; }
void Json::setWriteMeta(bool v) { writeMeta_ = v; }
void Json::setUseCheckpointData(bool v) { useCheckpointData_ = v; }

//! assign the correct values to firstScalarName and firstVectorName, only if
//! properties has been set
void Json::Piece::setVTKValues() {
  // set values for firstScalarName and firstVectorName from the values in
  // pointDataArrays
  for (auto pointDataArray : properties.pointDataArrays) {
    if (firstScalarName == "" && pointDataArray.nComponents != 3)
      firstScalarName = pointDataArray.name;
    if (firstVectorName == "" && pointDataArray.nComponents == 3)
      firstVectorName = pointDataArray.name;
  }
}

void Json::writeCombinedTypesVector(JsonUtils::Group &group, uint64_t nValues,
                                    bool output3DMeshes, const char *dsname) {
  if (output3DMeshes) {
    std::vector<int32_t> values(nValues, 12);
    return group.writeSimpleVec<int32_t>(values, dsname);
  } else {
    std::vector<int32_t> values(nValues, 9);
    return group.writeSimpleVec<int32_t>(values, dsname);
  }
}

namespace JsonUtils {
File::File(const char *filename, bool mpiio)
    : filename_(filename), mpiio_(mpiio), content_(),
      ownRank_(DihuContext::ownRankNoCommWorld()),
      worldSize_(DihuContext::nRanksCommWorld()) {}

File::~File() {
  if (mpiio_) {
    // TODO: check if its possible to write using mpiio
    if (ownRank_ == 0) {
      std::ofstream file(filename_);
      file << content_;
    }
  } else {
    std::ofstream file(filename_);
    file << content_;
  }
}

nlohmann::json &File::getFileContent() { return content_; }
int32_t File::getOwnRank() const { return ownRank_; }
int32_t File::getWorldSize() const { return worldSize_; }
bool File::isMPIIO() const { return mpiio_; }

Group File::newGroup(const char *name) { return Group(this, name); }

Group::Group(File *file, const char *name)
    : file_(file), groupContent_(file->getFileContent()[name]) {}

Group::Group(File *file, nlohmann::json &obj, const char *name)
    : file_(file), groupContent_(obj[name]) {}

Group Group::newGroup(const char *name) {
  return Group(file_, groupContent_, name);
}
} // namespace JsonUtils
} // namespace OutputWriter
