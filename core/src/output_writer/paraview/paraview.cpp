#include "output_writer/paraview/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include <utility/python_utility.h>
#include <utility/petsc_utility.h>
#include <mesh/structured_regular_fixed.h>
#include <mesh/structured_deformable.h>
#include <mesh/unstructured_deformable.h>
#include <mesh/mesh.h>
#include "output_writer/paraview/series_writer.h"

namespace OutputWriter {

SeriesWriter Paraview::seriesWriter_;

Paraview::Paraview(DihuContext context, PythonConfig settings,
                   std::shared_ptr<Partition::RankSubset> rankSubset)
    : Generic(context, settings, rankSubset) {
  binaryOutput_ = settings.getOptionBool("binary", true);
  fixedFormat_ = settings.getOptionBool("fixedFormat", true);
  combineFiles_ = settings.getOptionBool("combineFiles", false);
}

std::string Paraview::encodeBase64Vec(const Vec &vector,
                                      bool withEncodedSizePrefix) {
  PetscInt vectorSize = 0;
  VecGetSize(vector, &vectorSize);

  std::vector<dof_no_t> indices(vectorSize);
  std::iota(indices.begin(), indices.end(),
            0); // fill with increasing numbers: 0,1,2,...
  std::vector<double> values(vectorSize);
  VecGetValues(vector, vectorSize, indices.data(), values.data());

  return encodeBase64Float(values.begin(), values.end(), withEncodedSizePrefix);
}

std::string Paraview::convertToAscii(const Vec &vector, bool fixedFormat) {
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(vector, vectorValues);

  return convertToAscii(vectorValues, fixedFormat);
}

std::string Paraview::convertToAscii(const std::vector<double> &vector,
                                     bool fixedFormat) {
  std::stringstream result;
  for (int i = 0; i < vector.size(); i++) {
    if (fixedFormat) {
      result << std::setw(16) << std::scientific << vector[i] << " ";
    } else {
      result << vector[i] << " ";
    }
    /*if (i % 3 == 2)   // this breaks validity for paraview but creates better
    analysable results (1 point per row)
    {
      result << std::endl << std::string(5,'\t');
    }*/
  }
  return result.str();
}

std::string Paraview::convertToAscii(const std::vector<element_no_t> &vector,
                                     bool fixedFormat) {
  std::stringstream result;
  for (auto value : vector) {
    if (fixedFormat) {
      result << std::setw(16) << std::scientific << (float)(value) << " ";
    } else {
      result << value << " ";
    }
  }
  return result.str();
}

#ifdef PETSC_USE_64BIT_INDICES
std::string Paraview::convertToAscii(const std::vector<int> &vector,
                                     bool fixedFormat) {
  std::stringstream result;
  for (auto value : vector) {
    if (fixedFormat) {
      result << std::setw(16) << std::scientific << (float)(value) << " ";
    } else {
      result << value << " ";
    }
  }
  return result.str();
}
#endif

SeriesWriter &Paraview::seriesWriter() { return Paraview::seriesWriter_; }

} // namespace OutputWriter
