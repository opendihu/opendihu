#pragma once

#include <Python.h> // has to be the first included header

#include "precice/precice.hpp"
// #include "control/precice/initialize.h"

namespace Control {

class ReadWriteDataBase {
public:
  void preciceReadVolumeData(
      std::shared_ptr<precice::Participant> &preciceParticipant);
  void preciceWriteVolumeData(
      std::shared_ptr<precice::Participant> &preciceParticipant);

  void preciceReadSurfaceData(
      std::shared_ptr<precice::Participant> &preciceParticipant);
  void preciceWriteSurfaceData(
      std::shared_ptr<precice::Participant> &preciceParticipant);

  // void setNeumannBoundaryConditions(typename
  // PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData &preciceData);
  // void setDirichletBoundaryConditions(typename
  // PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData &preciceData);

  std::vector<double> displacementValues_;
  std::vector<double> velocityValues_;
  std::vector<double> tractionValues_;
  std::vector<Vec3> displacementVectors_;
  std::vector<Vec3> velocityVectors_;
  std::vector<Vec3> tractionVectors_;
  std::vector<double> scalarValues_;
  std::vector<double> scalarValuesOfMesh_;
  std::vector<Vec3> geometryValues_;
};

} // namespace Control

#include "control/precice/read_write.tpp"
