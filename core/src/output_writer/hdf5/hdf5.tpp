#include "output_writer/hdf5/hdf5.h"

#include "output_writer/loop_count_n_field_variables_of_mesh.h"

namespace OutputWriter {
template <typename DataType>
void HDF5::write(DataType &data, int timeStepNo, double currentTime,
                 int callCountIncrement) {
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime,
                             callCountIncrement)) {
    return;
  }
  Control::PerformanceMeasurement::start("durationHDF5Output");

  std::set<std::string> combined1DMeshes;
  std::set<std::string> combined2DMeshes;
  std::set<std::string> combined3DMeshes;

  // determine filename, broadcast from rank 0
  std::stringstream filename;
  filename << this->filenameBaseWithNo_ << ".h5";
  int filenameLength = filename.str().length();

  // broadcast length of filename
  MPIUtility::handleReturnValue(MPI_Bcast(&filenameLength, 1, MPI_INT, 0,
                                          this->rankSubset_->mpiCommunicator()),
                                "MPI_Bcast");

  std::vector<char> receiveBuffer(filenameLength + 1, char(0));
  strcpy(receiveBuffer.data(), filename.str().c_str());
  MPIUtility::handleReturnValue(MPI_Bcast(receiveBuffer.data(), filenameLength,
                                          MPI_CHAR, 0,
                                          this->rankSubset_->mpiCommunicator()),
                                "MPI_Bcast");
  LOG(DEBUG) << "open HDF5 file using MPI IO \"" << receiveBuffer.data()
             << "\".";
  hid_t fileID = openHDF5File(receiveBuffer.data());

  herr_t err;
  err = this->writeAttrString(fileID, "version", DihuContext::versionText());
  assert(err >= 0);

  err = this->writeAttrString(fileID, "meta", DihuContext::metaText());
  assert(err >= 0);

  err = this->writeAttrDouble(fileID, "currentTime", this->currentTime_);
  assert(err >= 0);

  err = this->writeAttrInt(fileID, "timeStepNo", this->timeStepNo_);
  assert(err >= 0);

  if (combineFiles_) {
    Control::PerformanceMeasurement::start("durationHDF51D");

    LOG(DEBUG)
        << "FieldVariablesForOutputWriter: "
        << StringUtility::demangle(
               typeid(typename DataType::FieldVariablesForOutputWriter).name());

    // create a PolyData file that combines all 1D meshes into one file
    writePolyDataFile<typename DataType::FieldVariablesForOutputWriter>(
        fileID, data.getFieldVariablesForOutputWriter(), combined1DMeshes);

    Control::PerformanceMeasurement::stop("durationHDF51D");
    Control::PerformanceMeasurement::start("durationHDF53D");

    // create an UnstructuredMesh file that combines all 3D meshes into one file
    writeCombinedUnstructuredGridFile<
        typename DataType::FieldVariablesForOutputWriter>(
        fileID, data.getFieldVariablesForOutputWriter(), combined3DMeshes,
        true);

    Control::PerformanceMeasurement::stop("durationHDF53D");
    Control::PerformanceMeasurement::start("durationHDF52D");

    // create an UnstructuredMesh file that combines all 2D meshes into one file
    writeCombinedUnstructuredGridFile<
        typename DataType::FieldVariablesForOutputWriter>(
        fileID, data.getFieldVariablesForOutputWriter(), combined2DMeshes,
        false);

    Control::PerformanceMeasurement::stop("durationHDF52D");
  }

  err = H5Fclose(fileID);
  assert(err >= 0);
  Control::PerformanceMeasurement::stop("durationHDF5Output");
}
} // namespace OutputWriter
