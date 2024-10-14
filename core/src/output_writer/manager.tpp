#include "output_writer/manager.h"

#include "output_writer/python_callback/python_callback.h"
#include "output_writer/python_file/python_file.h"
#include "output_writer/paraview/paraview.h"
#include "output_writer/exfile/exfile.h"
#include "output_writer/megamol/megamol.h"
#include "output_writer/hdf5/hdf5.h"
#include "control/diagnostic_tool/performance_measurement.h"

namespace OutputWriter {

template <typename DataType>
void Manager::writeOutput(DataType &problemData, int timeStepNo,
                          double currentTime, int callCountIncrement) const {
  // start duration measurement
  Control::PerformanceMeasurement::start("durationWriteOutput");

  for (auto &outputWriter : this->outputWriter_) {
    if (std::dynamic_pointer_cast<Exfile>(outputWriter) != nullptr) {
      LogScope s("WriteOutputExfile");
      Control::PerformanceMeasurement::start("durationWriteOutputExfile");

      std::shared_ptr<Exfile> writer =
          std::static_pointer_cast<Exfile>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime,
                              callCountIncrement);

      Control::PerformanceMeasurement::stop("durationWriteOutputExfile");
    } else if (std::dynamic_pointer_cast<Paraview>(outputWriter) != nullptr) {
      LogScope s("WriteOutputParaview");
      Control::PerformanceMeasurement::start("durationWriteOutputParaview");

      std::shared_ptr<Paraview> writer =
          std::static_pointer_cast<Paraview>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime,
                              callCountIncrement);

      Control::PerformanceMeasurement::stop("durationWriteOutputParaview");
    } else if (std::dynamic_pointer_cast<PythonCallback>(outputWriter) !=
               nullptr) {
      LogScope s("WriteOutputPythonCallback");
      Control::PerformanceMeasurement::start(
          "durationWriteOutputPythonCallback");

      std::shared_ptr<PythonCallback> writer =
          std::static_pointer_cast<PythonCallback>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime,
                              callCountIncrement);

      Control::PerformanceMeasurement::stop(
          "durationWriteOutputPythonCallback");
    } else if (std::dynamic_pointer_cast<PythonFile>(outputWriter) != nullptr) {
      LogScope s("WriteOutputPythonFile");
      Control::PerformanceMeasurement::start("durationWriteOutputPythonFile");

      std::shared_ptr<PythonFile> writer =
          std::static_pointer_cast<PythonFile>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime,
                              callCountIncrement);

      Control::PerformanceMeasurement::stop("durationWriteOutputPythonFile");
    } else if (std::dynamic_pointer_cast<MegaMol>(outputWriter) != nullptr) {
      LogScope s("WriteOutputMegamol");
      Control::PerformanceMeasurement::start("durationWriteOutputMegamol");

      std::shared_ptr<MegaMol> writer =
          std::static_pointer_cast<MegaMol>(outputWriter);
      writer->write<DataType>(problemData, timeStepNo, currentTime,
                              callCountIncrement);

      Control::PerformanceMeasurement::stop("durationWriteOutputMegamol");
    } else if (std::dynamic_pointer_cast<HDF5>(outputWriter) != nullptr) {
      LogScope s("WriteOutputHDF5");
      Control::PerformanceMeasurement::start("durationWriteOutputHDF5");

      std::shared_ptr<HDF5> writer =
          std::static_pointer_cast<HDF5>(outputWriter);
      writer->write<DataType>(problemData, nullptr, timeStepNo, currentTime,
                              callCountIncrement);

      Control::PerformanceMeasurement::stop("durationWriteOutputHDF5");
    }
  }

  // stop duration measurement
  Control::PerformanceMeasurement::stop("durationWriteOutput");
}

} // namespace OutputWriter
