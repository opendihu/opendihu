#pragma once

#include <Python.h> // has to be the first included header
#include <list>
#include <memory>

#include "control/types.h"
#include "data_management/data.h"
#include "output_writer/generic.h"

namespace OutputWriter {

class Manager {
public:
  /** call all output writers to write output, timeStepNo of -1 means no time
   * step number in output filename
   * @param timeStepNo -1 means no time step number in output filename
   * @param callCountIncrement if the output should be written now, not
   * considering the set output frequency, this is only needed for
   * FastMonodomainSolver
   */
  void writeOutput(int timeStepNo, double currentTime);

  template <typename DataType>
  void writeOutput(DataType &problemData, int timeStepNo = -1,
                   double currentTime = 0.0, int callCountIncrement = 1) const;

  //! parse settings and create output writers from specification in
  //! "OutputWriter" list, if rankSubset is not given, use the
  //! rankSubsetForCollectiveOperations, which is the collection of all
  //! available ranks
  void initialize(DihuContext context, PythonConfig settings,
                  std::shared_ptr<Partition::RankSubset> rankSubset = nullptr);

  //! if this manager contains any output writers
  bool hasOutputWriters();

  //! get the filename of the first output writer
  std::string filename();

  //! set the filename for the first output writer
  void setFilename(std::string filename);

protected:
  //! helper function that creates an outputWriter
  void createOutputWriterFromSettings(
      DihuContext context, PythonConfig settings,
      std::shared_ptr<Partition::RankSubset> rankSubset);

  std::list<std::shared_ptr<Generic>>
      outputWriter_; //< list of active output writers
};

} // namespace OutputWriter

#include "output_writer/manager.tpp"
