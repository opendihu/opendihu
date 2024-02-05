#include "output_writer/python_callback/python_callback.h"

#include <Python.h> // has to be the first included header
#include <iostream>

#include "easylogging++.h"
#include "output_writer/python_callback/python_callback_writer.h"

namespace OutputWriter {

template <typename DataType>
void PythonCallback::write(DataType &data, int timeStepNo, double currentTime,
                           int callCountIncrement) {
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime,
                             callCountIncrement)) {
    return;
  }

  LOG(TRACE) << "PythonCallback::write timeStepNo=" << timeStepNo
             << ", currentTime=" << currentTime;

  // call implementation specific for FunctionSpace type
  PythonCallbackWriter<typename DataType::FunctionSpace,
                       typename DataType::FieldVariablesForOutputWriter>::
      callCallback(callback_, data.getFieldVariablesForOutputWriter(),
                   this->timeStepNo_, this->currentTime_,
                   this->onlyNodalValues_);
}

} // namespace OutputWriter
