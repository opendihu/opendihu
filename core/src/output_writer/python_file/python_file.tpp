#include "output_writer/python_file/python_file.h"

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <cstdio>

#include "easylogging++.h"
#include "output_writer/python/python.h"
#include "output_writer/python_file/python_stiffness_matrix_writer.h"

namespace OutputWriter
{

template<typename DataType>
void PythonFile::write(DataType& data, int timeStepNo, double currentTime)
{ 
  LOG(TRACE) << "PythonFile::write ";
 
  // check if output should be written in this timestep and prepare filename
  if (!Generic::prepareWrite(data, timeStepNo, currentTime))
  {
    return;
  }
  
  // build python object for data
  PyObject *pyData = Python<typename DataType::BasisOnMesh, typename DataType::OutputFieldVariables>::
    buildPyDataObject(data.getOutputFieldVariables(), timeStepNo, currentTime, this->onlyNodalValues_);
  
  // determine file name
  std::stringstream s;
  s<<this->filename_<<".py";
  std::string filename = s.str();
  
  LOG(DEBUG) << "filename is [" << filename << "]";
  
  // open file, to see if directory needs to be created
  std::ofstream ofile = openFile(filename);
  if(ofile.is_open())
    ofile.close();
  
#if PY_MAJOR_VERSION >= 3  
  FILE *fileC = fopen(filename.c_str(), "w");
  int fileDescriptorC = fileno(fileC);
  PyObject *file = PyFile_FromFd(fileDescriptorC, NULL, "w", -1, NULL, NULL, NULL, true);
#else 
  // python 2.7
  char filenameC[filename.size()+1];
  std::strcpy(filenameC, filename.c_str());
  
  char mode[3];
  std::strcpy(mode, "wb");

  PyObject *file = PyFile_FromString(filenameC, mode);
#endif
 
  if (!file) 
  {
    LOG(ERROR) << "Could not open file \""<<filename<<"\" for output of python object.";
    return;
  }
  
  // pickle is the python library to serialize objects
  bool usePickle = PythonUtility::getOptionBool(specificSettings_, "binary", false);
  
  if (usePickle)
  {
    // load pickle module if is was not loaded in an earlier call
    static PyObject *module = NULL;
    if (module == NULL)
    {
      //module = PyImport_ImportModuleNoBlock("cpickle");
      if (module == NULL)
        module = PyImport_ImportModuleNoBlock("pickle");
    }
    
    if (module == NULL)
    {
      LOG(ERROR) << "Could not load python pickle package";
    }
    else
    {
      PyObject *pickle;
      std::string methodNameStr("dump");
      char methodName[methodNameStr.size()+1];
      std::strcpy(methodName, methodNameStr.c_str());
      
      std::string formatStr("(O O i)");
      char format[formatStr.size()+1];
      std::strcpy(format, formatStr.c_str());
      
      pickle = PyObject_CallMethod(module, methodName, format, pyData, file, 1);
      Py_XDECREF(pickle);
    }
  }
  else
  {
    // output python object
    PyFile_WriteObject(pyData, file, 0);
    LOG(DEBUG) << "PyFile_WriteObject done";
  }
  
  // decrement reference counters for python objects
  Py_XDECREF(pyData);
  Py_XDECREF(file);
  
  LOG(DEBUG) << "writeNumpySolution";
  
  // for regular fixed also output stiffness matrix 
  PythonStiffnessMatrixWriter<DataType>::writeNumpySolution(data, this->filename_);
}

};