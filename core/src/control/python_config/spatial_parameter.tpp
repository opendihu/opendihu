#include "control/python_config/spatial_parameter.h"

//! initialize the spatial parameter from the settings, at the given parameter
//! name "keyString"
template <typename FunctionSpaceType, typename ValueType>
void SpatialParameterBase<FunctionSpaceType, ValueType>::initialize(
    PythonConfig specificSettings, std::string keyString,
    ValueType defaultValue, std::shared_ptr<FunctionSpaceType> functionSpace) {
  inputMeshIsGlobal_ =
      specificSettings.getOptionBool("inputMeshIsGlobal", true);

  std::stringstream path;
  path << specificSettings << "[\"" << keyString << "\"]";

  VLOG(1) << "initialize option " << path.str();

  // get option
  PyObject *pyObject = specificSettings.getOptionPyObject(keyString);

  // default value
  if (pyObject == Py_None) {
    LOG(WARNING) << path.str() << " not set in \"" << Control::settingsFileName
                 << "\"."
                 << " Assuming default value \"" << defaultValue << "\".";
    values_.resize(1, defaultValue);
    valueIndices_.resize(functionSpace->nElementsLocal(), 0);
  } else {
    if (PyList_Check(pyObject) || PyDict_Check(pyObject)) {

      // if the type is Matrix or array and only one list is given, this is
      // interpreted as constant matrix for all elements (instead of different
      // scalar values for the elements)
      if (PyList_Check(pyObject) &&
          MathUtility::isArrayOrMatrix<ValueType>::value) {

        if (PyList_Size(pyObject) > 0) {
          // if the first item is scalar
          PyObject *firstItem = PyList_GetItem(pyObject, 0);
          if (!PyList_Check(firstItem)) {
            LOG(DEBUG) << path.str()
                       << " interpreting list as a single matrix entry of type "
                       << StringUtility::demangle(typeid(ValueType).name());

            values_.push_back(
                PythonUtility::convertFromPython<ValueType>::get(pyObject));
          }
        }
      }

      // if the above was not the case, interpret given values as list of
      // different values for the elements
      if (values_.empty()) {
        values_ = PythonUtility::convertFromPython<std::vector<ValueType>>::get(
            pyObject);
      }

      // set indices in valueIndices_ according to the number of values that
      // were present in the config
      SetValueNos<FunctionSpaceType, ValueType>::set(
          functionSpace, inputMeshIsGlobal_, path.str(), values_,
          valueIndices_);
    } else {
      // if it is a single value, use it for all elements
      values_.resize(
          1, PythonUtility::convertFromPython<ValueType>::get(pyObject));
      valueIndices_.resize(functionSpace->nElementsLocal(), 0);
    }
  }
}

//! get the value of the parameter in the specified element
template <typename FunctionSpaceType, typename ValueType>
void SpatialParameterBase<FunctionSpaceType, ValueType>::getValue(
    element_no_t elementNoLocal, ValueType &value) const {
  assert(elementNoLocal >= 0 && elementNoLocal < valueIndices_.size());
  int index = valueIndices_[elementNoLocal];
  value = values_[index];
}

//! the value of the parameter in the given element
/*template <typename FunctionSpaceType, typename ValueType>
ValueType SpatialParameterBase<FunctionSpaceType,ValueType>::
value(element_no_t elementNoLocal) const
{
  assert (elementNoLocal >= 0 && elementNoLocal < valueIndices_.size());
  int index = valueIndices_[elementNoLocal];
  return values_[index];
}*/

//! the value of the parameter in the given element
template <typename FunctionSpaceType, typename ValueType>
const ValueType &SpatialParameterBase<FunctionSpaceType, ValueType>::value(
    element_no_t elementNoLocal) const {
  assert(elementNoLocal >= 0 && elementNoLocal < valueIndices_.size());
  int index = valueIndices_[elementNoLocal];
  return values_[index];
}

//! the value of the parameter in the given element
template <typename FunctionSpaceType>
Vc::double_v SpatialParameter<FunctionSpaceType, double>::value(
    Vc::int_v elementNoLocal) const {
  assert(elementNoLocal[0] >= 0 &&
         elementNoLocal[0] < this->valueIndices_.size());
  Vc::double_v result;

  for (int vcComponentNo = 0; vcComponentNo < Vc::double_v::size();
       vcComponentNo++) {
    int elementNo = elementNoLocal[vcComponentNo];
    if (elementNo != -1) {
      int index = this->valueIndices_[elementNo];
      result[vcComponentNo] = (double)this->values_[index];
    }
  }

  return result;
}

//! the value of the parameter in the given element
template <typename FunctionSpaceType, int nRows, int nColumns>
MathUtility::Matrix<nRows, nColumns, Vc::double_v>
SpatialParameter<FunctionSpaceType, MathUtility::Matrix<nRows, nColumns>>::
    value(Vc::int_v elementNoLocal) const {
  assert(elementNoLocal[0] >= 0 &&
         elementNoLocal[0] < this->valueIndices_.size());
  MathUtility::Matrix<nRows, nColumns, Vc::double_v> result;

  for (int vcComponentNo = 0; vcComponentNo < Vc::double_v::size();
       vcComponentNo++) {
    int elementNo = elementNoLocal[vcComponentNo];
    if (elementNo != -1) {
      int index = this->valueIndices_[elementNo];
      for (int rowNo = 0; rowNo < nRows; rowNo++) {
        for (int columnNo = 0; columnNo < nColumns; columnNo++) {
          result(rowNo, columnNo)[vcComponentNo] =
              this->values_[index](rowNo, columnNo);
        }
      }
    }
  }

  return result;
}

template <typename FunctionSpaceType, typename ValueType>
void SetValueNos<FunctionSpaceType, ValueType>::set(
    std::shared_ptr<FunctionSpaceType> functionSpace, bool inputMeshIsGlobal,
    std::string path, std::vector<ValueType> &values,
    std::vector<int> &valueIndices) {
  element_no_t nElementsLocal = functionSpace->nElementsLocal();
  global_no_t nElementsGlobal = functionSpace->nElementsGlobal();

  VLOG(1) << path
          << ", set value indices for normal mesh, n values: " << values.size()
          << ", nElementsLocal: " << nElementsLocal
          << ", nElementsGlobal: " << nElementsGlobal
          << ", inputMeshIsGlobal: " << inputMeshIsGlobal;

  // fill with natural integers starting from 0
  valueIndices.resize(values.size());
  std::iota(valueIndices.begin(), valueIndices.end(), 0);

  // set last indices for which no entries were given to point to the 0th value
  valueIndices.resize(nElementsLocal, 0);

  // if number of given values in the config equals the number of elements in
  // the function space
  if (inputMeshIsGlobal) {
    if (values.size() != nElementsGlobal && values.size() != 1) {
      LOG(ERROR) << path << ": The number of entries given in the list ("
                 << values.size() << ") does not match "
                 << "the global number of elements (" << nElementsGlobal
                 << "). \"inputMeshIsGlobal\" is True. \n"
                 << "Note, if you only specify one entry, it will be used for "
                    "all elements. "
                 << " If you specify more than one, every element gets its own "
                    "value and, thus, the number of entries must match the "
                    "number of elements.";
    } else if (values.size() != 1) {
      std::vector<ValueType> localValues;
      localValues.reserve(nElementsLocal);

      // extract local element values
      for (element_no_t elementNoLocal = 0; elementNoLocal < nElementsLocal;
           elementNoLocal++) {
        global_no_t elementNoGlobalNatural =
            functionSpace->meshPartition()->getElementNoGlobalNatural(
                elementNoLocal);

        assert(elementNoGlobalNatural >= 0 &&
               elementNoGlobalNatural < values.size());
        localValues.push_back(values[elementNoGlobalNatural]);
      }
      values.assign(localValues.begin(), localValues.end());
    }
  } else {
    if (values.size() != nElementsLocal && values.size() != 1) {
      LOG(ERROR) << path << ": The number of entries given in the list ("
                 << values.size() << ") does not match "
                 << "the local number of elements (" << nElementsLocal
                 << "). \"inputMeshIsGlobal\" is False. \n"
                 << "Note, if you only specify one entry, it will be used for "
                    "all elements. "
                 << " If you specify more than one, every element gets its own "
                    "value and, thus, the number of entries must match the "
                    "number of elements.";
    }
  }

  VLOG(1) << "values: " << values << ", nElementsLocal: " << nElementsLocal
          << ", nElementsGlobal: " << nElementsGlobal
          << ", valueIndices: " << valueIndices;
}

template <int D, typename BasisFunctionType, typename ValueType>
void SetValueNos<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<D>,
                                              BasisFunctionType>,
                 ValueType>::
    set(std::shared_ptr<FunctionSpace::FunctionSpace<
            Mesh::CompositeOfDimension<D>, BasisFunctionType>>
            functionSpace,
        bool inputMeshIsGlobal, std::string path,
        std::vector<ValueType> &values, std::vector<int> &valueIndices) {
  element_no_t nElementsLocal = functionSpace->nElementsLocal();
  global_no_t nElementsGlobal = functionSpace->nElementsGlobal();

  VLOG(1) << path << ", set value indices for composite mesh, n values: "
          << values.size() << ", nElementsLocal: " << nElementsLocal
          << ", nElementsGlobal: " << nElementsGlobal
          << ", inputMeshIsGlobal: " << inputMeshIsGlobal
          << ", nSubMeshes: " << functionSpace->meshPartition()->nSubMeshes();

  if (values.size() == functionSpace->meshPartition()->nSubMeshes()) {
    LOG(DEBUG) << path << " contains a value for each submesh: " << values;

    // loop over local elements
    for (element_no_t elementNoLocal = 0; elementNoLocal < nElementsLocal;
         elementNoLocal++) {
      int subMeshNo = 0;
      element_no_t elementOnMeshNoLocal;
      functionSpace->meshPartition()->getSubMeshNoAndElementNoLocal(
          elementNoLocal, subMeshNo, elementOnMeshNoLocal);

      valueIndices.push_back(subMeshNo);
    }
  } else {
    valueIndices.resize(values.size());
    std::iota(valueIndices.begin(), valueIndices.end(), 0);
    valueIndices.resize(nElementsLocal, 0);

    // if number of given values in the config equals the number of elements in
    // the function space
    if (inputMeshIsGlobal) {
      if (values.size() != nElementsGlobal && values.size() != 1) {
        LOG(ERROR) << path << ": The number of entries given (" << values.size()
                   << ") does not match "
                   << "the global number of elements (" << nElementsGlobal
                   << "). \"inputMeshIsGlobal\" is True.";
      } else {
        std::vector<ValueType> localValues;
        localValues.reserve(nElementsLocal);

        // extract local element values
        for (element_no_t elementNoLocal = 0; elementNoLocal < nElementsLocal;
             elementNoLocal++) {
          global_no_t elementNoGlobalNatural =
              functionSpace->meshPartition()->getElementNoGlobalNatural(
                  elementNoLocal);

          assert(elementNoGlobalNatural > 0 &&
                 elementNoGlobalNatural < values.size());
          localValues.push_back(values[elementNoGlobalNatural]);
        }
        values.assign(localValues.begin(), localValues.end());
      }
    } else {
      if (values.size() != nElementsLocal && values.size() != 1) {
        LOG(ERROR) << path << ": The number of entries given (" << values.size()
                   << ") does not match "
                   << "the local number of elements (" << nElementsLocal
                   << "). \"inputMeshIsGlobal\" is False.";
      }
    }
  }

  VLOG(1) << "values: " << values << ", nElementsLocal: " << nElementsLocal
          << ", nElementsGlobal: " << nElementsGlobal
          << ", valueIndices: " << valueIndices;
}
