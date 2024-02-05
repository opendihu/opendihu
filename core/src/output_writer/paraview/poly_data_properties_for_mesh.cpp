#include "output_writer/paraview/poly_data_properties_for_mesh.h"

#include <Python.h> // has to be the first included header
#include <iostream>
#include <vector>

#include "utility/vector_operators.h"
#include "control/types.h"

namespace OutputWriter {

// output operator
std::ostream &operator<<(std::ostream &stream, PolyDataPropertiesForMesh rhs) {
  stream << "(D: " << rhs.dimensionality
         << ", nPointsLocal: " << rhs.nPointsLocal
         << ", nCellsLocal: " << rhs.nCellsLocal
         << ", nPointsGlobal: " << rhs.nPointsGlobal
         << ", nCellsGlobal: " << rhs.nCellsGlobal << ", pointDataArrays:";
  for (auto pointDataArray : rhs.pointDataArrays) {
    stream << " [" << pointDataArray.name << "," << pointDataArray.nComponents
           << ": ";
    for (int componentNo = 0;
         componentNo < pointDataArray.componentNames.size(); componentNo++)
      stream << "\"" << pointDataArray.componentNames[componentNo] << "\" ";
    stream << "]";
  }
  stream << ", " << rhs.unstructuredMeshConnectivityValues.size()
         << ", unstructuredMeshConnectivityValues: ";

  for (int i = 0; i < rhs.unstructuredMeshConnectivityValues.size(); i++) {
    stream << rhs.unstructuredMeshConnectivityValues[i] << " ";
  }

  stream << ")";
  return stream;
}

} // namespace OutputWriter
