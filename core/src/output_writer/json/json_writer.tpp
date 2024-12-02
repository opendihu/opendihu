#include "output_writer/json/json_writer.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "output_writer/json/json.h"
#include "output_writer/json/loop_output_point_data.h"
#include "field_variable/field_variable.h"

namespace OutputWriter {
// regular fixed
template <int D, typename BasisFunctionType,
          typename FieldVariablesForOutputWriterType>
void JsonWriter<
    FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,
                                 BasisFunctionType>,
    FieldVariablesForOutputWriterType>::
    outputFile(
        JsonUtils::Group &group,
        FieldVariablesForOutputWriterType fieldVariables,
        const std::string &meshName,
        std::shared_ptr<FunctionSpace::FunctionSpace<
            Mesh::StructuredRegularFixedOfDimension<D>, BasisFunctionType>>
            mesh,
        int nFieldVariablesOfMesh, const PythonConfig &specificSettings,
        double currentTime) {
  // write a RectilinearGrid

  // get type of geometry field
  typedef FieldVariable::FieldVariable<
      FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<D>,
                                   BasisFunctionType>,
      3>
      GeometryFieldType;

  LOG(DEBUG) << "Write RectilinearGrid";
  // coordinates of grid
  std::array<std::vector<double>, 3> coordinates;
  int dimensionNo = 0;
  for (; dimensionNo < D; dimensionNo++) {
    double meshWidth = mesh->meshWidth();
    node_no_t nNodes =
        mesh->meshPartition()->nNodesLocalWithGhosts(dimensionNo);

    LOG(DEBUG) << "dimension " << dimensionNo << ", meshWidth: " << meshWidth;

    coordinates[dimensionNo].resize(nNodes);

    for (node_no_t nodeNo = 0; nodeNo < nNodes; nodeNo++) {
      double coordinate =
          (mesh->meshPartition()->beginNodeGlobalNatural(dimensionNo) +
           nodeNo) *
          meshWidth;
      VLOG(1) << "coordinate: " << coordinate << ", nodeNo=" << nodeNo;
      coordinates[dimensionNo][nodeNo] = coordinate;
    }
  }

  // set other coordinates to 0
  for (; dimensionNo < 3; dimensionNo++) {
    coordinates[dimensionNo].resize(1);
    coordinates[dimensionNo][0] = 0.0;
  }

  JsonLoopOverTuple::loopOutputPointData(group, fieldVariables, meshName,
                                         false);
  JsonUtils::writePartitionFieldVariable<GeometryFieldType>(
      group, mesh->geometryField());
  group.writeSimpleVec(coordinates[0], "coordinates_0");
  group.writeSimpleVec(coordinates[1], "coordinates_1");
  group.writeSimpleVec(coordinates[2], "coordinates_2");
}

// structured deformable
template <int D, typename BasisFunctionType,
          typename FieldVariablesForOutputWriterType>
void JsonWriter<
    FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    FieldVariablesForOutputWriterType>::
    outputFile(JsonUtils::Group &group,
               FieldVariablesForOutputWriterType fieldVariables,
               const std::string &meshName,
               std::shared_ptr<FunctionSpace::FunctionSpace<
                   Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType>>
                   mesh,
               int nFieldVariablesOfMesh, const PythonConfig &specificSettings,
               double currentTime) {
  // write a StructuredGrid

  // get type of geometry field
  typedef FieldVariable::FieldVariable<
      FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,
                                   BasisFunctionType>,
      3>
      GeometryFieldType;

  LOG(DEBUG) << "Write StructuredGrid";

  JsonLoopOverTuple::loopOutputPointData(group, fieldVariables, meshName,
                                         false);
  JsonUtils::writePartitionFieldVariable<GeometryFieldType>(
      group, mesh->geometryField());
  JsonUtils::writeFieldVariable<GeometryFieldType>(group,
                                                   mesh->geometryField());
}

// unstructured deformable
template <int D, typename BasisFunctionType,
          typename FieldVariablesForOutputWriterType>
void JsonWriter<
    FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,
                                 BasisFunctionType>,
    FieldVariablesForOutputWriterType>::
    outputFile(
        JsonUtils::Group &group,
        FieldVariablesForOutputWriterType fieldVariables,
        const std::string &meshName,
        std::shared_ptr<FunctionSpace::FunctionSpace<
            Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>
            mesh,
        int nFieldVariablesOfMesh, const PythonConfig &specificSettings,
        double currentTime) {
  LOG(DEBUG) << "Write UnstructuredGrid";

  // get type of geometry field
  typedef FunctionSpace::FunctionSpace<
      Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>
      FunctionSpace;
  typedef FieldVariable::FieldVariable<FunctionSpace, 3> GeometryFieldType;

  JsonLoopOverTuple::loopOutputPointData(group, fieldVariables, meshName,
                                         false);
  JsonUtils::writePartitionFieldVariable<GeometryFieldType>(
      group, mesh->geometryField());
  JsonUtils::writeFieldVariable<GeometryFieldType>(group,
                                                   mesh->geometryField());

  // get the elements point lists
  std::vector<node_no_t> values;
  values.reserve(mesh->nElementsLocal() *
                 FunctionSpace::averageNNodesPerElement());

  // loop over elements and collect point numbers of the element
  for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal();
       elementNo++) {
    std::array<dof_no_t, FunctionSpace::nDofsPerElement()> dofsOfElement =
        mesh->getElementDofNosLocal(elementNo);
    for (typename std::array<
             dof_no_t, FunctionSpace::nDofsPerElement()>::const_iterator iter =
             dofsOfElement.begin();
         iter != dofsOfElement.end(); iter++) {
      dof_no_t dofNo = *iter;
      if (dofNo % FunctionSpace::nDofsPerNode() == 0) {
        node_no_t nodeNo = dofNo / FunctionSpace::nDofsPerNode();
        values.push_back(nodeNo);
      }
    }
  }
  group.writeSimpleVec(values, "connectivity");

  // offsets
  values.clear();
  values.resize(mesh->nElementsLocal());
  for (element_no_t elementNo = 0; elementNo < mesh->nElementsLocal();
       elementNo++) {
    values[elementNo] = (elementNo + 1) * FunctionSpace::nNodesPerElement();
  }
  group.writeSimpleVec(values, "offsets");

  // cell types
  int cellType = 0;
  switch (D) {
  case 1:
    cellType = 3; // VTK_LINE
    break;
  case 2:
    cellType = 8; // VTK_PIXEL
    break;
  case 3:
    cellType = 11; // VTK_VOXEL
    break;
  }
  values.clear();
  values.resize(mesh->nElementsLocal(), cellType);
  group.writeSimpleVec(values, "types");
}
} // namespace OutputWriter
