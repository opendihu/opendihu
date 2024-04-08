#pragma once

#include <hdf5.h>

#include "output_writer/generic.h"
#include "output_writer/poly_data_properties_for_mesh.h"

namespace OutputWriter {
class HDF5 : public Generic {
public:
  HDF5(DihuContext context, PythonConfig specificSettings,
       std::shared_ptr<Partition::RankSubset> rankSubset = nullptr);

  //! write out solution to file, if timeStepNo is not -1, this value will be
  //! part of the filename
  template <typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1,
             int callCountIncrement = 1);

  template <typename FieldVariablesForOutputWriterType>
  void
  writePolyDataFile(hid_t fileID,
                    const FieldVariablesForOutputWriterType &fieldVariables,
                    std::set<std::string> &combinedMeshesOut);

  template <typename FieldVariablesForOutputWriterType>
  void writeCombinedUnstructuredGridFile(
      hid_t fileID, const FieldVariablesForOutputWriterType &fieldVariables,
      std::set<std::string> &combinedMeshesOut, bool output3DMeshes);

protected:
  /** one VTKPiece is the XML element that will be output as <Piece></Piece>. It
   * is created from one or multiple opendihu meshes
   */
  struct VTKPiece {
    std::set<std::string>
        meshNamesCombinedMeshes; //< the meshNames of the combined meshes, or
                                 // only one meshName if it is not a merged mesh
    std::vector<std::string>
        meshNamesCombinedMeshesVector; //< the same as meshNamesCombinedMeshes,
                                       // but as vector that preserves the
                                       // order, this is important for the
                                       // output file
    PolyDataPropertiesForMesh properties; //< the properties of the merged mesh

    std::string firstScalarName; //< name of the first scalar field variable of
                                 // the mesh (this is for Paraview such that it
                                 // selects this as the default scalar field)
    std::string
        firstVectorName; //< name of the first vector field variable with 3
                         // components of the mesh (this is for Paraview such
                         // that it selects this as the default vector field)

    //! constructor, initialize nPoints and nCells to 0
    VTKPiece();

    //! assign the correct values to firstScalarName and firstVectorName, only
    //! if properties has been set
    void setVTKValues();
  };

  //! write the values vector combined to the file, correctly encoded,
  //! identifier is an id to access cached values
  template <typename T>
  void writeCombinedValuesVector(hid_t fileID, const std::vector<T> &values,
                                 const char *dsname,
                                 bool writeFloatsAsInt = false);

  //! write a vector containing nValues "12" (if output3DMeshes) or "9" (if
  //! !output3DMeshes) values for the types for an unstructured grid
  void writeCombinedTypesVector(hid_t fileHandle, int nValues,
                                bool output3DMeshes, const char *dsname);

private:
  hid_t openHDF5File(const char *filename, bool mpiio);

  // TODO: make into template function
  herr_t writeAttrInt(hid_t fileID, const char *key, int value);
  herr_t writeAttrDouble(hid_t fileID, const char *key, double value);
  herr_t writeAttrString(hid_t fileID, const char *key,
                         const std::string &value);

  //! helper method that writes the unstructured grid file
  template <typename FieldVariablesForOutputWriterType>
  void writeCombinedUnstructuredGridFile(
      hid_t fileID, const FieldVariablesForOutputWriterType &fieldVariables,
      PolyDataPropertiesForMesh &polyDataPropertiesForMesh,
      const std::map<std::string, PolyDataPropertiesForMesh>
          &meshPropertiesUnstructuredGridFile,
      std::vector<std::string> meshNames, bool meshPropertiesInitialized);

  bool combineFiles_;

  std::map<std::string, PolyDataPropertiesForMesh>
      meshPropertiesPolyDataFile_; //< mesh information for a poly data file
                                   //(*.vtp), for 1D data
  std::map<std::string, PolyDataPropertiesForMesh>
      meshPropertiesUnstructuredGridFile2D_; //< mesh information for a combined
                                             // unstructured grid file (*.vtu),
                                             // for 2D data
  std::map<std::string, PolyDataPropertiesForMesh>
      meshPropertiesUnstructuredGridFile3D_; //< mesh information for a combined
                                             // unstructured grid file (*.vtu),
                                             // for 3D data
  VTKPiece
      vtkPiece1D_; //< the VTKPiece data structure used for PolyDataFile, 1D
  VTKPiece vtkPiece3D_; //< the VTKPiece data structure used for

  int nCellsPreviousRanks1D_ = 0; //< sum of number of cells on other processes
                                  // with lower rank no., for vtp file
  int nPointsPreviousRanks1D_ =
      0; //< sum of number of points on other
         // processes with lower rank no., for vtp file

  std::map<std::string, int>
      nCellsPreviousRanks3D_; //< sum of number of cells on other processes with
                              // lower rank no., for vtu file
  std::map<std::string, int>
      nPointsPreviousRanks3D_; //< sum of number of points on other processes
                               // with lower rank no., for vtu file
};

namespace HDF5Utils {
template <typename T>
static void writeSimpleVec(hid_t fileID, const std::vector<T> &data,
                           const char *dsname);
template <>
void writeSimpleVec<int32_t>(hid_t fileID, const std::vector<int32_t> &data,
                             const char *dsname);
template <>
void writeSimpleVec<double>(hid_t fileID, const std::vector<double> &data,
                            const char *dsname);

//! write the given field variable as VTK <DataArray> element to file, if
//! onlyParallelDatasetElement write the <PDataArray> element
template <typename FieldVariableType>
static void writeFieldVariable(hid_t fileID, FieldVariableType &fieldVariable);

//! write the a field variable indicating which ranks own which portion of the
//! domain as DataSet element to file
template <typename FieldVariableType>
static void writePartitionFieldVariable(hid_t fileID,
                                        FieldVariableType &geometryField);
} // namespace HDF5Utils
} // namespace OutputWriter

#include "output_writer/hdf5/hdf5.tpp"
#include "output_writer/hdf5/hdf5_write_combined_values.tpp"
#include "output_writer/hdf5/hdf5_write_combined_file_1D.tpp"
#include "output_writer/hdf5/hdf5_write_combined_file_2D3D.tpp"
