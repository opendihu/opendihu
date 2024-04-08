#pragma once

#include <hdf5.h>

#include "output_writer/generic.h"
#include "output_writer/poly_data_properties_for_mesh.h"

namespace OutputWriter {
class HDF5 : public Generic {
public:
  //! constructor
  HDF5(DihuContext context, PythonConfig specificSettings,
       std::shared_ptr<Partition::RankSubset> rankSubset = nullptr);

  //! write out solution to file, if timeStepNo is not -1, this value will be
  //! part of the filename
  template <typename DataType>
  void write(DataType &data, int timeStepNo = -1, double currentTime = -1,
             int callCountIncrement = 1);

  //! write all 1D field variables into a group. This is uses MPI IO.
  //! It can be enabled with the "combineFiles" option. on return,
  //! combinedMeshesOut contains the 1D mesh names that were written to the hdf5
  //! file.
  template <typename FieldVariablesForOutputWriterType>
  void
  writePolyDataFile(hid_t fileID,
                    const FieldVariablesForOutputWriterType &fieldVariables,
                    std::set<std::string> &combinedMeshesOut);

  //! write all data of all 3D or 2D field ! variables (depending on
  //! output3DMeshes) into a group. This is uses MPI IO. It can be enabled with
  //! the "combineFiles" option. on return, combinedMeshesOut contains the 3D or
  //! 2D mesh names that were written to the vtu file.
  template <typename FieldVariablesForOutputWriterType>
  void writeCombinedUnstructuredGridFile(
      hid_t fileID, const FieldVariablesForOutputWriterType &fieldVariables,
      std::set<std::string> &combinedMeshesOut, bool output3DMeshes);

protected:
  /** one Piece is the file output. It is created from one or multiple opendihu
   * meshes
   */
  struct Piece {
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
                                 // the mesh
    std::string
        firstVectorName; //< name of the first vector field variable with 3
                         // components of the mesh

    //! constructor, initialize nPoints and nCells to 0
    Piece();

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
  //! open a HDF5 file with a given filename
  hid_t openHDF5File(const char *filename, bool mpiio);

  //! helper method that writes the unstructured grid file
  template <typename FieldVariablesForOutputWriterType>
  void writeCombinedUnstructuredGridFile(
      hid_t fileID, const FieldVariablesForOutputWriterType &fieldVariables,
      PolyDataPropertiesForMesh &polyDataPropertiesForMesh,
      const std::map<std::string, PolyDataPropertiesForMesh>
          &meshPropertiesUnstructuredGridFile,
      std::vector<std::string> meshNames, bool meshPropertiesInitialized);

  bool combineFiles_; //< if set everything is combined into a single file

  std::map<std::string, PolyDataPropertiesForMesh>
      meshPropertiesPolyDataFile_; //< mesh information for a data file, for 1D
                                   // data
  std::map<std::string, PolyDataPropertiesForMesh>
      meshPropertiesUnstructuredGridFile2D_; //< mesh information for a combined
                                             // unstructured grid file,
                                             // for 2D data
  std::map<std::string, PolyDataPropertiesForMesh>
      meshPropertiesUnstructuredGridFile3D_; //< mesh information for a combined
                                             // unstructured grid file
                                             // for 3D data
  Piece piece1D_; //< the Piece data structure used for PolyDataFile, 1D
  Piece piece3D_; //< the Piece data structure used for

  int nCellsPreviousRanks1D_ = 0;  //< sum of number of cells on other processes
                                   // with lower rank no.
  int nPointsPreviousRanks1D_ = 0; //< sum of number of points on other
                                   // processes with lower rank no.

  std::map<std::string, int>
      nCellsPreviousRanks3D_; //< sum of number of cells on other processes with
                              // lower rank no.
  std::map<std::string, int>
      nPointsPreviousRanks3D_; //< sum of number of points on other processes
                               // with lower rank no.
};

namespace HDF5Utils {
//! write a key(string) value pair to a given fileID
template <typename T>
static herr_t writeAttr(hid_t fileID, const char *key, T value);
//! write a key(string) value (int32_t) pair to a given fileID
template <>
herr_t writeAttr<int32_t>(hid_t fileID, const char *key, int32_t value);
//! write a key(string) value (double) pair to a given fileID
template <>
herr_t writeAttr<double>(hid_t fileID, const char *key, double value);
//! write a key(string) value (string) pair to a given fileID
template <>
herr_t writeAttr<const std::string &>(hid_t fileID, const char *key,
                                      const std::string &value);

//! write a dataset with a specific name to a given fileID
template <typename T>
static herr_t writeSimpleVec(hid_t fileID, const std::vector<T> &data,
                             const char *dsname);
//! write a dataset(vec<int32_t>) with a specific name to a given fileID
template <>
herr_t writeSimpleVec<int32_t>(hid_t fileID, const std::vector<int32_t> &data,
                               const char *dsname);
//! write a dataset(vec<double>) with a specific name to a given fileID
template <>
herr_t writeSimpleVec<double>(hid_t fileID, const std::vector<double> &data,
                              const char *dsname);

//! write the given field variable to a given fileID
template <typename FieldVariableType>
static herr_t writeFieldVariable(hid_t fileID,
                                 FieldVariableType &fieldVariable);

//! write the a field variable indicating which ranks own which portion of the
//! domain as DataSet element to a given fileID
template <typename FieldVariableType>
static herr_t writePartitionFieldVariable(hid_t fileID,
                                          FieldVariableType &geometryField);
} // namespace HDF5Utils
} // namespace OutputWriter

#include "output_writer/hdf5/hdf5.tpp"
#include "output_writer/hdf5/hdf5_write_combined_values.tpp"
#include "output_writer/hdf5/hdf5_write_combined_file_1D.tpp"
#include "output_writer/hdf5/hdf5_write_combined_file_2D3D.tpp"
