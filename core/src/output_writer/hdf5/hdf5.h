#pragma once

#include <hdf5.h>

#include "output_writer/generic.h"
#include "output_writer/poly_data_properties_for_mesh.h"

namespace OutputWriter {

// forward declaration
namespace HDF5Utils {
class Group;
class File;
} // namespace HDF5Utils

class HDF5 : public Generic {
public:
  //! constructor
  HDF5(DihuContext context, PythonConfig specificSettings,
       std::shared_ptr<Partition::RankSubset> rankSubset = nullptr);

  //! write out solution to file, if timeStepNo is not -1, this value will be
  //! part of the filename
  template <typename DataType>
  void write(DataType &data, const char *filename = nullptr,
             int timeStepNo = -1, double currentTime = -1,
             int callCountIncrement = 1);

  //! write all 1D field variables into a group. This is uses MPI IO.
  //! It can be enabled with the "combineFiles" option. on return,
  //! combinedMeshesOut contains the 1D mesh names that were written to the hdf5
  //! file.
  template <typename FieldVariablesForOutputWriterType>
  void
  writePolyDataFile(HDF5Utils::Group &group,
                    const FieldVariablesForOutputWriterType &fieldVariables,
                    std::set<std::string> &combinedMeshesOut);

  //! write all data of all 3D or 2D field ! variables (depending on
  //! output3DMeshes) into a group. This is uses MPI IO. It can be enabled with
  //! the "combineFiles" option. on return, combinedMeshesOut contains the 3D or
  //! 2D mesh names that were written to the vtu file.
  template <typename FieldVariablesForOutputWriterType>
  void writeCombinedUnstructuredGridFile(
      HDF5Utils::Group &group,
      const FieldVariablesForOutputWriterType &fieldVariables,
      std::set<std::string> &combinedMeshesOut, bool output3DMeshes);

  //! Enable or disable combine files option
  void setCombineFiles(bool v);
  //! Enable or disable write meta option
  void setWriteMeta(bool v);

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

  //! write a vector containing nValues "12" (if output3DMeshes) or "9" (if
  //! !output3DMeshes) values for the types for an unstructured grid
  herr_t writeCombinedTypesVector(HDF5Utils::Group &group, uint64_t nValues,
                                  bool output3DMeshes, const char *dsname);

private:
  //! helper method, to write a specific variables, called by the write function
  template <typename FieldVariablesForOutputWriterType>
  void innerWrite(const FieldVariablesForOutputWriterType &variables,
                  const char *filename = nullptr, int timeStepNo = -1,
                  double currentTime = -1, int callCountIncrement = 1);

  //! helper method that writes the unstructured grid file
  template <typename FieldVariablesForOutputWriterType>
  void writeCombinedUnstructuredGridFile(
      HDF5Utils::Group &group,
      const FieldVariablesForOutputWriterType &fieldVariables,
      PolyDataPropertiesForMesh &polyDataPropertiesForMesh,
      const std::map<std::string, PolyDataPropertiesForMesh>
          &meshPropertiesUnstructuredGridFile,
      std::vector<std::string> meshNames, bool meshPropertiesInitialized);

  bool combineFiles_; //< if set everything is combined into a single file
  bool writeMeta_;    //< if set, additional metadata is written to attributes

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

//! HDF5 File abstraction, automatically opens a File and closes the File on
//! destructor call. Provides helper function for creating new Groups and
//! writing Attributes to the HDF5 file
class File {
public:
  //! Constructor opens a new file with a given filename, either with mpiio or
  //! not
  File(const char *filename, bool mpiio);
  //! Close the file and cleanup everything else still open
  ~File();

  //! Get the HDF5 File ID
  hid_t getFileID() const;
  //! Get the own rank, this value is cached
  int32_t getOwnRank() const;
  //! Get the world size, this value is cached
  int32_t getWorldSize() const;
  //! Returns true if the file was opened with MPIIO, false if not
  bool isMPIIO() const;

  //! Create a new group with a given name
  Group newGroup(const char *name) const;

  //! Write a integer attribute to the root node with a given key
  herr_t writeAttrInt(const char *key, int32_t value) const;
  //! Write a double attribute to the root node with a given key
  herr_t writeAttrDouble(const char *key, double value) const;
  //! Write a string attribute to the root node with a given key
  herr_t writeAttrStr(const char *key, const std::string &value) const;

private:
  std::string filename_; //< filename in which is been written
  const bool mpiio_;     //< stored value if the file was opened with mpiio
  hid_t fileID_;         //< stored fileID
  int32_t ownRank_;      //< own rank cached
  int32_t worldSize_;    //< world size cached
};

//! HDF5 Group abstraction that provides helper function for create another
//! nested group or writing datasets to this group
class Group {
public:
  //! Constructor for creating a new group for a given file with a given name
  Group(const File *f, const char *name);
  //! Destructor that closes the group
  ~Group();

  //! Create a new Group inside the current group, nested groups
  Group newGroup(const char *name) const;

  //! write a dataset with a specific name to the current Group
  template <typename T>
  herr_t writeSimpleVec(const std::vector<T> &data, const std::string &dsname);

private:
  //! private constructor for creating a group for a given hid_t id, so
  //! basically a already existing group. This constructor only gets called by
  //! Group::newGroup, so it can be a private constructor which will then not
  //! lead to missuse of this constructor.
  Group(const File *f, hid_t id, const char *name);

  //! inner write method, that writes a dataset with a specific name to the
  //! current group with a specific typeId and memTypeId to a mpiio file
  herr_t writeVectorMPIIO(const void *data, const std::string &dsname,
                          const size_t dataSize, hid_t typeId, hid_t memTypeId,
                          size_t dsize);

  //! inner write method, that writes a dataset with a specific name to the
  //! current group with a specific typeId and memTypeId to a regular file
  herr_t writeVector(const void *data, const std::string &dsname,
                     const size_t dataSize, hid_t typeId, hid_t memTypeId,
                     size_t dsize);

  //! write a vector of integers to a given dataset ID as an attribute with a
  //! given key
  herr_t writeAttrArray(hid_t dsetID, const char *key,
                        const std::vector<int64_t> &data) const;

  const File
      *file_; //< Const handle to a File, weak pointer because we dont cleanup
              // that file, because we dont want to call the destructor.
  hid_t groupID_; //< HDF5 ID of this Group
};

//! write the given field variable to a given group
template <typename FieldVariableType>
static herr_t writeFieldVariable(Group &group,
                                 FieldVariableType &fieldVariable);

//! write the a field variable indicating which ranks own which portion of the
//! domain as DataSet element to a given group
template <typename FieldVariableType>
static herr_t writePartitionFieldVariable(Group &group,
                                          FieldVariableType &geometryField);
} // namespace HDF5Utils
} // namespace OutputWriter

#include "output_writer/hdf5/hdf5.tpp"
#include "output_writer/hdf5/hdf5_write_combined_file_1D.tpp"
#include "output_writer/hdf5/hdf5_write_combined_file_2D3D.tpp"
