#pragma once

#include "output_writer/generic.h"
#include "output_writer/poly_data_properties_for_mesh.h"

#include <nlohmann/json.hpp>

namespace OutputWriter {

// forward declaration
namespace JsonUtils {
class Group;
class File;
} // namespace JsonUtils

class Json : public Generic {
public:
  //! constructor
  Json(DihuContext context, PythonConfig specificSettings,
       std::shared_ptr<Partition::RankSubset> rankSubset = nullptr);

  //! write out solution to file, if timeStepNo is not -1, this value will be
  //! part of the filename
  template <typename DataType>
  void write(DataType &data, const char *filename = nullptr,
             int timeStepNo = -1, double currentTime = -1,
             int callCountIncrement = 1);

  //! write all 1D field variables into a group. This is uses MPI IO.
  //! It can be enabled with the "combineFiles" option. on return,
  //! combinedMeshesOut contains the 1D mesh names that were written to the json
  //! file.
  template <typename FieldVariablesForOutputWriterType>
  void
  writePolyDataFile(JsonUtils::Group &group,
                    const FieldVariablesForOutputWriterType &fieldVariables,
                    std::set<std::string> &combinedMeshesOut,
                    bool useUniqueName = false);

  //! write all data of all 3D or 2D field ! variables (depending on
  //! output3DMeshes) into a group. This is uses MPI IO. It can be enabled with
  //! the "combineFiles" option. on return, combinedMeshesOut contains the 3D or
  //! 2D mesh names that were written to the vtu file.
  template <typename FieldVariablesForOutputWriterType>
  void writeCombinedUnstructuredGridFile(
      JsonUtils::Group &group,
      const FieldVariablesForOutputWriterType &fieldVariables,
      std::set<std::string> &combinedMeshesOut, bool output3DMeshes,
      bool useUniqueName = false);

  //! Enable or disable combine files option
  void setCombineFiles(bool v);
  //! Enable or disable write meta option
  void setWriteMeta(bool v);
  //! Enable or disable use checkpointing data option
  void setUseCheckpointData(bool v);

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
  void writeCombinedTypesVector(JsonUtils::Group &group, uint64_t nValues,
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
      JsonUtils::Group &group,
      const FieldVariablesForOutputWriterType &fieldVariables,
      PolyDataPropertiesForMesh &polyDataPropertiesForMesh,
      const std::map<std::string, PolyDataPropertiesForMesh>
          &meshPropertiesUnstructuredGridFile,
      std::vector<std::string> meshNames, bool meshPropertiesInitialized,
      bool useUniqueName = false);

  bool combineFiles_; //< if set everything is combined into a single file
  bool writeMeta_;    //< if set, additional metadata is written to attributes
  bool useCheckpointData_; //< if set the output writer uses
                           // getFieldVariablesForCheckpointing rather than
                           // getFieldVariablesForOutputWriter

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

namespace JsonUtils {

//! Json File abstraction, automatically opens a File and closes the File on
//! destructor call. Provides helper function for creating new Groups and
//! writing Attributes to the Json file
class File {
public:
  //! Constructor opens a new file with a given filename, either with mpiio or
  //! not
  File(const char *filename, bool mpiio);
  //! Close the file and cleanup everything else still open
  ~File();

  //! Get the Json File Content
  nlohmann::json &getFileContent();
  //! Get the own rank, this value is cached
  int32_t getOwnRank() const;
  //! Get the world size, this value is cached
  int32_t getWorldSize() const;
  //! Returns true if the file was opened with MPIIO, false if not
  bool isMPIIO() const;

  //! Create a new group with a given name
  Group newGroup(const char *name);

  //! Write a integer attribute to the root node with a given key
  template <typename T,
            std::enable_if_t<std::is_same<T, int32_t>::value, bool> = true>
  void writeAttr(const char *key, const T &v);
  //! Write a double attribute to the root node with a given key
  template <typename T,
            std::enable_if_t<std::is_same<T, double>::value, bool> = true>
  void writeAttr(const char *key, const T &v);
  //! Write a string attribute to the root node with a given key
  template <typename T,
            std::enable_if_t<std::is_same<T, std::string>::value, bool> = true>
  void writeAttr(const char *key, const T &v);

private:
  std::string filename_;   //< filename in which is been written
  const bool mpiio_;       //< stored value if the file was opened with mpiio
  nlohmann::json content_; //< json content
  int32_t ownRank_;        //< own rank cached
  int32_t worldSize_;      //< world size cached
};

//! Json Group abstraction that provides helper function for create another
//! nested group or writing datasets to this group
class Group {
public:
  //! Constructor for creating a new group for a given file with a given name
  Group(File *file, const char *name);
  //! Destructor that closes the group
  ~Group() = default;

  //! Create a new Group inside the current group, nested groups
  Group newGroup(const char *name);

  //! write a dataset with a specific name to the current Group
  template <typename T>
  void writeSimpleVec(const std::vector<T> &data, const std::string &dsname);

private:
  //! Constructor for creating a new group for a given file with a given name
  //! and a specific json object
  Group(File *file, nlohmann::json &obj, const char *name);

  //! inner write method, that writes a dataset with a specific name to the
  //! current group with a specific typeId and memTypeId to a mpiio file
  template <typename T>
  void writeVectorMPIIO(const std::vector<T> &data, const std::string &dsname);

  //! inner write method, that writes a dataset with a specific name to the
  //! current group with a specific typeId and memTypeId to a regular file
  template <typename T>
  void writeVector(const std::vector<T> &data, const std::string &dsname);

  File *file_; //< Handle to a File, weak pointer because we dont cleanup
               // that file, because we dont want to call the destructor.
  nlohmann::json &groupContent_; //< json obj of this Group, needs to be a
                                 // reference to link back to original obj
};

//! write the given field variable to a given group
template <typename FieldVariableType>
static void writeFieldVariable(Group &group, FieldVariableType &fieldVariable);

//! write the a field variable indicating which ranks own which portion of the
//! domain as DataSet element to a given group
template <typename FieldVariableType>
static void writePartitionFieldVariable(Group &group,
                                        FieldVariableType &geometryField);
} // namespace JsonUtils
} // namespace OutputWriter

#include "output_writer/json/json.tpp"
#include "output_writer/json/json_write_combined_file_1D.tpp"
#include "output_writer/json/json_write_combined_file_2D3D.tpp"
