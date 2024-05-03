#pragma once

#include <Python.h> // has to be the first included header

#include "specialized_solver/muscle_contraction_solver.h"

#include "time_stepping_scheme/crank_nicolson.h"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver.h"

#include "precice/precice.hpp"
#include "control/precice/initialize.h"

namespace Control {

// Forward declaration, required for the anonymous enums
template <typename NestedSolver> class PreciceAdapterInitialize;

// template <typename NestedSolver>
// struct PreciceAdapterInitialize<NestedSolver>::PreciceSurfaceData;
/** This is a base class of the precice adapter that contains functionality that
 * depends on the type of the nested solver. All solvers that should be able to
 * use precice surface coupling have to implement this interface.
 */

template <typename NestedSolver>
class PreciceAdapterNestedSolver : public Runnable {};

/** Partial specialization for tendon or pure mechanics solver in a coupling
 * scheme, muscle contraction solver (nonlinear elasticity with active stress)
 */
template <typename T1>
class PreciceAdapterNestedSolver<FastMonodomainSolver<T1>> {
public:
  std::vector<double> displacementValues_;
  std::vector<double> velocityValues_;
  std::vector<double> tractionValues_;
  std::vector<Vec3> displacementVectors_;
  std::vector<Vec3> velocityVectors_;
  std::vector<Vec3> tractionVectors_;
  std::vector<double> scalarValues_;
  std::vector<double> scalarValuesOfMesh_;
  std::vector<Vec3> geometryValues_;
  //! define the type of the nested solver
  typedef FastMonodomainSolver<T1> NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolverType::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<FunctionSpace> functionSpace(NestedSolverType &nestedSolver);

  //! parse the options in "preciceVolumeData" and initialize all variables in
  //! precice, store in variable preciceVolumeData_
  template <typename PreciceVolumeData, typename PreciceVolumeMesh>
  void initializePreciceVolumeData(
      PythonConfig &specificSettings, NestedSolverType &nestedSolver,
      std::shared_ptr<precice::Participant> &preciceParticipant,
      std::vector<PreciceVolumeData> &preciceVolumeData,
      std::vector<std::shared_ptr<PreciceVolumeMesh>> &preciceVolumeMeshes) {
    // parse settings for coupling participants
    // loop over items of the key "preciceData"
    std::string settingsKey("preciceVolumeData");
    PyObject *listPy = specificSettings.getOptionPyObject(settingsKey);
    std::vector<PyObject *> list =
        PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);
    PythonConfig preciceDataConfig(specificSettings, settingsKey);

    // get the slot connector data from the nested solver
    int slotNo = 0;
    using SlotConnectorDataType =
        typename NestedSolverType::SlotConnectorDataType;
    std::shared_ptr<SlotConnectorDataType> slotConnectorData =
        nestedSolver.getSlotConnectorData();

    // get all present slot names
    std::vector<std::string> slotNames;
    SlotConnectorDataHelper<SlotConnectorDataType>::getSlotNames(
        slotConnectorData, slotNames);

    // loop over items of the list under "preciceData", the order of the items
    // is also the order of the data connector slots
    for (int listIndex = 0; listIndex < list.size(); listIndex++) {
      // extract the list item as PythonConfig {"mode": "write", "variableName":
      // ...}
      PythonConfig currentPreciceData(preciceDataConfig, listIndex);

      // create a new preciceData instance that will be added to
      // preciceVolumeData_ at the end of this loop
      PreciceVolumeData preciceData;

      // parse options
      preciceData.isGeometryField =
          currentPreciceData.getOptionBool("isGeometryField", false);
      preciceData.opendihuMeshName =
          currentPreciceData.getOptionString("opendihuMeshName", "");

      // parse slot name
      preciceData.slotName = currentPreciceData.getOptionString("slotName", "");
      preciceData.slotNo = slotNo;

      LOG(DEBUG) << "isGeometryField: " << preciceData.isGeometryField
                 << ", slotName: \"" << preciceData.slotName << "\".";

      if (preciceData.slotName == "") {
        preciceData.slotNo = slotNo;
        LOG(DEBUG) << "Using slotNo " << slotNo
                   << " because slotName was not given.";
      } else {
        LOG(DEBUG) << "The slotName \"" << preciceData.slotName
                   << "\" is given.";
        std::vector<std::string>::iterator slotNameIter =
            std::find(slotNames.begin(), slotNames.end(), preciceData.slotName);

        // if the "from" and "to" slot names match to slot names in the current
        // slotConnectorData
        if (slotNameIter != slotNames.end()) {
          LOG(DEBUG) << "Slot found";
          // determine slot nos
          preciceData.slotNo = std::distance(slotNames.begin(), slotNameIter);
          LOG(DEBUG) << "Slot \"" << preciceData.slotName << "\" is slot no "
                     << preciceData.slotNo;
        } else {
          LOG(DEBUG) << "Slot not found";
          preciceData.slotNo = slotNo;
          LOG(WARNING)
              << "A slot with name \"" << preciceData.slotName
              << "\" does not exist. Available slots: " << slotNames
              << ". Use `None` or the empty string to specify the slot by its "
                 "number (referenced slots will be in the order of "
              << "the items under \"preciceData\").\n Using slot No " << slotNo;
        }
      }

      // parse the mesh
      std::string currentMeshName =
          currentPreciceData.getOptionString("preciceMeshName", "");

      // find the mesh in the already parsed precice meshes
      typename std::vector<std::shared_ptr<PreciceVolumeMesh>>::iterator iter =
          std::find_if(preciceVolumeMeshes.begin(), preciceVolumeMeshes.end(),
                       [&currentMeshName](
                           std::shared_ptr<PreciceVolumeMesh> preciceMesh) {
                         return preciceMesh->preciceMeshName == currentMeshName;
                       });

      // if the mesh is not in preciceVolumeMeshes, create it and add it to
      // preciceVolumeMeshes
      if (iter == preciceVolumeMeshes.end()) {
        //   LOG(DEBUG) << "mesh was not yet initialized.";

        // create new precice mesh object
        std::shared_ptr<PreciceVolumeMesh> preciceMesh =
            std::make_shared<PreciceVolumeMesh>();
        preciceMesh->preciceMeshName = currentMeshName;

        std::vector<Vec3> geometryValues;
        int nArrayItems = 1;

        // get the function space that is used to initialize the mapping
        std::shared_ptr<typename NestedSolverType::FunctionSpace>
            functionSpace = nullptr;
        if (preciceData.opendihuMeshName != "") {
          functionSpace =
              DihuContext::meshManager()
                  ->functionSpace<typename NestedSolverType::FunctionSpace>(
                      preciceData.opendihuMeshName);
          LOG(DEBUG) << "Using opendihu mesh with name \""
                     << preciceData.opendihuMeshName
                     << "\" to initialize precice mapping.";

          // get the node positions of the opendihu mesh to initialize the
          // mapping
          preciceMesh->nNodesLocal =
              functionSpace->geometryField().nDofsLocalWithoutGhosts();
          functionSpace->geometryField().getValuesWithoutGhosts(geometryValues);

          // store opendihu mesh name
          preciceMesh->opendihuMeshName = functionSpace->meshName();

          // the following cannot happen
          if (preciceMesh->opendihuMeshName != preciceData.opendihuMeshName)
            LOG(FATAL)
                << "Initializing precice mesh from settings opendihuMeshName=\""
                << preciceData.opendihuMeshName
                << "\" and the resulting mesh is "
                << preciceMesh->opendihuMeshName;
        } else {
          LOG(DEBUG) << "get mesh from mesh partition, slot No "
                     << preciceData.slotNo;
          LOG(DEBUG) << "slot connector data type: "
                     << StringUtility::demangle(
                            typeid(SlotConnectorDataType).name());

          // get the mesh partition
          std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
              SlotConnectorDataHelper<SlotConnectorDataType>::
                  getMeshPartitionBase(slotConnectorData, preciceData.slotNo,
                                       0);

          if (!meshPartitionBase) {
            LOG(FATAL) << "Could not get mesh for slot No "
                       << preciceData.slotNo;
          } else
            LOG(DEBUG) << "got mesh partition for slot No "
                       << preciceData.slotNo;

          // get opendihu mesh name
          preciceMesh->opendihuMeshName =
              SlotConnectorDataHelper<SlotConnectorDataType>::getMeshName(
                  slotConnectorData, preciceData.slotNo);

          int nDofsLocalWithoutGhosts =
              meshPartitionBase->nDofsLocalWithoutGhosts();
          nArrayItems =
              SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
                  slotConnectorData,
                  preciceData.slotNo); // number of fibers if there are fibers

          preciceMesh->nNodesLocal = nDofsLocalWithoutGhosts * nArrayItems;

          // get the vector of values [0,1,...,nDofsLocalWithGhosts]
          const std::vector<PetscInt> &dofNosLocalWithGhosts =
              meshPartitionBase->dofNosLocal();
          std::vector<PetscInt> dofNosLocalWithoutGhosts(
              dofNosLocalWithGhosts.begin(),
              dofNosLocalWithGhosts.begin() + nDofsLocalWithoutGhosts);

          // loop over fibers if there are any
          for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
            static std::vector<Vec3> nodePositionsFiber;
            nodePositionsFiber.clear();
            SlotConnectorDataHelper<SlotConnectorDataType>::
                slotGetGeometryValues(slotConnectorData, preciceData.slotNo,
                                      arrayIndex, dofNosLocalWithoutGhosts,
                                      nodePositionsFiber);
            geometryValues.insert(geometryValues.end(),
                                  nodePositionsFiber.begin(),
                                  nodePositionsFiber.end());
          }

          LOG(DEBUG) << "collected " << geometryValues.size()
                     << " node positions from the " << nArrayItems << " fibers";
        }

        // transform to contiguous memory layout for precice
        std::vector<double> geometryValuesContiguous(3 * geometryValues.size());
        for (int entryNo = 0; entryNo < geometryValues.size(); entryNo++)
          for (int componentNo = 0; componentNo < 3; componentNo++)
            geometryValuesContiguous[3 * entryNo + componentNo] =
                geometryValues[entryNo][componentNo];

        if (geometryValuesContiguous.size() != preciceMesh->nNodesLocal * 3)
          LOG(FATAL) << "size mismatch: " << geometryValuesContiguous.size()
                     << "!=" << preciceMesh->nNodesLocal * 3;

        // resize buffer for vertex ids
        preciceMesh->preciceVertexIds.resize(preciceMesh->nNodesLocal);

        // give the node positions to precice and get the vertex ids
        preciceParticipant->setMeshVertices(preciceMesh->preciceMeshName,
                                            geometryValuesContiguous,
                                            preciceMesh->preciceVertexIds);

        // store newly created mesh in preciceData
        preciceData.preciceMesh = preciceMesh;

        // store the precice mesh to the vector of meshes
        preciceVolumeMeshes.push_back(preciceMesh);

        // output message about mesh names
        std::stringstream message;
        message << "Initialized precice mesh \"" << preciceMesh->preciceMeshName
                << "\" from opendihu mesh \"" << preciceMesh->opendihuMeshName
                << "\" with " << preciceMesh->nNodesLocal << " local nodes";
        if (nArrayItems > 1)
          message << " (and other corresponding meshes, in total "
                  << nArrayItems << " fibers/compartments)";
        LOG(INFO) << message.str() << ".";
      } else {
        LOG(DEBUG) << "Use existing precice mesh " << currentMeshName;
        preciceData.preciceMesh = *iter;
      }

      // parse mode
      std::string mode = currentPreciceData.getOptionString("mode", "");
      if (mode == "read") {
        preciceData.ioType = PreciceAdapterInitialize<
            NestedSolverType>::PreciceVolumeData::ioRead;
      } else if (mode == "write") {
        preciceData.ioType = PreciceAdapterInitialize<
            NestedSolverType>::PreciceVolumeData::ioWrite;
      } else {
        LOG(FATAL) << currentPreciceData << "[\"mode\"] is \"" << mode << "\", "
                   << "possible values are: \"read\", \"write\".";
      }

      // parse variable name
      preciceData.preciceDataName =
          currentPreciceData.getOptionString("preciceDataName", "variable");

      // increment slotNo, this value is used for slots that are not identified
      // by slotName in the config
      if (!preciceData.isGeometryField)
        slotNo++;

      // check if mesh size matches the specified slot
      // get the mesh partition
      std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
          SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase(
              slotConnectorData, preciceData.slotNo, 0);

      int nDofsLocalWithoutGhosts =
          meshPartitionBase->nDofsLocalWithoutGhosts();
      int nArrayItems =
          SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
              slotConnectorData,
              preciceData.slotNo); // number of fibers if there are fibers

      if (nDofsLocalWithoutGhosts * nArrayItems !=
          preciceData.preciceMesh->nNodesLocal) {
        LOG(DEBUG) << ", all available slots: "
                   << SlotConnectorDataHelper<SlotConnectorDataType>::getString(
                          slotConnectorData);
        LOG(FATAL)
            << currentPreciceData
            << ": Mesh does not match slot in "
               "PreciceAdapterVolumeCoupling.\n\n "
            << "The " << (listIndex + 1)
            << (listIndex == 0
                    ? "st"
                    : (listIndex == 1 ? "nd" : (listIndex == 2 ? "rd" : "th")))
            << " list item under \"preciceData\" uses slot "
            << preciceData.slotNo << " (slotName \"" << preciceData.slotName
            << "\")"
            << " and preciceMesh \"" << preciceData.preciceMesh->preciceMeshName
            << "\".\nThe slot has " << nDofsLocalWithoutGhosts << " dofs * "
            << nArrayItems << " array items (fibers) = "
            << nDofsLocalWithoutGhosts * nArrayItems
            << " dofs in total. The mesh was initialized from opendihu mesh \""
            << preciceData.preciceMesh->opendihuMeshName << "\" with "
            << preciceData.preciceMesh->nNodesLocal << " nodes. \n"
            << "(isGeometryField: " << std::boolalpha
            << preciceData.isGeometryField << ")\n"
            << "You can do the following:\n"
            << "- Rename the preciceMesh to a unique name and update the "
               "precice "
               "config xml file accordingly.\n"
            << "- Check the solver structure file to find out the meshes that "
               "are associated with the slots.\n"
            << "  Check that the correct mesh is specified and that the "
               "precice "
               "mesh was not initialized earlier with a different opendihu "
               "mesh.\n"
            << "  Every precice mesh is initialized the first time in appears "
               "under \"preciceData\", sometimes reordering the entries can "
               "help.\n"
            << "  If the desired mesh is not available at the current solver, "
               "maybe insert a MapDofs class.";
      }

      LOG(INFO) << "Precice data \"" << preciceData.preciceDataName
                << "\" maps to "
                << (preciceData.isGeometryField ? "the geometry field of " : "")
                << "slot " << preciceData.slotNo << " (\""
                << preciceData.slotName << "\") and uses precice mesh \""
                << preciceData.preciceMesh->preciceMeshName
                << "\", which is opendihu mesh \""
                << preciceData.preciceMesh->opendihuMeshName << "\" with "
                << preciceData.preciceMesh->nNodesLocal << " local nodes.";

      // store preciceData to vector
      preciceVolumeData.push_back(preciceData);
    }
  };

  template <typename SurfaceDataVector, typename VolumeDataVector>
  void
  preciceReadData(NestedSolverType &nestedSolver,
                  std::shared_ptr<precice::Participant> &preciceParticipant,
                  SurfaceDataVector &preciceSurfaceData,
                  VolumeDataVector &preciceVolumeData, DihuContext &context) {
    LOG(DEBUG) << "read volume data from precice";
    double preciceDt = preciceParticipant->getMaxTimeStepSize();

    using SlotConnectorDataType =
        typename NestedSolverType::SlotConnectorDataType;
    std::shared_ptr<SlotConnectorDataType> slotConnectorData =
        nestedSolver.getSlotConnectorData();

    // loop over data
    for (auto &preciceData : preciceVolumeData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceVolumeData::ioRead) {
        int nEntries = preciceData.preciceMesh->nNodesLocal;

        if (preciceData.isGeometryField) {
          nEntries = preciceData.preciceMesh->nNodesLocal * 3;
        }

        // allocate temporary memory
        scalarValues_.resize(nEntries);

        // get all data at once
        preciceParticipant->readData(preciceData.preciceMesh->preciceMeshName,
                                     preciceData.preciceDataName,
                                     preciceData.preciceMesh->preciceVertexIds,
                                     preciceDt, scalarValues_);

        // get the mesh partition
        std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
            SlotConnectorDataHelper<
                SlotConnectorDataType>::getMeshPartitionBase(slotConnectorData,
                                                             preciceData.slotNo,
                                                             0);

        int nDofsLocalWithoutGhosts =
            meshPartitionBase->nDofsLocalWithoutGhosts();

        // get the vector of values [0,1,...,nDofsLocalWithGhosts]
        const std::vector<PetscInt> &dofNosLocalWithGhosts =
            meshPartitionBase->dofNosLocal();
        std::vector<PetscInt> dofNosLocalWithoutGhosts(
            dofNosLocalWithGhosts.begin(),
            dofNosLocalWithGhosts.begin() + nDofsLocalWithoutGhosts);

        int nArrayItems =
            SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
                slotConnectorData,
                preciceData.slotNo); // number of fibers if there are fibers

        // store received data in field variable
        if (preciceData.isGeometryField) {
          // loop over fibers if there are any
          for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
            // fill the vector geometryValues_ with the geometry values of the
            // current fiber or mesh
            geometryValues_.resize(nDofsLocalWithoutGhosts);
            for (int dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts;
                 dofNoLocal++) {
              for (int componentNo = 0; componentNo < 3; componentNo++) {
                geometryValues_[dofNoLocal][componentNo] =
                    scalarValues_[3 * (arrayIndex * nDofsLocalWithoutGhosts +
                                       dofNoLocal) +
                                  componentNo];
                //
                std::cout << geometryValues_[dofNoLocal][componentNo]
                          << std::endl;
              }
            }

            SlotConnectorDataHelper<SlotConnectorDataType>::
                slotSetGeometryValues(slotConnectorData, preciceData.slotNo,
                                      arrayIndex, dofNosLocalWithoutGhosts,
                                      geometryValues_);
          }
        } else {
          // loop over fibers if there are any
          for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
            // fill the vector geometryValues_ with the geometry values of the
            // current fiber or mesh
            scalarValuesOfMesh_.assign(
                scalarValues_.begin() + arrayIndex * nDofsLocalWithoutGhosts,
                scalarValues_.begin() +
                    (arrayIndex + 1) * nDofsLocalWithoutGhosts);

            SlotConnectorDataHelper<SlotConnectorDataType>::slotSetValues(
                slotConnectorData, preciceData.slotNo, arrayIndex,
                dofNosLocalWithoutGhosts, scalarValuesOfMesh_);
          }
        }
      }
    }
  }

  template <typename SurfaceDataVector, typename VolumeDataVector>
  void
  preciceWriteData(NestedSolverType &nestedSolver,
                   std::shared_ptr<precice::Participant> &preciceParticipant,
                   SurfaceDataVector &preciceSurfaceData,
                   VolumeDataVector &preciceVolumeData, double scalingFactor) {
    // write data to precice
    LOG(DEBUG) << "write volume data to precice";

    using SlotConnectorDataType =
        typename NestedSolverType::SlotConnectorDataType;
    std::shared_ptr<SlotConnectorDataType> slotConnectorData =
        nestedSolver.getSlotConnectorData();

    // loop over data
    for (auto &preciceData : preciceVolumeData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceVolumeData::ioWrite) {
        scalarValues_.clear();

        // get the mesh partition
        std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
            SlotConnectorDataHelper<
                SlotConnectorDataType>::getMeshPartitionBase(slotConnectorData,
                                                             preciceData.slotNo,
                                                             0);

        int nDofsLocalWithoutGhosts =
            meshPartitionBase->nDofsLocalWithoutGhosts();

        // get the vector of values [0,1,...,nDofsLocalWithGhosts]
        const std::vector<PetscInt> &dofNosLocalWithGhosts =
            meshPartitionBase->dofNosLocal();
        std::vector<PetscInt> dofNosLocalWithoutGhosts(
            dofNosLocalWithGhosts.begin(),
            dofNosLocalWithGhosts.begin() + nDofsLocalWithoutGhosts);

        int nArrayItems =
            SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
                slotConnectorData,
                preciceData.slotNo); // number of fibers if there are fibers

        // if it is a geometry field, get the node positions of a mesh
        if (preciceData.isGeometryField) {
          geometryValues_.clear();

          // loop over fibers if there are any
          for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
            static std::vector<Vec3> geometryValuesFiber;
            geometryValuesFiber.clear();
            SlotConnectorDataHelper<SlotConnectorDataType>::
                slotGetGeometryValues(slotConnectorData, preciceData.slotNo,
                                      arrayIndex, dofNosLocalWithoutGhosts,
                                      geometryValuesFiber);

            geometryValues_.insert(geometryValues_.end(),
                                   geometryValuesFiber.begin(),
                                   geometryValuesFiber.end());
          }
          LOG(DEBUG) << "Using geometry field of opendihu meshes of "
                     << nArrayItems << " array items (e.g., fibers).";

          // transform to contiguous memory layout for precice
          scalarValues_.resize(3 * geometryValues_.size());

          for (int entryNo = 0; entryNo < geometryValues_.size(); entryNo++)
            for (int componentNo = 0; componentNo < 3; componentNo++)
              scalarValues_[3 * entryNo + componentNo] =
                  geometryValues_[entryNo][componentNo];

        } else {
          // the data is a normal slot, no geometry field

          // loop over fibers if there are any
          for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
            static std::vector<double> values;
            values.clear();
            // static void slotGetValues(std::shared_ptr<SlotConnectorDataType>
            // slotConnectorData,
            //   int slotNo, int arrayIndex, const std::vector<dof_no_t>
            //   &dofNosLocal, std::vector<double> &values);
            SlotConnectorDataHelper<SlotConnectorDataType>::slotGetValues(
                slotConnectorData, preciceData.slotNo, arrayIndex,
                dofNosLocalWithoutGhosts, values);
            scalarValues_.insert(scalarValues_.end(), values.begin(),
                                 values.end());

            LOG(DEBUG) << "arrayIndex " << arrayIndex << ", add "
                       << values.size()
                       << " values, now number: " << scalarValues_.size();
          }
        }

        // scale the values by a factor given in the config
        for (double &value : scalarValues_)
          value *= scalingFactor;

        // check dim for scalarValues_
        if (preciceData.isGeometryField) {
          if (scalarValues_.size() != 3 * preciceData.preciceMesh->nNodesLocal)
            LOG(FATAL) << "Wrong number of scalar values (isGeometryField): "
                       << scalarValues_.size()
                       << " != " << preciceData.preciceMesh->nNodesLocal
                       << ", nArrayItems: " << nArrayItems;
        } else {
          if (scalarValues_.size() != preciceData.preciceMesh->nNodesLocal)
            LOG(FATAL) << "Wrong number of scalar values: "
                       << scalarValues_.size()
                       << " != " << preciceData.preciceMesh->nNodesLocal
                       << ", nArrayItems: " << nArrayItems;
        }

        // write values in precice
        preciceParticipant->writeData(preciceData.preciceMesh->preciceMeshName,
                                      preciceData.preciceDataName,
                                      preciceData.preciceMesh->preciceVertexIds,
                                      scalarValues_);
      }
    }
    LOG(DEBUG) << "write volume data to precice complete";
  }

  //! initialize dirichlet boundary conditions by adding new dofs and prescribed
  //! values for all bottom or top nodes
  void addDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<
          typename SpatialDiscretization::DirichletBoundaryConditionsBase<
              FunctionSpace, 6>::ElementWithNodes>
          &dirichletBoundaryConditionElements);

  //! update existing boundary conditions with new values
  void updateDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<std::pair<global_no_t, std::array<double, 6>>>
          newDirichletBoundaryConditionValues);

  //! update the neumann boundary conditions by replacing the complete object
  void updateNeumannBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>>
          neumannBoundaryConditions);

  //! get the displacement and velocity vectors of the given local dof nos
  void getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                     const std::vector<dof_no_t> &dofNosLocal,
                                     std::vector<double> &displacementValues,
                                     std::vector<double> &velocityValues);

  //! get the traction vectors of the given local dof nos
  void getTractionValues(NestedSolverType &nestedSolver,
                         const std::vector<dof_no_t> &dofNosLocal,
                         std::vector<double> &tractionValues);

  //! get at Petsc Vec that stores all values of the current state, to be used
  //! to store and restore checkpoints
  Vec currentState(NestedSolverType &nestedSolver);

  //! get the field variable of the deformation gradient
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 9>>
  deformationGradientField(NestedSolverType &nestedSolver);

  void reset(NestedSolverType &nestedSolver);

protected:
  void saveFiberData(NestedSolverType &nestedSolver);

  void loadFiberData(NestedSolverType &nestedSolver);
};

/** Partial specialization for tendon or pure mechanics solver in a coupling
 * scheme, muscle contraction solver (nonlinear elasticity with active stress)
 */
template <typename T1, typename T2>
class PreciceAdapterNestedSolver<MuscleContractionSolver<T1, T2>> {
public:
  std::vector<double> displacementValues_;
  std::vector<double> velocityValues_;
  std::vector<double> tractionValues_;
  std::vector<Vec3> displacementVectors_;
  std::vector<Vec3> velocityVectors_;
  std::vector<Vec3> tractionVectors_;
  std::vector<double> scalarValues_;
  std::vector<double> scalarValuesOfMesh_;
  std::vector<Vec3> geometryValues_;
  //! define the type of the nested solver
  typedef MuscleContractionSolver<T1, T2> NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolverType::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<FunctionSpace> functionSpace(NestedSolverType &nestedSolver);
  //! parse the options in "preciceVolumeData" and initialize all variables in
  //! precice, store in variable preciceVolumeData_
  template <typename PreciceVolumeData, typename PreciceVolumeMesh>
  void initializePreciceVolumeData(
      PythonConfig &specificSettings, NestedSolverType &nestedSolver,
      std::shared_ptr<precice::Participant> &preciceParticipant,
      std::vector<PreciceVolumeData> &preciceVolumeData,
      std::vector<std::shared_ptr<PreciceVolumeMesh>> &preciceVolumeMeshes) {
    // parse settings for coupling participants
    // loop over items of the key "preciceData"
    std::string settingsKey("preciceVolumeData");
    PyObject *listPy = specificSettings.getOptionPyObject(settingsKey);
    std::vector<PyObject *> list =
        PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);
    PythonConfig preciceDataConfig(specificSettings, settingsKey);

    // get the slot connector data from the nested solver
    int slotNo = 0;
    using SlotConnectorDataType =
        typename NestedSolverType::SlotConnectorDataType;
    std::shared_ptr<SlotConnectorDataType> slotConnectorData =
        nestedSolver.getSlotConnectorData();

    // get all present slot names
    std::vector<std::string> slotNames;
    SlotConnectorDataHelper<SlotConnectorDataType>::getSlotNames(
        slotConnectorData, slotNames);

    // loop over items of the list under "preciceData", the order of the items
    // is also the order of the data connector slots
    for (int listIndex = 0; listIndex < list.size(); listIndex++) {
      // extract the list item as PythonConfig {"mode": "write", "variableName":
      // ...}
      PythonConfig currentPreciceData(preciceDataConfig, listIndex);

      // create a new preciceData instance that will be added to
      // preciceVolumeData_ at the end of this loop
      PreciceVolumeData preciceData;

      // parse options
      preciceData.isGeometryField =
          currentPreciceData.getOptionBool("isGeometryField", false);
      preciceData.opendihuMeshName =
          currentPreciceData.getOptionString("opendihuMeshName", "");

      // parse slot name
      preciceData.slotName = currentPreciceData.getOptionString("slotName", "");
      preciceData.slotNo = slotNo;

      LOG(DEBUG) << "isGeometryField: " << preciceData.isGeometryField
                 << ", slotName: \"" << preciceData.slotName << "\".";

      if (preciceData.slotName == "") {
        preciceData.slotNo = slotNo;
        LOG(DEBUG) << "Using slotNo " << slotNo
                   << " because slotName was not given.";
      } else {
        LOG(DEBUG) << "The slotName \"" << preciceData.slotName
                   << "\" is given.";
        std::vector<std::string>::iterator slotNameIter =
            std::find(slotNames.begin(), slotNames.end(), preciceData.slotName);

        // if the "from" and "to" slot names match to slot names in the current
        // slotConnectorData
        if (slotNameIter != slotNames.end()) {
          LOG(DEBUG) << "Slot found";
          // determine slot nos
          preciceData.slotNo = std::distance(slotNames.begin(), slotNameIter);
          LOG(DEBUG) << "Slot \"" << preciceData.slotName << "\" is slot no "
                     << preciceData.slotNo;
        } else {
          LOG(DEBUG) << "Slot not found";
          preciceData.slotNo = slotNo;
          LOG(WARNING)
              << "A slot with name \"" << preciceData.slotName
              << "\" does not exist. Available slots: " << slotNames
              << ". Use `None` or the empty string to specify the slot by its "
                 "number (referenced slots will be in the order of "
              << "the items under \"preciceData\").\n Using slot No " << slotNo;
        }
      }

      // parse the mesh
      std::string currentMeshName =
          currentPreciceData.getOptionString("preciceMeshName", "");

      // find the mesh in the already parsed precice meshes
      typename std::vector<std::shared_ptr<PreciceVolumeMesh>>::iterator iter =
          std::find_if(preciceVolumeMeshes.begin(), preciceVolumeMeshes.end(),
                       [&currentMeshName](
                           std::shared_ptr<PreciceVolumeMesh> preciceMesh) {
                         return preciceMesh->preciceMeshName == currentMeshName;
                       });

      // if the mesh is not in preciceVolumeMeshes, create it and add it to
      // preciceVolumeMeshes
      if (iter == preciceVolumeMeshes.end()) {
        //   LOG(DEBUG) << "mesh was not yet initialized.";

        // create new precice mesh object
        std::shared_ptr<PreciceVolumeMesh> preciceMesh =
            std::make_shared<PreciceVolumeMesh>();
        preciceMesh->preciceMeshName = currentMeshName;

        std::vector<Vec3> geometryValues;
        int nArrayItems = 1;

        // get the function space that is used to initialize the mapping
        std::shared_ptr<typename NestedSolverType::FunctionSpace>
            functionSpace = nullptr;
        if (preciceData.opendihuMeshName != "") {
          functionSpace =
              DihuContext::meshManager()
                  ->functionSpace<typename NestedSolverType::FunctionSpace>(
                      preciceData.opendihuMeshName);
          LOG(DEBUG) << "Using opendihu mesh with name \""
                     << preciceData.opendihuMeshName
                     << "\" to initialize precice mapping.";

          // get the node positions of the opendihu mesh to initialize the
          // mapping
          preciceMesh->nNodesLocal =
              functionSpace->geometryField().nDofsLocalWithoutGhosts();
          functionSpace->geometryField().getValuesWithoutGhosts(geometryValues);

          // store opendihu mesh name
          preciceMesh->opendihuMeshName = functionSpace->meshName();

          // the following cannot happen
          if (preciceMesh->opendihuMeshName != preciceData.opendihuMeshName)
            LOG(FATAL)
                << "Initializing precice mesh from settings opendihuMeshName=\""
                << preciceData.opendihuMeshName
                << "\" and the resulting mesh is "
                << preciceMesh->opendihuMeshName;
        } else {
          LOG(DEBUG) << "get mesh from mesh partition, slot No "
                     << preciceData.slotNo;
          LOG(DEBUG) << "slot connector data type: "
                     << StringUtility::demangle(
                            typeid(SlotConnectorDataType).name());

          // get the mesh partition
          std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
              SlotConnectorDataHelper<SlotConnectorDataType>::
                  getMeshPartitionBase(slotConnectorData, preciceData.slotNo,
                                       0);

          if (!meshPartitionBase) {
            LOG(FATAL) << "Could not get mesh for slot No "
                       << preciceData.slotNo;
          } else
            LOG(DEBUG) << "got mesh partition for slot No "
                       << preciceData.slotNo;

          // get opendihu mesh name
          preciceMesh->opendihuMeshName =
              SlotConnectorDataHelper<SlotConnectorDataType>::getMeshName(
                  slotConnectorData, preciceData.slotNo);

          int nDofsLocalWithoutGhosts =
              meshPartitionBase->nDofsLocalWithoutGhosts();
          nArrayItems =
              SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
                  slotConnectorData,
                  preciceData.slotNo); // number of fibers if there are fibers

          preciceMesh->nNodesLocal = nDofsLocalWithoutGhosts * nArrayItems;

          // get the vector of values [0,1,...,nDofsLocalWithGhosts]
          const std::vector<PetscInt> &dofNosLocalWithGhosts =
              meshPartitionBase->dofNosLocal();
          std::vector<PetscInt> dofNosLocalWithoutGhosts(
              dofNosLocalWithGhosts.begin(),
              dofNosLocalWithGhosts.begin() + nDofsLocalWithoutGhosts);

          // loop over fibers if there are any
          for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
            static std::vector<Vec3> nodePositionsFiber;
            nodePositionsFiber.clear();
            SlotConnectorDataHelper<SlotConnectorDataType>::
                slotGetGeometryValues(slotConnectorData, preciceData.slotNo,
                                      arrayIndex, dofNosLocalWithoutGhosts,
                                      nodePositionsFiber);
            geometryValues.insert(geometryValues.end(),
                                  nodePositionsFiber.begin(),
                                  nodePositionsFiber.end());
          }

          LOG(DEBUG) << "collected " << geometryValues.size()
                     << " node positions from the " << nArrayItems << " fibers";
        }

        // transform to contiguous memory layout for precice
        std::vector<double> geometryValuesContiguous(3 * geometryValues.size());
        for (int entryNo = 0; entryNo < geometryValues.size(); entryNo++)
          for (int componentNo = 0; componentNo < 3; componentNo++)
            geometryValuesContiguous[3 * entryNo + componentNo] =
                geometryValues[entryNo][componentNo];

        if (geometryValuesContiguous.size() != preciceMesh->nNodesLocal * 3)
          LOG(FATAL) << "size mismatch: " << geometryValuesContiguous.size()
                     << "!=" << preciceMesh->nNodesLocal * 3;

        // resize buffer for vertex ids
        preciceMesh->preciceVertexIds.resize(preciceMesh->nNodesLocal);

        // give the node positions to precice and get the vertex ids
        preciceParticipant->setMeshVertices(preciceMesh->preciceMeshName,
                                            geometryValuesContiguous,
                                            preciceMesh->preciceVertexIds);

        // store newly created mesh in preciceData
        preciceData.preciceMesh = preciceMesh;

        // store the precice mesh to the vector of meshes
        preciceVolumeMeshes.push_back(preciceMesh);

        // output message about mesh names
        std::stringstream message;
        message << "Initialized precice mesh \"" << preciceMesh->preciceMeshName
                << "\" from opendihu mesh \"" << preciceMesh->opendihuMeshName
                << "\" with " << preciceMesh->nNodesLocal << " local nodes";
        if (nArrayItems > 1)
          message << " (and other corresponding meshes, in total "
                  << nArrayItems << " fibers/compartments)";
        LOG(INFO) << message.str() << ".";
      } else {
        LOG(DEBUG) << "Use existing precice mesh " << currentMeshName;
        preciceData.preciceMesh = *iter;
      }

      // parse mode
      std::string mode = currentPreciceData.getOptionString("mode", "");
      if (mode == "read") {
        preciceData.ioType = PreciceAdapterInitialize<
            NestedSolverType>::PreciceVolumeData::ioRead;
      } else if (mode == "write") {
        preciceData.ioType = PreciceAdapterInitialize<
            NestedSolverType>::PreciceVolumeData::ioWrite;
      } else {
        LOG(FATAL) << currentPreciceData << "[\"mode\"] is \"" << mode << "\", "
                   << "possible values are: \"read\", \"write\".";
      }

      // parse variable name
      preciceData.preciceDataName =
          currentPreciceData.getOptionString("preciceDataName", "variable");

      // increment slotNo, this value is used for slots that are not identified
      // by slotName in the config
      if (!preciceData.isGeometryField)
        slotNo++;

      // check if mesh size matches the specified slot
      // get the mesh partition
      std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
          SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase(
              slotConnectorData, preciceData.slotNo, 0);

      int nDofsLocalWithoutGhosts =
          meshPartitionBase->nDofsLocalWithoutGhosts();
      int nArrayItems =
          SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
              slotConnectorData,
              preciceData.slotNo); // number of fibers if there are fibers

      if (nDofsLocalWithoutGhosts * nArrayItems !=
          preciceData.preciceMesh->nNodesLocal) {
        LOG(DEBUG) << ", all available slots: "
                   << SlotConnectorDataHelper<SlotConnectorDataType>::getString(
                          slotConnectorData);
        LOG(FATAL)
            << currentPreciceData
            << ": Mesh does not match slot in "
               "PreciceAdapterVolumeCoupling.\n\n "
            << "The " << (listIndex + 1)
            << (listIndex == 0
                    ? "st"
                    : (listIndex == 1 ? "nd" : (listIndex == 2 ? "rd" : "th")))
            << " list item under \"preciceData\" uses slot "
            << preciceData.slotNo << " (slotName \"" << preciceData.slotName
            << "\")"
            << " and preciceMesh \"" << preciceData.preciceMesh->preciceMeshName
            << "\".\nThe slot has " << nDofsLocalWithoutGhosts << " dofs * "
            << nArrayItems << " array items (fibers) = "
            << nDofsLocalWithoutGhosts * nArrayItems
            << " dofs in total. The mesh was initialized from opendihu mesh \""
            << preciceData.preciceMesh->opendihuMeshName << "\" with "
            << preciceData.preciceMesh->nNodesLocal << " nodes. \n"
            << "(isGeometryField: " << std::boolalpha
            << preciceData.isGeometryField << ")\n"
            << "You can do the following:\n"
            << "- Rename the preciceMesh to a unique name and update the "
               "precice "
               "config xml file accordingly.\n"
            << "- Check the solver structure file to find out the meshes that "
               "are associated with the slots.\n"
            << "  Check that the correct mesh is specified and that the "
               "precice "
               "mesh was not initialized earlier with a different opendihu "
               "mesh.\n"
            << "  Every precice mesh is initialized the first time in appears "
               "under \"preciceData\", sometimes reordering the entries can "
               "help.\n"
            << "  If the desired mesh is not available at the current solver, "
               "maybe insert a MapDofs class.";
      }

      LOG(INFO) << "Precice data \"" << preciceData.preciceDataName
                << "\" maps to "
                << (preciceData.isGeometryField ? "the geometry field of " : "")
                << "slot " << preciceData.slotNo << " (\""
                << preciceData.slotName << "\") and uses precice mesh \""
                << preciceData.preciceMesh->preciceMeshName
                << "\", which is opendihu mesh \""
                << preciceData.preciceMesh->opendihuMeshName << "\" with "
                << preciceData.preciceMesh->nNodesLocal << " local nodes.";

      // store preciceData to vector
      preciceVolumeData.push_back(preciceData);
    }
  };

  template <typename SurfaceDataVector, typename VolumeDataVector>
  void
  preciceReadData(NestedSolverType &nestedSolver,
                  std::shared_ptr<precice::Participant> &preciceParticipant,
                  SurfaceDataVector &preciceSurfaceData,
                  VolumeDataVector &preciceVolumeData, DihuContext &context) {

    LOG(DEBUG) << "read volume data from precice";
    double preciceDt = preciceParticipant->getMaxTimeStepSize();

    using SlotConnectorDataType =
        typename NestedSolverType::SlotConnectorDataType;
    std::shared_ptr<SlotConnectorDataType> slotConnectorData =
        nestedSolver.getSlotConnectorData();

    // loop over data
    for (auto &preciceData : preciceVolumeData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceVolumeData::ioRead) {
        int nEntries = preciceData.preciceMesh->nNodesLocal;

        if (preciceData.isGeometryField) {
          nEntries = preciceData.preciceMesh->nNodesLocal * 3;
        }

        // allocate temporary memory
        scalarValues_.resize(nEntries);

        // get all data at once
        preciceParticipant->readData(preciceData.preciceMesh->preciceMeshName,
                                     preciceData.preciceDataName,
                                     preciceData.preciceMesh->preciceVertexIds,
                                     preciceDt, scalarValues_);

        // get the mesh partition
        std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
            SlotConnectorDataHelper<
                SlotConnectorDataType>::getMeshPartitionBase(slotConnectorData,
                                                             preciceData.slotNo,
                                                             0);

        int nDofsLocalWithoutGhosts =
            meshPartitionBase->nDofsLocalWithoutGhosts();

        // get the vector of values [0,1,...,nDofsLocalWithGhosts]
        const std::vector<PetscInt> &dofNosLocalWithGhosts =
            meshPartitionBase->dofNosLocal();
        std::vector<PetscInt> dofNosLocalWithoutGhosts(
            dofNosLocalWithGhosts.begin(),
            dofNosLocalWithGhosts.begin() + nDofsLocalWithoutGhosts);

        int nArrayItems =
            SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
                slotConnectorData,
                preciceData.slotNo); // number of fibers if there are fibers

        // store received data in field variable
        if (preciceData.isGeometryField) {
          // loop over fibers if there are any
          for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
            // fill the vector geometryValues_ with the geometry values of the
            // current fiber or mesh
            geometryValues_.resize(nDofsLocalWithoutGhosts);
            for (int dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts;
                 dofNoLocal++) {
              for (int componentNo = 0; componentNo < 3; componentNo++) {
                geometryValues_[dofNoLocal][componentNo] =
                    scalarValues_[3 * (arrayIndex * nDofsLocalWithoutGhosts +
                                       dofNoLocal) +
                                  componentNo];
                //
                std::cout << geometryValues_[dofNoLocal][componentNo]
                          << std::endl;
              }
            }

            SlotConnectorDataHelper<SlotConnectorDataType>::
                slotSetGeometryValues(slotConnectorData, preciceData.slotNo,
                                      arrayIndex, dofNosLocalWithoutGhosts,
                                      geometryValues_);
          }
        } else {
          // loop over fibers if there are any
          for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
            // fill the vector geometryValues_ with the geometry values of the
            // current fiber or mesh
            scalarValuesOfMesh_.assign(
                scalarValues_.begin() + arrayIndex * nDofsLocalWithoutGhosts,
                scalarValues_.begin() +
                    (arrayIndex + 1) * nDofsLocalWithoutGhosts);

            SlotConnectorDataHelper<SlotConnectorDataType>::slotSetValues(
                slotConnectorData, preciceData.slotNo, arrayIndex,
                dofNosLocalWithoutGhosts, scalarValuesOfMesh_);
          }
        }
      }
    }

    LOG(DEBUG) << "read surface data from precice";
    // loop over data
    for (auto &preciceData : preciceSurfaceData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceData::ioRead) {
        // allocate memory
        int nEntries = preciceData.preciceMesh->nNodesLocal * 3;

        // if the data is displacements and velocities
        if (!preciceData.displacementsName.empty()) {
          displacementValues_.resize(nEntries);
          velocityValues_.resize(nEntries);

          // get all data at once
          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.displacementsName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              displacementValues_);

          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.velocitiesName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              velocityValues_);

          setDirichletBoundaryConditions(preciceData, nestedSolver);
        }
        // if the data is traction
        else if (!preciceData.tractionName.empty()) {
          tractionValues_.resize(nEntries);
          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              tractionValues_);

          setNeumannBoundaryConditions(preciceData, nestedSolver, context);
        } else {
          LOG(FATAL) << "Unknown precice data (read), none of displacements, "
                        "velocities or traction is set.";
        }
      }
    };
  }

  template <typename SurfaceDataVector, typename VolumeDataVector>
  void
  preciceWriteData(NestedSolverType &nestedSolver,
                   std::shared_ptr<precice::Participant> &preciceParticipant,
                   SurfaceDataVector &preciceSurfaceData,
                   VolumeDataVector &preciceVolumeData, double scalingFactor) {
    // write data to precice
    LOG(DEBUG) << "write volume data to precice";

    using SlotConnectorDataType =
        typename NestedSolverType::SlotConnectorDataType;
    std::shared_ptr<SlotConnectorDataType> slotConnectorData =
        nestedSolver.getSlotConnectorData();

    // loop over data
    for (auto &preciceData : preciceVolumeData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceVolumeData::ioWrite) {
        scalarValues_.clear();

        // get the mesh partition
        std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase =
            SlotConnectorDataHelper<
                SlotConnectorDataType>::getMeshPartitionBase(slotConnectorData,
                                                             preciceData.slotNo,
                                                             0);

        int nDofsLocalWithoutGhosts =
            meshPartitionBase->nDofsLocalWithoutGhosts();

        // get the vector of values [0,1,...,nDofsLocalWithGhosts]
        const std::vector<PetscInt> &dofNosLocalWithGhosts =
            meshPartitionBase->dofNosLocal();
        std::vector<PetscInt> dofNosLocalWithoutGhosts(
            dofNosLocalWithGhosts.begin(),
            dofNosLocalWithGhosts.begin() + nDofsLocalWithoutGhosts);

        int nArrayItems =
            SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(
                slotConnectorData,
                preciceData.slotNo); // number of fibers if there are fibers

        // if it is a geometry field, get the node positions of a mesh
        if (preciceData.isGeometryField) {
          geometryValues_.clear();

          // loop over fibers if there are any
          for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
            static std::vector<Vec3> geometryValuesFiber;
            geometryValuesFiber.clear();
            SlotConnectorDataHelper<SlotConnectorDataType>::
                slotGetGeometryValues(slotConnectorData, preciceData.slotNo,
                                      arrayIndex, dofNosLocalWithoutGhosts,
                                      geometryValuesFiber);

            geometryValues_.insert(geometryValues_.end(),
                                   geometryValuesFiber.begin(),
                                   geometryValuesFiber.end());
          }
          LOG(DEBUG) << "Using geometry field of opendihu meshes of "
                     << nArrayItems << " array items (e.g., fibers).";

          // transform to contiguous memory layout for precice
          scalarValues_.resize(3 * geometryValues_.size());

          for (int entryNo = 0; entryNo < geometryValues_.size(); entryNo++)
            for (int componentNo = 0; componentNo < 3; componentNo++)
              scalarValues_[3 * entryNo + componentNo] =
                  geometryValues_[entryNo][componentNo];

        } else {
          // the data is a normal slot, no geometry field

          // loop over fibers if there are any
          for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++) {
            static std::vector<double> values;
            values.clear();
            // static void slotGetValues(std::shared_ptr<SlotConnectorDataType>
            // slotConnectorData,
            //   int slotNo, int arrayIndex, const std::vector<dof_no_t>
            //   &dofNosLocal, std::vector<double> &values);
            SlotConnectorDataHelper<SlotConnectorDataType>::slotGetValues(
                slotConnectorData, preciceData.slotNo, arrayIndex,
                dofNosLocalWithoutGhosts, values);
            scalarValues_.insert(scalarValues_.end(), values.begin(),
                                 values.end());

            LOG(DEBUG) << "arrayIndex " << arrayIndex << ", add "
                       << values.size()
                       << " values, now number: " << scalarValues_.size();
          }
        }

        // scale the values by a factor given in the config
        for (double &value : scalarValues_)
          value *= scalingFactor;

        // check dim for scalarValues_
        if (preciceData.isGeometryField) {
          if (scalarValues_.size() != 3 * preciceData.preciceMesh->nNodesLocal)
            LOG(FATAL) << "Wrong number of scalar values (isGeometryField): "
                       << scalarValues_.size()
                       << " != " << preciceData.preciceMesh->nNodesLocal
                       << ", nArrayItems: " << nArrayItems;
        } else {
          if (scalarValues_.size() != preciceData.preciceMesh->nNodesLocal)
            LOG(FATAL) << "Wrong number of scalar values: "
                       << scalarValues_.size()
                       << " != " << preciceData.preciceMesh->nNodesLocal
                       << ", nArrayItems: " << nArrayItems;
        }

        // write values in precice
        preciceParticipant->writeData(preciceData.preciceMesh->preciceMeshName,
                                      preciceData.preciceDataName,
                                      preciceData.preciceMesh->preciceVertexIds,
                                      scalarValues_);
      }
    }
    LOG(DEBUG) << "write volume data to precice complete";

    LOG(DEBUG) << "write surface data to precice";

    // loop over data
    for (auto &preciceData : preciceSurfaceData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceData::ioWrite) {
        // if the data is displacements and velocities
        if (!preciceData.displacementsName.empty()) {
          // convert geometry values to precice data layout
          displacementValues_.clear();
          velocityValues_.clear();

          this->getDisplacementVelocityValues(
              nestedSolver, preciceData.preciceMesh->dofNosLocal,
              displacementValues_, velocityValues_);

#ifndef NDEBUG
          LOG(DEBUG) << "write displacements data to precice: "
                     << displacementValues_;
#endif
          // scale displacement and velocity values
          for (double &value : displacementValues_)
            value *= scalingFactor;

          for (double &value : velocityValues_)
            value *= scalingFactor;

          // write displacement values in precice
          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.displacementsName,
              preciceData.preciceMesh->preciceVertexIds, displacementValues_);

          // write velocity values in precice
          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.velocitiesName,
              preciceData.preciceMesh->preciceVertexIds, velocityValues_);
        }
        // if the data is traction
        else if (!preciceData.tractionName.empty() &&
                 preciceData.average == true) {
          LOG(INFO) << "Write averaged traction";
          // convert geometry values to precice data layout
          tractionValues_.clear();
          this->getTractionValues(nestedSolver,
                                  preciceData.preciceMesh->dofNosLocal,
                                  tractionValues_);
          // average z-values of traction
          double average_traction = 0.0;
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            average_traction += tractionValues_[i];
          }
          average_traction /= (tractionValues_.size() / 3);
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            tractionValues_[i] = average_traction;
          }
#ifndef NDEBUG
          LOG(DEBUG) << "Write averaged-traction data to precice: "
                     << tractionValues_[2];
#endif
          // scale traction values, they are always scaled by the factor of -1
          for (double &value : tractionValues_) {
            value *= scalingFactor;
            value *= -1;
          }

          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, tractionValues_);
        } else if (!preciceData.tractionName.empty()) {

          LOG(INFO) << "Write non-averaged traction";

          // convert geometry values to precice data layout
          tractionValues_.clear();
          this->getTractionValues(nestedSolver,
                                  preciceData.preciceMesh->dofNosLocal,
                                  tractionValues_);

#ifndef NDEBUG
          LOG(DEBUG) << "write traction data to precice: " << tractionValues_;
          std::stringstream s;
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            s << " " << tractionValues_[i];
          }
          LOG(DEBUG) << "z values of traction: " << s.str();
#endif
          // scale traction values, they are always scaled by the factor of -1
          for (double &value : tractionValues_) {
            value *= scalingFactor;
            value *= -1;
          }

          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, tractionValues_);
        } else {
          LOG(FATAL) << "Unknown precice data (write), none of displacements, "
                        "velocities or traction is set.";
        }
      }
    }

    LOG(DEBUG) << "write traction data to precice complete";
  }

  template <typename SurfaceData>
  void setDirichletBoundaryConditions(SurfaceData &preciceData,
                                      NestedSolverType &nestedSolver) {
    std::vector<std::pair<global_no_t, std::array<double, 6>>>
        newDirichletBCValues;

    auto functionSpace = this->functionSpace(nestedSolver);

    // if this rank has no data, do not set any boundary conditions
    if (preciceData.preciceMesh->nNodesLocal != 0) {

      // loop over selected nodes of surface mesh
      for (const int &surfaceDofNoLocal :
           preciceData.preciceMesh->selectedDofNosLocal) {
        LOG(DEBUG) << "In mesh: " << preciceData.preciceMesh->preciceMeshName
                   << ", surfaceDofNosLocal = " << surfaceDofNoLocal
                   << "and dofNosLocal ="
                   << preciceData.preciceMesh->dofNosLocal[surfaceDofNoLocal];

        dof_no_t dofNoLocal =
            preciceData.preciceMesh->dofNosLocal[surfaceDofNoLocal];
        global_no_t dofNoGlobal =
            functionSpace->meshPartition()->getDofNoGlobalPetsc(dofNoLocal);
        // assign received values to dirichlet bc vector of size 6
        std::array<double, 6> newDirichletBCValue;

        for (int i = 0; i < 3; i++) {
          newDirichletBCValue[i] =
              displacementValues_[3 * surfaceDofNoLocal + i];
          newDirichletBCValue[3 + i] =
              velocityValues_[3 * surfaceDofNoLocal + i];
        }

        newDirichletBCValues.push_back(
            std::pair<global_no_t, std::array<double, 6>>(dofNoGlobal,
                                                          newDirichletBCValue));
      }

      LOG(DEBUG) << "read data from precice complete, displacement values: "
                 << displacementValues_
                 << ", velocityValues: " << velocityValues_;
      LOG(DEBUG) << "read and set Dirichlet BC: " << newDirichletBCValues;
    }

    //! set new dirichlet boundary condition values
    this->updateDirichletBoundaryConditions(nestedSolver, newDirichletBCValues);
  }

  template <typename SurfaceData>
  void setNeumannBoundaryConditions(SurfaceData &preciceData,
                                    NestedSolverType &nestedSolver,
                                    DihuContext &context) {
    // set traction values as neumann boundary conditions
    using FunctionSpace =
        typename PreciceAdapterNestedSolver<NestedSolverType>::FunctionSpace;
    using ElementWithFacesType =
        typename SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>::ElementWithFaces;

    std::vector<ElementWithFacesType> neumannBoundaryConditionElements;

    auto functionSpace = this->functionSpace(nestedSolver);

    // if this rank has no data, do not set any boundary conditions
    if (preciceData.preciceMesh->nNodesLocal != 0) {

      const int nNodesX = functionSpace->nNodesLocalWithoutGhosts(0);
      const int nNodesY = functionSpace->nNodesLocalWithoutGhosts(1);
      const int nNodesZ = functionSpace->nNodesLocalWithoutGhosts(2);

      int nElementsX = functionSpace->meshPartition()->nElementsLocal(0);
      int nElementsY = functionSpace->meshPartition()->nElementsLocal(1);
      int nElementsZ = functionSpace->meshPartition()->nElementsLocal(2);

      // set node and element indics in z direction for bottom surface
      int nodeIndexZ = 0;
      int elementalNodeIndexZ = 0;
      int elementIndexZ = 0;

      // for top surface
      if (preciceData.preciceMesh->face ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
        nodeIndexZ = nNodesZ - 1;
        elementIndexZ = nElementsZ - 1;
        elementalNodeIndexZ = 2;
      }

      // loop over elements
      for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++) {
        for (int elementIndexX = 0; elementIndexX < nElementsX;
             elementIndexX++) {
          ElementWithFacesType elementWithFaces;
          element_no_t elementNoLocal =
              elementIndexZ * nElementsX * nElementsY +
              elementIndexY * nElementsX + elementIndexX;
          elementWithFaces.elementNoLocal = elementNoLocal;

          // set surface dofs
          Mesh::face_t face = Mesh::face_t::face2Minus;
          if (preciceData.preciceMesh->face ==
              PreciceAdapterInitialize<
                  NestedSolverType>::PreciceSurfaceMesh::face2Plus)
            face = Mesh::face_t::face2Plus;

          elementWithFaces.face = face;

          // get dofs indices within the numbering of the volume element that
          // correspond to the selected face
          const int nDofsPerNode = FunctionSpace::nDofsPerNode();
          const int nSurfaceDofs =
              ::FunctionSpace::FunctionSpaceBaseDim<
                  2,
                  typename FunctionSpace::BasisFunction>::nNodesPerElement() *
              nDofsPerNode;
          std::array<dof_no_t, nSurfaceDofs> surfaceDofs;
          FunctionSpace::getFaceDofs(face, surfaceDofs);

          elementWithFaces.surfaceDofs.assign(surfaceDofs.begin(),
                                              surfaceDofs.end());

          // loop over the nodes of the element
          for (int elementalNodeIndexY = 0; elementalNodeIndexY < 3;
               elementalNodeIndexY++) {
            for (int elementalNodeIndexX = 0; elementalNodeIndexX < 3;
                 elementalNodeIndexX++) {
              int elementalDofIndex = elementalNodeIndexZ * 9 +
                                      elementalNodeIndexY * 3 +
                                      elementalNodeIndexX;

              dof_no_t dofNoLocal =
                  functionSpace->getDofNo(elementNoLocal, elementalDofIndex);

              // only iterate over local dofs, the ghost dofs are not considered
              // here (although it would be correct to communicate them with
              // precice)
              if (dofNoLocal >= functionSpace->nDofsLocalWithoutGhosts())
                continue;

              int valueIndex = dofNoLocal - nodeIndexZ * nNodesX * nNodesY;

              LOG(DEBUG) << "(x,y,z)=(" << elementalNodeIndexX << ","
                         << elementalNodeIndexY << "," << elementalNodeIndexZ
                         << ") dofNoLocal " << dofNoLocal
                         << ", valueIndex: " << valueIndex << "/"
                         << tractionValues_.size() / 3;

              assert(valueIndex >= 0);
              if (valueIndex >= preciceData.preciceMesh->nNodesLocal) {
                LOG(ERROR) << "valueIndex: " << valueIndex
                           << ", dofNoLocal: " << dofNoLocal
                           << ", nNodes: " << nNodesX << "," << nNodesY
                           << ", nodeIndexZ: " << nodeIndexZ
                           << ", nNodesLocal: "
                           << preciceData.preciceMesh->nNodesLocal
                           << " elementNoLocal: " << elementNoLocal
                           << ", elementalDofIndex: " << elementalDofIndex
                           << ", elementalNodeIndexZ: " << elementalNodeIndexZ
                           << ", meshPartition: "
                           << *functionSpace->meshPartition();
              }
              assert(valueIndex < preciceData.preciceMesh->nNodesLocal);

              Vec3 traction;
              for (int i = 0; i < 3; i++) {
                traction[i] = tractionValues_[3 * valueIndex + i];
              }

              dof_no_t surfaceDof = elementalDofIndex;

              // for top surface
              if (preciceData.preciceMesh->face ==
                  PreciceAdapterInitialize<
                      NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
                surfaceDof = 18 + elementalDofIndex;
              }
              elementWithFaces.dofVectors.push_back(
                  std::pair<dof_no_t, Vec3>(surfaceDof, traction));

              // LOG(INFO) << "dofVectors: " << elementWithFaces.dofVectors <<
              // ", traction: " << traction;
            }
          }
          neumannBoundaryConditionElements.push_back(elementWithFaces);
        }
      }
    }
    // create new Neumann BC object
    using NeumannBoundaryConditionsType =
        SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>;
    std::shared_ptr<NeumannBoundaryConditionsType> neumannBoundaryConditions =
        std::make_shared<NeumannBoundaryConditionsType>(context);
    neumannBoundaryConditions->initialize(functionSpace,
                                          neumannBoundaryConditionElements);
    neumannBoundaryConditions->setDeformationGradientField(
        this->deformationGradientField(nestedSolver));

#ifndef NDEBUG
    std::stringstream s;
    for (int i = 0; i < neumannBoundaryConditionElements.size(); i++) {
      s << "{el. " << neumannBoundaryConditionElements[i].elementNoLocal
        << ", \"" << Mesh::getString(neumannBoundaryConditionElements[i].face)
        << "\", dofVectors: [";
      for (int j = 0; j < neumannBoundaryConditionElements[i].dofVectors.size();
           j++) {
        s << "(" << neumannBoundaryConditionElements[i].dofVectors[j].first
          << ": (";
        for (int k = 0;
             k <
             neumannBoundaryConditionElements[i].dofVectors[j].second.size();
             k++) {
          if (k != 0)
            s << ",";
          s << neumannBoundaryConditionElements[i].dofVectors[j].second[k];
        }
        s << ")),";
      }
      s << "], surfaceDofs: ";
      for (int j = 0;
           j < neumannBoundaryConditionElements[i].surfaceDofs.size(); j++)
        s << neumannBoundaryConditionElements[i].surfaceDofs[j] << ",";
      s << "} ";
    }
    LOG(DEBUG) << "read and set Neumann BC:\n" << s.str();
#endif

    // set Neumann BCs in the static hyperelasticity of the
    // TimeSteppingScheme::DynamicHyperelasticitySolver solver
    this->updateNeumannBoundaryConditions(nestedSolver,
                                          neumannBoundaryConditions);

    LOG(DEBUG) << "read data from precice complete, traction values: "
               << tractionValues_;
  }
  //! initialize dirichlet boundary conditions by adding new dofs and prescribed
  //! values for all bottom or top nodes
  void addDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<
          typename SpatialDiscretization::DirichletBoundaryConditionsBase<
              FunctionSpace, 6>::ElementWithNodes>
          &dirichletBoundaryConditionElements);

  //! update existing boundary conditions with new values
  void updateDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<std::pair<global_no_t, std::array<double, 6>>>
          newDirichletBoundaryConditionValues);

  //! update the neumann boundary conditions by replacing the complete object
  void updateNeumannBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>>
          neumannBoundaryConditions);

  //! get the displacement and velocity vectors of the given local dof nos
  void getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                     const std::vector<dof_no_t> &dofNosLocal,
                                     std::vector<double> &displacementValues,
                                     std::vector<double> &velocityValues);

  //! get the traction vectors of the given local dof nos
  void getTractionValues(NestedSolverType &nestedSolver,
                         const std::vector<dof_no_t> &dofNosLocal,
                         std::vector<double> &tractionValues);

  //! get at Petsc Vec that stores all values of the current state, to be used
  //! to store and restore checkpoints
  Vec currentState(NestedSolverType &nestedSolver);

  //! get the field variable of the deformation gradient
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 9>>
  deformationGradientField(NestedSolverType &nestedSolver);

  void reset(NestedSolverType &nestedSolver);

  void saveFiberData(NestedSolverType &nestedSolver);

  void loadFiberData(NestedSolverType &nestedSolver);
};

/** Partial specialization for tendon or pure mechanics solver in a coupling
 * scheme, muscle contraction solver (nonlinear elasticity with active stress)
 */
template <typename T1, typename T2, typename T3>
class PreciceAdapterNestedSolver<
    Control::Coupling<T1, MuscleContractionSolver<T2, T3>>> {
public:
  std::vector<double> displacementValues_;
  std::vector<double> velocityValues_;
  std::vector<double> tractionValues_;
  std::vector<Vec3> displacementVectors_;
  std::vector<Vec3> velocityVectors_;
  std::vector<Vec3> tractionVectors_;
  std::vector<double> scalarValues_;
  std::vector<double> scalarValuesOfMesh_;
  std::vector<Vec3> geometryValues_;
  //! define the type of the nested solver
  typedef Control::Coupling<T1, MuscleContractionSolver<T2, T3>>
      NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef
      typename NestedSolverType::TimeStepping2Type::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<FunctionSpace> functionSpace(NestedSolverType &nestedSolver);

  //! parse the options in "preciceVolumeData" and initialize all variables in
  //! precice, store in variable preciceVolumeData_
  template <typename PreciceVolumeData, typename PreciceVolumeMesh>
  void initializePreciceVolumeData(
      PythonConfig &specificSettings, NestedSolverType &nestedSolver,
      std::shared_ptr<precice::Participant> &preciceParticipant,
      std::vector<PreciceVolumeData> &preciceVolumeData,
      std::vector<std::shared_ptr<PreciceVolumeMesh>> &preciceVolumeMeshes){};

  template <typename SurfaceDataVector, typename VolumeDataVector>
  void
  preciceReadData(NestedSolverType &nestedSolver,
                  std::shared_ptr<precice::Participant> &preciceParticipant,
                  SurfaceDataVector &preciceSurfaceData,
                  VolumeDataVector &preciceVolumeData, DihuContext &context) {
    LOG(DEBUG) << "read surface data from precice";
    double preciceDt = preciceParticipant->getMaxTimeStepSize();
    // loop over data
    for (auto &preciceData : preciceSurfaceData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceData::ioRead) {
        // allocate memory
        int nEntries = preciceData.preciceMesh->nNodesLocal * 3;

        // if the data is displacements and velocities
        if (!preciceData.displacementsName.empty()) {
          displacementValues_.resize(nEntries);
          velocityValues_.resize(nEntries);

          // get all data at once
          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.displacementsName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              displacementValues_);

          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.velocitiesName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              velocityValues_);

          setDirichletBoundaryConditions(preciceData, nestedSolver);
        }
        // if the data is traction
        else if (!preciceData.tractionName.empty()) {
          tractionValues_.resize(nEntries);
          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              tractionValues_);

          setNeumannBoundaryConditions(preciceData, nestedSolver, context);
        } else {
          LOG(FATAL) << "Unknown precice data (read), none of displacements, "
                        "velocities or traction is set.";
        }
      }
    };
  }

  template <typename SurfaceDataVector, typename VolumeDataVector>
  void
  preciceWriteData(NestedSolverType &nestedSolver,
                   std::shared_ptr<precice::Participant> &preciceParticipant,
                   SurfaceDataVector &preciceSurfaceData,
                   VolumeDataVector &preciceVolumeData, double scalingFactor) {
    // write data to precice
    LOG(DEBUG) << "write surface data to precice";

    // loop over data
    for (auto &preciceData : preciceSurfaceData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceData::ioWrite) {
        // if the data is displacements and velocities
        if (!preciceData.displacementsName.empty()) {
          // convert geometry values to precice data layout
          displacementValues_.clear();
          velocityValues_.clear();

          this->getDisplacementVelocityValues(
              nestedSolver, preciceData.preciceMesh->dofNosLocal,
              displacementValues_, velocityValues_);

#ifndef NDEBUG
          LOG(DEBUG) << "write displacements data to precice: "
                     << displacementValues_;
#endif
          // scale displacement and velocity values
          for (double &value : displacementValues_)
            value *= scalingFactor;

          for (double &value : velocityValues_)
            value *= scalingFactor;

          // write displacement values in precice
          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.displacementsName,
              preciceData.preciceMesh->preciceVertexIds, displacementValues_);

          // write velocity values in precice
          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.velocitiesName,
              preciceData.preciceMesh->preciceVertexIds, velocityValues_);
        }
        // if the data is traction
        else if (!preciceData.tractionName.empty() &&
                 preciceData.average == true) {
          LOG(INFO) << "Write averaged traction";
          // convert geometry values to precice data layout
          tractionValues_.clear();
          this->getTractionValues(nestedSolver,
                                  preciceData.preciceMesh->dofNosLocal,
                                  tractionValues_);
          // average z-values of traction
          double average_traction = 0.0;
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            average_traction += tractionValues_[i];
          }
          average_traction /= (tractionValues_.size() / 3);
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            tractionValues_[i] = average_traction;
          }
#ifndef NDEBUG
          LOG(DEBUG) << "Write averaged-traction data to precice: "
                     << tractionValues_[2];
#endif
          // scale traction values, they are always scaled by the factor of -1
          for (double &value : tractionValues_) {
            value *= scalingFactor;
            value *= -1;
          }

          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, tractionValues_);
        } else if (!preciceData.tractionName.empty()) {

          LOG(INFO) << "Write non-averaged traction";

          // convert geometry values to precice data layout
          tractionValues_.clear();
          this->getTractionValues(nestedSolver,
                                  preciceData.preciceMesh->dofNosLocal,
                                  tractionValues_);

#ifndef NDEBUG
          LOG(DEBUG) << "write traction data to precice: " << tractionValues_;
          std::stringstream s;
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            s << " " << tractionValues_[i];
          }
          LOG(DEBUG) << "z values of traction: " << s.str();
#endif
          // scale traction values, they are always scaled by the factor of -1
          for (double &value : tractionValues_) {
            value *= scalingFactor;
            value *= -1;
          }

          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, tractionValues_);
        } else {
          LOG(FATAL) << "Unknown precice data (write), none of displacements, "
                        "velocities or traction is set.";
        }
      }
    }

    LOG(DEBUG) << "write traction data to precice complete";
  }

  template <typename SurfaceData>
  void setDirichletBoundaryConditions(SurfaceData &preciceData,
                                      NestedSolverType &nestedSolver) {
    std::vector<std::pair<global_no_t, std::array<double, 6>>>
        newDirichletBCValues;

    auto functionSpace = this->functionSpace(nestedSolver);

    // if this rank has no data, do not set any boundary conditions
    if (preciceData.preciceMesh->nNodesLocal != 0) {

      // loop over selected nodes of surface mesh
      for (const int &surfaceDofNoLocal :
           preciceData.preciceMesh->selectedDofNosLocal) {
        LOG(DEBUG) << "In mesh: " << preciceData.preciceMesh->preciceMeshName
                   << ", surfaceDofNosLocal = " << surfaceDofNoLocal
                   << "and dofNosLocal ="
                   << preciceData.preciceMesh->dofNosLocal[surfaceDofNoLocal];

        dof_no_t dofNoLocal =
            preciceData.preciceMesh->dofNosLocal[surfaceDofNoLocal];
        global_no_t dofNoGlobal =
            functionSpace->meshPartition()->getDofNoGlobalPetsc(dofNoLocal);
        // assign received values to dirichlet bc vector of size 6
        std::array<double, 6> newDirichletBCValue;

        for (int i = 0; i < 3; i++) {
          newDirichletBCValue[i] =
              displacementValues_[3 * surfaceDofNoLocal + i];
          newDirichletBCValue[3 + i] =
              velocityValues_[3 * surfaceDofNoLocal + i];
        }

        newDirichletBCValues.push_back(
            std::pair<global_no_t, std::array<double, 6>>(dofNoGlobal,
                                                          newDirichletBCValue));
      }

      LOG(DEBUG) << "read data from precice complete, displacement values: "
                 << displacementValues_
                 << ", velocityValues: " << velocityValues_;
      LOG(DEBUG) << "read and set Dirichlet BC: " << newDirichletBCValues;
    }

    //! set new dirichlet boundary condition values
    this->updateDirichletBoundaryConditions(nestedSolver, newDirichletBCValues);
  }

  template <typename SurfaceData>
  void setNeumannBoundaryConditions(SurfaceData &preciceData,
                                    NestedSolverType &nestedSolver,
                                    DihuContext &context) {
    // set traction values as neumann boundary conditions
    using FunctionSpace =
        typename PreciceAdapterNestedSolver<NestedSolverType>::FunctionSpace;
    using ElementWithFacesType =
        typename SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>::ElementWithFaces;

    std::vector<ElementWithFacesType> neumannBoundaryConditionElements;

    auto functionSpace = this->functionSpace(nestedSolver);

    // if this rank has no data, do not set any boundary conditions
    if (preciceData.preciceMesh->nNodesLocal != 0) {

      const int nNodesX = functionSpace->nNodesLocalWithoutGhosts(0);
      const int nNodesY = functionSpace->nNodesLocalWithoutGhosts(1);
      const int nNodesZ = functionSpace->nNodesLocalWithoutGhosts(2);

      int nElementsX = functionSpace->meshPartition()->nElementsLocal(0);
      int nElementsY = functionSpace->meshPartition()->nElementsLocal(1);
      int nElementsZ = functionSpace->meshPartition()->nElementsLocal(2);

      // set node and element indics in z direction for bottom surface
      int nodeIndexZ = 0;
      int elementalNodeIndexZ = 0;
      int elementIndexZ = 0;

      // for top surface
      if (preciceData.preciceMesh->face ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
        nodeIndexZ = nNodesZ - 1;
        elementIndexZ = nElementsZ - 1;
        elementalNodeIndexZ = 2;
      }

      // loop over elements
      for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++) {
        for (int elementIndexX = 0; elementIndexX < nElementsX;
             elementIndexX++) {
          ElementWithFacesType elementWithFaces;
          element_no_t elementNoLocal =
              elementIndexZ * nElementsX * nElementsY +
              elementIndexY * nElementsX + elementIndexX;
          elementWithFaces.elementNoLocal = elementNoLocal;

          // set surface dofs
          Mesh::face_t face = Mesh::face_t::face2Minus;
          if (preciceData.preciceMesh->face ==
              PreciceAdapterInitialize<
                  NestedSolverType>::PreciceSurfaceMesh::face2Plus)
            face = Mesh::face_t::face2Plus;

          elementWithFaces.face = face;

          // get dofs indices within the numbering of the volume element that
          // correspond to the selected face
          const int nDofsPerNode = FunctionSpace::nDofsPerNode();
          const int nSurfaceDofs =
              ::FunctionSpace::FunctionSpaceBaseDim<
                  2,
                  typename FunctionSpace::BasisFunction>::nNodesPerElement() *
              nDofsPerNode;
          std::array<dof_no_t, nSurfaceDofs> surfaceDofs;
          FunctionSpace::getFaceDofs(face, surfaceDofs);

          elementWithFaces.surfaceDofs.assign(surfaceDofs.begin(),
                                              surfaceDofs.end());

          // loop over the nodes of the element
          for (int elementalNodeIndexY = 0; elementalNodeIndexY < 3;
               elementalNodeIndexY++) {
            for (int elementalNodeIndexX = 0; elementalNodeIndexX < 3;
                 elementalNodeIndexX++) {
              int elementalDofIndex = elementalNodeIndexZ * 9 +
                                      elementalNodeIndexY * 3 +
                                      elementalNodeIndexX;

              dof_no_t dofNoLocal =
                  functionSpace->getDofNo(elementNoLocal, elementalDofIndex);

              // only iterate over local dofs, the ghost dofs are not considered
              // here (although it would be correct to communicate them with
              // precice)
              if (dofNoLocal >= functionSpace->nDofsLocalWithoutGhosts())
                continue;

              int valueIndex = dofNoLocal - nodeIndexZ * nNodesX * nNodesY;

              LOG(DEBUG) << "(x,y,z)=(" << elementalNodeIndexX << ","
                         << elementalNodeIndexY << "," << elementalNodeIndexZ
                         << ") dofNoLocal " << dofNoLocal
                         << ", valueIndex: " << valueIndex << "/"
                         << tractionValues_.size() / 3;

              assert(valueIndex >= 0);
              if (valueIndex >= preciceData.preciceMesh->nNodesLocal) {
                LOG(ERROR) << "valueIndex: " << valueIndex
                           << ", dofNoLocal: " << dofNoLocal
                           << ", nNodes: " << nNodesX << "," << nNodesY
                           << ", nodeIndexZ: " << nodeIndexZ
                           << ", nNodesLocal: "
                           << preciceData.preciceMesh->nNodesLocal
                           << " elementNoLocal: " << elementNoLocal
                           << ", elementalDofIndex: " << elementalDofIndex
                           << ", elementalNodeIndexZ: " << elementalNodeIndexZ
                           << ", meshPartition: "
                           << *functionSpace->meshPartition();
              }
              assert(valueIndex < preciceData.preciceMesh->nNodesLocal);

              Vec3 traction;
              for (int i = 0; i < 3; i++) {
                traction[i] = tractionValues_[3 * valueIndex + i];
              }

              dof_no_t surfaceDof = elementalDofIndex;

              // for top surface
              if (preciceData.preciceMesh->face ==
                  PreciceAdapterInitialize<
                      NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
                surfaceDof = 18 + elementalDofIndex;
              }
              elementWithFaces.dofVectors.push_back(
                  std::pair<dof_no_t, Vec3>(surfaceDof, traction));

              // LOG(INFO) << "dofVectors: " << elementWithFaces.dofVectors <<
              // ", traction: " << traction;
            }
          }
          neumannBoundaryConditionElements.push_back(elementWithFaces);
        }
      }
    }
    // create new Neumann BC object
    using NeumannBoundaryConditionsType =
        SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>;
    std::shared_ptr<NeumannBoundaryConditionsType> neumannBoundaryConditions =
        std::make_shared<NeumannBoundaryConditionsType>(context);
    neumannBoundaryConditions->initialize(functionSpace,
                                          neumannBoundaryConditionElements);
    neumannBoundaryConditions->setDeformationGradientField(
        this->deformationGradientField(nestedSolver));

#ifndef NDEBUG
    std::stringstream s;
    for (int i = 0; i < neumannBoundaryConditionElements.size(); i++) {
      s << "{el. " << neumannBoundaryConditionElements[i].elementNoLocal
        << ", \"" << Mesh::getString(neumannBoundaryConditionElements[i].face)
        << "\", dofVectors: [";
      for (int j = 0; j < neumannBoundaryConditionElements[i].dofVectors.size();
           j++) {
        s << "(" << neumannBoundaryConditionElements[i].dofVectors[j].first
          << ": (";
        for (int k = 0;
             k <
             neumannBoundaryConditionElements[i].dofVectors[j].second.size();
             k++) {
          if (k != 0)
            s << ",";
          s << neumannBoundaryConditionElements[i].dofVectors[j].second[k];
        }
        s << ")),";
      }
      s << "], surfaceDofs: ";
      for (int j = 0;
           j < neumannBoundaryConditionElements[i].surfaceDofs.size(); j++)
        s << neumannBoundaryConditionElements[i].surfaceDofs[j] << ",";
      s << "} ";
    }
    LOG(DEBUG) << "read and set Neumann BC:\n" << s.str();
#endif

    // set Neumann BCs in the static hyperelasticity of the
    // TimeSteppingScheme::DynamicHyperelasticitySolver solver
    this->updateNeumannBoundaryConditions(nestedSolver,
                                          neumannBoundaryConditions);

    LOG(DEBUG) << "read data from precice complete, traction values: "
               << tractionValues_;
  }
  //! initialize dirichlet boundary conditions by adding new dofs and prescribed
  //! values for all bottom or top nodes
  void addDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<
          typename SpatialDiscretization::DirichletBoundaryConditionsBase<
              FunctionSpace, 6>::ElementWithNodes>
          &dirichletBoundaryConditionElements);

  //! update existing boundary conditions with new values
  void updateDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<std::pair<global_no_t, std::array<double, 6>>>
          newDirichletBoundaryConditionValues);

  //! update the neumann boundary conditions by replacing the complete object
  void updateNeumannBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>>
          neumannBoundaryConditions);

  //! get the displacement and velocity vectors of the given local dof nos
  void getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                     const std::vector<dof_no_t> &dofNosLocal,
                                     std::vector<double> &displacementValues,
                                     std::vector<double> &velocityValues);

  //! get the traction vectors of the given local dof nos
  void getTractionValues(NestedSolverType &nestedSolver,
                         const std::vector<dof_no_t> &dofNosLocal,
                         std::vector<double> &tractionValues);

  //! get at Petsc Vec that stores all values of the current state, to be used
  //! to store and restore checkpoints
  Vec currentState(NestedSolverType &nestedSolver);

  //! get the field variable of the deformation gradient
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 9>>
  deformationGradientField(NestedSolverType &nestedSolver);

  void reset(NestedSolverType &nestedSolver);

  void saveFiberData(NestedSolverType &nestedSolver);

  void loadFiberData(NestedSolverType &nestedSolver);
};

/** Partial specialization for tendon or pure mechanics solver, dynamic
 * nonlinear elasticity
 */
template <typename Material>
class PreciceAdapterNestedSolver<
    TimeSteppingScheme::DynamicHyperelasticitySolver<Material>> {
public:
  std::vector<double> displacementValues_;
  std::vector<double> velocityValues_;
  std::vector<double> tractionValues_;
  std::vector<Vec3> displacementVectors_;
  std::vector<Vec3> velocityVectors_;
  std::vector<Vec3> tractionVectors_;
  std::vector<double> scalarValues_;
  std::vector<double> scalarValuesOfMesh_;
  std::vector<Vec3> geometryValues_;
  //! define the type of the nested solver
  typedef TimeSteppingScheme::DynamicHyperelasticitySolver<Material>
      NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolverType::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<typename TimeSteppingScheme::DynamicHyperelasticitySolver<
      Material>::FunctionSpace>
  functionSpace(NestedSolverType &nestedSolver);

  //! parse the options in "preciceVolumeData" and initialize all variables in
  //! precice, store in variable preciceVolumeData_
  template <typename PreciceVolumeData, typename PreciceVolumeMesh>
  void initializePreciceVolumeData(
      PythonConfig &specificSettings, NestedSolverType &nestedSolver,
      std::shared_ptr<precice::Participant> &preciceParticipant,
      std::vector<PreciceVolumeData> &preciceVolumeData,
      std::vector<std::shared_ptr<PreciceVolumeMesh>> &preciceVolumeMeshes){};

  template <typename SurfaceDataVector, typename VolumeDataVector>
  void
  preciceReadData(NestedSolverType &nestedSolver,
                  std::shared_ptr<precice::Participant> &preciceParticipant,
                  SurfaceDataVector &preciceSurfaceData,
                  VolumeDataVector &preciceVolumeData, DihuContext &context) {
    LOG(DEBUG) << "read surface data from precice";
    double preciceDt = preciceParticipant->getMaxTimeStepSize();
    // loop over data
    for (auto &preciceData : preciceSurfaceData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceData::ioRead) {
        // allocate memory
        int nEntries = preciceData.preciceMesh->nNodesLocal * 3;

        // if the data is displacements and velocities
        if (!preciceData.displacementsName.empty()) {
          displacementValues_.resize(nEntries);
          velocityValues_.resize(nEntries);

          // get all data at once
          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.displacementsName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              displacementValues_);

          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.velocitiesName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              velocityValues_);

          setDirichletBoundaryConditions(preciceData, nestedSolver);
        }
        // if the data is traction
        else if (!preciceData.tractionName.empty()) {
          tractionValues_.resize(nEntries);
          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              tractionValues_);

          setNeumannBoundaryConditions(preciceData, nestedSolver, context);
        } else {
          LOG(FATAL) << "Unknown precice data (read), none of displacements, "
                        "velocities or traction is set.";
        }
      }
    };
  }

  template <typename SurfaceDataVector, typename VolumeDataVector>
  void
  preciceWriteData(NestedSolverType &nestedSolver,
                   std::shared_ptr<precice::Participant> &preciceParticipant,
                   SurfaceDataVector &preciceSurfaceData,
                   VolumeDataVector &preciceVolumeData, double scalingFactor) {
    // write data to precice
    LOG(DEBUG) << "write surface data to precice";

    // loop over data
    for (auto &preciceData : preciceSurfaceData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceData::ioWrite) {
        // if the data is displacements and velocities
        if (!preciceData.displacementsName.empty()) {
          // convert geometry values to precice data layout
          displacementValues_.clear();
          velocityValues_.clear();

          this->getDisplacementVelocityValues(
              nestedSolver, preciceData.preciceMesh->dofNosLocal,
              displacementValues_, velocityValues_);

#ifndef NDEBUG
          LOG(DEBUG) << "write displacements data to precice: "
                     << displacementValues_;
#endif
          // scale displacement and velocity values
          for (double &value : displacementValues_)
            value *= scalingFactor;

          for (double &value : velocityValues_)
            value *= scalingFactor;

          // write displacement values in precice
          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.displacementsName,
              preciceData.preciceMesh->preciceVertexIds, displacementValues_);

          // write velocity values in precice
          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.velocitiesName,
              preciceData.preciceMesh->preciceVertexIds, velocityValues_);
        }
        // if the data is traction
        else if (!preciceData.tractionName.empty() &&
                 preciceData.average == true) {
          LOG(INFO) << "Write averaged traction";
          // convert geometry values to precice data layout
          tractionValues_.clear();
          this->getTractionValues(nestedSolver,
                                  preciceData.preciceMesh->dofNosLocal,
                                  tractionValues_);
          // average z-values of traction
          double average_traction = 0.0;
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            average_traction += tractionValues_[i];
          }
          average_traction /= (tractionValues_.size() / 3);
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            tractionValues_[i] = average_traction;
          }
#ifndef NDEBUG
          LOG(DEBUG) << "Write averaged-traction data to precice: "
                     << tractionValues_[2];
#endif
          // scale traction values, they are always scaled by the factor of -1
          for (double &value : tractionValues_) {
            value *= scalingFactor;
            value *= -1;
          }

          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, tractionValues_);
        } else if (!preciceData.tractionName.empty()) {

          LOG(INFO) << "Write non-averaged traction";

          // convert geometry values to precice data layout
          tractionValues_.clear();
          this->getTractionValues(nestedSolver,
                                  preciceData.preciceMesh->dofNosLocal,
                                  tractionValues_);

#ifndef NDEBUG
          LOG(DEBUG) << "write traction data to precice: " << tractionValues_;
          std::stringstream s;
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            s << " " << tractionValues_[i];
          }
          LOG(DEBUG) << "z values of traction: " << s.str();
#endif
          // scale traction values, they are always scaled by the factor of -1
          for (double &value : tractionValues_) {
            value *= scalingFactor;
            value *= -1;
          }

          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, tractionValues_);
        } else {
          LOG(FATAL) << "Unknown precice data (write), none of displacements, "
                        "velocities or traction is set.";
        }
      }
    }

    LOG(DEBUG) << "write traction data to precice complete";
  }

  template <typename SurfaceData>
  void setDirichletBoundaryConditions(SurfaceData &preciceData,
                                      NestedSolverType &nestedSolver) {
    std::vector<std::pair<global_no_t, std::array<double, 6>>>
        newDirichletBCValues;

    auto functionSpace = this->functionSpace(nestedSolver);

    // if this rank has no data, do not set any boundary conditions
    if (preciceData.preciceMesh->nNodesLocal != 0) {
      // loop over nodes
      const int nNodesX = functionSpace->nNodesLocalWithoutGhosts(0);
      const int nNodesY = functionSpace->nNodesLocalWithoutGhosts(1);
      const int nNodesZ = functionSpace->nNodesLocalWithoutGhosts(2);

      // set node index in z direction for bottom surface
      int nodeIndexZ = 0;

      // for top surface
      if (preciceData.preciceMesh->face ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
        nodeIndexZ = nNodesZ - 1;
      }

      // loop over nodes to set the received values
      newDirichletBCValues.reserve(nNodesX * nNodesY);
      int valueIndex = 0;

      // loop over nodes of surface mesh
      for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++) {
        for (int nodeIndexX = 0; nodeIndexX < nNodesX;
             nodeIndexX++, valueIndex++) {
          node_no_t nodeNoLocal = nodeIndexZ * nNodesX * nNodesY +
                                  nodeIndexY * nNodesX + nodeIndexX;

          dof_no_t dofNoLocal = nodeNoLocal;
          global_no_t dofNoGlobal =
              functionSpace->meshPartition()->getDofNoGlobalPetsc(dofNoLocal);

          // assign received values to dirichlet bc vector of size 6
          std::array<double, 6> newDirichletBCValue;

          for (int i = 0; i < 3; i++) {
            newDirichletBCValue[i] = displacementValues_[3 * valueIndex + i];
            newDirichletBCValue[3 + i] = velocityValues_[3 * valueIndex + i];
          }

          newDirichletBCValues.push_back(
              std::pair<global_no_t, std::array<double, 6>>(
                  dofNoGlobal, newDirichletBCValue));
        }
      }

      LOG(DEBUG) << "read data from precice complete, displacement values: "
                 << displacementValues_
                 << ", velocityValues: " << velocityValues_;
      LOG(DEBUG) << "read and set Dirichlet BC: " << newDirichletBCValues;
    }

    //! set new dirichlet boundary condition values
    this->updateDirichletBoundaryConditions(nestedSolver, newDirichletBCValues);
  }

  template <typename SurfaceData>
  void setNeumannBoundaryConditions(SurfaceData &preciceData,
                                    NestedSolverType &nestedSolver,
                                    DihuContext &context) {
    // set traction values as neumann boundary conditions
    using FunctionSpace =
        typename PreciceAdapterNestedSolver<NestedSolverType>::FunctionSpace;
    using ElementWithFacesType =
        typename SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>::ElementWithFaces;

    std::vector<ElementWithFacesType> neumannBoundaryConditionElements;

    auto functionSpace = this->functionSpace(nestedSolver);

    // if this rank has no data, do not set any boundary conditions
    if (preciceData.preciceMesh->nNodesLocal != 0) {

      const int nNodesX = functionSpace->nNodesLocalWithoutGhosts(0);
      const int nNodesY = functionSpace->nNodesLocalWithoutGhosts(1);
      const int nNodesZ = functionSpace->nNodesLocalWithoutGhosts(2);

      int nElementsX = functionSpace->meshPartition()->nElementsLocal(0);
      int nElementsY = functionSpace->meshPartition()->nElementsLocal(1);
      int nElementsZ = functionSpace->meshPartition()->nElementsLocal(2);

      // set node and element indics in z direction for bottom surface
      int nodeIndexZ = 0;
      int elementalNodeIndexZ = 0;
      int elementIndexZ = 0;

      // for top surface
      if (preciceData.preciceMesh->face ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
        nodeIndexZ = nNodesZ - 1;
        elementIndexZ = nElementsZ - 1;
        elementalNodeIndexZ = 2;
      }

      // loop over elements
      for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++) {
        for (int elementIndexX = 0; elementIndexX < nElementsX;
             elementIndexX++) {
          ElementWithFacesType elementWithFaces;
          element_no_t elementNoLocal =
              elementIndexZ * nElementsX * nElementsY +
              elementIndexY * nElementsX + elementIndexX;
          elementWithFaces.elementNoLocal = elementNoLocal;

          // set surface dofs
          Mesh::face_t face = Mesh::face_t::face2Minus;
          if (preciceData.preciceMesh->face ==
              PreciceAdapterInitialize<
                  NestedSolverType>::PreciceSurfaceMesh::face2Plus)
            face = Mesh::face_t::face2Plus;

          elementWithFaces.face = face;

          // get dofs indices within the numbering of the volume element that
          // correspond to the selected face
          const int nDofsPerNode = FunctionSpace::nDofsPerNode();
          const int nSurfaceDofs =
              ::FunctionSpace::FunctionSpaceBaseDim<
                  2,
                  typename FunctionSpace::BasisFunction>::nNodesPerElement() *
              nDofsPerNode;
          std::array<dof_no_t, nSurfaceDofs> surfaceDofs;
          FunctionSpace::getFaceDofs(face, surfaceDofs);

          elementWithFaces.surfaceDofs.assign(surfaceDofs.begin(),
                                              surfaceDofs.end());

          // loop over the nodes of the element
          for (int elementalNodeIndexY = 0; elementalNodeIndexY < 3;
               elementalNodeIndexY++) {
            for (int elementalNodeIndexX = 0; elementalNodeIndexX < 3;
                 elementalNodeIndexX++) {
              int elementalDofIndex = elementalNodeIndexZ * 9 +
                                      elementalNodeIndexY * 3 +
                                      elementalNodeIndexX;

              dof_no_t dofNoLocal =
                  functionSpace->getDofNo(elementNoLocal, elementalDofIndex);

              // only iterate over local dofs, the ghost dofs are not considered
              // here (although it would be correct to communicate them with
              // precice)
              if (dofNoLocal >= functionSpace->nDofsLocalWithoutGhosts())
                continue;

              int valueIndex = dofNoLocal - nodeIndexZ * nNodesX * nNodesY;

              LOG(DEBUG) << "(x,y,z)=(" << elementalNodeIndexX << ","
                         << elementalNodeIndexY << "," << elementalNodeIndexZ
                         << ") dofNoLocal " << dofNoLocal
                         << ", valueIndex: " << valueIndex << "/"
                         << tractionValues_.size() / 3;

              assert(valueIndex >= 0);
              if (valueIndex >= preciceData.preciceMesh->nNodesLocal) {
                LOG(ERROR) << "valueIndex: " << valueIndex
                           << ", dofNoLocal: " << dofNoLocal
                           << ", nNodes: " << nNodesX << "," << nNodesY
                           << ", nodeIndexZ: " << nodeIndexZ
                           << ", nNodesLocal: "
                           << preciceData.preciceMesh->nNodesLocal
                           << " elementNoLocal: " << elementNoLocal
                           << ", elementalDofIndex: " << elementalDofIndex
                           << ", elementalNodeIndexZ: " << elementalNodeIndexZ
                           << ", meshPartition: "
                           << *functionSpace->meshPartition();
              }
              assert(valueIndex < preciceData.preciceMesh->nNodesLocal);

              Vec3 traction;
              for (int i = 0; i < 3; i++) {
                traction[i] = tractionValues_[3 * valueIndex + i];
              }

              dof_no_t surfaceDof = elementalDofIndex;

              // for top surface
              if (preciceData.preciceMesh->face ==
                  PreciceAdapterInitialize<
                      NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
                surfaceDof = 18 + elementalDofIndex;
              }
              elementWithFaces.dofVectors.push_back(
                  std::pair<dof_no_t, Vec3>(surfaceDof, traction));

              // LOG(INFO) << "dofVectors: " << elementWithFaces.dofVectors <<
              // ", traction: " << traction;
            }
          }
          neumannBoundaryConditionElements.push_back(elementWithFaces);
        }
      }
    }
    // create new Neumann BC object
    using NeumannBoundaryConditionsType =
        SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>;
    std::shared_ptr<NeumannBoundaryConditionsType> neumannBoundaryConditions =
        std::make_shared<NeumannBoundaryConditionsType>(context);
    neumannBoundaryConditions->initialize(functionSpace,
                                          neumannBoundaryConditionElements);
    neumannBoundaryConditions->setDeformationGradientField(
        this->deformationGradientField(nestedSolver));

#ifndef NDEBUG
    std::stringstream s;
    for (int i = 0; i < neumannBoundaryConditionElements.size(); i++) {
      s << "{el. " << neumannBoundaryConditionElements[i].elementNoLocal
        << ", \"" << Mesh::getString(neumannBoundaryConditionElements[i].face)
        << "\", dofVectors: [";
      for (int j = 0; j < neumannBoundaryConditionElements[i].dofVectors.size();
           j++) {
        s << "(" << neumannBoundaryConditionElements[i].dofVectors[j].first
          << ": (";
        for (int k = 0;
             k <
             neumannBoundaryConditionElements[i].dofVectors[j].second.size();
             k++) {
          if (k != 0)
            s << ",";
          s << neumannBoundaryConditionElements[i].dofVectors[j].second[k];
        }
        s << ")),";
      }
      s << "], surfaceDofs: ";
      for (int j = 0;
           j < neumannBoundaryConditionElements[i].surfaceDofs.size(); j++)
        s << neumannBoundaryConditionElements[i].surfaceDofs[j] << ",";
      s << "} ";
    }
    LOG(DEBUG) << "read and set Neumann BC:\n" << s.str();
#endif

    // set Neumann BCs in the static hyperelasticity of the
    // TimeSteppingScheme::DynamicHyperelasticitySolver solver
    this->updateNeumannBoundaryConditions(nestedSolver,
                                          neumannBoundaryConditions);

    LOG(DEBUG) << "read data from precice complete, traction values: "
               << tractionValues_;
  }
  //! initialize dirichlet boundary conditions by adding prescribed values for
  //! all bottom or top nodes
  void addDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<
          typename SpatialDiscretization::DirichletBoundaryConditionsBase<
              FunctionSpace, 6>::ElementWithNodes>
          &dirichletBoundaryConditionElements);

  //! update existing boundary conditions with new values
  void updateDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<std::pair<global_no_t, std::array<double, 6>>>
          newDirichletBoundaryConditionValues);

  //! update the neumann boundary conditions by replacing the complete object
  void updateNeumannBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>>
          neumannBoundaryConditions);

  //! get the displacement and velocity vectors of the given local dof nos
  void getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                     const std::vector<dof_no_t> &dofNosLocal,
                                     std::vector<double> &displacementValues,
                                     std::vector<double> &velocityValues);

  //! get the traction vectors of the given local dof nos
  void getTractionValues(NestedSolverType &nestedSolver,
                         const std::vector<dof_no_t> &dofNosLocal,
                         std::vector<double> &tractionValues);

  //! get at Petsc Vec that stores all values of the current state, to be used
  //! to store and restore checkpoints
  Vec currentState(NestedSolverType &nestedSolver);

  //! get the field variable of the deformation gradient
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 9>>
  deformationGradientField(NestedSolverType &nestedSolver);

  void reset(NestedSolverType &nestedSolver);

  void saveFiberData(NestedSolverType &nestedSolver);

  void loadFiberData(NestedSolverType &nestedSolver);
};

/** Partial specialization for tendon or pure mechanics solver, static nonlinear
 * elasticity
 */
template <typename Material>
class PreciceAdapterNestedSolver<
    SpatialDiscretization::HyperelasticitySolver<Material>> {
public:
  std::vector<double> displacementValues_;
  std::vector<double> velocityValues_;
  std::vector<double> tractionValues_;
  std::vector<Vec3> displacementVectors_;
  std::vector<Vec3> velocityVectors_;
  std::vector<Vec3> tractionVectors_;
  std::vector<double> scalarValues_;
  std::vector<double> scalarValuesOfMesh_;
  std::vector<Vec3> geometryValues_;
  //! define the type of the nested solver
  typedef SpatialDiscretization::HyperelasticitySolver<Material>
      NestedSolverType;

  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolverType::FunctionSpace FunctionSpace;

  typedef typename SpatialDiscretization::DirichletBoundaryConditionsBase<
      FunctionSpace, 6>::ElementWithNodes ElementWithNodes;

  //! get the function space of the nested solver, after it has been initialized
  std::shared_ptr<typename SpatialDiscretization::HyperelasticitySolver<
      Material>::FunctionSpace>
  functionSpace(NestedSolverType &nestedSolver);

  //! parse the options in "preciceVolumeData" and initialize all variables in
  //! precice, store in variable preciceVolumeData_
  template <typename PreciceVolumeData, typename PreciceVolumeMesh>
  void initializePreciceVolumeData(
      PythonConfig &specificSettings, NestedSolverType &nestedSolver,
      std::shared_ptr<precice::Participant> &preciceParticipant,
      std::vector<PreciceVolumeData> &preciceVolumeData,
      std::vector<std::shared_ptr<PreciceVolumeMesh>> &preciceVolumeMeshes){};

  template <typename SurfaceDataVector, typename VolumeDataVector>
  void
  preciceReadData(NestedSolverType &nestedSolver,
                  std::shared_ptr<precice::Participant> &preciceParticipant,
                  SurfaceDataVector &preciceSurfaceData,
                  VolumeDataVector &preciceVolumeData, DihuContext &context) {
    LOG(DEBUG) << "read surface data from precice";
    double preciceDt = preciceParticipant->getMaxTimeStepSize();
    // loop over data
    for (auto &preciceData : preciceSurfaceData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceData::ioRead) {
        // allocate memory
        int nEntries = preciceData.preciceMesh->nNodesLocal * 3;

        // if the data is displacements and velocities
        if (!preciceData.displacementsName.empty()) {
          displacementValues_.resize(nEntries);
          velocityValues_.resize(nEntries);

          // get all data at once
          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.displacementsName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              displacementValues_);

          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.velocitiesName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              velocityValues_);

          setDirichletBoundaryConditions(preciceData, nestedSolver);
        }
        // if the data is traction
        else if (!preciceData.tractionName.empty()) {
          tractionValues_.resize(nEntries);
          preciceParticipant->readData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, preciceDt,
              tractionValues_);

          setNeumannBoundaryConditions(preciceData, nestedSolver, context);
        } else {
          LOG(FATAL) << "Unknown precice data (read), none of displacements, "
                        "velocities or traction is set.";
        }
      }
    };
  }

  template <typename SurfaceDataVector, typename VolumeDataVector>
  void
  preciceWriteData(NestedSolverType &nestedSolver,
                   std::shared_ptr<precice::Participant> &preciceParticipant,
                   SurfaceDataVector &preciceSurfaceData,
                   VolumeDataVector &preciceVolumeData, double scalingFactor) {
    // write data to precice
    LOG(DEBUG) << "write surface data to precice";

    // loop over data
    for (auto &preciceData : preciceSurfaceData) {
      if (preciceData.ioType ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceData::ioWrite) {
        // if the data is displacements and velocities
        if (!preciceData.displacementsName.empty()) {
          // convert geometry values to precice data layout
          displacementValues_.clear();
          velocityValues_.clear();

          this->getDisplacementVelocityValues(
              nestedSolver, preciceData.preciceMesh->dofNosLocal,
              displacementValues_, velocityValues_);

#ifndef NDEBUG
          LOG(DEBUG) << "write displacements data to precice: "
                     << displacementValues_;
#endif
          // scale displacement and velocity values
          for (double &value : displacementValues_)
            value *= scalingFactor;

          for (double &value : velocityValues_)
            value *= scalingFactor;

          // write displacement values in precice
          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.displacementsName,
              preciceData.preciceMesh->preciceVertexIds, displacementValues_);

          // write velocity values in precice
          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.velocitiesName,
              preciceData.preciceMesh->preciceVertexIds, velocityValues_);
        }
        // if the data is traction
        else if (!preciceData.tractionName.empty() &&
                 preciceData.average == true) {
          LOG(INFO) << "Write averaged traction";
          // convert geometry values to precice data layout
          tractionValues_.clear();
          this->getTractionValues(nestedSolver,
                                  preciceData.preciceMesh->dofNosLocal,
                                  tractionValues_);
          // average z-values of traction
          double average_traction = 0.0;
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            average_traction += tractionValues_[i];
          }
          average_traction /= (tractionValues_.size() / 3);
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            tractionValues_[i] = average_traction;
          }
#ifndef NDEBUG
          LOG(DEBUG) << "Write averaged-traction data to precice: "
                     << tractionValues_[2];
#endif
          // scale traction values, they are always scaled by the factor of -1
          for (double &value : tractionValues_) {
            value *= scalingFactor;
            value *= -1;
          }

          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, tractionValues_);
        } else if (!preciceData.tractionName.empty()) {

          LOG(INFO) << "Write non-averaged traction";

          // convert geometry values to precice data layout
          tractionValues_.clear();
          this->getTractionValues(nestedSolver,
                                  preciceData.preciceMesh->dofNosLocal,
                                  tractionValues_);

#ifndef NDEBUG
          LOG(DEBUG) << "write traction data to precice: " << tractionValues_;
          std::stringstream s;
          for (int i = 2; i < tractionValues_.size(); i += 3) {
            s << " " << tractionValues_[i];
          }
          LOG(DEBUG) << "z values of traction: " << s.str();
#endif
          // scale traction values, they are always scaled by the factor of -1
          for (double &value : tractionValues_) {
            value *= scalingFactor;
            value *= -1;
          }

          preciceParticipant->writeData(
              preciceData.preciceMesh->preciceMeshName,
              preciceData.tractionName,
              preciceData.preciceMesh->preciceVertexIds, tractionValues_);
        } else {
          LOG(FATAL) << "Unknown precice data (write), none of displacements, "
                        "velocities or traction is set.";
        }
      }
    }

    LOG(DEBUG) << "write traction data to precice complete";
  }

  template <typename SurfaceData>
  void setDirichletBoundaryConditions(SurfaceData &preciceData,
                                      NestedSolverType &nestedSolver) {
    std::vector<std::pair<global_no_t, std::array<double, 6>>>
        newDirichletBCValues;

    auto functionSpace = this->functionSpace(nestedSolver);

    // if this rank has no data, do not set any boundary conditions
    if (preciceData.preciceMesh->nNodesLocal != 0) {
      // loop over nodes
      const int nNodesX = functionSpace->nNodesLocalWithoutGhosts(0);
      const int nNodesY = functionSpace->nNodesLocalWithoutGhosts(1);
      const int nNodesZ = functionSpace->nNodesLocalWithoutGhosts(2);

      // set node index in z direction for bottom surface
      int nodeIndexZ = 0;

      // for top surface
      if (preciceData.preciceMesh->face ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
        nodeIndexZ = nNodesZ - 1;
      }

      // loop over nodes to set the received values
      newDirichletBCValues.reserve(nNodesX * nNodesY);
      int valueIndex = 0;

      // loop over nodes of surface mesh
      for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++) {
        for (int nodeIndexX = 0; nodeIndexX < nNodesX;
             nodeIndexX++, valueIndex++) {
          node_no_t nodeNoLocal = nodeIndexZ * nNodesX * nNodesY +
                                  nodeIndexY * nNodesX + nodeIndexX;

          dof_no_t dofNoLocal = nodeNoLocal;
          global_no_t dofNoGlobal =
              functionSpace->meshPartition()->getDofNoGlobalPetsc(dofNoLocal);

          // assign received values to dirichlet bc vector of size 6
          std::array<double, 6> newDirichletBCValue;

          for (int i = 0; i < 3; i++) {
            newDirichletBCValue[i] = displacementValues_[3 * valueIndex + i];
            newDirichletBCValue[3 + i] = velocityValues_[3 * valueIndex + i];
          }

          newDirichletBCValues.push_back(
              std::pair<global_no_t, std::array<double, 6>>(
                  dofNoGlobal, newDirichletBCValue));
        }
      }

      LOG(DEBUG) << "read data from precice complete, displacement values: "
                 << displacementValues_
                 << ", velocityValues: " << velocityValues_;
      LOG(DEBUG) << "read and set Dirichlet BC: " << newDirichletBCValues;
    }

    //! set new dirichlet boundary condition values
    this->updateDirichletBoundaryConditions(nestedSolver, newDirichletBCValues);
  }

  template <typename SurfaceData>
  void setNeumannBoundaryConditions(SurfaceData &preciceData,
                                    NestedSolverType &nestedSolver,
                                    DihuContext &context) {
    // set traction values as neumann boundary conditions
    using FunctionSpace =
        typename PreciceAdapterNestedSolver<NestedSolverType>::FunctionSpace;
    using ElementWithFacesType =
        typename SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>::ElementWithFaces;

    std::vector<ElementWithFacesType> neumannBoundaryConditionElements;

    auto functionSpace = this->functionSpace(nestedSolver);

    // if this rank has no data, do not set any boundary conditions
    if (preciceData.preciceMesh->nNodesLocal != 0) {

      const int nNodesX = functionSpace->nNodesLocalWithoutGhosts(0);
      const int nNodesY = functionSpace->nNodesLocalWithoutGhosts(1);
      const int nNodesZ = functionSpace->nNodesLocalWithoutGhosts(2);

      int nElementsX = functionSpace->meshPartition()->nElementsLocal(0);
      int nElementsY = functionSpace->meshPartition()->nElementsLocal(1);
      int nElementsZ = functionSpace->meshPartition()->nElementsLocal(2);

      // set node and element indics in z direction for bottom surface
      int nodeIndexZ = 0;
      int elementalNodeIndexZ = 0;
      int elementIndexZ = 0;

      // for top surface
      if (preciceData.preciceMesh->face ==
          PreciceAdapterInitialize<
              NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
        nodeIndexZ = nNodesZ - 1;
        elementIndexZ = nElementsZ - 1;
        elementalNodeIndexZ = 2;
      }

      // loop over elements
      for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++) {
        for (int elementIndexX = 0; elementIndexX < nElementsX;
             elementIndexX++) {
          ElementWithFacesType elementWithFaces;
          element_no_t elementNoLocal =
              elementIndexZ * nElementsX * nElementsY +
              elementIndexY * nElementsX + elementIndexX;
          elementWithFaces.elementNoLocal = elementNoLocal;

          // set surface dofs
          Mesh::face_t face = Mesh::face_t::face2Minus;
          if (preciceData.preciceMesh->face ==
              PreciceAdapterInitialize<
                  NestedSolverType>::PreciceSurfaceMesh::face2Plus)
            face = Mesh::face_t::face2Plus;

          elementWithFaces.face = face;

          // get dofs indices within the numbering of the volume element that
          // correspond to the selected face
          const int nDofsPerNode = FunctionSpace::nDofsPerNode();
          const int nSurfaceDofs =
              ::FunctionSpace::FunctionSpaceBaseDim<
                  2,
                  typename FunctionSpace::BasisFunction>::nNodesPerElement() *
              nDofsPerNode;
          std::array<dof_no_t, nSurfaceDofs> surfaceDofs;
          FunctionSpace::getFaceDofs(face, surfaceDofs);

          elementWithFaces.surfaceDofs.assign(surfaceDofs.begin(),
                                              surfaceDofs.end());

          // loop over the nodes of the element
          for (int elementalNodeIndexY = 0; elementalNodeIndexY < 3;
               elementalNodeIndexY++) {
            for (int elementalNodeIndexX = 0; elementalNodeIndexX < 3;
                 elementalNodeIndexX++) {
              int elementalDofIndex = elementalNodeIndexZ * 9 +
                                      elementalNodeIndexY * 3 +
                                      elementalNodeIndexX;

              dof_no_t dofNoLocal =
                  functionSpace->getDofNo(elementNoLocal, elementalDofIndex);

              // only iterate over local dofs, the ghost dofs are not considered
              // here (although it would be correct to communicate them with
              // precice)
              if (dofNoLocal >= functionSpace->nDofsLocalWithoutGhosts())
                continue;

              int valueIndex = dofNoLocal - nodeIndexZ * nNodesX * nNodesY;

              LOG(DEBUG) << "(x,y,z)=(" << elementalNodeIndexX << ","
                         << elementalNodeIndexY << "," << elementalNodeIndexZ
                         << ") dofNoLocal " << dofNoLocal
                         << ", valueIndex: " << valueIndex << "/"
                         << tractionValues_.size() / 3;

              assert(valueIndex >= 0);
              if (valueIndex >= preciceData.preciceMesh->nNodesLocal) {
                LOG(ERROR) << "valueIndex: " << valueIndex
                           << ", dofNoLocal: " << dofNoLocal
                           << ", nNodes: " << nNodesX << "," << nNodesY
                           << ", nodeIndexZ: " << nodeIndexZ
                           << ", nNodesLocal: "
                           << preciceData.preciceMesh->nNodesLocal
                           << " elementNoLocal: " << elementNoLocal
                           << ", elementalDofIndex: " << elementalDofIndex
                           << ", elementalNodeIndexZ: " << elementalNodeIndexZ
                           << ", meshPartition: "
                           << *functionSpace->meshPartition();
              }
              assert(valueIndex < preciceData.preciceMesh->nNodesLocal);

              Vec3 traction;
              for (int i = 0; i < 3; i++) {
                traction[i] = tractionValues_[3 * valueIndex + i];
              }

              dof_no_t surfaceDof = elementalDofIndex;

              // for top surface
              if (preciceData.preciceMesh->face ==
                  PreciceAdapterInitialize<
                      NestedSolverType>::PreciceSurfaceMesh::face2Plus) {
                surfaceDof = 18 + elementalDofIndex;
              }
              elementWithFaces.dofVectors.push_back(
                  std::pair<dof_no_t, Vec3>(surfaceDof, traction));

              // LOG(INFO) << "dofVectors: " << elementWithFaces.dofVectors <<
              // ", traction: " << traction;
            }
          }
          neumannBoundaryConditionElements.push_back(elementWithFaces);
        }
      }
    }
    // create new Neumann BC object
    using NeumannBoundaryConditionsType =
        SpatialDiscretization::NeumannBoundaryConditions<
            FunctionSpace, Quadrature::Gauss<3>, 3>;
    std::shared_ptr<NeumannBoundaryConditionsType> neumannBoundaryConditions =
        std::make_shared<NeumannBoundaryConditionsType>(context);
    neumannBoundaryConditions->initialize(functionSpace,
                                          neumannBoundaryConditionElements);
    neumannBoundaryConditions->setDeformationGradientField(
        this->deformationGradientField(nestedSolver));

#ifndef NDEBUG
    std::stringstream s;
    for (int i = 0; i < neumannBoundaryConditionElements.size(); i++) {
      s << "{el. " << neumannBoundaryConditionElements[i].elementNoLocal
        << ", \"" << Mesh::getString(neumannBoundaryConditionElements[i].face)
        << "\", dofVectors: [";
      for (int j = 0; j < neumannBoundaryConditionElements[i].dofVectors.size();
           j++) {
        s << "(" << neumannBoundaryConditionElements[i].dofVectors[j].first
          << ": (";
        for (int k = 0;
             k <
             neumannBoundaryConditionElements[i].dofVectors[j].second.size();
             k++) {
          if (k != 0)
            s << ",";
          s << neumannBoundaryConditionElements[i].dofVectors[j].second[k];
        }
        s << ")),";
      }
      s << "], surfaceDofs: ";
      for (int j = 0;
           j < neumannBoundaryConditionElements[i].surfaceDofs.size(); j++)
        s << neumannBoundaryConditionElements[i].surfaceDofs[j] << ",";
      s << "} ";
    }
    LOG(DEBUG) << "read and set Neumann BC:\n" << s.str();
#endif

    // set Neumann BCs in the static hyperelasticity of the
    // TimeSteppingScheme::DynamicHyperelasticitySolver solver
    this->updateNeumannBoundaryConditions(nestedSolver,
                                          neumannBoundaryConditions);

    LOG(DEBUG) << "read data from precice complete, traction values: "
               << tractionValues_;
  }
  //! initialize dirichlet boundary conditions by adding prescribed values for
  //! all bottom or top nodes
  void addDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<
          typename SpatialDiscretization::DirichletBoundaryConditionsBase<
              FunctionSpace, 6>::ElementWithNodes>
          &dirichletBoundaryConditionElements);

  //! update existing boundary conditions with new values
  void updateDirichletBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::vector<std::pair<global_no_t, std::array<double, 6>>>
          newDirichletBoundaryConditionValues);

  //! update the neumann boundary conditions by replacing the complete object
  void updateNeumannBoundaryConditions(
      NestedSolverType &nestedSolver,
      std::shared_ptr<SpatialDiscretization::NeumannBoundaryConditions<
          FunctionSpace, Quadrature::Gauss<3>, 3>>
          neumannBoundaryConditions);

  //! get the displacement and velocity vectors of the given local dof nos
  void getDisplacementVelocityValues(NestedSolverType &nestedSolver,
                                     const std::vector<dof_no_t> &dofNosLocal,
                                     std::vector<double> &displacementValues,
                                     std::vector<double> &velocityValues);

  //! get the traction vectors of the given local dof nos
  void getTractionValues(NestedSolverType &nestedSolver,
                         const std::vector<dof_no_t> &dofNosLocal,
                         std::vector<double> &tractionValues);

  //! get at Petsc Vec that stores all values of the current state, to be used
  //! to store and restore checkpoints
  Vec currentState(NestedSolverType &nestedSolver);

  //! get the field variable of the deformation gradient
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace, 9>>
  deformationGradientField(NestedSolverType &nestedSolver);

  void reset(NestedSolverType &nestedSolver);

  void saveFiberData(NestedSolverType &nestedSolver);

  void loadFiberData(NestedSolverType &nestedSolver);
};

} // namespace Control

#include "control/precice/nested_solver.tpp"