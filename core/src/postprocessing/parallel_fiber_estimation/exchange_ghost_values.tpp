#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
exchangeGhostValues(const std::array<bool,4> &subdomainIsAtBorder)
{
  struct GhostValues
  {
    std::vector<double> nodePositionValues;   ///< the values of the node positions of the ghost elements, in sequential order like [x,x,x,x, ... y,y,y,y ... z,z,z,...]
    std::vector<double> solutionValues;      ///< values of solution field variable inside the ghost elements
    std::vector<double> gradientValues;      ///< values of the gradient field variable, consecutive like nodePositionValues
    std::array<element_no_t,3> nElementsPerCoordinateDirection;   ///< size of the ghost mesh
  };
  std::array<GhostValues,4> ghostValuesBuffer;  ///< [face], data for meshes containing ghost elements for the sides, face0Minus, face0Plus, face1Minus, face1Plus
  std::array<GhostValues,4> boundaryValues;  ///< [face], data to be send

  // communicate ghost elements to neighbouring subdomains
  // determine elements on own domain that are to be send to the neighbouring domain
  // loop over faces: face0Minus = 0, face0Plus, face1Minus, face1Plus
  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    if (!subdomainIsAtBorder[face])
    {
      LOG(DEBUG) << "determine ghost elements for face " << Mesh::getString((Mesh::face_t)face);

      // get information about neighbouring rank and boundary elements for face
      int neighbourRankNo;
      std::vector<dof_no_t> dofNos;
      meshPartition_->getBoundaryElements((Mesh::face_t)face, neighbourRankNo, ghostValuesBuffer[face].nElementsPerCoordinateDirection, dofNos);

      LOG(DEBUG) << "getBoundaryElements for face " << Mesh::getString((Mesh::face_t)face) << " returned neighbourRankNo=" << neighbourRankNo
        << ", nElementsPerCoordinateDirection: " << ghostValuesBuffer[face].nElementsPerCoordinateDirection << ", " << dofNos.size();

      // get relevant values in own domain that will be send to the neighbouring domain
      problem_->data().functionSpace()->geometryField().setRepresentationGlobal();
      problem_->data().solution()->setRepresentationGlobal();
      data_.gradient()->setRepresentationGlobal();
      problem_->data().functionSpace()->geometryField().startGhostManipulation();
      problem_->data().solution()->startGhostManipulation();
      data_.gradient()->startGhostManipulation();

      problem_->data().functionSpace()->geometryField().getValues(dofNos, boundaryValues[face].nodePositionValues);
      problem_->data().solution()->getValues(dofNos, boundaryValues[face].solutionValues);
      data_.gradient()->getValues(dofNos, boundaryValues[face].gradientValues);

      int nNodePositionValues = boundaryValues[face].nodePositionValues.size();
      int nSolutionValues = boundaryValues[face].solutionValues.size();
      int nGradientValues = boundaryValues[face].gradientValues.size();
      assert(nSolutionValues*3 == nGradientValues);
      assert(nNodePositionValues == nGradientValues);

#if 1
      if (neighbourRankNo != -1)
      {
        std::stringstream s;
        s << "04_ghost_elements_face_" << Mesh::getString((Mesh::face_t)face);
        PyObject_CallFunction(functionOutputGhostElements_, "s i O O f", s.str().c_str(), currentRankSubset_->ownRankNo(),
                              PythonUtility::convertToPython<std::vector<double>>::get(boundaryValues[face].nodePositionValues),
                              PythonUtility::convertToPython<std::array<element_no_t,3>>::get(ghostValuesBuffer[face].nElementsPerCoordinateDirection), 0.05);
        PythonUtility::checkForError();
      }
#endif

      LOG(DEBUG) << "exchange ghosts with neighbour " << neighbourRankNo << ": nNodePositionValues=" << nNodePositionValues << ", nSolutionValues=" << nSolutionValues << ", nGradientValues=" << nGradientValues;

#if !defined(WRITE_CHECKPOINT_GHOST_MESH) && !defined(USE_CHECKPOINT_GHOST_MESH)
      for (int i = 0; i < 2; i++)
      {

        if (((face == Mesh::face_t::face0Minus || face == Mesh::face_t::face1Minus) && i == 0)
          ||((face == Mesh::face_t::face0Plus || face == Mesh::face_t::face1Plus) && i == 1))
        {
          LOG(DEBUG) << "receive from rank " << neighbourRankNo;

          // receive from neighbouring process
          // blocking receive call to receive node position values
          ghostValuesBuffer[face].nodePositionValues.resize(nNodePositionValues);
          MPIUtility::handleReturnValue(MPI_Recv(ghostValuesBuffer[face].nodePositionValues.data(), nNodePositionValues, MPI_DOUBLE,
                                                 neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

          // blocking receive call to receive solution values
          ghostValuesBuffer[face].solutionValues.resize(nSolutionValues);
          MPIUtility::handleReturnValue(MPI_Recv(ghostValuesBuffer[face].solutionValues.data(), nSolutionValues, MPI_DOUBLE,
                                                 neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

          // blocking receive call to receive gradient values
          ghostValuesBuffer[face].gradientValues.resize(nGradientValues);
          MPIUtility::handleReturnValue(MPI_Recv(ghostValuesBuffer[face].gradientValues.data(), nGradientValues, MPI_DOUBLE,
                                                 neighbourRankNo, 0, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

          LOG(DEBUG) << "receive from rank " << neighbourRankNo << " completed";
        }
        else
        {
          LOG(DEBUG) << "send to rank " << neighbourRankNo;

          // send values to neighbouring process
          // blocking send call to send solution values
          MPIUtility::handleReturnValue(MPI_Send(boundaryValues[face].nodePositionValues.data(), nNodePositionValues, MPI_DOUBLE,
                                                 neighbourRankNo, 0, currentRankSubset_->mpiCommunicator()), "MPI_Send");

          // blocking send call to send solution values
          MPIUtility::handleReturnValue(MPI_Send(boundaryValues[face].solutionValues.data(), nSolutionValues, MPI_DOUBLE,
                                                 neighbourRankNo, 0, currentRankSubset_->mpiCommunicator()), "MPI_Send");

          // blocking send call to send gradient values
          MPIUtility::handleReturnValue(MPI_Send(boundaryValues[face].gradientValues.data(), nGradientValues, MPI_DOUBLE,
                                                 neighbourRankNo, 0, currentRankSubset_->mpiCommunicator()), "MPI_Send");

          LOG(DEBUG) << "send to rank " << neighbourRankNo << " completed";
        }
      }
#endif

#ifdef WRITE_CHECKPOINT_GHOST_MESH
      std::stringstream filenameOut;
      filenameOut << "checkpoint_ghost_mesh_" << neighbourRankNo << "_" << Mesh::getString((Mesh::face_t)face) << ".txt";
      std::ofstream fileOut;
      fileOut.open(filenameOut.str().c_str(), std::ios::out | std::ios::trunc);

      assert(fileOut.is_open());

      fileOut << ghostValuesBuffer[face].nodePositionValues.size() << " "
        << ghostValuesBuffer[face].solutionValues.size() << " "
        << ghostValuesBuffer[face].gradientValues.size() << " ";


      for (int i = 0; i < ghostValuesBuffer[face].nodePositionValues.size(); i++)
        fileOut << ghostValuesBuffer[face].nodePositionValues[i] << " ";

      for (int i = 0; ghostValuesBuffer[face].solutionValues.size(); i++)
        fileOut << ghostValuesBuffer[face].solutionValues[i] << " ";

      for (int i = 0; i < ghostValuesBuffer[face].gradientValues.size(); i++)
        fileOut << ghostValuesBuffer[face].gradientValues[i] << " ";

      fileOut.close();
      LOG(DEBUG) << "Wrote ghost meshes (" << ghostValuesBuffer[face].nodePositionValues.size() << " node position values, "
        << ghostValuesBuffer[face].solutionValues.size() << " solution values, " << ghostValuesBuffer[face].gradientValues.size()
        << " gradient values) to checkpoint file " << filenameOut.str();
#endif

#ifdef USE_CHECKPOINT_GHOST_MESH

      std::stringstream filenameIn;
      filenameIn << "checkpoint_ghost_mesh_" << neighbourRankNo << "_" << Mesh::getString((Mesh::face_t)face) << ".txt";
      std::ifstream fileIn;
      fileIn.open(filenameIn.str().c_str(), std::ios::in);

      assert(fileIn.is_open());

      int size1, size2, size3;
      fileIn >> size1 >> size2 >> size3;
      ghostValuesBuffer[face].nodePositionValues.resize(size1);
      ghostValuesBuffer[face].solutionValues.resize(size2);
      ghostValuesBuffer[face].gradientValues.resize(size3);

      for (int i = 0; i < ghostValuesBuffer[face].nodePositionValues.size(); i++)
        fileIn >> ghostValuesBuffer[face].nodePositionValues[i];

      for (int i = 0; i < ghostValuesBuffer[face].solutionValues.size(); i++)
        fileIn >> ghostValuesBuffer[face].solutionValues[i];

      for (int i = 0; i < ghostValuesBuffer[face].gradientValues.size(); i++)
        fileIn >> ghostValuesBuffer[face].gradientValues[i];

      fileIn.close();
      LOG(DEBUG) << "Loaded ghost meshes (" << ghostValuesBuffer[face].nodePositionValues.size() << " node position values, "
        << ghostValuesBuffer[face].solutionValues.size() << " solution values, " << ghostValuesBuffer[face].gradientValues.size()
        << " gradient values) from checkpoint file " << filenameIn.str();
#endif

#if 1
      if (neighbourRankNo != -1)
      {
        std::stringstream s;
        s << "05_received_ghost_elements_face_" << Mesh::getString((Mesh::face_t)face);
        PyObject_CallFunction(functionOutputGhostElements_, "s i O O f", s.str().c_str(), currentRankSubset_->ownRankNo(),
                              PythonUtility::convertToPython<std::vector<double>>::get(ghostValuesBuffer[face].nodePositionValues),
                              PythonUtility::convertToPython<std::array<element_no_t,3>>::get(ghostValuesBuffer[face].nElementsPerCoordinateDirection), 0.05);
        PythonUtility::checkForError();
      }
#endif
    }

  }

  // create ghost element meshes from received data

  LOG(DEBUG) << "split rank subset, own rank: " << currentRankSubset_->ownRankNo();
  std::shared_ptr<Partition::RankSubset> rankSubsetSingleRank = std::make_shared<Partition::RankSubset>(currentRankSubset_->ownRankNo(), currentRankSubset_->mpiCommunicator());

  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    if (!ghostValuesBuffer[face].nodePositionValues.empty())
    {
      std::stringstream meshName;
      meshName << this->functionSpace_->meshName() << "_ghost_" << Mesh::getString((Mesh::face_t)face);

      // transform the node position data from vector of double to vector of Vec3
      int nNodes = ghostValuesBuffer[face].nodePositionValues.size() / 3;
      std::vector<Vec3> nodePositions(nNodes);

      for (node_no_t nodeIndex = 0; nodeIndex < nNodes; nodeIndex++)
      {
        nodePositions[nodeIndex][0] = ghostValuesBuffer[face].nodePositionValues[nodeIndex*3 + 0];
        nodePositions[nodeIndex][1] = ghostValuesBuffer[face].nodePositionValues[nodeIndex*3 + 1];
        nodePositions[nodeIndex][2] = ghostValuesBuffer[face].nodePositionValues[nodeIndex*3 + 2];
      }

      LOG(DEBUG) << "create ghost mesh with nElementsPerCoordinateDirection: " << ghostValuesBuffer[face].nElementsPerCoordinateDirection;

      // create ghost mesh
      std::array<int,3> nRanks({1,1,1});
      context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetSingleRank);
      this->ghostMesh_[face] = context_.meshManager()->template createFunctionSpace<FunctionSpaceType>(
        meshName.str(), nodePositions, ghostValuesBuffer[face].nElementsPerCoordinateDirection, nRanks);


      LOG(DEBUG) << "created ghost mesh for face " << Mesh::getString((Mesh::face_t)face) << ".";

      // create solution field variable for ghost mesh
      LOG(DEBUG) << "create solution field variable on ghost mesh";
      this->ghostMeshSolution_[face] = this->ghostMesh_[face]->template createFieldVariable<1>("solution");
      this->ghostMeshSolution_[face]->setValuesWithGhosts(ghostValuesBuffer[face].solutionValues);

      LOG(DEBUG) << "create gradient field variable on ghost mesh";
      // create gradient field variable for ghost mesh
      this->ghostMeshGradient_[face] = this->ghostMesh_[face]->template createFieldVariable<3>("gradient");
      int nGradientValues = ghostValuesBuffer[face].gradientValues.size()/3;

      std::vector<Vec3> gradientValues(nGradientValues);
      for (int i = 0; i < nGradientValues; i++)
      {
        gradientValues[i][0] = ghostValuesBuffer[face].gradientValues[3*i + 0];
        gradientValues[i][1] = ghostValuesBuffer[face].gradientValues[3*i + 1];
        gradientValues[i][2] = ghostValuesBuffer[face].gradientValues[3*i + 2];
      }
      this->ghostMeshGradient_[face]->setValuesWithGhosts(gradientValues);
    }
    else
    {
      this->ghostMesh_[face] = nullptr;
      LOG(DEBUG) << "ghost mesh for face " << Mesh::getString((Mesh::face_t)face) << " is null, because ghostValuesBuffer[face].nodePositionValues is empty()";
    }

    this->functionSpace_->setGhostMesh((Mesh::face_t)face, this->ghostMesh_[face]);
    LOG(DEBUG) << "after settings ghost mesh for face " << Mesh::getString((Mesh::face_t)face) << ": ";
    this->functionSpace_->debugOutputGhostMeshSet();
  }

}

};  // namespace
