#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include "partition/rank_subset.h"
#include "control/diagnostic_tool/stimulation_logging.h"

//! get element lengths and vmValues from the other ranks
template <int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates, nAlgebraics,
                              DiffusionTimeSteppingScheme>::fetchFiberData() {
  VLOG(1) << "fetchFiberData";
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances =
      nestedSolvers_.instancesLocal();

  CellmlAdapterType &cellmlAdapter =
      instances[0].timeStepping1().instancesLocal()[0].discretizableInTime();

  int nInstancesLocalCellml;
  int nAlgebraicsLocalCellml;
  int nParametersPerInstance;
  cellmlAdapter.getNumbers(nInstancesLocalCellml, nAlgebraicsLocalCellml,
                           nParametersPerInstance);

  // loop over fibers and communicate element lengths and initial values to the
  // ranks that participate in computing
  int fiberNo = 0;
  int fiberDataNo = 0;
  for (int i = 0; i < instances.size(); i++) {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances =
        instances[i]
            .timeStepping1()
            .instancesLocal(); // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++) {
      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace =
          innerInstances[j].data().functionSpace();
      LOG(DEBUG) << "instance (inner,outer)=(" << i << "," << j
                 << "), fiberNo: " << fiberNo << ", functionSpace "
                 << fiberFunctionSpace->meshName() << ","
                 << fiberFunctionSpace->meshPartition()->rankSubset()->size()
                 << " ranks ("
                 << *fiberFunctionSpace->meshPartition()->rankSubset() << ")"
                 << ", " << innerInstances.size() << " inner instances";

      // communicate element lengths
      std::vector<double> localLengths(fiberFunctionSpace->nElementsLocal());

      // loop over local elements and compute element lengths
      for (element_no_t elementNoLocal = 0;
           elementNoLocal < fiberFunctionSpace->nElementsLocal();
           elementNoLocal++) {
        std::array<Vec3, FiberFunctionSpace::nDofsPerElement()>
            geometryElementValues;
        fiberFunctionSpace->geometryField().getElementValues(
            elementNoLocal, geometryElementValues);
        double elementLength = MathUtility::distance<3>(
            geometryElementValues[0], geometryElementValues[1]);
        localLengths[elementNoLocal] = elementLength;
      }

      std::shared_ptr<Partition::RankSubset> rankSubset =
          fiberFunctionSpace->meshPartition()->rankSubset();
      MPI_Comm mpiCommunicator = rankSubset->mpiCommunicator();
      int computingRank = fiberNo % rankSubset->size();

      std::vector<int> nElementsOnRanks(rankSubset->size());
      std::vector<int> nDofsOnRanks(rankSubset->size());
      std::vector<int> offsetsOnRanks(rankSubset->size());
      std::vector<int> nParametersOnRanks(rankSubset->size());
      std::vector<int> parameterOffsetsOnRanks(rankSubset->size());

      double *elementLengthsReceiveBuffer = nullptr;
      double *vmValuesReceiveBuffer = nullptr;
      std::vector<double> parametersReceiveBuffer(
          fiberFunctionSpace->nDofsGlobal() * nParametersPerInstance);

      if (computingRank == rankSubset->ownRankNo()) {
        // allocate buffers
        fiberData_[fiberDataNo].elementLengths.resize(
            fiberFunctionSpace->nElementsGlobal());
        fiberData_[fiberDataNo].vmValues.resize(
            fiberFunctionSpace->nDofsGlobal());

        // resize buffer of further data that will be transferred back in
        // updateFiberData()
        int nStatesAndAlgebraicsValues = statesForTransferIndices_.size() +
                                         algebraicsForTransferIndices_.size() -
                                         1;

        // if also the computeStateInformation should be output, increase buffer
        // size by one entry per node
        if (setComputeStateInformation_)
          nStatesAndAlgebraicsValues++;

        fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues.resize(
            fiberFunctionSpace->nDofsGlobal() * nStatesAndAlgebraicsValues);

        elementLengthsReceiveBuffer =
            fiberData_[fiberDataNo].elementLengths.data();
        vmValuesReceiveBuffer = fiberData_[fiberDataNo].vmValues.data();
      }

      for (int rankNo = 0; rankNo < rankSubset->size(); rankNo++) {
        nElementsOnRanks[rankNo] =
            fiberFunctionSpace->meshPartition()->nNodesLocalWithGhosts(0,
                                                                       rankNo) -
            1;
        offsetsOnRanks[rankNo] =
            fiberFunctionSpace->meshPartition()->beginNodeGlobalNatural(0,
                                                                        rankNo);
        nDofsOnRanks[rankNo] =
            fiberFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(
                0, rankNo);
        parameterOffsetsOnRanks[rankNo] =
            fiberFunctionSpace->meshPartition()->beginNodeGlobalNatural(
                0, rankNo) *
            nParametersPerInstance;
        nParametersOnRanks[rankNo] =
            fiberFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(
                0, rankNo) *
            nParametersPerInstance;
      }

      // int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype
      // sendtype,
      //          void *recvbuf, const int *recvcounts, const int *displs,
      //          MPI_Datatype recvtype, int root, MPI_Comm comm)
      //
      VLOG(1) << "Gatherv of element lengths to rank " << computingRank
              << ", values " << localLengths << ", sizes: " << nElementsOnRanks
              << ", offsets: " << offsetsOnRanks;

      MPI_Gatherv(localLengths.data(), fiberFunctionSpace->nElementsLocal(),
                  MPI_DOUBLE, elementLengthsReceiveBuffer,
                  nElementsOnRanks.data(), offsetsOnRanks.data(), MPI_DOUBLE,
                  computingRank, mpiCommunicator);

      // get own vm values
      std::vector<double> vmValuesLocal;
      innerInstances[j].data().solution()->getValuesWithoutGhosts(
          0, vmValuesLocal);

      // communicate Vm values
      LOG(DEBUG) << "Gatherv of values to rank " << computingRank
                 << ", sizes: " << nDofsOnRanks
                 << ", offsets: " << offsetsOnRanks << ", local values "
                 << vmValuesLocal;

      MPI_Gatherv(
          vmValuesLocal.data(), fiberFunctionSpace->nDofsLocalWithoutGhosts(),
          MPI_DOUBLE, vmValuesReceiveBuffer, nDofsOnRanks.data(),
          offsetsOnRanks.data(), MPI_DOUBLE, computingRank, mpiCommunicator);

      // communicate parameter values
      // get own parameter values

      // get the data_.parameters() raw pointer
      innerInstances[j].discretizableInTime().data().prepareParameterValues();

      double *parameterValuesLocal =
          innerInstances[j].discretizableInTime().data().parameterValues();
      // size of this array is fiberFunctionSpace->nDofsLocalWithoutGhosts() *
      // nAlgebraics parameterValuesLocal has struct of array memory layout with
      // space for a total of nAlgebraics_ parameters [i0p0, i1p0, i2p0, ...
      // i0p1, i1p1, i2p1, ...]

      // only the actual parameter values should be sent, not the rest of the
      // parameters buffer therefore allocate a send buffer with the according
      // size
      int nParametersLocal = nParametersPerInstance *
                             fiberFunctionSpace->nDofsLocalWithoutGhosts();
      std::vector<double> parametersSendBuffer(nParametersLocal);

      // loop over the actual parameter values for every dof
      for (int dofNoLocal = 0;
           dofNoLocal < fiberFunctionSpace->nDofsLocalWithoutGhosts();
           dofNoLocal++) {
        for (int parameterNo = 0; parameterNo < nParametersPerInstance;
             parameterNo++) {
          // store parameter values to send buffer
          parametersSendBuffer[dofNoLocal * nParametersPerInstance +
                               parameterNo] = parameterValuesLocal
              [parameterNo * fiberFunctionSpace->nDofsLocalWithoutGhosts() +
               dofNoLocal];
        }
      }

      if (VLOG_IS_ON(1)) {
        VLOG(1) << "Gatherv of parameters to rank " << computingRank
                << ", send buffer: " << parametersSendBuffer << " contains "
                << nParametersLocal << " parameters "
                << " " << nParametersPerInstance << " per instances with "
                << fiberFunctionSpace->nDofsLocalWithoutGhosts()
                << " local instances.";
      }

      // receive data
      MPI_Gatherv(parametersSendBuffer.data(), nParametersLocal, MPI_DOUBLE,
                  parametersReceiveBuffer.data(), nParametersOnRanks.data(),
                  parameterOffsetsOnRanks.data(), MPI_DOUBLE, computingRank,
                  mpiCommunicator);

      // store result from parametersReceiveBuffer (for current fiber) to
      // fiberPointBuffersParameters_ (for a vc vector) loop over number of
      // instances of the problem on the current fiber
      if (computingRank == rankSubset->ownRankNo()) {
        int nInstancesOnFiber = fiberData_[fiberDataNo].vmValues.size();

        for (int instanceNo = 0; instanceNo < nInstancesOnFiber; instanceNo++) {
          // compute indices for fiberPointBuffersParameters_
          global_no_t valueIndexAllFibers =
              fiberData_[fiberDataNo].valuesOffset + instanceNo;

          if (useVc_) {
            global_no_t pointBuffersNo =
                valueIndexAllFibers / Vc::double_v::size();
            int entryNo = valueIndexAllFibers % Vc::double_v::size();

            // LOG(DEBUG) << "valueIndexAllFibers: " << valueIndexAllFibers <<
            // ", (" << pointBuffersNo << "," << entryNo << ")";

            // set all received parameter values for the current instance in the
            // correct slot in the vc vector of the current pointBuffer compute
            // buffer
            for (int parameterNo = 0; parameterNo < nParametersPerInstance;
                 parameterNo++) {
              fiberPointBuffersParameters_
                  [pointBuffersNo][parameterNo][entryNo] =
                      parametersReceiveBuffer[instanceNo *
                                                  nParametersPerInstance +
                                              parameterNo];
            }

            if (VLOG_IS_ON(1)) {
              if (entryNo == Vc::double_v::size() - 1) {
                VLOG(1) << "stored " << nParametersPerInstance
                        << " parameters in buffer no " << pointBuffersNo << ": "
                        << fiberPointBuffersParameters_[pointBuffersNo];
              }
            }
          } else {
            int instanceNoToCompute =
                fiberDataNo * nInstancesOnFiber + instanceNo;

            // set all received parameter values for the current instance in the
            // correct slot in the vc vector of the current pointBuffer compute
            // buffer
            for (int parameterNo = 0; parameterNo < nParametersPerInstance;
                 parameterNo++) {
              // gpuParameters_[parameterNo*nInstances + instanceNo]
              gpuParameters_[parameterNo * nInstancesToCompute_ +
                             instanceNoToCompute] =
                  parametersReceiveBuffer[instanceNo * nParametersPerInstance +
                                          parameterNo];
            }
          }
        }

        if (!useVc_) {
          int nElementsOnFiber = fiberFunctionSpace->nElementsGlobal();
          for (int elementNo = 0; elementNo < nElementsOnFiber; elementNo++) {
            gpuElementLengths_[fiberDataNo * nElementsOnFiber + elementNo] =
                elementLengthsReceiveBuffer[elementNo];
          }
        }
      }

      // get the data_.parameters() raw pointer
      innerInstances[j].discretizableInTime().data().restoreParameterValues();

      // increase index for fiberData_ struct
      if (computingRank == rankSubset->ownRankNo())
        fiberDataNo++;
    }
  }

  if (useVc_) {
    // copy Vm values to compute buffers
    for (int fiberDataNo = 0; fiberDataNo < fiberData_.size(); fiberDataNo++) {
      int nValues = fiberData_[fiberDataNo].vmValues.size();

      for (int valueNo = 0; valueNo < nValues; valueNo++) {
        global_no_t valueIndexAllFibers =
            fiberData_[fiberDataNo].valuesOffset + valueNo;

        global_no_t pointBuffersNo = valueIndexAllFibers / Vc::double_v::size();
        int entryNo = valueIndexAllFibers % Vc::double_v::size();

        fiberPointBuffers_[pointBuffersNo].states[0][entryNo] =
            fiberData_[fiberDataNo].vmValues[valueNo];
      }
    }
  }
}

//! send vmValues data from fiberData_ back to the fibers where it belongs to
//! and set in the respective field variable
template <int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates, nAlgebraics,
                              DiffusionTimeSteppingScheme>::updateFiberData() {
  // copy Vm and other states/algebraics from compute buffers to fiberData_
  if (useVc_) {
    // loop over vector of fibers that have been computed locally
    for (int fiberDataNo = 0; fiberDataNo < fiberData_.size(); fiberDataNo++) {
      int nValues = fiberData_[fiberDataNo].vmValues.size();

      // loop over all nodes on the current entire fiber
      for (int valueNo = 0; valueNo < nValues; valueNo++) {
        // compute indices to access fiberPointBuffers_ variable
        global_no_t valueIndexAllFibers =
            fiberData_[fiberDataNo].valuesOffset + valueNo;

        global_no_t pointBuffersNo = valueIndexAllFibers / Vc::double_v::size();
        int entryNo = valueIndexAllFibers % Vc::double_v::size();

        assert(statesForTransferIndices_.size() > 0);
        const int stateToTransfer =
            statesForTransferIndices_[0]; // transfer the first state value

        // collect first values of first state for transfer, which is the Vm
        // values, store under vmValues
        fiberData_[fiberDataNo].vmValues[valueNo] =
            fiberPointBuffers_[pointBuffersNo].states[stateToTransfer][entryNo];

        // loop over further states to transfer
        int furtherDataIndex = 0;
        for (int i = 1; i < statesForTransferIndices_.size();
             i++, furtherDataIndex++) {
          const int stateToTransfer = statesForTransferIndices_[i];

          // store further states to transfer under
          // furtherStatesAndAlgebraicsValues
          fiberData_[fiberDataNo]
              .furtherStatesAndAlgebraicsValues[furtherDataIndex * nValues +
                                                valueNo] =
              fiberPointBuffers_[pointBuffersNo]
                  .states[stateToTransfer][entryNo];
        }

        // loop over algebraics to transfer
        for (int i = 0; i < algebraicsForTransferIndices_.size();
             i++, furtherDataIndex++) {
          // store further algebraics to transfer under
          // furtherStatesAndAlgebraicsValues
          fiberData_[fiberDataNo]
              .furtherStatesAndAlgebraicsValues[furtherDataIndex * nValues +
                                                valueNo] =
              fiberPointBuffersAlgebraicsForTransfer_[pointBuffersNo][i]
                                                     [entryNo];
        }

        // add the information about whether the point is constant or
        // not_constant or neighbour_not_constant
        if (setComputeStateInformation_) {
          // also store under furtherStatesAndAlgebraicsValues
          fiberData_[fiberDataNo]
              .furtherStatesAndAlgebraicsValues[furtherDataIndex * nValues +
                                                valueNo] =
              fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo];

          // if fiber has not been stimulated, set value to -1
          if (onlyComputeIfHasBeenStimulated_ &&
              !fiberHasBeenStimulated_[fiberDataNo]) {
            fiberData_[fiberDataNo]
                .furtherStatesAndAlgebraicsValues[furtherDataIndex * nValues +
                                                  valueNo] = -1;
          }
          // the same value is set for all Vc::double_v::size() entries of the
          // Vc::double_v vector (different entryNo's)
        }
      }
      LOG(DEBUG) << "states and algebraics for transfer at fiberDataNo="
                 << fiberDataNo << ": "
                 << fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues;
      LOG(DEBUG)
          << "size: "
          << fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues.size()
          << ", nValues: " << nValues;
    }
  }

  LOG(TRACE) << "updateFiberData";
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances =
      nestedSolvers_.instancesLocal();

  // loop over fibers and communicate element lengths and initial values to the
  // ranks that participate in computing

  int fiberNo = 0;
  int fiberDataNo = 0;
  for (int i = 0; i < instances.size(); i++) {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances =
        instances[i]
            .timeStepping1()
            .instancesLocal(); // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++) {
      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace =
          innerInstances[j].data().functionSpace();
      // CellmlAdapterType &cellmlAdapter =
      // innerInstances[j].discretizableInTime();

      // prepare helper variables for Scatterv
      std::shared_ptr<Partition::RankSubset> rankSubset =
          fiberFunctionSpace->meshPartition()->rankSubset();
      MPI_Comm mpiCommunicator = rankSubset->mpiCommunicator();
      int computingRank =
          fiberNo % rankSubset->size(); // rank which computes the current fiber

      std::vector<int> nDofsOnRanks(rankSubset->size());
      std::vector<int> offsetsOnRanks(rankSubset->size());

      for (int rankNo = 0; rankNo < rankSubset->size(); rankNo++) {
        offsetsOnRanks[rankNo] =
            fiberFunctionSpace->meshPartition()->beginNodeGlobalNatural(0,
                                                                        rankNo);
        nDofsOnRanks[rankNo] =
            fiberFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(
                0, rankNo);
      }

      double *sendBufferVmValues = nullptr;
      if (computingRank == rankSubset->ownRankNo()) {
        sendBufferVmValues = fiberData_[fiberDataNo].vmValues.data();
      }

      //  int MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int
      //  *displs, MPI_Datatype sendtype,
      //                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
      //                  int root, MPI_Comm comm)
      // communicate Vm values
      std::vector<double> vmValuesLocal(
          fiberFunctionSpace->nDofsLocalWithoutGhosts());
      MPI_Scatterv(sendBufferVmValues, nDofsOnRanks.data(),
                   offsetsOnRanks.data(), MPI_DOUBLE, vmValuesLocal.data(),
                   fiberFunctionSpace->nDofsLocalWithoutGhosts(), MPI_DOUBLE,
                   computingRank, mpiCommunicator);

      VLOG(1) << "Scatterv from rank " << computingRank
              << ", sizes: " << nDofsOnRanks << ", offsets: " << offsetsOnRanks
              << ", received local values " << vmValuesLocal;

      // store Vm values in CellmlAdapter and diffusion FiniteElementMethod
      VLOG(1) << "fiber " << fiberDataNo << ", set values " << vmValuesLocal;
      innerInstances[j].data().solution()->setValuesWithoutGhosts(
          0, vmValuesLocal);
      instances[i]
          .timeStepping2()
          .instancesLocal()[j]
          .data()
          .solution()
          ->setValuesWithoutGhosts(0, vmValuesLocal);

      // ----------------------
      // communicate further states and algebraics that are selected by the
      // options "statesForTransfer" and "algebraicsForTransfer"

      std::vector<int> nValuesOnRanks(rankSubset->size());
      int nStatesAndAlgebraicsValues = statesForTransferIndices_.size() +
                                       algebraicsForTransferIndices_.size() - 1;

      // if also the computeStateInformation should be communicated, the buffer
      // has entry more per node
      if (setComputeStateInformation_)
        nStatesAndAlgebraicsValues++;

      // fill offset and size helper variables for MPI_Scatterv
      for (int rankNo = 0; rankNo < rankSubset->size(); rankNo++) {
        offsetsOnRanks[rankNo] =
            fiberFunctionSpace->meshPartition()->beginNodeGlobalNatural(0,
                                                                        rankNo);
        nValuesOnRanks[rankNo] =
            fiberFunctionSpace->meshPartition()->nNodesLocalWithoutGhosts(
                0, rankNo);
      }

      // receive buffer valuesLocal[furtherDataIndex * nValues + valueNo]
      std::vector<double> valuesLocal(
          fiberFunctionSpace->nDofsLocalWithoutGhosts() *
          nStatesAndAlgebraicsValues);
      std::vector<MPI_Request> scatterRequests(nStatesAndAlgebraicsValues);

      // loop over variable to transfer, because of the memory layout it is not
      // possible to do this with a single MPI_Scatterv
      for (int variableNo = 0; variableNo < nStatesAndAlgebraicsValues;
           variableNo++) {
        // get send buffer for MPI_Scatterv
        double *sendBuffer = nullptr;
        if (computingRank == rankSubset->ownRankNo()) {
          sendBuffer =
              fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues.data() +
              variableNo * fiberFunctionSpace->nDofsGlobal();
        }
        double *receiveBuffer =
            valuesLocal.data() +
            variableNo * fiberFunctionSpace->nDofsLocalWithoutGhosts();

        MPI_Iscatterv(sendBuffer, nValuesOnRanks.data(), offsetsOnRanks.data(),
                      MPI_DOUBLE, receiveBuffer,
                      fiberFunctionSpace->nDofsLocalWithoutGhosts(), MPI_DOUBLE,
                      computingRank, mpiCommunicator,
                      &scatterRequests[variableNo]);

        // debugging output
        if (VLOG_IS_ON(1)) {
          std::stringstream s;
          for (int i = 0; i < fiberFunctionSpace->nDofsLocalWithoutGhosts();
               i++)
            s << ","
              << valuesLocal[variableNo *
                                 fiberFunctionSpace->nDofsLocalWithoutGhosts() +
                             i];
          VLOG(1) << "Scatterv furtherStatesAndAlgebraicsValues from rank "
                  << computingRank << ", sizes: " << nValuesOnRanks
                  << ", offsets: " << offsetsOnRanks
                  << ", variableNo: " << variableNo
                  << " sendBuffer: " << sendBuffer
                  << ", received local values: (" << s.str() << ")";
        }
      }
      MPI_Waitall(nStatesAndAlgebraicsValues, scatterRequests.data(),
                  MPI_STATUSES_IGNORE);

      // store received states and algebraics values in diffusion
      // slotConnectorData loop over further states to transfer
      int furtherDataIndex = 0;
      for (int stateIndex = 1; stateIndex < statesForTransferIndices_.size();
           stateIndex++, furtherDataIndex++) {
        // store in diffusion

        // get field variable
        std::vector<::Data::ComponentOfFieldVariable<FiberFunctionSpace, 1>>
            &variable1 = instances[i]
                             .timeStepping2()
                             .instancesLocal()[j]
                             .getSlotConnectorData()
                             ->variable1;

        if (stateIndex >= variable1.size()) {
          continue;
        }
        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace, 1>>
            fieldVariableStates = variable1[stateIndex].values;

        int nValues = fiberFunctionSpace->nDofsLocalWithoutGhosts();
        double *values = valuesLocal.data() + furtherDataIndex * nValues;

        // int componentNo, int nValues, const dof_no_t *dofNosLocal, const
        // double *values
        fieldVariableStates->setValues(
            0, nValues,
            fiberFunctionSpace->meshPartition()->dofNosLocal().data(), values);

        // store in cellmlAdapter
        std::shared_ptr<
            FieldVariable::FieldVariable<FiberFunctionSpace, nStates>>
            fieldVariableStatesCellML = instances[i]
                                            .timeStepping1()
                                            .instancesLocal()[j]
                                            .getSlotConnectorData()
                                            ->variable1[stateIndex]
                                            .values;

        const int componentNo = statesForTransferIndices_[stateIndex];

        // int componentNo, int nValues, const dof_no_t *dofNosLocal, const
        // double *values
        fieldVariableStatesCellML->setValues(
            componentNo, nValues,
            fiberFunctionSpace->meshPartition()->dofNosLocal().data(), values);

        VLOG(1) << "store " << nValues << " values for additional state "
                << statesForTransferIndices_[stateIndex];
      }

      // loop over algebraics to transfer
      for (int algebraicIndex = 0;
           algebraicIndex < algebraicsForTransferIndices_.size();
           algebraicIndex++, furtherDataIndex++) {
        // store in diffusion

        // get field variable
        std::vector<::Data::ComponentOfFieldVariable<FiberFunctionSpace, 1>>
            &variable2 = instances[i]
                             .timeStepping2()
                             .instancesLocal()[j]
                             .getSlotConnectorData()
                             ->variable2;

        if (algebraicIndex >= variable2.size()) {
          continue;
        }

        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace, 1>>
            fieldVariableAlgebraics = variable2[algebraicIndex].values;

        int nValues = fiberFunctionSpace->nDofsLocalWithoutGhosts();
        double *values = valuesLocal.data() + furtherDataIndex * nValues;

        // int componentNo, int nValues, const dof_no_t *dofNosLocal, const
        // double *values
        fieldVariableAlgebraics->setValues(
            0, nValues,
            fiberFunctionSpace->meshPartition()->dofNosLocal().data(), values);

        // store in CellmlAdapter
        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace, 1>>
            fieldVariableAlgebraicsCellML = instances[i]
                                                .timeStepping1()
                                                .instancesLocal()[j]
                                                .getSlotConnectorData()
                                                ->variable2[algebraicIndex]
                                                .values;

        // const int componentNo =
        // algebraicsForTransferIndices_[algebraicIndex];

        // int componentNo, int nValues, const dof_no_t *dofNosLocal, const
        // double *values
        fieldVariableAlgebraicsCellML->setValues(
            0, nValues,
            fiberFunctionSpace->meshPartition()->dofNosLocal().data(), values);

        LOG(DEBUG) << "store " << nValues << " values for algebraic "
                   << algebraicsForTransferIndices_[algebraicIndex];
        LOG(DEBUG) << *fieldVariableAlgebraics;
      }

      // store the information about whether the point is constant or
      // not_constant or neighbour_not_constant
      if (setComputeStateInformation_) {
        // get field variable
        std::vector<::Data::ComponentOfFieldVariable<FiberFunctionSpace, 1>>
            &variable2 = instances[i]
                             .timeStepping2()
                             .instancesLocal()[j]
                             .getSlotConnectorData()
                             ->variable2;

        int algebraicIndex = algebraicsForTransferIndices_.size();
        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace, 1>>
            fieldVariableAlgebraics = variable2[algebraicIndex].values;

        int nValues = fiberFunctionSpace->nDofsLocalWithoutGhosts();
        double *values = valuesLocal.data() + furtherDataIndex * nValues;

        // int componentNo, int nValues, const dof_no_t *dofNosLocal, const
        // double *values
        fieldVariableAlgebraics->setValues(
            0, nValues,
            fiberFunctionSpace->meshPartition()->dofNosLocal().data(), values);
      }

      // increase index for fiberData_ struct
      if (computingRank == rankSubset->ownRankNo())
        fiberDataNo++;
    }
  }
}

template <int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<
    nStates, nAlgebraics,
    DiffusionTimeSteppingScheme>::restoreFiberDataCheckpoint() {
  fiberPointBuffers_ = fiberPointBuffersLastCheckpoint_;
  updateFiberData();
}

template <int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<
    nStates, nAlgebraics,
    DiffusionTimeSteppingScheme>::saveFiberDataCheckpoint() {
  fiberPointBuffersLastCheckpoint_ = fiberPointBuffers_;
}