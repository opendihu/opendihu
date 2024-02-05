#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"

#include "utility/mpi_utility.h"

// ---- incompressible case ----

template <typename DisplacementsFunctionSpaceType,
          typename PressureFunctionSpaceType, typename Term, int nComponents>
std::string PartitionedPetscVecForHyperelasticity<
    DisplacementsFunctionSpaceType, PressureFunctionSpaceType, Term,
    nComponents, std::enable_if_t<Term::isIncompressible, Term>>::
    getString(bool horizontal, std::string vectorName) const {
  // do not assemble a horizontal string for console in release mode, because
  // this is only needed for debugging output
  //#ifdef NDEBUG
  //  if (horizontal)
  //    return std::string("");
  //#endif

#ifndef NDEBUG

  std::stringstream result;
  int ownRankNo = this->meshPartition_->ownRankNo();

  // the master rank collects all data and writes the file
  if (ownRankNo == 0) {
    // handle displacement values
    std::vector<std::vector<global_no_t>> dofNosGlobalNatural(
        this->meshPartition_->nRanks());
    std::vector<std::vector<double>> values(this->meshPartition_->nRanks());
    std::vector<int> nDofsOnRank(this->meshPartition_->nRanks());
    std::vector<int> nPressureDofsOnRank(this->meshPartition_->nRanks());

    VLOG(1) << "values: " << values;

    // handle own rank
    // get dof nos in global natural ordering
    int nDofsLocalWithoutGhosts =
        this->meshPartition_->nDofsLocalWithoutGhosts();
    nDofsOnRank[0] = nDofsLocalWithoutGhosts;

    std::vector<global_no_t> dofNosGlobalNaturalOwn;
    this->meshPartition_->getDofNosGlobalNatural(dofNosGlobalNaturalOwn);

    dofNosGlobalNatural[0].resize(nDofsLocalWithoutGhosts);
    std::copy(dofNosGlobalNaturalOwn.begin(), dofNosGlobalNaturalOwn.end(),
              dofNosGlobalNatural[0].begin());

    // get displacement values
    values[0].resize(nComponents * nDofsLocalWithoutGhosts);
    for (int i = 0; i < nComponents; i++) {
      this->getValues(i, nDofsLocalWithoutGhosts,
                      this->meshPartition_->dofNosLocal().data(),
                      values[0].data() + i * nDofsLocalWithoutGhosts);
    }

    VLOG(1) << "get own displacement values: " << values;

    for (int rankNo = 0; rankNo < this->meshPartition_->nRanks(); rankNo++) {
      VLOG(1) << "rank " << rankNo;

      if (rankNo == 0)
        continue;

      // receive number of dofs on rank
      MPIUtility::handleReturnValue(
          MPI_Recv(&nDofsOnRank[rankNo], 1, MPI_INT, rankNo, 0,
                   this->meshPartition_->rankSubset()->mpiCommunicator(),
                   MPI_STATUS_IGNORE),
          "MPI_Recv");

      VLOG(1) << ", nDofsOnRank: " << nDofsOnRank[rankNo];

      VLOG(1) << "recv from " << rankNo << " " << nDofsOnRank[rankNo]
              << " dofs and displacements values";

      // receive dof nos
      dofNosGlobalNatural[rankNo].resize(nDofsOnRank[rankNo]);
      MPIUtility::handleReturnValue(
          MPI_Recv(dofNosGlobalNatural[rankNo].data(), nDofsOnRank[rankNo],
                   MPI_UNSIGNED_LONG_LONG, rankNo, 0,
                   this->meshPartition_->rankSubset()->mpiCommunicator(),
                   MPI_STATUS_IGNORE),
          "MPI_Recv");

      VLOG(1) << "received displacements dofs: " << dofNosGlobalNatural[rankNo];

      // receive values
      values[rankNo].resize(nComponents * nDofsOnRank[rankNo]);
      MPIUtility::handleReturnValue(
          MPI_Recv(values[rankNo].data(), nComponents * nDofsOnRank[rankNo],
                   MPI_DOUBLE, rankNo, 0,
                   this->meshPartition_->rankSubset()->mpiCommunicator(),
                   MPI_STATUS_IGNORE),
          "MPI_Recv");

      VLOG(1) << "received displacements values: " << values[rankNo];
    }

    VLOG(1) << "values: " << values;

    // handle pressure values
    std::vector<std::vector<global_no_t>> dofNosGlobalNaturalPressure(
        this->meshPartition_->nRanks());
    std::vector<std::vector<double>> valuesPressure(
        this->meshPartition_->nRanks());

    // handle own rank
    // get global natural dof nos for pressure
    nDofsLocalWithoutGhosts =
        this->meshPartitionPressure_->nDofsLocalWithoutGhosts();
    nPressureDofsOnRank[0] = nDofsLocalWithoutGhosts;

    VLOG(1) << "nRanks: " << this->meshPartition_->nRanks()
            << ", nDofsLocalWithoutGhosts: " << nDofsLocalWithoutGhosts;

    std::vector<global_no_t> dofNosGlobalNaturalPressureOwn;
    this->meshPartitionPressure_->getDofNosGlobalNatural(
        dofNosGlobalNaturalPressureOwn);

    dofNosGlobalNaturalPressure[0].resize(nDofsLocalWithoutGhosts);
    std::copy(dofNosGlobalNaturalPressureOwn.begin(),
              dofNosGlobalNaturalPressureOwn.end(),
              dofNosGlobalNaturalPressure[0].begin());

    // get pressure values
    valuesPressure[0].resize(nDofsLocalWithoutGhosts);
    this->getValues(componentNoPressure_, nDofsLocalWithoutGhosts,
                    this->meshPartitionPressure_->dofNosLocal().data(),
                    valuesPressure[0].data());

    // loop over other ranks
    for (int rankNo = 0; rankNo < this->meshPartitionPressure_->nRanks();
         rankNo++) {
      if (rankNo == 0)
        continue;

      // receive number of dofs on rank
      VLOG(1) << "pressure rank " << rankNo
              << " n dofs: " << nPressureDofsOnRank[rankNo];

      // receive number of dofs on rank
      MPIUtility::handleReturnValue(
          MPI_Recv(&nPressureDofsOnRank[rankNo], 1, MPI_INT, rankNo, 0,
                   this->meshPartition_->rankSubset()->mpiCommunicator(),
                   MPI_STATUS_IGNORE),
          "MPI_Recv");

      VLOG(1) << "recv from " << rankNo << " " << nPressureDofsOnRank[rankNo]
              << " dofs and pressure values";

      // receive dof nos
      dofNosGlobalNaturalPressure[rankNo].resize(nPressureDofsOnRank[rankNo]);
      MPIUtility::handleReturnValue(
          MPI_Recv(
              dofNosGlobalNaturalPressure[rankNo].data(),
              nPressureDofsOnRank[rankNo], MPI_UNSIGNED_LONG_LONG, rankNo, 0,
              this->meshPartitionPressure_->rankSubset()->mpiCommunicator(),
              MPI_STATUS_IGNORE),
          "MPI_Recv");

      VLOG(1) << "received pressureDofs: "
              << dofNosGlobalNaturalPressure[rankNo];

      // receive values
      valuesPressure[rankNo].resize(nPressureDofsOnRank[rankNo]);
      MPIUtility::handleReturnValue(
          MPI_Recv(
              valuesPressure[rankNo].data(), nPressureDofsOnRank[rankNo],
              MPI_DOUBLE, rankNo, 0,
              this->meshPartitionPressure_->rankSubset()->mpiCommunicator(),
              MPI_STATUS_IGNORE),
          "MPI_Recv");

      VLOG(1) << "received pressure values: " << valuesPressure[rankNo];
    }

    // sort displacement values according to global natural dof no
    std::vector<std::pair<global_no_t, std::array<double, nComponents>>>
        displacementEntries;
    displacementEntries.reserve(this->nEntriesGlobal_);

    for (int rankNo = 0; rankNo < this->meshPartition_->nRanks(); rankNo++) {
      assert(nDofsOnRank[rankNo] == dofNosGlobalNatural[rankNo].size());
      for (int i = 0; i < nDofsOnRank[rankNo]; i++) {
        std::pair<global_no_t, std::array<double, nComponents>>
            displacementEntry;
        displacementEntry.first = dofNosGlobalNatural[rankNo][i];

        std::array<double, nComponents> valuesDof;
        for (int componentNo = 0; componentNo < nComponents; componentNo++) {
          valuesDof[componentNo] =
              values[rankNo][componentNo * nDofsOnRank[rankNo] + i];
        }
        displacementEntry.second = valuesDof;

        displacementEntries.push_back(displacementEntry);
      }
    }

    // sort list according to dof no.s
    std::sort(displacementEntries.begin(), displacementEntries.end(),
              [&](std::pair<global_no_t, std::array<double, nComponents>> a,
                  std::pair<global_no_t, std::array<double, nComponents>> b) {
                return a.first < b.first;
              });

    // sort pressure values according to global natural dof no
    std::vector<std::pair<global_no_t, double>> pressureEntries;
    pressureEntries.reserve(this->nEntriesGlobal_);

    for (int rankNo = 0; rankNo < this->meshPartitionPressure_->nRanks();
         rankNo++) {
      int nDofsOnRank = dofNosGlobalNaturalPressure[rankNo].size();
      for (int i = 0; i < nDofsOnRank; i++) {
        pressureEntries.push_back(std::pair<global_no_t, double>(
            dofNosGlobalNaturalPressure[rankNo][i], valuesPressure[rankNo][i]));
      }
    }

    // sort list according to dof no.s
    std::sort(
        pressureEntries.begin(), pressureEntries.end(),
        [&](std::pair<global_no_t, double> a,
            std::pair<global_no_t, double> b) { return a.first < b.first; });

    if (VLOG_IS_ON(1)) {
      VLOG(1) << "dofNosGlobalNatural: " << dofNosGlobalNatural;
      VLOG(1) << "values: " << values;
      VLOG(1) << "nDofsOnRank: " << nDofsOnRank;
      VLOG(1) << "nPressureDofsOnRank: " << nPressureDofsOnRank;
      VLOG(1) << "valuesPressure: " << valuesPressure;
      VLOG(1) << "displacementEntries: " << displacementEntries;
      VLOG(1) << "pressureEntries: " << pressureEntries;
    }

    // write file
    std::string newline = "\n";
    std::string separator = ", ";
    std::array<std::string, nComponents> componentNames;

    componentNames[0] = "ux";
    componentNames[1] = "uy";
    componentNames[2] = "uz";
    if (nComponents == 6) {
      componentNames[3] = "vx";
      componentNames[4] = "vy";
      componentNames[5] = "vz";
    }

    if (horizontal) {
      newline = "";
      separator = ", ";
      result << std::endl;
    } else {
      result << vectorName << "r" << this->meshPartitionPressure_->nRanks()
             << " = [";
    }

    // loop over not-pressure components (u and possibly v)
    for (int componentNo = 0; componentNo < nComponents; componentNo++) {
      // start of component
      if (horizontal) {
        result << componentNames[componentNo] << " = ["
               << newline; // print e.g. "ux = ["
      }

      // write displacement values
      for (int i = 0; i < displacementEntries.size(); i++) {
        if (i != 0)
          result << separator;
        result << displacementEntries[i].first << ":"
               << (fabs(displacementEntries[i].second[componentNo]) < 1e-13
                       ? 0
                       : displacementEntries[i].second[componentNo]);
      }

      // end of component
      if (horizontal) {
        result << newline << "]; " << std::endl;
      } else {
        result << ", ...\n ";
      }
    }

    if (horizontal) {
      result << " p = [" << newline;
    } else {
      result << ", ...\n ";
    }

    for (int i = 0; i < pressureEntries.size(); i++) {
      if (i != 0)
        result << separator;
      result << pressureEntries[i].first << ":"
             << (fabs(pressureEntries[i].second) < 1e-13
                     ? 0
                     : pressureEntries[i].second);
    }

    if (horizontal) {
      result << newline << "];" << std::endl;
    } else {
      result << "]; " << std::endl;
    }
  } else {
    // all other ranks send the data to rank 0
    int nDofsLocalWithoutGhosts =
        this->meshPartition_->nDofsLocalWithoutGhosts();

    // send number of local dofs
    MPIUtility::handleReturnValue(
        MPI_Send(&nDofsLocalWithoutGhosts, 1, MPI_INT, 0, 0,
                 this->meshPartition_->rankSubset()->mpiCommunicator()),
        "MPI_Send");

    // send global natural dof nos for displacements
    std::vector<global_no_t> dofNosGlobalNatural;
    this->meshPartition_->getDofNosGlobalNatural(dofNosGlobalNatural);

    assert(dofNosGlobalNatural.size() == nDofsLocalWithoutGhosts);

    VLOG(1) << "send to 0 " << nDofsLocalWithoutGhosts
            << " dofs and displacements values";

    MPIUtility::handleReturnValue(
        MPI_Send(dofNosGlobalNatural.data(), nDofsLocalWithoutGhosts,
                 MPI_UNSIGNED_LONG_LONG, 0, 0,
                 this->meshPartition_->rankSubset()->mpiCommunicator()),
        "MPI_Send");

    VLOG(1) << "sent displacements dofs: " << dofNosGlobalNatural;

    // send displacement values
    std::vector<double> values(nComponents * nDofsLocalWithoutGhosts);

    for (int componentNo = 0; componentNo < nComponents; componentNo++) {
      this->getValues(componentNo, nDofsLocalWithoutGhosts,
                      this->meshPartition_->dofNosLocal().data(),
                      values.data() + componentNo * nDofsLocalWithoutGhosts);
    }

    MPIUtility::handleReturnValue(
        MPI_Send(values.data(), nComponents * nDofsLocalWithoutGhosts,
                 MPI_DOUBLE, 0, 0,
                 this->meshPartition_->rankSubset()->mpiCommunicator()),
        "MPI_Send");

    VLOG(1) << "sent displacements values: " << values;

    nDofsLocalWithoutGhosts =
        this->meshPartitionPressure_->nDofsLocalWithoutGhosts();

    // send number of local pressure dofs
    MPIUtility::handleReturnValue(
        MPI_Send(&nDofsLocalWithoutGhosts, 1, MPI_INT, 0, 0,
                 this->meshPartition_->rankSubset()->mpiCommunicator()),
        "MPI_Send");

    // send global natural dof nos for pressure
    VLOG(1) << "send to 0 " << nDofsLocalWithoutGhosts
            << " pressure dofs and values";

    std::vector<global_no_t> dofNosGlobalNaturalPressure;
    this->meshPartitionPressure_->getDofNosGlobalNatural(
        dofNosGlobalNaturalPressure);

    assert(dofNosGlobalNaturalPressure.size() == nDofsLocalWithoutGhosts);

    MPIUtility::handleReturnValue(
        MPI_Send(dofNosGlobalNaturalPressure.data(), nDofsLocalWithoutGhosts,
                 MPI_UNSIGNED_LONG_LONG, 0, 0,
                 this->meshPartitionPressure_->rankSubset()->mpiCommunicator()),
        "MPI_Send");

    VLOG(1) << "sent pressure dofs: " << dofNosGlobalNaturalPressure;

    // send pressure values
    std::vector<double> valuesPressure(nDofsLocalWithoutGhosts);

    this->getValues(componentNoPressure_, nDofsLocalWithoutGhosts,
                    this->meshPartitionPressure_->dofNosLocal().data(),
                    valuesPressure.data());

    MPIUtility::handleReturnValue(
        MPI_Send(valuesPressure.data(), nDofsLocalWithoutGhosts, MPI_DOUBLE, 0,
                 0,
                 this->meshPartitionPressure_->rankSubset()->mpiCommunicator()),
        "MPI_Send");

    VLOG(1) << "sent pressure values: " << valuesPressure;
  }

  return result.str();
#else

  // in release mode, do not gather all the data from all processes
  return std::string("");
#endif
}
