#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_for_hyperelasticity.h"

template <typename DisplacementsFunctionSpaceType,
          typename PressureFunctionSpaceType, typename Term,
          int nDisplacementComponents>
PartitionedPetscMatForHyperelasticityBase<DisplacementsFunctionSpaceType,
                                          PressureFunctionSpaceType, Term,
                                          nDisplacementComponents>::
    PartitionedPetscMatForHyperelasticityBase(
        std::shared_ptr<PartitionedPetscVecForHyperelasticity<
            DisplacementsFunctionSpaceType, PressureFunctionSpaceType, Term,
            nDisplacementComponents>>
            partitionedPetscVecForHyperelasticity,
        int nNonZerosDiagonal, int nNonZerosOffdiagonal, std::string name)
    : PartitionedPetscMatOneComponent<FunctionSpace::Generic>(

          // create generic function space with number of entries as in the
          // given vector
          DihuContext::meshManager()
              ->createGenericFunctionSpace(
                  partitionedPetscVecForHyperelasticity->nEntriesLocal(),
                  partitionedPetscVecForHyperelasticity->meshPartition(),
                  std::string("genericMeshForMatrix") + name)
              ->meshPartition(),

          nNonZerosDiagonal, nNonZerosOffdiagonal, name),
      partitionedPetscVecForHyperelasticity_(
          partitionedPetscVecForHyperelasticity) {

  /*
   PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>>
   meshPartition, int nNonZerosDiagonal, int nNonZerosOffdiagonal, std::string
   name);
   */
}

template <typename DisplacementsFunctionSpaceType,
          typename PressureFunctionSpaceType, typename Term,
          int nDisplacementComponents>
void PartitionedPetscMatForHyperelasticityBase<
    DisplacementsFunctionSpaceType, PressureFunctionSpaceType, Term,
    nDisplacementComponents>::setValue(int componentNoRow, PetscInt row,
                                       int componentNoColumn, PetscInt column,
                                       PetscScalar value, InsertMode mode) {
  if (VLOG_IS_ON(2)) {
    std::stringstream stream;
    stream << "\"" << this->name_ << "\" setValue "
           << (mode == INSERT_VALUES ? "(insert)" : "(add)") << ", row " << row
           << ", column " << column << ", value " << value;
    VLOG(2) << stream.str();
  }

  // replace dirichlet BC values with the prescribed values
  if (partitionedPetscVecForHyperelasticity_->isPrescribed(componentNoRow,
                                                           row) ||
      partitionedPetscVecForHyperelasticity_->isPrescribed(componentNoColumn,
                                                           column)) {
    return;
  }

  if (Term::isIncompressible) {
    assert(componentNoRow < nDisplacementComponents + 1);
    assert(componentNoColumn < nDisplacementComponents + 1); //  < 4
  } else {
    assert(componentNoRow < nDisplacementComponents);
    assert(componentNoColumn < nDisplacementComponents); //  < 3
  }

  // assert(row < this->meshPartitionRows_->nDofsLocalWithGhosts());
  // assert(column < this->meshPartitionColumns_->nDofsLocalWithGhosts());

  // determine new indices
  row = partitionedPetscVecForHyperelasticity_->nonBCDofNoGlobal(componentNoRow,
                                                                 row);
  column = partitionedPetscVecForHyperelasticity_->nonBCDofNoGlobal(
      componentNoColumn, column);

  // this wraps the standard PETSc MatSetValue on the global matrix
  PetscErrorCode ierr;
  ierr = MatSetValues(this->globalMatrix_, 1, &row, 1, &column, &value, mode);
  CHKERRV(ierr);
}

template <typename DisplacementsFunctionSpaceType,
          typename PressureFunctionSpaceType, typename Term,
          int nDisplacementComponents>
void PartitionedPetscMatForHyperelasticityBase<
    DisplacementsFunctionSpaceType, PressureFunctionSpaceType, Term,
    nDisplacementComponents>::setValue(int componentNoRow, Vc::int_v row,
                                       int componentNoColumn, Vc::int_v column,
                                       PetscScalar value, InsertMode mode) {
  // loop over rows
  for (int vcComponentNo = 0; vcComponentNo < Vc::double_v::size();
       vcComponentNo++) {
    if (row[vcComponentNo] == -1 || column[vcComponentNo] == -1)
      break;

    // call the normal, scalar setValue
    this->setValue(componentNoRow, int(row[vcComponentNo]), componentNoColumn,
                   int(column[vcComponentNo]), value, mode);
  }
}

template <typename DisplacementsFunctionSpaceType,
          typename PressureFunctionSpaceType, typename Term,
          int nDisplacementComponents>
void PartitionedPetscMatForHyperelasticityBase<
    DisplacementsFunctionSpaceType, PressureFunctionSpaceType, Term,
    nDisplacementComponents>::setValue(int componentNoRow, Vc::int_v row,
                                       int componentNoColumn, Vc::int_v column,
                                       Vc::double_v value, InsertMode mode) {
  // loop over rows
  for (int vcComponentNo = 0; vcComponentNo < Vc::double_v::size();
       vcComponentNo++) {
    if (row[vcComponentNo] == -1 || column[vcComponentNo] == -1)
      break;

    // call the normal, scalar setValue
    this->setValue(componentNoRow, (int)(row[vcComponentNo]), componentNoColumn,
                   (int)(column[vcComponentNo]), (double)(value[vcComponentNo]),
                   mode);
  }
}

template <typename PressureFunctionSpaceType, typename Term,
          int nDisplacementComponents>
void PartitionedPetscMatForHyperelasticity<
    FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,
                                 BasisFunction::LagrangeOfOrder<2>>,
    PressureFunctionSpaceType, Term,
    nDisplacementComponents>::dumpMatrixGlobalNatural(std::string filename) {
  VLOG(1) << "dumpMatrixGlobalNatural, name: " << this->name_;

  // the local matrices have ordering
  // r0         r1         r2
  // ux uy uz p ux uy uz p ...

  // the final global matrix has ordering:
  // ux_0 ux_1 ux_2 ... uy_0 uy_1 uy_2 ... ... p_0 p_1 ...
  // \-----v-----/
  //    global natural

  std::shared_ptr<Partition::MeshPartition<DisplacementsFunctionSpaceType>>
      meshPartition =
          this->partitionedPetscVecForHyperelasticity_->meshPartition();

  std::shared_ptr<Partition::MeshPartition<PressureFunctionSpaceType>>
      meshPartitionPressure;

  if (Term::isIncompressible) {
    meshPartitionPressure =
        this->partitionedPetscVecForHyperelasticity_->meshPartitionPressure();
  }

  Mat matrix = this->globalMatrix_;

  std::stringstream result;
  int ownRankNo = meshPartition->ownRankNo();
  int nRanks = meshPartition->nRanks();

  PetscInt nRowsGlobal = 0;
  PetscInt nColumnsGlobal = 0;
  PetscErrorCode ierr;

  global_no_t nDisplacementDofsGlobal = meshPartition->nDofsGlobal();
  dof_no_t nDisplacementDofsLocal = meshPartition->nDofsLocalWithoutGhosts();

  global_no_t nPressureDofsGlobal = 0;
  dof_no_t nPressureDofsLocal = 0;

  if (Term::isIncompressible) {
    nPressureDofsGlobal = meshPartitionPressure->nDofsGlobal();
    nPressureDofsLocal = meshPartitionPressure->nDofsLocalWithoutGhosts();
  }

  // get global matrix sizes
  nRowsGlobal =
      nDisplacementComponents * nDisplacementDofsGlobal + nPressureDofsGlobal;
  nColumnsGlobal = nRowsGlobal;

  PetscInt nRowsMatrixNonBc = 0;
  PetscInt nColumnsMatrixNonBc = 0;
  ierr = MatGetSize(matrix, &nRowsMatrixNonBc, &nColumnsMatrixNonBc);
  CHKERRV(ierr);

  VLOG(1) << "nRowsGlobal: " << nRowsGlobal
          << ", nRows actual matrix: " << nRowsMatrixNonBc;
  VLOG(1) << "nColumnsGlobal: " << nColumnsGlobal
          << ", nColumns actual matrix: " << nColumnsMatrixNonBc;

  // get displacement value columns
  const PetscInt *ownershipRanges;
  ierr = MatGetOwnershipRanges(matrix, &ownershipRanges);
  CHKERRV(ierr);

  PetscInt ownershipBegin = 0;
  PetscInt ownershipEnd = 0;
  ierr = MatGetOwnershipRange(matrix, &ownershipBegin, &ownershipEnd);
  CHKERRV(ierr);

  PetscInt nRowsMatrixNonBcLocal = ownershipEnd - ownershipBegin;
  VLOG(1) << "nRowsMatrixNonBcLocal: " << nRowsMatrixNonBcLocal
          << ", nDisplacementDofsLocal: " << nDisplacementDofsLocal << "/"
          << nDisplacementDofsGlobal
          << ", nPressureDofsLocal: " << nPressureDofsLocal << "/"
          << nPressureDofsGlobal;

  // get global matrix values

  // get values in row-major format with global indexing
  std::vector<double> matrixValuesLocal(
      nRowsMatrixNonBcLocal * nColumnsMatrixNonBc, 0.0);

  std::vector<PetscInt> rowIndices(nRowsMatrixNonBcLocal);
  std::vector<PetscInt> columnIndices(nColumnsMatrixNonBc);

  std::iota(rowIndices.begin(), rowIndices.end(), ownershipBegin);
  std::iota(columnIndices.begin(), columnIndices.end(), 0);

  ierr = MatGetValues(matrix, nRowsMatrixNonBcLocal, rowIndices.data(),
                      nColumnsMatrixNonBc, columnIndices.data(),
                      matrixValuesLocal.data());

  std::vector<double> matrixGlobalNonBc;
  if (ownRankNo == 0)
    matrixGlobalNonBc.resize(nRowsMatrixNonBc * nColumnsMatrixNonBc);

  std::vector<int> sizesOnRanks(nRanks);
  std::vector<int> offsets(nRanks);

  for (int rankNo = 0; rankNo < nRanks; rankNo++) {
    int nRowsOnRank = ownershipRanges[rankNo + 1] - ownershipRanges[rankNo];
    sizesOnRanks[rankNo] = nRowsOnRank * nColumnsMatrixNonBc;

    // setup offsets for MPI_Gatherv
    if (rankNo == 0) {
      offsets[rankNo] = 0;
    } else {
      offsets[rankNo] = offsets[rankNo - 1] + sizesOnRanks[rankNo - 1];
    }
  }

  MPI_Gatherv(matrixValuesLocal.data(),
              nRowsMatrixNonBcLocal * nColumnsMatrixNonBc, MPI_DOUBLE,
              matrixGlobalNonBc.data(), sizesOnRanks.data(), offsets.data(),
              MPI_DOUBLE, 0, meshPartition->mpiCommunicator());

  // gather number of local dofs
  const int nComponents =
      nDisplacementComponents + (Term::isIncompressible ? 1 : 0);

  std::vector<PetscInt> nDofsLocalRanks(
      nRanks * 2); // nDisplacementDofs, nPressureDofs for every rank

  std::vector<PetscInt> nDofsLocal{nDisplacementDofsLocal, nPressureDofsLocal};

  MPI_Gather(nDofsLocal.data(), 2, MPI_INT, nDofsLocalRanks.data(), 2, MPI_INT,
             0, meshPartition->mpiCommunicator());

  VLOG(1) << "nDofsLocal: " << nDofsLocal
          << ", nDofsLocalRanks: " << nDofsLocalRanks;

  // gather dofNoLocalToDofNoNonBcGlobal_
  std::array<std::vector<std::vector<PetscInt>>, nComponents>
      dofNoLocalToDofNoNonBcGlobalRanks;
  std::array<std::vector<PetscInt>, nComponents>
      dofNoLocalToDofNoNonBcGlobalRanksBuffer;

  // and gather localDofNo to globalNaturalDofNo mappings
  std::array<std::vector<std::vector<global_no_t>>, nComponents>
      dofNosGlobalNaturalRanks; // dofNosGlobalNaturalRanks[componentNo][rankNo][dofNoLocal];
  std::array<std::vector<global_no_t>, nComponents>
      dofNosGlobalNaturalRanksBuffer;

  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    // setup offsets for MPI_Gatherv
    std::vector<int> offsets(nRanks);
    std::vector<int> sizes(nRanks);

    PetscInt totalSize = 0;
    if (ownRankNo == 0) {
      offsets[0] = 0;

      if (componentNo < nDisplacementComponents) {
        sizes[0] = nDisplacementDofsLocal;
      } else {
        sizes[0] = nPressureDofsLocal;
      }
      totalSize = sizes[0];

      for (int rankNo = 1; rankNo < nRanks; rankNo++) {
        if (componentNo < nDisplacementComponents) {
          sizes[rankNo] = nDofsLocalRanks[2 * rankNo + 0];
        } else {
          sizes[rankNo] = nDofsLocalRanks[2 * rankNo + 1];
        }

        offsets[rankNo] = offsets[rankNo - 1] + sizes[rankNo - 1];
        totalSize += sizes[rankNo];
      }
    } else {
      if (componentNo < nDisplacementComponents) {
        sizes[ownRankNo] = nDisplacementDofsLocal;
      } else {
        sizes[ownRankNo] = nPressureDofsLocal;
      }
    }

    VLOG(1) << "component " << componentNo << ", sizes: " << sizes
            << ", offsets: " << offsets << ", totalSize: " << totalSize;

    // gather dofNoLocalToDofNoNonBcGlobal_
    assert(this->partitionedPetscVecForHyperelasticity_
               ->dofNoLocalToDofNoNonBcGlobal()[componentNo]
               .size() >= sizes[ownRankNo]); // dofNoLocalToDofNoNonBcGlobal_
                                             // may contain ghost values

    dofNoLocalToDofNoNonBcGlobalRanksBuffer[componentNo].resize(totalSize);

    VLOG(1) << "MPI_Gatherv, send " << sizes[ownRankNo]
            << " ints, sizes: " << sizes << ", offsets: " << offsets
            << ", totalSize of recv buffer: " << totalSize;
    MPI_Gatherv(this->partitionedPetscVecForHyperelasticity_
                    ->dofNoLocalToDofNoNonBcGlobal()[componentNo]
                    .data(),
                sizes[ownRankNo], MPI_INT,
                dofNoLocalToDofNoNonBcGlobalRanksBuffer[componentNo].data(),
                sizes.data(), offsets.data(), MPI_INT, 0,
                meshPartition->mpiCommunicator());

    std::vector<global_no_t> dofNosGlobalNatural;

    if (componentNo < nDisplacementComponents) {
      meshPartition->getDofNosGlobalNatural(dofNosGlobalNatural);
    } else {
      meshPartitionPressure->getDofNosGlobalNatural(dofNosGlobalNatural);
    }
    assert(dofNosGlobalNatural.size() == sizes[ownRankNo]);

    dofNosGlobalNaturalRanksBuffer[componentNo].resize(totalSize);

    MPI_Gatherv(dofNosGlobalNatural.data(), sizes[ownRankNo],
                MPI_UNSIGNED_LONG_LONG,
                dofNosGlobalNaturalRanksBuffer[componentNo].data(),
                sizes.data(), offsets.data(), MPI_UNSIGNED_LONG_LONG, 0,
                meshPartition->mpiCommunicator());

    // copy the data, split by rankNo
    if (ownRankNo == 0) {
      dofNosGlobalNaturalRanks[componentNo].resize(nRanks);
      dofNoLocalToDofNoNonBcGlobalRanks[componentNo].resize(nRanks);

      for (int rankNo = 0; rankNo < nRanks; rankNo++) {
        dofNosGlobalNaturalRanks[componentNo][rankNo].resize(sizes[rankNo]);
        dofNoLocalToDofNoNonBcGlobalRanks[componentNo][rankNo].resize(
            sizes[rankNo]);

        for (dof_no_t dofNoLocal = 0; dofNoLocal < sizes[rankNo];
             dofNoLocal++) {
          dofNosGlobalNaturalRanks[componentNo][rankNo][dofNoLocal] =
              dofNosGlobalNaturalRanksBuffer[componentNo]
                                            [offsets[rankNo] + dofNoLocal];
          dofNoLocalToDofNoNonBcGlobalRanks[componentNo][rankNo][dofNoLocal] =
              dofNoLocalToDofNoNonBcGlobalRanksBuffer[componentNo]
                                                     [offsets[rankNo] +
                                                      dofNoLocal];
        }
      }
    }
  }

  if (ownRankNo == 0) {
    if (VLOG_IS_ON(1)) {
      VLOG(1) << "dofNosGlobalNaturalRanks: " << dofNosGlobalNaturalRanks;
      VLOG(1) << "dofNoLocalToDofNoNonBcGlobalRanks: "
              << dofNoLocalToDofNoNonBcGlobalRanks;

      if (nRowsMatrixNonBc < 200) {
        std::stringstream s;
        std::stringstream s2;
        for (int j = 0; j < nRowsMatrixNonBc; j++) {
          for (int i = 0; i < nColumnsMatrixNonBc; i++) {
            s << matrixGlobalNonBc[j * nColumnsMatrixNonBc + i] << ", ";
            if (fabs(matrixGlobalNonBc[j * nColumnsMatrixNonBc + i]) < 1e-7)
              s2 << ".";
            else if (fabs(matrixGlobalNonBc[j * nColumnsMatrixNonBc + i] -
                          1.0) < 1e-7)
              s2 << "1";
            else if (fabs(matrixGlobalNonBc[j * nColumnsMatrixNonBc + i] -
                          matrixGlobalNonBc[i * nColumnsMatrixNonBc + j]) <
                     2e-4)
              s2 << "s";
            else
              s2 << "x";
          }
          s << std::endl;
          s2 << std::endl;
        }
        VLOG(1) << "matrixGlobalNonBc: (" << nRowsMatrixNonBc << "x"
                << nColumnsMatrixNonBc << "):" << std::endl
                << s.str();
        // LOG(DEBUG) << "matrixGlobalNonBc: (" << nRowsMatrixNonBc << "x" <<
        // nColumnsMatrixNonBc << ") nonzeros:" << std::endl << s2.str();
      }
    }

    // copy matrixGlobalNonBc values to right place
    // old values in matrixGlobalNonBc are in dofNoNonBcGlobal indexing, new
    // values in matrixGlobalNatural will be in global natural ordering
    std::vector<double> matrixGlobalNatural(nRowsGlobal * nColumnsGlobal, 0.0);

    // set to identity
    for (int i = 0; i < nRowsGlobal; i++)
      matrixGlobalNatural[i * nColumnsGlobal + i] = 1.0;

    // copy entries of old matrix to new one
    int rowNo = 0;
    for (int rankNoRow = 0; rankNoRow < nRanks; rankNoRow++) {
      for (int componentNoRow = 0; componentNoRow < nComponents;
           componentNoRow++) {
        int nDofsLocalRankRow = nDofsLocalRanks[2 * rankNoRow + 0];
        if (componentNoRow == nDisplacementComponents) // pressure value
        {
          nDofsLocalRankRow = nDofsLocalRanks[2 * rankNoRow + 1];
        }

        for (int dofNoLocalRow = 0; dofNoLocalRow < nDofsLocalRankRow;
             dofNoLocalRow++, rowNo++) {
          int dofNoNonBcGlobalRow =
              dofNoLocalToDofNoNonBcGlobalRanks[componentNoRow][rankNoRow]
                                               [dofNoLocalRow];

          VLOG(2) << "row [rank,comp,dof] = [" << rankNoRow << ","
                  << componentNoRow << "," << dofNoLocalRow
                  << "], nonBcGlobal: " << dofNoNonBcGlobalRow;

          // if dof is dirichlet bc, do not set any values (row is already set
          // to identity matrix row)
          if (dofNoNonBcGlobalRow == -1)
            continue;

          int dofNoGlobalNaturalRow =
              componentNoRow * nDisplacementDofsGlobal +
              dofNosGlobalNaturalRanks[componentNoRow][rankNoRow]
                                      [dofNoLocalRow];

          VLOG(2) << "    globalNatural: " << dofNoGlobalNaturalRow;

          int columnNo = 0;
          for (int rankNoColumn = 0; rankNoColumn < nRanks; rankNoColumn++) {
            for (int componentNoColumn = 0; componentNoColumn < nComponents;
                 componentNoColumn++) {
              int nDofsLocalRankColumn = nDofsLocalRanks[2 * rankNoColumn + 0];
              if (componentNoColumn ==
                  nDisplacementComponents) // pressure value
              {
                nDofsLocalRankColumn = nDofsLocalRanks[2 * rankNoColumn + 1];
              }

              for (int dofNoLocalColumn = 0;
                   dofNoLocalColumn < nDofsLocalRankColumn;
                   dofNoLocalColumn++, columnNo++) {
                int dofNoNonBcGlobalColumn =
                    dofNoLocalToDofNoNonBcGlobalRanks[componentNoColumn]
                                                     [rankNoColumn]
                                                     [dofNoLocalColumn];

                VLOG(2) << "column [rank,comp,dof] = [" << rankNoColumn << ","
                        << componentNoColumn << "," << dofNoLocalColumn
                        << "], nonBcGlobal: " << dofNoNonBcGlobalColumn;

                if (dofNoNonBcGlobalColumn == -1)
                  continue;

                int dofNoGlobalNaturalColumn =
                    componentNoColumn * nDisplacementDofsGlobal +
                    dofNosGlobalNaturalRanks[componentNoColumn][rankNoColumn]
                                            [dofNoLocalColumn];

                VLOG(2) << "    globalNatural: " << dofNoGlobalNaturalColumn;
                VLOG(2) << "    value[" << dofNoNonBcGlobalRow << ","
                        << dofNoNonBcGlobalColumn << "]="
                        << matrixGlobalNonBc[dofNoNonBcGlobalRow *
                                                 nColumnsMatrixNonBc +
                                             dofNoNonBcGlobalColumn]
                        << " at matrix[" << dofNoGlobalNaturalRow << ","
                        << dofNoGlobalNaturalColumn << "]";

                matrixGlobalNatural[dofNoGlobalNaturalRow * nColumnsGlobal +
                                    dofNoGlobalNaturalColumn] =
                    matrixGlobalNonBc[dofNoNonBcGlobalRow *
                                          nColumnsMatrixNonBc +
                                      dofNoNonBcGlobalColumn];
              }
            }
          }
        }
      }
    }

    if (nRowsGlobal < 200 && VLOG_IS_ON(1)) {
      std::stringstream s2;
      for (int j = 0; j < nRowsGlobal; j++) {
        for (int i = 0; i < nColumnsGlobal; i++) {
          if (fabs(matrixGlobalNatural[j * nColumnsGlobal + i]) < 1e-7)
            s2 << ".";
          else if (fabs(matrixGlobalNatural[j * nColumnsGlobal + i] - 1.0) <
                   1e-7)
            s2 << "1";
          else if (fabs(matrixGlobalNatural[j * nColumnsGlobal + i] -
                        matrixGlobalNatural[i * nColumnsGlobal + j]) < 2e-4)
            s2 << "s";
          else
            s2 << "x";
        }
        s2 << std::endl;
      }

      VLOG(1) << "final matrix " << nRowsGlobal << "x" << nColumnsGlobal
              << " entries, " << nDisplacementDofsGlobal << " dofs_u, "
              << nPressureDofsGlobal << " dofs_p, nonzeros: " << std::endl
              << s2.str();
    }

    // write file
    std::ofstream file;
    std::stringstream matrixName;

    if (filename.find("/") != std::string::npos) {
      matrixName << filename.substr(filename.rfind("/") + 1);
    } else {
      matrixName << filename;
    }
    // matrixName << "r" << nRanks;

    filename += std::string(".m");
    OutputWriter::Generic::openFile(file, filename);

    // write header
    file << "% " << nRowsGlobal << "x" << nColumnsGlobal << " entries, "
         << nDisplacementDofsGlobal << " dofs_u, " << nPressureDofsGlobal
         << " dofs_p, " << meshPartition->nRanks() << " MPI ranks" << std::endl
         << matrixName.str() << " = ..." << std::endl
         << "[";

    for (int j = 0; j < nRowsGlobal; j++) {
      if (j != 0)
        file << "; ..." << std::endl << " ";

      for (int i = 0; i < nColumnsGlobal; i++) {
        if (i != 0)
          file << ", ";

        file << matrixGlobalNatural[j * nColumnsGlobal + i];
      }
    }
    file << "];" << std::endl;

    file.close();

    LOG(INFO) << "Matrix written to \"" << filename << "\".";
  }
}

//! get a submatrix of the upper left part (only displacements)
template <typename DisplacementsFunctionSpaceType,
          typename PressureFunctionSpaceType, typename Term,
          int nDisplacementComponents>
Mat PartitionedPetscMatForHyperelasticityBase<
    DisplacementsFunctionSpaceType, PressureFunctionSpaceType, Term,
    nDisplacementComponents>::getSubmatrix(int rowVariableNo,
                                           int columnVariableNo) {
  // assert that rowVariableNo and columnVariableNo have valid values
  if (nDisplacementComponents ==
      3) // for static problem, there is a submatrix for u (and potentially p)
  {
    if (Term::isIncompressible) {
      assert(rowVariableNo >= 0 && rowVariableNo < 2);
      assert(columnVariableNo >= 0 && columnVariableNo < 2);
    } else {
      assert(rowVariableNo == 0);
      assert(columnVariableNo == 0);
    }
  } else if (nDisplacementComponents ==
             6) // for dynamic problem, there are submatrices for u,v (and
                // potentially p)
  {
    if (Term::isIncompressible) {
      assert(rowVariableNo >= 0 && rowVariableNo < 3);
      assert(columnVariableNo >= 0 && columnVariableNo < 3);
    } else {
      assert(rowVariableNo >= 0 && rowVariableNo < 2);
      assert(columnVariableNo >= 0 && columnVariableNo < 2);
    }
  }

  MPI_Comm mpiCommunicator =
      this->partitionedPetscVecForHyperelasticity_->meshPartition()
          ->mpiCommunicator();

  IS indexSetRows = nullptr;
  IS indexSetColumns = nullptr;

  // no velocity components
  if (nDisplacementComponents == 3) {
    if (rowVariableNo == 0) {
      indexSetRows = this->partitionedPetscVecForHyperelasticity_
                         ->displacementDofsGlobal();
    } else if (rowVariableNo == 1) {
      indexSetRows =
          this->partitionedPetscVecForHyperelasticity_->pressureDofsGlobal();
    }

    if (columnVariableNo == 0) {
      indexSetColumns = this->partitionedPetscVecForHyperelasticity_
                            ->displacementDofsGlobal();
    } else if (columnVariableNo == 1) {
      indexSetColumns =
          this->partitionedPetscVecForHyperelasticity_->pressureDofsGlobal();
    }
  } else // with velocity components
  {
    if (rowVariableNo == 0) {
      indexSetRows = this->partitionedPetscVecForHyperelasticity_
                         ->displacementDofsGlobal();
    } else if (rowVariableNo == 1) {
      indexSetRows =
          this->partitionedPetscVecForHyperelasticity_->velocityDofsGlobal();
    } else if (rowVariableNo == 2) {
      indexSetRows =
          this->partitionedPetscVecForHyperelasticity_->pressureDofsGlobal();
    }

    if (columnVariableNo == 0) {
      indexSetColumns = this->partitionedPetscVecForHyperelasticity_
                            ->displacementDofsGlobal();
    } else if (columnVariableNo == 1) {
      indexSetColumns =
          this->partitionedPetscVecForHyperelasticity_->velocityDofsGlobal();
    } else if (columnVariableNo == 2) {
      indexSetColumns =
          this->partitionedPetscVecForHyperelasticity_->pressureDofsGlobal();
    }
  }

  assert(indexSetRows);
  assert(indexSetColumns);

  Mat submatrix;
  PetscErrorCode ierr;

  ierr = MatCreateSubMatrix(this->globalMatrix_, indexSetRows, indexSetColumns,
                            MAT_INITIAL_MATRIX, &submatrix);
  CHKERRABORT(mpiCommunicator, ierr);

  return submatrix;
}
