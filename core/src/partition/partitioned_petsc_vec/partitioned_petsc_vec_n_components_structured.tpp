#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"

//! constructor
template <typename MeshType, typename BasisFunctionType, int nComponents>
PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                         nComponents>::
    PartitionedPetscVecNComponentsStructured(
        std::shared_ptr<Partition::MeshPartition<
            FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>,
            MeshType>>
            meshPartition,
        std::string name)
    : PartitionedPetscVecBase<
          FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>>(
          meshPartition, name),
      nExtractedComponents_(0) {
  // dm_ = meshPartition->dmElements();

  createVector();
}

//! constructor, copy from existing vector
template <typename MeshType, typename BasisFunctionType, int nComponents>
template <int nComponents2>
PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                         nComponents>::
    PartitionedPetscVecNComponentsStructured(
        PartitionedPetscVec<
            FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>,
            nComponents2> &rhs,
        std::string name, bool reuseData, int rhsComponentNoBegin)
    : PartitionedPetscVecBase<
          FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>>(
          rhs.meshPartition(), name) {
  if (reuseData && nComponents > nComponents2) {
    reuseData = false;
    LOG(WARNING) << "Could not create PartitionedPetscVec \"" << name
                 << "\" as copy of rhs (\"" << rhs.name()
                 << "\") with reuseData=true, nComponents=" << nComponents
                 << " > " << nComponents2 << ". Setting reuseData to false.";
  }

  if (reuseData) {
    if (rhs.currentRepresentation_ ==
        Partition::values_representation_t::representationContiguous) {
      LOG(FATAL) << "Constructor of PartitionedPetscVec with reuseData=true is "
                    "not for contiguous representation of rhs! Call "
                    "rhs.setRepresentationGlobal beforehand.";
    }

    // reuse the Petsc Vec's of the rhs PartitionedPetscVec

    // loop over the components of this field variable
    for (int componentNo = 0; componentNo < std::min(nComponents, nComponents2);
         componentNo++) {
      if (rhsComponentNoBegin + componentNo >= rhs.vectorGlobal_.size()) {
        LOG(FATAL) << "Trying to create a partitioned petsc vec with "
                   << nComponents << " from another one with " << nComponents2
                   << ", starting at component " << rhsComponentNoBegin;
      }
      vectorGlobal_[componentNo] =
          rhs.vectorGlobal_[rhsComponentNoBegin + componentNo];
      vectorLocal_[componentNo] =
          rhs.vectorLocal_[rhsComponentNoBegin + componentNo];
    }

    valuesContiguous_ = rhs.valuesContiguous_;
    extractedData_ = rhs.extractedData_;
    savedVectorLocal_ = rhs.savedVectorLocal_;
    savedVectorGlobal_ = rhs.savedVectorGlobal_;

    this->currentRepresentation_ = rhs.currentRepresentation_;

    // create VecNest object, if number of components > 1
    if (nComponents > 1) {
      VLOG(2) << "\"" << this->name_ << "\" has " << nComponents
              << " components, VecCreateNest";
      PetscErrorCode ierr;
      ierr = VecCreateNest(this->meshPartition_->mpiCommunicator(), nComponents,
                           NULL, vectorGlobal_.data(), &vectorNestedGlobal_);
      CHKERRV(ierr);
    }
  } else {
    // create new Petsc Vec's
    createVector();

    LOG(DEBUG)
        << "\"" << this->name_ << "\" contruct empty vector (local size: "
        << this->meshPartition_->nNodesLocalWithoutGhosts() << ") from rhs \""
        << rhs.name() << "\", representation: "
        << Partition::valuesRepresentationString[rhs.currentRepresentation()];

    // copy the values of rhs
    setValues(rhs);
  }
}

//! create a distributed Petsc vector, according to partition
template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                              nComponents>::createVector() {
  VLOG(2) << "\"" << this->name_ << "\" createVector with " << nComponents
          << " components, size local: "
          << this->meshPartition_->nNodesLocalWithoutGhosts()
          << ", global: " << this->meshPartition_->nNodesGlobal()
          << ", ghost dof nos global/petsc: "
          << this->meshPartition_->ghostDofNosGlobalPetsc();
  PetscErrorCode ierr;

  // The local vector contains the nodal/dof values for the local portion of the
  // current rank. This includes ghost nodes. The global vector manages the
  // whole data and also only stores the local portion of it on the current
  // rank. It shares the memory with the local vector. The local vector uses
  // local indices whereas the global vector is accessed using globalPetsc
  // indexing. If we want to manipulate data, we have to ensure consistency in
  // the parallel global vector (VecGhostUpdateBegin,VecGhostUpdateEnd) and then
  // fetch the local portion of the global vector together with ghost values
  // into a local vector (VecGhostGetLocalForm), then manipulate the values
  // (using standard VecSetValues/VecGetValues routines) and commit the results
  // (VecGhostRestoreLocalForm, VecGhostUpdateBegin, VecGhostUpdateEnd).

  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    const dof_no_t nGhostDofs = this->meshPartition_->nDofsLocalWithGhosts() -
                                this->meshPartition_->nDofsLocalWithoutGhosts();
    VLOG(2) << "\"" << this->name_ << "\" component " << componentNo
            << ", VecCreateGhost " << *this->meshPartition_->rankSubset()
            << " size local: "
            << this->meshPartition_->nDofsLocalWithoutGhosts()
            << ", global: " << this->meshPartition_->nDofsGlobal()
            << ", n dofs: " << nGhostDofs;
    ierr = VecCreateGhost(this->meshPartition_->mpiCommunicator(),
                          this->meshPartition_->nDofsLocalWithoutGhosts(),
                          this->meshPartition_->nDofsGlobal(), nGhostDofs,
                          this->meshPartition_->ghostDofNosGlobalPetsc().data(),
                          &vectorGlobal_[componentNo]);
    CHKERRV(ierr);

#if 0
    // debugging tests, learn how ghost value communcation works
    if (componentNo == 0)
    {
      Vec debugLocal_;
      Vec debugGlobal_;

      int ownRankNo = this->meshPartition_->rankSubset()->ownRankNo();
      std::vector<int> ghostDofNosGlobalDebug({ownRankNo+1});
      std::vector<int> ghostDofNosLocalDebug({1});
      std::vector<int> dofNosLocalDebug({0,1});
      std::vector<int> dofNosGlobalDebug({ownRankNo,ownRankNo+1});

      int nGhostDofsDebug = 1;
      int nValuesLocalDebug = 1;
      int nValuesGlobalDebug = 9;

      if (ownRankNo == 7)
      {
        nGhostDofsDebug = 0;
        nValuesLocalDebug = 2;
      }

      // step 1: construct geometry field
      // createVector()
      ierr = VecCreateGhost(this->meshPartition_->mpiCommunicator(), nValuesLocalDebug,
                            nValuesGlobalDebug, nGhostDofsDebug, ghostDofNosGlobalDebug.data(), &debugGlobal_); CHKERRV(ierr);

      //ierr = PetscObjectSetName((PetscObject) debugGlobal_, "debug"); CHKERRV(ierr);

      //ierr = VecGhostUpdateBegin(debugGlobal_, INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);
      //ierr = VecGhostUpdateEnd(debugGlobal_, INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);

      //ierr = VecSetFromOptions(debugGlobal_); CHKERRV(ierr);
      ierr = VecGhostGetLocalForm(debugGlobal_, &debugLocal_); CHKERRV(ierr);
      // now representation is local
      // --

      // step 2: this->geometryField_->setValuesWithoutGhosts(geometryValues)
      // field_variable: this->values_->setValues()
      // setValues()
      std::vector<int> dofNoDebug({0,1});
      std::vector<double> valuesDebug({(double)10+ownRankNo, (double)10+ownRankNo+1});
      valuesDebug.resize(nValuesLocalDebug);
      ierr = VecSetValues(debugLocal_, nValuesLocalDebug, dofNoDebug.data(), valuesDebug.data(), INSERT_VALUES); CHKERRV(ierr);
      // now representation is local
      // --

      // set representation global
      //ierr = VecGhostRestoreLocalForm(debugGlobal_, &debugLocal_); CHKERRV(ierr);
      // --

      // step 3: this->geometryField_->finishGhostManipulation()
      // finishGhostManipulation()
      // representation must be local
      //ierr = VecGhostRestoreLocalForm(debugGlobal_, &debugLocal_); CHKERRV(ierr);
      //ierr = VecGhostUpdateBegin(debugGlobal_, ADD_VALUES, SCATTER_REVERSE); CHKERRV(ierr);
      //ierr = VecGhostUpdateEnd(debugGlobal_, ADD_VALUES, SCATTER_REVERSE); CHKERRV(ierr);
      // now representation is global
      // --

      // now representation is global

      // startGhostManipulation
      //ierr = VecGhostUpdateBegin(debugGlobal_, INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);
      //ierr = VecGhostUpdateEnd(debugGlobal_, INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);
      //ierr = VecGhostGetLocalForm(debugGlobal_, &debugLocal_); CHKERRV(ierr);
      // --

      // now representation is local
      // zeroGhostBuffer
      //std::vector<double> zeroValues(nGhostDofsDebug, 0.0);
      //ierr = VecSetValues(debugLocal_, nGhostDofsDebug, ghostDofNosLocalDebug.data(), zeroValues.data(), INSERT_VALUES); CHKERRV(ierr);
      //ierr = VecSetValues(debugGlobal_, nGhostDofsDebug, ghostDofNosGlobalDebug.data(), zeroValues.data(), INSERT_VALUES); CHKERRV(ierr);
      // --

      // now representation is global

      // get global values

      // set representation global
      //ierr = VecGhostRestoreLocalForm(debugGlobal_, &debugLocal_); CHKERRV(ierr);
      // --

      // set representation local
      //ierr = VecGhostGetLocalForm(debugGlobal_, &debugLocal_); CHKERRV(ierr);
      // --

      ierr = VecGhostUpdateBegin(debugGlobal_, INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);
      ierr = VecGhostUpdateEnd(debugGlobal_, INSERT_VALUES, SCATTER_FORWARD); CHKERRV(ierr);


      // get values from global vector
      std::vector<double> globalValues(nValuesLocalDebug+nGhostDofsDebug);
      ierr = VecGetValues(debugGlobal_, nValuesLocalDebug+nGhostDofsDebug, dofNosGlobalDebug.data(), globalValues.data()); CHKERRV(ierr);
      LOG(INFO) << "local + ghost values (from global vector): " << globalValues;

      // get values from local vector
      ierr = VecGetValues(debugLocal_, nValuesLocalDebug, dofNoDebug.data(), valuesDebug.data()); CHKERRV(ierr);
      LOG(DEBUG) << "local values (from local vector): " << valuesDebug;

      // get ghost values (from local vector)
      std::vector<double> localGhostValues(nGhostDofsDebug);
      ierr = VecGetValues(debugLocal_, nGhostDofsDebug, ghostDofNosLocalDebug.data(), localGhostValues.data()); CHKERRV(ierr);

      LOG(DEBUG) << "ghost values (from local vector): " << localGhostValues;

      // get ghost values (from global vector)
      ierr = VecGetValues(debugGlobal_, nGhostDofsDebug, ghostDofNosGlobalDebug.data(), localGhostValues.data()); CHKERRV(ierr);

      LOG(DEBUG) << "ghost values (from global vector): " << localGhostValues;

      //ierr = VecGetValues();
    }
#endif

    // initialize PETSc vector object
    // ierr = VecCreate(this->meshPartition_->mpiCommunicator(),
    // &vectorLocal_[componentNo]); CHKERRV(ierr);
    ierr = PetscObjectSetName((PetscObject)vectorGlobal_[componentNo],
                              this->name_.c_str());
    CHKERRV(ierr);

    // set sparsity type and other options
    ierr = VecSetFromOptions(vectorGlobal_[componentNo]);
    CHKERRV(ierr);
    ierr = VecGhostGetLocalForm(vectorGlobal_[componentNo],
                                &vectorLocal_[componentNo]);
    CHKERRV(ierr);

    // ignore negative indices. This is needed when Vc::int_v contains the
    // indices
    ierr = VecSetOption(vectorLocal_[componentNo], VEC_IGNORE_NEGATIVE_INDICES,
                        PETSC_TRUE);
    CHKERRV(ierr);
  }

  // createVector acts like startGhostManipulation as it also gets the local
  // vector (VecGhostGetLocalForm) to work on.
  this->currentRepresentation_ =
      Partition::values_representation_t::representationLocal;
  nExtractedComponents_ = 0;

  // create VecNest object, if number of components > 1
  if (nComponents > 1) {
    VLOG(2) << "\"" << this->name_ << "\" has " << nComponents
            << " components, VecCreateNest";
    ierr = VecCreateNest(this->meshPartition_->mpiCommunicator(), nComponents,
                         NULL, vectorGlobal_.data(), &vectorNestedGlobal_);
    CHKERRV(ierr);
  }
  VLOG(1) << "created field variable \"" << this->name_ << "\", " << this;
}

// set the internal representation to be global, i.e. using the global vectors
template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::setRepresentationGlobal() {
  VLOG(2) << "\"" << this->name_
          << "\" setRepresentationGlobal, previous representation: "
          << this->getCurrentRepresentationString();
  // set the internal representation to be global, i.e. using the global vectors

  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationLocal) {
    // restore local vector back into global vector, does not consider ghost
    // dofs use finishGhostManipulation to handle ghost dofs, that requires
    // communication

    PetscErrorCode ierr;
    // loop over the components of this field variable
    for (int componentNo = 0; componentNo < nComponents; componentNo++) {
      ierr = VecGhostRestoreLocalForm(vectorGlobal_[componentNo],
                                      &vectorLocal_[componentNo]);
      CHKERRV(ierr);
    }
    this->currentRepresentation_ =
        Partition::values_representation_t::representationGlobal;
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationContiguous) {
    if (nComponents == 1) {
      LOG(FATAL) << "nComponents is 1 and representation is contiguous. This "
                    "should not happen. It would lead to infinite recursion.";
    }

    // this sets the representation to local
    restoreValuesContiguous();

    setRepresentationGlobal();
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationGlobal) {
    // already global, do nothing
  } else {
    LOG(FATAL) << "Cannot set vector representation from \""
               << this->getCurrentRepresentationString() << "\" to \"global\".";
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::setRepresentationLocal() {
  VLOG(2) << "\"" << this->name_
          << "\" setRepresentationLocal, previous representation: "
          << this->getCurrentRepresentationString();

  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationGlobal) {
    // retrieve the local vector from the global vector, does not consider ghost
    // dofs use startGhostManipulation to also set the correct values for the
    // ghost dofs, that requires communication

    PetscErrorCode ierr;
    // loop over the components of this field variable
    for (int componentNo = 0; componentNo < nComponents; componentNo++) {
      ierr = VecGhostGetLocalForm(vectorGlobal_[componentNo],
                                  &vectorLocal_[componentNo]);
      CHKERRV(ierr);
    }
    this->currentRepresentation_ =
        Partition::values_representation_t::representationLocal;
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationContiguous) {
    // this sets the representation to local
    restoreValuesContiguous();
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationInvalid) {
    LOG(FATAL) << "\"" << this->name_
               << "\" setRepresentationLocal, previous representation: "
               << this->getCurrentRepresentationString()
               << ". This is not directly possible, call "
                  "restoreExtractedComponent instead.";
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationLocal) {
    // already local, do nothing
  } else {
    LOG(FATAL) << "Cannot set vector representation from \""
               << this->getCurrentRepresentationString() << "\" to \"local\".";
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::startGhostManipulation() {
  VLOG(2) << "\"" << this->name_ << "\" startGhostManipulation";

  if (this->currentRepresentation_ !=
      Partition::values_representation_t::representationGlobal) {
    LOG(FATAL)
        << "\"" << this->name_
        << "\", startGhostManipulation called when representation is not "
           "global (but "
        << this->getCurrentRepresentationString()
        << "), this overwrites the previous values and fetches the last from "
           "the global vectors!"
        << std::endl
        << "Call setRepresentationGlobal() before startGhostManipulation() or "
           "check if startGhostManipulation() "
        << "is even necessary (because the representation is already local).";
  }

  // copy the global values into the local vectors, distributing ghost values
  PetscErrorCode ierr;

  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    ierr = VecGhostUpdateBegin(vectorGlobal_[componentNo], INSERT_VALUES,
                               SCATTER_FORWARD);
    CHKERRV(ierr);
    // ierr = VecZeroEntries(vectorLocal_[componentNo]); CHKERRV(ierr);
    // ierr = DMGlobalToLocalBegin(*dm_, vectorGlobal_[componentNo],
    // INSERT_VALUES, vectorLocal_[componentNo]); CHKERRV(ierr);
  }

  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    ierr = VecGhostUpdateEnd(vectorGlobal_[componentNo], INSERT_VALUES,
                             SCATTER_FORWARD);
    CHKERRV(ierr);
    // ierr = DMGlobalToLocalEnd(*dm_, vectorGlobal_[componentNo],
    // INSERT_VALUES, vectorLocal_[componentNo]); CHKERRV(ierr);
    ierr = VecGhostGetLocalForm(vectorGlobal_[componentNo],
                                &vectorLocal_[componentNo]);
    CHKERRV(ierr);
  }

  this->currentRepresentation_ =
      Partition::values_representation_t::representationLocal;
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::finishGhostManipulation() {
  VLOG(2) << "\"" << this->name_ << "\" finishGhostManipulation";

  if (this->currentRepresentation_ !=
      Partition::values_representation_t::representationLocal) {
    LOG(ERROR) << "\"" << this->name_
               << "\", finishGhostManipulation called when representation is "
                  "not local (it is "
               << this->getCurrentRepresentationString()
               << "), (probably no previous startGhostManipulation)";
  }

  // Copy the local values vectors into the global vector. ADD_VALUES means that
  // ghost values are reduced (summed up)
  PetscErrorCode ierr;
  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    ierr = VecGhostRestoreLocalForm(vectorGlobal_[componentNo],
                                    &vectorLocal_[componentNo]);
    CHKERRV(ierr);

    ierr = VecGhostUpdateBegin(vectorGlobal_[componentNo], ADD_VALUES,
                               SCATTER_REVERSE);
    CHKERRV(ierr);
  }

  // loop over the components of this field variable
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    ierr = VecGhostUpdateEnd(vectorGlobal_[componentNo], ADD_VALUES,
                             SCATTER_REVERSE);
    CHKERRV(ierr);
  }
  this->currentRepresentation_ =
      Partition::values_representation_t::representationGlobal;
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                              nComponents>::zeroGhostBuffer() {
  VLOG(2) << "\"" << this->name_ << "\" zeroGhostBuffer";

  // set local ghost values to 0
  int nValues = this->meshPartition_->nDofsLocalWithGhosts() -
                this->meshPartition_->nDofsLocalWithoutGhosts();
  const PetscInt *indices = this->meshPartition_->dofNosLocal().data() +
                            this->meshPartition_->nDofsLocalWithoutGhosts();
  std::vector<double> values(nValues, 0.0);

  PetscErrorCode ierr;
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    ierr = VecSetValues(vectorLocal_[componentNo], nValues, indices,
                        values.data(), INSERT_VALUES);
    CHKERRV(ierr);
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::getValues(int componentNo,
                                                         PetscInt ni,
                                                         const PetscInt ix[],
                                                         PetscScalar y[]) {
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationGlobal) {
    VLOG(1) << "getValues called in global vector representation, must be "
               "local, now set to local";
    setRepresentationLocal();
  }

  // depending on which data representation is active, use vectorLocal or
  // valuesContiguous
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationContiguous) {
    if (componentNo > 0) {
      // shift indices
      std::vector<PetscInt> &indices = temporaryIndicesVector_;
      if (indices.size() < ni)
        indices.resize(ni);
      for (int i = 0; i < ni; i++) {
        indices[i] =
            ix[i] +
            componentNo * this->meshPartition_->nDofsLocalWithoutGhosts();
      }

      // this wraps the standard PETSc VecGetValues on the local vector
      PetscErrorCode ierr;
      ierr = VecGetValues(valuesContiguous_, ni, indices.data(), y);
      CHKERRV(ierr);
    } else {
      // this wraps the standard PETSc VecGetValues on the local vector
      PetscErrorCode ierr;
      ierr = VecGetValues(valuesContiguous_, ni, ix, y);
      CHKERRV(ierr);
    }
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationLocal) {
    // this wraps the standard PETSc VecGetValues on the local vector
    PetscErrorCode ierr;
    ierr = VecGetValues(vectorLocal_[componentNo], ni, ix, y);
    CHKERRV(ierr);
  }

  // debugging output
  if (VLOG_IS_ON(3)) {
    std::stringstream str;
    str << "\"" << this->name_ << "\" getValues(componentNo=" << componentNo
        << ", indices=";
    for (int i = 0; i < ni; i++) {
      str << ix[i] << " ";
    }
    str << ") [representation=" << this->getCurrentRepresentationString()
        << "]: ";
    for (int i = 0; i < ni; i++) {
      if (fabs(y[i]) > 1e-10)
        str << y[i] << ", ";
      else
        str << "0, ";
    }
    str << "]";
    VLOG(3) << str.str();
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType,
    nComponents>::getValuesGlobalPetscIndexing(int componentNo, PetscInt ni,
                                               const PetscInt ix[],
                                               PetscScalar y[]) {
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationGlobal) {
    VLOG(1) << "getValuesGlobalPetscIndexing called in global vector "
               "representation, must be local, now set to local";
    setRepresentationLocal();
  }

  assert(this->currentRepresentation_ ==
         Partition::values_representation_t::representationLocal);

  // this wraps the standard PETSc VecGetValues on the global vector
  PetscErrorCode ierr;
  ierr = VecGhostRestoreLocalForm(vectorGlobal_[componentNo],
                                  &vectorLocal_[componentNo]);
  CHKERRV(ierr);
  ierr = VecGetValues(vectorGlobal_[componentNo], ni, ix, y);
  CHKERRV(ierr);
  ierr = VecGhostGetLocalForm(vectorGlobal_[componentNo],
                              &vectorLocal_[componentNo]);
  CHKERRV(ierr);

  // debugging output
  if (VLOG_IS_ON(3)) {
    std::stringstream str;
    str << "\"" << this->name_
        << "\" getValuesGlobalPetscIndexing(componentNo=" << componentNo
        << ", indices=";
    for (int i = 0; i < ni; i++) {
      str << ix[i] << " ";
    }
    str << "): ";
    for (int i = 0; i < ni; i++) {
      str << y[i] << " ";
    }
    VLOG(3) << str.str();
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::setValues(int componentNo,
                                                         PetscInt ni,
                                                         const PetscInt ix[],
                                                         const PetscScalar y[],
                                                         InsertMode iora) {
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationGlobal) {
    VLOG(1) << "setValues called in global vector representation, must be "
               "local, now set to local";
    setRepresentationLocal();
  }

  // debugging output
  if (VLOG_IS_ON(3)) {
    std::stringstream str;
    str << "\"" << this->name_ << "\" setValues(componentNo=" << componentNo
        << ", indices=";
    for (int i = 0; i < ni; i++) {
      str << ix[i] << " ";
    }
    str << ", values=";
    for (int i = 0; i < ni; i++) {
      str << y[i] << " ";
    }
    str << (iora == INSERT_VALUES
                ? "INSERT_VALUES"
                : (iora == ADD_VALUES ? "ADD_VALUES" : "unknown"));
    str << ") [representation=" << this->getCurrentRepresentationString()
        << "]: ";
    VLOG(3) << str.str();
  }

  // depending on which data representation is active, use vectorLocal or
  // valuesContiguous
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationContiguous) {
    if (componentNo > 0) {
      // shift indices
      std::vector<PetscInt> &indices = temporaryIndicesVector_;
      indices.resize(ni);
      for (int i = 0; i < ni; i++) {
        indices[i] =
            ix[i] +
            componentNo * this->meshPartition_->nDofsLocalWithoutGhosts();
      }

      // this wraps the standard PETSc VecSetValues on the local vector
      PetscErrorCode ierr;
      ierr = VecSetValues(valuesContiguous_, ni, indices.data(), y, iora);
      CHKERRV(ierr);
    } else {
      // this wraps the standard PETSc VecSetValues on the local vector
      PetscErrorCode ierr;
      ierr = VecSetValues(valuesContiguous_, ni, ix, y, iora);
      CHKERRV(ierr);
    }
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationLocal) {

#ifndef NDEBUG

    for (int i = 0; i < ni; i++) {
      // local indices can be -1, which means masked out (VecSetOption(x,
      // VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE) is set)
      if (ix[i] < -1 || ix[i] >= this->meshPartition_->nDofsLocalWithGhosts()) {
        LOG(ERROR) << "\"" << this->name_ << "\" setValues ix[" << i
                   << "]=" << ix[i] << ", nDofsLocalWithGhosts: "
                   << this->meshPartition_->nDofsLocalWithGhosts();
      }
      assert(ix[i] >= -1 &&
             ix[i] < this->meshPartition_->nDofsLocalWithGhosts());
    }
#endif

    // this wraps the standard PETSc VecSetValues on the local vector
    PetscErrorCode ierr;
    ierr = VecSetValues(vectorLocal_[componentNo], ni, ix, y, iora);
    CHKERRV(ierr);

    // Note, there is also VecSetValuesLocal which acts on the global vector but
    // the indices must be provided in the local ordering. For this to work one
    // has to provide the mapping in advance (VecSetLocalToGlobalMapping)
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::setValue(int componentNo,
                                                        PetscInt row,
                                                        PetscScalar value,
                                                        InsertMode mode) {
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationGlobal) {
    VLOG(1) << "setValue called in global vector representation, must be "
               "local, now set to local";
    setRepresentationLocal();
  }

  VLOG(3) << "\"" << this->name_ << "\" setValue(componentNo=" << componentNo
          << ", row=" << row << ", value=" << value
          << (mode == INSERT_VALUES
                  ? "INSERT_VALUES"
                  : (mode == ADD_VALUES ? "ADD_VALUES" : "unknown"));

  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationLocal) {
    PetscErrorCode ierr;
    ierr = VecSetValue(vectorLocal_[componentNo], row, value, mode);
    CHKERRV(ierr);
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationContiguous) {
    PetscErrorCode ierr;
    PetscInt index =
        row + componentNo * this->meshPartition_->nDofsLocalWithoutGhosts();
    ierr = VecSetValue(valuesContiguous_, index, value, mode);
    CHKERRV(ierr);
  }
}

//! set values from another vector
template <typename MeshType, typename BasisFunctionType, int nComponents>
template <int nComponents2>
void PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                              nComponents>::
    setValues(PartitionedPetscVec<
              FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>,
              nComponents2> &rhs) {
  VLOG(3)
      << "\"" << this->name_ << "\" setValues(rhs vector \"" << rhs.name()
      << "\"), rhs representation: "
      << Partition::valuesRepresentationString[rhs.currentRepresentation()]
      << ", own representation: "
      << Partition::valuesRepresentationString[this->currentRepresentation()];

  // copy existing values from rhs PartitionedPetscVec, depending on the rhs
  // representation
  PetscErrorCode ierr;
  if (rhs.currentRepresentation() ==
      Partition::values_representation_t::representationGlobal) {
    setRepresentationGlobal();

    // copy from global vector
    for (int componentNo = 0; componentNo < std::min(nComponents, nComponents2);
         componentNo++) {
      LOG(DEBUG) << "copy component " << componentNo << " from \"" << rhs.name()
                 << "\" to \"" << this->name() << "\".";
      ierr = VecCopy(rhs.valuesGlobal(componentNo), vectorGlobal_[componentNo]);
      CHKERRV(ierr);

      // this sometimes makes Petsc hang, sometimes even later in solve, use
      // with caution!
    }
  } else if (rhs.currentRepresentation() ==
             Partition::values_representation_t::representationLocal) {
    setRepresentationLocal();

    // copy from local vector
    for (int componentNo = 0; componentNo < std::min(nComponents, nComponents2);
         componentNo++) {
      ierr = VecCopy(rhs.valuesLocal(componentNo), vectorLocal_[componentNo]);
      CHKERRV(ierr);
    }
  } else if (rhs.currentRepresentation() ==
             Partition::values_representation_t::representationContiguous) {
    // copy from contiguous vector
    ierr = VecCopy(rhs.getValuesContiguous(), valuesContiguous_);
    CHKERRV(ierr);
    this->currentRepresentation_ =
        Partition::values_representation_t::representationContiguous;
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                              nComponents>::zeroEntries() {
  VLOG(3) << "\"" << this->name_ << "\" zeroEntries";

  PetscErrorCode ierr;
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationLocal) {
    for (int componentNo = 0; componentNo < nComponents; componentNo++) {
      ierr = VecZeroEntries(vectorLocal_[componentNo]);
      CHKERRV(ierr);
    }
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationGlobal) {
    for (int componentNo = 0; componentNo < nComponents; componentNo++) {
      ierr = VecZeroEntries(vectorGlobal_[componentNo]);
      CHKERRV(ierr);
    }
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationContiguous) {
    ierr = VecZeroEntries(valuesContiguous_);
    CHKERRV(ierr);
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
Vec &PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::valuesLocal(int componentNo) {
  assert(componentNo >= 0 && componentNo < nComponents);

  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationInvalid) {
    LOG(FATAL) << "Cannot use field variable with representation set to "
                  "invalid. This happens because extractComponentShared was "
               << "called and the field variable was used afterwards. You can "
                  "only access the field variable again after "
               << "restoreExtractedComponent has been called.";
  }

  if (this->currentRepresentation_ !=
      Partition::values_representation_t::representationLocal) {
    VLOG(1) << "valuesLocal called in not local vector representation ("
            << this->getCurrentRepresentationString()
            << "), now set to local (without considering ghost dofs, call "
               "startGhostManipulation if ghosts are needed!)";
    setRepresentationLocal();
  }

  VLOG(2) << "\"" << this->name_ << "\" valuesLocal(component=" << componentNo
          << ")";

  return vectorLocal_[componentNo];
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
Vec &PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::valuesGlobal(int componentNo) {
  assert(componentNo >= 0 && componentNo < nComponents);

  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationInvalid) {
    LOG(FATAL) << "Cannot use field variable with representation set to "
                  "invalid. This happens because extractComponentShared was "
               << "called and the field variable was used afterwards. You can "
                  "only access the field variable again after "
               << "restoreExtractedComponent has been called.";
  }

  if (this->currentRepresentation_ !=
      Partition::values_representation_t::representationGlobal) {
    VLOG(1) << "valuesGlobal called in not global vector representation ("
            << this->getCurrentRepresentationString()
            << "), now set to global (without considering ghost dofs, call "
               "finishGhostManipulation if ghosts are needed!)";
    setRepresentationGlobal();
  }

  VLOG(2) << "\"" << this->name_ << "\" valuesGlobal(component=" << componentNo
          << ")";

  return vectorGlobal_[componentNo];
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
Vec &PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                              nComponents>::valuesGlobal() {
  if (this->currentRepresentation_ !=
      Partition::values_representation_t::representationGlobal)
    setRepresentationGlobal();

  if (nComponents == 1)
    return valuesGlobal(0);

  return vectorNestedGlobal_;
}

//! fill a contiguous vector with all components after each other, "struct of
//! array"-type data layout. after manipulation of the vector has finished one
//! has to call restoreValuesContiguous
template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::setRepresentationContiguous() {
  VLOG(2) << "\"" << this->name_ << "\" setRepresentationContiguous()";

  // if the contiguous representation is already being used, do nothing, also if
  // the is only one component, do not use the contiguous representation
  if (this->currentRepresentation_ ==
          Partition::values_representation_t::representationContiguous ||
      nComponents == 1) {
    return;
  }

  // if the representation is global, set to local without considering ghosts,
  // because in contiguous values we do not have ghosts
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationGlobal) {
    setRepresentationLocal();
  }

  // here the representation is local

  PetscErrorCode ierr;

  // create contiguos vector if it does not exist yet
  if (valuesContiguous_ == PETSC_NULL) {
    // ierr = VecCreate(this->meshPartition_->mpiCommunicator(),
    // &valuesContiguous_);
    // CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);
    ierr = VecCreate(MPI_COMM_SELF, &valuesContiguous_);
    CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);
    ierr =
        PetscObjectSetName((PetscObject)valuesContiguous_, this->name_.c_str());
    CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);

    // initialize size of vector
    int nEntriesLocal =
        this->meshPartition_->nDofsLocalWithoutGhosts() * nComponents;
    // int nEntriesGlobal = this->meshPartition_->nDofsGlobal() * nComponents;
    // // this could also be set the nEntriesLocal, but the the communicator
    // would have to be different (MPI_COMM_SELF)
    int nEntriesGlobal =
        nEntriesLocal; // this could also be set the nEntriesLocal, but the the
                       // communicator would have to be different
                       // (MPI_COMM_SELF)
    ierr = VecSetSizes(valuesContiguous_, nEntriesLocal, nEntriesGlobal);
    CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);

    // set sparsity type and other options
    ierr = VecSetFromOptions(valuesContiguous_);
    CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);

    LOG(DEBUG) << "\"" << this->name_
               << "\" (structured) create valuesContiguous_, nComponents = "
               << nComponents << ", nEntriesLocal = " << nEntriesLocal
               << ", nEntriesGlobal = " << nEntriesGlobal
               << ", rank subset: " << *this->meshPartition_->rankSubset()
               << ", but using MPI_COMM_SELF because valuesContiguous_ is "
                  "completely local";
  }

  if (VLOG_IS_ON(3)) {
    std::stringstream s;
    output(s);
    VLOG(3) << "before copy: " << s.str();
  }

  // implementation with get and set values
#if 0
  int nEntries = this->meshPartition_->nDofsLocalWithoutGhosts();

  std::vector<PetscInt> getIndices(nEntries);
  std::iota(getIndices.begin(), getIndices.end(), 0);

  std::vector<PetscInt> setIndices(nEntries);

  std::vector<double> values(nEntries);

  // copy values from component vectors to contiguous vector
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    // retrive values from local vector
    ierr = VecGetValues(vectorLocal_[componentNo], nEntries, getIndices.data(), values.data()); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);

    LOG(DEBUG) << "get values " << values;

    // set values in contiguous vector
    std::iota(setIndices.begin(), setIndices.end(), this->meshPartition_->nDofsLocalWithoutGhosts() * componentNo);
    ierr = VecSetValues(valuesContiguous_, nEntries, setIndices.data(), values.data(), INSERT_VALUES); CHKERRABORT(this->meshPartition_->mpiCommunicator(),ierr);
  }

#endif
// more efficient implementation that avoids data copy
#if 1
  double *valuesDataContiguous;
  ierr = VecGetArray(valuesContiguous_, &valuesDataContiguous);
  CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);

  // copy values from component vectors to contiguous vector
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    const double *valuesDataComponent;
    ierr = VecGetArrayRead(vectorLocal_[componentNo], &valuesDataComponent);
    CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);

    VLOG(3) << "  copy "
            << this->meshPartition_->nDofsLocalWithoutGhosts() * sizeof(double)
            << " bytes to contiguous array, first two values: "
            << valuesDataComponent[0] << "," << valuesDataComponent[1];
    memcpy(valuesDataContiguous +
               componentNo * this->meshPartition_->nDofsLocalWithoutGhosts(),
           valuesDataComponent,
           this->meshPartition_->nDofsLocalWithoutGhosts() * sizeof(double));

    ierr = VecRestoreArrayRead(vectorLocal_[componentNo], &valuesDataComponent);
    CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);
  }

  ierr = VecRestoreArray(valuesContiguous_, &valuesDataContiguous);
  CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);
#endif
  this->currentRepresentation_ =
      Partition::values_representation_t::representationContiguous;

  if (VLOG_IS_ON(3)) {
    std::stringstream s;
    output(s);
    VLOG(3) << "after copy: " << s.str();
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                              nComponents>::
    setRepresentation(Partition::values_representation_t representation,
                      values_modified_t values) {
  VLOG(2) << "\"" << this->name_ << "\" setRepresentation("
          << Partition::valuesRepresentationString[representation] << ", "
          << (values == values_modified_t::values_unchanged ? "values unchanged"
                                                            : "values changed")
          << "), previous representation: "
          << this->getCurrentRepresentationString();

  if (this->currentRepresentation_ == representation) {
    return;
  }

  switch (representation) {
  case Partition::values_representation_t::representationLocal:
    if (values == values_modified_t::values_unchanged) {
      // avoid restoreValuesContiguous() in implementation above
      this->currentRepresentation_ = representation;
    } else {
      setRepresentationLocal();
    }
    break;
  case Partition::values_representation_t::representationGlobal:
    if (values == values_modified_t::values_unchanged) {
      // avoid restoreValuesContiguous() in implementation above
      this->currentRepresentation_ = representation;
    } else {
      setRepresentationGlobal();
    }
    break;
  case Partition::values_representation_t::representationContiguous:
    if (values == values_modified_t::values_unchanged) {
      this->currentRepresentation_ = representation;
    } else {
      setRepresentationContiguous();
    }
    break;
  default:
    LOG(FATAL) << "The representation cannot be set to " << representation
               << ": '" << Partition::valuesRepresentationString[representation]
               << "'";
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
bool PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType,
    nComponents>::isExtractComponentSharedPossible(int componentNo) {
  if (nComponents >= 1 && componentNo == componentNo - 1) {
    return false;
  }

  // if (nExtractedComponents_ != 0)
  //   return false;
  return true;
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                              nComponents>::
    extractComponentShared(
        int componentNo,
        std::shared_ptr<PartitionedPetscVec<
            FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>, 1>>
            extractedPartitionedPetscVec) {
  VLOG(2) << "\"" << this->name_
          << "\" extractComponentShared(componentNo=" << componentNo
          << ") nComponents = " << nComponents;

  if (this->currentRepresentation_ !=
      Partition::values_representation_t::representationContiguous) {
    VLOG(1) << "Called extractComponentShared with "
            << this->getCurrentRepresentationString()
            << " representation, representation needs to be contiguous, set to "
               "contiguous";
    setRepresentationContiguous();
  }

  if (nComponents == 1) {
    setRepresentationLocal();
    extractedPartitionedPetscVec->setRepresentationLocal();

    // store the original vectorLocal and vectorGlobal of the field variable to
    // be extracted
    savedVectorLocal_ = extractedPartitionedPetscVec->valuesLocal();
    savedVectorGlobal_ = extractedPartitionedPetscVec->valuesGlobal();

    // set the local and global vector of the foreign field variable
    extractedPartitionedPetscVec->valuesLocal() = vectorLocal_[0];
    extractedPartitionedPetscVec->valuesGlobal() = vectorGlobal_[0];
  } else {
    LOG(DEBUG) << "\"" << this->name() << "\": get array read";

    if (nExtractedComponents_ == 0) {
      // get data array from valuesContiguous;
      PetscErrorCode ierr;
      ierr = VecGetArrayRead(this->valuesContiguous_, &extractedData_);
      CHKERRV(ierr);

      int nDofsLocalWithGhosts = this->meshPartition()->nDofsLocalWithGhosts();
      int nDofsLocalWithoutGhosts =
          this->meshPartition()->nDofsLocalWithoutGhosts();
      int nGhostValuesExtractedFieldVariable =
          nDofsLocalWithGhosts - nDofsLocalWithoutGhosts;

      int nValuesFollowingExtractedComponent =
          nDofsLocalWithoutGhosts * (nComponents - 1 - componentNo);
      if (nGhostValuesExtractedFieldVariable >
          nValuesFollowingExtractedComponent) {
        LOG(ERROR)
            << "Getting array of component " << componentNo << "/"
            << nComponents << "."
            << " This may lead to usage of unallocated memory for the ghosts "
            << "in the extracted field variable. (read the comment below in "
               "the code why)";
      }

      // Save the values following the normal range in extractedData_.
      // This is because after VecPlaceArray the extractedPartitionedPetscVec
      // uses the space immediately after the local values to store its ghost
      // values. This is only possible because valuesContiguous_ is assumed to
      // be long enough after the component to be extracted (dangerous
      // heuristic). When the array is later restored back to valuesContiguous_
      // (by restoreExtractedComponent), these following values have also
      // changed.
      savedValues_.resize(nGhostValuesExtractedFieldVariable);
      std::copy(extractedData_ + nDofsLocalWithoutGhosts,
                extractedData_ + nDofsLocalWithGhosts, savedValues_.data());

      VLOG(1) << "\"" << this->name_ << "\": " << this
              << " extract values for the first time";
    }

    // set array in field variable
    extractedPartitionedPetscVec->setRepresentationGlobal();

    LOG(DEBUG) << "\"" << extractedPartitionedPetscVec->name()
               << "\": place array";
    // extractedPartitionedPetscVec->currentRepresentation_ =
    // Partition::values_representation_t::representationGlobal;
    PetscErrorCode ierr;
    ierr = VecPlaceArray(
        extractedPartitionedPetscVec->valuesGlobal(0),
        extractedData_ +
            componentNo * this->meshPartition_->nDofsLocalWithoutGhosts());
    CHKERRV(ierr);

    nExtractedComponents_++;
    VLOG(1) << "\"" << this->name_ << "\",\""
            << extractedPartitionedPetscVec->name() << "\": " << this
            << " now after get-array-read incremented, nExtractedComponents_="
            << nExtractedComponents_;
  }

  this->currentRepresentation_ =
      Partition::values_representation_t::representationInvalid;
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
template <int nComponents2>
void PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                              nComponents>::
    restoreExtractedComponent(
        std::shared_ptr<PartitionedPetscVec<
            FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>,
            nComponents2>>
            extractedPartitionedPetscVec,
        int componentNo) {
  VLOG(2) << "\"" << this->name_
          << "\" restoreExtractedComponent() nComponents = " << nComponents
          << ", current representation: "
          << this->getCurrentRepresentationString();

  if (this->currentRepresentation_ !=
      Partition::values_representation_t::representationInvalid) {
    LOG(ERROR) << "restoreExtractedComponent was called on a vector with "
                  "representation "
               << this->getCurrentRepresentationString()
               << ", should be representationInvalid. Check that "
                  "extractComponentShared was called previously.";
    return;
  }

  if (nComponents == 1) {
    // set the local and global vector of the foreign field variable
    extractedPartitionedPetscVec->valuesLocal() = savedVectorLocal_;
    extractedPartitionedPetscVec->valuesGlobal() = savedVectorGlobal_;

    this->currentRepresentation_ =
        Partition::values_representation_t::representationLocal;
  } else {
    nExtractedComponents_--;
    VLOG(1) << "\"" << this->name_ << "\",\""
            << extractedPartitionedPetscVec->name() << "\": " << this
            << " in restoreExtractedComponent(" << componentNo
            << ") after decrement, nExtractedComponents_="
            << nExtractedComponents_;

    if (nExtractedComponents_ == 0) {
      // only restore the array if no component is extracted (it could be that
      // extractComponentShared has been called multiple times on multiple
      // components) assert that the extracted data array is set
      assert(extractedData_);

      // restore the saved values
      if (!savedValues_.empty()) {
        int nDofsLocalWithoutGhosts =
            this->meshPartition()->nDofsLocalWithoutGhosts();
        std::copy(savedValues_.data(),
                  savedValues_.data() + savedValues_.size(),
                  (double *)extractedData_ + nDofsLocalWithoutGhosts);
      }

      LOG(DEBUG) << "\"" << this->name() << "\": restore array read";

      // restore the data array to the valuesContiguous Vec
      PetscErrorCode ierr;
      ierr = VecRestoreArrayRead(this->valuesContiguous_, &extractedData_);
      CHKERRV(ierr);

      this->currentRepresentation_ =
          Partition::values_representation_t::representationContiguous;
    } else {
      VLOG(1) << "\"" << this->name_ << "\",\""
              << extractedPartitionedPetscVec->name() << "\": " << this
              << " do not restore array read, because nExtractedComponents_="
              << nExtractedComponents_ << " > 0";
    }

    LOG(DEBUG) << "\"" << extractedPartitionedPetscVec->name()
               << "\": reset array";

    // restore the data array in the extracted vector
    PetscErrorCode ierr;
    ierr = VecResetArray(extractedPartitionedPetscVec->valuesGlobal(0));
    CHKERRV(ierr);
    // ierr =
    // VecResetArray(extractedPartitionedPetscVec->getValuesContiguous());
    // CHKERRV(ierr);
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                              nComponents>::
    extractComponentCopy(
        int componentNo,
        std::shared_ptr<PartitionedPetscVec<
            FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>, 1>>
            extractedPartitionedPetscVec) {
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationGlobal) {
    VLOG(1) << "Called extractComponents for global representation, "
               "representation needs to be local or contiguous, set to local";
    setRepresentationLocal();
  }

  VLOG(2) << "\"" << this->name_
          << "\" extractComponentCopy(componentNo=" << componentNo << ")";

  // prepare source vector
  Vec vectorSource;
  dof_no_t dofStart = 0;
  // if the contiguous representation is already being used, return contiguous
  // vector
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationContiguous) {
    vectorSource = valuesContiguous_;
    dofStart = componentNo * this->meshPartition_->nDofsLocalWithoutGhosts();
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationLocal) {
    vectorSource = vectorLocal_[componentNo];
  } else
    assert(false);

  const double *valuesSource;
  PetscErrorCode ierr;
  ierr = VecGetArrayRead(vectorSource, &valuesSource);
  CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);

  extractedPartitionedPetscVec->setRepresentationLocal();
  double *valuesTarget;
  ierr =
      VecGetArray(extractedPartitionedPetscVec->valuesLocal(0), &valuesTarget);
  CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);

  VLOG(1) << "  copy "
          << this->meshPartition_->nDofsLocalWithoutGhosts() * sizeof(double)
          << " bytes (\"" << this->name_ << "\" component " << componentNo
          << ") to \"" << extractedPartitionedPetscVec->name() << "\"";
  memcpy(valuesTarget, valuesSource + dofStart,
         this->meshPartition_->nDofsLocalWithoutGhosts() * sizeof(double));

  // restore memory
  ierr = VecRestoreArray(extractedPartitionedPetscVec->valuesLocal(0),
                         &valuesTarget);
  CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);
  ierr = VecRestoreArrayRead(vectorSource, &valuesSource);
  CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);
}

//! set the values for the given component from a petsc Vec
template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::setValues(int componentNo,
                                                         Vec petscVector,
                                                         std::string name) {
  VLOG(3) << "\"" << this->name_
          << "\" setValues from petscVector, representation: "
          << this->getCurrentRepresentationString() << ".";

  assert(componentNo >= 0);
  assert(componentNo < vectorGlobal_.size());

  PetscErrorCode ierr;

  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationContiguous) {
    // prepare source vector
    const double *valuesSource;
    ierr = VecGetArrayRead(petscVector, &valuesSource);
    CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);

    // prepare target vector
    // if the contiguous representation is already being used, return contiguous
    // vector
    dof_no_t dofStart =
        componentNo * this->meshPartition_->nDofsLocalWithoutGhosts();

    double *valuesTarget;
    ierr = VecGetArray(valuesContiguous_, &valuesTarget);
    CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);

    VLOG(1) << "  copy "
            << this->meshPartition_->nDofsLocalWithoutGhosts() * sizeof(double)
            << " bytes to "
            << "(\"" << this->name_ << "\" component " << componentNo
            << ") from \"" << name << "\"";
    memcpy(valuesTarget + dofStart, valuesSource,
           this->meshPartition_->nDofsLocalWithoutGhosts() * sizeof(double));

    // restore memory
    ierr = VecRestoreArrayRead(petscVector, &valuesSource);
    CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);
    ierr = VecRestoreArray(valuesContiguous_, &valuesTarget);
    CHKERRABORT(this->meshPartition_->mpiCommunicator(), ierr);
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationLocal) {
    VLOG(1) << "input vector: " << PetscUtility::getStringVector(petscVector);
    ierr = VecCopy(petscVector, vectorLocal_[componentNo]);
    CHKERRV(ierr);
    VLOG(1) << "now local vector is: "
            << PetscUtility::getStringVector(vectorLocal_[componentNo]);
  } else if (this->currentRepresentation_ ==
             Partition::values_representation_t::representationGlobal) {
    ierr = VecCopy(petscVector, vectorGlobal_[componentNo]);
    CHKERRV(ierr);
  }
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<MeshType, BasisFunctionType,
                                              nComponents>::
    setValues(int componentNo,
              std::shared_ptr<PartitionedPetscVec<
                  FunctionSpace::FunctionSpace<MeshType, BasisFunctionType>, 1>>
                  fieldVariable) {
  setValues(componentNo, fieldVariable->valuesLocal(0), fieldVariable->name());
}

//! get a vector of local dof nos (from meshPartition), without ghost dofs
template <typename MeshType, typename BasisFunctionType, int nComponents>
std::vector<PetscInt> &PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::localDofNosWithoutGhosts() {
  assert(this->meshPartition_);
  return this->meshPartition_->localDofNosWithoutGhosts();
}

template <typename MeshType, typename BasisFunctionType, int nComponents>
void PartitionedPetscVecNComponentsStructured<
    MeshType, BasisFunctionType, nComponents>::output(std::ostream &stream) {
#ifndef NDEBUG
  // this method gets all local non-ghost values and outputs them to stream,
  // only on rank 0
  PetscMPIInt ownRankNo = this->meshPartition_->ownRankNo();
  PetscMPIInt nRanks = this->meshPartition_->nRanks();

  int componentNo = 0;
  Vec vector = vectorLocal_[componentNo];
  if (this->currentRepresentation_ ==
      Partition::values_representation_t::representationContiguous)
    vector = valuesContiguous_;
  else if (this->currentRepresentation_ ==
           Partition::values_representation_t::representationGlobal)
    vector = vectorGlobal_[componentNo];

  // get global size of vector
  PetscInt nEntries, nEntriesLocal;
  PetscErrorCode ierr;
  ierr = VecGetSize(vector, &nEntries);
  CHKERRV(ierr);
  ierr = VecGetLocalSize(vector, &nEntriesLocal);
  CHKERRV(ierr);

  // ident 3 tabs + 1 space to get a better alignment with
  // slot_connector_data.tpp operator<<
  stream << "\t\t\t "
         << "vector \"" << this->name_ << "\" (" << nEntries << " global, "
         << nEntriesLocal << " local entries (per component), representation: "
         << this->getCurrentRepresentationString()
         << ", vectors local[0]: " << vectorLocal_[0]
         << ", global[0]: " << vectorGlobal_[0]
         << ", contiguous: " << valuesContiguous_ << ")";

  // loop over components
  for (int componentNo = 0; componentNo < nComponents; componentNo++) {
    Vec vector = vectorLocal_[componentNo];
    if (this->currentRepresentation_ ==
        Partition::values_representation_t::representationContiguous) {
      vector = valuesContiguous_;
    } else if (this->currentRepresentation_ ==
               Partition::values_representation_t::representationGlobal) {
      vector = vectorGlobal_[componentNo];
    }

    // get global size of vector
    PetscInt nEntries, nEntriesLocal;
    PetscErrorCode ierr;
    ierr = VecGetSize(vector, &nEntries);
    CHKERRV(ierr);
    ierr = VecGetLocalSize(vector, &nEntriesLocal);
    CHKERRV(ierr);

    /*if (this->currentRepresentation_ ==
    Partition::values_representation_t::representationGlobal)
    {
      if (nRanks > 1)
        continue;
    }*/

    // retrieve local values
    int nDofsLocal = this->meshPartition_->nDofsLocalWithoutGhosts();
    std::vector<PetscInt> &indices = temporaryIndicesVector_;
    indices.resize(nDofsLocal);

    // get the local values
    std::vector<double> localValues(nDofsLocal);

    if (this->currentRepresentation_ ==
        Partition::values_representation_t::representationLocal) {
      // local representation
      vector = vectorLocal_[componentNo];

      // prepare indices
      for (int i = 0; i < nDofsLocal; i++) {
        indices[i] = this->meshPartition_->dofNosLocal()[i];
      }

      ierr =
          VecGetValues(vector, nDofsLocal, indices.data(), localValues.data());
      CHKERRV(ierr);
    } else if (this->currentRepresentation_ ==
               Partition::values_representation_t::representationContiguous) {
      // contiguous representation
      vector = valuesContiguous_;

      // prepare indices
      for (int i = 0; i < nDofsLocal; i++) {
        indices[i] =
            this->meshPartition_->dofNosLocal()[i] + componentNo * nDofsLocal;
      }

      ierr =
          VecGetValues(vector, nDofsLocal, indices.data(), localValues.data());
      CHKERRV(ierr);
    } else if (this->currentRepresentation_ ==
               Partition::values_representation_t::representationGlobal) {
      // global representation
      vector = vectorGlobal_[componentNo];

      // prepare indices
      for (int i = 0; i < nDofsLocal; i++) {
        // transfer the local nos in global dof nos, using the PETSc
        // localToGlobal mapping for the dofs
        this->meshPartition_->getDofNoGlobalPetsc(
            this->meshPartition_->dofNosLocal(), indices);
      }

      ierr =
          VecGetValues(vector, nDofsLocal, indices.data(), localValues.data());
      CHKERRV(ierr);
    }

    std::vector<int> localSizesSend(nRanks, 0);
    localSizesSend[ownRankNo] = nDofsLocal;

    std::vector<int> localSizes(nRanks);
    // MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    // void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm
    // comm) Note that th recvcount argument at the root indicates the number of
    // items it receives from each process, not the total number of items it
    // receives.

    // spread the sizes on the processes
    MPIUtility::handleReturnValue(
        MPI_Allgather(&nDofsLocal, 1, MPI_INT, localSizes.data(), 1, MPI_INT,
                      this->meshPartition_->mpiCommunicator()),
        "MPI_Gather (5)");

    int maxLocalSize;
    MPIUtility::handleReturnValue(
        MPI_Allreduce(localSizes.data() + ownRankNo, &maxLocalSize, 1, MPI_INT,
                      MPI_MAX, this->meshPartition_->mpiCommunicator()),
        "MPI_Allreduce");

    // LOG(DEBUG) << "localSizesSend: " << localSizesSend << ", localSizes: " <<
    // localSizes << ", maxLocalSize: " << maxLocalSize << ", nRanks: " <<
    // nRanks << ", rankSubset: " << *this->meshPartition_->rankSubset();

    // gather all local values to rank 0
    std::vector<double> recvBuffer(maxLocalSize * nRanks);
    std::vector<double> sendBuffer(maxLocalSize, 0.0);
    std::copy(localValues.begin(), localValues.end(), sendBuffer.begin());

    // VLOG(1) << " sendBuffer: " << sendBuffer;

    MPIUtility::handleReturnValue(
        MPI_Gather(sendBuffer.data(), maxLocalSize, MPI_DOUBLE,
                   recvBuffer.data(), maxLocalSize, MPI_DOUBLE, 0,
                   this->meshPartition_->mpiCommunicator()),
        "MPI_Gather (6)");

    if (ownRankNo == 0) {
      // ident 4 tabs to get a better alignment with slot_connector_data.tpp
      // operator<<
      stream << std::endl
             << "\t\t\t\t\"" << this->name_ << "\" component " << componentNo
             << ": local values on ranks: [";
      std::vector<double> globalValues(this->meshPartition_->nDofsGlobal());
      for (int rankNo = 0; rankNo < nRanks; rankNo++) {
        if (rankNo != 0)
          stream << ";";

        int dofNoLocalEnd = localSizes[rankNo];
        if (!VLOG_IS_ON(1)) {
          dofNoLocalEnd = std::min((int)400, (int)localSizes[rankNo]);
        }

        for (dof_no_t dofNoLocal = 0; dofNoLocal < localSizes[rankNo];
             dofNoLocal++) {
          if (dofNoLocal == dofNoLocalEnd) {
            stream << "... (" << localSizes[rankNo]
                   << " entries, only showing the first " << dofNoLocalEnd
                   << ")";
            break;
          }

          double value = recvBuffer[rankNo * maxLocalSize + dofNoLocal];
          stream << " " << value;
        }
      }
      stream << "], ";
      // if (vector == valuesContiguous_)
      //   stream << std::endl;
    } else {
      stream << "\"" << this->name_ << "\" component " << componentNo
             << ": local values (dof global:value) [";

      // VLOG(1) << "localValues: " << localValues;
      dof_no_t dofNoLocalEnd = this->meshPartition_->nDofsLocalWithoutGhosts();
      if (!VLOG_IS_ON(1)) {
        dofNoLocalEnd = std::min((int)300, (int)dofNoLocalEnd);
      }
      for (dof_no_t dofNoLocal = 0; dofNoLocal < dofNoLocalEnd; dofNoLocal++) {
        if (dofNoLocal != 0)
          stream << "  ";

        double value = localValues[dofNoLocal];
        global_no_t dofNoGlobal =
            this->meshPartition_->getDofNoGlobalPetsc(dofNoLocal);

        stream << dofNoGlobal << ":" << value;
        if (dofNoLocal == 299 && !VLOG_IS_ON(1)) {
          stream << " (" << this->meshPartition_->nDofsLocalWithoutGhosts()
                 << " entries total, only showing the first 300 (call with "
                    "-vmodule=partitioned_petsc_vec_n*=1 to show all))";
        }
      }

      if (vector == valuesContiguous_) {
        stream << "]";
      } else {
        stream << "], ";
      }
    }

    // show ghost values
    if (vector != valuesContiguous_) {
      int nGhostValues = this->meshPartition_->nDofsLocalWithGhosts() -
                         this->meshPartition_->nDofsLocalWithoutGhosts();

      stream << nGhostValues << " local ghosts";

      // retrieving ghost values only works if the current representation is
      // local
      if (nGhostValues != 0) {
        setRepresentationLocal();

        stream << " (dof global:value): [";

        // get ghost values
        std::vector<double> localGhostValues(nGhostValues);
        indices.resize(nGhostValues);

        const std::vector<PetscInt> &ghostDofNosGlobalPetsc =
            this->meshPartition_->ghostDofNosGlobalPetsc();

        if (this->currentRepresentation_ ==
            Partition::values_representation_t::representationLocal) {
          // local representation
          vector = vectorLocal_[componentNo];

          // prepare indices
          for (int i = nDofsLocal;
               i < this->meshPartition_->nDofsLocalWithGhosts(); i++) {
            indices[i - nDofsLocal] = this->meshPartition_->dofNosLocal()[i];
          }

          ierr = VecGetValues(vector, nGhostValues, indices.data(),
                              localGhostValues.data());
          CHKERRV(ierr);
        } else if (this->currentRepresentation_ ==
                   Partition::values_representation_t::representationGlobal) {
          // global representation
          vector = vectorGlobal_[componentNo];

          ierr =
              VecGetValues(vector, nGhostValues, ghostDofNosGlobalPetsc.data(),
                           localGhostValues.data());
          CHKERRV(ierr);
        }

        int nGhostValuesShow = nGhostValues;
        if (!VLOG_IS_ON(1)) {
          nGhostValuesShow = std::min((int)100, (int)nGhostValues);
        }
        for (int i = 0; i < nGhostValuesShow; i++) {
          if (i != 0)
            stream << "  ";
          stream << ghostDofNosGlobalPetsc[i] << ":" << localGhostValues[i];
        }
        stream << "]" << std::endl;

        // also output vector using Petsc viewer (not so nice and fails after
        // 1000 files)
#if 0
        PetscViewer viewer;
        static int counter = 0;
        std::stringstream vectorOutputFilename;
        vectorOutputFilename << "vector_" << counter++ << ".txt";
        ierr = PetscViewerASCIIOpen(this->meshPartition_->mpiCommunicator(), vectorOutputFilename.str().c_str(), &viewer); CHKERRV(ierr);
        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_INDEX); CHKERRV(ierr);
        ierr = VecView(vectorGlobal_[componentNo], viewer); CHKERRV(ierr);

        if (ownRankNo == 0)
        {
          stream << "(Vector also written to \"" << vectorOutputFilename.str() << "\".)";
        }
#endif
      }
    } // componentNo
  }
#endif
}
