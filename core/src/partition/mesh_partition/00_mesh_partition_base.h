#pragma once

#include <memory>
#include <petscdmda.h>

#include "control/types.h"
#include "partition/rank_subset.h"

namespace Partition {

/** base class for mesh partition */
class MeshPartitionBase {
public:
  //! constructor, store the rankSubset
  MeshPartitionBase(std::shared_ptr<RankSubset> rankSubset);

  //! virtual destructor
  virtual ~MeshPartitionBase();

  //! number of ranks
  int nRanks() const;

  //! own rank id in the current communcator
  int ownRankNo();

  //! own rank id in the current communcator
  std::shared_ptr<RankSubset> rankSubset() const;

  //! number of entries in the current partition (this usually refers to the
  //! elements)
  virtual element_no_t nElementsLocal() const = 0;

  //! number of nodes in total
  virtual global_no_t nElementsGlobal() const = 0;

  //! remove all dofs from the vector that are not handled in the local
  //! partition
  virtual void
  extractLocalDofsWithoutGhosts(std::vector<double> &values) const = 0;

  //! get the MPI communicator that is needed for the work portion
  MPI_Comm mpiCommunicator() const;

  //! fill the dofNosLocal vector
  void createLocalDofOrderings();

  //! get a vector of local dof nos, range [0,nDofsLocalWithoutGhosts] are the
  //! dofs without ghost dofs, the whole vector are the dofs with ghost dofs
  //! (only for structured mesh)e
  const std::vector<PetscInt> &dofNosLocal() const;

  //! get the number of local dofs, with ghosts
  virtual dof_no_t nDofsLocalWithGhosts() const = 0;

  //! get the number of local dofs, without ghosts
  virtual dof_no_t nDofsLocalWithoutGhosts() const = 0;

  //! number of dofs in total
  virtual global_no_t nDofsGlobal() const = 0;

  //! number of dofs in total for the numbering used with boundary conditions
  virtual global_no_t nDofsGlobalForBoundaryConditions() const;

  //! get a vector of local dof nos including ghost dofs, in the natural
  //! ordering
  virtual const std::vector<dof_no_t> &dofNosLocalNaturalOrdering() const;

  //! get the node no in global petsc ordering from a local node no
  virtual global_no_t getNodeNoGlobalPetsc(node_no_t nodeNoLocal) const = 0;

  //! get the local dof no for a global petsc dof no, does not work for ghost
  //! nodes
  virtual dof_no_t getDofNoLocal(global_no_t dofNoGlobalPetsc,
                                 bool &isLocal) const = 0;

  //! get the rank on which the global natural node is located
  virtual int
  getRankOfDofNoGlobalNatural(global_no_t dofNoGlobalNatural) const = 0;

  //! transform the global natural numbering to the local numbering
  virtual node_no_t
  getNodeNoLocalFromGlobalNatural(global_no_t nodeNoGlobalNatural,
                                  bool &isOnLocalDomain) const = 0;

  //! get the local element no from the global element no, isOnLocalDomain is
  //! true
  virtual element_no_t getElementNoLocal(global_no_t elementNoGlobalPetsc,
                                         bool &isOnLocalDomain) const = 0;

  //! get the global natural element no for a local element no
  virtual global_no_t
  getElementNoGlobalNatural(element_no_t elementNoLocal) const = 0;

  //! get a PETSc IS (index set) with the same information as dofNosLocal_
  const IS &dofNosLocalIS() const;

  //! get a PETSc IS (index set) with the same information as dofNosLocal_, but
  //! without ghost dofs
  const IS &dofNosLocalNonGhostIS() const;

protected:
  std::shared_ptr<RankSubset>
      rankSubset_; //< the set of ranks that compute something where this
                   //partition is a part of, also holds the MPI communciator

  std::vector<dof_no_t> dofNosLocal_; //< vector of all local nos of non-ghost
                                      //dofs followed by the ghost dofs
  IS dofNosLocalIS_; //< index set (IS) with the indices of the local dof nos
                     //(including ghosts)
  IS dofNosLocalNonGhostIS_; //< index set (IS) with the indices of the local
                             //dof nos (without ghosts)
};

} // namespace Partition
