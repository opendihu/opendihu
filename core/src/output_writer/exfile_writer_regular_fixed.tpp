#include "output_writer/exfile_writer.h"


namespace OutputWriter
{   

//! write exnode file to given stream
template<int D, typename BasisFunctionType>
void ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
outputExelem(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables)
{
  
}

//! write exnode file to given stream
template<int D, typename BasisFunctionType>
void ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
outputExnode(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables)
{
  
}

  
};  //namespace