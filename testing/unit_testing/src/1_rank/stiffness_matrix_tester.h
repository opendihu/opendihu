#pragma once

#include <Python.h> // this has to be the first included header
#include <vector>
#include <map>
#include "utility/petsc_utility.h"

namespace SpatialDiscretization {

class StiffnessMatrixTester {
public:
  template <typename MeshType, typename BasisFunctionType,
            typename QuadratureType, typename EquationType>
  static void
  compareMatrix(FiniteElementMethod<MeshType, BasisFunctionType, QuadratureType,
                                    EquationType> &finiteElementMethod,
                std::vector<double> &referenceMatrix) {
    Mat &stiffnessMatrix =
        finiteElementMethod.data_.stiffnessMatrix()->valuesGlobal();
    std::vector<double> matrix;
    PetscUtility::getMatrixEntries(stiffnessMatrix, matrix);

    ASSERT_EQ(matrix.size(), referenceMatrix.size())
        << "Matrix has wrong number of entries";
    for (unsigned int i = 0; i < matrix.size(); i++) {
      double difference = fabs(matrix[i] - referenceMatrix[i]);
      EXPECT_LE(difference, 1e-14)
          << "Matrix entry no. " << i << " differs by " << difference
          << ", should be " << referenceMatrix[i] << ", but is " << matrix[i];
    }
  }

  template <typename MeshType, typename BasisFunctionType,
            typename QuadratureType, typename EquationType>
  static void
  compareRhs(FiniteElementMethod<MeshType, BasisFunctionType, QuadratureType,
                                 EquationType> &finiteElementMethod,
             std::vector<double> &referenceRhs) {
    Vec &rhsVector = finiteElementMethod.data_.rightHandSide()->valuesLocal();
    std::vector<double> rhs;
    PetscUtility::getVectorEntries(rhsVector, rhs);

    ASSERT_EQ(rhs.size(), referenceRhs.size())
        << "Rhs has wrong number of entries";
    for (unsigned int i = 0; i < rhs.size(); i++) {
      double difference = fabs(rhs[i] - referenceRhs[i]);
      EXPECT_LE(difference, 1e-14)
          << "Rhs entry no. " << i << " differs by " << difference
          << ", should be " << referenceRhs[i] << ", but is " << rhs[i];
    }
  }

  template <typename MeshType, typename BasisFunctionType,
            typename QuadratureType, typename EquationType>
  static void compareSolution(
      FiniteElementMethod<MeshType, BasisFunctionType, QuadratureType,
                          EquationType> &finiteElementMethod,
      std::vector<double> &referenceSolution, double tolerance = 2e-14) {
    ASSERT_TRUE(finiteElementMethod.data_.solution() != nullptr)
        << "Solution vector is not set in finite element method.";
    Vec &solutionVector = finiteElementMethod.data_.solution()->valuesLocal();
    std::vector<double> solution;
    PetscUtility::getVectorEntries(solutionVector, solution);

    ASSERT_EQ(solution.size(), referenceSolution.size())
        << "Solution has wrong number of entries";
    for (unsigned int i = 0; i < solution.size(); i++) {
      double difference = fabs(solution[i] - referenceSolution[i]);
      EXPECT_LE(difference, tolerance)
          << "Solution entry no. " << i << " differs by " << difference
          << ", should be " << referenceSolution[i] << ", but is "
          << solution[i] << " (tolerance: " << tolerance << ")";
    }
  }

  template <typename MeshType, typename BasisFunctionType,
            typename QuadratureType, typename EquationType>
  static void checkDirichletBCInSolution(
      FiniteElementMethod<MeshType, BasisFunctionType, QuadratureType,
                          EquationType> &finiteElementMethod,
      std::map<int, double> &dirichletBC) {
    Vec &solutionVector = finiteElementMethod.data_.solution()->valuesLocal();
    std::vector<double> solution;
    PetscUtility::getVectorEntries(solutionVector, solution);

    for (std::pair<int, double> entry : dirichletBC) {
      double difference = fabs(solution[entry.first] - entry.second);
      EXPECT_LE(difference, 1e-13)
          << "Dirichlet BC on node " << entry.first << ", value "
          << entry.second << " is not met! "
          << "Actual value: " << solution[entry.first]
          << ", Difference: " << difference;
    }
  }

  template <typename MeshType1, typename BasisFunctionType1,
            typename QuadratureType1, typename EquationType1,
            typename MeshType2, typename BasisFunctionType2,
            typename QuadratureType2, typename EquationType2>
  static void
  checkEqual(FiniteElementMethod<MeshType1, BasisFunctionType1, QuadratureType1,
                                 EquationType1> &finiteElementMethod1,
             FiniteElementMethod<MeshType2, BasisFunctionType2, QuadratureType2,
                                 EquationType2> &finiteElementMethod2) {
    // rhsVector
    Vec &rhsVector1 = finiteElementMethod1.data_.rightHandSide()->valuesLocal();
    std::vector<double> rhs1;
    PetscUtility::getVectorEntries(rhsVector1, rhs1);
    Vec &rhsVector2 = finiteElementMethod2.data_.rightHandSide()->valuesLocal();
    std::vector<double> rhs2;
    PetscUtility::getVectorEntries(rhsVector2, rhs2);

    ASSERT_EQ(rhs1.size(), rhs2.size()) << "Rhs has wrong number of entries";
    for (unsigned int i = 0; i < rhs1.size(); i++) {
      double difference = fabs(rhs1[i] - rhs2[i]);
      EXPECT_LE(difference, 1e-14)
          << "Rhs entry no. " << i << " differs by " << difference << ", "
          << rhs1[i] << " != " << rhs2[i];
    }

    // stiffness matrix
    Mat &stiffnessMatrix1 =
        finiteElementMethod1.data_.stiffnessMatrix()->valuesGlobal();
    std::vector<double> matrix1;
    PetscUtility::getMatrixEntries(stiffnessMatrix1, matrix1);
    Mat &stiffnessMatrix2 =
        finiteElementMethod2.data_.stiffnessMatrix()->valuesGlobal();
    std::vector<double> matrix2;
    PetscUtility::getMatrixEntries(stiffnessMatrix2, matrix2);

    ASSERT_EQ(matrix1.size(), matrix2.size())
        << "Matrix has wrong number of entries";
    for (unsigned int i = 0; i < matrix1.size(); i++) {
      double difference = fabs(matrix1[i] - matrix2[i]);
      EXPECT_LE(difference, 2e-14)
          << "Matrix entry no. " << i << " differs by " << difference
          << ", entry1: " << matrix1[i] << " != " << matrix2[i];
    }
  }

  template <typename T1, typename T2>
  static void testMassMatrix(T1 &finiteElementMethod1, T2 &finiteElementMethod2,
                             std::vector<double> &rhsValues) {
    // create the discretization matrix if it does not already exist
    finiteElementMethod2.setMassMatrix();
    Mat &massMatrix = finiteElementMethod2.data_.massMatrix()->valuesGlobal();

    PetscInt n, m;
    MatGetSize(massMatrix, &n, &m);
    LOG(DEBUG) << "matrix size: " << n << "x" << m
               << ", rhsValues size: " << rhsValues.size() << std::endl;
    ASSERT_EQ(n, rhsValues.size());
    Vec rhsStrong, rhsWeak;

    PetscUtility::createVector(rhsStrong, n, "rhs strong");
    PetscUtility::createVector(rhsWeak, n, "rhs weak");

    PetscUtility::setVector(rhsValues, rhsStrong);

    LOG(DEBUG) << "massMatrix: "
               << PetscUtility::getStringMatrixVector(massMatrix, rhsStrong);
    MatMult(massMatrix, rhsStrong, rhsWeak);

    LOG(DEBUG) << "massMatrix * rhsStrong = rhsWeak: "
               << PetscUtility::getStringVector(rhsWeak);

    // massMatrix * f_strong = rhs_weak
    Vec &rhs = finiteElementMethod1.data_.rightHandSide()
                   ->valuesLocal(); // rhs in weak formulation

    LOG(DEBUG) << "using stencil: " << PetscUtility::getStringVector(rhs);

    std::vector<double> rhsWeakDMatrix, rhsWeakStencil;
    PetscUtility::getVectorEntries(rhsWeak, rhsWeakDMatrix);
    PetscUtility::getVectorEntries(rhs, rhsWeakStencil);

    // compare vectors
    ASSERT_EQ(rhsWeakDMatrix.size(), rhsWeakStencil.size());

    for (unsigned int i = 0; i < rhsWeakDMatrix.size(); i++) {
      double difference = fabs(rhsWeakDMatrix[i] - rhsWeakStencil[i]);
      EXPECT_LE(difference, 1e-14)
          << "Rhs entry number " << i
          << " is different. Using massMatrix: " << rhsWeakDMatrix[i]
          << ", using stencil: " << rhsWeakStencil[i]
          << ", Difference: " << difference;
    }
  }
};

} // namespace SpatialDiscretization
