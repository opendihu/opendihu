#include <Python.h> // this has to be the first included header

#include <iostream>
#include <cstdlib>
#include <fstream>

#include "gtest/gtest.h"
#include "opendihu.h"
#include "arg.h"
#include "stiffness_matrix_tester.h"
#include "node_positions_tester.h"

namespace SpatialDiscretization {

TEST(LaplaceTest, MatrixIsCorrect3DStencils) {
  std::string pythonConfig = R"(
# Laplace 3D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [4, 4, 4],
    "physicalExtent": [4.0, 4.0, 4.0],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  FiniteElementMethod<Mesh::StructuredRegularFixedOfDimension<3>,
                      BasisFunction::LagrangeOfOrder<>, Quadrature::None,
                      Equation::Static::Laplace>
      equationDiscretized(settings);

  equationDiscretized.run();

  std::vector<double> referenceMatrix(15625, 0.0);

  const double stencil[2][2][2] = {
      {{-4. / 12, 0. / 12}, {0. / 12, 1. / 12}}, // center
      {{0. / 12, 1. / 12}, {1. / 12, 1. / 12}},  // bottom
  };

  // fill with stencil values
  auto matrixIndex = [](int i, int j) { return j * 125 + i; };

  // loop over elements
  for (int x = 0; x < 4; x++) {
    for (int y = 0; y < 4; y++) {
      for (int z = 0; z < 4; z++) {
        // nodes:
        // 6--7
        // 4--5
        //
        // 2--3
        // 0--1
        int node0 = z * 25 + y * 5 + x;
        int node1 = z * 25 + y * 5 + x + 1;
        int node2 = z * 25 + (y + 1) * 5 + x;
        int node3 = z * 25 + (y + 1) * 5 + x + 1;
        int node4 = (z + 1) * 25 + y * 5 + x;
        int node5 = (z + 1) * 25 + y * 5 + x + 1;
        int node6 = (z + 1) * 25 + (y + 1) * 5 + x;
        int node7 = (z + 1) * 25 + (y + 1) * 5 + x + 1;

        // add contribution from node 0 to all nodes          x  y  z
        referenceMatrix[matrixIndex(node0, node0)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node0, node1)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node0, node2)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node0, node3)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node0, node4)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node0, node5)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node0, node6)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node0, node7)] += stencil[1][1][1];

        // add contribution from node 1 to all nodes
        referenceMatrix[matrixIndex(node1, node0)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node1, node1)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node1, node2)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node1, node3)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node1, node4)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node1, node5)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node1, node6)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node1, node7)] += stencil[1][0][1];

        // add contribution from node 2 to all nodes
        referenceMatrix[matrixIndex(node2, node0)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node2, node1)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node2, node2)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node2, node3)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node2, node4)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node2, node5)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node2, node6)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node2, node7)] += stencil[0][1][1];

        // add contribution from node 3 to all nodes
        referenceMatrix[matrixIndex(node3, node0)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node3, node1)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node3, node2)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node3, node3)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node3, node4)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node3, node5)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node3, node6)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node3, node7)] += stencil[0][0][1];

        // add contribution from node 4 to all nodes          x  y  z
        referenceMatrix[matrixIndex(node4, node0)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node4, node1)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node4, node2)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node4, node3)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node4, node4)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node4, node5)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node4, node6)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node4, node7)] += stencil[1][1][0];

        // add contribution from node 5 to all nodes
        referenceMatrix[matrixIndex(node5, node0)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node5, node1)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node5, node2)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node5, node3)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node5, node4)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node5, node5)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node5, node6)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node5, node7)] += stencil[1][0][0];

        // add contribution from node 6 to all nodes
        referenceMatrix[matrixIndex(node6, node0)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node6, node1)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node6, node2)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node6, node3)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node6, node4)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node6, node5)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node6, node6)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node6, node7)] += stencil[0][1][0];

        // add contribution from node 7 to all nodes
        referenceMatrix[matrixIndex(node7, node0)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node7, node1)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node7, node2)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node7, node3)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node7, node4)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node7, node5)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node7, node6)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node7, node7)] += stencil[0][0][0];
      }
    }
  }

  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
}

TEST(LaplaceTest, StructuredDeformableMatrixIsCorrect3DStencils) {
  std::string pythonConfig = R"(
# Laplace 3D
config = {
  "disablePrinting": False,
  "disableMatrixPrinting": True,
  "FiniteElementMethod" : {
    "nElements": [4, 4, 4],
    "physicalExtent": [4.0, 4.0, 4.0],
    "relativeTolerance": 1e-15,
  },
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  FiniteElementMethod<Mesh::StructuredDeformableOfDimension<3>,
                      BasisFunction::LagrangeOfOrder<>, Quadrature::Gauss<2>,
                      Equation::Static::Laplace>
      equationDiscretized(settings);

  equationDiscretized.run();

  std::vector<double> referenceMatrix(15625, 0.0);

  const double stencil[2][2][2] = {
      {{-4. / 12, 0. / 12}, {0. / 12, 1. / 12}}, // center
      {{0. / 12, 1. / 12}, {1. / 12, 1. / 12}},  // bottom
  };

  // fill with stencil values
  auto matrixIndex = [](int i, int j) { return j * 125 + i; };

  // loop over elements
  for (int x = 0; x < 4; x++) {
    for (int y = 0; y < 4; y++) {
      for (int z = 0; z < 4; z++) {
        // nodes:
        // 6--7
        // 4--5
        //
        // 2--3
        // 0--1
        int node0 = z * 25 + y * 5 + x;
        int node1 = z * 25 + y * 5 + x + 1;
        int node2 = z * 25 + (y + 1) * 5 + x;
        int node3 = z * 25 + (y + 1) * 5 + x + 1;
        int node4 = (z + 1) * 25 + y * 5 + x;
        int node5 = (z + 1) * 25 + y * 5 + x + 1;
        int node6 = (z + 1) * 25 + (y + 1) * 5 + x;
        int node7 = (z + 1) * 25 + (y + 1) * 5 + x + 1;

        // add contribution from node 0 to all nodes          x  y  z
        referenceMatrix[matrixIndex(node0, node0)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node0, node1)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node0, node2)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node0, node3)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node0, node4)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node0, node5)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node0, node6)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node0, node7)] += stencil[1][1][1];

        // add contribution from node 1 to all nodes
        referenceMatrix[matrixIndex(node1, node0)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node1, node1)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node1, node2)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node1, node3)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node1, node4)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node1, node5)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node1, node6)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node1, node7)] += stencil[1][0][1];

        // add contribution from node 2 to all nodes
        referenceMatrix[matrixIndex(node2, node0)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node2, node1)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node2, node2)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node2, node3)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node2, node4)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node2, node5)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node2, node6)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node2, node7)] += stencil[0][1][1];

        // add contribution from node 3 to all nodes
        referenceMatrix[matrixIndex(node3, node0)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node3, node1)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node3, node2)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node3, node3)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node3, node4)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node3, node5)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node3, node6)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node3, node7)] += stencil[0][0][1];

        // add contribution from node 4 to all nodes          x  y  z
        referenceMatrix[matrixIndex(node4, node0)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node4, node1)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node4, node2)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node4, node3)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node4, node4)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node4, node5)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node4, node6)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node4, node7)] += stencil[1][1][0];

        // add contribution from node 5 to all nodes
        referenceMatrix[matrixIndex(node5, node0)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node5, node1)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node5, node2)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node5, node3)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node5, node4)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node5, node5)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node5, node6)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node5, node7)] += stencil[1][0][0];

        // add contribution from node 6 to all nodes
        referenceMatrix[matrixIndex(node6, node0)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node6, node1)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node6, node2)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node6, node3)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node6, node4)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node6, node5)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node6, node6)] += stencil[0][0][0];
        referenceMatrix[matrixIndex(node6, node7)] += stencil[0][1][0];

        // add contribution from node 7 to all nodes
        referenceMatrix[matrixIndex(node7, node0)] += stencil[1][1][1];
        referenceMatrix[matrixIndex(node7, node1)] += stencil[1][0][1];
        referenceMatrix[matrixIndex(node7, node2)] += stencil[0][1][1];
        referenceMatrix[matrixIndex(node7, node3)] += stencil[0][0][1];
        referenceMatrix[matrixIndex(node7, node4)] += stencil[1][1][0];
        referenceMatrix[matrixIndex(node7, node5)] += stencil[1][0][0];
        referenceMatrix[matrixIndex(node7, node6)] += stencil[0][1][0];
        referenceMatrix[matrixIndex(node7, node7)] += stencil[0][0][0];
      }
    }
  }

  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
}

TEST(LaplaceTest, SolverManagerWorks) {
  std::string pythonConfig = R"(
# Laplace 1D
n = 5
    
# boundary conditions
bc = {}
bc[0] = 1.0
bc[n] = 0.0

config = {
  "disablePrinting": False,
  "disableMatrixPrinting": False,
  "Solvers" : {
    "linearSolver": {
      "relativeTolerance": 1e-15
    },
  },
  "FiniteElementMethod" : {
    "nElements": n,
    "physicalExtent": 4.0,
    "dirichletBoundaryConditions": bc,
    "solverName": "linearSolver"
  }
}
)";

  DihuContext settings(argc, argv, pythonConfig);

  FiniteElementMethod<Mesh::StructuredRegularFixedOfDimension<1>,
                      BasisFunction::LagrangeOfOrder<>, Quadrature::None,
                      Equation::Static::Laplace>
      equationDiscretized(settings);

  equationDiscretized.run();

  std::vector<double> referenceMatrix = {
      1, 0,    0,    0,    0,    0, 0, -2.5, 1.25, 0,    0,    0,
      0, 1.25, -2.5, 1.25, 0,    0, 0, 0,    1.25, -2.5, 1.25, 0,
      0, 0,    0,    1.25, -2.5, 0, 0, 0,    0,    0,    0,    1};
  std::vector<double> referenceRhs = {1, -1.25, 0, 0, 0, 0};
  std::map<int, double> dirichletBC = {{0, 1.0}, {5, 0.0}};

  StiffnessMatrixTester::compareMatrix(equationDiscretized, referenceMatrix);
  StiffnessMatrixTester::compareRhs(equationDiscretized, referenceRhs);
  StiffnessMatrixTester::checkDirichletBCInSolution(equationDiscretized,
                                                    dirichletBC);
}

} // namespace SpatialDiscretization
