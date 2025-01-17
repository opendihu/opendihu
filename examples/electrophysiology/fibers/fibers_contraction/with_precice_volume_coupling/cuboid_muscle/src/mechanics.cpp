#include <Python.h>
#include <iostream>
#include <cstdlib>
#include <iostream>
#include "easylogging++.h"
#include "opendihu.h"

int main(int argc, char *argv[]) {
  // initialize OpenDihu
  DihuContext settings(argc, argv);

  // define MuscleContractionSolver for 3D mechanics
<<<<<<< HEAD
  Control::PreciceAdapter<MuscleContractionSolver<
=======
  Control::PreciceAdapterVolumeCoupling<MuscleContractionSolver<
>>>>>>> develop
      Mesh::StructuredDeformableOfDimension<3>,
      Equation::SolidMechanics::
          TransverselyIsotropicMooneyRivlinIncompressibleActive3D>>
      problem(settings);

  problem.run();

  return EXIT_SUCCESS;
}