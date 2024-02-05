#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[]) {
  // 3D Laplace equation

  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  Postprocessing::ParallelFiberEstimation<BasisFunction::LagrangeOfOrder<1>>
      problem(settings);

  problem.run();

  return EXIT_SUCCESS;
}
