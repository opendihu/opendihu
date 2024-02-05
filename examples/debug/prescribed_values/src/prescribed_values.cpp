#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[]) {
  // 1D diffusion equation du/dt = c du^2/dx^2 + r(x,t)

  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  PrescribedValues<
      FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<2>,
                                   BasisFunction::LagrangeOfOrder<2>>>
      problem(settings);

  problem.run();

  return EXIT_SUCCESS;
}
