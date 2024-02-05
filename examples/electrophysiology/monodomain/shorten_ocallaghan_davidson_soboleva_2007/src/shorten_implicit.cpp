#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[]) {
  // 1D reaction-diffusion equation du/dt = c du^2/dx^2 + R(t), R is from cellml
  // file

  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  LOG(DEBUG) << std::string(80, '=');

  OperatorSplitting::Strang<
      TimeSteppingScheme::Heun<CellmlAdapter<58, 77> // nStates,nAlgebraics:
                                                     // 58,77 = Shorten, 4,9 =
                                                     // Hodgkin Huxley
                               >,
      TimeSteppingScheme::CrankNicolson<
          SpatialDiscretization::FiniteElementMethod<
              Mesh::StructuredRegularFixedOfDimension<1>,
              BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>,
              Equation::Dynamic::IsotropicDiffusion>>>
      problem(settings);
  problem.run();

  return EXIT_SUCCESS;
}
