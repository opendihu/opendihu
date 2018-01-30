#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 1D reaction-diffusion equation du/dt = c du^2/dx^2 + R(t), R is from cellml file
  
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  
  LOG(DEBUG)<<std::string(80, '=');
  
  OperatorSplitting::Godunov<
    TimeSteppingScheme::ExplicitEuler<
      CellmlAdapter
    >,
    TimeSteppingScheme::ExplicitEuler<
      SpatialDiscretization::FiniteElementMethod<
        Mesh::RegularFixed<1>,
        BasisFunction::Lagrange<>,
        Integrator::Gauss<2>,
        Equation::Dynamic::Diffusion
      >
    >
  >
  operatorSplitting(settings);
       
  Computation computation(settings, operatorSplitting);
  computation.run();
  
  return EXIT_SUCCESS;
}