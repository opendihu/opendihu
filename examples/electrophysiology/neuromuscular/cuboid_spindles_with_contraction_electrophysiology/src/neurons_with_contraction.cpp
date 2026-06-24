#include <Python.h>
#include <iostream>
#include <cstdlib>

#include <iostream>
#include "easylogging++.h"

#include "opendihu.h"

int main(int argc, char *argv[]) {
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  // define helper function space for various activation signals, this is
  // actually a vector space
  using HelperFunctionSpace =
      FunctionSpace::FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,
                                   BasisFunction::LagrangeOfOrder<1>>;

  // define overall structure of solvers
    Control::Coupling<

                           // muscle spindles solver
        TimeSteppingScheme::Heun<
            CellmlAdapter<4, 9, HelperFunctionSpace>
        >,

        // motoneuron with electro-mechanics
        // Control::Coupling<
        // mapping motor neuron signals + cortical input to actual inputs
        Control::MapDofs<
            HelperFunctionSpace,

            // electro-mechanics solver
            Control::Coupling<
                // electrophysiology solver (mockup)
                FastMonodomainSolver<Control::MultipleInstances< // subdomains in xy-plane
                OperatorSplitting::Strang<
                        Control::MultipleInstances< // fiber reaction term
                            TimeSteppingScheme::Heun<CellmlAdapter<
                                9, 19, // nStates, nAlgebraics
                                FunctionSpace::FunctionSpace<
                                    Mesh::StructuredDeformableOfDimension<1>,
                                    BasisFunction::LagrangeOfOrder<1>>>>>,
                        Control::MultipleInstances<            // fiber diffusion
                            TimeSteppingScheme::ImplicitEuler< // note that implicit euler
                                                                // gives lower error in
                                                                // this case than crank
                                                                // nicolson
                                SpatialDiscretization::FiniteElementMethod<
                                    Mesh::StructuredDeformableOfDimension<1>,
                                    BasisFunction::LagrangeOfOrder<1>,
                                    Quadrature::Gauss<2>,
                                    Equation::Dynamic::IsotropicDiffusion>>>>>>,
                // mechanics solver
                MuscleContractionSolver<
                    Mesh::StructuredDeformableOfDimension<3>
                    // BasisFunction::LagrangeOfOrder<2>
                >
            >
        >
        // >  
    >
  problem(settings);

  // run problem
  problem.run();

  return EXIT_SUCCESS;
}
