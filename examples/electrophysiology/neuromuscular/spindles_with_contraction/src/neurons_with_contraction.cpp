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
                Control::MultipleInstances<
                    PrescribedValues<FunctionSpace::FunctionSpace<
                        Mesh::StructuredDeformableOfDimension<3>,
                        BasisFunction::LagrangeOfOrder<1>
                    >>
                >,
                // mechanics solver
                MuscleContractionSolver<
                    Mesh::CompositeOfDimension<3>
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
