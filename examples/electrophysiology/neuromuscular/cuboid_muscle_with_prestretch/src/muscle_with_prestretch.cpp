#include <iostream>
#include <cstdlib>

#include "opendihu.h"

////// Possible options

// #define Shorten
// #define HodgkinHuxley
#define HodgkinHuxlexRazumova

// #define FiberDiffusionSolver TimeSteppingScheme::ImplicitEuler
// #define FiberDiffusionSolver TimeSteppingScheme::CrankNicolson
template <typename T>
using FiberDiffusionSolver = TimeSteppingScheme::ImplicitEuler<T>;
// template<typename T> using FiberDiffusionSolver =
// TimeSteppingScheme::CrankNicolson<T>;

////// ^^^^^^^^^^^^^^^^

#ifdef Shorten
#define N_STATES 57
#define N_ALGEBRAICS 1
#endif
#ifdef HodgkinHuxley
#define N_STATES 4
#define N_ALGEBRAICS 9
#endif
#ifdef HodgkinHuxlexRazumova
#define N_STATES 9
#define N_ALGEBRAICS 19
#endif

// define material
struct Material : Equation::SolidMechanics::HyperelasticityBase {
  static constexpr bool isIncompressible =
      true; //< if the formulation is incompressible, then,
            // strainEnergyDensityFunctionVolumetric will not be considered
  static constexpr bool usesFiberDirection =
      false; //< if the decoupled form uses the 4th or 5th invariants, Ibar4,
             // Ibar2, this means it is an anisotropic material
  static constexpr bool usesActiveStress =
      false; //< if the value of an active stress term will be added to the
             // stress

  // material parameters
  static constexpr auto c1 = PARAM(0); //< material parameter
  static constexpr auto c2 = PARAM(1); //< material parameter

  static constexpr int nMaterialParameters =
      2; //< number of material parameters

  //! the isochoric part of the decoupled strain energy function,
  //! Ψ_iso(Ibar1,Ibar2,Ibar4,Ibar5), in terms of the reduced invariants
  static const auto constexpr strainEnergyDensityFunctionIsochoric =
      c1 * (Ibar1 - INT(3)) + c2 * (Ibar2 - INT(3));

  //! the volumetric part of the decoupled strain energy function, Ψ_vol(J),
  //! only used for compressible formulation (isIncompressible == false)
  static const auto constexpr strainEnergyDensityFunctionVolumetric = INT(0);

  //! coupled form of the strain energy function, Ψ(I1,I2,I3), as alternative to
  //! the two decoupled functions
  static const auto constexpr strainEnergyDensityFunctionCoupled = INT(0);

  //! another coupled form of the strain energy function, Ψ(C), dependent on
  //! right Cauchy Green tensor, C. it must only depend on variables C11, C12,
  //! C13, C22, C23, C33.
  static const auto constexpr strainEnergyDensityFunctionCoupledDependentOnC =
      INT(0);
};


// Fast monodomain: 0D / 1D
using MonodomainSolver = FastMonodomainSolver< // a wrapper that improves
                                               // performance of multidomain
    Control::MultipleInstances<                // fibers
        OperatorSplitting::Strang<
            Control::MultipleInstances<TimeSteppingScheme::Heun< // fiber
                                                                 // reaction
                                                                 // term
                CellmlAdapter<N_STATES, N_ALGEBRAICS, // depends on the cellml
                                                      // model
                              FunctionSpace::FunctionSpace<
                                  Mesh::StructuredDeformableOfDimension<1>,
                                  BasisFunction::LagrangeOfOrder<1>>>>>,
            Control::MultipleInstances<
                FiberDiffusionSolver< // fiber diffusion, note that implicit
                                      // euler gives lower error in this case
                                      // than crank nicolson
                    SpatialDiscretization::FiniteElementMethod<
                        Mesh::StructuredDeformableOfDimension<1>,
                        BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2>,
                        Equation::Dynamic::IsotropicDiffusion>>>>>>;

int main(int argc, char *argv[]) {
  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);
  Control::Coupling<
      // static trans-iso material for "prestretch"
      SpatialDiscretization::HyperelasticitySolver<Material>,
      // actual simlation
      Control::Coupling<
          // Term1
          MonodomainSolver,
          // Term2
          MuscleContractionSolver<Mesh::StructuredDeformableOfDimension<3>>>

      >
      problem(settings);

  problem.run();

  return EXIT_SUCCESS;
}
