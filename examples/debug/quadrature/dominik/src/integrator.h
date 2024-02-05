#pragma once
#define numNodes 5
namespace Quadrature {

class Quadrature {
public:
};

/**
 * Quadrature class that can be used if no integrator is necessary, because
 * the integral is solve analytically. This is the case if Mesh::RegularFixed is
 * used. Those meshes use stencils.
 */
class None {
public:
  static constexpr int numberEvaluations() { return 0; }
};

}; // namespace Quadrature
