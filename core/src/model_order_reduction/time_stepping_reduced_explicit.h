#pragma once

#include "control/dihu_context.h"
#include "function_space/function_space.h"
#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"

namespace ModelOrderReduction {
template <typename TimeSteppingExplicitType>
class TimeSteppingSchemeOdeReducedExplicit
    : public TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType> {
public:
  //! constructor
  TimeSteppingSchemeOdeReducedExplicit(DihuContext context, std::string name);

  //! destructor
  virtual ~TimeSteppingSchemeOdeReducedExplicit(){};

  //! initialize timestepping member
  void initialize();

  //! evaluates the right hand side function
  void evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output,
                                                 int timeStepNo,
                                                 double currentTime);

protected:
};

} // namespace ModelOrderReduction

#include "model_order_reduction/time_stepping_reduced_explicit.tpp"
