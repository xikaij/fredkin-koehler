#include "BoundaryIntegralFMM.h"

template<>
InputParameters validParams<BoundaryIntegralFMM>()
{
  InputParameters params = validParams<GeneralUserObject>();

  params.set<std::string>("built_by_action") = "add_user_object";

  return params;
}

BoundaryIntegralFMM::BoundaryIntegralFMM(InputParameters parameters) :
    GeneralUserObject(parameters)
{
}

void
BoundaryIntegralFMM::initialize()
{

}

void
BoundaryIntegralFMM::execute()
{

}

void
BoundaryIntegralFMM::finalize()
{

}
