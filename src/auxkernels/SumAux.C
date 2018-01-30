//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SumAux.h"

template <>
InputParameters
validParams<SumAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("coupled_var", "Coupled variable");
  return params;
}

SumAux::SumAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    // We can couple in a value from one of our kernels with a call to coupledValueAux
    _coupled_val(coupledValue("coupled_var"))
{
}
 
Real
SumAux::computeValue()
{
  return _coupled_val[_qp] + _u[_qp];
}
