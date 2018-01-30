/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "CoupledDirichletBC.h"

template <>
InputParameters
validParams<CoupledDirichletBC>()
{
  InputParameters params = validParams<NodalBC>();

  params.addRequiredCoupledVar("coupled_var", "Value on the Boundary");
  return params;
}

CoupledDirichletBC::CoupledDirichletBC(const InputParameters & parameters)
  : NodalBC(parameters),

    /**
     * Get a reference to the coupled variable's values.
     */
    _coupled_var_val(coupledValue("coupled_var"))
{
}

Real
CoupledDirichletBC::computeQpResidual()
{
  return _u[_qp] - _coupled_var_val[_qp];
}
