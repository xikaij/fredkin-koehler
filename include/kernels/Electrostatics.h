/*
 * From Ferret by J. Mangeri <mangeri@fzu.cz>
 */

#ifndef ELECTROSTATICS_H
#define ELECTROSTATICS_H

#include "Kernel.h"

class Electrostatics;

template<>
InputParameters validParams<Electrostatics>();

class Electrostatics: public Kernel
{
public:

  Electrostatics(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const Real _permittivity;
  const Real _len_scale;

};

#endif
