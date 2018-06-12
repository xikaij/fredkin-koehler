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

#ifndef MULTIAPPADDTRANSFER_H
#define MULTIAPPADDTRANSFER_H

#include "MultiAppTransfer.h"

// Forward declarations
class MultiAppAddTransfer;
class MooseVariableFEBase;

template <>
InputParameters validParams<MultiAppAddTransfer>();

/**
 * Copy the value to the target domain from the nearest node in the source domain.
 */
class MultiAppAddTransfer : public MultiAppTransfer
{
public:
  MultiAppAddTransfer(const InputParameters & parameters);

  /**
   * Performs basic error checking that the variable exists on MultiApp.
   */
  virtual void initialSetup() override;

  /**
   * Performs the transfer of a variable (Nonlinear or Auxiliary) to/from the Multiapp.
   */
  virtual void execute() override;

protected:
  /**
   * Performs the transfer of a variable between two problems.
   */
  void transfer(FEProblemBase & to_problem, FEProblemBase & from_problem);

  /**
   * Performs the transfer of values between a node or element.
   */
  void transferDofObject(libMesh::DofObject * to_object,
                         libMesh::DofObject * from_object,
                         MooseVariableFEBase & to_var,
                         MooseVariableFEBase & from_var);

  /// The name of the variable to transfer to
  const VariableName & _to_var_name;

  /// Name of variable transfering from
  const VariableName & _from_var_name;
};

#endif // MULTIAPPADDTRANSFER_H
