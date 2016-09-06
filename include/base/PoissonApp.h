#ifndef POISSONAPP_H
#define POISSONAPP_H

#include "MooseApp.h"

class PoissonApp;

template<>
InputParameters validParams<PoissonApp>();

class PoissonApp : public MooseApp
{
public:
  PoissonApp(InputParameters parameters);
  virtual ~PoissonApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* POISSONAPP_H */
