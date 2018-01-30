#ifndef FREDKINKOEHLERAPP_H
#define FREDKINKOEHLERAPP_H

#include "MooseApp.h"

class FredkinKoehlerApp;

template<>
InputParameters validParams<FredkinKoehlerApp>();

class FredkinKoehlerApp : public MooseApp
{
public:
  FredkinKoehlerApp(InputParameters parameters);
  virtual ~FredkinKoehlerApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* FREDKINKOEHLERAPP_H */
