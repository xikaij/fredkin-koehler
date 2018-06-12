#ifndef KOALAAPP_H
#define KOALAAPP_H

#include "MooseApp.h"

class KoalaApp;

template<>
InputParameters validParams<KoalaApp>();

class KoalaApp : public MooseApp
{
public:
  KoalaApp(InputParameters parameters);
  virtual ~KoalaApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* FREDKINKOEHLERAPP_H */
