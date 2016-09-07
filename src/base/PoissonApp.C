#include "PoissonApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

// Poisson Includes
#include "BoundaryIntegralFMM.h"

template<>
InputParameters validParams<PoissonApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

PoissonApp::PoissonApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  PoissonApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  PoissonApp::associateSyntax(_syntax, _action_factory);
}

PoissonApp::~PoissonApp()
{
}

// External entry point for dynamic application loading
extern "C" void PoissonApp__registerApps() { PoissonApp::registerApps(); }
void
PoissonApp::registerApps()
{
  registerApp(PoissonApp);
}

// External entry point for dynamic object registration
extern "C" void PoissonApp__registerObjects(Factory & factory) { PoissonApp::registerObjects(factory); }
void
PoissonApp::registerObjects(Factory & factory)
{
}

// External entry point for dynamic syntax association
extern "C" void PoissonApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { PoissonApp::associateSyntax(syntax, action_factory); }
void
PoissonApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
