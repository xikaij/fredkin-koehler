#include "KoalaApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

// Kernels
#include "Electrostatics.h"
#include "PolarElectricEStrong.h"

// UserObjects
#include "BoundaryIntegralFMM.h"

// BoundaryCondition
#include "CoupledDirichletBC.h"

// Transfer
#include "MultiAppAddTransfer.h"

template<>
InputParameters validParams<KoalaApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

KoalaApp::KoalaApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  KoalaApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  KoalaApp::associateSyntax(_syntax, _action_factory);
}

KoalaApp::~KoalaApp()
{
}

// External entry point for dynamic application loading
extern "C" void KoalaApp__registerApps() { KoalaApp::registerApps(); }
void
KoalaApp::registerApps()
{
  registerApp(KoalaApp);
}

// External entry point for dynamic object registration
extern "C" void KoalaApp__registerObjects(Factory & factory) { KoalaApp::registerObjects(factory); }
void
KoalaApp::registerObjects(Factory & factory)
{
  // Kernels
  registerKernel(Electrostatics);
  registerKernel(PolarElectricEStrong);

  // UserObjects
  registerUserObject(BoundaryIntegralFMM);

  // BoundaryCondition
  registerBoundaryCondition(CoupledDirichletBC);
}

// External entry point for dynamic syntax association
extern "C" void KoalaApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { KoalaApp::associateSyntax(syntax, action_factory); }
void
KoalaApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
