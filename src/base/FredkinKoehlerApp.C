#include "FredkinKoehlerApp.h"
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
InputParameters validParams<FredkinKoehlerApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

FredkinKoehlerApp::FredkinKoehlerApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  FredkinKoehlerApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  FredkinKoehlerApp::associateSyntax(_syntax, _action_factory);
}

FredkinKoehlerApp::~FredkinKoehlerApp()
{
}

// External entry point for dynamic application loading
extern "C" void FredkinKoehlerApp__registerApps() { FredkinKoehlerApp::registerApps(); }
void
FredkinKoehlerApp::registerApps()
{
  registerApp(FredkinKoehlerApp);
}

// External entry point for dynamic object registration
extern "C" void FredkinKoehlerApp__registerObjects(Factory & factory) { FredkinKoehlerApp::registerObjects(factory); }
void
FredkinKoehlerApp::registerObjects(Factory & factory)
{
  // Kernels
  registerKernel(Electrostatics);
  registerKernel(PolarElectricEStrong);

  // UserObjects
  registerUserObject(BoundaryIntegralFMM);

  // BoundaryCondition
  registerBoundaryCondition(CoupledDirichletBC);

  // Transfers
  registerTransfer(MultiAppAddTransfer);
}

// External entry point for dynamic syntax association
extern "C" void FredkinKoehlerApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { FredkinKoehlerApp::associateSyntax(syntax, action_factory); }
void
FredkinKoehlerApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
