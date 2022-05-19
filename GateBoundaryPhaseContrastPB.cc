


#include "GateBoundaryPhaseContrastPB.hh"

#include "GateEMStandardProcessMessenger.hh"

//-----------------------------------------------------------------------------
GateBoundaryPhaseContrastPB::GateBoundaryPhaseContrastPB():GateVProcess("BoundaryPhaseContrast")
{
  SetDefaultParticle("gamma");
  SetProcessInfo("Boundary process for X-ray");
// pMessenger = new GateEMStandardProcessMessenger(this) ;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
G4VProcess* GateBoundaryPhaseContrastPB::CreateProcess(G4ParticleDefinition *)
{
  return new G4BoundaryPhaseContrast(GetG4ProcessName());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void GateBoundaryPhaseContrastPB::ConstructProcess(G4ProcessManager * manager)
{
  manager->AddDiscreteProcess(GetProcess());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
bool GateBoundaryPhaseContrastPB::IsApplicable(G4ParticleDefinition * par)
{
  return ( par == G4Gamma::GammaDefinition() );
}
//-----------------------------------------------------------------------------

MAKE_PROCESS_AUTO_CREATOR_CC(GateBoundaryPhaseContrastPB)
