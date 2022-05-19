//____________________G4PhaseContrast.hh_________________

#ifndef G4BoundaryPhaseContrast_h
#define G4BoundaryPhaseContrast_h 1

// _________________________includes___________________________

#include "G4Step.hh"                           
#include "G4VDiscreteProcess.hh"                              
#include "G4EmProcessSubType.hh"                              
#include "G4PhysicalConstants.hh"                             
#include "G4ParallelWorldProcess.hh"                          
#include "G4Gamma.hh"                                         
#include "G4TransportationManager.hh"                         
#include "G4GeometryTolerance.hh"                             
#include "G4SystemOfUnits.hh"                                
#include "GateConfiguration.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "G4Material.hh"

using namespace std;

#ifdef GATE_USE_XRAYLIB
#include <xraylib.h>
#endif

//____________________________class Definition__________________

class G4BoundaryPhaseContrast : public G4VDiscreteProcess {

public:

//________________________________FUNCTIONS_________________________

//1st function = dot product 
static G4double dot(G4double d[],G4double n[]);

//2nd function = normal 
static G4double diff (G4double a[],G4double x,G4double y,G4double z, G4int n);

//3rd function = find boundary
static G4double boundary(G4double a[],G4double x,G4double y ,G4double z);

//4th function = which fitted function?
static G4double fit(G4double x, G4double y , G4double z);
static void index(G4double f,G4double a[]);

//5th function = material
static G4double Mat(G4double a[],G4double i,G4double j,G4double k);

//____________________________________________________________________

    G4BoundaryPhaseContrast(const G4String &processName = "BoundaryPhaseContrast", G4ProcessType type = fElectromagnetic);
    ~G4BoundaryPhaseContrast();

private:

    G4BoundaryPhaseContrast(const G4BoundaryPhaseContrast &right);


    G4BoundaryPhaseContrast &operator=(const G4BoundaryPhaseContrast &right);


public:

    G4bool IsApplicable(const G4ParticleDefinition &aParticleType);

    void DoReflection();

    G4double GetRindex(G4Material *Material, G4double Energy);

    G4double GetMeanFreePath(const G4Track &aTrack, G4double , G4ForceCondition *condition);

    G4VParticleChange *PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);


private:

    G4double GetIncidentAngle();

private:
    G4double TotalMomentum;
    G4ThreeVector OldMomentum;
    G4ThreeVector OldPolarization;

    G4ThreeVector NewMomentum;
    G4ThreeVector NewPolarization;

    G4ThreeVector theGlobalNormal;

    G4double Rindex1;
    G4double Rindex2;

    G4double cost1, cost2, sint1, sint2;

    G4Material *Material1;
    G4Material *Material2;

    G4double kCarTolerance;

};

//_____________________________inline____________________________-

inline
G4bool G4BoundaryPhaseContrast::IsApplicable(const G4ParticleDefinition &aParticleType) {  
    return ( &aParticleType == G4Gamma::Gamma() );                                     
}                                                                                      

inline                                                                           
void G4BoundaryPhaseContrast::DoReflection() {                                     
    G4double PdotN = OldMomentum * theGlobalNormal;                              
    NewMomentum = OldMomentum - (2.*PdotN) * theGlobalNormal;                     
}  

inline

#ifdef GATE_USE_XRAYLIB
G4double G4BoundaryPhaseContrast::GetRindex(G4Material *Material, G4double Energy) {
#else
G4double G4BoundaryPhaseContrast::GetRindex(G4Material *Material, G4double) {
#endif
    G4double delta = 0.0;
    G4double Density = Material->GetDensity() / (g / cm3);

#ifdef GATE_USE_XRAYLIB
#if XRAYLIB_MAJOR > 3
    xrl_error *error = NULL;
    for (unsigned int i = 0; i < Material->GetElementVector()->size(); ++i) {
      delta += (1 - Refractive_Index_Re(Material->GetElementVector()->at(i)->GetSymbol(), Energy/(keV), 1.0,&error)) * Material->GetFractionVector()[i];
      if (error != NULL) {
	G4cerr << "error message: " << error->message << "\n";
	xrl_clear_error(&error);
      }
    }
#else
    for (unsigned int i = 0; i < Material->GetElementVector()->size(); ++i)
      delta += ( 1 - Refractive_Index_Re(Material->GetElementVector()->at(i)->GetSymbol(), Energy/(keV), 1.0)) * Material->GetFractionVector()[i]; //it was (1- Refractive_Ind....)
#endif    
#else
    G4Exception( "G4XrayBoundaryProcess::GetRindex", "GetRindex", FatalException, "Xraylib is not available\n");
#endif

    return 1 - delta * Density;
}           
                                           
#endif 



