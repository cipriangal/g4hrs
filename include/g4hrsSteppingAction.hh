
#ifndef __REMOLLSTEPPINGACTION_HH
#define __REMOLLSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "g4hrsEventAction.hh"
#include "globals.hh"
#include "TFile.h"
#include "TH2D.h"

class g4hrsTransportFunction; 
class g4hrsEventAction;
class g4hrsTune;

class g4hrsSteppingAction : public G4UserSteppingAction
{
  public:
    g4hrsSteppingAction(g4hrsEventAction*);
~g4hrsSteppingAction();
void Write();

    virtual void UserSteppingAction(const G4Step*);

    void SetEnableKryptonite(G4bool k){ fEnableKryptonite = k; }

  private:
    G4bool drawFlag;
	G4double rad;
    G4bool fEnableKryptonite;
	g4hrsTransportFunction* fTransportFunction;

  g4hrsEventAction *fEvAct;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };

	int nelements;

	G4double fSeptumAngle;		// septum angle obtained from messenger (constant)
	G4double fHRSAngle;		// HRS angle obtained from messenger (constant)
	G4double septum_angle;		// local septum angle (changes sign depending on L/R HRS)
	G4double hrs_angle;		// local HRS angle (changes sign depending on L/R HRS)
	G4double fHRSMomentum; 		// HRS central momentum
	bool goodParticle;
	double sign; 	// y/phi sign flip for using transport function on left arm

	G4double fMinEKill;

	g4hrsTune* fTune;
	
	G4int fLHRS;
	G4int fRHRS;
	
	G4double fX0;
	G4double fY0;
	G4double fZ0;
	G4double fTh0;
	G4double fPh0;
	G4double fP0;		
	G4double fX0_tr;
	G4double fY0_tr;
	G4double fZ0_tr;
	G4double fTh0_tr;
	G4double fPh0_tr;
	G4double fP0_tr;

	float r0[5];
	
	int numTF;
	int numTFvar;
	G4double TFdata[4][12];

	int numVB;
	int numVar;
	G4double VBdata[14][12];
	G4String VBnames[14];
	
	int numZCrit;
        int numZCritVar;
        G4double ZCritData[23][5];
        G4String ZCritNames[23];

TFile *fout;
TH2D *all,*accept;

};

#endif//__REMOLLSTEPPINGACTION_HH
