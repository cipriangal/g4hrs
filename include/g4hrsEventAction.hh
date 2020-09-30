
#ifndef g4hrsEventAction_h
#define g4hrsEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Timer.hh"

#include "globals.hh"

class G4Event;
class g4hrsIO;

class g4hrsEventAction : public G4UserEventAction
{
  public:
    g4hrsEventAction();
    virtual ~g4hrsEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
	
    void SetIO( g4hrsIO *io ){ fIO = io; }

  void SetFill(G4bool val){fFillTree=val;}
  private:
  //  G4int gemCollID, hcalCollID, bbcalCollID;

    double fGEMres;

    g4hrsIO *fIO;

    // Timer for benchmarking of simulation time per event
    G4Timer fTimer;
  G4bool fFillTree;
  public:
};

#endif

    
