#include "g4hrsEventAction.hh"
#include "g4hrsGenericDetectorHit.hh"
#include "g4hrsGenericDetectorSum.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

#include "g4hrsIO.hh"


g4hrsEventAction::g4hrsEventAction() {
  fFillTree=true;
}

g4hrsEventAction::~g4hrsEventAction(){
}


void g4hrsEventAction::BeginOfEventAction(const G4Event* ev){

	fIO->ClearVirtualBoundaryData();
 
	if(ev->GetEventID() % 100000 == 0) {
		G4cout << ev->GetEventID() << "\n";
	}
	fFillTree=false;
  // Start timer at event 0
  if (ev->GetEventID() == 0) fTimer.Start();
  // Pretty ongoing status
  if ((ev->GetEventID() % 1) == 0) {
    // Stop timer (running timer cannot be read)
    fTimer.Stop();
    // Print event number
//    G4cout << "Event " << ev->GetEventID();
    // Only print duration per event when meaningful (avoid division by zero)
    if (ev->GetEventID() > 0)
//      G4cout << " (" << std::setprecision(3) << std::fixed
  //      << 1000.*fTimer.GetRealElapsed()/1000.0 << " ms/event)";
    // Carriage return without newline
//    G4cout << "\r" << std::flush;
    // Start timer again
    fTimer.Start();
  }
}

void g4hrsEventAction::EndOfEventAction(const G4Event* evt ) {
  //G4SDManager   *SDman = G4SDManager::GetSDMpointer();
  G4HCofThisEvent *HCE = evt->GetHCofThisEvent();

  G4VHitsCollection *thiscol;

  if( HCE ){

      // Traverse all hit collections, sort by output type
      for( int hcidx = 0; hcidx < HCE->GetCapacity(); hcidx++ ){
          thiscol = HCE->GetHC(hcidx);
          if(thiscol){ // This is NULL if nothing is stored
              // Dyanmic cast to test types, process however see fit and feed to IO

              ////  Generic Detector Hits ///////////////////////////////////
              if( g4hrsGenericDetectorHitsCollection *thiscast = 
                      dynamic_cast<g4hrsGenericDetectorHitsCollection *>(thiscol)){
                  for( unsigned int hidx = 0; hidx < thiscast->GetSize(); hidx++ ){
                      fIO->AddGenericDetectorHit((g4hrsGenericDetectorHit *) 
                              thiscast->GetHit(hidx) );	  
                  }
              }

              ////  Generic Detector Sum ////////////////////////////////////
              if( g4hrsGenericDetectorSumCollection *thiscast = 
                      dynamic_cast<g4hrsGenericDetectorSumCollection *>(thiscol)){
                  for( unsigned int hidx = 0; hidx < thiscast->GetSize(); hidx++ ){
                      fIO->AddGenericDetectorSum((g4hrsGenericDetectorSum *) 
                              thiscast->GetHit(hidx) );
                  }
              }

          }
      }

  }

  fIO->SetVirtualBoundaryData();

  //  G4cout<<evt->GetEventID()<<G4endl;

  // Fill tree and reset buffers
  if(fFillTree){
    //G4cout<<"saved"<<G4endl;
    fIO->FillTree();
    fIO->Flush();
  }
  //std::cin.ignore();
  return;
}

