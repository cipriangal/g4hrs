#include "g4hrsVEventGen.hh"

#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

#include "g4hrsBeamTarget.hh"
#include "g4hrsVertex.hh"
#include "g4hrsEvent.hh"
#include "g4hrsRun.hh"
#include "g4hrsRunData.hh"

g4hrsVEventGen::g4hrsVEventGen() {
    fBeamTarg = g4hrsBeamTarget::GetBeamTarget();
    fRunData  = g4hrsRun::GetRun()->GetData();

    fIsVPosHCSSet = false;
    fIsVPosTCSSet = false;
    fIsVMomHCSSet = false;
    fIsVMomTCSSet = false;
    
    //fSampType       = kMainTarget;
    fSampType       = kFullTarget;//FIXME
    
    fApplyMultScatt = true;
    fSeptumAngle = 5.*deg;
}

g4hrsVEventGen::~g4hrsVEventGen() {
}

g4hrsEvent *g4hrsVEventGen::GenerateEvent() {
    // Set up beam/target vertex
    g4hrsVertex vert   = fBeamTarg->SampleVertex(fSampType);

    /////////////////////////////////////////////////////////////////////
    // Create and initialize values for event
    g4hrsEvent *thisev = new g4hrsEvent();
    thisev->fVertexPos    = fBeamTarg->fVer;
    if( fApplyMultScatt ) {
        thisev->fBeamMomentum = fBeamTarg->fSampE*(fBeamTarg->fDir.unit());
    } else {
        thisev->fBeamMomentum = fBeamTarg->fSampE*G4ThreeVector(0.0, 0.0, 1.0);
    }
    /////////////////////////////////////////////////////////////////////

    SamplePhysics(&vert, thisev);

    PolishEvent(thisev);

    return thisev;
}
	
void g4hrsVEventGen::PolishEvent(g4hrsEvent *ev) {
    /*!
       Here it's our job to:
          Make sure the event is sane
          Apply multiple scattering effects to the final
        products if applicable
      Calculate rates from our given luminosity
      Calculate measured asymmetry from polarization
      Calculate vertex offsets
     */

    if( !ev->EventIsSane() ) {
        G4cerr << __FILE__ << " line " << __LINE__ << ":  Event check failed for generator " << fName << ".  Aborting" << G4endl;
        ev->Print();
        exit(1);
    }

    G4ThreeVector rotax      = (fBeamTarg->fDir.cross(G4ThreeVector(0.0, 0.0, 1.0))).unit();
    G4RotationMatrix msrot;
    msrot.rotate(fBeamTarg->fDir.theta(), rotax);

    std::vector<G4ThreeVector>::iterator iter;

    if( fApplyMultScatt ) {
        for( iter = ev->fPartRealMom.begin(); iter != ev->fPartRealMom.end(); iter++ ) {
            //  rotate direction vectors based on multiple scattering
            (*iter) *= msrot;
        }

        // Rotate position offsets due to multiple scattering
        for( iter = ev->fPartPos.begin(); iter != ev->fPartPos.end(); iter++ ) {
            //  rotate direction vectors based on multiple scattering
            (*iter) *= msrot;
        }
    }

    // Add base vertex
    for( iter = ev->fPartPos.begin(); iter != ev->fPartPos.end(); iter++ ) {
        (*iter) += ev->fVertexPos;
    }
    
    if ( ev->fRate == 0 ){// If the rate is set to 0 then calculate it using the cross section
    	ev->fRate  = ev->fEffXs*fBeamTarg->GetEffLumin()/((G4double) fRunData->GetNthrown());
    }
    else { // For LUND - calculate rate and cross section	
    	ev->fEffXs = ev->fRate*((G4double) fRunData->GetNthrown())/(fBeamTarg->GetEffLumin());
    	ev->fRate = ev->fRate/((G4double) fRunData->GetNthrown());
    }

    ev->fmAsym = ev->fAsym*fBeamTarg->fBeamPol;

	if(fIsVPosHCSSet && fIsVPosTCSSet) {
        	G4cerr << __FILE__ << " line " << __LINE__ << ": Cannot specify vertex position in both HCS and TCS" << fName << ".  Aborting" << G4endl;
		exit(1);				
	}	
	if(fIsVMomHCSSet && fIsVMomTCSSet) {
        	G4cerr << __FILE__ << " line " << __LINE__ << ": Cannot specify vertex momentum in both HCS and TCS" << fName << ".  Aborting" << G4endl;
		exit(1);				
	}	

		
	if(fIsVPosHCSSet || fIsVMomHCSSet || fIsVPosTCSSet || fIsVMomTCSSet) {

		
		G4RotationMatrix rotate_hall;
		rotate_hall.rotateZ(-90.*deg);
		rotate_hall.rotateY(fSeptumAngle);
		
		G4AffineTransform hallAxis_targ = G4AffineTransform(rotate_hall);
		G4AffineTransform hall_targ = hallAxis_targ.Inverse();

		if(fIsVPosHCSSet) {
			G4ThreeVector fSetVPosHCS;
			fSetVPosHCS = fSetVPos;
			for( iter = ev->fPartPos.begin(); iter != ev->fPartPos.end(); iter++ ) {
				(*iter) = fSetVPosHCS;
			}
		}


		if(fIsVPosTCSSet) {	
			G4ThreeVector fSetVPosHCS;
			fSetVPosHCS = hall_targ.TransformPoint(fSetVPos);
			for( iter = ev->fPartPos.begin(); iter != ev->fPartPos.end(); iter++ ) {
				(*iter) = fSetVPosHCS;
			}		
		}

		if(fIsVMomHCSSet) {
			G4double theta = fSetVMom[0]*deg;
			G4double phi = fSetVMom[1]*deg;
			G4double mom = (fBeamTarg->fBeamE)*(1.+fSetVMom[2]);
			G4ThreeVector fSetVMomHCS = G4ThreeVector(mom*sin(theta)*cos(phi),mom*sin(theta)*sin(phi),mom*cos(theta));		
			for( iter = ev->fPartMom.begin(); iter != ev->fPartMom.end(); iter++ ) {
				(*iter) = fSetVMomHCS;
			}
			for( iter = ev->fPartRealMom.begin(); iter != ev->fPartRealMom.end(); iter++ ) {
				(*iter) = fSetVMomHCS;
			}
		}


		if(fIsVMomTCSSet) {
			G4double theta = fSetVMom[0];
			G4double phi = fSetVMom[1];
			G4double mom = (fBeamTarg->fBeamE)*(1.+fSetVMom[2]);
			G4double pztr = mom/sqrt(theta*theta + phi*phi + 1.*1.);  
			G4double pxtr = pztr*theta;
			G4double pytr = pztr*phi;
			G4ThreeVector fSetVMomTCS = G4ThreeVector(pxtr, pytr, pztr);
			G4ThreeVector fSetVMomHCS = hall_targ.TransformPoint(fSetVMomTCS);
			for( iter = ev->fPartMom.begin(); iter != ev->fPartMom.end(); iter++ ) {
				(*iter) = fSetVMomHCS;
			}
			for( iter = ev->fPartRealMom.begin(); iter != ev->fPartRealMom.end(); iter++ ) {
				(*iter) = fSetVMomHCS;
			}
		}

	}

    return;
}

