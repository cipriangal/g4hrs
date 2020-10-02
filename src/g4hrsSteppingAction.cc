#include "g4hrsSteppingAction.hh"
#include "g4hrsTransportFunction.hh"
#include "g4hrsTune.hh"
//#include "g4hrsSteppingActionMessenger.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SteppingManager.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4TouchableHandle.hh"
#include "G4TouchableHistory.hh"
#include "G4AffineTransform.hh"
#include "G4SystemOfUnits.hh"

g4hrsSteppingAction::g4hrsSteppingAction()
:drawFlag(false)
{

	numTF = 12;
	numTFvar = 4;

	numVB = 14;
	numVar = 6;
	VBnames[0] = "virtualBoundaryPhys_sen";
	VBnames[1] = "virtualBoundaryPhys_sm";
	VBnames[2] = "virtualBoundaryPhys_sex";
	VBnames[3] = "virtualBoundaryPhys_col";
	VBnames[4] = "virtualBoundaryPhys_q1en";
	VBnames[5] = "virtualBoundaryPhys_q1ex";
	VBnames[6] = "virtualBoundaryPhys_q2en";
	VBnames[7] = "virtualBoundaryPhys_q2ex";
	VBnames[8] = "virtualBoundaryPhys_den";
	VBnames[9] = "virtualBoundaryPhys_dex";
	VBnames[10] = "virtualBoundaryPhys_q3en";
	VBnames[11] = "virtualBoundaryPhys_q3ex";
	VBnames[12] = "virtualBoundaryPhys_vdc";
	VBnames[13] = "virtualBoundaryPhys_fp";	
	
	numZCrit = 24;
	numZCritVar = 6;
	ZCritNames[0] =  "virtualBoundaryPhys_zpinch1";
	ZCritNames[1] =  "virtualBoundaryPhys_zpinch2";
	ZCritNames[2] =  "virtualBoundaryPhys_zpinch3";
	ZCritNames[3] =  "virtualBoundaryPhys_ztarg";
	ZCritNames[4] =  "virtualBoundaryPhys_zfield0";
	ZCritNames[5] =  "virtualBoundaryPhys_zfield1";
	ZCritNames[6] =  "virtualBoundaryPhys_zfield2";
	ZCritNames[7] =  "virtualBoundaryPhys_zmidtosep";
	ZCritNames[8] =  "virtualBoundaryPhys_zsep1";
	ZCritNames[9] =  "virtualBoundaryPhys_zsep2";
	ZCritNames[10] = "virtualBoundaryPhys_zsep3";
	ZCritNames[11] = "virtualBoundaryPhys_zsep4";
        ZCritNames[12] = "virtualBoundaryPhys_zup1";
        ZCritNames[13] = "virtualBoundaryPhys_zup2";
        ZCritNames[14] = "virtualBoundaryPhys_zdown1";
        ZCritNames[15] = "virtualBoundaryPhys_zdown2";
        ZCritNames[16] = "virtualBoundaryPhys_zdown3";
        ZCritNames[17] = "virtualBoundaryPhys_zdown4";
        ZCritNames[18] = "virtualBoundaryPhys_zdown5";
        ZCritNames[19] = "virtualBoundaryPhys_zdown6";
        ZCritNames[20] = "virtualBoundaryPhys_zdown7";
        ZCritNames[21] = "virtualBoundaryPhys_zdown8";
        ZCritNames[22] = "virtualBoundaryPhys_zdown9";
        ZCritNames[23] = "virtualBoundaryPhys_zsieve";


///  new g4hrsSteppingActionMessenger(this);
	rad = 1.;
    fEnableKryptonite = true;
	fTransportFunction = new g4hrsTransportFunction();

	fTune = g4hrsTune::GetTune();	
	
	fSeptumAngle = fTune->septumAngle;
	fHRSAngle = fTune->HRSAngle;
	fHRSMomentum = fTune->HRSMomentum;
	nelements = 5;
	fMinEKill = 0.850*GeV;

}

void g4hrsSteppingAction::UserSteppingAction(const G4Step *aStep) {
    G4Track* fTrack = aStep->GetTrack();
    G4Material* material = fTrack->GetMaterial();

	// Don't continue in these materials
	// Don't track secondaries
	// Don't track low energy (<0.9 GeV) particles

	if(	((  material->GetName()=="Tungsten"
		||  material->GetName()=="Copper" ) 
		&&  fEnableKryptonite )
		|| fTrack->GetParentID() != 0
		|| fTrack->GetTotalEnergy() < fMinEKill ) {

			fTrack->SetTrackStatus(fStopAndKill);

		} 
		
//	G4cout << "Step " << fTrack->GetCurrentStepNumber() << " in " << material->GetName() << G4endl;

	////////////////////////////////////////////////
	// Virtual boundaries and transport functions //
	////////////////////////////////////////////////
	
	if(fTrack->GetParentID() == 0) { 
		G4StepPoint* prePoint = aStep->GetPreStepPoint(); 
		// Get hall coordinates in meters to output to ROOT	
		G4double x = fTrack->GetPosition().x()/1000.;
		G4double y = fTrack->GetPosition().y()/1000.;
		G4double z = fTrack->GetPosition().z()/1000.;
		// Get position 3-vector in millimeters for transforms
		G4ThreeVector position = fTrack->GetPosition();
		// Get momentum 3-vector
		G4ThreeVector momentum = fTrack->GetMomentum();

		G4String realVolName = fTrack->GetVolume()->GetName();

//		G4cout << std::setw(20) << realVolName << "\t" << material->GetName() << G4endl;		

		if(fTrack->GetCurrentStepNumber()<=1) {
		
			G4ThreeVector position0 = prePoint->GetPosition()/1000.;	
			G4ThreeVector momentum0 = prePoint->GetMomentum();

			fX0 = position0.x();
			fY0 = position0.y();
			fZ0 = position0.z();
			fTh0 = momentum0.getTheta();
			fPh0 = momentum0.getPhi();	
			fP0 = momentum0.mag();

	/*
			// Test transform with special position/angle	
			fX0 = 0.*mm;
			fY0 = 0.*mm;
			fZ0 = 0.;
			fP0 = 1.063*GeV;
			fTh0 = 5.*deg;
			fPh0 = 0.*deg;
			G4ThreeVector position0 = G4ThreeVector(fX0, fY0, fZ0);
			G4ThreeVector momentum0 = G4ThreeVector(fP0*sin(fTh0)*cos(fPh0),fP0*sin(fTh0)*sin(fPh0),fP0*cos(fTh0));
	*/


			if(momentum0.x() > 0.){
				fLHRS = 1;
				septum_angle = -fSeptumAngle; 
				hrs_angle = -fHRSAngle;
				sign = -1.;
			}
			if(momentum0.x() < 0.){
				fRHRS = 1;
				septum_angle = fSeptumAngle;
				hrs_angle = fHRSAngle;
				sign = +1.;
			}

			G4RotationMatrix rotate_targ;
			rotate_targ.rotateY(septum_angle);	
			rotate_targ.rotateZ(90.*deg);
							
			G4AffineTransform transportAxis_targ = G4AffineTransform(rotate_targ);
			G4AffineTransform transport_targ = transportAxis_targ.Inverse();

			//Transform
			G4ThreeVector position0_tr = transport_targ.TransformPoint(position0);
			G4ThreeVector momentum0_tr = transport_targ.TransformPoint(momentum0);


			fX0_tr = position0_tr.x();
			fY0_tr = position0_tr.y();
			fZ0_tr = position0_tr.z();
			fTh0_tr = momentum0_tr.x()/momentum0_tr.z();
			fPh0_tr = momentum0_tr.y()/momentum0_tr.z();
			fP0_tr = momentum0_tr.mag();
			r0[0] = (float)fX0_tr;
			r0[1] = (float)(fTh0_tr);
			r0[2] = (float)(fY0_tr*sign);			
			r0[3] = (float)(fPh0_tr*sign);
			r0[4] = (float)((fP0-fHRSMomentum)/fHRSMomentum);	

			// Only need to call transport functions once for event
			goodParticle = fTransportFunction->CallTransportFunction(r0, TFdata[0], TFdata[1], TFdata[2], TFdata[3]); 	
		
		}

	/* TROUBLESHOOTING 
		// Center of LHRS focal plane in HCS for coordinate transform test
		G4ThreeVector fphc = G4ThreeVector(21770.815*sin(fHRSAngle),8330.421,1103.79+21770.815*cos(fHRSAngle));
		// Unit momentum or central trajectory exiting dipole in HCS for coordinate transform test
		G4ThreeVector cthc = G4ThreeVector(cos(45.*deg)*sin(fHRSAngle),sin(45.*deg),cos(45.*deg)*cos(fHRSAngle));
	*/
		
		
		// Navigator of parallel world to get volume name (to identify virtual boundaries) and hall to transport coordinate transformation
		G4Navigator* fParallelNavigator = G4TransportationManager::GetTransportationManager()->GetNavigator("g4hrsparallel");
		G4String volName = fParallelNavigator->LocateGlobalPointAndSetup(position)->GetName();

		if(volName.find("virtualBoundaryPhys") != G4String::npos) { 

			// Must explicitly define transform for septum, as these virtual boundaries are in hall coordinate system
			//Rotations
			G4RotationMatrix rotate_sen, rotate_sm, rotate_sex;
			rotate_sen.rotateY(septum_angle);
			rotate_sen.rotateZ(90.*deg);
			rotate_sm.rotateY((septum_angle + hrs_angle)/2.);
			rotate_sm.rotateZ(90.*deg);
			rotate_sex.rotateY(hrs_angle);
			rotate_sex.rotateZ(90.*deg);
			//Translations
			double pivotZOffset = 115.12*cm;// RR match survey
			double septumZPosition = 70.0*cm; //match survey
			septumZPosition+=pivotZOffset; //put origin at target center, not Hall A pivot
			double septumLength = 74.*cm;
			double z_sen = septumZPosition - septumLength/2.;
			double x_sen = z_sen*tan(fSeptumAngle);
			//Length of chord connecting central trajectory at septum entrance and septum exit
			double chord = septumLength/cos((fHRSAngle+fSeptumAngle)/2.);
			//Radius of circle for central trajectory in ideal septum
			double R = sqrt((chord*chord)/(2.*(1-cos(fHRSAngle-fSeptumAngle))));
			//Center of aforementioned circle
			double xc = R*cos(fSeptumAngle) + x_sen;
			double zc = z_sen - R*sin(fSeptumAngle);
			//x,z position of central trajectory at septum middle
			double z_sm = septumZPosition;
			double x_sm = xc - sqrt(R*R - (z_sm - zc)*(z_sm - zc));
			//x,z position of central trajectory at septum exit
			double z_sex = septumZPosition + septumLength/2.;
			double x_sex = xc - sqrt(R*R - (z_sex - zc)*(z_sex - zc));
		
			if(fRHRS) {
				x_sen*=-1.;
				x_sm*=-1.;
				x_sex*=-1.;
			}

			G4ThreeVector translate_sen = G4ThreeVector(x_sen,0.,z_sen);
			G4ThreeVector translate_sm = G4ThreeVector(x_sm,0.,z_sm);
			G4ThreeVector translate_sex = G4ThreeVector(x_sex,0.,z_sex);
	/*
			//TROUBLESHOOTING
			G4ThreeVector senhc = G4ThreeVector(x_sen,-2.,z_sen);
			G4ThreeVector semhc = G4ThreeVector(x_sm,0.,z_sm);
			G4ThreeVector sexhc = G4ThreeVector(x_sex,0.,z_sex);
			G4ThreeVector senct = G4ThreeVector(sin(fSeptumAngle),0.,cos(fSeptumAngle));
			G4ThreeVector semct = G4ThreeVector(sin((fSeptumAngle+fHRSAngle)/2.),0.,cos((fSeptumAngle+fHRSAngle)/2.));
			G4ThreeVector sexct = G4ThreeVector(sin(fHRSAngle),0.,cos(fHRSAngle));
			G4ThreeVector sentr = position_transform_sen.TransformPoint(senhc);
			G4ThreeVector sencttr = momentum_transform_sen.TransformPoint(senct);
			G4ThreeVector semtr = position_transform_sm.TransformPoint(semhc);
			G4ThreeVector semcttr = momentum_transform_sm.TransformPoint(semct);
			G4ThreeVector sextr = position_transform_sex.TransformPoint(sexhc);
			G4ThreeVector sexcttr = momentum_transform_sex.TransformPoint(sexct);
	*/

			// Transform and record virtual detector data			
			for(int i = 0; i < numVB; i++) {

				G4AffineTransform position_transform;
				G4AffineTransform momentum_transform;

				if(volName.find(VBnames[i]) != G4String::npos ) {
					if(volName.find("virtualBoundaryPhys_sen") != G4String::npos ) {
						position_transform.SetNetTranslation(translate_sen);
						position_transform.SetNetRotation(rotate_sen);
						momentum_transform.SetNetRotation(rotate_sen);
						position_transform.Invert();
						momentum_transform.Invert();
					} else if(volName.find("virtualBoundaryPhys_sm") != G4String::npos )  {		
						position_transform.SetNetTranslation(translate_sm);
						position_transform.SetNetRotation(rotate_sm);
						momentum_transform.SetNetRotation(rotate_sm);
						position_transform.Invert();
						momentum_transform.Invert();
					} else if(volName.find("virtualBoundaryPhys_sex") != G4String::npos )  {		
						position_transform.SetNetTranslation(translate_sex);
						position_transform.SetNetRotation(rotate_sex);
						momentum_transform.SetNetRotation(rotate_sex);
						position_transform.Invert();
						momentum_transform.Invert();
					} else {
						//This gets the full transform (rotation + translation) that is required for transforming POSITION		
						position_transform = fParallelNavigator->GetGlobalToLocalTransform();
						//However, the MOMENTUM transform requires rotation ONLY!! 
						G4RotationMatrix momentum_rotation = position_transform.NetRotation();
						momentum_transform = G4AffineTransform(momentum_rotation);
					}

					// Hall coordinates
					VBdata[i][0] = x;
					VBdata[i][1] = y;
					VBdata[i][2] = z;
					VBdata[i][3] = momentum.theta()/rad;
					VBdata[i][4] = momentum.phi()/rad;
					VBdata[i][5] = momentum.mag();	
					// Transport coordinates
					G4ThreeVector pos_tr = position_transform.TransformPoint(position)/1000.;
					G4ThreeVector mom_tr = momentum_transform.TransformPoint(momentum);
					VBdata[i][6] = pos_tr.x();
					VBdata[i][7] = pos_tr.y();
					VBdata[i][8] = pos_tr.z();
					VBdata[i][9] = mom_tr.x()/mom_tr.z();
					VBdata[i][10] = mom_tr.y()/mom_tr.z();
					VBdata[i][11] = mom_tr.mag();		


				}
			}

			// Record tracking detector data
			for(int i = 0; i < numZCrit; i++) {
				if(volName.find(ZCritNames[i]) != G4String::npos) {
					ZCritData[i][0] = x;
					ZCritData[i][1] = y;
					ZCritData[i][2] = z;
					ZCritData[i][3] = momentum.theta()/rad;
					ZCritData[i][4] = momentum.phi()/rad;
                                        ZCritData[i][5] = momentum.mag();		
		}
			}

		} //end if volName contains virtualBoundaryPhys
	} //end if ParentID == 0
}
