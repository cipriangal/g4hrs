#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "g4hrsBeamTarget.hh"
#include "g4hrsMultScatt.hh"
#include "g4hrsMaterial.hh"

#include <math.h>

#define __MAX_MAT 100
#define Euler 0.5772157

g4hrsBeamTarget *g4hrsBeamTarget::gSingleton = NULL;

g4hrsBeamTarget::g4hrsBeamTarget(){
    gSingleton = this;
    fTargVol = NULL;

    UpdateInfo();

    fOldRaster = true;
    fRasterX = fRasterY = 5.0*mm;
    fX0 = fY0 = fTh0 = fPh0 = fdTh = fdPh = 0.0;

    fCorrTh = fCorrPh = 0.0;

    fMS = new g4hrsMultScatt();

    fBeamE   = gDefaultBeamE;
    fBeamPol = gDefaultBeamPol;

    fBeamCurr = gDefaultBeamCur;

    fEcut = 1e-6*MeV;

    fDefaultMat = new G4Material("Default_proton"   , 1., 1.0, 1e-19*g/mole);

    fAlreadyWarned = false;

    fTargetMaterial = "Pb208";

}

g4hrsBeamTarget::~g4hrsBeamTarget(){
}

g4hrsBeamTarget *g4hrsBeamTarget::GetBeamTarget() {
    if( gSingleton == NULL ){
	gSingleton = new g4hrsBeamTarget();
    }
    return gSingleton;
}

G4String g4hrsBeamTarget::GetTargetMaterial() {

 	return fTargetMaterial;
	
}

G4double g4hrsBeamTarget::GetEffLumin(){
	// Beam current in e-/sec
	G4double electronCurrent = (fBeamCurr/ampere)/(e_SI);
	// Effective material length in 1/mm^2 (see calculation below)
	// Yields luminosity in e-/(s mm^2)
	G4double lumin = fEffMatLen*electronCurrent;
    	return lumin;
}

void g4hrsBeamTarget::UpdateInfo(){
    std::vector<G4VPhysicalVolume *>::iterator it;

    fTargLength   = 0.0;
    fTotalLength   = 0.0;
    fAllVols.clear();

    double thiszlen = -1e9;

    if( fTargVol ){
        if( !dynamic_cast<G4Tubs *>( fTargVol->GetLogicalVolume()->GetSolid() )  &&
            !dynamic_cast<G4Box *>( fTargVol->GetLogicalVolume()->GetSolid() ) 
          ){
            G4cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
                ":  Target volume not made of G4Tubs or G4Box" << G4endl; 
            exit(1);
        } 

        if( dynamic_cast<G4Tubs *>( fTargVol->GetLogicalVolume()->GetSolid() ) ){
            thiszlen = ((G4Tubs *) fTargVol->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0;
        }
        if( dynamic_cast<G4Box *>( fTargVol->GetLogicalVolume()->GetSolid() ) ){
            thiszlen = ((G4Box *) fTargVol->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0;
        }

        printf("[g4hrsBeamTarget::UpdateInfo] Target volume material %s\n", fTargVol->GetLogicalVolume()->GetMaterial()->GetName().c_str() );

        fTargLength  = thiszlen*fTargVol->GetLogicalVolume()->GetMaterial()->GetDensity();
        fTotalLength = thiszlen*fTargVol->GetLogicalVolume()->GetMaterial()->GetDensity();
    }

    // First upstream
    for(it = fUpstreamVols.begin(); it != fUpstreamVols.end(); it++ ){

	// Assume everything is non-nested tubes
	if( !dynamic_cast<G4Tubs *>( (*it)->GetLogicalVolume()->GetSolid() )  &&
            !dynamic_cast<G4Box *>( (*it)->GetLogicalVolume()->GetSolid() ) 
            ){
	    G4cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		":  Target volume not made of G4Tubs or G4Box" << G4endl; 
	    exit(1);
	}

	if( dynamic_cast<G4Tubs *>( (*it)->GetLogicalVolume()->GetSolid() ) ){
            thiszlen = ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0;
        }
	if( dynamic_cast<G4Box *>( (*it)->GetLogicalVolume()->GetSolid() ) ){
            thiszlen = ((G4Box *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0;
        }

        printf("[g4hrsBeamTarget::UpdateInfo] Target volume material %s\n", (*it)->GetLogicalVolume()->GetMaterial()->GetName().c_str() );

        fTotalLength  += thiszlen*(*it)->GetLogicalVolume()->GetMaterial()->GetDensity();

        fAllVols.push_back(*it);
    }

    if( fTargVol ){
        fAllVols.push_back(fTargVol);
    }

    // Finally Downstream
    for(it = fDownstreamVols.begin(); it != fDownstreamVols.end(); it++ ){

	// Assume everything is non-nested tubes
	if( !dynamic_cast<G4Tubs *>( (*it)->GetLogicalVolume()->GetSolid() )  &&
            !dynamic_cast<G4Box *>( (*it)->GetLogicalVolume()->GetSolid() ) 
            ){
	    G4cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		":  Target volume not made of G4Tubs or G4Box" << G4endl; 
	    exit(1);
	}

	if( dynamic_cast<G4Tubs *>( (*it)->GetLogicalVolume()->GetSolid() ) ){
            thiszlen = ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0;
        }
	if( dynamic_cast<G4Box *>( (*it)->GetLogicalVolume()->GetSolid() ) ){
            thiszlen = ((G4Box *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0;
        }

        printf("[g4hrsBeamTarget::UpdateInfo] Target volume material %s\n", (*it)->GetLogicalVolume()->GetMaterial()->GetName().c_str() );

        fTotalLength  += thiszlen*(*it)->GetLogicalVolume()->GetMaterial()->GetDensity();
        
        fAllVols.push_back(*it);
    }

    return;
}


void g4hrsBeamTarget::SetTargetLen(G4double z){
    std::vector<G4VPhysicalVolume *>::iterator it;

    if( fTargVol ){
        if( !dynamic_cast<G4Tubs *>( fTargVol->GetLogicalVolume()->GetSolid() )  &&
                !dynamic_cast<G4Box *>( fTargVol->GetLogicalVolume()->GetSolid() ) 
          ){
            G4cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
                ":  Target volume not made of G4Tubs or G4Box" << G4endl; 
            exit(1);
        } 


	G4GeometryManager::GetInstance()->OpenGeometry(fTargVol);
        if( dynamic_cast<G4Tubs *>( fTargVol->GetLogicalVolume()->GetSolid() ) ){
            ((G4Tubs *) fTargVol->GetLogicalVolume()->GetSolid())->SetZHalfLength(z/2.0);
        }
        if( dynamic_cast<G4Box *>( fTargVol->GetLogicalVolume()->GetSolid() ) ){
            ((G4Box *) fTargVol->GetLogicalVolume()->GetSolid())->SetZHalfLength(z/2.0);
        }
	G4GeometryManager::GetInstance()->CloseGeometry(true, false, fTargVol);
    }

    // Move positions of upstream and downstream volumes
    for(it = fUpstreamVols.begin(); it != fUpstreamVols.end(); it++ ){
	G4GeometryManager::GetInstance()->OpenGeometry((*it));
        double zpos = (*it)->GetFrameTranslation().z();
        (*it)->SetTranslation( G4ThreeVector(0.0, 0.0, zpos-z/2.0) );
	G4GeometryManager::GetInstance()->CloseGeometry(true, false, (*it));
    }
    for(it = fDownstreamVols.begin(); it != fDownstreamVols.end(); it++ ){
	G4GeometryManager::GetInstance()->OpenGeometry((*it));

        double zpos = (*it)->GetFrameTranslation().z();
        (*it)->SetTranslation( G4ThreeVector(0.0, 0.0, zpos+z/2.0) );
	G4GeometryManager::GetInstance()->CloseGeometry(true, false, (*it));
    }

    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->GeometryHasBeenModified();

    UpdateInfo();
}

void g4hrsBeamTarget::SetTargetPos(G4double z){
    std::vector<G4VPhysicalVolume *>::iterator it;

    double targz = -1e9;
    if( fTargVol ){
        targz = fTargVol->GetFrameTranslation().z();
        fTargVol->SetTranslation( G4ThreeVector(0.0, 0.0, z ) );
    }
	
    for(it = fUpstreamVols.begin(); it != fUpstreamVols.end(); it++ ){
	G4GeometryManager::GetInstance()->OpenGeometry((*it));
	G4String volname = (*it)->GetName();
        double zpos = (*it)->GetFrameTranslation().z();
        (*it)->SetTranslation( G4ThreeVector(0.0, 0.0, z+(zpos-targz) ) );
	G4GeometryManager::GetInstance()->CloseGeometry(true, false, (*it));
    }
    for(it = fDownstreamVols.begin(); it != fDownstreamVols.end(); it++ ){
	G4GeometryManager::GetInstance()->OpenGeometry((*it));

	G4String volname = (*it)->GetName();
        double zpos = (*it)->GetFrameTranslation().z();
        (*it)->SetTranslation( G4ThreeVector(0.0, 0.0, z+(zpos-targz)) );
	G4GeometryManager::GetInstance()->CloseGeometry(true, false, (*it));
    }


    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->GeometryHasBeenModified();

    UpdateInfo();
}

void g4hrsBeamTarget::SetTargetMaterial(G4String targMat) {

	fTargetMaterial = targMat;
	G4Material* targ_material;
    	mMaterialManager = g4hrsMaterial::GetHRSMaterialManager();
	if(fTargetMaterial == "Pb208") {
		targ_material = mMaterialManager->lead208;	
	} else if(fTargetMaterial == "C12") {
		targ_material = mMaterialManager->carbon12;
	} else if(fTargetMaterial == "Ca40") {
		targ_material = mMaterialManager->calcium40;	
	} else if(fTargetMaterial == "Ca48") {
		targ_material = mMaterialManager->calcium48;	
	} else if(fTargetMaterial == "Vacuum") {
		targ_material = mMaterialManager->vacuum;
	} else if(fTargetMaterial == "Galaxy") {
                targ_material = mMaterialManager->galaxy;
        } else if(fTargetMaterial == "Diamond") { 
                targ_material = mMaterialManager->diamond;
        }  else {
		G4cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		":  Invalid target material selected" << G4endl; 
		exit(1);
	}

	G4GeometryManager::GetInstance()->OpenGeometry(fTargVol);
	((G4LogicalVolume *) fTargVol->GetLogicalVolume())->SetMaterial(targ_material);
	G4GeometryManager::GetInstance()->CloseGeometry(true, false, fTargVol);
	G4RunManager* runManager = G4RunManager::GetRunManager();
	runManager->GeometryHasBeenModified();
	UpdateInfo();
}


////////////////////////////////////////////////////////////////////////////////////////////
//  Sampling functions

g4hrsVertex g4hrsBeamTarget::SampleVertex(SampType_t samp){
    g4hrsVertex thisvert;

    G4double rasx = CLHEP::RandFlat::shoot( fX0 - fRasterX/2.0, fX0 + fRasterX/2.0);
    G4double rasy = CLHEP::RandFlat::shoot( fY0 - fRasterY/2.0, fY0 + fRasterY/2.0);
    G4double ztrav, len;

    // Sample where along target weighted by density (which roughly corresponds to A
    // or the number of electrons, which is probably good enough for this

    // Figure out how far along the target we got
    switch( samp ){
	case kMainTarget: 
	    fSampLen = fTargLength;
	    break;

    case kWalls:
	    G4cerr << "ERROR" << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		": scattering from cell walls has been specified, but handling not implemented" << G4endl;
	    exit(1);
	    break;

	    /*
	case kWalls:
	    fSampLen = fTotalLength-fLH2Length;
	    break;
	    */
	case kFullTarget:
	    fSampLen = fTotalLength;
	    break;
    }
    const bool debugPrint=false;
    //debug block
    if(debugPrint){
       G4cout<<__PRETTY_FUNCTION__<<" "<<__LINE__<<G4endl
    	  <<"\t tgt L*density "<<fTargLength<<" total L*density "<<fTotalLength<<G4endl
    	  <<"\tsampling length*density (should include diamond foils) "<<fSampLen<<G4endl;
    }

    std::vector<G4VPhysicalVolume *>::iterator it;

    //calculate total z length of the target
    G4double zzSampLen(0);
    //G4cout<<"Volumes "<<fAllVols.size()<<G4endl;
    for(it = fAllVols.begin(); it != fAllVols.end(); it++ )
      zzSampLen += ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0;
    //sample flat in z
    G4double zzSelectedZ = CLHEP::RandFlat::shoot(0.0, zzSampLen);

    G4Material *mat;
    G4double rhoLen(0);
    G4double zzLen(0);
    bool foundIt(false);

    //calculate the length*density up to the point where the sampling was done
    for(it = fAllVols.begin(); it != fAllVols.end() && !foundIt; it++ ){
      G4double volLen = ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0;
      mat = (*it)->GetLogicalVolume()->GetMaterial();
      if(debugPrint){
	G4cout<<"\t"<<mat->GetName()<<" "<< (*it)->GetLogicalVolume()->GetSolid()->GetName()<<" "<<mat->GetDensity()<<G4endl;
	G4cout<<"\t B "<<zzSelectedZ<<" "<<zzLen<<" "<<volLen<<" "<<rhoLen<<" "<<foundIt<<G4endl;
      }
      if( zzSelectedZ - zzLen <= volLen ){
	rhoLen += (zzSelectedZ-zzLen) * mat->GetDensity();
	foundIt=true;
      }else{
	rhoLen += ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0*mat->GetDensity();
	zzLen  += volLen;
      }
      if(debugPrint)
	G4cout<<"\t A "<<zzSelectedZ<<" "<<zzLen<<" "<<volLen<<" "<<rhoLen<<" "<<foundIt<<G4endl;
    }

    //assign that to the ztrav instead of sampling
    ztrav = rhoLen;
    //ztrav = CLHEP::RandFlat::shoot(0.0, fSampLen);
    if(debugPrint)
      G4cout<<"&& "<<zzSampLen<<" "<<zzSelectedZ<<" "<<rhoLen<<" "<<fSampLen<<G4endl;

    G4double zinvol(0);

    G4double cumz   = 0.0;
    G4double radsum = 0.0;

    int      nmsmat = 0;
    double   msthick[__MAX_MAT];
    double   msA[__MAX_MAT];
    double   msZ[__MAX_MAT];

    // Figure out the material we are in and the radiation length we traversed
    bool foundvol = false;

    for(it = fAllVols.begin(); it != fAllVols.end() && !foundvol; it++ ){
        mat = (*it)->GetLogicalVolume()->GetMaterial();

        len = ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0*mat->GetDensity();
        switch( samp ){
            case kMainTarget: 
                /*
                   if( !isLH2 ){
                   radsum += len/mat->GetDensity()/mat->GetRadlen();
                   } else {
                   */
                foundvol = true;
                zinvol = ztrav/mat->GetDensity();
                radsum += zinvol/mat->GetRadlen();
                //}
                break;

            case kWalls:

                G4cerr << "WARNING " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
                    ": scattering from cell walls has been specified, but handling not implemented" << G4endl;
                /*
                   if( isLH2 ){
                   radsum += len/mat->GetDensity()/mat->GetRadlen();
                   } else {
                   if( ztrav - cumz < len ){
                   foundvol = true;
                   zinvol = (ztrav - cumz)/mat->GetDensity();
                   radsum += zinvol/mat->GetRadlen();
                   } else {
                   radsum += len/mat->GetDensity()/mat->GetRadlen();
                   cumz   += len;
                   }
                   }
                   */
                break;

            case kFullTarget:
                if( ztrav - cumz < len ){
                    foundvol = true;
                    zinvol = (ztrav - cumz)/mat->GetDensity();
                    radsum += zinvol/mat->GetRadlen();
                } else {
                    radsum += len/mat->GetDensity()/mat->GetRadlen();
                    cumz   += len;
                }
                break;
        }
	if(debugPrint)
	  G4cout<<"\t"<<mat->GetName()<<" "<<foundvol<<" "<<ztrav<<" "<<cumz<<" "<<zinvol<<" "<<len<<G4endl;

        if( mat->GetBaseMaterial() ){
            G4cerr << __FILE__ << " " << __PRETTY_FUNCTION__ << ":  The material you're using isn't" <<
                " defined in a way we can use for multiple scattering calculations" << G4endl;
            G4cerr << "Aborting" << G4endl; 
            exit(1);
        }

        if( foundvol ){
            // For our vertex
            thisvert.fMaterial = mat;
            thisvert.fRadLen   = radsum;

            // For our own info
            fTravLen = zinvol;
            fRadLen = radsum;

            fVer    = G4ThreeVector( rasx, rasy, 
                    zinvol - (*it)->GetFrameTranslation().z()  
                    - ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength() );
	    if(debugPrint){
	      G4cout<<"\t"<<(*it)->GetLogicalVolume()->GetSolid()->GetName()<<G4endl;
	      G4cout<<zinvol<<" "<<(*it)->GetFrameTranslation().z()
		    <<" "<<((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()
		    <<" "<<zinvol - (*it)->GetFrameTranslation().z() - ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()
		    <<G4endl;
	      std::cin.ignore();
	      // std::ofstream fout("z_pos.txt",std::ifstream::app);
	      // //fout<<zzSelectedZ<<std::endl;
	      // fout<<zinvol<<std::endl;
	      // fout.close();	      
	    }

            G4double masssum = 0.0;
            const G4int *atomvec = mat->GetAtomsVector();
            const G4ElementVector *elvec = mat->GetElementVector();
            const G4double *fracvec = mat->GetFractionVector();

            for( unsigned int i = 0; i < elvec->size(); i++ ){
                // FIXME:  Not sure why AtomsVector would ever return null
                // but it does - SPR 2/5/13.  Just going to assume unit
                // weighting for now if that is the case
                if( atomvec ){
                    masssum += (*elvec)[i]->GetA()*atomvec[i];
                } else {
                    masssum += (*elvec)[i]->GetA();
                }
                msthick[nmsmat] = mat->GetDensity()*zinvol*fracvec[i];
                msA[nmsmat] = (*elvec)[i]->GetA()*mole/g;
                msZ[nmsmat] = (*elvec)[i]->GetZ();
		//G4cout<<"~~ "<<nmsmat<<" "<<msthick[nmsmat]<<" "<<msA[nmsmat]<<" "<<msZ[nmsmat]<<G4endl;
                nmsmat++;
            }

		// target density in grams per cubic millimeter
		G4double rho = mat->GetDensity()/(g/mm3);
		// target length in millimeters
		G4double t = ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0/mm;
		// target atomic mass number in grams per mole
		G4double Amass = mat->GetA()/(g/mole);

		// Effective material length in 1/mm^2
            	fEffMatLen = rho*t*(Avogadro/Amass);

        } else {
            const G4ElementVector *elvec = mat->GetElementVector();
            const G4double *fracvec = mat->GetFractionVector();
            for( unsigned int i = 0; i < elvec->size(); i++ ){

                msthick[nmsmat] = len*fracvec[i];
                msA[nmsmat] = (*elvec)[i]->GetA()*mole/g;
                msZ[nmsmat] = (*elvec)[i]->GetZ();
		//G4cout<<"@@ "<<nmsmat<<" "<<msthick[nmsmat]<<" "<<msA[nmsmat]<<" "<<msZ[nmsmat]<<G4endl;
                nmsmat++;
            }
        }
    }


    if( !foundvol ){
        if( !fAlreadyWarned ){
            G4cerr << "WARNING: " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ": Could not find sampling volume" << G4endl;
            fAlreadyWarned = true;
        }

        thisvert.fMaterial = fDefaultMat;
        thisvert.fRadLen   = 0.0;
    }

    // Sample multiple scattering + angles
    G4double msth, msph;
    G4double bmth, bmph;

    if( nmsmat > 0 ){
      //G4cout<<"?? "<<nmsmat<<" "<<msthick[0]<<" "<<msA[0]<<" "<<msZ[0]<<G4endl;
        fMS->Init( fBeamE, nmsmat, msthick, msA, msZ );
        msth = fMS->GenerateMSPlane();
        msph = fMS->GenerateMSPlane();
    } else {
        msth = 0.0;
        msph = 0.0;
    }

    assert( !std::isnan(msth) && !std::isnan(msph) );

    if(fOldRaster){
        bmth = CLHEP::RandGauss::shoot(fTh0, fdTh);
        bmph = CLHEP::RandGauss::shoot(fPh0, fdPh);

        if( fRasterX > 0 ){ bmth += fCorrTh*(rasx-fX0)/fRasterX/2; }
        if( fRasterY > 0 ){ bmph += fCorrPh*(rasy-fY0)/fRasterY/2; }

        // Initial direction
        fDir = G4ThreeVector(0.0, 0.0, 1.0);

        fDir.rotateY( bmth); // Positive th pushes to positive X (around Y-axis)
        fDir.rotateX(-bmph); // Positive ph pushes to positive Y (around X-axis)
    } else{
        G4ThreeVector bmVec = G4ThreeVector(fVer.x(),fVer.y(),-1*(-8000.0*mm-fVer.z())); // in mm
        fDir = G4ThreeVector(bmVec.unit());
    }

    fDir.rotateY(msth);
    fDir.rotateX(msph);

    // Sample beam energy based on radiation
    // We do this so it doesn't affect the weighting
    //
    // This can be ignored and done in a generator by itself

    G4double  Ekin = fBeamE - electron_mass_c2;
    G4double  bt   = fRadLen*4.0/3.0;
    G4double  prob_sample, eloss, sample, env, value, ref;

    G4double prob = 1.- pow(fEcut/Ekin,bt) - bt/(bt+1.)*(1.- pow(fEcut/Ekin,bt+1.))
        + 0.75*bt/(2.+bt)*(1.- pow(fEcut/Ekin,bt+2.));
    prob = prob/(1.- bt*Euler + bt*bt/2.*(Euler*Euler+pi*pi/6.)); /* Gamma function */

    prob_sample = G4UniformRand();

    if (prob_sample <= prob) {
        do {
            sample = G4UniformRand();
            eloss = fEcut*pow(Ekin/fEcut,sample);
            env = 1./eloss;
            value = 1./eloss*(1.-eloss/Ekin+0.75*pow(eloss/Ekin,2))*pow(eloss/Ekin,bt);

            sample = G4UniformRand();
            ref = value/env;
        } while (sample > ref);

        fSampE = fBeamE - eloss;
        assert( fSampE > electron_mass_c2 );
    } else {
        fSampE = fBeamE;
    }


    thisvert.fBeamE = fSampE;

    assert( fBeamE >= electron_mass_c2 );

    //RR -- from Hanjie. She calculates ionization losses here since this class has all the materials to be looped

    fIonLoss = 0.0;
    G4double ion_mostprob = 0.0;  // the ionization most probable energy loss dE/dx (MeV/(g/cm2)) 
    bool pathend = false;

    G4double kai = 0.15355; // 2*pi*Na*re^2*me*c^2*z;
    G4double vgamma = thisvert.fBeamE/electron_mass_c2;
    G4double vbeta = sqrt(1. - 1./pow(vgamma,2));
    G4double vlen = 0;
    G4double zcum = 0;
    G4double zztrav = ztrav/(g/cm2);

    for(it = fAllVols.begin(); it != fAllVols.end() && !pathend; it++ ){
        mat = (*it)->GetLogicalVolume()->GetMaterial();
        vlen = ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0*mat->GetDensity();
        G4double vvlen = vlen / (g/cm2);

   // calculate the most probable energy loss for ionization dE/dx based on the "Particle Detectors C.Grupen B.Shwartz" eq(1.22)
   // ionization parameters
     G4IonisParamMat *ion_par = mat->GetIonisation();
     G4double I = ion_par->GetMeanExcitationEnergy();
     G4double X0 = ion_par->GetX0density();
     G4double X1 = ion_par->GetX1density();
     G4double C = ion_par->GetCdensity();
     G4double Aden = ion_par->GetAdensity();
     G4double Mden = ion_par->GetMdensity();
  
     G4double del_den = 0.0; // density effect PDG 2-18 Eq. (33.7) same as what is in G4hBetheBlochModel::BetheBlochFormula
         G4double X = log(vgamma*vbeta)/log(10);
        if ( X < X0 ) {
             del_den = 0.0 ;
        } else {
             del_den = 2*log(10)*X - C ;
             if ( X < X1 ) del_den += Aden*pow((X1-X),Mden) ;
        }

        G4double term1 = log(2.*electron_mass_c2*vgamma*vgamma*vbeta*vbeta/I);

        if( (zztrav - zcum) < vlen){
          G4double AA = mat->GetA() *mole/g;
          G4double Kai = kai * mat->GetZ()/AA/(vbeta*vbeta)*(zztrav-zcum);

          G4double term2 = log(Kai/I);
          ion_mostprob += Kai*( term1 + term2 +0.2 - vbeta*vbeta - del_den );

   // method 2: geant4 ionisation, should be mean energy not the most probably energy
   // inclue the shell correction, which should be needed for low energy beam (so that the shell electron is not in rest anymore)
   //
   //          G4ParticleTable* aparticleTable = G4ParticleTable::GetParticleTable();
   //          G4ParticleDefinition* aparticle = aparticleTable->FindParticle("e-");

   //        G4hBetheBlochModel *betheBlochModel = new G4hBetheBlochModel("Bethe-Bloch") ;
  //        G4double ionlossBB = betheBlochModel->TheValue(aparticle,mat, thisvert.fBeamE);

  // method 3: ionization energy loss (mean energy lossi dE/dx get from bob)
  //    double pb_x0 = 1.123 ; // MeV/(g/cm2) 
  //    ion_eloss = pb_x0 * fRadLen *6.39358;  

 //    double c12_x0 = 1.745 ; // MeV/(g/cm2)
 //    ion_eloss = c12_x0 * fRadLen * 42.70;

         pathend = true;
      } else { 
          G4double AA = mat->GetA() *mole/g;
          G4double Kai = kai * mat->GetZ()/AA/(vbeta*vbeta)*vvlen;
   
          G4double term2 = log(Kai/I);
          ion_mostprob += Kai*( term1 + term2 + 0.2 - vbeta*vbeta - del_den );
          zcum += vvlen;
      }

      }
      fIonLoss = ion_mostprob;
   

    return thisvert;
}


