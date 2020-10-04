#include "g4hrsGenNuclElastic.hh"

#include "CLHEP/Random/RandFlat.h"

#include "Randomize.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"

#include "g4hrsEvent.hh"
#include "g4hrsVertex.hh"
#include "g4hrsBeamTarget.hh"
#include "g4hrsMultScatt.hh"
#include "g4hrstypes.hh"
#include "g4hrsDatabase.hh"

#include <math.h>

#define Euler 0.5772157
#define NINTERVAL 3

g4hrsGenNuclElastic::g4hrsGenNuclElastic(){
    fTh_min =	1.*deg;
    fTh_max =	10.*deg;
    fPh_min = 	0.*deg;
    fPh_max = 	360.*deg;	  

    fE_min = 80.0*MeV; // Absolute minimum of electron energy
                            // to generate

    fApplyMultScatt = true;
    fBeamTarg = g4hrsBeamTarget::GetBeamTarget();

    fDatabase = new g4hrsDatabase(fBeamTarg->GetTargetMaterial());
    fDiamondDB = new g4hrsDatabase("Diamond");


}

g4hrsGenNuclElastic::~g4hrsGenNuclElastic(){
    delete fDatabase;
    fDatabase = NULL;
    delete fDiamondDB;
    fDiamondDB = NULL;
}

void g4hrsGenNuclElastic::SamplePhysics(g4hrsVertex *vert, g4hrsEvent *evt){
    // Generate ep event
    
    //  Crazy weighting for brem because ep cross section blows up at low Q2

    // Get initial beam energy instead of using other sampling
    double beamE = fBeamTarg->fBeamE;
    double Ekin  = beamE - electron_mass_c2;

    bool bypass_target = false;

    double bremcut = fBeamTarg->fEcut;

    // Approximation for Q2, just needs to be order of magnitude
    double effQ2 = 2.0*beamE*beamE*(1.0-cos(5.0*deg));

    // Internal radiation
    double int_bt = 0.75*(alpha/pi)*( log( effQ2/(electron_mass_c2*electron_mass_c2) ) - 1.0 );
    // External radiation
    double radlen = fBeamTarg->fRadLen;

    double bt;
    if( !bypass_target ){
	bt = (4.0/3.0)*(int_bt + radlen);
    } else {
	bt = 0.0;
    }

    double prob, prob_sample, sample, eloss, value;
    value = 1.0;
    prob = 1.- pow(bremcut/Ekin,bt) - bt/(bt+1.)*(1.- pow(bremcut/Ekin,bt+1.))
	+ 0.75*bt/(2.+bt)*(1.- pow(bremcut/Ekin,bt+2.));
    prob = prob/(1.- bt*Euler + bt*bt/2.*(Euler*Euler+pi*pi/6.)); // Gamma function 
    prob_sample = G4UniformRand();        // Random sampling 

    double env, ref;

    if (prob_sample <= prob) {//Bremsstrahlung has taken place!
	do {
	    sample = G4UniformRand();
	    eloss = fBeamTarg->fEcut*pow(Ekin/fBeamTarg->fEcut,sample);
	    env = 1./eloss;
	    value = 1./eloss*(1.-eloss/Ekin+0.75*pow(eloss/Ekin,2))*pow(eloss/Ekin,bt);

	    sample = G4UniformRand();
	    ref = value/env;
	} while (sample > ref);


	beamE -= eloss;
    }
   
    //RR -- from Hanjie
    beamE -= fBeamTarg->fIonLoss;


    if( beamE < electron_mass_c2 ){ 
	evt->SetEffCrossSection(0.0);
	evt->SetAsymmetry(0.0);
	return; 
    }

    // Set event information to our new sampling
    evt->fBeamE = beamE;
    evt->fBeamMomentum = evt->fBeamMomentum.unit()*sqrt(beamE*beamE - electron_mass_c2*electron_mass_c2);;

    ////////////////////////////////////////////////////////////////////////////////////////////

		


    // sample with 1.0/(1-cos)^2

	double inv_cthmax = 1./(1.-cos(fTh_max));
	double inv_cthmin = 1./(1.-cos(fTh_min));
	
	double sampv = 1./CLHEP::RandFlat::shoot(inv_cthmax, inv_cthmin);
	assert( -1.0 < sampv && sampv < 1.0 );

    	double th = acos(1.0-sampv);
    	double ph = CLHEP::RandFlat::shoot(fPh_min, fPh_max);

	// Weight for cross section to account for non-uniform sampling
	double dcosth = cos(fTh_min) - cos(fTh_max);
	double dphi = fPh_max - fPh_min;
	double norm = (inv_cthmin - inv_cthmax)/(cos(fTh_min) - cos(fTh_max));	
	double V = sampv*sampv*norm*dcosth*dphi;		

	double thisZ = vert->GetMaterial()->GetZ();
	double thisA = vert->GetMaterial()->GetA()/(g/mole);
	double M_nucl = thisA*(amu_c2/GeV)*GeV;			// GeV/GeV is redundant but the units should be explicit

	//final energy	
	double ef = (M_nucl*beamE)/(M_nucl + beamE*(1. - cos(th)));	

	// G4cout<<"\t A= "<<thisA<<" "<<thisZ<<G4endl;
	// G4cout<<"\t A= "<<int(thisA)<<" "<<int(thisZ)<<G4endl;

    double Q2  = 2.0*beamE*ef*(1.0-cos(th));
    evt->SetQ2( Q2 );

     //RR - For debugging purposes to make sure the flag worked
    //bool table = fDatabase->Table; G4cout << table << G4endl;


	// Get cross section from database (units of millibarns)
	// Multiply by millibarn, now it's in mm^2 (Geant's favorite length unit)
	double sigma = fDatabase->Interpolate(beamE,th,0,0)*millibarn;  
	if(int(thisA)==12){
	  sigma = fDiamondDB->Interpolate(beamE,th,0,0)*millibarn;  	
	}
    // Suppress too low angles from being generated
    // If we're in the multiple-scattering regime
    // the cross sections are senseless.  We'll define this 
    // as anything less than three sigma of the characteristic
    // width
    
    if( th < 3.0*fBeamTarg->fMS->GetPDGTh() ){
	sigma = 0.0;
    }

	evt->SetEffCrossSection(sigma*V);
	evt->SetWeight(V);

    if( vert->GetMaterial()->GetNumberOfElements() != 1 ){
	G4cerr << __FILE__ << " line " << __LINE__ << 
	    ": Error!  Some lazy programmer didn't account for complex materials in the moller process!" << G4endl;
	exit(1);
    }
	
	G4double APV = fDatabase->Interpolate(beamE,th,0,1);
	G4double APV1 = fDatabase->Interpolate(beamE,th,1,1);
	if(int(thisA)==12){
	  APV  = fDiamondDB->Interpolate(beamE,th,0,1); 
	  APV1 = fDiamondDB->Interpolate(beamE,th,1,1); 
	}

	G4double sens;
	if(APV > 0.) {
		sens = (APV - APV1)/APV;
	} else {
		sens = 0.;
	}
    
	evt->SetAsymmetry(APV);
	evt->SetSensitivity(sens);
    	evt->SetW2( M_nucl*M_nucl + 2.*M_nucl*(beamE - ef) - Q2 );

    // REradiate////////////////////////////////////////////////////////////////////////////
    // We're going to use the new kinematics for this guy

    int_bt = (alpha/pi)*( log( Q2/(electron_mass_c2*electron_mass_c2) ) - 1.0 );
    Ekin = ef - electron_mass_c2;;

    prob = 1.- pow(bremcut/Ekin, int_bt) - int_bt/(int_bt+1.)*(1.- pow(bremcut/Ekin,int_bt+1.))
	+ 0.75*int_bt/(2.+int_bt)*(1.- pow(bremcut/Ekin,int_bt+2.));
    prob = prob/(1.- int_bt*Euler + int_bt*int_bt/2.*(Euler*Euler+pi*pi/6.)); /* Gamma function */
    prob_sample = G4UniformRand();        /* Random sampling */

    if (prob_sample <= prob) {//Bremsstrahlung has taken place!
	do {
	    sample = G4UniformRand();
	    eloss = fBeamTarg->fEcut*pow(Ekin/fBeamTarg->fEcut,sample);
	    env = 1./eloss;
	    value = 1./eloss*(1.-eloss/Ekin+0.75*pow(eloss/Ekin,2))*pow(eloss/Ekin,bt);

	    sample = G4UniformRand();
	    ref = value/env;
	} while (sample > ref);

	ef = Ekin-eloss+electron_mass_c2;
	assert( ef > electron_mass_c2 );
    }

    ///////////////////////////////////////////////////////////////////////////////////////

    evt->ProduceNewParticle( G4ThreeVector(0.0, 0.0, 0.0), 
	    G4ThreeVector(ef*cos(ph)*sin(th), ef*sin(ph)*sin(th), ef*cos(th) ), 
	    "e-" );

    return;

}

G4double g4hrsGenNuclElastic::RadProfile(G4double eloss, G4double btt){
     double Ekin = fBeamTarg->fBeamE - electron_mass_c2;
     double retval = 1./eloss*(1.-eloss/Ekin+0.75*pow(eloss/Ekin,2))*pow(eloss/Ekin,btt);

     if( std::isnan(retval) || std::isinf(retval) ){
	 G4cerr << __FILE__ << " line " << __LINE__ << ": ERROR" << G4endl;
	 G4cerr << "Ekin " << Ekin/GeV << " GeV   btt = " << btt << " retval = " << retval << G4endl;
	 fprintf(stderr, "eloss = %e GeV\n", eloss/GeV);
     }

     assert( !std::isnan(retval) && !std::isinf(retval) );

     return retval;
}

G4double g4hrsGenNuclElastic::EnergNumInt(G4double btt, G4double a0, G4double b0){
    const int nbin = 1000;
    double sum = 0.0;
    double bremcut = fBeamTarg->fEcut;

    int j;
    double boolc[5] = {7.0, 32.0, 12.0, 32.0, 7.0};

    double a, b, thissum;

    for(int i =0;i<nbin;i++) {
	// Integrate over sample spacings that are logarithmic
	a = bremcut*pow(b0/a0, (((double) i)/nbin));
	b = bremcut*pow(b0/a0, (((double) i+1.0)/nbin));

	// Boole's rule
	thissum = 0.0;
	for( j = 0; j < 5; j++ ){
	    thissum +=  boolc[j]*RadProfile( (b-a)*j*0.25 + a, btt);
	}
	sum += thissum*(b-a)/90.0;
    }

     assert( !std::isnan(sum) && !std::isinf(sum) );

    return sum;
}















