#include "g4hrsGenBeam.hh"

#include "CLHEP/Random/RandFlat.h"

#include "g4hrsEvent.hh"
#include "g4hrsVertex.hh"
#include "g4hrsBeamTarget.hh"

#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include "g4hrstypes.hh"

#include <math.h>

g4hrsGenBeam::g4hrsGenBeam(){
    fApplyMultScatt = true;
    fBeamTarg = g4hrsBeamTarget::GetBeamTarget();

    fZpos = -5.0*m;
}

g4hrsGenBeam::~g4hrsGenBeam(){
}

void g4hrsGenBeam::SamplePhysics(g4hrsVertex *vert, g4hrsEvent *evt){
    // Get initial beam energy instead of using other sampling
    double beamE = fBeamTarg->fBeamE;
    evt->fBeamE = beamE;
    evt->fBeamMomentum = evt->fBeamMomentum.unit()*sqrt(beamE*beamE - electron_mass_c2*electron_mass_c2);;

    // Override target sampling z
    evt->fVertexPos.setZ( fZpos );

    G4double rasx = CLHEP::RandFlat::shoot( -2.0*mm, 2.0*mm);
    G4double rasy = CLHEP::RandFlat::shoot( -3.0*mm, 3.0*mm);

    evt->ProduceNewParticle( G4ThreeVector(rasx, rasy, -10.0*mm), 
	    evt->fBeamMomentum, 
	    "e-" );

    evt->SetEffCrossSection(0.0);
    evt->SetAsymmetry(0.0);

    evt->SetQ2(0.0);
    evt->SetW2(0.0);

    return;

}
