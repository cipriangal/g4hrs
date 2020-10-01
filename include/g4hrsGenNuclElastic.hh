#ifndef __REMOLLGENPELASTIC_HH 
#define __REMOLLGENPELASTIC_HH 
/*!
 * ep elastic event generator
 *
 * Seamus Riordan
 * February 12, 2013
 *
 * Based heavily on previous work from mollersim
*/

#include "g4hrsVEventGen.hh"

class g4hrsBeamTarget;
class g4hrsDatabase;

class g4hrsGenNuclElastic : public g4hrsVEventGen {
    public:
	g4hrsGenNuclElastic();
	~g4hrsGenNuclElastic();

    private:
	void SamplePhysics(g4hrsVertex *, g4hrsEvent *);

	G4double RadProfile(G4double,G4double);
	G4double EnergNumInt(G4double,G4double,G4double);

	g4hrsBeamTarget *fBeamTarg;

        g4hrsDatabase *fDatabase;
        g4hrsDatabase *fDiamondDB;
};

#endif//__REMOLLGENPELASTIC_HH 
