#ifndef __REMOLLBEAMTARGET_HH
#define __REMOLLBEAMTARGET_HH

#include "g4hrstypes.hh"
#include "g4hrsglobs.hh"
#include "g4hrsVertex.hh"
#include "G4ThreeVector.hh"
#include <vector>

/*!
     Class that contains information on 
     the beam and target.  It needs to be
     aware of and consistant with what is
     in the geometry.

     This is responsible for:
         Rastering, arbitrary beam angle
	 Sampling along the target
	 Pre-vertex multiple scattering
	 External radiative effects
	 Luminosity calculations
 
     This is implemented in the singleton model

*/

class G4VPhysicalVolume;
class G4Material;
class g4hrsMultScatt;
class g4hrsMaterial;

class g4hrsBeamTarget {
    private: 
	static g4hrsBeamTarget *gSingleton;
	g4hrsBeamTarget();


	std::vector <G4VPhysicalVolume *> fUpstreamVols;
	std::vector <G4VPhysicalVolume *> fDownstreamVols;
	std::vector <G4VPhysicalVolume *> fAllVols;
	G4VPhysicalVolume *fTargVol;
	
      	g4hrsMaterial *mMaterialManager;


    public:
	static g4hrsBeamTarget *GetBeamTarget();
	~g4hrsBeamTarget();

	G4double GetEffLumin();

	void Reset(){ fUpstreamVols.clear(); fDownstreamVols.clear(); fTargVol = NULL; UpdateInfo(); }
	void AddUpstreamWall( G4VPhysicalVolume *v ){ fUpstreamVols.push_back(v);  UpdateInfo(); }
	void AddDownstreamWall( G4VPhysicalVolume *v ){ fDownstreamVols.push_back(v);  UpdateInfo(); }
	void SetTargetVolume( G4VPhysicalVolume *v ){ fTargVol = v;  UpdateInfo(); }
	void SetTargetPos(G4double z);
	void SetTargetLen(G4double l);
	void SetTargetMaterial(G4String);
	G4String GetTargetMaterial(); 
	void SetBeamCurrent(G4double I){ fBeamCurr = I; }

	g4hrsVertex SampleVertex(SampType_t);

	G4double fBeamE;
	G4double fBeamCurr;
	G4double fBeamPol;
	G4String fTargetMaterial;
	
	G4VPhysicalVolume * GetTargVol(){ return fTargVol; }

	g4hrsMultScatt *fMS;
    private:
	void UpdateInfo();

	G4double fTargLength;
	G4double fTotalLength;

	G4double fTargRad;
	G4Material *fDefaultMat;

	bool fAlreadyWarned;

	
    public:
	// Base position, angle *sampled* info
	G4ThreeVector fVer, fDir;
	G4double fSampE, fRadLen, fSampLen;
	G4double fTravLen;
	G4double fEcut, fEffMatLen;

	// Base position/angle sampling info
        G4bool fOldRaster;
	G4double fRasterX, fRasterY;
	G4double fX0, fY0;
	G4double fTh0, fPh0;
	G4double fdTh, fdPh, fCorrTh, fCorrPh;

        //ionization loss
        G4double fIonLoss;
};


#endif//__REMOLLBEAMTARGET_HH

