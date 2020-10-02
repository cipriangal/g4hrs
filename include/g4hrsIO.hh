#ifndef g4hrsIO_HH
#define g4hrsIO_HH

#include "TROOT.h"
#include "TObject.h"
#include "G4Run.hh"
#include "g4hrstypes.hh"

#include "G4String.hh"

class TFile;
class TTree;

class g4hrsGenericDetectorHit;
class g4hrsGenericDetectorSum;
class g4hrsEvent;
class g4hrsSteppingAction;

#include <xercesc/dom/DOMElement.hpp>


#define __IO_MAXHIT 10000
#define __FILENAMELEN 255

// Units for output
#define __E_UNIT GeV
#define __L_UNIT m
#define __T_UNIT ns
#define __ANG_UNIT rad
#define __ASYMM_SCALE 1e-9 // ppb

class g4hrsIO {
    public:
	 g4hrsIO();
	~g4hrsIO();

	void SetFilename(G4String  fn);
	G4String GetFilename(){return fFilename;}

	void FillTree();
	void Flush();
	void WriteTree();

	void WriteRun();

	void InitializeTree();

	void GrabGDMLFiles( G4String fn );

    private:
	TFile *fFile;
	TTree *fTree;

	char fFilename[__FILENAMELEN];

	std::vector<G4String>       fGDMLFileNames;

	void SearchGDMLforFiles(G4String );
	void TraverseChildren(  xercesc::DOMElement * );

	//  Interfaces and buffers to the tree
	//  This is potentially going to get very long...

	// Event data
    public:
	void SetEventData(g4hrsEvent *);
    private:
	Int_t fNEvPart;

	Double_t fEvRate;
	Double_t fEvEffXS;
	Double_t fEvAsym;
	Double_t fEvmAsym;
	Double_t fEvSens;
	Double_t fEvWeight;
	Double_t fEvBeamP;
	Double_t fEvQ2;
	Double_t fEvW2;
	Double_t fEvThCoM;

	Double_t fBmX;
	Double_t fBmY;
	Double_t fBmZ;
	Double_t fBmdX;
	Double_t fBmdY;
	Double_t fBmdZ;
	Double_t fBmth;
	Double_t fBmph;

	Int_t fEvPID[__IO_MAXHIT];

	Double_t fEvPart_X[__IO_MAXHIT];
	Double_t fEvPart_Y[__IO_MAXHIT];
	Double_t fEvPart_Z[__IO_MAXHIT];
	Double_t fEvPart_Px[__IO_MAXHIT];
	Double_t fEvPart_Py[__IO_MAXHIT];
	Double_t fEvPart_Pz[__IO_MAXHIT];
	Double_t fEvPart_Th[__IO_MAXHIT];
	Double_t fEvPart_Ph[__IO_MAXHIT];
	Double_t fEvPart_P[__IO_MAXHIT];
	Double_t fEvPart_tPx[__IO_MAXHIT];
	Double_t fEvPart_tPy[__IO_MAXHIT];
	Double_t fEvPart_tPz[__IO_MAXHIT];


	//  GenericDetectorHit
    public:
	void AddGenericDetectorHit(g4hrsGenericDetectorHit *);
    private:
	Int_t fNGenDetHit;
	Int_t fGenDetHit_det[__IO_MAXHIT];
	Int_t fGenDetHit_id[__IO_MAXHIT];

	Int_t fGenDetHit_trid[__IO_MAXHIT];
	Int_t fGenDetHit_pid[__IO_MAXHIT];
	Int_t fGenDetHit_gen[__IO_MAXHIT];
	Int_t fGenDetHit_mtrid[__IO_MAXHIT];

	Double_t fGenDetHit_X[__IO_MAXHIT];
	Double_t fGenDetHit_Y[__IO_MAXHIT];
	Double_t fGenDetHit_Z[__IO_MAXHIT];
	Double_t fGenDetHit_R[__IO_MAXHIT];
	Double_t fGenDetHit_Ph[__IO_MAXHIT];

	Double_t fGenDetHit_Px[__IO_MAXHIT];
	Double_t fGenDetHit_Py[__IO_MAXHIT];
	Double_t fGenDetHit_Pz[__IO_MAXHIT];
	Double_t fGenDetHit_P[__IO_MAXHIT];
	Double_t fGenDetHit_E[__IO_MAXHIT];
	Double_t fGenDetHit_M[__IO_MAXHIT];

	Double_t fGenDetHit_Vx[__IO_MAXHIT];
	Double_t fGenDetHit_Vy[__IO_MAXHIT];
	Double_t fGenDetHit_Vz[__IO_MAXHIT];

	Int_t fCollCut;

	//  GenericDetectorSum
    public:
	void AddGenericDetectorSum(g4hrsGenericDetectorSum *);
    private:
	Int_t fNGenDetSum;
	Int_t fGenDetSum_det[__IO_MAXHIT];
	Int_t fGenDetSum_id[__IO_MAXHIT];
	Double_t fGenDetSum_edep[__IO_MAXHIT];

	// Virtual boundary data 
    public:
	void SetSteppingAction(g4hrsSteppingAction *stepper) { fSteppingAction = stepper;}
	void ClearVirtualBoundaryData();
	void SetVirtualBoundaryData();
    private:

	g4hrsSteppingAction* fSteppingAction;

	G4int fLHRS;
	G4int fRHRS;

	G4double fX0;
	G4double fY0;
	G4double fZ0;
	G4double fP0;
	G4double fTh0;
	G4double fPh0;
	G4double fX0_tr;
	G4double fY0_tr;
	G4double fZ0_tr;
	G4double fP0_tr;
	G4double fTh0_tr;
	G4double fPh0_tr;

	int numTF;
	int numTFvar;
	G4double TFdata[4][12];
	char TFName[12][6];
	char TFVarName[4][5];

	int numVB;
	int numVar;
	G4double VBdata[14][12];
	char VBName[14][6];
	char VarName[12][5];

	int numZCrit;
	int numZCritVar;
	G4double ZCritData[24][5];
	char ZCritName[24][10];
	char ZCritVarName[5][5];

};

#endif//g4hrsIO_HH
