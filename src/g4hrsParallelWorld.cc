#include "g4hrsParallelWorld.hh"
#include "g4hrsGenericDetector.hh"
#include "g4hrsBeamTarget.hh"
#include "g4hrsRun.hh"
#include "g4hrsRunData.hh"
#include "g4hrsIO.hh"
#include "g4hrsMaterial.hh"
#include "g4hrsEMFieldSetup.hh"
#include "g4hrsTune.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"

#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "globals.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"

#include "G4ios.hh"

#include "G4UnitsTable.hh"

#include "G4Polycone.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"

//visual
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Box.hh"

#define __DET_STRLEN 200
#define __MAX_DETS 5000

g4hrsParallelWorld::g4hrsParallelWorld(G4String parallelWorldName) :G4VUserParallelWorld(parallelWorldName) {

    // Default geometry file

	fTune = g4hrsTune::GetTune();
    fHRSAngle=fTune->HRSAngle;
    fSeptumAngle=fTune->septumAngle;

    fTargetW = 2.0*2.54*cm;
    fTargetH = 2.0*2.54*cm;
    fTargetL = 5*mm;

    fTargetX = 0.0;
    fTargetY = 0.0;
    fTargetZ = 0.0;

    fSeptum_UseUniformB=0;

    fPivot2SieveFace = 800*mm;
    fPivotXOffset =  0.0*deg;
    fPivotYOffset =  0.0*deg;
    fPivotZOffset =  1151.2*mm;  // RR new target position (-10 cm)

    fSetupSieveSlit = false;

    fUseSeptumPlusStdHRS=0;
    fSnakeModel = 49;

    fSetupStdScatChamber = 0;

    //#SetupHRS have the following candidates: 
    //# 0: do nothing; 1: will build septum, sieve and VB; 
    //# 2: add Q1; 3: add Q2 ; 4: add Dipole and Q3  
    fSetupHRS = 4;


    fScatChamberRin   =  484.025*mm;
    fScatChamberRout  =  509.425*mm;
    fScatChamberL     =  692.15*mm;

    fSetupCREXGeometry = true;

    //fEMFieldSetup = NULL;

    fFieldX = 45000*mm;
    fFieldY = 45000*mm;
    fFieldZ = 45000*mm;

}

g4hrsParallelWorld::~g4hrsParallelWorld() {
    }

void g4hrsParallelWorld::Construct() {
    
	G4VPhysicalVolume *worldVolume = GetWorld();
	ConstructSD(worldVolume->GetLogicalVolume());

	return;
}



void g4hrsParallelWorld::ConstructSD(G4LogicalVolume* motherLogical) {

    // Add detectors here
/*	
	G4SDManager *SDMan = G4SDManager::GetSDMpointer();
	G4VSensitiveDetector* virtual45 = new g4hrsGenericDetector("virtual45",45);
	SDMan->AddNewDetector(virtual45);

	G4VSolid* parallelDetSolid = new G4Box("parallelDetSolid", 500., 500., 0.05);
        G4LogicalVolume* parallelDetLogical = new G4LogicalVolume(parallelDetSolid,NULL,"parallelDetLogical",0,0,0);		
	parallelDetLogical->SetSensitiveDetector(virtual45);
	G4VPhysicalVolume* parallelDetPhys = new G4PVPlacement(0,G4ThreeVector(fTargetX,fTargetY,fTargetZ + 69.99937 - 74./2. + fPivotZOffset - 10.*cm),
		parallelDetLogical,"parallelDetPhys",motherLogical,0,0);
*/

	double pSeptumBoxX      = 130.0  * cm;
	double pSeptumBoxY      = 80.   * cm;
	double pSeptumZ      = 74.0   * cm;
	
	///////////////////////////////////
	// Solids for virtual boundaries //
	///////////////////////////////////
	G4double vb_thickness = 1. * mm;
	G4VSolid* FPSolid = new G4Tubs("FPTub",0,75.*cm,vb_thickness,0.0*deg,360.0*deg);
	G4VSolid* SeptumSolid = new G4Box("SeptumPlane",pSeptumBoxX/2.,pSeptumBoxY/2.,vb_thickness/2.);
	G4VSolid* PlaneSolid1 = new G4Tubs("PlaneTub",0,35.*cm,vb_thickness,0.0*deg,360.0*deg); //circles
	G4VSolid* PlaneSolid2 = new G4Tubs("PlaneTub",0,75.*cm,vb_thickness,0.0*deg,360.0*deg); //circles
	G4VSolid* CollimatorSolid = new G4Tubs("ColTub",0,20.*cm,vb_thickness,0.0*deg,360.0*deg); //circles
	G4VSolid* LocalAxis = new G4Cons("LocalAxis",20.*cm, 30.*cm, 0.0*cm, 0.1*cm, 40.*cm, 5.0*deg, 350*deg);


	G4LogicalVolume* LPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,0,"LPlaneLogical1",0,0,0);
	G4LogicalVolume* RPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,0,"RPlaneLogical1",0,0,0);
	G4LogicalVolume* LPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,0,"LPlaneLogical1",0,0,0);
	G4LogicalVolume* SepPlaneLogical = new G4LogicalVolume(SeptumSolid,0,"SepPlaneLogical",0,0,0);
	G4LogicalVolume* RPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,0,"RPlaneLogical1",0,0,0);
	G4LogicalVolume* CollimatorLogical = new G4LogicalVolume(CollimatorSolid,0,"CollimatorLogical",0,0,0);
	G4LogicalVolume* LocalAxisLog = new G4LogicalVolume(LocalAxis, 0, "LocalAxisLog",0,0,0);

      
        G4LogicalVolume* LFPLogical = new G4LogicalVolume(FPSolid,0,"LFPLogical",0,0,0);
        G4LogicalVolume* RFPLogical = new G4LogicalVolume(FPSolid,0,"RFPLogical",0,0,0);

	G4VisAttributes* virtualBoundaryVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));	


	LPlaneLogical1->SetVisAttributes(virtualBoundaryVisAtt);	
	LPlaneLogical2->SetVisAttributes(virtualBoundaryVisAtt);
	RPlaneLogical1->SetVisAttributes(virtualBoundaryVisAtt);	
	RPlaneLogical2->SetVisAttributes(virtualBoundaryVisAtt);
	SepPlaneLogical->SetVisAttributes(virtualBoundaryVisAtt);
	LFPLogical->SetVisAttributes(virtualBoundaryVisAtt);
	RFPLogical->SetVisAttributes(virtualBoundaryVisAtt);
	CollimatorLogical->SetVisAttributes(virtualBoundaryVisAtt);
	LocalAxisLog->SetVisAttributes(virtualBoundaryVisAtt);


	///////////////////////////////
	// Septum virtual boundaries //
	///////////////////////////////
	
	// Septum virtual boundaries will be placed directly in the (parallel) world volume
	
	//double pSeptumPlaceZ = 70.414 * cm;
	double pSeptumPlaceZ = 70.0 * cm;
	new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ - 0.5 * pSeptumZ + fPivotZOffset),
		SepPlaneLogical,"virtualBoundaryPhys_sen",motherLogical,0,0);//sen
	new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ + fPivotZOffset),
		SepPlaneLogical,"virtualBoundaryPhys_sm",motherLogical,0,0);//sm
	new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ + 0.5 * pSeptumZ + fPivotZOffset),
	    	SepPlaneLogical,"virtualBoundaryPhys_sex",motherLogical,0,0);//sex

	///////////////////////////////
	// Septum tracking detectors //
	///////////////////////////////
	
	// Septum tracking detectors will be placed directly in the (parallel) world volume
	// at engineering-critical z locations ("pinch points")

	// pinch points are known from CAD in HCS, so the pivot offset is required
	double zpinch1 =  -7.8*cm + fPivotZOffset;
	double zpinch2 =  21.3*cm + fPivotZOffset;
	double zpinch3 = 109.5*cm + fPivotZOffset;

	// other tracking detectors in G4 coordinates
	double ztarg = 1.*cm;
	double zfield0 = 99.378*cm;
	double zfield1 = 100.378*cm;
	double zfield2 = 101.378*cm;
	double zmidtosep = 69.189*cm;

	// tracking in the septum	
	double zsep1 = (pSeptumPlaceZ - 0.5 * pSeptumZ + fPivotZOffset) + 0.1 * pSeptumZ; 	
	double zsep2 = (pSeptumPlaceZ - 0.5 * pSeptumZ + fPivotZOffset) + 0.2 * pSeptumZ; 	
	double zsep3 = (pSeptumPlaceZ - 0.5 * pSeptumZ + fPivotZOffset) + 0.3 * pSeptumZ; 	
	double zsep4 = (pSeptumPlaceZ - 0.5 * pSeptumZ + fPivotZOffset) + 0.4 * pSeptumZ; 	

         //tracking the full path of the acceptance
         //RR upstream pivot
         double zup[2] = { fPivotZOffset-506*mm ,fPivotZOffset-459*mm };
         //RR downstream pivot
         double zdown[9] = { fPivotZOffset + 7*mm,  fPivotZOffset + 65*mm ,  fPivotZOffset + 236*mm,  fPivotZOffset + 353*mm,
         fPivotZOffset + 471*mm,  fPivotZOffset + 510*mm,  fPivotZOffset + 590*mm,  fPivotZOffset + 1012*mm,  fPivotZOffset + 1114*mm };

         //RR -- adding the sieve boundary
         double sieve = fPivotZOffset - 155.6*mm; 
         new G4PVPlacement(0,G4ThreeVector(0,0,sieve),SepPlaneLogical,"virtualBoundaryPhys_zsieve",motherLogical,0,0);



	new G4PVPlacement(0,G4ThreeVector(0,0,zpinch1),
		SepPlaneLogical,"virtualBoundaryPhys_zpinch1",motherLogical,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,zpinch2),
		SepPlaneLogical,"virtualBoundaryPhys_zpinch2",motherLogical,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,zpinch3),
		SepPlaneLogical,"virtualBoundaryPhys_zpinch3",motherLogical,0,0);

	new G4PVPlacement(0,G4ThreeVector(0,0,ztarg),
		SepPlaneLogical,"virtualBoundaryPhys_ztarg",motherLogical,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,zfield0),
		SepPlaneLogical,"virtualBoundaryPhys_zfield0",motherLogical,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,zfield1),
		SepPlaneLogical,"virtualBoundaryPhys_zfield1",motherLogical,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,zfield2),
		SepPlaneLogical,"virtualBoundaryPhys_zfield2",motherLogical,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,zmidtosep),
		SepPlaneLogical,"virtualBoundaryPhys_zmidtosep",motherLogical,0,0);

	new G4PVPlacement(0,G4ThreeVector(0,0,zsep1),
		SepPlaneLogical,"virtualBoundaryPhys_zsep1",motherLogical,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,zsep2),
		SepPlaneLogical,"virtualBoundaryPhys_zsep2",motherLogical,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,zsep3),
		SepPlaneLogical,"virtualBoundaryPhys_zsep3",motherLogical,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,zsep4),
		SepPlaneLogical,"virtualBoundaryPhys_zsep4",motherLogical,0,0);

        //RR adding instances of these volumes 
        new G4PVPlacement(0,G4ThreeVector(0,0,zup[0]),
               SepPlaneLogical,"virtualBoundaryPhys_zup1",motherLogical,0,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,zup[1]),
            SepPlaneLogical,"virtualBoundaryPhys_zup2",motherLogical,0,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,zdown[0]),
               SepPlaneLogical,"virtualBoundaryPhys_zdown1",motherLogical,0,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,zdown[1]),
               SepPlaneLogical,"virtualBoundaryPhys_zdown2",motherLogical,0,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,zdown[2]),
               SepPlaneLogical,"virtualBoundaryPhys_zdown3",motherLogical,0,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,zdown[3]),
               SepPlaneLogical,"virtualBoundaryPhys_zdown4",motherLogical,0,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,zdown[4]),
               SepPlaneLogical,"virtualBoundaryPhys_zdown5",motherLogical,0,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,zdown[5]),
               SepPlaneLogical,"virtualBoundaryPhys_zdown6",motherLogical,0,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,zdown[6]),
               SepPlaneLogical,"virtualBoundaryPhys_zdown7",motherLogical,0,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,zdown[7]),
               SepPlaneLogical,"virtualBoundaryPhys_zdown8",motherLogical,0,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,zdown[8]),
               SepPlaneLogical,"virtualBoundaryPhys_zdown9",motherLogical,0,0);



	////////////////////////////
	// HRS virtual boundaries //
	////////////////////////////
	
	// LHRS (RHRS) virutal boundaries will be placed in a logical LHRS container (RHRS container) which will be placed in the world volume
	// Note that the HRS containers are created along -y axis, rotated to along the z axis, then rotated +/- HRS angle
	// Therefore the placement of virtual boundaries in the HRS containers will not be with respect to the hall axes 


	////////////////////
	// HRS containers //
	////////////////////

	double pHRSContainerRin=1.35*m,pHRSContainerRout=25*m;//NickieMode it just has to be big enough to fit my VB at FP
	const int kNPlane_HRSContainer=7;
	double rInner_HRSContainer[] = {	pHRSContainerRin,
						pHRSContainerRin,
						2.5*m,
						3.7*m,
						9.0*m,
						pHRSContainerRout-3.0*m,
						pHRSContainerRout};
	
	double rOuter_HRSContainer[] = {	pHRSContainerRout,
						pHRSContainerRout,
						pHRSContainerRout,
						pHRSContainerRout,
						pHRSContainerRout,
						pHRSContainerRout,
						pHRSContainerRout};
	
	double zPlane_HRSContainer[] = {	-2.0*m,
						1.0*m,
						1.0*m,
						2.0*m,
						2.0*m,
						11.0*m,
						11.0*m};
	
	
	G4Polycone* HRSContainerSolid = new G4Polycone("HRSContainer",257.625*deg,24.75*deg,
		kNPlane_HRSContainer,zPlane_HRSContainer,rInner_HRSContainer,rOuter_HRSContainer);

	G4LogicalVolume* LHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		0,"LHRSContainerLogical",0,0);
	G4LogicalVolume* RHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		0,"RHRSContainerLogical",0,0);

	
	G4VisAttributes *HallVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
	HallVisAtt->SetVisibility(false);

	LHRSContainerLogical->SetVisAttributes(HallVisAtt); 
	RHRSContainerLogical->SetVisAttributes(HallVisAtt);
	
	G4RotationMatrix *pRotLHRSContainer=new G4RotationMatrix();
	pRotLHRSContainer->rotateX(-270*deg);
	pRotLHRSContainer->rotateZ(-fHRSAngle);  

	G4RotationMatrix *pRotRHRSContainer=new G4RotationMatrix();
	pRotRHRSContainer->rotateX(-270*deg);
	pRotRHRSContainer->rotateZ(fHRSAngle);  

	new G4PVPlacement(pRotLHRSContainer,G4ThreeVector(0.,0.,0.+fPivotZOffset),
		LHRSContainerLogical,"LHRSContainerPhys",motherLogical,0,0,0);
	new G4PVPlacement(pRotRHRSContainer,G4ThreeVector(0,0,0.+fPivotZOffset),
		RHRSContainerLogical,"RHRSContainerPhys",motherLogical,0,0,0);

	
	G4RotationMatrix* LHRSRot = new G4RotationMatrix();
	LHRSRot->rotateY(fHRSAngle); 
	G4RotationMatrix* RHRSRot = new G4RotationMatrix();
	RHRSRot->rotateY(-fHRSAngle); 
	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 
	pRotX90deg->rotateZ(90*deg);
	G4RotationMatrix *pRotX45deg=new G4RotationMatrix();
	pRotX45deg->rotateX(-45*deg);
	pRotX45deg->rotateZ(90*deg); 
	G4RotationMatrix *pRotVDC = new G4RotationMatrix();
	pRotVDC->rotateX(-45*deg);// RR this to be in transport coordinates
        pRotVDC->rotateZ(90*deg);
	G4RotationMatrix *pRotFP = new G4RotationMatrix();
	pRotFP->rotateX(-45*deg);
	pRotFP->rotateZ(90*deg);

//	new G4PVPlacement(0,G4ThreeVector(0,0,0),LocalAxisLog,"localaxis_lhrs",motherLogical,0,0,0);

	// Positions from SNAKE
	double pivotOffset = 115.21*cm;
	double col = 1.38*m;
        //Updating these positions to place the parallel planes at the entrance and exit of the quads
	double q1en = 172.*cm;
    	double q1ex = q1en + 70.*cm;
	double q2en = q1ex + 127.71*cm;
	double q2ex = q2en + 182.66*cm;
	double dpen = q2ex + 443.74*cm;
	double dpex_z = 15.9006973*m;
	double dpex_x = 2.4603032*m;
	double q3en_z = 17.0257042*m;
	double q3en_x = 3.5853101*m;
	double q3l = 1.8268*m;
	double q3ex_z = q3en_z + (q3l / sqrt(2.));
	double q3ex_x = q3en_x + (q3l / sqrt(2.)); 
	double vdc_drift = 3.4538*m;
	double vdc_z = q3ex_z + (vdc_drift / sqrt(2.));
	double vdc_x = q3ex_x + (vdc_drift / sqrt(2.));
	double fp_drift = 0.9*m;
	double fp_z = vdc_z + (fp_drift / sqrt(2.));
	double fp_x = vdc_x + (fp_drift / sqrt(2.));

	///////////////////////////////////
	// Collimator virtual boundaries //	
	///////////////////////////////////
	//LHRS
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-col,0),
		CollimatorLogical,"virtualBoundaryPhys_col_LHRS",LHRSContainerLogical,0,0,0);//col
	//RHRS
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-col,0),
		CollimatorLogical,"virtualBoundaryPhys_col_RHRS",RHRSContainerLogical,0,0,0);//col
	
	///////////////////////////
	// Q1 virtual boundaries //
	///////////////////////////
	//LHRS
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-q1en,0),
		LPlaneLogical1,"virtualBoundaryPhys_q1en_LHRS",LHRSContainerLogical,0,0,0);//q1en	
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-q1ex,0),
		LPlaneLogical1,"virtualBoundaryPhys_q1ex_LHRS",LHRSContainerLogical,0,0,0);//q1ex
	//RHRS
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-q1en,0),
		RPlaneLogical1,"virtualBoundaryPhys_q1en_RHRS",RHRSContainerLogical,0,0,0);//q1en	
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-q1ex,0),
		RPlaneLogical1,"virtualBoundaryPhys_q1ex_RHRS",RHRSContainerLogical,0,0,0);//q1ex

	///////////////////////////
	// Q2 virtual boundaries //
	///////////////////////////
	//LHRS
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-q2en,0),		    
		LPlaneLogical2,"virtualBoundaryPhys_q2en_LHRS",LHRSContainerLogical,0,0,0);//q2en
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-q2ex,0),
		LPlaneLogical2,"virtualBoundaryPhys_q2ex_LHRS",LHRSContainerLogical,0,0,0);//q2ex
	//RHRS
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-q2en,0),
		RPlaneLogical2,"virtualBoundaryPhys_q2en_RHRS",RHRSContainerLogical,0,0,0);//q2en
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-q2ex,0),
		RPlaneLogical2,"virtualBoundaryPhys_q2ex_RHRS",RHRSContainerLogical,0,0,0);//q2ex

	//////////////////////////
	// D virtual boundaries //
	//////////////////////////
	//LHRS
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-dpen,0),
		LPlaneLogical2,"virtualBoundaryPhys_den_LHRS",LHRSContainerLogical,0,0,0);//den
	new G4PVPlacement(pRotX45deg,G4ThreeVector(0, -dpex_z,dpex_x),
		LPlaneLogical2,"virtualBoundaryPhys_dex_LHRS",LHRSContainerLogical,0,0,0);//dex
	//RHRS
	new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-dpen,0),
		RPlaneLogical2,"virtualBoundaryPhys_den_RHRS",RHRSContainerLogical,0,0,0);//den
	new G4PVPlacement(pRotX45deg,G4ThreeVector(0, -dpex_z,dpex_x),
		RPlaneLogical2,"virtualBoundaryPhys_dex_RHRS",RHRSContainerLogical,0,0,0);//dex
	  
  	///////////////////////////
	// Q3 virtual boundaries //
	///////////////////////////
	//LHRS
	new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-q3en_z,q3en_x),
		LPlaneLogical2,"virtualBoundaryPhys_q3en_LHRS",LHRSContainerLogical,0,0,0);//q3en
	new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-q3ex_z,q3ex_x),
		LPlaneLogical2,"virtualBoundaryPhys_q3ex_LHRS",LHRSContainerLogical,0,0,0);//q3ex
	//RHRS
	new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-q3en_z,q3en_x),
		RPlaneLogical2,"virtualBoundaryPhys_q3en_RHRS",RHRSContainerLogical,0,0,0);//q3en
	new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-q3ex_z,q3ex_x),
		RPlaneLogical2,"virtualBoundaryPhys_q3ex_RHRS",RHRSContainerLogical,0,0,0);//q3ex

	////////////////////////////
	// VDC virtual boundaries //	
	////////////////////////////
	//LHRS
	new G4PVPlacement(pRotVDC,G4ThreeVector(0,-vdc_z,vdc_x),
		LFPLogical,"virtualBoundaryPhys_vdc_LHRS",LHRSContainerLogical,0,0,0);//vdc	
	//RHRS
	new G4PVPlacement(pRotVDC,G4ThreeVector(0,-vdc_z,vdc_x),
		RFPLogical,"virtualBoundaryPhys_vdc_RHRS",RHRSContainerLogical,0,0,0);//vdc	

	///////////////////////////
	// FP virtual boundaries //	
	///////////////////////////
	//LHRS
	new G4PVPlacement(pRotFP,G4ThreeVector(0,-fp_z,fp_x),
		LFPLogical,"virtualBoundaryPhys_fp_LHRS",LHRSContainerLogical,0,0,0);//vdc
	//RHRS	
	new G4PVPlacement(pRotFP,G4ThreeVector(0,-fp_z,fp_x),
		RFPLogical,"virtualBoundaryPhys_fp_RHRS",RHRSContainerLogical,0,0,0);//vdc
	
	return;

}





















