#include "g4hrsDetectorConstruction.hh"
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


//visual
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#define __DET_STRLEN 200
#define __MAX_DETS 5000

g4hrsDetectorConstruction::g4hrsDetectorConstruction() {

    // Default geometry file

    fWorldVolume = NULL;

    fHRSAngle=12.5*deg;
    fSeptumAngle=5.0*deg;

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
    fPivotZOffset =  1151.2*mm;   //RR new target position (-10 cm)


    fSetupSieveSlit = false;

    fUseSeptumPlusStdHRS=0;
    fSnakeModel = 49;
	fTune = g4hrsTune::GetTune();
    fSetupStdScatChamber = 0;

	// create septum and HRS fields
//    	fEMFieldSetup = new g4hrsEMFieldSetup();
    
    //#SetupHRS have the following candidates: 
    //# 0: do nothing; 1: will build septum, sieve and VB; 
    //# 2: add Q1; 3: add Q2 ; 4: add Dipole and Q3  
    fSetupHRS = 4;


    fScatChamberRin   =  484.025*mm;
    fScatChamberRout  =  509.425*mm;
    fScatChamberL     =  692.15*mm;

    fScatChamberXOffset = 0.0;
    fScatChamberYOffset = 0.0;
    fScatChamberZOffset = -fPivotZOffset;
    fSetupCREXGeometry = true;

    //fEMFieldSetup = NULL;

    fFieldX = 45000*mm;
    fFieldY = 45000*mm;
    fFieldZ = 45000*mm;

}

g4hrsDetectorConstruction::~g4hrsDetectorConstruction() {
        if(fEMFieldSetup){ 
            delete fEMFieldSetup;
            fEMFieldSetup = NULL;
        }

}

G4VPhysicalVolume* g4hrsDetectorConstruction::Construct() {
    
    G4VPhysicalVolume *worldVolume;

    fEMFieldSetup = new g4hrsEMFieldSetup();
 
	//    g4hrsMaterial::g4hrsMaterial();
    mMaterialManager = g4hrsMaterial::GetHRSMaterialManager();

    
    int z, nelements;
    double a, density;

    // Air
    G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
    G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

    G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
    Air->AddElement(N, 70.*perCent);
    Air->AddElement(O, 30.*perCent);

    G4VSolid *worldSolid = new G4Box("worldbox", 50.0*m, 50.0*m, 50.0*m ); 
//    G4LogicalVolume *worldLog = new G4LogicalVolume(worldSolid, Air, "worldLog", 0, 0, 0);
    G4LogicalVolume *worldLog = new G4LogicalVolume(worldSolid, mMaterialManager->vacuum, "worldLog", 0, 0, 0);

    fWorldVolume = new G4PVPlacement(0, G4ThreeVector(), worldLog, "world", 0, false, 0);
    
    CreateTarget(worldLog);

    G4NistManager *nist = G4NistManager::Instance();
    G4VSolid *acceptLimiter_sv = new G4Box("acceptLimiter_solid", 130*cm/2, 80*cm/2, 0.1*mm);
    G4LogicalVolume *acceptLimiter_lv = new G4LogicalVolume(acceptLimiter_sv, 
							    nist->FindOrBuildMaterial("G4_Galactic"),
							    "acceptLimiter_lv");
    new G4PVPlacement(0,G4ThreeVector(0,0,10*cm),
		      acceptLimiter_lv,
		      "acceptLimiter",
		      worldLog,
		      false,0,true);

    CreateSeptum(worldLog);
//    CreateTargetChamber(worldLog);
    CreateHRS(worldLog);

    return fWorldVolume;
}

void g4hrsDetectorConstruction::CreateTarget(G4LogicalVolume *pMotherLogVol){

    g4hrsBeamTarget *beamtarg = g4hrsBeamTarget::GetBeamTarget();
    beamtarg->Reset();

	// Make lead the default target
	// Can be changed in macro
	G4Material* targ_material = mMaterialManager->lead208;

    G4VSolid* targetSolid  = new G4Tubs("targetBox", 0.0, fTargetW, fTargetL / 2.0, 0, 360*deg );
    G4LogicalVolume* targetLogical = new G4LogicalVolume(targetSolid,targ_material,"targetLogical",0,0,0);

    G4VPhysicalVolume *phystarg = new G4PVPlacement(0,G4ThreeVector(fTargetX, fTargetY, fTargetZ),
        targetLogical,"targetPhys",pMotherLogVol,0,0);
    
    beamtarg->SetTargetVolume(phystarg);

    return;

}

void g4hrsDetectorConstruction::CreateSeptum(G4LogicalVolume *pMotherLogVol){
	const double inch=2.54*cm;
	double startphi,deltaphi;
	G4String SDname;

	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4VSensitiveDetector* sieveSlitSD=new g4hrsGenericDetector(SDname="sieveSlit", 1);
	G4VSensitiveDetector* septumWindowSD=new g4hrsGenericDetector(SDname="septumWindow", 2);

	G4SDManager* SDman = G4SDManager::GetSDMpointer();

        //g4hrsMaterial* mMaterialManager = g4hrsMaterial::GetHRSMaterialManager();


	/////////////////////////////////////////////////
	//From Hall A NIM Paper, the standard sieve slit
	//Each spectrometer is equipped with a set of collimators, positioned 1:109 +/- 0.005 
	//and 1:101+/-0.005 m from the target on the left and right spectrometers, respectively. 
	//There is a large collimator, made of 80 mm thick tungsten, with a 121.8 mm vertical 
	//opening and a 62:9 mm horizontal opening at the entrance face. The opening in this 
	//collimator expands to 129.7 by 66.8 mm at the exit face. A second smaller collimator, 
	//made of the same material, is 50.0 by 21.3 mm at the entrance face and 53.2 by 22:6 mm 
	//at the exit face. The third collimator is the sieve slit, which is used to study the 
	//optical properties of the spectro- meters. The sieve is a 5 mm thick stainless steel 
	//sheet with a pattern of 49 holes 7 x 7, spaced 25 mm apart vertically and 12:5 mm apart 
	//horizontally. Two of the holes, one in the center and one displaced two rows vertically 
	//and one horizontally, are 4 mm in diameter. The rest of the holes are 2 mm in diameter. 
	//The sieve slits are positioned 75 mm further from the target than the other collimators.

       
	G4RotationMatrix *pRotLHRS=new G4RotationMatrix();
	pRotLHRS->rotateY(fHRSAngle); 
	G4RotationMatrix *pRotRHRS=new G4RotationMatrix();
	pRotRHRS->rotateY(-fHRSAngle); 

	G4RotationMatrix *pRotLSeptum=new G4RotationMatrix();
	pRotLSeptum->rotateY(fSeptumAngle); 
	G4RotationMatrix *pRotRSeptum=new G4RotationMatrix();
	pRotRSeptum->rotateY(-fSeptumAngle); 
	G4RotationMatrix *pRotNone=new G4RotationMatrix();
	pRotNone->rotateY(0); 

	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 
	/////////////////////////////////////////////////
        //
	///////////////////////////////////////
	//Sieve Slit for HRS-Angle=5.65
	///////////////////////////////////////

	////////////////the following is for 6 degree sieve//////////////

	//these 2 will be read from ini file, it is 80.0cm for 6 deg and 120 cm for 12.5 deg

	double pSieveSlitX=2.205*inch; //33.13*mm
	double pSieveSlitY=5.134*inch; //130.40*mm
	double pSieveSlitZ=0.2*inch;

	double pSieveSlitHoleR=0.6985*mm;           //radius of small hole 0.055/2 inch
	double pSieveSlitLargeHoleR=1.3462*mm;      //radius of large hole 0.106/2 inch
	double pSieveSlitHoldL=pSieveSlitZ+0.1*mm;  //need to make it longer to avoid round off in the subtraction

	//the big center hole horizontal and vertical offset, From servey
	double pSieveSlitLargeHoleH=0*mm;    //positive means shift away from the beam line
	double pSieveSlitLargeHoleV=0*mm;

	//please note that the following constants are for the sieve in the right arm only
	//need to mirror(flip) it in order to put in the left arm  
	double pSieveSlitDeltaH[8]={0.537*inch, 0.188*inch, 0.188*inch, 0.188*inch,
		0.241*inch, 0.241*inch, 0.241*inch, 0.381*inch}; //in inch, from left to right
	double pSieveSlitDeltaV[8]={0.496*inch, 0.524*inch, 0.524*inch, 0.524*inch,
		0.524*inch, 0.524*inch, 0.524*inch, 1.494*inch}; //in inch, from top to buttom

	//the whole position relative to the slit center 
	double pSieveSlitHolePosH[7], pSieveSlitHolePosV[7];
	for(int ii=0;ii<7;ii++)
	{
		pSieveSlitHolePosH[ii] = (ii==0)?pSieveSlitX/2.0:pSieveSlitHolePosH[ii-1];
		pSieveSlitHolePosH[ii] -= pSieveSlitDeltaH[ii];

		pSieveSlitHolePosV[ii] = (ii==0)?pSieveSlitY/2.0:pSieveSlitHolePosV[ii-1];
		pSieveSlitHolePosV[ii] -= pSieveSlitDeltaV[ii];
	}


	//now start to build box then subtract 49 holes 
	G4VSolid* sieveSlitWholeSolid=new G4Box("sieveSlitWholeBox",pSieveSlitX/2.0,
		pSieveSlitY/2.0,pSieveSlitZ/2.0);
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* sieveSlitHoleSolid=new G4Tubs("sieveSlitHoleTubs",0,pSieveSlitHoleR,
		pSieveSlitHoldL/2.0,startphi,deltaphi); 
	G4VSolid* sieveSlitLargeHoleSolid=new G4Tubs("sieveSlitLargeHoleTubs",0,
		pSieveSlitLargeHoleR,pSieveSlitHoldL/2.0,startphi,deltaphi); 
	G4SubtractionSolid* sieveSlitSolid=(G4SubtractionSolid*)sieveSlitWholeSolid;
	char strName[100];
	for(int ih=0;ih<7;ih++)
	{
		for(int iv=0;iv<7;iv++)
		{
			sprintf(strName,"sieveSlitHole_H%d_V%d",ih,iv);
			if((ih==3 && iv==3) || (ih==4 && iv==1)) 
			{
				//now dig large holes in the block
				sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
					sieveSlitLargeHoleSolid,0,
					G4ThreeVector(pSieveSlitHolePosH[ih],pSieveSlitHolePosV[iv],0));
			}
			else
			{
				//now dig small holes in the block
				sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
					sieveSlitHoleSolid,0,
					G4ThreeVector(pSieveSlitHolePosH[ih],pSieveSlitHolePosV[iv],0));
			}
		}
	}
	sieveSlitSolid->SetName("sieveSlitSolid");

        G4VisAttributes *LeadVisAtt = new G4VisAttributes(G4Colour(204/255.,204./255.,255./255.));

	G4LogicalVolume* sieveSlitLogical = new G4LogicalVolume(sieveSlitSolid,
		mMaterialManager->tungsten,"sieveSlitLogical",0,0,0);
	sieveSlitLogical->SetVisAttributes(LeadVisAtt); 

	SDman->AddNewDetector(sieveSlitSD);
	sieveSlitLogical->SetSensitiveDetector(sieveSlitSD);

	//calculate the center position in the Lab frame
	double pSieveSlitCenterHOffset=pSieveSlitLargeHoleH-pSieveSlitHolePosH[3];
	double pSieveSlitCenterVOffset=pSieveSlitLargeHoleV-pSieveSlitHolePosV[3];

	//place the sieve slits in the hall
	double pLSieveSlitPos_X=(fPivot2SieveFace+pSieveSlitZ/2.0)*sin(fSeptumAngle)+
		pSieveSlitCenterHOffset+fPivotXOffset;
	double pLSieveSlitPos_Y=pSieveSlitCenterVOffset+fPivotYOffset;
	double pLSieveSlitPos_Z=(fPivot2SieveFace+pSieveSlitZ/2.0)*cos(fSeptumAngle)+fPivotZOffset;
	double pRSieveSlitPos_X=(fPivot2SieveFace+pSieveSlitZ/2.0)*sin(-fSeptumAngle)-
		pSieveSlitCenterHOffset+fPivotXOffset;
	double pRSieveSlitPos_Y=pSieveSlitCenterVOffset+fPivotYOffset;
	double pRSieveSlitPos_Z=(fPivot2SieveFace+pSieveSlitZ/2.0)*cos(-fSeptumAngle)+fPivotZOffset;

	if(fSetupSieveSlit)
	{ 
		G4RotationMatrix *pRotLSieve=new G4RotationMatrix();
		pRotLSieve->rotateY(-fSeptumAngle-180*deg);
		new G4PVPlacement(pRotLSieve,
		G4ThreeVector(pLSieveSlitPos_X,pLSieveSlitPos_Y,pLSieveSlitPos_Z),
		sieveSlitLogical,"leftSieveSlitPhys",pMotherLogVol,0,0);
	}
	if(fSetupSieveSlit)
	{
	  new G4PVPlacement(pRotRSeptum,
	  G4ThreeVector(pRSieveSlitPos_X,pRSieveSlitPos_Y,pRSieveSlitPos_Z),
	  sieveSlitLogical,"rightSieveSlitPhys",pMotherLogVol,0,0);
	}

	/////////////////////////
	// Septum block 
	/////////////////////////
	//by Jixie: Allow one to setup septum without HRS
	double col_distance = 1.38*m;//or is it 1.39?
        ////////////////////////////////////////////////////////////////////
        //Septum block, 140 cm width, 84.4 cm height and 74 cm in length, silicon steel
        //Tunnel size: started from x=8.4cm, 30.4cm wide and 24.4 cm in height
        //located at z=700 mm, no rotation
        double pSeptumX=140.0*cm;
        double pSeptumY=84.4*cm;
        double pSeptumZ=74.0*cm;
        double pSeptumTunnelX=30.55*cm;
        double pSeptumTunnelY=24.4*cm;//-2.0*inch;  //By Jixie @20120205: Add 2 inches of iron
        double pSeptumBeamHoleX=7.8*cm;
        double pSeptumBeamHoleY=8.0*cm;



        double pSeptumTunnelPos_X=8.425*cm+pSeptumTunnelX/2.0;
        double pSeptumPos_Z=70.0*cm;		




        G4VSolid* septumBlockSolid = new G4Box("septumBlockBox",pSeptumX/2.0,
                pSeptumY/2.0,pSeptumZ/2.0);

        //Left and right tunnels, treat the cu coils as part of the block
        //By Jixie: I reduced this by 0.2cm for the Helium bag
        G4VSolid* septumTunnelSolid = new G4Box("septumTunnelBox",pSeptumTunnelX/2.0-0.2*cm,
                pSeptumTunnelY/2.0-0.2*cm,pSeptumZ/2.0+1.0*mm);

        //beam pine hole
        G4VSolid* septumBeamHoleSolid = new G4Box("septumBeamHoleBox",pSeptumBeamHoleX/2.0,
                pSeptumBeamHoleY/2.0,pSeptumZ/2.0+1.0*mm);

        //dig 3 holes, left, right tunnel and beam hole
        G4SubtractionSolid* septumBlockSubLSolid=new G4SubtractionSolid("septumBlockSubL",
                septumBlockSolid,septumTunnelSolid,0,G4ThreeVector(pSeptumTunnelPos_X,0,0));
        G4SubtractionSolid* septumBlockSubLRSolid=new G4SubtractionSolid("septumBlockSubLR",
                septumBlockSubLSolid,septumTunnelSolid,0,G4ThreeVector(-pSeptumTunnelPos_X,0,0));
        G4SubtractionSolid* septumBlockSubLRCSolid=new G4SubtractionSolid("septumBlockSubLRC",
                septumBlockSubLRSolid,septumBeamHoleSolid);


        G4VisAttributes *IronVisAtt = new G4VisAttributes(G4Colour(100./255.,149./255.,237./255.));

//        G4LogicalVolume* septumLogical = new G4LogicalVolume(septumBlockSubLRSolid,
        G4LogicalVolume* septumLogical = new G4LogicalVolume(septumBlockSubLRCSolid,
                mMaterialManager->siliconsteel,"septumLogical",0,0,0);
        septumLogical->SetVisAttributes(IronVisAtt);

	G4FieldManager* septumFieldManager = fEMFieldSetup->GetFieldManager();
	pMotherLogVol->SetFieldManager(septumFieldManager,true);

        G4int placeseptum = 1;
        //put it in the hall, no rotation
        if( placeseptum ){
          new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPos_Z+fPivotZOffset),
                            septumLogical,"septumPhys",pMotherLogVol,0,0,0);
        }

        double pSeptumCoilRadius = 22.06*cm;
        double pSeptumCoilCenterX = 2.83*cm;
        double pSeptumCoilCenterY = -10.25*cm;
        double pSeptumCoilThickness = 4.5*cm;
        G4VSolid* septumCoilCylinderSolid = new G4Tubs("septumCoilCylinderTub",
                0,pSeptumCoilRadius,pSeptumCoilThickness/2.0,0,360*deg);
  //8
        double septumCoilRecX = 25.0*cm;
        double septumCoilRecY = 15.0*cm;
        G4VSolid* septumCoilRecSolid = new G4Box("septumCoilRecBox",
                septumCoilRecX/2.0,septumCoilRecY/2.0,pSeptumCoilThickness/2.0);

        G4IntersectionSolid* septumCoilSolid = new G4IntersectionSolid("septumCoilSolid",
                septumCoilCylinderSolid,septumCoilRecSolid,0,
                G4ThreeVector(-pSeptumCoilCenterX-septumCoilRecX/2.0,-pSeptumCoilCenterY+septumCoilRecY/2.0,0));

        G4LogicalVolume* septumCoilLogical = new G4LogicalVolume(septumCoilSolid,
                mMaterialManager->copper,"septumCoilLogical",0,0,0);

        G4VisAttributes *CuBrownVisAtt = new G4VisAttributes(G4Colour(1.0,0.5,0.5));

        septumCoilLogical->SetVisAttributes(CuBrownVisAtt);

        //place 16 copies into the septum container
        double pSeptumCoilPos_X_in   = pSeptumTunnelPos_X-pSeptumTunnelX/2.0-pSeptumCoilThickness/2.0;
        double pSeptumCoilPos_X_out  = pSeptumTunnelPos_X+pSeptumTunnelX/2.0+pSeptumCoilThickness/2.0;
        double pSeptumCoilPos_Y      = pSeptumCoilCenterY-pSeptumTunnelY/2.0;
        double pSeptumCoilPos_Z_up   = pSeptumPos_Z-pSeptumZ/2.0+pSeptumCoilCenterX;
        double pSeptumCoilPos_Z_down = pSeptumPos_Z+pSeptumZ/2.0-pSeptumCoilCenterX;


        //place up stream coils in the following order (looking downstream)
        //#####2######1###7######6#####
        //######      |###|      |#####
        //######      |###|      |#####
        //#####3######0###4######5#####
        G4RotationMatrix* pSeptumCoilRotFrontDown = new G4RotationMatrix();
        pSeptumCoilRotFrontDown->rotateY(90*deg);
        G4RotationMatrix* pSeptumCoilRotFrontUp = new G4RotationMatrix();
        pSeptumCoilRotFrontUp->rotateY(90*deg);
        pSeptumCoilRotFrontUp->rotateX(180*deg);

        if( placeseptum ){
        new G4PVPlacement(pSeptumCoilRotFrontDown,
                G4ThreeVector(pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_up+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,0,0);
        new G4PVPlacement(pSeptumCoilRotFrontUp,
                G4ThreeVector(pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_up+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,1,0);
        new G4PVPlacement(pSeptumCoilRotFrontUp,
                G4ThreeVector(pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_up+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,2,0);
        new G4PVPlacement(pSeptumCoilRotFrontDown,
                G4ThreeVector(pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_up+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,3,0);

        new G4PVPlacement(pSeptumCoilRotFrontDown,
                G4ThreeVector(-pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_up+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,4,0);
        new G4PVPlacement(pSeptumCoilRotFrontUp,
                G4ThreeVector(-pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_up+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,5,0);
        new G4PVPlacement(pSeptumCoilRotFrontUp,
                G4ThreeVector(-pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_up+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,6,0);
        new G4PVPlacement(pSeptumCoilRotFrontDown,
                G4ThreeVector(-pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_up+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,7,0);
        }

        //place down stream coils in the following order (looking downstream)
        //####10######9###15#####14####
        //######      |###|      |#####
        //######      |###|      |#####
        //####11######8###12#####13####
        G4RotationMatrix* pSeptumCoilRotBackDown = new G4RotationMatrix();
        pSeptumCoilRotBackDown->rotateY(270*deg);
        G4RotationMatrix* pSeptumCoilRotBackUp = new G4RotationMatrix();
        pSeptumCoilRotBackUp->rotateY(270*deg);
        pSeptumCoilRotBackUp->rotateX(180*deg);

        if( placeseptum ){
        new G4PVPlacement(pSeptumCoilRotBackDown,
                G4ThreeVector(pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,8,0);
        new G4PVPlacement(pSeptumCoilRotBackUp,
                G4ThreeVector(pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,9,0);
        new G4PVPlacement(pSeptumCoilRotBackUp,
                G4ThreeVector(pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,10,0);
        new G4PVPlacement(pSeptumCoilRotBackDown,
                G4ThreeVector(pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,11,0);

        new G4PVPlacement(pSeptumCoilRotBackDown,
                G4ThreeVector(-pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,12,0);
        new G4PVPlacement(pSeptumCoilRotBackUp,
                G4ThreeVector(-pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,13,0);
        new G4PVPlacement(pSeptumCoilRotBackUp,
                G4ThreeVector(-pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,14,0);
        new G4PVPlacement(pSeptumCoilRotBackDown,
                G4ThreeVector(-pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down+fPivotZOffset),
                septumCoilLogical,"septumCoilPhys",pMotherLogVol,true,15,0);
        }
    

        if( fSnakeModel == 49 || fSnakeModel == 48 || fSnakeModel == 51 || fSnakeModel == 53){

	  /////////////////////////
	  // HRS Q1              // Nickie is also adding the Q1 collimator here
	  /////////////////////////
	  double pQ1Rin=15.0*cm;
	  //double pQ1Length=(1.698M-1.36*m)*2+0.8*m;

	  //Nickie puts in Paul's collimator////////////////////////////////////////////////////////////////
	  //Paul's collimator///////////////////////////////////////////////////////////////////////////////
	  double PaulColT  = 0.01 * m; //This is obviously just a placeholder.
	  double buffer2   = 4.5 * cm;
	  double box2h     = 11.7 * cm - 20 * cm * sin(acos(18.9/20.0));
	  double box2w     = 4.0 * cm;
	  double box3h     = 2.0 * cm;
	  double box3w     = 10.0 * cm;
	  //double box3h     = 10.0 * cm;
	  //double box3w     = 20.0 * cm;
	  double box3x     = 4.0 * cm;
	  double box3y     = -1.88 * box3x + 14.74 * cm;
	  double box3deg   = 28.0 * deg;

	  double pHallCenter2Col=col_distance;//This is Nickie's correction: 30 cm from Q1 entrance
	  double pHRSVBPos_X=(pHallCenter2Col+PaulColT/2) * cos(fHRSAngle);
	  double pHRSVBPos_Y=(pHallCenter2Col+PaulColT/2) * sin(fHRSAngle);

	  G4RotationMatrix *rot28=new G4RotationMatrix();
	  rot28->rotateY( 90 * deg);
	  rot28->rotateX( 90 * deg + 28 * deg);
	  G4RotationMatrix *rot28m=new G4RotationMatrix();
	  rot28m->rotateY( 90 * deg);
	  rot28m->rotateX( - 90 * deg - 28 * deg);

	  G4VSolid* collCircle1 = new G4Tubs("collCirlce1", 0,              2. * pQ1Rin,        PaulColT, 0.,               360.0 * deg);
	  G4VSolid* collCircle2 = new G4Tubs("collCirlce2", 20.5*cm-buffer2,20.5*cm,      2*PaulColT,
					     -asin(11.7/20.5), 2 * asin(11.7/20.5));
	  G4VSolid* collCircle3 = new G4Tubs("collCirlce3", 20.0*cm,      20.0*cm+buffer2,2*PaulColT, -acos(18.9/20.0), 2 * acos(18.9/20.0));
	  G4Box*    collBox1    = new G4Box ("collBox1"   , 2 * PaulColT, 11.7 * cm,   (2.33+2.9)/2.0 * cm);
	  G4Box*    collBox2    = new G4Box ("collBox2"   , 2 * PaulColT, box2h / 2, box2w / 2);
	  G4Box*    collBox3    = new G4Box ("collBox3"   ,     PaulColT, box3h / 2, box3w / 2);

	  G4SubtractionSolid* subtraction1 = new G4SubtractionSolid("subtraction1", collCircle1,  collCircle2,
								    0, G4ThreeVector(-14.5 * cm, 0., 0.));
	  //pRotX90deg, G4ThreeVector(14.5, 0., 0.));
	  G4SubtractionSolid* subtraction2 = new G4SubtractionSolid("subtraction2", subtraction1, collCircle3,
								    0, G4ThreeVector(-22.9 * cm, 0., 0.));
	  G4SubtractionSolid* subtraction3 = new G4SubtractionSolid("subtraction3", subtraction2, collBox1,
								    0, G4ThreeVector(.285 * cm, 0., 0.));
	  G4SubtractionSolid* subtraction4
	    = new G4SubtractionSolid("subtraction4", subtraction3, collBox2,
				     0, G4ThreeVector(-0.5 * box2w,
						      20 * cm * sin(acos(18.9/20.0)) + 0.5 * box2h,
						      0.));
	  G4SubtractionSolid* subtraction5
	    = new G4SubtractionSolid("subtraction5", subtraction4, collBox2,
				     0, G4ThreeVector(-0.5 * box2w,
						      -20 * cm * sin(acos(18.9/20.0)) - 0.5 * box2h,
						      0.));
	  G4UnionSolid* subtraction6
	    = new G4UnionSolid("subtraction6", subtraction5, collBox3,
			       //rot28, G4ThreeVector(- box3x - 0.5 * box3w * cos(box3deg),
			       rot28, G4ThreeVector(- box3x - 0.5 * box3h * cos(box3deg),
						    //box3y + 0.5 * box3w * sin(box3deg),
						    box3y + 0.5 * box3h * sin(box3deg),
						    0.));
	  G4UnionSolid* subtraction7
	    = new G4UnionSolid("subtraction7", subtraction6, collBox3,
			       //rot28m, G4ThreeVector(- box3x - 0.5 * box3w * cos(box3deg),
			       rot28m, G4ThreeVector(- box3x - 0.5 * box3h * cos(box3deg),
						     //- box3y - 0.5 * box3w * sin(box3deg),
						     - box3y - 0.5 * box3h * sin(box3deg),
						     0.));
	  
	  G4LogicalVolume* PaulLogical = new G4LogicalVolume(subtraction7, mMaterialManager->tungsten, "PaulLogical", 0, 0, 0);
	  G4RotationMatrix *pRotLHRScol=new G4RotationMatrix();
	  //pRotLHRScol->rotateZ(180 * deg);
	  pRotLHRScol->rotateY(-fHRSAngle); 
	  G4RotationMatrix *pRotRHRScol=new G4RotationMatrix();
	  
	  pRotRHRScol->rotateZ(180 * deg);
	  pRotRHRScol->rotateY(-fHRSAngle); 
	  //pRotRHRScol->rotateX(90. * deg); 

	  PaulLogical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
	  
	  //new G4PVPlacement(pRotRHRScol,G4ThreeVector(-pHRSVBPos_Y,0,pHRSVBPos_X),
	  //PaulLogical,"PaulPhys",motherLogical,0,0,0);
	  //new G4PVPlacement(pRotLHRScol,G4ThreeVector(pHRSVBPos_Y,0,pHRSVBPos_X),
	  //PaulLogical,"PaulPhys",motherLogical,0,0,0);
	  //End of Paul's collimator////////////////////////////////////////////////////////////////////////   	  
	  
	}
	else if( fSnakeModel != 49 || fSnakeModel > 50)
	{
		//place VB @ Septum entrance window, 10.4cm width and 24.4cm height, 
		//The following declared as module variales already
		//double mHRSVBWidth=104*mm;
		//double mHRSVBHeight=244*mm;
		//double mHRSVBThick=0.0508*mm; 
		double mHRSVBWidth=800*mm;
		double mHRSVBHeight=800*mm;
		double mHRSVBThick=0.0508*mm; 
		//acceptance is 20 mrad

                G4VisAttributes *LightYellowVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.2)); 
	  
		G4VSolid* septumWindowSolid = new G4Box("septumWindowBox",mHRSVBWidth/2.0,
			mHRSVBHeight/2.0,mHRSVBThick/2.0);
		G4LogicalVolume* septumWindowLogical = new G4LogicalVolume(septumWindowSolid,
		        mMaterialManager->mylar,"septumWindowLogical",0,0,0);
		SDman->AddNewDetector(septumWindowSD);
		septumWindowLogical->SetSensitiveDetector(septumWindowSD);
		septumWindowLogical->SetVisAttributes(LightYellowVisAtt); 

		G4VSolid* TargetSolid = new G4Tubs("TargetTub",0.0,50*cm,
			mHRSVBThick/2.0,0.0,360.0*deg);
		G4LogicalVolume* TargetLogical = new G4LogicalVolume(TargetSolid,
			mMaterialManager->mylar,"TargetLogical",0,0,0);
		SDman->AddNewDetector(septumWindowSD);
		TargetLogical->SetSensitiveDetector(septumWindowSD);
		TargetLogical->SetVisAttributes(LightYellowVisAtt); 

                double mPivot2HRSVBFace = 1347.17*mm;

		double pTunnel2Beam_X=8.4*cm;
		//put both left and right septum entrance window, which should not less than pTunnel2Beam_X
		if(fSetupHRS){
		  double pLSeptumWindowPos_X=(mPivot2HRSVBFace+mHRSVBThick/2.0)*
		    sin(fSeptumAngle)+fPivotXOffset;
		  if(mPivot2HRSVBFace>1190*mm){
		    //need to shift it in X to make it barely touch the septum tunnel		
		    pLSeptumWindowPos_X=mHRSVBWidth/2*cos(fSeptumAngle)+pTunnel2Beam_X;
		  }
		  double pLSeptumWindowPos_Y=fPivotYOffset;
		  double pLSeptumWindowPos_Z=(mPivot2HRSVBFace+mHRSVBThick/2.0)*
		    cos(fSeptumAngle)+fPivotZOffset;
		  //new G4PVPlacement(pRotLSeptum,
		  //G4ThreeVector(pLSeptumWindowPos_X,pLSeptumWindowPos_Y,pLSeptumWindowPos_Z),
		  //septumWindowLogical,"virtualBoundaryPhys_LHRS",motherLogical,0,0,0);
		}if(fSetupHRS){
		  double pRSeptumWindowPos_X=(mPivot2HRSVBFace+mHRSVBThick/2.)*
		    sin(-fSeptumAngle)+fPivotXOffset;			
		  if(mPivot2HRSVBFace>1190*mm){
		    //need to shift it in X to make it barely touch the septum tunnel		
		    pRSeptumWindowPos_X=-mHRSVBWidth/2*cos(-fSeptumAngle)-pTunnel2Beam_X;
		  }
		  double pRSeptumWindowPos_Y=fPivotYOffset;
		  double pRSeptumWindowPos_Z=(mPivot2HRSVBFace+mHRSVBThick/2.0)*
		    cos(-fSeptumAngle)+fPivotZOffset;
		  //new G4PVPlacement(pRotRSeptum,
		  //G4ThreeVector(pRSeptumWindowPos_X,pRSeptumWindowPos_Y,pRSeptumWindowPos_Z),
		  //septumWindowLogical,"virtualBoundaryPhys_RHRS",motherLogical,0,0,0);
		}
	}


	return;
}

void g4hrsDetectorConstruction::CreateTargetChamber(G4LogicalVolume *pMotherLogVol){
	const double inch=2.54*cm;

	G4RotationMatrix *pRotScatInHall=new G4RotationMatrix();
	pRotScatInHall->rotateX(90.*deg);

//        g4hrsMaterial* mMaterialManager = g4hrsMaterial::GetHRSMaterialManager();
	double startphi=0.*deg, deltaphi=360.*deg;
	/////////////////////////////////////////////////////////
	//scattering chamber container
	/////////////////////////////////////////////////////////
	//This is just a container to enclose the taraget chamber, 10 mm larger than the 
	//scattering chamber itself.
	//With this container, all stuff inside do not need a rotation
	


	double pScatChamberContainerRin=fScatChamberRin-1*cm;      
	double pScatChamberContainerRout=fScatChamberRout+1*cm;
	double pScatChamberContainerL=fScatChamberL+(3.50+17.0+1.25)*inch*2+10.0*mm;
	G4VSolid* scatChamberContainerExtendedSolid = new G4Tubs("scatChamberContainerExtendedTubs",
		pScatChamberContainerRin,pScatChamberContainerRout,
		pScatChamberContainerL/2.0,0.,360.*deg);
	G4VSolid* scatChamberContainerExtraSolid = new G4Tubs("scatChamberContainerExtraTubs",
		0,pScatChamberContainerRout+1*mm,17.25*inch/2.0,0.,360.*deg);
	G4SubtractionSolid* scatChamberContainerSolid=new G4SubtractionSolid("scatChamberContainerSolid",
		scatChamberContainerExtendedSolid,scatChamberContainerExtraSolid,
		0,G4ThreeVector(0,0,-fScatChamberL/2-17.25*inch/2.0));


        G4VisAttributes *HallVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
        HallVisAtt->SetVisibility(false);


	G4LogicalVolume* scatChamberContainerLogical = new G4LogicalVolume(scatChamberContainerSolid,
		mMaterialManager->vacuum,"scatChamberContainerLogical",0,0,0);
	scatChamberContainerLogical->SetVisAttributes(HallVisAtt); 

	G4VPhysicalVolume* scatChamberContainerPhys=new G4PVPlacement(pRotScatInHall,
		G4ThreeVector(fScatChamberXOffset,fScatChamberYOffset,fScatChamberZOffset),
		scatChamberContainerLogical,"scatChamberContainerPhys",pMotherLogVol,0,0);

	//////////////////////////
	// build the target chamber.
	//////////////////////////
	//Build the simplified scatter chamber,it contains 2 windows of rectangles 
	//The following already defined in the config file
	//double fScatChamberRin=17.875*inch,fScatChamberRout=18.875*inch,fScatChamberL=27.25*inch;

	G4VSolid* scatChamberWholeSolid=0;
	//If mSetupScatChamber==1, setup the body only, 
	//If mSetupScatChamber==2, setup the body plus top flange and bottom flange, this
	//will make the program slower


	if(fSetupStdScatChamber==1)
	{
	  scatChamberWholeSolid = new G4Tubs("scatChamberWholeTubs",
					     fScatChamberRin,fScatChamberRout,fScatChamberL/2.0,0.,360.*deg);
	}
	else if(fSetupStdScatChamber>=2)
	{
		startphi=0.*deg; deltaphi=360.*deg;
		const int kNPlane_SC=11;
		double rInner_SC[] = {0,0,fScatChamberRin,
			fScatChamberRin,fScatChamberRin,fScatChamberRin,
			fScatChamberRin,fScatChamberRin,fScatChamberRin,
			0,0};
		double rOuter_SC[] = {
			fScatChamberRout+1.0*inch,fScatChamberRout+1.0*inch,fScatChamberRout+1.0*inch,
			fScatChamberRout,fScatChamberRout,fScatChamberRout+1.0*inch,
			fScatChamberRout+1.0*inch,fScatChamberRout,fScatChamberRout,
			fScatChamberRout+1*inch,fScatChamberRout+1*inch
		};
		double zPlane_SC[] = {
			-fScatChamberL/2-4.50*inch,-fScatChamberL/2-3.25*inch,-fScatChamberL/2-1.0*inch,
			-fScatChamberL/2,fScatChamberL/2+0.25*inch,fScatChamberL/2+1.25*inch,
			fScatChamberL/2+3.50*inch,fScatChamberL/2+3.50*inch,fScatChamberL/2+20.5*inch,
			fScatChamberL/2+20.5*inch,fScatChamberL/2+21.75*inch
		};

		G4Polycone* SCWholeSolid = new G4Polycone("SCPolycone",startphi,deltaphi,
			kNPlane_SC,zPlane_SC,rInner_SC,rOuter_SC);

		scatChamberWholeSolid = SCWholeSolid;
	}


	//these are the subtraction part, not the scatter chamber itself
	double pSCWindowRin=fScatChamberRin-1*mm;
	double pSCWindowRout=fScatChamberRout+1*mm;
	double pSCEntranceWindowH=6.44*inch;
	double pSCDownCapH=15.0*inch;
	
	//rectangle EntranceWindow covering 80 to 100 degrees
	startphi=80*deg;deltaphi=20*deg;
	G4VSolid* SCEntranceWindowSolid = new G4Tubs("SCEntranceWindowTubs",
		pSCWindowRin,pSCWindowRout,pSCEntranceWindowH/2.0,startphi,deltaphi);
	//rectangle DownCap covering -225 to 45 degrees
	startphi=-225*deg;deltaphi=270*deg;
	G4VSolid* SCDownCapSolid = new G4Tubs("SCDownCapTubs",
		pSCWindowRin,pSCWindowRout,pSCDownCapH/2.0,startphi,deltaphi);
	// subtract the Entrance window 
	G4SubtractionSolid* SCSubtractEntranceSolid=new G4SubtractionSolid(
		"SCSubtractEntrance",scatChamberWholeSolid,SCEntranceWindowSolid);

	// subtract the Exit window 
	G4SubtractionSolid* SCSubtractEntranceNExitSolid=new G4SubtractionSolid(
		"SCSubtractEntranceNExit",SCSubtractEntranceSolid,SCDownCapSolid);

        G4VisAttributes *WhiteVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0)); 

	//setup the scatter chamber
	G4LogicalVolume* scatChamberLogical = new G4LogicalVolume(
		SCSubtractEntranceNExitSolid,mMaterialManager->aluminum,"scatChamberLogical",0,0,0);
	scatChamberLogical->SetVisAttributes(WhiteVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),
		scatChamberLogical,"scatChamberPhys",scatChamberContainerLogical,0,0);

	
	/////////////////////////
	// target chamber window covers 
	/////////////////////////
	//Covers for EntranceWindow 

        double mScatChamberEntranceWindowThick = 0.254*mm;
        double mScatChamberExitWindowThick = 0.508*mm;
        G4VisAttributes *LightYellowVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.2)); 

	//EntranceWindowCover
	double pSCEntranceWindowCoverH=pSCEntranceWindowH+0.8*inch;
	double pSCEntranceWindowCoverRin=fScatChamberRout;
	double pSCEntranceWindowCoverRout=pSCEntranceWindowCoverRin+mScatChamberEntranceWindowThick;

	startphi=78*deg;deltaphi=24*deg;
	G4VSolid* SCEntranceWindowCoverSolid = new G4Tubs("SCEntranceWindowCoverTubs",
		pSCEntranceWindowCoverRin,pSCEntranceWindowCoverRout,
		pSCEntranceWindowCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCEntranceWindowCoverLogical = new G4LogicalVolume(
		SCEntranceWindowCoverSolid,mMaterialManager->aluminum,"SCEntranceWindowCoverLogical",0,0,0);
	SCEntranceWindowCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),SCEntranceWindowCoverLogical,
		"SCEntranceWindowCoverPhys",scatChamberContainerLogical,false,0);

	//DownCapCover
	double pSCDownCapCoverH=pSCDownCapH+0.8*inch;
	double pSCDownCapCoverRin=fScatChamberRout;
	double pSCDownCapCoverRout=fScatChamberRout+mScatChamberExitWindowThick;

	startphi=-227*deg;deltaphi=274*deg;
	G4VSolid* SCDownCapCoverSolid = new G4Tubs("SCDownCapCoverTubs",
		pSCDownCapCoverRin,pSCDownCapCoverRout,pSCDownCapCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCDownCapCoverLogical = new G4LogicalVolume(
		SCDownCapCoverSolid,mMaterialManager->aluminum,"SCDownCapCoverLogical",0,0,0);
	SCDownCapCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),SCDownCapCoverLogical,
		"SCDownCapCoverPhys",scatChamberContainerLogical,false,0);




	return;

}


void g4hrsDetectorConstruction::CreateHRS(G4LogicalVolume* motherLogical)
{
//        g4hrsMaterial* mMaterialManager = g4hrsMaterial::GetHRSMaterialManager();

	int IS_NIM = 0;//if it's 0, then we go by SNAKE, if it is 1, we go by NIM.
	//As it turns out, I think we want to go by SNAKE, so that is why it is zero right now.
	G4VPhysicalVolume* theHRSPhys=0;

	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 
	G4RotationMatrix *pRotX45deg=new G4RotationMatrix();
	pRotX45deg->rotateX(-45*deg); 
	G4RotationMatrix *pRotX30deg=new G4RotationMatrix();
	pRotX30deg->rotateX(60*deg); 
	G4RotationMatrix *pRotX105deg=new G4RotationMatrix();
	pRotX105deg->rotateX(165*deg); 
	G4RotationMatrix *pRotXLHRSdeg=new G4RotationMatrix();
	pRotXLHRSdeg->rotateY(fHRSAngle); 
	G4RotationMatrix *pRotXRHRSdeg=new G4RotationMatrix();
	pRotXRHRSdeg->rotateY(-fHRSAngle); 

	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4String SDname;
	G4VSensitiveDetector* Q1WindowSD=new g4hrsGenericDetector(SDname="Q1Window", 2);
	//Build the HRS, positions were taken fron the NIM paper
	//The Apertures are defined as:
	//Q1 exit is a circle of radius 0.15m
	//Q3 entrance and exit are circles of radius 0.30m
	//The dipole entrance and exit are trapezoids:
	//-0.40 m<x<0.40 m (+x is down) 
	//y=+-(0.125*(1-(1.25*x/8.40)) (smallest gap at positive x, x in m)

	/////////////////////////
	// HRS QQDQ Containner
	/////////////////////////
	//Build a container using polycone, covering 270+/-12 degrees, 1.46m to 20.76m
	//3.05m below the beam line is the ground, need to subtract everything beneath that
	//looks like this:
	/*
	HRS container:covering 24 degrees in X-Z plane, 10.05 m height.
	Stuff inside this containner will also use the pivot as the origin, but do not 
	need to worry about the rotation of the HRS.
	//                         7.0m above beam line 
	//                       -----------------------------| 20.76m from pivot,
	//                      /                             |
	//                     /                        Q3    |
	//                    /                               |
	//                   /                       E        /
	//                  /                     L           |
	//          --------                   O              |
	//         /                        P                 |
	//    -----                     I                     |
	//----|----- Q1 --- Q2 ---- D     --------------------|------ beam line -----
	//    |                                               |
	//    ------------------------------------------------|
	//    1.46m from povot, 3.05 m below beam line

	*/

	//double pHRSContainerRin=1.46*m,pHRSContainerRout=20.76*m;//JixieMode
	double pHRSContainerRin=1.46*m,pHRSContainerRout=25*m;//NickieMode it just has to be big enough to fit my VB at FP
	//double pBeamLine2Ground;
	//pBeamLine2Ground=-3.05*m;
	//build the container with polycone

	const int kNPlane_HRSContainer=7;
	double rInner_HRSContainer[] = {pHRSContainerRin,pHRSContainerRin,2.5*m,
		3.7*m,9.0*m,pHRSContainerRout-3.0*m,pHRSContainerRout};
	double rOuter_HRSContainer[] = {pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,
		pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,pHRSContainerRout};
	//double zPlane_HRSContainer[] = {-2.0*m,1.0*m,1.0*m,
	//2.0*m,2.0*m,7.0*m,7.0*m};
	double zPlane_HRSContainer[] = {-2.0*m,1.0*m,1.0*m,
					2.0*m,2.0*m,11.0*m,11.0*m};
	G4Polycone* HRSContainerSolid = new G4Polycone("HRSContainer",258.0*deg,24.0*deg,
		kNPlane_HRSContainer,zPlane_HRSContainer,rInner_HRSContainer,rOuter_HRSContainer);

	////build the container using tube
	//double pHRSContainerHeight=14.0*m;
	//G4VSolid* HRSContainerTub = new G4Tubs("HRSContainerTub",pHRSContainerRin,
	//	pHRSContainerRout,pHRSContainerHeight/2.0,258.0*deg,24.0*deg);
	////ground
	//double pGoundHeight=10*m;
	//G4VSolid* GroundTub = new G4Tubs("GroundTub",0,30*m,pGoundHeight/2,0*deg,360.0*deg);
	////tube subtract the ground
	//G4SubtractionSolid* HRSContainerSolid = new G4SubtractionSolid("HRSContainer",
	//		HRSContainerTub,GroundTub,
	//		0,G4ThreeVector(0,0,pBeamLine2Ground-pGoundHeight/2));

	G4LogicalVolume* LHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		mMaterialManager->vacuum,"LHRSContainerLogical",0,0,0);
	G4LogicalVolume* RHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		mMaterialManager->vacuum,"RHRSContainerLogical",0,0,0);

	G4VisAttributes* MagFieldVisAtt = new G4VisAttributes(G4Colour(1., 0., 1.));	

        G4VisAttributes *HallVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
        HallVisAtt->SetVisibility(false);


	LHRSContainerLogical->SetVisAttributes(HallVisAtt); 
	RHRSContainerLogical->SetVisAttributes(HallVisAtt);
	//RHRSContainerLogical->SetVisAttributes(MagFieldVisAtt); 

	G4RotationMatrix *pRotLHRSContainer=new G4RotationMatrix();
	pRotLHRSContainer->rotateX(-270*deg);
	pRotLHRSContainer->rotateZ(-fHRSAngle);  

	G4RotationMatrix *pRotRHRSContainer=new G4RotationMatrix();
	pRotRHRSContainer->rotateX(-270*deg);
	pRotRHRSContainer->rotateZ(fHRSAngle);  

	if(fSetupHRS>=2)
	{
		new G4PVPlacement(pRotLHRSContainer,G4ThreeVector(0,0,0.+fPivotZOffset),
			LHRSContainerLogical,"LHRSContainerPhys",motherLogical,0,0,0);
	}
	if(fSetupHRS>=2)
	{
		new G4PVPlacement(pRotRHRSContainer,G4ThreeVector(0,0,0.+fPivotZOffset),
			RHRSContainerLogical,"RHRSContainerPhys",motherLogical,0,0,0);
	}


	//////////////////////////
	//Q1 Entrance Window //
	//////////////////////////
	bool pSetupQ1EntranceWindow=true;
	//G2P and CREX have their own VB defined at a different location
	if( fSetupCREXGeometry)  pSetupQ1EntranceWindow=false;
	if(pSetupQ1EntranceWindow) 
	{

		//this part is trying to place a virtual boundary at the Q1 entrance

		//Place both left and right VB for HRS, which is pHRSContainerRin+4*mm away from the 
		//hall center(1.462m). This aperture is a round disk of 29.8 cm diameter and 2 mm thick
		//The real Q1 vacumn entrance to hall center is 1.312m, 

		double pHRSQ1WinThick = 2*mm;
		G4VSolid* HRSQ1WinSolid = new G4Tubs("HRSQ1WinTub",0.0,14.9*cm,
			pHRSQ1WinThick/2,0.0,360.0*deg);
		G4LogicalVolume* HRSQ1WinLogical = new G4LogicalVolume(HRSQ1WinSolid,
			mMaterialManager->mylar,"HRSQ1WinLogical",0,0,0);
		SDman->AddNewDetector(Q1WindowSD);
		HRSQ1WinLogical->SetSensitiveDetector(Q1WindowSD);

                G4VisAttributes *LightYellowVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.2)); 
		HRSQ1WinLogical->SetVisAttributes(LightYellowVisAtt); 

		//since the container has been rotated by 90 deg about x axis,
		//y'=z  z'=-y ==> I have to put this offset as -y
		double pHallCenter2Q1Win=pHRSContainerRin+4*mm;  //place it at the first 1.464 m

	}
	//BELOW ARE SNAKE VALUES, NOT NIM VALUES
	double pTarget =             0.0  * cm;
	double pQ1en   = pTarget + 160.0  * cm;
	double pQ1ex   = pQ1en   +  94.13 * cm;
	double pQ2en   = pQ1ex   + 115.58 * cm;
	double pQ2ex   = pQ2en   + 182.66 * cm;
	//ABOVE ARE SNAKE VALUES, NOT NIM VALUES
	
	/////////////////////////
	// HRS Q1              // Nickie is also adding the Q1 collimator here
	/////////////////////////
	double pHallCenter2LQ1Face;//=1.69*m;//NIM
	double pHallCenter2RQ1Face;//=1.69*m;//NIM
	//pHallCenter2LQ1Face = ( IS_NIM == 1 ) ? 169.0 * cm : pQ1en;
	//pHallCenter2RQ1Face = ( IS_NIM == 1 ) ? 169.0 * cm : pQ1en;
	pHallCenter2LQ1Face = pQ1en;
	pHallCenter2RQ1Face = pQ1en;
	//double pHallCenter2LQ1FaceMag=1.600*m;//SNAKE
	//double pHallCenter2RQ1FaceMag=1.600*m;//SNAKE
	//double fringe_extension = 1.25;
	//double fringe_extension = 1.2;
	double fringe_extension = 1.00;

	//flag to trigger sos quad instead of old quad for the proper runs
	//K. Allada and B. Schmookler:
	//Magnetic Length = 70 cm
	//Radius to Pole Tip = 12.827  cm
	
	int sos = fTune->sosQuad;

	double q1shift = sos ? 0.0 * m : 0.0 * m;
	double pQ1Rin  = sos ? 12.827 * cm : 15.0*cm;
	double pQ1Rout = sos ? 35.0   * cm : 35.0*cm;//for now, keep the outer same for either
	double pQ1Length = sos ?  70. * cm : pQ1ex - pQ1en;
	double pQ1PosAct = pQ1ex - pQ1en;
	double pQ1LengthMag=pQ1Length * fringe_extension;
	//double pQ1LengthMag=94.1*cm;
	//double pQ1Length=(1.698M-1.36*m)*2+0.8*m;


	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q1Solid = new G4Tubs("Q1Tub",pQ1Rin,pQ1Rout,pQ1Length/2.0,0.0,360.0*deg);

	//build 2 copy since there are differennt fields involved in it
	G4LogicalVolume* LQ1Logical = new G4LogicalVolume(Q1Solid,
		mMaterialManager->siliconsteel,"LQ1Logical",0,0,0);
	G4LogicalVolume* RQ1Logical = new G4LogicalVolume(Q1Solid,
		mMaterialManager->siliconsteel,"RQ1Logical",0,0,0);

        G4VisAttributes *IronVisAtt = new G4VisAttributes(G4Colour(100./255.,149./255.,237./255.));
	LQ1Logical->SetVisAttributes(IronVisAtt); 
	RQ1Logical->SetVisAttributes(IronVisAtt); 

	if(fSetupHRS>=2){
	  //put it in the container, which also center at the hall center
	  //therefore only the z_at_lab position need to be considered
	  double pLQ1Pos_Z=(pHallCenter2LQ1Face+pQ1PosAct/2.0) + q1shift;
	  //double pLQ1Pos_Z=pQ1PosAct + q1shift;
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z,0),
			    LQ1Logical,"LQ1Phys",LHRSContainerLogical,0,0,0);
	}
	if(fSetupHRS>=2){
	  double pRQ1Pos_Z=(pHallCenter2RQ1Face+pQ1PosAct/2.0) + q1shift;
	  //double pRQ1Pos_Z=pQ1PosAct + q1shift;
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z,0),
			    RQ1Logical,"RQ1Phys",RHRSContainerLogical,0,0,0);
	}

	//vac
	double pQ1vacRin  = pQ1Rin;
	double pQ1vacRout = pQ1Rin + 5. * cm;
	double pQ1vacLength = 610.4 * mm; 
	double pQ1vacCenterZ = pHallCenter2LQ1Face + pQ1PosAct + pQ1vacLength / 2.;
	double pQ1vacCenterY = 0;
	G4VSolid* Q1vacSolid = new G4Tubs("Q1vacTub",pQ1vacRin,pQ1vacRout,pQ1vacLength/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ1vacLogical = new G4LogicalVolume(Q1vacSolid,
		mMaterialManager->siliconsteel,"LQ1vacLogical",0,0,0);
	G4LogicalVolume* RQ1vacLogical = new G4LogicalVolume(Q1vacSolid,
		mMaterialManager->siliconsteel,"RQ1vacLogical",0,0,0);

	LQ1vacLogical->SetVisAttributes(IronVisAtt); 
	RQ1vacLogical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ1vacInContainer=new G4RotationMatrix();
	pRotQ1vacInContainer->rotateX(90*deg); 
	
	if(fSetupHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ1vacInContainer,G4ThreeVector(0, -pQ1vacCenterZ,pQ1vacCenterY),
			LQ1vacLogical,"LQ1vacPhys",LHRSContainerLogical,0,0,0);
	}
	if(fSetupHRS>=4)
	{
		new G4PVPlacement(pRotQ1vacInContainer,G4ThreeVector(0, -pQ1vacCenterZ,pQ1vacCenterY),
			RQ1vacLogical,"RQ1vacPhys",RHRSContainerLogical,0,0,0);
	}
	

	/////////////////////////
	// HRS Q2 
	/////////////////////////
	double pHallCenter2LQ2Face;//=3.74*m;//NIM
	double pHallCenter2RQ2Face;//=3.74*m;//NIM
	//pHallCenter2LQ2Face = ( IS_NIM == 1 ) ? 374.0 * cm : pQ2en;
	//pHallCenter2RQ2Face = ( IS_NIM == 1 ) ? 374.0 * cm : pQ2en;
	pHallCenter2LQ2Face = pQ2en;
	pHallCenter2RQ2Face = pQ2en;
	//double pHallCenter2LQ2FaceMag=3.696*m;//SNAKE
	//double pHallCenter2RQ2FaceMag=3.696*m;//SNAKE
	double pQ2Rin=30.0*cm;
	double pQ2Rout=75.0*cm;
	double pQ2Length;//=180*cm;//NIM
	//pQ2Length = ( IS_NIM == 1 ) ? 180.0 * cm : pQ2ex - pQ2en;
	pQ2Length = pQ2ex - pQ2en;
	double pQ2LengthMag=pQ2Length * fringe_extension;
	//double pQ2LengthMag=182.6*cm;//SNAKE

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q2Solid = new G4Tubs("Q2Tub",pQ2Rin,pQ2Rout,pQ2Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ2Logical = new G4LogicalVolume(Q2Solid,
		mMaterialManager->siliconsteel,"LQ2Logical",0,0,0);
	G4LogicalVolume* RQ2Logical = new G4LogicalVolume(Q2Solid,
		mMaterialManager->siliconsteel,"RQ2Logical",0,0,0);

	LQ2Logical->SetVisAttributes(IronVisAtt); 
	RQ2Logical->SetVisAttributes(IronVisAtt); 

	if(fSetupHRS>=3)
	{
		//put it in the container, which also center at the hall center
		double pLQ2Pos_Z=(pHallCenter2LQ2Face+pQ2Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z,0),
			LQ2Logical,"LQ2Phys",LHRSContainerLogical,0,0,0);
	}
	if(fSetupHRS>=3)
	{
		double pRQ2Pos_Z=(pHallCenter2RQ2Face+pQ2Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z,0),
			RQ2Logical,"RQ2Phys",RHRSContainerLogical,0,0,0);
	}

	//vac
	
	double pQ2vacRin  = pQ2Rin;
	double pQ2vacRout = pQ2Rin + 5. * cm;
	double pQ2vacLength = ( 1166.06 - 675. ) * mm;
	double pQ2vacCenterZ = pHallCenter2LQ2Face - pQ2vacLength / 2.;
	double pQ2vacCenterY = 0;
	G4VSolid* Q2vacSolid = new G4Tubs("Q2vacTub",pQ2vacRin,pQ2vacRout,pQ2vacLength/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ2vacLogical = new G4LogicalVolume(Q2vacSolid,
		mMaterialManager->siliconsteel,"LQ2vacLogical",0,0,0);
	G4LogicalVolume* RQ2vacLogical = new G4LogicalVolume(Q2vacSolid,
		mMaterialManager->siliconsteel,"RQ2vacLogical",0,0,0);

	LQ2vacLogical->SetVisAttributes(IronVisAtt); 
	RQ2vacLogical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ2vacInContainer=new G4RotationMatrix();
	pRotQ2vacInContainer->rotateX(90*deg); 
	if(fSetupHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ2vacInContainer,G4ThreeVector(0, -pQ2vacCenterZ,pQ2vacCenterY),
			LQ2vacLogical,"LQ2vacPhys",LHRSContainerLogical,0,0,0);
	}
	if(fSetupHRS>=4)
	{
		new G4PVPlacement(pRotQ2vacInContainer,G4ThreeVector(0, -pQ2vacCenterZ,pQ2vacCenterY),
			RQ2vacLogical,"RQ2vacPhys",RHRSContainerLogical,0,0,0);
	}
	

	//vac2
	double pQ2vac2Rin     = pQ2Rin;
	double pQ2vac2Rout    = pQ2Rin + 1. * mm;
	//double pQ2vac2Length  = 561.7 * mm;
	double pQ2vac2Length  = 226.7 * mm; // SPR Shortening to avoid overlaps with dipole
	double pQ2vac2CenterZ = pHallCenter2RQ2Face + pQ2Length + pQ2vac2Length / 2.;
	double pQ2vac2CenterY = 0;
	G4VSolid* Q2vac2Solid = new G4Tubs("Q2vac2Tub",pQ2vac2Rin,pQ2vac2Rout,pQ2vac2Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ2vac2Logical = new G4LogicalVolume(Q2vac2Solid,
		mMaterialManager->siliconsteel,"LQ2vac2Logical",0,0,0);
	G4LogicalVolume* RQ2vac2Logical = new G4LogicalVolume(Q2vac2Solid,
		mMaterialManager->siliconsteel,"RQ2vac2Logical",0,0,0);

	LQ2vac2Logical->SetVisAttributes(IronVisAtt); 
	RQ2vac2Logical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ2vac2InContainer=new G4RotationMatrix();
	pRotQ2vac2InContainer->rotateX(90*deg); 
	if(fSetupHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ2vac2InContainer,G4ThreeVector(0, -pQ2vac2CenterZ,pQ2vac2CenterY),
			LQ2vac2Logical,"LQ2vac2Phys",LHRSContainerLogical,0,0,0);
	}
	if(fSetupHRS>=4)
	{
		new G4PVPlacement(pRotQ2vac2InContainer,G4ThreeVector(0, -pQ2vac2CenterZ,pQ2vac2CenterY),
			RQ2vac2Logical,"RQ2vac2Phys",RHRSContainerLogical,0,0,0);
	}

	/////////////////////////
	// HRS Dipole 
	/////////////////////////
	//The dipole is built as a disc subtraced subtract another disc to get the side 
	//then subtract by the tunnel disc
	double pDipoleBendAngle=45.*deg, pDipoleFaceAngle=30.*deg;
	double pDipoleR=8.4*m;
	double pDipoleRprime=pDipoleR*sin(pDipoleBendAngle/2.)/
		sin(180.*deg-pDipoleBendAngle/2.-pDipoleFaceAngle);
	//double pDipoleRprime=4.0518*m;

	double pDipoleRCenterY=pDipoleR;
	double pDipoleRCenterZ;//=9.96*m;
	//pDipoleRCenterZ = ( IS_NIM == 1 ) ? 9.96 * m : 9.961 * m;//snake
	pDipoleRCenterZ = 9.961 * m;//snake
	//pDipoleRCenterZ = ( IS_NIM == 1 ) ? 9.96 * m : 9.9547 * m;//seamus
	G4cout << pDipoleRCenterZ << G4endl;
	//double pDipoleRprimeCenterY=pDipoleRprime*cos(pDipoleFaceAngle);			// =2.865*m;
	//double pDipoleRprimeCenterZ=9.96*m+pDipoleRprime*sin(pDipoleFaceAngle);	//12.825*m;

	//the center of Rprime relative to R 
	double pRprime2R_Y=pDipoleRprime*cos(pDipoleFaceAngle)-pDipoleR;	// =-5.535*m;
	double pRprime2R_Z=pDipoleRprime*sin(pDipoleFaceAngle);				// =1.865*m;

	//the original disc
	G4VSolid* DipoleWholeTub = new G4Tubs("DipoleWholeTub",
		pDipoleR-0.8*m,pDipoleR+0.8*m,0.4*m,
		172.*deg,pDipoleBendAngle+16.*deg);
	//the disc to be subtracted, musu be thicker 
	G4VSolid* DipolePrimeTub = new G4Tubs("DipolePrimeTub",
		0,pDipoleR,0.5*m,
		180.*deg+pDipoleFaceAngle+pDipoleBendAngle,360.*deg-pDipoleFaceAngle*2.-pDipoleBendAngle);
	//subtract the small tube to form the shape of the sides
	G4SubtractionSolid* DipoleWithSides = new G4SubtractionSolid("DipoleWithSides",
		DipoleWholeTub,DipolePrimeTub,
		0,G4ThreeVector(pRprime2R_Y,-pRprime2R_Z,0));

	//the tunnel disc, I use a rectangle shape here
	//G4VSolid* DipoleTunnelTub = new G4Tubs("DipoleTunnelTub",
	//	pDipoleR-0.4*m,pDipoleR+0.4*m,0.125*m,
	//	170*deg,pDipoleBendAngle+20*deg);
	//The shape of the dipole tunnel is a trapzoid, x=+/-0.4; y=+/-(0.125*(1-(1.25*x/8.40)) Jixie
	//To build this shape, I have to use polycone Jixie
	// -5.22008 < x < -4.98099 J.L.R.
	//  -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784 J.L.R.
	double dy=0.125*(1.25*0.40/8.40)*m;
	const int kNPlane_DipoleTunnel=4;
	double rInner_DipoleTunnel[] = {pDipoleR-0.4*m,pDipoleR-0.4*m,pDipoleR-0.4*m,pDipoleR-0.4*m};
	double rOuter_DipoleTunnel[] = {pDipoleR-0.4*m,pDipoleR+0.4*m,pDipoleR+0.4*m,pDipoleR-0.4*m};
	double zPlane_DipoleTunnel[] = {-0.125*m-dy,-0.125*m+dy,0.125*m-dy,0.125*m+dy};
	G4Polycone* DipoleTunnelCone = new G4Polycone("DipoleTunnelCone",
		    170.0*deg,pDipoleBendAngle+20.*deg,
		    kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);//Jixie's original
	
	//subtract the tunnel disc
	G4SubtractionSolid* DipoleSolid = new G4SubtractionSolid("DipoleSolid",
		DipoleWithSides,DipoleTunnelCone);

	G4LogicalVolume* LDipoleLogical = new G4LogicalVolume(DipoleSolid,
		mMaterialManager->siliconsteel,"LDipoleLogical",0,0,0);
	G4LogicalVolume* RDipoleLogical = new G4LogicalVolume(DipoleSolid,
		mMaterialManager->siliconsteel,"RDipoleLogical",0,0,0);

        G4VisAttributes *OrangeVisAtt = new G4VisAttributes(G4Colour(1.0,0.5,0.0));//orange


	LDipoleLogical->SetVisAttributes(OrangeVisAtt); 
	RDipoleLogical->SetVisAttributes(OrangeVisAtt); 	

	G4RotationMatrix *pRotDipoleInContainer=new G4RotationMatrix();
	G4RotationMatrix *pRotDipoleFringe1InContainer=new G4RotationMatrix();
	G4RotationMatrix *pRotDipoleFringe2InContainer=new G4RotationMatrix();
	
	pRotDipoleInContainer->rotateY(90*deg);
	//pRotDipoleFringeInContainer->rotateY(-90*deg);
	pRotDipoleFringe1InContainer->rotateY(180*deg);
	pRotDipoleFringe1InContainer->rotateZ(90*deg);
	pRotDipoleFringe1InContainer->rotateX(90*deg); 
	pRotDipoleFringe2InContainer->rotateX(90*deg);
	pRotDipoleFringe2InContainer->rotateY(90*deg + 45.*deg); 

	G4RotationMatrix *pRot_den=new G4RotationMatrix();
	G4RotationMatrix *pRot_dex=new G4RotationMatrix();	
	pRot_den->rotateY(90*deg);
	pRot_den->rotateX(-60*deg);
	pRot_dex->rotateY(90*deg);
	pRot_dex->rotateX(-90*deg + 105 * deg);

	//pRotDipoleFringe2InContainer->rotateX(90*deg); 
	//if(0)
	if(fSetupHRS>=4)
	{
	  //put it in the container, which also center at the hall center
            // SPR:  FIXME  This is problematic out of the box
	  //new G4PVPlacement(pRotDipoleInContainer,
	//		    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
	//		    LDipoleLogical,"LDipolePhys",LHRSContainerLogical,0,0,0);
	}
	//if(0)
	if(fSetupHRS>=4)
	{
            // SPR:  FIXME This is problematic out of the box
	 // new G4PVPlacement(pRotDipoleInContainer,
	//		    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
	//		    RDipoleLogical,"RDipolePhys",RHRSContainerLogical,0,0,0);
	}

	


	/////////////////////////
	// HRS Q3 
	/////////////////////////
	double pQ3Rin=30.0*cm;
	double pQ3Rout=75.0*cm;
	double pQ3Length;//=180*cm;//NIM
	//pQ3Length = ( IS_NIM == 1 ) ? 180.0 * cm : 182.68 * cm; 
	pQ3Length = 182.68 * cm; 
	double pQ3CenterY;//=pDipoleR*(1-cos(pDipoleBendAngle))+2.4*m*sin(pDipoleBendAngle);//NIM
	double pQ3CenterZ;//=9.96*m+pDipoleR*sin(pDipoleBendAngle)+2.4*m*cos(pDipoleBendAngle);//NIM
	//pQ3CenterY = ( IS_NIM == 1 ) ? pDipoleR * ( 1 - cos( pDipoleBendAngle ) ) + 2.4 * m * sin( pDipoleBendAngle ) : 3.5853101 * m  + pQ3Length / sqrt(2.) / 2.;
	//pQ3CenterZ = ( IS_NIM == 1 ) ? 9.96 * m + pDipoleR * sin( pDipoleBendAngle ) + 2.4 * m * cos( pDipoleBendAngle ) : 17.0257042 * m  + pQ3Length / sqrt(2) / 2;
	pQ3CenterY =  3.5853101 * m  + pQ3Length / sqrt(2.) / 2.;
	pQ3CenterZ = 17.0257042 * m  + pQ3Length / sqrt(2.) / 2;

	double pQ3LengthMag=pQ3Length * fringe_extension;
	//double pQ3LengthMag=182.68*cm;//SNAKE

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q3Solid = new G4Tubs("Q3Tub",pQ3Rin,pQ3Rout,pQ3Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ3Logical = new G4LogicalVolume(Q3Solid,
		mMaterialManager->siliconsteel,"LQ3Logical",0,0,0);
	G4LogicalVolume* RQ3Logical = new G4LogicalVolume(Q3Solid,
		mMaterialManager->siliconsteel,"RQ3Logical",0,0,0);

	LQ3Logical->SetVisAttributes(IronVisAtt); 
	RQ3Logical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ3InContainer=new G4RotationMatrix();
	pRotQ3InContainer->rotateX(-45*deg); 
	if(fSetupHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),
			LQ3Logical,"LQ3Phys",LHRSContainerLogical,0,0,0);
	}
	if(fSetupHRS>=4)
	{
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),
			RQ3Logical,"RQ3Phys",RHRSContainerLogical,0,0,0);
	}


	//vacuum pipe
	double pQ3vacRin  = 30.0*cm;
	double pQ3vacRout = 35.0*cm;
//	double pQ3vacLength = 575. * mm; 
	double pQ3vacLength = 530. * mm; // SPR:  Shortening this so that we don't overlap with fringe volume
	double pQ3vacCenterY = pQ3CenterY + pQ3Length / 2. / sqrt( 2 ) + pQ3vacLength / 2. / sqrt( 2 );
	double pQ3vacCenterZ = pQ3CenterZ + pQ3Length / 2. / sqrt( 2 ) + pQ3vacLength / 2. / sqrt( 2 );
	G4VSolid* Q3vacSolid = new G4Tubs("Q3vacTub",pQ3vacRin,pQ3vacRout,pQ3vacLength/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ3vacLogical = new G4LogicalVolume(Q3vacSolid,
		mMaterialManager->siliconsteel,"LQ3vacLogical",0,0,0);
	G4LogicalVolume* RQ3vacLogical = new G4LogicalVolume(Q3vacSolid,
		mMaterialManager->siliconsteel,"RQ3vacLogical",0,0,0);

	LQ3vacLogical->SetVisAttributes(IronVisAtt); 
	RQ3vacLogical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ3vacInContainer=new G4RotationMatrix();
	pRotQ3vacInContainer->rotateX(-45*deg); 
	if(fSetupHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3vacInContainer,G4ThreeVector(0,-pQ3vacCenterZ,pQ3vacCenterY),
			LQ3vacLogical,"LQ3vacPhys",LHRSContainerLogical,0,0,0);
	}
	if(fSetupHRS>=4)
	{
		new G4PVPlacement(pRotQ3vacInContainer,G4ThreeVector(0,-pQ3vacCenterZ,pQ3vacCenterY),
			RQ3vacLogical,"RQ3vacPhys",RHRSContainerLogical,0,0,0);
	}
	//vacuum pipe
	double pQ3vac2CenterY = pQ3CenterY - pQ3Length / 2. / sqrt( 2 ) - pQ3vacLength / 2. / sqrt( 2 );
	double pQ3vac2CenterZ = pQ3CenterZ - pQ3Length / 2. / sqrt( 2 ) - pQ3vacLength / 2. / sqrt( 2 );
	G4VSolid* Q3vac2Solid = new G4Tubs("Q3vac2Tub",pQ3vacRin,pQ3vacRout,pQ3vacLength/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ3vac2Logical = new G4LogicalVolume(Q3vac2Solid,
		mMaterialManager->siliconsteel,"LQ3vac2Logical",0,0,0);
	G4LogicalVolume* RQ3vac2Logical = new G4LogicalVolume(Q3vac2Solid,
		mMaterialManager->siliconsteel,"RQ3vac2Logical",0,0,0);

	LQ3vac2Logical->SetVisAttributes(IronVisAtt); 
	RQ3vac2Logical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ3vac2InContainer=new G4RotationMatrix();
	pRotQ3vac2InContainer->rotateX(-45*deg); 
	if(fSetupHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3vac2InContainer,G4ThreeVector(0,-pQ3vac2CenterZ,pQ3vac2CenterY),
			LQ3vac2Logical,"LQ3vac2Phys",LHRSContainerLogical,0,0,0);
	}
	if(fSetupHRS>=4)
	{
		new G4PVPlacement(pRotQ3vac2InContainer,G4ThreeVector(0,-pQ3vac2CenterZ,pQ3vac2CenterY),
			RQ3vac2Logical,"RQ3vac2Phys",RHRSContainerLogical,0,0,0);
	}



	//////////////////////////////////////////////////////////
	//#ifdef G4DEBUG_GEOMETRY
	//plot the tunnel to view
	//////////////////////////
	//Q1Q2 tunnel
	//////////////////////////

	const int kNPlane_Q1Q2Tunnel=4;
	double rInner_Q1Q2Tunnel[] = {0,0,0,0};
	double rOuter_Q1Q2Tunnel[] = {pQ1Rin-0.1*mm,pQ1Rin-0.1*mm,pQ2Rin-0.1*mm,pQ2Rin-0.1*mm};
	double zPlane_Q1Q2Tunnel[] = {pHRSContainerRin+6*mm,pHallCenter2RQ1Face+pQ1PosAct+10.0*cm,
		pHallCenter2RQ1Face+pQ1PosAct+30*cm,9.7*m};
	G4Polycone* Q1Q2TunnelSolid = new G4Polycone("Q1Q2TunnelPolycone",0.0,360.0*deg,
		kNPlane_Q1Q2Tunnel,zPlane_Q1Q2Tunnel,rInner_Q1Q2Tunnel,rOuter_Q1Q2Tunnel);

	G4LogicalVolume* LQ1Q2TunnelLogical = new G4LogicalVolume(Q1Q2TunnelSolid,
		mMaterialManager->vacuum,"LQ1Q2TunnelLogical",0,0,0);
	G4LogicalVolume* RQ1Q2TunnelLogical = new G4LogicalVolume(Q1Q2TunnelSolid,
		mMaterialManager->vacuum,"RQ1Q2TunnelLogical",0,0,0);

        G4VisAttributes *LightYellowVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.2)); 

	LQ1Q2TunnelLogical->SetVisAttributes(LightYellowVisAtt); 
	RQ1Q2TunnelLogical->SetVisAttributes(LightYellowVisAtt); 
	/*
	if(fSetupHRS>=3){//This is the old Jixie way
	  //put it in the container, which also center at the hall center
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,0,0),
			    LQ1Q2TunnelLogical,"LQ1Q2TunnelPhys",LHRSContainerLogical,0,0,0);
	}
	if(fSetupHRS>=3){//This is the old Jixie way
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,0,0),
			    RQ1Q2TunnelLogical,"RQ1Q2TunnelPhys",RHRSContainerLogical,0,0,0);
	}
	*/
 	//This is a copy of Q1, but of smaller radius.//This is the new Nickie way
	//G4VSolid* Q1SolidMag = new G4Tubs("Q1TubMag",0,pQ1Rin,pQ1Length/2.0,0.0,360.0*deg);//NIM
	G4VSolid* Q1SolidMag = new G4Tubs("Q1TubMag",0,pQ1Rin,pQ1LengthMag/2.0,0.0,360.0*deg);//for fringe fields, note longer length

	//build 2 copy since there are different fields involved in it
	G4LogicalVolume* LQ1MagLogical = new G4LogicalVolume(Q1SolidMag,
		mMaterialManager->vacuum,"LQ1MagLogical",0,0,0);
	G4LogicalVolume* RQ1MagLogical = new G4LogicalVolume(Q1SolidMag,
		mMaterialManager->vacuum,"RQ1MagLogical",0,0,0);

	//Nickie's Q1 field
	G4FieldManager* LQ1FieldManager = fEMFieldSetup->GetFieldManagerFZBL1();
	G4FieldManager* RQ1FieldManager = fEMFieldSetup->GetFieldManagerFZBR1();
	G4bool allLocal = true;
	//G4bool allLocal = false;
	LQ1MagLogical->SetFieldManager(LQ1FieldManager,allLocal);
	RQ1MagLogical->SetFieldManager(RQ1FieldManager,allLocal);

	LQ1MagLogical->SetVisAttributes(LightYellowVisAtt); 
	RQ1MagLogical->SetVisAttributes(LightYellowVisAtt); 


	//Nickie's old way
	if(fSetupHRS>=2){
	  //put it in the container, which also center at the hall center
	  //therefore only the z_at_lab position need to be considered
	  double pLQ1Pos_Z=(pHallCenter2LQ1Face+pQ1PosAct/2.0 + q1shift);//NIM
	  //double pLQ1Pos_Z=(pQ1PosAct + q1shift);//NIM
	  //double pLQ1Pos_Z=(pHallCenter2LQ1FaceMag+pQ1LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z,0),
			    LQ1MagLogical,"LQ1MagPhys",LHRSContainerLogical,0,0,0);
	}if(fSetupHRS>=2){
	  double pRQ1Pos_Z=(pHallCenter2RQ1Face+pQ1PosAct/2.0 + q1shift);//NIM
	  //double pRQ1Pos_Z=(pQ1PosAct + q1shift);//NIM
	  //double pRQ1Pos_Z=(pHallCenter2RQ1FaceMag+pQ1LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z,0),
			    RQ1MagLogical,"RQ1MagPhys",RHRSContainerLogical,0,0,0);
	}
	
	//G4VSolid* Q2SolidMag = new G4Tubs("Q2TubMag",0,pQ2Rin,pQ2Length/2.0,0.0,360.0*deg);//NIM
	G4VSolid* Q2SolidMag = new G4Tubs("Q2TubMag",0,pQ2Rin,pQ2LengthMag/2.0,0.0,360.0*deg);//SNAKE

	G4LogicalVolume* LQ2MagLogical = new G4LogicalVolume(Q2SolidMag,
		mMaterialManager->vacuum,"LQ2MagLogical",0,0,0);
	G4LogicalVolume* RQ2MagLogical = new G4LogicalVolume(Q2SolidMag,
		mMaterialManager->vacuum,"RQ2MagLogical",0,0,0);

	//Nickie's Q2 field
	G4FieldManager* LQ2FieldManager = fEMFieldSetup->GetFieldManagerFZBL2();
	G4FieldManager* RQ2FieldManager = fEMFieldSetup->GetFieldManagerFZBR2();
	
	LQ2MagLogical->SetFieldManager(LQ2FieldManager,allLocal);
	RQ2MagLogical->SetFieldManager(RQ2FieldManager,allLocal);

	LQ2MagLogical->SetVisAttributes(LightYellowVisAtt); 
	RQ2MagLogical->SetVisAttributes(LightYellowVisAtt); 

	 //Nickie's old way
	if(fSetupHRS>=3){
	  //put it in the container, which also center at the hall center
	  double pLQ2Pos_Z=(pHallCenter2LQ2Face+pQ2Length/2.0);//NIM
	  //double pLQ2Pos_Z=(pHallCenter2LQ2FaceMag+pQ2LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z,0),
			    LQ2MagLogical,"LQ2MagPhys",LHRSContainerLogical,0,0,0);
	}if(fSetupHRS>=3){
	  double pRQ2Pos_Z=(pHallCenter2RQ2Face+pQ2Length/2.0);//NIM
	  //double pRQ2Pos_Z=(pHallCenter2RQ2FaceMag+pQ2LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z,0),
			    RQ2MagLogical,"RQ2MagPhys",RHRSContainerLogical,0,0,0);
	}
	

	//////////////////////////
	//Dipole tunnel
	//////////////////////////
	//G4VSolid* DipoleTunnelSolid = new G4Tubs("DipoleTunnelSolid",pDipoleR-0.4*m+0.1*mm,
	//	pDipoleR+0.4*m-0.1*mm,0.125*m-0.1*mm,180*deg,pDipoleBendAngle);
	double pLQ1Pos_Z_en2=(pHallCenter2LQ1Face);//NIM
	double pLQ1Pos_Z_ex2=(pHallCenter2LQ1Face + pQ1PosAct);//NIM
	double pLQ2Pos_Z_en2=(pHallCenter2LQ2Face);//NIM
	double pLQ2Pos_Z_ex2=(pHallCenter2LQ2Face + pQ2Length);//NIM
	double pRQ1Pos_Z_en2=(pHallCenter2RQ1Face);//NIM
	double pRQ1Pos_Z_ex2=(pHallCenter2RQ1Face + pQ1PosAct);//NIM
	double pRQ2Pos_Z_en2=(pHallCenter2RQ2Face);//NIM
	double pRQ2Pos_Z_ex2=(pHallCenter2RQ2Face + pQ2Length);//NIM
	double pLDPos_Z_en2 = pLQ2Pos_Z_ex2 + 4.42 * m;
	double pRDPos_Z_en2 = pRQ2Pos_Z_ex2 + 4.42 * m;

	double rInner_DipoleTunnel2[] = {0.,0.,0.,0.};
	double rOuter_DipoleTunnel2[] = {pDipoleR-0.4*m,pDipoleR+0.4*m,pDipoleR+0.4*m,pDipoleR-0.4*m};

	G4Polycone* DipoleTunnelCone1 = new G4Polycone("DipoleTunnelCone1",
	            150.0*deg,105.0*deg,
		    kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel2,rOuter_DipoleTunnel2);
	G4Polycone* DipoleTunnelCone2 = new G4Polycone("DipoleTunnelCone2",
	            175.0*deg,pDipoleBendAngle + 10 * deg,
		    kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);

	G4Polycone* DipoleTunnelSolid = new G4Polycone("DipoleTunnelSolid",
		180.0*deg,pDipoleBendAngle,
		kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);

	G4double psiadj = 180. * deg - 22.5 * deg - 30. * deg;
	G4double badj = 8.4 * m * sin( 30.0 * deg ) / sin( psiadj );
	G4double xadj = badj * sin( 22.5 * deg );
	G4double yadj = badj * cos( 22.5 * deg );
	G4IntersectionSolid* DipoleTunnelCone3 = new G4IntersectionSolid("DipoleTunnelCone3",
			 DipoleTunnelCone2, DipoleTunnelCone1, 0, G4ThreeVector(-yadj, -xadj, 0));
	

	//DipoleTunnelSolid is Jixie with flat entrance
	//DipoleTunnelCone3 is Nickie with 30 deg entrance

	//But that is not good enough. We also have to add extensions to have fringe fields.
	//This is just going to be a trapezoidal extensionpRotDipoleInContaine
	//G4Trd* DipoleFringeSolid = new G4Trd("DipoleFringeSolid", 0.125*m+dy, 0.125*m-dy, 1.4 * m + 0.4 * m * 2. * tan( pi / 6 ), 1.4 * m, 0.4 * m);
	//G4cout << 0.125*m+dy << " and " << 0.125*m-dy << " better be greater than zero!" << G4endl;
	G4Trap* DipoleFringeSolid1 = new G4Trap("DipoleFringeSolid1",
					       0.4 * m , -pi / 12. ,
					       0.0     , 0.125*m-dy,
					       2.05 * m + 0.4 * m * 2. * tan( pi / 12. ), 2.05 * m + 0.4 * m * 2. * tan( pi / 12. ),
					       0.0     , 0.125*m+dy,
					       2.05 * m, 2.05 * m  ,
					       0.0   );
	G4Trap* DipoleFringeSolid2 = new G4Trap("DipoleFringeSolid2",
					       0.4 * m , -pi / 12. ,
					       0.0     , 0.125*m-dy,
					       0.65 * m + 0.4 * m * 2. * tan( pi / 12. ), 0.65 * m + 0.4 * m * 2. * tan( pi / 12. ),
					       0.0     , 0.125*m+dy,
					       0.65 * m, 0.65 * m  ,
					       0.0    );
	
	//pDz  , pTheta,
	//pPhi , pDy1,
	//pDx1 , pDx2,
	//pAlp1, pDy2,
	//pDx3 , pDx4,
	//pAlp2,

	//G4UnionSolid* DipoleTunnelCone4 = new G4UnionSolid("DipoleTunnelCone4",
	//DipoleTunnelCone3, DipoleFringeSolid, pRotDipoleInContainer, G4ThreeVector(-8.4 * m, 0., 0.));
	G4UnionSolid* DipoleTunnelCone4 = new G4UnionSolid("DipoleTunnelCone4",
		   DipoleTunnelCone3, DipoleFringeSolid1, pRotDipoleFringe1InContainer, G4ThreeVector(-8.4 * m, +2.05 * m, 0.));
	G4UnionSolid* DipoleTunnelCone5 = new G4UnionSolid("DipoleTunnelCone4",
		   DipoleTunnelCone4, DipoleFringeSolid2, pRotDipoleFringe2InContainer,
							   G4ThreeVector(-8.4 * m * cos( pi / 4. ) + 0.3 * m / sqrt( 2. ),
									 -8.4 * m * sin( pi / 4. ) - 0.3 * m / sqrt( 2. ), 0.));

	G4LogicalVolume* LDipoleFringe1Logical = new G4LogicalVolume(DipoleFringeSolid1,
	        mMaterialManager->vacuum,"LDipoleFringe1Logical",0,0,0);
	G4LogicalVolume* RDipoleFringe1Logical = new G4LogicalVolume(DipoleFringeSolid1,
	        mMaterialManager->vacuum,"RDipoleFringe1Logical",0,0,0);
	G4LogicalVolume* LDipoleFringe2Logical = new G4LogicalVolume(DipoleFringeSolid2,
	        mMaterialManager->vacuum,"LDipoleFringe2Logical",0,0,0);
	G4LogicalVolume* RDipoleFringe2Logical = new G4LogicalVolume(DipoleFringeSolid2,
	        mMaterialManager->vacuum,"RDipoleFringe2Logical",0,0,0);

	G4LogicalVolume* LDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelCone5,
		mMaterialManager->vacuum,"LDipoleTunnelLogical",0,0,0);
	G4LogicalVolume* RDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelCone5,
		mMaterialManager->vacuum,"RDipoleTunnelLogical",0,0,0);
        /*
           // This utilizes Nickie's which are too short 
	G4LogicalVolume* LDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelCone3,
		mMaterialManager->vacuum,"LDipoleTunnelLogical",0,0,0);
	G4LogicalVolume* RDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelCone3,
		mMaterialManager->vacuum,"RDipoleTunnelLogical",0,0,0);
        */

        G4VisAttributes *YellowVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));      //yellow, color for Ultem

	LDipoleTunnelLogical->SetVisAttributes(YellowVisAtt); 
	RDipoleTunnelLogical->SetVisAttributes(YellowVisAtt); 
	LDipoleFringe1Logical->SetVisAttributes(YellowVisAtt);
	RDipoleFringe1Logical->SetVisAttributes(YellowVisAtt); 
	LDipoleFringe2Logical->SetVisAttributes(YellowVisAtt);
	RDipoleFringe2Logical->SetVisAttributes(YellowVisAtt); 
	//Nickie's dipole field
	G4FieldManager* LdipoleFieldManager = fEMFieldSetup->GetFieldManagerFZBL3();
	G4FieldManager* RdipoleFieldManager = fEMFieldSetup->GetFieldManagerFZBR3();
	//G4double minEps= 1.0e-9;  //   Minimum & value for smallest steps
	//G4double maxEps= 0.25*mm;  //   Maximum & value for largest steps
	
	//dipoleFieldManager->SetMinimumEpsilonStep( minEps );
	//LdipoleFieldManager->SetMaximumEpsilonStep( maxEps );
	//RdipoleFieldManager->SetMaximumEpsilonStep( maxEps );
	//dipoleFieldManager->SetDeltaOneStep( 0.5e-9 * mm );
	//dipoleFieldManager->GetChordFinder()->SetDeltaChord( 0.001 * mm );
	LDipoleTunnelLogical->SetFieldManager(LdipoleFieldManager,allLocal);
	RDipoleTunnelLogical->SetFieldManager(RdipoleFieldManager,allLocal);

        //double pFringeX_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) + 0.75 * m / sqrt( 2 ) :   2.4603032 * m + 0.75 * m / sqrt( 2 );
        //double pFringeZ_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) - 0.75 * m / sqrt( 2 ) : -15.9006973 * m - 0.75 * m / sqrt( 2 );

	if(fSetupHRS>=4)
	{
	  //put it in the container, which also center at the hall center
	  new G4PVPlacement(pRotDipoleInContainer,
			    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			    LDipoleTunnelLogical,"LDipoleTunnelPhys",LHRSContainerLogical,0,0,0);
	  /*
	  new G4PVPlacement(pRotDipoleFringe1InContainer,
			    G4ThreeVector(0., -7.75 * m, 0.),
			    LDipoleFringe1Logical,"LDipoleFringe1Phys",LHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotDipoleFringe2InContainer,
			    G4ThreeVector(0., pFringeZ_ex, pFringeX_ex),
			    LDipoleFringe2Logical,"LDipoleFringe2Phys",LHRSContainerLogical,0,0,0);
	  */
	}
	if(fSetupHRS>=4)
	{
	  new G4PVPlacement(pRotDipoleInContainer,
			    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			    RDipoleTunnelLogical,"RDipoleTunnelPhys",RHRSContainerLogical,0,0,0);
	  /*
	  new G4PVPlacement(pRotDipoleFringe1InContainer,
			    G4ThreeVector(0., -7.75 * m, 0.),
			    RDipoleFringe1Logical,"RDipoleFringe1Phys",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotDipoleFringe2InContainer,
			    G4ThreeVector(0., pFringeZ_ex, pFringeX_ex),
			    RDipoleFringe2Logical,"RDipoleFringe2Phys",RHRSContainerLogical,0,0,0);

	  */
	}


	//////////////////////////
	//Q3 tunnel
	//////////////////////////

	G4VSolid* Q3TunnelSolid = new G4Tubs("Q3TunnelTub",0,pQ3Rin-0.1*mm,
		pQ3Length,0.0,360.0*deg);
	G4LogicalVolume* LQ3TunnelLogical = new G4LogicalVolume(Q3TunnelSolid,
		mMaterialManager->vacuum,"LQ3TunnelLogical",0,0,0);
	G4LogicalVolume* RQ3TunnelLogical = new G4LogicalVolume(Q3TunnelSolid,
		mMaterialManager->vacuum,"RQ3TunnelLogical",0,0,0);

	LQ3TunnelLogical->SetVisAttributes(YellowVisAtt); 
	RQ3TunnelLogical->SetVisAttributes(YellowVisAtt); 

	//double pQ3TunnelPos_Y=pDipoleR*(1-cos(pDipoleBendAngle))+2.0*m*sin(pDipoleBendAngle);
	//double pQ3TunnelPos_Z=9.96*m+pDipoleR*sin(pDipoleBendAngle)+2.0*m*cos(pDipoleBendAngle);
	//double pQ3TunnelPos_Z=pDipoleR*sin(pDipoleBendAngle)+2.0*m*cos(pDipoleBendAngle);
	//pQ3TunnelPos_Z += ( IS_NIM == 1 ) ? 9.96*m : pDipoleRCenterZ; 
	//if(fSetupHRS>=4)
	//{
		//put it in the container, which also center at the hall center
		//new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3TunnelPos_Z,pQ3TunnelPos_Y),
				  //LQ3TunnelLogical,"LQ1Q2TunnelPhys",LHRSContainerLogical,0,0,0);
	//}
	//if(fSetupHRS>=4)
	//{
	  //new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3TunnelPos_Z,pQ3TunnelPos_Y),
	  //RQ3TunnelLogical,"RQ1Q2TunnelPhys",RHRSContainerLogical,0,0,0);
	//}

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	//G4VSolid* Q3MagSolid = new G4Tubs("Q3MagTub",0,pQ3Rin,pQ3Length/2.0,0.0,360.0*deg);//NIM
	G4VSolid* Q3MagSolid = new G4Tubs("Q3MagTub",0,pQ3Rin,pQ3LengthMag/2.0,0.0,360.0*deg);//SNAKE

	G4LogicalVolume* LQ3MagLogical = new G4LogicalVolume(Q3MagSolid,
		mMaterialManager->vacuum,"LQ3MagLogical",0,0,0);
	G4LogicalVolume* RQ3MagLogical = new G4LogicalVolume(Q3MagSolid,
		mMaterialManager->vacuum,"RQ3MagLogical",0,0,0);

	//Nickie's Q3 field
	G4FieldManager* LQ3FieldManager = fEMFieldSetup->GetFieldManagerFZBL4();
	G4FieldManager* RQ3FieldManager = fEMFieldSetup->GetFieldManagerFZBR4();
	LQ3MagLogical->SetFieldManager(LQ3FieldManager,allLocal);
	RQ3MagLogical->SetFieldManager(RQ3FieldManager,allLocal);

	LQ3MagLogical->SetVisAttributes(MagFieldVisAtt); 
	RQ3MagLogical->SetVisAttributes(MagFieldVisAtt); 

	G4RotationMatrix *pRotQ3MagInContainer=new G4RotationMatrix();
	pRotQ3MagInContainer->rotateX(-45*deg); 
	if(fSetupHRS>=4)
	{
	  //put it in the container, which also center at the hall center
	  //new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZMag,pQ3CenterYMag),//SNAKE
	  //LQ3MagLogical,"LQ3MagPhys",LHRSContainerLogical,0,0,0);//SNAKE
	  new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),//NIM
			    LQ3MagLogical,"LQ3MagPhys",LHRSContainerLogical,0,0,0);//NIM
	}
	if(fSetupHRS>=4)
	{
	  //new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZMag,pQ3CenterYMag),//SNAKE
	  //RQ3MagLogical,"RQ3MagPhys",RHRSContainerLogical,0,0,0);//SNAKE
	  new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),//NIM
	    		    RQ3MagLogical,"RQ3MagPhys",RHRSContainerLogical,0,0,0);//NIM
	}


	//Nickie's sensitive detector, at the focal plane, virtual boundary, vb, fp
	//double pFPR      = pDipoleR; //radius of curvature of dipole
	//double pFPA      = 9.96 * m; //distance from pivot to dipole entrance
	//double pFPB      = 2.4  * m; //distance from center of dipole to center of Q3
	//
	//double pFPB      = 2.4  * m + pQ3Length / 2.0 + 3.57 * m + 1.43 * m;
	//double pFPB      = 2.4  * m + pQ3Length / 2.0 + 3.57 * m + 0.7 * m;
	//double pFPCenterY=       pFPR * ( 1 - cos( pDipoleBendAngle ) ) + pFPB * sin( pDipoleBendAngle );
	//double pFPCenterZ=pFPA + pFPR *       sin( pDipoleBendAngle )   + pFPB * cos( pDipoleBendAngle );
        //double pFPA      = 9.96 * m;//distance from pivot to entrance of dipole
        //double pFPCenterX=( pFPA + pFPH + ( pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2) ) *-sin( mRHRSAngle );
        //double pFPCenterZ=( pFPA + pFPH + ( pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2) ) * cos( mRHRSAngle );

	double pFPCenterZ , pFPCenterY ;
	double pVDCCenterZ, pVDCCenterY;
	double pQZ1CenterZ, pQZ1CenterY;
	double pQZ2CenterZ, pQZ2CenterY;
	//if( IS_NIM == 1 ){
	//double pFPR      = 8.4  * m;//radius of curvature of dipole
	//double pFPA      = ( IS_NIM == 1 ) ? 9.96 * m : pDipoleRCenterZ;//distance from pivot to entrance of dipole
	//double pFPH      = pFPR * tan ( 45. / 2. * deg );
	//pFPCenterZ=( pFPA + pFPH + ( pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2) );
	//pFPCenterY=(                 pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2);
	//}else{
	pVDCCenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) ;
	pVDCCenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) ;
	pQZ1CenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 40.0 * cm;
	pQZ1CenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 40.0 * cm;
	pQZ2CenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 52.0 * cm;
	pQZ2CenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 52.0 * cm;
	pFPCenterZ  = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m + 1.43 * m ) / sqrt(2) ;
	pFPCenterY  = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m + 1.43 * m ) / sqrt(2) ;
	//}
	
	G4double vb_thickness = .00001 * mm;
	G4VSolid* FPSolid = new G4Tubs("FPTub",0,pQ3Rout * 2,vb_thickness,0.0,360.0*deg);
	G4VSolid* PlaneSolid1 = new G4Tubs("PlaneTub",0,pQ1Rout,vb_thickness,0.0,360.0*deg); //circles
	G4VSolid* PlaneSolid2 = new G4Tubs("PlaneTub",0,pQ2Rout,vb_thickness,0.0,360.0*deg); //circles

	double X_S0 = 29.5 * cm;
	double Y_S0 = 35.5 * cm;
	double Z_S0 = vb_thickness;

	double X_S1 = 29.5 * cm;
	double Y_S1 = 35.5 * cm;
	double Z_S1 = vb_thickness;

	double X_Q1 = 29.5 * cm;
	double Y_Q1 = 35.5 * cm;
	double Z_Q1 = vb_thickness;

	double X_Q2 = 29.5 * cm;
	double Y_Q2 = 35.5 * cm;
	double Z_Q2 = vb_thickness;

	G4VSolid* Scint0  = new G4Box("Scint0" , X_S0 / 2.0, Y_S0 / 2.0, Z_S0 / 2.0);
	G4VSolid* Scint1  = new G4Box("Scint1" , X_S1 / 2.0, Y_S1 / 2.0, Z_S1 / 2.0);
	G4VSolid* Quartz1 = new G4Box("Quartz1", X_Q1 / 2.0, Y_Q1 / 2.0, Z_Q1 / 2.0);
	G4VSolid* Quartz2 = new G4Box("Quartz2", X_Q2 / 2.0, Y_Q2 / 2.0, Z_Q2 / 2.0);

	G4VSolid* magneticSolid = new G4Box("magneticBox",fFieldX/2.0,fFieldY/2.0,fFieldZ/2.0);

	G4LogicalVolume* LFPLogical = new G4LogicalVolume(FPSolid,
		mMaterialManager->vacuum,"LFPLogical",0,0,0);
	G4LogicalVolume* RFPLogical = new G4LogicalVolume(FPSolid,
		mMaterialManager->vacuum,"RFPLogical",0,0,0);

	G4LogicalVolume* LPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,
		mMaterialManager->vacuum,"LPlaneLogical1",0,0,0);
	G4LogicalVolume* RPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,
		mMaterialManager->vacuum,"RPlaneLogical1",0,0,0);
	G4LogicalVolume* LPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,
		mMaterialManager->vacuum,"LPlaneLogical1",0,0,0);
	G4LogicalVolume* RPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,
		mMaterialManager->vacuum,"RPlaneLogical1",0,0,0);

	LFPLogical->SetVisAttributes(MagFieldVisAtt); 
	RFPLogical->SetVisAttributes(MagFieldVisAtt); 
	LPlaneLogical1->SetVisAttributes(MagFieldVisAtt); 
	RPlaneLogical1->SetVisAttributes(MagFieldVisAtt); 
	LPlaneLogical2->SetVisAttributes(MagFieldVisAtt); 
	RPlaneLogical2->SetVisAttributes(MagFieldVisAtt); 

	G4RotationMatrix *pRotVDCInContainer=new G4RotationMatrix();
	pRotVDCInContainer->rotateX(0.*deg); 
	G4RotationMatrix *pRotFPInContainer=new G4RotationMatrix();
	pRotFPInContainer->rotateX(-45*deg); 
	if(fSnakeModel == 49 || fSnakeModel == 48 || fSnakeModel > 51 ){
	  double pSeptumX      = 140.0  * cm;
	  double pSeptumY      = 84.4   * cm;
	  double pSeptumZ      = 74.0   * cm;
	  //double pSeptumPlaceZ = 70.414 * cm;
	  double pSeptumPlaceZ = 69.99937 * cm;
	  
	}
	//double pLQ1Pos_Z_en=(pHallCenter2LQ1Face);//NIM
	//double pLQ1Pos_Z_ex=(pHallCenter2LQ1Face + pQ1Length);//NIM
	//double pLQ2Pos_Z_en=(pHallCenter2LQ2Face);//NIM
	//double pLQ2Pos_Z_ex=(pHallCenter2LQ2Face + pQ2Length);//NIM
	//double pRQ1Pos_Z_en=(pHallCenter2RQ1Face);//NIM
	//double pRQ1Pos_Z_ex=(pHallCenter2RQ1Face + pQ1Length);//NIM
	//double pRQ2Pos_Z_en=(pHallCenter2RQ2Face);//NIM
	//double pRQ2Pos_Z_ex=(pHallCenter2RQ2Face + pQ2Length);//NIM
	//double pLDPos_Z_en = pLQ2Pos_Z_ex + 4.42 * m;
	//double pRDPos_Z_en = pRQ2Pos_Z_ex + 4.42 * m;
	//double pLDPos_X_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) : 2.4603032 * m;
	//double pRDPos_X_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) : 2.4603032 * m;
	//double pLDPos_Z_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) : -15.9006973 * m;
	//double pRDPos_Z_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) : -15.9006973 * m;

	double pHallCenter2Col = 1.38 * m;
	double pPaulColT        = 0.01 * m;
	double pPaulX = ( - pHallCenter2Col - pPaulColT * 2. ) * cos(fHRSAngle) ;
	double pPaulY = ( - pHallCenter2Col - pPaulColT * 2. ) * sin(fHRSAngle);
	if(fSnakeModel == 49 || fSnakeModel == 48 || fSnakeModel > 50 ){
	  //double pLQ1Pos_Z=(pHallCenter2LQ1FaceMag+pQ1LengthMag/1.0);//SNAKE

	}
	if(fSnakeModel == 49 || fSnakeModel == 48 || fSnakeModel > 50 ){
	}

	//#endif
	//////////////////////////////////////////////////////////

	return;

}


g4hrsEMField* g4hrsDetectorConstruction::GetEMFieldFromSetup(){
	
	return fEMFieldSetup->GetEMField();
	
}



