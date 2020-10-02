// ********************************************************************
//
// $Id: g4hrsMaterial.cc,v 1.0, 2010/12/26 HRS Exp $
//This class must be launched from HRSDetectorConstruction and should
//be instanced only once
//
// ********************************************************************
#include <stdio.h>
#include <math.h>
#include <string>
using namespace std;

#include "g4hrsMaterial.hh"


//#define G4DEBUG_MATERIAL 1

/////////////////////////////////////////////////////////////////////

g4hrsMaterial* g4hrsMaterial::fHRSMaterialManager=0;
g4hrsMaterial* g4hrsMaterial::GetHRSMaterialManager()
{ 
	if(!fHRSMaterialManager) 
	{
            fHRSMaterialManager = new g4hrsMaterial();
	  G4cout<<"Error! No instance of g4hrsMaterial found! \n"
			<<"Use g4hrsMaterial::g4hrsMaterial() to build one in HRSDetectorConstruction::Construct()!"
			<<G4endl;
	}
	return fHRSMaterialManager; 
}

/////////////////////////////////////////////////////////////////////

g4hrsMaterial::g4hrsMaterial()
{
	
	//default 55% NH3 and 45% He4 mixture. In volumn ratio 
	mSolidNH3D=817*mg/cm3;
	mLiquidHeD=145*mg/cm3;
	mNH3WeightRatio=0.8732;

	ConstructMaterials();
	fHRSMaterialManager=this;
	//std::cout<<"g4hrsMaterial() construction done!"<<std::endl;
}

g4hrsMaterial::~g4hrsMaterial()
{
	DestroyMaterials(); 
	//std::cout<<"delete g4hrsMaterial ...done!"<<std::endl;
}


/////////////////////////////////////////////////////////////////////
void g4hrsMaterial::ConstructMaterials()
{
	G4double a;
	G4double z;
	G4double density;
	G4double weightRatio;
	G4String name;
	G4String symbol;
	G4int nElem,nComponent,nAtoms;
	G4double pressure;
	G4double temperature;

	G4NistManager* pMatMan = G4NistManager::Instance();

	// elements for mixtures and compounds
	/*
	a = 1.01*g/mole;
	G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);
	a = 9.012*g/mole;
	G4Element* elBe = new G4Element(name="Beryllium", symbol="Be", z=4., a);
	a = 12.01*g/mole;
	G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);
	a = 14.01*g/mole;
	G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7., a);
	a = 16.00*g/mole;
	G4Element* elO = new G4Element(name="Oxigen", symbol="O", z=8., a);
	a=18.9984032*g/mole;
	G4Element* elF = new G4Element(name="Fluorine", symbol="F", z=9., a);
	a = 28.09*g/mole;
	G4Element* elSi= new G4Element(name="Silicon", symbol="Si", z=14., a);
	a = 35.45*g/mole;
	G4Element* elCl= new G4Element(name="Chlorine", symbol="Cl", z=17., a);
	a = 51.9961*g/mole;
	G4Element* elCr= new G4Element(name="Chromium", symbol="Cr", z=24., a);
	a = 54.94*g/mole;
	G4Element* elMn = new G4Element(name="Manganese", symbol="Mn", z=25., a);
	a = 55.845*g/mole;
	G4Element* elFe= new G4Element(name="Iron", symbol="Fe", z=26., a);
	a = 58.6934*g/mole;
	G4Element* elNi= new G4Element(name="Nickel", symbol="Ni", z=28., a);
	a = 63.55*g/mole;
	G4Element* elCu= new G4Element(name="Copper", symbol="Cu", z=29., a);
	a = 95.96*g/mole;
	G4Element* elMo = new G4Element(name="Molydb", symbol="Mo", z=42., a);
	*/
	int zz=0;
	G4Element* elH  = pMatMan->FindOrBuildElement(zz=1);
	G4Element* elBe = pMatMan->FindOrBuildElement(zz=4);
	G4Element* elC  = pMatMan->FindOrBuildElement(zz=6);
	G4Element* elN  = pMatMan->FindOrBuildElement(zz=7);
	G4Element* elO  = pMatMan->FindOrBuildElement(zz=8);
	G4Element* elF  = pMatMan->FindOrBuildElement(zz=9);
	G4Element* elAl = pMatMan->FindOrBuildElement(zz=13);
	G4Element* elSi = pMatMan->FindOrBuildElement(zz=14);
	G4Element* elP  = pMatMan->FindOrBuildElement(zz=15);
	G4Element* elS  = pMatMan->FindOrBuildElement(zz=16);
	G4Element* elCl = pMatMan->FindOrBuildElement(zz=17);
	//G4Element* elCa = pMatMan->FindOrBuildElement(zz=20);
	G4Element* elCr = pMatMan->FindOrBuildElement(zz=24);
	G4Element* elMn = pMatMan->FindOrBuildElement(zz=25);
	G4Element* elFe = pMatMan->FindOrBuildElement(zz=26);
	G4Element* elCo = pMatMan->FindOrBuildElement(zz=27);
	G4Element* elNi = pMatMan->FindOrBuildElement(zz=28);
	G4Element* elCu = pMatMan->FindOrBuildElement(zz=29);
	G4Element* elMo = pMatMan->FindOrBuildElement(zz=42);


	//outspace vacuum
	density = universe_mean_density;  //1.e-25*g/cm3;
	pressure = 3.0e-18 *pascal;
	temperature = 2.7 *kelvin;  //cosmic background temperature
	a = 1.3e-7 *g/mole;
	galaxy = new G4Material(name="Galaxy",z=1,a,density,kStateGas,temperature,pressure);

	//A pressure that causes the Hg column to rise 1 millimeter is called ONE torr (1mm Hg)
	// Vaccum of 1.e-6 torr at room temperature,  1 atmosphere = 101325*pascal = 760 *torr
	density = 1.e-6/760.0 * 1.29*mg/cm3; //0.001 of air density
	pressure = 1.e-6/760.0 *atmosphere;
	temperature = 293.15 *kelvin;  //room temperature
	a = 28.97 *g/mole;
	vacuum = new G4Material(name="Vacuum",z=1,a,density,kStateGas,temperature,pressure);

	//G4cout << "Vacuum: " << z << " " << density << " " << kStateGas << " " << temperature << " " << pressure << G4endl;
	//Vacuum: 1 1.05941e+07 3 293.15 0.832133
	// Air
	density = 1.29*mg/cm3;
	air = new G4Material(name="Air", density, nElem=2);
	air->AddElement(elN, weightRatio=0.7);
	air->AddElement(elO, weightRatio=0.3);
/*
	// Argon gas
	a = 39.95*g/mole;
	density = 1.782*mg/cm3;
	argonGas = new G4Material(name="ArgonGas",z=18.,a,density,kStateGas,temperature,pressure);

	//kapton
	density = 1.42*g/cm3;
	kapton = new G4Material(name="Kapton", density, nElem=4);
	kapton->AddElement(elH, nAtoms=10);
	kapton->AddElement(elC, nAtoms=22);
	kapton->AddElement(elO, nAtoms=5);
	kapton->AddElement(elN, nAtoms=2);

	/////////////////////////////////////////////////////////
	//ultem C37H24O6N2
	density = 1.27*g/cm3;
	ultem = new G4Material(name="Ultem", density, nElem=4);
	ultem->AddElement(elH, nAtoms=37);
	ultem->AddElement(elC, nAtoms=24);
	ultem->AddElement(elO, nAtoms=6);
	ultem->AddElement(elN, nAtoms=2);

	//fused quartz SiO2
	density = 2.20*g/cm3;
	SiO2 = new G4Material(name="Fused-Quartz", density, nElem=2);
	SiO2->AddElement(elSi, nAtoms=1);
	SiO2->AddElement(elO, nAtoms=2);

	//epoxy-resin , C11H12O3
	//density = 0.95*g/cm3;     //volumn ratio 60%:40%
	density = 1.268*g/cm3;  //weight ratio 60%:40%
	epoxy = new G4Material(name="Epoxy-Resin", density, nElem=3);
	epoxy->AddElement(elH, nAtoms=12);
	epoxy->AddElement(elC, nAtoms=11);
	epoxy->AddElement(elO, nAtoms=3);

	//G10FR4 is 60% of SiO2 and 40% epoxy
	density = 1.7*g/cm3;
	G10FR4 = new G4Material(name="G10FR4", density, nComponent=2);
	G10FR4->AddMaterial(SiO2, weightRatio=0.6);
	G10FR4->AddMaterial(epoxy, weightRatio=0.4);

	//CH2: 0.94g/cm^3 ,H-[CH2]n-H
	CH2 = new G4Material(name="CH2", density=0.94*g/cm3, nElem=2);
	CH2->AddElement(elH, nAtoms=2);
	CH2->AddElement(elC, nAtoms=1);

	//scintillator
	density = 1.032*g/cm3;
	scintillator = new G4Material(name="EJ204", density, nElem=2);
	scintillator->AddElement(elH, nAtoms=521);
	scintillator->AddElement(elC, nAtoms=474);

	//pcbNchip Cu9Si5C460H506O138 epoxy + silicon + copper
	//weight ratio 57.7:4:1
	density = 2.31*g/cm3;
	pcbNchip = new G4Material(name="Pcb&Chip", density, nElem=5);
	pcbNchip->AddElement(elH, nAtoms=506);
	pcbNchip->AddElement(elC, nAtoms=460);
	pcbNchip->AddElement(elO, nAtoms=138);
	pcbNchip->AddElement(elSi, nAtoms=5);
	pcbNchip->AddElement(elCu, nAtoms=9);

	//cable Cu4C6H9Cl3
	density = 30.27/(2.15*0.04*100.) *g/cm3;
	cable = new G4Material(name="Cable", density, nElem=4);
	cable->AddElement(elH, nAtoms=9);
	cable->AddElement(elC, nAtoms=6);
	cable->AddElement(elCl, nAtoms=3);
	cable->AddElement(elCu, nAtoms=4);

	/////////////////////////////////////////////////////////
	//at 0 degree, 1 atm, density = 0.1786 mg/cm^3
	//Helium Gas at 1 atm & room temprature, density=0.163 *mg/cm3;
	a = 4.0026*g/mole;
	pressure = 1.0 *atmosphere ;
	temperature = 293.0 *kelvin ;
	density = 0.163 *mg/cm3 ;
	density = (pressure/temperature)/(1.0*atmosphere/273*kelvin) * 0.1786 * mg/cm3;
	heliumGas = new G4Material(name="HeliumGas", z=2., a, density, kStateGas, temperature, pressure);

	//Helium Gas at 1 atm & 4.22k temprature, 
	a = 4.0026*g/mole;
	pressure = 1.0 *atmosphere ;
	temperature = 4.22 *kelvin ;
	density = (pressure/temperature)/(1.0*atmosphere/273*kelvin) * 0.1786 * mg/cm3;
	heliumGas4k = new G4Material(name="HeliumGas4k", z=2., a, density, kStateGas, temperature, pressure);

	//double mLiquidHeD=0.145*mg/cm3;
	//liquid He 0.145 *g/cm3
	a = 4.0026*g/mole;
	density = mLiquidHeD;
	liquidHe = new G4Material(name="LiquidHe", z=2., a, density);

	//liquid He3 0.08243 *g/cm3
	a = 30.069*g/mole;
	density = 0.08243 *g/cm3;
	liquidHe3 = new G4Material(name="LiquidHe3", z=2., a, density);

	//Ethane gas 1.252 *mg/cm3
	a = 3.01603*g/mole;
	pressure = 1.0 *atmosphere ;
	temperature = 293.0 *kelvin ;
	density = 1.252 *mg/cm3;
	ethane = new G4Material(name="ethane", density, 2, kStateGas, temperature, pressure);
	ethane -> AddElement(elH,6);
	ethane -> AddElement(elC,2);

	//bigbite wire chamber gas is a mixture of He4 and ethane 
	EthaneHe = new G4Material("Ethane+He",0.00056*g/cm3,2,
		kStateGas,temperature,temperature);
	EthaneHe -> AddMaterial(heliumGas,0.21);		// % in relative mass
	EthaneHe -> AddMaterial(ethane,0.79);			// % in relative mass

	//liquid H2 0.07099 *g/cm3
	a = 1.00794*g/mole;
	density = 0.07099 *g/cm3;
	liquidH2 = new G4Material(name="LiquidH2", z=1., a, density);

	//liquid D2 0.180 *g/cm3
	a = 2.01410*g/mole;
	density = 0.180 *g/cm3;
	liquidD2 = new G4Material(name="LiquidD2", z=1., a, density);
*/

	//Define diamond
	a = 12.01 * g/mole;
	density = 3.515*g/cm3;
	diamond = new G4Material("Diamond", z=6, a, density);

        //mylar
	density = 1.39*g/cm3;
	mylar = new G4Material(name="Mylar", density, nElem=2);
	mylar->AddElement(elH, nAtoms=8);
	mylar->AddElement(elC, nAtoms=10);

	//Aluminum
	a = 26.982*g/mole;
	density = 2.70*g/cm3;
	aluminum = new G4Material(name="Aluminum", z=13., a, density);

	//carbon12
	a = 12.01078*g/mole;
	density = 2.267*g/cm3;
	carbon12 = new G4Material(name="Carbon12", z=6., a, density);

	//Calcium40
	a = 40.078*g/mole;//natural
	density = 1.55*g/cm3;//natural
	calcium40 = new G4Material(name="Calcium40", z=20., a, density);
	
	//Calcium48
	a = 47.952534*g/mole;//48Ca
	density = 1.855*g/cm3;//48Ca
	calcium48 = new G4Material(name="Calcium48", z=20., a, density);
/*	
	//Iron
	a = 55.845*g/mole;
	density = 7.874*g/cm3;
	iron = new G4Material(name="Iron", z=26., a, density);

	//StainlessSteel, Material Names : stainless steel
	//Material : Fe-Cr-Ni-Mo,  93xx: 3.25% Ni,1.2%Cr,0.12%Mo,
	density = 7.85 *g/cm3;
	stainlesssteel = new G4Material(name="StainlessSteel", density, nElem=4);
	stainlesssteel->AddElement(elFe, weightRatio=0.9543);
	stainlesssteel->AddElement(elCr, weightRatio=0.0325);
	stainlesssteel->AddElement(elNi, weightRatio=0.0120);
	stainlesssteel->AddElement(elMo, weightRatio=0.0012);

	//304 stainless stell http://baike.baidu.com/view/973507.htm
	//  C      Si     Mn     P     S       Cr           Ni          Mo
	//<0.08   <1.00  <2.00 <0.05 <0.03 18.00-20.00   8.00-10.50    2.84
	//code [06Cr19Ni10], which means: other 6%, Cr 19%, Ni 10% 
	density = 7.93 *g/cm3;
	stainlesssteel304 = new G4Material(name="StainlessSteel304", density, nElem=9);
	stainlesssteel304->AddElement(elFe, weightRatio=0.65);
	stainlesssteel304->AddElement(elCr, weightRatio=0.19);
	stainlesssteel304->AddElement(elNi, weightRatio=0.10);
	stainlesssteel304->AddElement(elMn, weightRatio=0.02);
	stainlesssteel304->AddElement(elSi, weightRatio=0.01);
	stainlesssteel304->AddElement(elMo, weightRatio=0.0284);
	stainlesssteel304->AddElement(elC,  weightRatio=0.0008);
	stainlesssteel304->AddElement(elP,  weightRatio=0.0005);
	stainlesssteel304->AddElement(elS,  weightRatio=0.0003);
*/
	//Silicon Steel(electrical steel),
	//Electrical steel is an iron alloy which may have from zero to 6.5% silicon (Si:5Fe). 
	//Silicon significantly increases the electrical resistivity of the steel, which 
	//decreases the induced eddy currents and thus reduces the core loss. Manganese and 
	//aluminum can be added up to 0.5%.

	//There are two main types of electrical steel: grain-oriented and non-oriented.
	//Grain-oriented electrical steel usually has a silicon level of 3% (Si:11Fe). It is 
	//processed in such a way that the optimum properties are developed in the rolling 
	//direction, due to a tight control (proposed by Norman P. Goss) of the crystal 
	//orientation relative to the sheet. Due to the special orientation, the magnetic 
	//flux density is increased by 30% in the coil rolling direction, although its magnetic 
	//saturation is decreased by 5%. It is used for the cores of high-efficiency transformers, 
	//electric motor and generators. Cold Rolled Grain-oriented steel is often abbreviated to CRGO.
	//
	//Non-oriented electrical steel usually has a silicon level of 2 to 3.5% and has similar 
	//magnetic properties in all directions, which makes it isotropic. It is less expensive 
	//and is used in applications where the direction of magnetic flux is changing, such as 
	//electric motors and generators. It is also used when efficiency is less important or 
	//when there is insufficient space to correctly orient components to take advantage of 
	//the anisotropic properties of grain-oriented electrical steel. Cold Rolled Non Grain-oriented 
	//steel is often abbreviated to CRNGO.

	//I am using CRGO here : Fe11Si,  3% Si,
	density = 7.65 *g/cm3;
	siliconsteel = new G4Material(name="SiliconSteel", density, nElem=2);
	siliconsteel->AddElement(elFe, nAtoms=11);
	siliconsteel->AddElement(elSi, nAtoms=1);

	//Tantalum
	a = 180.95*g/mole;
	density = 16.69*g/cm3;
	tantalum = new G4Material(name="Tantalum", z=73., a, density);

	//Tungsten
	a = 183.84*g/mole;
	density = 19.25*g/cm3;
	tungsten = new G4Material(name="Tungsten", z=74., a, density);

	//lead
	a = 207.2*g/mole;
	density = 11.34*g/cm3;
	//G4cout << "A gram is " << g << G4endl;
	lead = new G4Material(name="Lead", z=82., a, density);

	//lead208
	a = 207.9766521*g/mole;
	density = 11.38*g/cm3;
	lead208 = new G4Material(name="Lead208", z=82., a, density);
	// G4cout<<lead208->GetDensity()<<" "<<density<<" "<<g<<" "<<cm<<" "<<cm3<<G4endl;
	// G4cout<<joule<<" "<<second<<" "<<meter<<G4endl;
	// G4cout<<lead208->GetTotNbOfElectPerVolume()<<G4endl;
	// std::cin.ignore();



/*
	//solid NH3 	
	//double mSolidNH3D = 0.817*g/cm3;
	density=mSolidNH3D;
	solidNH3 = new G4Material(name="SolidNH3", density, nElem=2);
	solidNH3->AddElement(elH, nAtoms=3);
	solidNH3->AddElement(elN, nAtoms=1);

	//SolidNH3(55%)+LiquidHe(45%) in volumn, if detector.ini provides volumn ratio
	//double mNH3VolumnRatio=0.55;
	//density = pLiquidHeD*(1.0-mNH3VolumnRatio)+pSolidNH3D*mNH3VolumnRatio;
	//double mNH3WeightRatio=pSolidNH3D*mNH3VolumnRatio/density;
	//When detector.ini provides the weight ratio
	//mNH3WeightRatio=0.8732;
	density = 1.0/((1.0-mNH3WeightRatio)/mLiquidHeD+mNH3WeightRatio/mSolidNH3D);
	NH3He = new G4Material(name="SolidNH3+LiquidHe", density, nComponent=2);
	NH3He->AddMaterial(solidNH3, weightRatio=mNH3WeightRatio);
	NH3He->AddMaterial(liquidHe, weightRatio=1.0-mNH3WeightRatio);
*/
	//Copper
	a = 63.546*g/mole;
	density = 8.96*g/cm3;
	copper = new G4Material(name="Copper", z=29., a, density);
/*
	//check here for tha name of the material
	//http://geant4.cern.ch/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/apas09.html
	beryllium = pMatMan->FindOrBuildMaterial("G4_Be");       // pure beryllium

	//Beryllium oxide
	BeO = new G4Material(name="BeO", density=3.02*g/cm3, nComponent=2); 
	BeO->AddElement(elBe, nAtoms=1);
	BeO->AddElement(elO, nAtoms=1);


	//G4_POLYETHYLENE  (C_2H_4)_N-Polyethylene  0.94 g/cm3,  
	//             1     0.143711
	//             6     0.856289
	boron   = pMatMan->FindOrBuildMaterial("G4_B");       // pure boron
	plastic = pMatMan->FindOrBuildMaterial("G4_POLYETHYLENE");  // pure plastic = SWX213
	//http://www.shieldwerx.com/poly-neutron.html
	// 5% borated polyethylene = SWX203
	boratedpoly05 = new G4Material(name="BoratedPoly05", density=1.06*g/cm3, nComponent=2); 
	boratedpoly05->AddMaterial(boron,0.05);
	boratedpoly05->AddMaterial(plastic,0.95);
	// 5% borated polyethylene = SWX210
	boratedpoly30 = new G4Material(name="BoratedPoly30", density=1.19*g/cm3, nComponent=2); 
	boratedpoly30->AddMaterial(boron,0.30);
	boratedpoly30->AddMaterial(plastic,0.70);

	//PCTFE (Polychlorotrifluoroethylene)  [CF2-CFCl]n
	PCTFE = new G4Material(name="PCTFE", density=2.12*g/cm3, nElem=3);
	PCTFE->AddElement(elC, nAtoms=2); 
	PCTFE->AddElement(elF, nAtoms=2); 
	PCTFE->AddElement(elCl, nAtoms=1); 

	//mu metal density=8.25 g/cm3
	///Typical Chemistry (Wt.%)
	//	C	0.015
	//	Mn	0.50
	//	P	0.005 max
	//	S	0.001 max
	//	Si	0.30
	//	Cr	0.02 max
	//	Ni	80.20
	//	Mo	4.85
	//	Al	0.01 max
	//	Co	0.02 max
	//	Fe	Balance
	const double percent=0.01;
	MuMetal = new G4Material(name="MuMetal", density=8.25*g/cm3, nElem=11);
	MuMetal->AddElement(elC, weightRatio=0.015*percent); 
	MuMetal->AddElement(elMn,weightRatio=0.500*percent); 
	MuMetal->AddElement(elP, weightRatio=0.005*percent); 
	MuMetal->AddElement(elS, weightRatio=0.001*percent); 
	MuMetal->AddElement(elSi,weightRatio=0.300*percent); 
	MuMetal->AddElement(elCr,weightRatio=0.020*percent); 
	MuMetal->AddElement(elNi,weightRatio=80.20*percent); 
	MuMetal->AddElement(elMo,weightRatio=4.850*percent);
	MuMetal->AddElement(elAl,weightRatio=0.010*percent); 
	MuMetal->AddElement(elCo,weightRatio=0.020*percent); 
	MuMetal->AddElement(elFe,weightRatio=14.079*percent);  


	//absorber, some thing made of iron but have huge desity so
	//it will absorb all particle 
	a = 55.845*g/mole;
	density = 7.874*g/cm3*1.0E10;
	absorber = new G4Material(name="Absorber", z=26., a, density);
*/

	BuildMaterialMap();

#ifdef G4DEBUG_MATERIAL
	if(G4DEBUG_MATERIAL>=0)
	{
		G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
		G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	}
#endif

}

G4Material* g4hrsMaterial::FindMaterial(G4String name)
{
	map<G4String, G4Material*>::iterator it = mMaterialMap.find(name);
	if(it == mMaterialMap.end())
	{
		G4cout<<"Material \""<<name<<"\" NOT found..."<<G4endl;
		return NULL;
	}
	return mMaterialMap[name];
}

void g4hrsMaterial::BuildMaterialMap()
{
	G4MaterialTable* matTable = (G4MaterialTable*)G4Material::GetMaterialTable();
	for(size_t i=0;i<matTable->size();i++) 
	{
		G4Material* pMat = (*(matTable))[i];
		mMaterialMap[pMat->GetName()] = pMat;
	}
}

/////////////////////////////////////////////////////////////////////
void g4hrsMaterial::DestroyMaterials()
{
	// Destroy all allocated elements and materials
	size_t i;
	G4ElementTable* elemTable = (G4ElementTable*)G4Element::GetElementTable();
	for(i=0;i<elemTable->size();i++) { if((*(elemTable))[i]) delete (*(elemTable))[i]; }
	elemTable->clear();
	G4MaterialTable* matTable = (G4MaterialTable*)G4Material::GetMaterialTable();
	for(i=0;i<matTable->size();i++) { if((*(matTable))[i]) delete (*(matTable))[i]; }
	matTable->clear();
}

