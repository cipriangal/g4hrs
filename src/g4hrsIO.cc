#include "g4hrsIO.hh"

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "G4ParticleDefinition.hh"

#include "g4hrsGenericDetectorHit.hh"
#include "g4hrsGenericDetectorSum.hh"
#include "g4hrsEvent.hh"
#include "g4hrsRun.hh"
#include "g4hrsRunData.hh"
#include "g4hrsBeamTarget.hh"
#include "g4hrsSteppingAction.hh"

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMNode.hpp>


g4hrsIO::g4hrsIO(){

	numTF = 12;
	numTFvar = 4;
	
	sprintf(TFName[0],"sen");
	sprintf(TFName[1],"sm");
	sprintf(TFName[2],"sex");
	sprintf(TFName[3],"col");
	sprintf(TFName[4],"q1ex");
	sprintf(TFName[5],"q2ex");
	sprintf(TFName[6],"den");
	sprintf(TFName[7],"dex");
	sprintf(TFName[8],"q3en");
	sprintf(TFName[9],"q3m");
	sprintf(TFName[10],"q3ex");
	sprintf(TFName[11],"fp");

	sprintf(TFVarName[0],"x");
	sprintf(TFVarName[1],"th");
	sprintf(TFVarName[2],"y");
	sprintf(TFVarName[3],"ph");


	numVB = 14;
	numVar = 6;

	sprintf(VBName[0],"sen");
	sprintf(VBName[1],"sm");
	sprintf(VBName[2],"sex");
	sprintf(VBName[3],"col");
	sprintf(VBName[4],"q1en");
	sprintf(VBName[5],"q1ex");
	sprintf(VBName[6],"q2en");
	sprintf(VBName[7],"q2ex");
	sprintf(VBName[8],"den");
	sprintf(VBName[9],"dex");
	sprintf(VBName[10],"q3en");
	sprintf(VBName[11],"q3ex");
	sprintf(VBName[12],"vdc");
	sprintf(VBName[13],"fp");

	sprintf(VarName[0],"x");
	sprintf(VarName[1],"y");
	sprintf(VarName[2],"z");
	sprintf(VarName[3],"th");
	sprintf(VarName[4],"ph");
	sprintf(VarName[5],"p");

	numZCrit = 24;
	numZCritVar = 6;
	
	sprintf(ZCritName[0], "zpinch1");
	sprintf(ZCritName[1], "zpinch2");
	sprintf(ZCritName[2], "zpinch3");
	sprintf(ZCritName[3], "ztarg");
	sprintf(ZCritName[4], "zfield0");
	sprintf(ZCritName[5], "zfield1");
	sprintf(ZCritName[6], "zfield2");
	sprintf(ZCritName[7], "zmidtosep");
	sprintf(ZCritName[8], "zsep1");
	sprintf(ZCritName[9], "zsep2");
	sprintf(ZCritName[10],"zsep3");
	sprintf(ZCritName[11],"zsep4");
        sprintf(ZCritName[12],"zup1");
        sprintf(ZCritName[13],"zup2");
        sprintf(ZCritName[14],"zdown1");
        sprintf(ZCritName[15],"zdown2");
        sprintf(ZCritName[16],"zdown3");
        sprintf(ZCritName[17],"zdown4");
        sprintf(ZCritName[18],"zdown5");
        sprintf(ZCritName[19],"zdown6");
        sprintf(ZCritName[20],"zdown7");
        sprintf(ZCritName[21],"zdown8");
        sprintf(ZCritName[22],"zdown9");
        sprintf(ZCritName[23],"zsieve");


	sprintf(ZCritVarName[0],"x");
	sprintf(ZCritVarName[1],"y");
	sprintf(ZCritVarName[2],"z");
	sprintf(ZCritVarName[3],"th");
	sprintf(ZCritVarName[4],"ph");
        sprintf(ZCritVarName[5],"p");

    fTree = NULL;
    InitializeTree();
    // Default filename
    strcpy(fFilename, "g4hrsout.root");

    fFile = NULL;
}

g4hrsIO::~g4hrsIO(){
    if( fTree ){ delete fTree; }
    fTree = NULL;
    if( fFile ){ delete fFile; }
    fFile = NULL;
}

void g4hrsIO::SetFilename(G4String fn){
    G4cout << "Setting output file to " << fn << G4endl;
    strcpy(fFilename, fn.data());
}

void g4hrsIO::InitializeTree(){
    if( fFile ){
	fFile->Close();
	delete fFile;
    }

    fFile = new TFile(fFilename, "RECREATE");

    if( fTree ){ delete fTree; }

    fTree = new TTree("T", "Geant4 Moller Simulation");

    fTree->SetMaxTreeSize(1900000000); // 1.9GB

    // Event information
    fTree->Branch("rate",     &fEvRate,   "rate/D");
    fTree->Branch("ev.A",     &fEvAsym,   "ev.A/D");
    fTree->Branch("ev.S",     &fEvSens,   "ev.S/D");
    fTree->Branch("ev.Am",    &fEvmAsym,  "ev.Am/D");
    fTree->Branch("ev.xs",    &fEvEffXS,  "ev.xs/D");
    fTree->Branch("ev.V",     &fEvWeight, "ev.V/D");
    fTree->Branch("ev.Q2",    &fEvQ2,     "ev.Q2/D");
    fTree->Branch("ev.W2",    &fEvW2,     "ev.W2/D");
    fTree->Branch("ev.thcom", &fEvThCoM,  "ev.thcom/D");
    fTree->Branch("ev.beamp",  &fEvBeamP,   "ev.beamp/D");

    fTree->Branch("ev.nuclA", &fEvNuclA, "ev.nuclA/I");

    fTree->Branch("ev.npart", &fNEvPart   ,     "ev.npart/I");
    fTree->Branch("ev.pid",   &fEvPID,      "ev.pid[ev.npart]/I");
    fTree->Branch("ev.vx",    &fEvPart_X,   "ev.vx[ev.npart]/D");
    fTree->Branch("ev.vy",    &fEvPart_Y,   "ev.vy[ev.npart]/D");
    fTree->Branch("ev.vz",    &fEvPart_Z,   "ev.vz[ev.npart]/D");
    fTree->Branch("ev.p",     &fEvPart_P,   "ev.p[ev.npart]/D");
    fTree->Branch("ev.px",    &fEvPart_Px,  "ev.px[ev.npart]/D");
    fTree->Branch("ev.py",    &fEvPart_Py,  "ev.py[ev.npart]/D");
    fTree->Branch("ev.pz",    &fEvPart_Pz,  "ev.pz[ev.npart]/D");
    fTree->Branch("ev.th",    &fEvPart_Th,     "ev.th[ev.npart]/D");
    fTree->Branch("ev.ph",    &fEvPart_Ph,     "ev.ph[ev.npart]/D");
    fTree->Branch("ev.tpx",    &fEvPart_tPx,  "ev.tpx[ev.npart]/D");
    fTree->Branch("ev.tpy",    &fEvPart_tPy,  "ev.tpy[ev.npart]/D");
    fTree->Branch("ev.tpz",    &fEvPart_tPz,  "ev.tpz[ev.npart]/D");

    fTree->Branch("bm.x",    &fBmX,  "bm.x/D");
    fTree->Branch("bm.y",    &fBmY,  "bm.y/D");
    fTree->Branch("bm.z",    &fBmZ,  "bm.z/D");
    fTree->Branch("bm.dx",    &fBmdX,  "bm.dx/D");
    fTree->Branch("bm.dy",    &fBmdY,  "bm.dy/D");
    fTree->Branch("bm.dz",    &fBmdZ,  "bm.dz/D");
    fTree->Branch("bm.th",    &fBmth,  "bm.th/D");
    fTree->Branch("bm.ph",    &fBmph,  "bm.ph/D");

    // GenericDetectorHit
    fTree->Branch("hit.n",    &fNGenDetHit,     "hit.n/I");
    fTree->Branch("hit.det",  &fGenDetHit_det,  "hit.det[hit.n]/I");
    fTree->Branch("hit.vid",  &fGenDetHit_id,   "hit.vid[hit.n]/I");

    fTree->Branch("hit.pid",  &fGenDetHit_pid,   "hit.pid[hit.n]/I");
    fTree->Branch("hit.trid", &fGenDetHit_trid,  "hit.trid[hit.n]/I");
    fTree->Branch("hit.mtrid",&fGenDetHit_mtrid, "hit.mtrid[hit.n]/I");
    fTree->Branch("hit.gen",  &fGenDetHit_gen,   "hit.gen[hit.n]/I");

    fTree->Branch("hit.x",    &fGenDetHit_X,   "hit.x[hit.n]/D");
    fTree->Branch("hit.y",    &fGenDetHit_Y,   "hit.y[hit.n]/D");
    fTree->Branch("hit.z",    &fGenDetHit_Z,   "hit.z[hit.n]/D");
    fTree->Branch("hit.r",    &fGenDetHit_R,   "hit.r[hit.n]/D");
    fTree->Branch("hit.ph",   &fGenDetHit_Ph,  "hit.ph[hit.n]/D");

    fTree->Branch("hit.px",   &fGenDetHit_Px,   "hit.px[hit.n]/D");
    fTree->Branch("hit.py",   &fGenDetHit_Py,   "hit.py[hit.n]/D");
    fTree->Branch("hit.pz",   &fGenDetHit_Pz,   "hit.pz[hit.n]/D");

    fTree->Branch("hit.vx",   &fGenDetHit_Vx,   "hit.vx[hit.n]/D");
    fTree->Branch("hit.vy",   &fGenDetHit_Vy,   "hit.vy[hit.n]/D");
    fTree->Branch("hit.vz",   &fGenDetHit_Vz,   "hit.vz[hit.n]/D");

    fTree->Branch("hit.p",    &fGenDetHit_P,   "hit.p[hit.n]/D");
    fTree->Branch("hit.e",    &fGenDetHit_E,   "hit.e[hit.n]/D");
    fTree->Branch("hit.m",    &fGenDetHit_M,   "hit.m[hit.n]/D");

    fTree->Branch("hit.colCut",    &fCollCut,     "hit.colCut/I");

    // GenericDetectorSum
    fTree->Branch("sum.n",    &fNGenDetSum,     "sum.n/I");
    fTree->Branch("sum.det",  &fGenDetSum_det,  "sum.det[sum.n]/I");
    fTree->Branch("sum.vid",  &fGenDetSum_id,   "sum.vid[sum.n]/I");
    fTree->Branch("sum.edep", &fGenDetSum_edep, "sum.edep[sum.n]/D");

	// Target data (used as input to transport functions)
	
	fTree->Branch("x_tg",	&fX0,	"fX0/D");
	fTree->Branch("y_tg",	&fY0,	"fY0/D");
	fTree->Branch("z_tg",	&fZ0,	"fZ0/D");
	fTree->Branch("th_tg",	&fTh0,	"fTh0/D");
	fTree->Branch("ph_tg",	&fPh0,	"fPh0/D");
	fTree->Branch("p_tg",	&fP0,	"fP0/D");
	fTree->Branch("x_tg_tr",	&fX0_tr,	"fX0_tr/D");
	fTree->Branch("y_tg_tr",	&fY0_tr,	"fY0_tr/D");
	fTree->Branch("z_tg_tr",	&fZ0_tr,	"fZ0_tr/D");
	fTree->Branch("th_tg_tr",	&fTh0_tr,	"fTh0_tr/D");
	fTree->Branch("ph_tg_tr",	&fPh0_tr,	"fPh0_tr/D");
	fTree->Branch("p_tg_tr",	&fP0_tr,	"fP0_tr/D");

	//Transport function data
	
	for(int i = 0; i < numTF; i++) {
		for(int j = 0; j < numTFvar; j++) {
			char thisVar[15];
			sprintf(thisVar,"%s_%s_tf",TFVarName[j],TFName[i]);
			fTree->Branch(thisVar, &TFdata[j][i], Form("%s/D",thisVar));
		}
	}


	//Virtual boundary data

	for(int i = 0; i < numVB; i++) {
		for(int j = 0; j < numVar; j++) {
			char thisVar[15];
			sprintf(thisVar,"%s_%s",VarName[j],VBName[i]);
			fTree->Branch(thisVar, &VBdata[i][j], Form("%s/D",thisVar));
			sprintf(thisVar,"%s_%s_tr",VarName[j],VBName[i]);
			fTree->Branch(thisVar, &VBdata[i][j+numVar], Form("%s/D",thisVar));
		}
	}

	//Critical z-position data
	
	for(int i = 0; i < numZCrit; i++) {
		for(int j = 0; j < numZCritVar; j++) {
			char thisVar[15];
			sprintf(thisVar,"%s_%s",ZCritVarName[j],ZCritName[i]);
			fTree->Branch(thisVar, &ZCritData[i][j], Form("%s/D",thisVar));
		}
	}


    return;
}

void g4hrsIO::FillTree(){
    if( !fTree ){ 
	fprintf(stderr, "Error %s: %s line %d - Trying to fill non-existant tree\n", __PRETTY_FUNCTION__, __FILE__, __LINE__ );
	return; 
    }

    fTree->Fill();
    fTree->GetCurrentFile();
}

void g4hrsIO::Flush(){
    //  Set arrays to 0
    fNGenDetHit = 0;
    fNGenDetSum = 0;
    fCollCut = 1; // default
}

void g4hrsIO::WriteTree(){
    assert( fFile );
    assert( fTree );

    if( !fFile->IsOpen() ){
	G4cerr << "ERROR: " << __FILE__ << " line " << __LINE__ << ": TFile not open" << G4endl;
	exit(1);
    }

    G4cout << "Writing output to " << fFile->GetName() << "... ";

    fFile->cd();

    fTree->Write("T", TObject::kOverwrite);
    g4hrsRun::GetRun()->GetData()->Write("run_data", TObject::kOverwrite); 

    fTree->ResetBranchAddresses();
    delete fTree;
    fTree = NULL;

    fFile->Close();

    delete fFile;
    fFile = NULL;

    G4cout << "written" << G4endl;

    return;
}

///////////////////////////////////////////////////////////////////////////////
// Interfaces to output section ///////////////////////////////////////////////

// Event Data

void g4hrsIO::SetEventData(g4hrsEvent *ev){
    int n = ev->fPartType.size();
    if( n > __IO_MAXHIT ){
	G4cerr << "WARNING: " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ":  Buffer size exceeded!" << G4endl;
	return;
    }

    fNEvPart = n;

    fEvRate   = ev->fRate; //already in Hz 
    fEvEffXS  = ev->fEffXs/barn; 		
    fEvAsym   = ev->fAsym/__ASYMM_SCALE;
    fEvSens   = ev->fSens;
    fEvWeight = ev->fWeight;
    fEvmAsym  = ev->fmAsym/__ASYMM_SCALE;
    fEvBeamP  = ev->fBeamMomentum.mag()/__E_UNIT;

    fEvQ2     = ev->fQ2/__E_UNIT/__E_UNIT;
    fEvW2     = ev->fW2/__E_UNIT/__E_UNIT;
    fEvThCoM  = ev->fThCoM/deg; // specify this in degrees over anything else

    fEvNuclA = ev->fNuclA;

    int idx;
    for( idx = 0; idx < n; idx++ ){
	fEvPID[idx] = ev->fPartType[idx]->GetPDGEncoding();

	fEvPart_X[idx] = ev->fPartPos[idx].x()/__L_UNIT;
	fEvPart_Y[idx] = ev->fPartPos[idx].y()/__L_UNIT;
	fEvPart_Z[idx] = ev->fPartPos[idx].z()/__L_UNIT;

	fEvPart_Px[idx] = ev->fPartRealMom[idx].x()/__E_UNIT;
	fEvPart_Py[idx] = ev->fPartRealMom[idx].y()/__E_UNIT;
	fEvPart_Pz[idx] = ev->fPartRealMom[idx].z()/__E_UNIT;
	fEvPart_Th[idx] = ev->fPartRealMom[idx].theta();
	fEvPart_Ph[idx] = ev->fPartRealMom[idx].phi();

	fEvPart_P[idx] = ev->fPartRealMom[idx].mag()/__E_UNIT;

	fEvPart_tPx[idx] = ev->fPartMom[idx].x()/__E_UNIT;
	fEvPart_tPy[idx] = ev->fPartMom[idx].y()/__E_UNIT;
	fEvPart_tPz[idx] = ev->fPartMom[idx].z()/__E_UNIT;
    }

    /////////////////////////////////////////////////
    //  Set beam data as well

    g4hrsBeamTarget *bt = g4hrsBeamTarget::GetBeamTarget();

    fBmX = bt->fVer.x()/__L_UNIT;
    fBmY = bt->fVer.y()/__L_UNIT;
    fBmZ = bt->fVer.z()/__L_UNIT;
    
    fBmdX = bt->fDir.x();
    fBmdY = bt->fDir.y();
    fBmdZ = bt->fDir.z();
    fBmth = bt->fDir.theta();
    fBmph = bt->fDir.phi()/deg;

    //    G4cout << "** fDir:: " << bt->fDir.x()/deg << "  " << bt->fDir.y()/deg << "  " << bt->fVer.z()/mm << G4endl;

    return;
}

// GenericDetectorHit

void g4hrsIO::AddGenericDetectorHit(g4hrsGenericDetectorHit *hit){
    int n = fNGenDetHit;
    if( n > __IO_MAXHIT ){
	G4cerr << "WARNING: " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ":  Buffer size exceeded!" << G4endl;
	return;
    }

    fGenDetHit_det[n]  = hit->fDetID;
    fGenDetHit_id[n]   = hit->fCopyID;

    fGenDetHit_trid[n] = hit->fTrID;
    fGenDetHit_mtrid[n]= hit->fmTrID;
    fGenDetHit_pid[n]  = hit->fPID;
    fGenDetHit_gen[n]  = hit->fGen;

    fGenDetHit_X[n]  = hit->f3X.x()/__L_UNIT;
    fGenDetHit_Y[n]  = hit->f3X.y()/__L_UNIT;
    fGenDetHit_Z[n]  = hit->f3X.z()/__L_UNIT;
    fGenDetHit_R[n]  = sqrt(hit->f3X.x()*hit->f3X.x()+hit->f3X.y()*hit->f3X.y())/__L_UNIT;
    fGenDetHit_Ph[n] = hit->f3X.phi()/deg;

    fGenDetHit_Px[n]  = hit->f3P.x()/__E_UNIT;
    fGenDetHit_Py[n]  = hit->f3P.y()/__E_UNIT;
    fGenDetHit_Pz[n]  = hit->f3P.z()/__E_UNIT;

    fGenDetHit_Vx[n]  = hit->f3V.x()/__L_UNIT;
    fGenDetHit_Vy[n]  = hit->f3V.y()/__L_UNIT;
    fGenDetHit_Vz[n]  = hit->f3V.z()/__L_UNIT;

    fGenDetHit_P[n]  = hit->fP/__E_UNIT;
    fGenDetHit_E[n]  = hit->fE/__E_UNIT;
    fGenDetHit_M[n]  = hit->fM/__E_UNIT;

    fNGenDetHit++;

    // for collimator cut
    if( (hit->fDetID==200 && hit->f3X.perp()/__L_UNIT < 0.03) || 
      	(hit->fDetID==201 && hit->f3X.perp()/__L_UNIT < 0.05) )
      fCollCut=0;
}


// GenericDetectorSum

void g4hrsIO::AddGenericDetectorSum(g4hrsGenericDetectorSum *hit){
    int n = fNGenDetSum;
    if( n > __IO_MAXHIT ){
	G4cerr << "WARNING: " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ":  Buffer size exceeded!" << G4endl;
	return;
    }

    fGenDetSum_edep[n] = hit->fEdep/__E_UNIT;
    fGenDetSum_det[n]  = hit->fDetID;
    fGenDetSum_id[n]   = hit->fCopyID;

    fNGenDetSum++;
}

// Virtual boundary data

void g4hrsIO::ClearVirtualBoundaryData() {

	fSteppingAction->fLHRS = 0;
	fSteppingAction->fRHRS = 0;
	
	fSteppingAction->fX0 = -333.;
	fSteppingAction->fY0 = -333.;
	fSteppingAction->fZ0 = -333.;
	fSteppingAction->fTh0 = -333.;
	fSteppingAction->fPh0 = -333.;
	fSteppingAction->fP0 = -333.;
	fSteppingAction->fX0_tr = -333.;
	fSteppingAction->fY0_tr = -333.;
	fSteppingAction->fZ0_tr = -333.;
	fSteppingAction->fTh0_tr = -333.;
	fSteppingAction->fPh0_tr = -333.;
	fSteppingAction->fP0_tr = -333.;

	for(int i = 0; i < numTF; i++) {
		for(int j = 0; j < numTFvar; j++) {
			fSteppingAction->TFdata[j][i] = -333.;
		}
	}	

	for(int i = 0; i<numVB; i++) {
		for(int j = 0; j<2*numVar; j++) {
			fSteppingAction->VBdata[i][j] = -333.;
		}
	}

	for(int i = 0; i<numZCrit; i++) {
		for(int j = 0; j< numZCritVar; j++) {
			fSteppingAction->ZCritData[i][j] = -333.;
		}
	}

}

void g4hrsIO::SetVirtualBoundaryData() {

	fLHRS = fSteppingAction->fLHRS;
	fRHRS = fSteppingAction->fRHRS;

	fX0 = fSteppingAction->fX0;
	fY0 = fSteppingAction->fY0;
	fZ0 = fSteppingAction->fZ0;
	fTh0 = fSteppingAction->fTh0;
	fPh0 = fSteppingAction->fPh0;
	fP0 = fSteppingAction->fP0;
	fX0_tr = fSteppingAction->fX0_tr;
	fY0_tr = fSteppingAction->fY0_tr;
	fZ0_tr = fSteppingAction->fZ0_tr;
	fTh0_tr = fSteppingAction->fTh0_tr;
	fPh0_tr = fSteppingAction->fPh0_tr;
	fP0_tr = fSteppingAction->fP0_tr;

	for(int i = 0; i < numTF; i++) {
		for(int j = 0; j < numTFvar; j++) {
			TFdata[j][i] = fSteppingAction->TFdata[j][i];
		}
	}	

	for(int i = 0; i < numVB; i++) {
		for(int j = 0; j < 2*numVar; j++) {
			VBdata[i][j] = fSteppingAction->VBdata[i][j];
		}
	}	
	
	for(int i = 0; i < numZCrit; i++) {
		for(int j = 0; j < numZCritVar; j++) {
			ZCritData[i][j] = fSteppingAction->ZCritData[i][j];
		}
	}
	
}

/*---------------------------------------------------------------------------------*/

void g4hrsIO::GrabGDMLFiles(G4String fn){
    // Reset list
    fGDMLFileNames.clear();

    g4hrsRunData *rundata = g4hrsRun::GetRun()->GetData();
    rundata->ClearGDMLFiles();

    xercesc::XMLPlatformUtils::Initialize();
    SearchGDMLforFiles(fn);
    xercesc::XMLPlatformUtils::Terminate();


    // Store filename

    unsigned int idx;

    // Copy into buffers
    for( idx = 0; idx < fGDMLFileNames.size(); idx++ ){
	G4cout << "Found GDML file " << fGDMLFileNames[idx] << G4endl;
	rundata->AddGDMLFile(fGDMLFileNames[idx]);
    }

    return;
}

void g4hrsIO::SearchGDMLforFiles(G4String fn){
    /*!  Chase down files to be included by GDML.
     *   Mainly look for file tags and perform recursively */

    struct stat thisfile;

    int ret = stat(fn.data(), &thisfile);

    if( ret != 0 ){
	G4cerr << "ERROR opening file " << fn <<  " in " << __PRETTY_FUNCTION__ << ": " << strerror(errno) << G4endl;
	exit(1);
    }

   xercesc::XercesDOMParser *xmlParser = new xercesc::XercesDOMParser();

   // Make sure file exists - otherwise freak out

   fGDMLFileNames.push_back(fn.data());

   xmlParser->parse( fn.data() );
   xercesc::DOMDocument* xmlDoc = xmlParser->getDocument();

   xercesc::DOMElement* elementRoot = xmlDoc->getDocumentElement();

   TraverseChildren( elementRoot );
   return;
}

void g4hrsIO::TraverseChildren( xercesc::DOMElement *thisel ){

   xercesc::DOMNodeList*      children = thisel->getChildNodes();
   const XMLSize_t nodeCount = children->getLength();

   for( XMLSize_t xx = 0; xx < nodeCount; ++xx ){
       xercesc::DOMNode* currentNode = children->item(xx);
       if( currentNode->getNodeType() ){   // true is not NULL

	   if( currentNode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE ){ // is element 
	       xercesc::DOMElement* currentElement
		   = dynamic_cast< xercesc::DOMElement* >( currentNode );
	       if( xercesc::XMLString::equals(currentElement->getTagName(), xercesc::XMLString::transcode("file"))){
		   SearchGDMLforFiles(G4String(xercesc::XMLString::transcode(currentElement->getAttribute(xercesc::XMLString::transcode("name")))));
	       }

	       if( currentElement->getChildNodes()->getLength() > 0 ){
		   TraverseChildren( currentElement );
	       }
	   }
       }
   }

}









