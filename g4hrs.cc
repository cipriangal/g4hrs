/*!
  g4hrs - PREX/CREX Optics Simluation

  Seamus Riordan, et al.
  riordan@jlab.org

*/

#include "CLHEP/Random/Random.h"

#include "g4hrsRunAction.hh"
#include "g4hrsRun.hh"
#include "g4hrsRunData.hh"
#include "g4hrsPrimaryGeneratorAction.hh"
#include "g4hrsEventAction.hh"
#include "g4hrsSteppingAction.hh"
#include "g4hrsDetectorConstruction.hh"
#include "g4hrsParallelWorld.hh"

#include "g4hrsIO.hh"
#include "g4hrsMessenger.hh"

//  Standard physics list
#include "G4Version.hh"
#include "G4PhysListFactory.hh"
#include "G4OpticalPhysics.hh"
#if G4VERSION_NUMBER < 1000
#include "LHEP.hh"
#endif

#include "G4RunManager.hh"

#include "G4UnitsTable.hh"

#include "G4RunManagerKernel.hh"

#include "G4ParallelWorldPhysics.hh"

//to make gui.mac work
#include <G4UImanager.hh>
#include <G4UIExecutive.hh>
#include <G4UIterminal.hh>

#ifdef G4UI_USE_QT
#include "G4UIQt.hh"
#endif

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4UI_USE_TCSH
#include "G4UItcsh.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include <sys/types.h>
#include <sys/stat.h>

#ifdef __APPLE__
#include <unistd.h>
#endif

#include <time.h>

int main(int argc, char** argv){

    clock_t tStart=clock();
    // Initialize the CLHEP random engine used by
    // "shoot" type functions
    unsigned int seed = time(0) + (int) getpid();

    unsigned int devrandseed = 0;
    //  /dev/urandom doens't block
    FILE *fdrand = fopen("/dev/urandom", "r");
    if( fdrand ){
	fread(&devrandseed, sizeof(int), 1, fdrand);
	seed += devrandseed;
	fclose(fdrand);
    }

    G4Random::createInstance();
    G4Random::setTheSeed(seed);

    g4hrsRun::GetRun()->GetData()->SetSeed(seed);

    g4hrsIO *io = new g4hrsIO();

    //-------------------------------
    // Initialization of Run manager
    //-------------------------------
    //G4cout << "RunManager construction starting...." << G4endl;
    G4RunManager * runManager = new G4RunManager;

    g4hrsMessenger *rmmess = new g4hrsMessenger();
    rmmess->SetIO(io);

    // Detector geometry
    G4VUserDetectorConstruction* detector = new g4hrsDetectorConstruction();
	G4String parallelWorldName = "g4hrsparallel";
    detector->RegisterParallelWorld(new g4hrsParallelWorld(parallelWorldName));
    runManager->SetUserInitialization(detector);
    rmmess->SetDetCon( ((g4hrsDetectorConstruction *) detector) );

    ((g4hrsDetectorConstruction *) detector)->SetIO(io);

//    rmmess->SetEmFieldSetup(((g4hrsDetectorConstruction *) detector)->GetEMFieldSetup());
//    rmmess->SetEMField(((g4hrsDetectorConstruction *) detector)->GetEMFieldFromSetup());

    // Physics we want to use
    G4int verbose = 0;
    G4PhysListFactory factory;
    #if G4VERSION_NUMBER < 1000
    G4VModularPhysicsList* physlist = factory.GetReferencePhysList("LHEP");
    #else
    //G4VModularPhysicsList* physlist = factory.GetReferencePhysList("FTFP_BERT_LIV");
    //G4VModularPhysicsList* physlist = factory.GetReferencePhysList("QGSP_BERT_HP");
    G4VModularPhysicsList* physlist = factory.GetReferencePhysList("FTFP_BERT");
    #endif
    physlist->RegisterPhysics(new G4OpticalPhysics());
    physlist->SetVerboseLevel(verbose);
    physlist->RegisterPhysics(new G4ParallelWorldPhysics(parallelWorldName));
    runManager->SetUserInitialization(physlist);

    //-------------------------------
    // UserAction classes
    //-------------------------------
    G4UserRunAction* run_action = new g4hrsRunAction;
    ((g4hrsRunAction *) run_action)->SetIO(io);
    runManager->SetUserAction(run_action);

    G4VUserPrimaryGeneratorAction* gen_action = new g4hrsPrimaryGeneratorAction;
    ((g4hrsPrimaryGeneratorAction *) gen_action)->SetIO(io);
    rmmess->SetPriGen((g4hrsPrimaryGeneratorAction *)gen_action);
    runManager->SetUserAction(gen_action);

    G4UserEventAction* event_action = new g4hrsEventAction;
    ((g4hrsEventAction *) event_action)->SetIO(io);

    runManager->SetUserAction(event_action);
    g4hrsSteppingAction* stepping_action = new g4hrsSteppingAction((g4hrsEventAction*)event_action);
    //G4UserSteppingAction* stepping_action = new g4hrsSteppingAction((g4hrsEventAction*)event_action);
    runManager->SetUserAction(stepping_action);
    rmmess->SetStepAct((g4hrsSteppingAction *) stepping_action);
	io->SetSteppingAction((g4hrsSteppingAction *) stepping_action);

    // New units

    G4UIsession* session = 0;

    //----------------
    // Visualization:
    //----------------

    if (argc==1)   // Define UI session for interactive mode.
    {

	// G4UIterminal is a (dumb) terminal.
#if defined(G4UI_USE_QT)
	session = new G4UIQt(argc,argv);
#elif defined(G4UI_USE_WIN32)
	session = new G4UIWin32();
#elif defined(G4UI_USE_XM)
	session = new G4UIXm(argc,argv);
#elif defined(G4UI_USE_TCSH)
	session = new G4UIterminal(new G4UItcsh);
#else
	session = new G4UIterminal();
#endif

    }

    g4hrsRunData *rundata = g4hrsRun::GetRun()->GetData();


//	Don't initialize here, initialize in macro		
//    	runManager->Initialize();

#ifdef G4VIS_USE
    // Visualization, if you choose to have it!
    //
    // Simple graded message scheme - give first letter or a digit:
    //  0) quiet,         // Nothing is printed.
    //  1) startup,       // Startup and endup messages are printed...
    //  2) errors,        // ...and errors...
    //  3) warnings,      // ...and warnings...
    //  4) confirmations, // ...and confirming messages...
    //  5) parameters,    // ...and parameters of scenes and views...
    //  6) all            // ...and everything available.

    //this is the initializing the run manager?? Right?
    G4VisManager* visManager = new G4VisExecutive;
    //visManager -> SetVerboseLevel (1);
    visManager ->Initialize();
#endif

    //get the pointer to the User Interface manager
    G4UImanager * UI = G4UImanager::GetUIpointer();

    if (session)   // Define UI session for interactive mode.
    {

        std::cout << "Executing gui "  << std::endl;
	// G4UIterminal is a (dumb) terminal.
	//UI->ApplyCommand("/control/execute myVis.mac");

#if defined(G4UI_USE_XM) || defined(G4UI_USE_WIN32) || defined(G4UI_USE_QT)
	// Customize the G4UIXm,Win32 menubar with a macro file :
	UI->ApplyCommand("/control/execute macros/gui.mac");
#endif

	session->SessionStart();
	delete session;
    }
    else           // Batch mode - not using the GUI
    {
#ifdef G4VIS_USE
	visManager->SetVerboseLevel("quiet");
#endif
	//these line will execute a macro without the GUI
	//in GEANT4 a macro is executed when it is passed to the command, /control/execute
	G4String command = "/control/execute ";
	G4String fileName = argv[1];

	/* Copy contents of macro into buffer to be written out
	 * into ROOT file
	 * */

        std::cout << "Executing " << argv[1] << std::endl;

	rundata->SetMacroFile(argv[1]);


	UI->ApplyCommand(command+fileName);
    }

    //if one used the GUI then delete it
#ifdef G4VIS_USE
    delete visManager;
#endif

    // Initialize Run manager
    // runManager->Initialize();
    stepping_action->Write();
    G4cout<<" Running time[s]: "<< (double) ((clock() - tStart)/CLOCKS_PER_SEC)<<G4endl;
    return 0;
}
