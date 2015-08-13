#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include "SampleHandler/DiskListEOS.h"
#include <TSystem.h>

// xAH Event Selection
#include "xAODAnaHelpers/BasicEventSelection.h"

// xAH Muons
#include "xAODAnaHelpers/MuonCalibrator.h"
#include "xAODAnaHelpers/MuonSelector.h"
#include "xAODAnaHelpers/MuonEfficiencyCorrector.h"

// xAH Electrons
#include "xAODAnaHelpers/ElectronCalibrator.h"
#include "xAODAnaHelpers/ElectronSelector.h"
#include "xAODAnaHelpers/ElectronEfficiencyCorrector.h"

// xAH Jets
#include "xAODAnaHelpers/JetCalibrator.h"
#include "xAODAnaHelpers/JetSelector.h"
#include "xAODAnaHelpers/BJetEfficiencyCorrector.h"

// xAH Overlap Removal
#include "xAODAnaHelpers/OverlapRemover.h"

// Our Balancing Algorithm
#include "ZJetBalance/BalanceAlgorithm.h"

#include <string>
#include <sys/stat.h>
#include <fstream>
#include <unistd.h>

#include "TEnv.h"
#include "TString.h"
#include "TSystem.h"

using namespace std;

int main( int argc, char* argv[] ) {

  std::string samplePath = ".";
  std::string inputTag = "";
  std::string outputTag = "";
  std::string submitDir = "submitDir";
  std::string leptonChoice = "Muon";
  int limitEvents = -1;

  std::string systName = "None";
  float systVal = 0;

  // True -> use Muons; False -> Use Electrons
  bool useMuons = true;

  /////////// Retrieve arguments //////////////////////////
  std::vector< std::string> options;
  for(int ii=1; ii < argc; ++ii){
    options.push_back( argv[ii] );
  }

  if (argc > 1 && options.at(0).compare("-h") == 0) {
    std::cout << std::endl
         << " runZJetBalance: ZJetBalance job submission" << std::endl
         << std::endl
         << " Optional arguments:" << std::endl
         << "  -h               Prints this menu" << std::endl
         << "  -inFile          Path to a folder, root file, or text file" << std::endl
         << "  -inputTag        A wildcarded file name to run on" << std::endl
         << "  -outputTag       Version string to be appended to job name" << std::endl
         << "  -submitDir       Name of output directory" << std::endl
         << "  -leptonChoice    Choose channel to run on, either Muon or Electron" << std::endl
         << "  -limitEvents     Limit the number of events processed for debugging" << std::endl
         << "  -syst            Name AND value for systematic" << std::endl
         << std::endl;
    exit(1);
  }

  int iArg = 0;
  while(iArg < argc-1) {
    if (options.at(iArg).compare("-h") == 0) {
       // Ignore if not first argument
       ++iArg;
    } else if (options.at(iArg).compare("-inFile") == 0) {
       if (iArg+1 == argc || iArg+1 == (int)options.size() || options.at(iArg+1)[0] == '-' ) {
         std::cout << " -inFile should be followed by a file or folder" << std::endl;
         return 1;
       } else {
         samplePath = options.at(iArg+1);
         iArg += 2;
       }
    } else if (options.at(iArg).compare("-inputTag") == 0) {
       if (iArg+1 == argc || iArg+1 == (int)options.size() || options.at(iArg+1)[0] == '-' ) {
         std::cout << " -inputTag is a wildcarded file name to run on" << std::endl;
         return 1;
       } else {
         inputTag = options.at(iArg+1);
         iArg += 2;
       }
    } else if (options.at(iArg).compare("-outputTag") == 0) {
       if (iArg+1 == argc || iArg+1 == (int)options.size() || options.at(iArg+1)[0] == '-' ) {
         std::cout << " -outputTag should be followed by a job version string" << std::endl;
         return 1;
       } else {
         outputTag = options.at(iArg+1);
         iArg += 2;
       }
    } else if (options.at(iArg).compare("-submitDir") == 0) {
       if (iArg+1 == argc || iArg+1 == (int)options.size() || options.at(iArg+1)[0] == '-' ) {
         std::cout << " -submitDir should be followed by a folder name" << std::endl;
         return 1;
       } else {
         submitDir = options.at(iArg+1);
         iArg += 2;
       }
    } else if (options.at(iArg).compare("-leptonChoice") == 0) {
       if (iArg+1 == argc || iArg+1 == (int)options.size() || options.at(iArg+1)[0] == '-' ) {
         std::cout << " -leptonChoice should be followed by either Muon or Electron" << std::endl;
         return 1;
       } else {
         leptonChoice = options.at(iArg+1);
         iArg += 2;
       }
    }  else if (options.at(iArg).compare("-limitEvents") == 0) {
       if (iArg+1 == argc || iArg+1 == (int)options.size() || options.at(iArg+1)[0] == '-' ) {
         std::cout << " -limitEvents should be followed by either Muon or Electron" << std::endl;
         return 1;
       } else {
         limitEvents = stoi(options.at(iArg+1));
         iArg += 2;
       }
    } else if (options.at(iArg).compare("-syst") == 0) {
       if (iArg+1 == argc || iArg+2 == argc || iArg+1 == (int)options.size() || options.at(iArg+1)[0] == '-' ) {
         std::cout << " -inFile should be followed by a systematic string and an integer" << std::endl;
         return 1;
       } else {
         std::stringstream ss;
         systName = options.at(iArg+1);
         ss << options.at(iArg+2);
         ss >> systVal;
         ss.str("");
         if(systVal == 0 && options.at(iArg+2)[0] != '0'){
           std::cout << " -inFile should be followed by a systematic string and an INTEGER" << std::endl;
           return 1;
         }
         iArg += 3;
       }
    }else{
      std::cout << "Couldn't understand argument " << options.at(iArg) << std::endl;
      return 1;
    }
  }//while arguments

  // parse lepton choice
  if ( leptonChoice.find( "Electron" ) != string::npos || leptonChoice.find( "electron" ) != string::npos ){
    useMuons = false;
    cout << "Using Electrons for Balance Algorithm" << endl;
  } else cout << "Using Muons for Balance Algorithm" << endl;


  //if grid job
  bool f_grid = false;
  bool f_lxbatch = false;
  bool f_samplelist = false;

  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // Construct the samples to run on:
  SH::SampleHandler sh;
  std::string containerName;
  std::string userName = getlogin();
  std::vector< std::string > outputContainerNames; //for grid only

  //Check if input is a directory or a file
  struct stat buf;
  stat(samplePath.c_str(), &buf);
  if( samplePath.substr(0, 4).find("eos") != std::string::npos){
    SH::DiskListEOS list(samplePath.c_str());
    if (inputTag.size() > 0){
      SH::scanDir (sh, list, inputTag); //Run on all files within dir containing inputTag
    }else{
      SH::scanDir (sh, list); //Run on all files within dir
    }
    std::cout << "Running on EOS directory " << samplePath << std::endl;
  }else if( S_ISDIR(buf.st_mode) ){ //if it is a local directory
    SH::DiskListLocal list (samplePath);
    if (inputTag.size() > 0){
      SH::scanDir (sh, list, inputTag); //Run on all files within dir containing inputTag
    }else{
      SH::scanDir (sh, list); //Run on all files within dir
    }
    std::cout << "Running Locally on directory  " << samplePath << std::endl;

  } else {  //if it is a file
    if( samplePath.substr( samplePath.size()-4 ).find(".txt") != std::string::npos){ //It is a text file of samples
      if( samplePath.find("grid") != string::npos || samplePath.find("Grid") != string::npos || samplePath.find("GRID") != string::npos ) //It is samples for the grid
        f_grid = true;
      if( samplePath.find("samplelist") != string::npos ) // It is a list of AOD sample files, don't treat them as individual containers
        f_samplelist = true;

      std::ifstream inFile( samplePath );
      if ( f_samplelist )
        SH::readFileList( sh , "sample", samplePath );
      while( !f_samplelist && std::getline(inFile, containerName) ){
        if (containerName.size() > 1 && containerName.find("#") != 0 ){
          std::cout << "Adding container " << containerName << std::endl;
          if(f_grid){
            SH::scanDQ2( sh, containerName);
            //Add output container name to file of containers
            //follows grid format: "user."+userName+".%in:name[1]%.%in:name[2]%.%in:name[3]%"+outputTag
            int startPosition = 0;
            int namePosition = 0;
            // add for JetETMiss "private" samples
            //startPosition = containerName.find_first_of(".", namePosition)+1;
            //namePosition = containerName.find_first_of(".", namePosition)+1;
            namePosition = containerName.find_first_of(".", namePosition)+1;
            namePosition = containerName.find_first_of(".", namePosition)+1;
            namePosition = containerName.find_first_of(".", namePosition)+1;
            std::string outstr = "user."+userName+"."+containerName.substr(startPosition, namePosition)+outputTag+"/";
            outputContainerNames.push_back( outstr );
          } else{
            //Get full path of file
            char fullPath[300];
            realpath( containerName.c_str(), fullPath );
            string thisPath = fullPath;
            //split into fileName and directory two levels above file
            string fileName = thisPath.substr(containerName.find_last_of("/")+1);
            thisPath = thisPath.substr(0, thisPath.find_last_of("/"));
            thisPath = thisPath.substr(0, thisPath.find_last_of("/"));
            std::cout << "path and filename are " << thisPath << " and " << fileName << std::endl;

            SH::DiskListLocal list (thisPath);
            //SH::SampleHandler sh_tmp;
            //SH::scanDir (sh_tmp, list);
            //sh.add( sh_tmp.findByName, ("*"+fileName).c_str() );
            SH::scanDir (sh, list, fileName); // specifying one particular file for testing
          }
        }
      }
    }else{ //It is a single root file to run on
      //Get full path of file
      char fullPath[300];
      realpath( samplePath.c_str(), fullPath );
      string thisPath = fullPath;
      //split into fileName and directory two levels above file
      string fileName = thisPath.substr(thisPath.find_last_of("/")+1);
      thisPath = thisPath.substr(0, thisPath.find_last_of("/"));
      thisPath = thisPath.substr(0, thisPath.find_last_of("/"));

      std::cout << "path and file " << thisPath << " and " << fileName << std::endl;
      SH::DiskListLocal list (thisPath);
      SH::scanDir (sh, list, fileName); // specifying one particular file for testing

    }
  }//it's a file

  ///////// Set output container name //////////////
  std::string outputName;
  if( outputTag.size() > 0)
    outputTag = "."+outputTag+"/";
  else
    outputTag = "/";

  if(f_grid)
    outputName = "user."+userName+".%in:name[1]%.%in:name[2]%.%in:name[3]%"+outputTag;
  // for JetETMiss "private" samples
  //  outputName = "user."+userName+".%in:name[3]%.%in:name[4]%.%in:name[5]%"+outputTag;
  else
    outputName = "%in:name%"+outputTag;

  // Set the name of the input TTree. It's always "CollectionTree" for xAOD files.
  sh.setMetaString( "nc_tree", "CollectionTree" );
  sh.setMetaString("nc_grid_filter", "*");  //Data files on grid to not end in .root
  sh.print();

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler( sh );

  // For debugging purposes, limit the amount of events that we loop over.
  if(limitEvents > -1) job.options()->setDouble (EL::Job::optMaxEvents, limitEvents);

  // To automatically delete submitDir
  job.options()->setDouble(EL::Job::optRemoveSubmitDir, 1);

  // Use TTreeCache to read-precache data files to speed up analysis
  job.options()->setDouble (EL::Job::optCacheSize, 10*1024*1024);
  job.options()->setDouble (EL::Job::optCacheLearnEntries, 20);

  // For Trigger
  job.options()->setString( EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_branch );
  
  // basic event selection : GRL, event cleaning, NPV
  BasicEventSelection* baseEventSel = new BasicEventSelection();
  if ( useMuons )
    baseEventSel->setName("baseEventSel")->setConfig( "$ROOTCOREBIN/data/ZJetBalance/baseEventForMuons.config" );
  else
    baseEventSel->setName("baseEventSel")->setConfig( "$ROOTCOREBIN/data/ZJetBalance/baseEventForElectrons.config" );

  // Declare all lepton operations first, initialization comes later. Slightly inconsistent with other algos, but saves system resources.
  /// MUONS ///
  MuonCalibrator* muonCalib = new MuonCalibrator();
  muonCalib->setName( "muonCalib" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/muonCalib.config")->setSyst( systName, systVal );

  MuonSelector* muonPreSelect = new MuonSelector();
  muonPreSelect->setName( "muonPreSelect" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/muonPreSelect.config");

  MuonSelector* muonSelect = new MuonSelector();
  muonSelect->setName( "muonSelect" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/muonSelect.config");

  MuonSelector* muonSelectForMuonInJetCorrection = new MuonSelector();
  muonSelectForMuonInJetCorrection->setName( "muonSelect" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/muonSelectForMuonInJetCorrection.config");

  MuonEfficiencyCorrector* muonCorrect = new MuonEfficiencyCorrector();
  muonCorrect->setName( "muonCorrect" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/muonCorrect.config");
  
  /// ELECTRONS ///
  ElectronCalibrator* electronCalib = new ElectronCalibrator();
  electronCalib->setName( "electronCalib" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/electronCalib.config")->setSyst( systName, systVal );

  ElectronSelector* electronPreSelect = new ElectronSelector();
  electronPreSelect->setName( "electronPreSelect" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/electronPreSelect.config");

  ElectronSelector* electronSelect = new ElectronSelector();
  electronSelect->setName( "electronSelect" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/electronSelect.config");

  ElectronEfficiencyCorrector* electronCorrect = new ElectronEfficiencyCorrector();
  electronCorrect->setName( "electronCorrect" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/electronCorrect.config");
    

  /// JETS ///
  // jet calibrator
  JetCalibrator* jetCalib = new JetCalibrator();
  jetCalib->setName( "jetCalib" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/jetCalib_AntiKt4EMTopo.config")->setSyst( systName, systVal );

  // jet selector
  JetSelector* jetSelect = new JetSelector();
  jetSelect->setName( "jetSelect" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/jetSelect.config" );

  // bjet efficiecny corrector
  BJetEfficiencyCorrector* bjetCorrectFix85 = new BJetEfficiencyCorrector();
  bjetCorrectFix85->setName( "bjetCorrectFix85" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/bjetCorrectFix85.config" );
  BJetEfficiencyCorrector* bjetCorrectFix77 = new BJetEfficiencyCorrector();
  bjetCorrectFix77->setName( "bjetCorrectFix77" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/bjetCorrectFix77.config" );
  BJetEfficiencyCorrector* bjetCorrectFix70 = new BJetEfficiencyCorrector();
  bjetCorrectFix70->setName( "bjetCorrectFix70" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/bjetCorrectFix70.config" );
  BJetEfficiencyCorrector* bjetCorrectFix60 = new BJetEfficiencyCorrector();
  bjetCorrectFix60->setName( "bjetCorrectFix60" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/bjetCorrectFix60.config" );
  
  BJetEfficiencyCorrector* bjetCorrectFlt70 = new BJetEfficiencyCorrector();
  bjetCorrectFlt70->setName( "bjetCorrectFlt70" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/bjetCorrectFlt70.config" );
  /// OVERLAP REMOVAL ///
  OverlapRemover* overlapRemover = new OverlapRemover();
  overlapRemover->setName( "overlapRemover" )->setConfig( "$ROOTCOREBIN/data/ZJetBalance/overlapRemoval.config" );

  /// BALANCING ALGORITHM ///
  BalanceAlgorithm* balAlg; 

//  muonCalib->m_debug    = true;
//  muonSelect->m_debug   = true;
//  muonCorrect->m_debug  = true;
//  jetCalib->m_debug     = true;
//  jetSelect->m_debug    = true;


  // ADD ALGOS TO JOB
  job.algsAdd( baseEventSel );

  job.algsAdd( muonCalib    );
  job.algsAdd( muonPreSelect);
  job.algsAdd( electronCalib     );
  job.algsAdd( electronPreSelect );
  job.algsAdd( jetCalib     );
  job.algsAdd( overlapRemover );

  job.algsAdd( jetSelect    );
  job.algsAdd( bjetCorrectFix60 );
  job.algsAdd( bjetCorrectFix70 );
  job.algsAdd( bjetCorrectFix77 );
  job.algsAdd( bjetCorrectFix85 );
  job.algsAdd( bjetCorrectFlt70 );
  job.algsAdd( muonSelectForMuonInJetCorrection   );

  if( useMuons ){
    job.algsAdd( muonSelect      );
    job.algsAdd( muonCorrect     ); 
    balAlg = new BalanceAlgorithm();
    balAlg->setName("ZJetBalanceAlgo")->setConfig( "$ROOTCOREBIN/data/ZJetBalance/zjetAlgo.config" );
  }else{
    job.algsAdd( electronSelect  );
    job.algsAdd( electronCorrect );
    balAlg = new BalanceAlgorithm();
    balAlg->setName("ZJetBalanceAlgo")->setConfig( "$ROOTCOREBIN/data/ZJetBalance/zeejetAlgo.config" );
  }

  job.algsAdd( balAlg );


  if(f_grid){
    EL::PrunDriver driver;
    driver.options()->setString("nc_outputSampleName", outputName);
    driver.options()->setString (EL::Job::optSubmitFlags, "--forceStaged");
    
    driver.options()->setDouble(EL::Job::optGridNFilesPerJob, 2);
    //driver.options()->setString(EL::Job::optRootVer, "5.34.25");
    //driver.options()->setString(EL::Job::optCmtConfig, "x86_64-slc6-gcc48-opt");
    //driver.options()->setDouble("nc_nGBPerJob", 1);
    //driver.options()->setString("nc_excludeSite", ???);
    //driver.options()->setString(EL::Job::optGridMergeOutput, "false");
    //driver.options()->setDouble(EL::Job::optGridMemory,10240); // 10 GB

    //driver.submit(job, submitDir); // with monitoring
    driver.submitOnly(job, submitDir); //without monitoring
  }else if( f_lxbatch){
    std::cout << "Currently not implemented! " << std::endl;
  }else{
    // Run the job using the local/direct driver:
    EL::DirectDriver driver;
    driver.options()->setString("nc_outputSampleName", outputName);
    driver.submit( job, submitDir );
  }

  ///// For grid, save list of ouput containers to the submission directory /////
  std::ofstream fileList((submitDir+"/outputContainers.txt"), std::ios_base::out);
  for( unsigned int iCont=0; iCont < outputContainerNames.size(); ++iCont){
    fileList << outputContainerNames.at(iCont)+"\n";
  }
  fileList.close();


  return 0;
}

