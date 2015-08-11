#include <xAODRootAccess/Init.h>
#include <SampleHandler/SampleHandler.h>
#include <SampleHandler/ScanDir.h>
#include <SampleHandler/DiskListEOS.h>
#include <SampleHandler/DiskListXRD.h>
#include <SampleHandler/DiskListLocal.h>
#include <SampleHandler/ToolsDiscovery.h>
#include <EventLoop/Job.h>
#include <EventLoop/DirectDriver.h>
#include <EventLoop/OutputStream.h>
#include <EventLoopGrid/PrunDriver.h>
#include <ZJetBalance/ZJetBalanceMiniTreeAnaSkeleton.h>
#include <ZJetBalance/ZJetBalanceMiniTree_GenBalanceHistograms.h>

#include <string>
#include <sys/stat.h>
#include <fstream>
#include <unistd.h>

#include <TEnv.h>
#include <TString.h>
#include <TSystem.h>


int main( int argc, char* argv[] ) {

  //
  // Create the EventLoop job:
  //
  EL::Job job;
  
  //
  // Init various job options
  //
  std::string configName = "$ROOTCOREBIN/data/ZJetBalance/ZJetBalanceMiniTreeAnaSkeleton.config";
  std::string treeName   = "outTree";
  std::string submitDir  = "submitDir";
  std::string outputName = "";
  
  std::string samplePath = ".";
  std::string inputTag;
  std::string outputTag = "";
  int algoId = -1;
  
  //
  // Set up various job options
  //
  
  /////////// Retrieve job arguments //////////////////////////
  std::vector< std::string> options;
  for(int ii=1; ii < argc; ++ii){
    options.push_back( std::string(argv[ii]) );
  }
  
  if (argc > 1 && options.at(0).compare("-h") == 0) {
    std::cout << std::endl
	      << " job submission" << std::endl
	      << std::endl
	      << " Optional arguments:" << std::endl
	      << "  -h               Prints this menu" << std::endl
	      << "  -inFile          Path to a folder, root file, or text file" << std::endl
	      << "  -algorithm       Algorithm ID ([1]=Skeleton [2]=GenHists)" << std::endl
	      << "  -outputTag       Version string to be appended to job name" << std::endl
	      << "  -submitDir       Name of output directory" << std::endl
	      << "  -configName      Path to config file" << std::endl
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
       char tmpChar = options.at(iArg+1)[0];
       if (iArg+1 == argc || tmpChar == '-' ) {
         std::cout << " -inFile should be followed by a file or folder" << std::endl;
         return 1;
       } else {
         samplePath = options.at(iArg+1);
         iArg += 2;
       }

    } else if (options.at(iArg).compare("-outputTag") == 0) {
       char tmpChar = options.at(iArg+1)[0];
       if (iArg+1 == argc || tmpChar == '-' ) {
         std::cout << " -outputTag should be followed by a job version string" << std::endl;
         return 1;
       } else {
         outputTag = options.at(iArg+1);
         iArg += 2;
       }

    } else if (options.at(iArg).compare("-submitDir") == 0) {
       char tmpChar = options.at(iArg+1)[0];
       if (iArg+1 == argc || tmpChar == '-' ) {
         std::cout << " -submitDir should be followed by a folder name" << std::endl;
         return 1;
       } else {
         submitDir = options.at(iArg+1);
         iArg += 2;
       }

    } else if (options.at(iArg).compare("-configName") == 0) {
       char tmpChar = options.at(iArg+1)[0];
       if (iArg+1 == argc || tmpChar == '-' ) {
         std::cout << " -configName should be followed by a config file" << std::endl;
         return 1;
       } else {
         configName = options.at(iArg+1);
         iArg += 2;
       }

    } else if (options.at(iArg).compare("-treeName") == 0) {
       char tmpChar = options.at(iArg+1)[0];
       if (iArg+1 == argc || tmpChar == '-' ) {
         std::cout << " -treeName should be followed by a tree name" << std::endl;
         return 1;
       } else {
         treeName = options.at(iArg+1);
         iArg += 2;
       }

    } else if (options.at(iArg).compare("-algorithm") == 0) {
       char tmpChar = options.at(iArg+1)[0];
       if (iArg+1 == argc || tmpChar == '-' ) {
         std::cout << " -treeName should be followed by a tree name" << std::endl;
         return 1;
       } else {
         algoId = strtol(options.at(iArg+1).c_str(), NULL, 0);
         iArg += 2;
       }
       
       
    } else{
      std::cout << "Couldn't understand argument " << options.at(iArg) << std::endl;
      return 1;
    }

  }//while arguments

  bool f_grid = false;
  bool f_lxbatch = false;  
  

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
      if( samplePath.find("grid") != std::string::npos ) //It is samples for the grid
        f_grid = true;

      std::ifstream inFile( samplePath );
      while(std::getline(inFile, containerName) ){
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
          }else{
            //Get full path of file
            char fullPath[300];
            realpath( containerName.c_str(), fullPath );
	    std::string thisPath = fullPath;
            //split into fileName and directory two levels above file
	    std::string fileName = thisPath.substr(containerName.find_last_of("/")+1);
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
      std::string thisPath = fullPath;
      //split into fileName and directory two levels above file
      std::string fileName = thisPath.substr(thisPath.find_last_of("/")+1);
      thisPath = thisPath.substr(0, thisPath.find_last_of("/"));
      
      std::cout << "path and file " << thisPath << " and " << fileName << std::endl;
      SH::ScanDir().sampleDepth(0).samplePattern(fileName).scan(sh, thisPath);
      
    }
  }//it's a file  


  ///////// Set output container name //////////////
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
  
  sh.setMetaString( "nc_tree", treeName );
  sh.setMetaString( "nc_grid_filter", "*");  //Data files on grid to not end in .root
  sh.print();
  
  job.sampleHandler( sh );
  // To automatically delete submitDir
  job.options()->setDouble(EL::Job::optRemoveSubmitDir, 1);

  //
  //  Set the number of events
  //
  TEnv* config = new TEnv(gSystem->ExpandPathName( configName.c_str() ));
  int nEvents = config->GetValue("MaxEvent",       -1);
  if(nEvents > 0)
    job.options()->setDouble(EL::Job::optMaxEvents, nEvents);
  
  //
  // Now the Event/Objection Selection
  //
  ZJetBalanceMiniTreeAnaSkeleton* analysisSkeleton = new ZJetBalanceMiniTreeAnaSkeleton();
  analysisSkeleton->setName("analysisSkeleton")->setConfig( "$ROOTCOREBIN/data/ZJetBalance/ZJetBalanceMiniTreeAnaSkeleton.config" );

  ZJetBalanceMiniTree_GenBalanceHistograms* analysisGenHistograms = new ZJetBalanceMiniTree_GenBalanceHistograms();
  analysisGenHistograms->setName("BalanceHistograms")->setConfig( "$ROOTCOREBIN/data/ZJetBalance/ZJetBalanceMiniTree_GenBalanceHistograms.config" );
  
  
  //
  // Add configured algos to event loop job
  //
  // !! UNFORTUNATELY, currently we cannot several algorithm in a singl job !!
  // (owing to how to access to TTree, to be modified in future)
  switch (algoId) {
  case 1:
    job.algsAdd( analysisSkeleton );
    break;
  case 2:
    job.algsAdd( analysisGenHistograms );
    break;
  default:
    Error("runZJetBalanceMiniTreeAna", "NO ALGORITHM IS SELECTED WITH ALGORITHM ID=%d", algoId);
    exit(EXIT_FAILURE);
  }  

  
  //
  // Submit the job
  //
  
  if(f_grid){
    EL::PrunDriver driver;
    driver.options()->setString("nc_outputSampleName", outputName);

    driver.options()->setDouble(EL::Job::optGridNFilesPerJob, 2);
    //driver.options()->setString(EL::Job::optRootVer, "5.34.25");
    //driver.options()->setString(EL::Job::optCmtConfig, "x86_64-slc6-gcc48-opt");
    //driver.options()->setDouble("nc_nGBPerJob", 1);
    //driver.options()->setString("nc_excludeSite", ???);
    //driver.options()->setString(EL::Job::optGridMergeOutput, "false");
    driver.options()->setDouble(EL::Job::optGridMemory,10240); // 10 GB

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
