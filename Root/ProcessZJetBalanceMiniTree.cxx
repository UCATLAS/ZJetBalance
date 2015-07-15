#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include <AthContainers/ConstDataVector.h>

#include <xAODTracking/VertexContainer.h>
#include <xAODJet/JetContainer.h>
#include <xAODEventInfo/EventInfo.h>
#include <ZJetBalance/ProcessZJetBalanceMiniTree.h>
#include <xAODAnaHelpers/HelperFunctions.h>
#include <xAODAnaHelpers/tools/ReturnCheck.h>

#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TEnv.h>
#include <TSystem.h>

#include <utility>      
#include <iostream>
#include <fstream>

using namespace std;

// this is needed to distribute the algorithm to the workers
ClassImp(ProcessZJetBalanceMiniTree)

ProcessZJetBalanceMiniTree :: ProcessZJetBalanceMiniTree () :
  m_debug(false)
{
  Info("ProcessZJetBalanceMiniTree()", "Calling constructor");
}


EL::StatusCode  ProcessZJetBalanceMiniTree :: configure ()
{
  m_configName = gSystem->ExpandPathName( m_configName.c_str() );
  Info("configure()", "Configuing ProcessZJetBalanceMiniTree Interface. User configuration read from : %s \n", m_configName.c_str());
  TEnv* config = new TEnv(m_configName.c_str());
  if( !config ) {
    Error("configure()", "Failed to read config file!");
    Error("configure()", "config name : %s",m_configName.c_str());
    return EL::StatusCode::FAILURE;
  }

  //
  // Read Input from .config file
  //
  m_debug                    = config->GetValue("Debug" ,          false );
  
  config->Print();
  Info("configure()", "ProcessZJetBalanceMiniTree Interface succesfully configured! \n");
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ProcessZJetBalanceMiniTree :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessZJetBalanceMiniTree :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  Info("histInitialize()", "Calling histInitialize \n");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessZJetBalanceMiniTree :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessZJetBalanceMiniTree :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  TTree *tree = wk()->tree();
  tree->SetBranchStatus ("*", 0);
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessZJetBalanceMiniTree :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  m_eventCounter = 0;
  
  Info("initialize()", "Succesfully initialized! \n");
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ProcessZJetBalanceMiniTree :: execute ()
{
  // Here you do everything that needs to be done on every single
  // event, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  if(m_debug) Info("execute()", "Processing Event");
  ++m_eventCounter;
  if (m_eventCounter%100==0) {
    Info("execute()", "%10d th event is been processed.", m_eventCounter);
  }
  
  wk()->tree()->GetEntry (wk()->treeEntry());
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ProcessZJetBalanceMiniTree :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessZJetBalanceMiniTree :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessZJetBalanceMiniTree :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}
