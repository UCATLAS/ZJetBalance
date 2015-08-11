#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include <AthContainers/ConstDataVector.h>

#include <xAODTracking/VertexContainer.h>
#include <xAODJet/JetContainer.h>
#include <xAODEventInfo/EventInfo.h>
#include <ZJetBalance/ZJetBalanceMiniTreeAnaSkeleton.h>
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

#include <stdlib.h>

using namespace std;

// this is needed to distribute the algorithm to the workers
ClassImp(ZJetBalanceMiniTreeAnaSkeleton)

ZJetBalanceMiniTreeAnaSkeleton :: ZJetBalanceMiniTreeAnaSkeleton ()
{
  Info("ZJetBalanceMiniTreeAnaSkeleton()", "Calling constructor");
}

EL::StatusCode  ZJetBalanceMiniTreeAnaSkeleton :: configure ()
{
  Info("configure()", "called");
  
  // load basic configuration
  // m_debug m_doPUreweighting m_lumiCalcFileNames m_PRWFileNames m_additional_weight m_btagJets m_btagOP
  if (this->LoadBasicConfiguration()==EL::StatusCode::FAILURE) {
    Error("configure()", "failed in LoadBasicConfiguration()");
    return EL::StatusCode::FAILURE;
  }
  
  
  
  //
  // Read Input from .config file
  //
  m_configName = gSystem->ExpandPathName( m_configName.c_str() );
  Info("configure()", "Configuing ZJetBalanceMiniTreeAnaSkeleton Interface. User configuration read from : %s \n", m_configName.c_str());
  TEnv* config = new TEnv(m_configName.c_str());
  if( !config ) {
    Error("configure()", "Failed to read config file!");
    Error("configure()", "config name : %s",m_configName.c_str());
    return EL::StatusCode::FAILURE;
  }
  
  if( m_btagJets ) {
    if( m_btagOP != "85" && m_btagOP != "77" && m_btagOP != "70" && m_btagOP != "60" ) {
      std::cout << "Invalid b-tag operating point " << m_btagOP << std::endl;
      return EL::StatusCode::FAILURE;
    }
    std::cout << "BTagging using the " << m_btagOP << "% operating point" << std::endl;
  }
  
  
  config->Print();
  Info("configure()", "ZJetBalanceMiniTreeAnaSkeleton Interface succesfully configured! \n");
  
  Info("configure()", "Ends");
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ZJetBalanceMiniTreeAnaSkeleton :: setupJob (EL::Job& job)
{
  Info("setupJob()", "called");
  
  Info("setupJob()", "ends");
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTreeAnaSkeleton :: histInitialize ()
{
  Info("histInitialize()", "called");
  
  // configuration from ENV file (note : histInistialize called before initialize)
  if ( this->configure() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  
  this->InitializeCutflows();
  
  // add your own histgrams
  m_h_ZM = this->book("ZM", "M_{Z} [GeV]", 60, 60, 120);
  
  Info("histInitialize()", "ends");
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTreeAnaSkeleton :: fileExecute ()
{
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ZJetBalanceMiniTreeAnaSkeleton :: initialize ()
{
  Info("initialize()", "Calling initialize");
  
  m_isMC          =  false;
  m_eventCounter  =  0;
  mcChannelNumber = -1; // needs for isMC decision
  
  // Pileup RW Tool //
  if (this->InitializePileupReweightingTool() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed in pileup RW tool initialization" );
    return EL::StatusCode::FAILURE;      
  }
  
  
  // pass all histograms instanted with the book function are passed
  // to the worker through 
  this->record( wk() );
  
  Info("initialize()", "Succesfully initialized! \n");
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ZJetBalanceMiniTreeAnaSkeleton :: execute ()
{
  this->LoadMiniTree();
  
  double weight_final=1.0;
  if (m_isMC) {
    double pileup_reweighting_factor = GetPileupReweightingFactor();
    weight_final = 
      mcEventWeight*weight_xs*pileup_reweighting_factor*m_additional_weight;
  }
  
  FillCutflowHistograms("Process NTuple", mcEventWeight, weight_final);
  m_h_ZM->Fill(ZM, weight_final);
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ZJetBalanceMiniTreeAnaSkeleton :: postExecute ()
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTreeAnaSkeleton :: finalize ()
{
  std::cout << "Finialize!" << std::endl;
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTreeAnaSkeleton :: histFinalize ()
{
  return EL::StatusCode::SUCCESS;
}

