#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include "EventLoop/OutputStream.h"
#include "AthContainers/ConstDataVector.h"

#include "xAODTracking/VertexContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEventInfo/EventInfo.h"
#include <ZJetBalance/BalanceAlgorithm.h>
#include <xAODAnaHelpers/HelperFunctions.h>
#include <xAODAnaHelpers/tools/ReturnCheck.h>

#include "TFile.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>

using namespace std;

// this is needed to distribute the algorithm to the workers
ClassImp(BalanceAlgorithm)

BalanceAlgorithm :: BalanceAlgorithm () :
  m_cutflowHist(0),
  m_cutflowHistW(0),
  m_GeV(1000.),
  m_treeStream("tree")
{
//  // Here you put any code for the base initialization of variables,
//  // e.g. initialize all pointers to 0.  Note that you should only put
//  // the most basic initialization here, since this method will be
//  // called on both the submission and the worker node.  Most of your
//  // initialization code will go into histInitialize() and
//  // initialize().

  Info("BalanceAlgorithm()", "Calling constructor");

}

EL::StatusCode  BalanceAlgorithm :: configure ()
{
  Info("configure()", "Configuing BalanceAlgorithm Interface. User configuration read from : %s \n", m_configName.c_str());
  m_configName = gSystem->ExpandPathName( m_configName.c_str() );
  TEnv* config = new TEnv(m_configName.c_str());
  if( !config ) {
    Error("configure()", "Failed to read config file!");
    Error("configure()", "config name : %s",m_configName.c_str());
    return EL::StatusCode::FAILURE;
  }

  // read debug flag from .config file
  m_useMuons                 = config->GetValue("UseMuons", true);
  m_inputJetContainerName    = config->GetValue("InputJetContainer",  "");
  m_inputJetAlgo             = config->GetValue("InputJetAlgo",       "");

  m_inputMuonForMuonInJetCorrectionContainerName   = config->GetValue("InputMuonForMuonInJetCorrectionContainer",  "");
  m_inputMuonForMuonInJetCorrectionAlgo            = config->GetValue("InputMuonForMuonInJetCorrectionAlgo",       "");
  m_debug                    = config->GetValue("Debug" ,      false );
  m_useCutFlow               = config->GetValue("UseCutFlow",  true);
  m_writeTree                = config->GetValue("WriteTree",  true);
  m_leadingJetPtCut          = config->GetValue("LeadingJetPtCut",  460000); //Default 400GeV
  m_truthLevelOnly           = config->GetValue("TruthLevelOnly",  false);
  m_eventDetailStr           = config->GetValue("EventDetailStr", "truth pileup");
  m_trigDetailStr            = config->GetValue("TrigDetailStr", "");
  m_jetDetailStr             = config->GetValue("JetDetailStr", "kinematic clean energy truth flavorTag layer");
  m_jetDetailStrSyst         = config->GetValue("JetDetailStrSyst", "kinematic clean energy");
  
  if( m_useMuons ){
    m_inputMuonContainerName   = config->GetValue("InputMuonContainer",  "");
    m_inputMuonAlgo            = config->GetValue("InputMuonAlgo",       "");
    m_muonDetailStr            = config->GetValue("MuonDetailStr", "kinematic");
  } else{
    m_inputElectronContainerName   = config->GetValue("InputElectronContainer",  "");
    m_inputElectronAlgo            = config->GetValue("InputElectronAlgo",       "");
    m_electronDetailStr        = config->GetValue("ElectronDetailStr", "kinematic");
  }

  config->Print();
  Info("configure()", "BalanceAlgorithm Interface succesfully configured! \n");

  if( m_inputJetContainerName.empty() ) {
    Error("configure()", "Jet InputContainer is empty!");
    return EL::StatusCode::FAILURE;
  }
  if( m_useMuons && m_inputMuonContainerName.empty() ) {
    Error("configure()", "Muon InputContainer is empty!");
    return EL::StatusCode::FAILURE;
  }
  if( (!m_useMuons) && m_inputElectronContainerName.empty() ) {
    Error("configure()", "Electron InputContainer is empty!");
    return EL::StatusCode::FAILURE;
  }
  if( m_inputMuonForMuonInJetCorrectionContainerName.empty() ) {
    Error("configure()", "Muon InputContainer is empty!");
    return EL::StatusCode::FAILURE;
  }

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode BalanceAlgorithm :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  job.useXAOD();
  xAOD::Init( "BalanceAlgorithm" ).ignore(); // call before opening first file

  EL::OutputStream outForTree( m_treeStream );
  job.outputAdd (outForTree);
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode BalanceAlgorithm :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  Info("histInitialize()", "Calling histInitialize \n");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode BalanceAlgorithm :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode BalanceAlgorithm :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode BalanceAlgorithm :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();
  m_eventCounter = -1;

  if ( this->configure() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("BalanceAlgorithm::initialize()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug), "");
  if( m_truthLevelOnly ) { m_isMC = true; }
  else { 
    m_isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) ? true : false;
  }
  
  Info("initialize()", "eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION )=%s", eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION )? "true" : "false");
  
  getLumiWeights(eventInfo);

  if(m_useCutFlow) {

    TFile *file = wk()->getOutputFile ("cutflow");
    m_cutflowHist  = (TH1D*)file->Get("cutflow");
    m_cutflowHistW = (TH1D*)file->Get("cutflow_weighted");

    m_cutflowFirst = m_cutflowHist->GetXaxis()->FindBin("LeadingJetPtCut");
    m_cutflowHistW->GetXaxis()->FindBin("LeadingJetPtCut");

    m_cutflowHist->GetXaxis()->FindBin("OppositeSignLeptonCut");
    m_cutflowHistW->GetXaxis()->FindBin("OppositeSignLeptonCut");

    m_cutflowHist->GetXaxis()->FindBin("ZMassLoose");
    m_cutflowHistW->GetXaxis()->FindBin("ZMassLoose");

  }

  Info("initialize()", "Succesfully initialized! \n");
  return EL::StatusCode::SUCCESS;
}

void BalanceAlgorithm::AddTree( std::string name ) {

  std::string treeName("outTree");
  // naming convention
  //if( ! name.empty() ) { treeName += "."; } // makes it hard to do command line stuff
  treeName += name; // add systematic
  TTree * outTree = new TTree(treeName.c_str(),treeName.c_str());
  if( !outTree ) {
    Error("AddTree()","Failed to get output tree!");
    // FIXME!! kill here
  }
  TFile* treeFile = wk()->getOutputFile( m_treeStream );
  outTree->SetDirectory( treeFile );

  MiniTree* miniTree = new MiniTree(m_event, outTree, treeFile);
  // only limited information available in truth xAODs
  if( m_truthLevelOnly ) {
    miniTree->AddEvent("truth");
    miniTree->AddJets("kinematic");
  } else { // reconstructed xAOD
    miniTree->AddEvent( m_eventDetailStr );
    miniTree->AddTrigger( m_trigDetailStr );
    if( !name.empty() ) { // save limited information for systematic variations
      miniTree->AddJets( m_jetDetailStrSyst );
      //miniTree->AddMuons( m_muonDetailStrSyst );
    } else {
      miniTree->AddJets ( m_jetDetailStr );
      if( m_useMuons )
        miniTree->AddMuons( m_muonDetailStr );
      else 
        miniTree->AddElectrons( m_electronDetailStr );
    }
  }
  m_myTrees[name] = miniTree;
  // see Worker.cxx line 134: the following function call takes ownership of the tree
  // from the treeFile; no output is written. so don't do that!
  // wk()->addOutput( outTree );

}



EL::StatusCode BalanceAlgorithm :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  if(m_debug) Info("execute()", "Applying selection");
  ++m_eventCounter;

  m_iCutflow = m_cutflowFirst;
  
  //----------------------------
  // Event information
  //---------------------------

  ///////////////////////////// Retrieve Containers /////////////////////////////////////////
  
  if(m_debug) Info("execute()", "Get Containers");
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("BalanceAlgorithm::execute()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug), "");

  if (m_eventCounter == 0) {
  }
  SG::AuxElement::ConstAccessor<float> NPVAccessor("NPV");
  const xAOD::VertexContainer* vertices = 0;
  if(!m_truthLevelOnly) {
    RETURN_CHECK("BalanceAlgorithm::execute()", HelperFunctions::retrieve(vertices, "PrimaryVertices", m_event, m_store), "");
  }
  if(!m_truthLevelOnly && !NPVAccessor.isAvailable( *eventInfo )) { // NPV might already be available
    // number of PVs with 2 or more tracks
    //eventInfo->auxdecor< int >( "NPV" ) = HelperFunctions::countPrimaryVertices(vertices, 2);
    // TMP for JetUncertainties uses the same variable
    eventInfo->auxdecor< float >( "NPV" ) = HelperFunctions::countPrimaryVertices(vertices, 2);
  }
  //Set this first, as it's needed by passCut()
  if(m_isMC)
    m_mcEventWeight = eventInfo->mcEventWeight();
  else
    m_mcEventWeight = 1;

  eventInfo->auxdecor< float >("weight_xs") = m_xs * m_filtEff;
  eventInfo->auxdecor< float >("weight") = m_mcEventWeight * m_xs * m_filtEff;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loop over Systematics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(m_debug) Info("execute()", "Systematic Loop");
  // did any collection pass the cuts?
  bool pass(false);
  bool doCutflow(m_useCutFlow); // will only stay true for nominal
  const xAOD::JetContainer*  signalJets  = 0;
  
  // Declare both containers, so we have them in the scope.
  const xAOD::MuonContainer* signalMuons = 0; 
  const xAOD::ElectronContainer* signalElectrons = 0;

  // if input comes from xAOD, or just running one collection,
  // then get the one collection and be done with it
  if( (m_inputJetAlgo.empty() && 
      (( m_useMuons && m_inputMuonAlgo.empty() ) || ( !m_useMuons && m_inputElectronAlgo.empty() ))) 
      || m_truthLevelOnly ) {
    RETURN_CHECK("BalanceAlgorithm::execute()", HelperFunctions::retrieve(signalJets, m_inputJetContainerName, m_event, m_store), "");
    
    if( m_useMuons ){
      RETURN_CHECK("BalanceAlgorithm::execute()", HelperFunctions::retrieve(signalMuons, m_inputMuonContainerName, m_event, m_store), "");
      pass = this->executeAnalysis( eventInfo, signalJets, signalMuons, vertices, doCutflow, "" );
    }else{
      RETURN_CHECK("BalanceAlgorithm::execute()", HelperFunctions::retrieve(signalElectrons, m_inputElectronContainerName, m_event, m_store), "");
      pass = this->executeAnalysis( eventInfo, signalJets, signalElectrons, vertices, doCutflow, "" );
    }
  }
  else { // get the list of systematics to run over

    // get vector of string giving the names
    std::vector<std::string>* systNames = 0;
    if ( m_store->contains< std::vector<std::string> >( m_inputJetAlgo ) ) {
      if(!m_store->retrieve( systNames, m_inputJetAlgo ).isSuccess()) {
        Info("execute()", "Cannot find vector from %s", m_inputJetAlgo.c_str());
        return StatusCode::FAILURE;
      }
    }
    
    // Add loop over muon systematics !!! TODO
    if( m_useMuons )
      RETURN_CHECK("BalanceAlgorithm::execute()", HelperFunctions::retrieve(signalMuons, m_inputMuonContainerName, m_event, m_store), "");
    else
      RETURN_CHECK("BalanceAlgorithm::execute()", HelperFunctions::retrieve(signalElectrons, m_inputElectronContainerName, m_event, m_store), "");

    // loop over systematics
    bool saveContainerNames(false);
    std::vector< std::string >* vecOutContainerNames = 0;
    if(saveContainerNames) { vecOutContainerNames = new std::vector< std::string >; }
    // shoudl only doCutflow for the nominal
    bool passOne(false);
    std::string inputJetContainerName("");
    for( auto systName : *systNames ) {

      inputJetContainerName = m_inputJetContainerName+systName;
      RETURN_CHECK("BalanceAlgorithm::execute()", HelperFunctions::retrieve(signalJets, inputJetContainerName, m_event, m_store), "");
      // allign with Dijet naming conventions
      if( systName.empty() ) { doCutflow = m_useCutFlow; } // only doCutflow for nominal
      else { doCutflow = false; }

      if( m_useMuons )
        passOne = this->executeAnalysis( eventInfo, signalJets, signalMuons, vertices, doCutflow, systName );
      else
        passOne = this->executeAnalysis( eventInfo, signalJets, signalElectrons, vertices, doCutflow, systName );
      
      // save the string if passing the selection
      if( saveContainerNames && passOne ) { vecOutContainerNames->push_back( systName ); }
      // the final decision - if at least one passes keep going!
      pass = pass || passOne;

    }

    // save list of systs that shoudl be considered down stream
    if( saveContainerNames ) {
      RETURN_CHECK( "execute()", m_store->record( vecOutContainerNames, m_name), "Failed to record vector of output container names.");
    }
  }
  if(!pass) {
    wk()->skipEvent();
  }
  return EL::StatusCode::SUCCESS;
  
}

// We overload the executeAnalysis function to use muons or electrons, depending on the input container
bool BalanceAlgorithm :: executeAnalysis ( const xAOD::EventInfo* eventInfo,
    const xAOD::JetContainer* signalJets,
    const xAOD::MuonContainer* signalMuons,
    const xAOD::VertexContainer* vertices,
    bool doCutflow,
    std::string systName) {

  /////////////////////////// Begin Selections  ///////////////////////////////


  ////  leadingJetPt trigger efficiency BEFORE MCCLEANING ////
  if( signalJets->at(0)->pt() < m_leadingJetPtCut){
      wk()->skipEvent();  return EL::StatusCode::SUCCESS;
  }
  if(doCutflow) passCut(); //Leading jet pT cut

  if( ( (signalMuons->at(0)->charge()) * (signalMuons->at(1)->charge()) ) > 0 ){
      wk()->skipEvent();  return EL::StatusCode::SUCCESS;
  }
  if(doCutflow) passCut(); // Opposite sign leptons cut.

  // create the Z-Object
  TLorentzVector Z = TLorentzVector();
  Z += signalMuons->at(0)->p4();
  Z += signalMuons->at(1)->p4();
  if( fabs(Z.M() - 90e3) > 35e3 ) {
      wk()->skipEvent();  return EL::StatusCode::SUCCESS;
  }
  if(doCutflow) passCut(); //Z mass cut

  if ( muonInJetCorrection (signalJets) == EL::StatusCode::FAILURE ) {
    Error("executeAnalysis()", "failure in muonInJetCorrection()");
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Selections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ///////////////////////// Add final variables ////////////////////////////////
  eventInfo->auxdecor< float >( "ZpT" )  = Z.Pt() / m_GeV;
  eventInfo->auxdecor< float >( "Zeta" ) = Z.Eta();
  eventInfo->auxdecor< float >( "Zphi" ) = Z.Phi();
  eventInfo->auxdecor< float >( "ZM" )   = Z.M() / m_GeV;
  eventInfo->auxdecor< float >( "dPhiZJet1" ) = Z.DeltaPhi( signalJets->at(0)->p4() );
  eventInfo->auxdecor< float >( "dEtaZJet1" ) = Z.Eta() - signalJets->at(0)->eta();
  eventInfo->auxdecor< float >( "pTRef1" ) = Z.Pt()*fabs(cos( Z.DeltaPhi( signalJets->at(0)->p4() ) )) / m_GeV;
  if( signalJets->size() > 1 ){ 
    eventInfo->auxdecor< float >( "dPhiZJet2" ) = Z.DeltaPhi( signalJets->at(1)->p4() );
    eventInfo->auxdecor< float >( "dEtaZJet2" ) = Z.Eta() - signalJets->at(1)->eta();
    eventInfo->auxdecor< float >( "pTRef2" ) = Z.Pt()*fabs(cos( Z.DeltaPhi( signalJets->at(1)->p4() ) )) / m_GeV; 
    eventInfo->auxdecor< float >( "jetDPhi" ) = signalJets->at(0)->p4().DeltaPhi( signalJets->at(1)->p4() );
    eventInfo->auxdecor< float >( "jetDEta" ) = signalJets->at(0)->eta() - signalJets->at(1)->eta();
    eventInfo->auxdecor< float >( "jetPtRatio" ) = signalJets->at(1)->pt() / signalJets->at(0)->pt();
  }

  // detector eta and punch-through-variable
  if( !m_truthLevelOnly ) {
    xAOD::JetFourMom_t fourVec;
    bool status(false);
    for( auto iJet : *signalJets ) {
      status = iJet->getAttribute<xAOD::JetFourMom_t>( "JetEMScaleMomentum", fourVec );
      if( status ) { iJet->auxdecor< float >( "emScaleEta") = fourVec.eta(); }
      else         { iJet->auxdecor< float >( "emScaleEta") = -999; }
      status = iJet->getAttribute<xAOD::JetFourMom_t>( "JetConstitScaleMomentum", fourVec );
      if( status ) { iJet->auxdecor< float >( "constitScaleEta") = fourVec.eta(); }
      else         { iJet->auxdecor< float >( "constitScaleEta") = -999; }
    }
  }


  if(m_debug){
    std::cout << "Event # " << m_eventCounter << std::endl;
  }

  /////////////////////////////////////// Output Plots ////////////////////////////////
  float weight(1);
  if( eventInfo->isAvailable< float >( "weight" ) ) {
    weight = eventInfo->auxdecor< float >( "weight" );
  }

  ///////////////////////////// fill the tree ////////////////////////////////////////////
  if(m_writeTree){
    if( m_myTrees.find( systName ) == m_myTrees.end() ) { AddTree( systName ); }
    if(eventInfo)   m_myTrees[systName]->FillEvent( eventInfo, m_event );
    if( m_truthLevelOnly ) {
      if(signalJets)  m_myTrees[systName]->FillJets( signalJets, -1 );
    } else {
      m_myTrees[systName]->FillTrigger( eventInfo );
      int pVLoc = HelperFunctions::getPrimaryVertexLocation( vertices );
      if(signalJets)  m_myTrees[systName]->FillJets( signalJets, pVLoc );
      if(signalMuons) m_myTrees[systName]->FillMuons( signalMuons, vertices->at(pVLoc) );
    }
    m_myTrees[systName]->Fill();
  }

  return true;
}

bool BalanceAlgorithm :: executeAnalysis ( const xAOD::EventInfo* eventInfo,
    const xAOD::JetContainer* signalJets,
    const xAOD::ElectronContainer* signalElectrons,
    const xAOD::VertexContainer* vertices,
    bool doCutflow,
    std::string systName) {

  /////////////////////////// Begin Selections  ///////////////////////////////


  ////  leadingJetPt trigger efficiency BEFORE MCCLEANING ////
  if( signalJets->at(0)->pt() < m_leadingJetPtCut){
      wk()->skipEvent();  return EL::StatusCode::SUCCESS;
  }
  if(doCutflow) passCut(); //Leading jet pT cut
  
  if( ( (signalElectrons->at(0)->charge()) * (signalElectrons->at(1)->charge()) ) > 0 ){
      wk()->skipEvent();  return EL::StatusCode::SUCCESS;
  }
  if(doCutflow) passCut(); // Opposite sign leptons cut.
  
  // create the Z-Object
  TLorentzVector Z = TLorentzVector();
  Z += signalElectrons->at(0)->p4();
  Z += signalElectrons->at(1)->p4();
  if( fabs(Z.M() - 90e3) > 35e3 ) {
      wk()->skipEvent();  return EL::StatusCode::SUCCESS;
  }
  if(doCutflow) passCut(); //Z mass cut

  if ( muonInJetCorrection (signalJets) == EL::StatusCode::FAILURE ) {
    Error("executeAnalysis()", "failure in muonInJetCorrection()");
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Selections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ///////////////////////// Add final variables ////////////////////////////////
  eventInfo->auxdecor< float >( "ZpT" )  = Z.Pt() / m_GeV;
  eventInfo->auxdecor< float >( "Zeta" ) = Z.Eta();
  eventInfo->auxdecor< float >( "Zphi" ) = Z.Phi();
  eventInfo->auxdecor< float >( "ZM" )   = Z.M() / m_GeV;
  eventInfo->auxdecor< float >( "dPhiZJet1" ) = Z.DeltaPhi( signalJets->at(0)->p4() );
  eventInfo->auxdecor< float >( "pTRef1" ) = Z.Pt()*fabs(cos( Z.DeltaPhi( signalJets->at(0)->p4() ) )) / m_GeV;
  if( signalJets->size() > 1 ){ 
    eventInfo->auxdecor< float >( "dPhiZJet2" ) = Z.DeltaPhi( signalJets->at(1)->p4() );
    eventInfo->auxdecor< float >( "pTRef2" ) = Z.Pt()*fabs(cos( Z.DeltaPhi( signalJets->at(1)->p4() ) )) / m_GeV; 
    eventInfo->auxdecor< float >( "jetDPhi" ) = signalJets->at(0)->p4().DeltaPhi( signalJets->at(1)->p4() );
    eventInfo->auxdecor< float >( "jetDEta" ) = signalJets->at(0)->eta() - signalJets->at(1)->eta();
    eventInfo->auxdecor< float >( "jetPtRatio" ) = signalJets->at(1)->pt() / signalJets->at(0)->pt();
  }

  // detector eta and punch-through-variable
//  if( !m_truthLevelOnly ) {
//    for( auto iJet : *signalJets ) {
//      xAOD::JetFourMom_t jetConstitScaleP4 = iJet->getAttribute<xAOD::JetFourMom_t>("JetConstitScaleMomentum");
//      xAOD::JetFourMom_t jetEMScaleP4 = iJet->getAttribute<xAOD::JetFourMom_t>("JetEMScaleMomentum");
//      iJet->auxdecor< float >( "constitScaleEta") = jetConstitScaleP4.eta();
//      iJet->auxdecor< float >( "emScaleEta")      = jetEMScaleP4.eta();
//    }
//  }


  if(m_debug){
    std::cout << "Event # " << m_eventCounter << std::endl;
  }

  /////////////////////////////////////// Output Plots ////////////////////////////////
  float weight(1);
  if( eventInfo->isAvailable< float >( "weight" ) ) {
    weight = eventInfo->auxdecor< float >( "weight" );
  }


  ///////////////////////////// fill the tree ////////////////////////////////////////////
  if(m_writeTree){
    if( m_myTrees.find( systName ) == m_myTrees.end() ) { AddTree( systName ); }
    if(eventInfo)   m_myTrees[systName]->FillEvent( eventInfo, m_event );
    if( m_truthLevelOnly ) {
      if(signalJets)  m_myTrees[systName]->FillJets( signalJets, -1 );
    } else {
      m_myTrees[systName]->FillTrigger( eventInfo );
      int pVLoc = HelperFunctions::getPrimaryVertexLocation( vertices );
      if(signalJets)  m_myTrees[systName]->FillJets( signalJets, pVLoc );
      if(signalElectrons) m_myTrees[systName]->FillElectrons( signalElectrons, vertices->at(pVLoc) );
    }
    m_myTrees[systName]->Fill();
  }

  return true;
}

//Easy method for automatically filling cutflow and incrementing counter
void BalanceAlgorithm::passCut() {
  m_cutflowHist->Fill(m_iCutflow, 1);
  m_cutflowHistW->Fill(m_iCutflow, m_mcEventWeight);
  m_iCutflow++;
}

//This grabs cross section, acceptance, and eventNumber information from the respective text file
//text format:     147915 2.3793E-01 5.0449E-03 499000
EL::StatusCode BalanceAlgorithm::getLumiWeights(const xAOD::EventInfo* eventInfo) {
  
  Info("getLumiWeights()", "m_isMC=%s", m_isMC? "true" : "false");
  
  if(!m_isMC){
    m_mcChannelNumber = eventInfo->runNumber();
    m_xs = 1;
    m_filtEff = 1;
    m_numAMIEvents = 0;
    return EL::StatusCode::SUCCESS;
  }

  m_mcChannelNumber = eventInfo->mcChannelNumber();
  //if mcChannelNumber = 0 need to retrieve from runNumber
  if(eventInfo->mcChannelNumber()==0) m_mcChannelNumber = eventInfo->runNumber();
  ifstream fileIn(  gSystem->ExpandPathName( "$ROOTCOREBIN/data/ZJetBalance/XsAcc_13TeV.txt" ) );
  std::string runNumStr = std::to_string( m_mcChannelNumber );
  std::string line;
  std::string subStr;
  while (getline(fileIn, line)){
    istringstream iss(line);
    iss >> subStr;
    if (subStr.find(runNumStr) != string::npos){
      iss >> subStr;
      sscanf(subStr.c_str(), "%e", &m_xs);
      iss >> subStr;
      sscanf(subStr.c_str(), "%e", &m_filtEff);
      iss >> subStr;
      sscanf(subStr.c_str(), "%i", &m_numAMIEvents);
      cout << "Setting xs / acceptance / numAMIEvents to " << m_xs << ":" << m_filtEff << ":" << m_numAMIEvents << endl;
      continue;
    }
  }
  if( m_numAMIEvents == 0){
    cerr << "ERROR: Could not find proper file information for file number " << runNumStr << endl;
    return EL::StatusCode::FAILURE;
  }
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode BalanceAlgorithm :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode BalanceAlgorithm :: finalize ()
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



EL::StatusCode BalanceAlgorithm :: histFinalize ()
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
  if( m_writeTree ) {
    std::string thisName;
    m_ss.str( std::string() );
    m_ss << m_mcChannelNumber;
    TFile * treeFile = wk()->getOutputFile( m_treeStream );
    if(m_useCutFlow) {
      TH1F* thisCutflowHist = (TH1F*) m_cutflowHist->Clone();
      thisName = thisCutflowHist->GetName();
      thisCutflowHist->SetName( (thisName+"_"+m_ss.str()).c_str() );
      thisCutflowHist->SetDirectory( treeFile );

      TH1F* thisCutflowHistW = (TH1F*) m_cutflowHistW->Clone();
      thisName = thisCutflowHistW->GetName();
      thisCutflowHistW->SetName( (thisName+"_"+m_ss.str()).c_str() );
      thisCutflowHistW->SetDirectory( treeFile );
    }
  }
  return EL::StatusCode::SUCCESS;
}

//===============================
EL::StatusCode BalanceAlgorithm :: muonInJetCorrection (const xAOD::JetContainer* signalJets)
{
  const xAOD::MuonContainer* muonsForCorr = 0;
  RETURN_CHECK("BalanceAlgorithm::muonInJetCorrection()", 
	       HelperFunctions::retrieve(muonsForCorr, m_inputMuonForMuonInJetCorrectionContainerName, m_event, m_store), "");
  
  for( auto signalJet : *signalJets ) {
    const TLorentzVector& jetP4 = signalJet->p4();
    
    double minimumDr = 0.4;
    const xAOD::Muon* closestMuon = 0;
    for ( auto muon : *muonsForCorr ) {
      const TLorentzVector& muonP4 = muon->p4();
      const double dR = jetP4.DeltaR(muonP4);
      if (dR<minimumDr) {
	minimumDr   = dR;
	closestMuon = muon;
      }
    }
    
    TLorentzVector muonInJetP4(0, 0, 0, 0);
    if (closestMuon) {
      muonInJetP4 = getFourMomentumOfMuonInJet (closestMuon);
    }
    
    const TLorentzVector correctedJetP4 = (jetP4+muonInJetP4);
    
    signalJet->auxdecor< float >("mucorrected_pt")  = correctedJetP4.Pt();
    signalJet->auxdecor< float >("mucorrected_phi") = correctedJetP4.Phi();
    signalJet->auxdecor< float >("mucorrected_eta") = correctedJetP4.Eta();
    signalJet->auxdecor< float >("mucorrected_m"  ) = correctedJetP4.M();
  }
  
  return EL::StatusCode::SUCCESS;
}

//===============================
TLorentzVector BalanceAlgorithm :: getFourMomentumOfMuonInJet (const xAOD::Muon* muon)  
{
  float eLoss=0.0;
  muon->parameter(eLoss,xAOD::Muon::EnergyLoss);
  const TLorentzVector& muonP4 = muon->p4();
  
  const double theta=muonP4.Theta();
  const double phi  =muonP4.Phi();
  
  const double eLossX=eLoss*sin(theta)*cos(phi);
  const double eLossY=eLoss*sin(theta)*sin(phi);
  const double eLossZ=eLoss*cos(theta);
  
  const TLorentzVector eLossP4(eLossX,eLossY,eLossZ,eLoss);  
  
  return muonP4-eLossP4;
}
