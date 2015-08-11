#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include <AthContainers/ConstDataVector.h>

#include <xAODTracking/VertexContainer.h>
#include <xAODJet/JetContainer.h>
#include <xAODEventInfo/EventInfo.h>
#include <ZJetBalance/ZJetBalanceMiniTreeAnaBase.h>
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
ClassImp(ZJetBalanceMiniTreeAnaBase)

ZJetBalanceMiniTreeAnaBase :: ZJetBalanceMiniTreeAnaBase ()
{
  Info("ZJetBalanceMiniTreeAnaBase()", "Calling constructor");
}

EL::StatusCode  ZJetBalanceMiniTreeAnaBase :: configure ()
{
  Info("configure()", "Called");
  if (this->LoadBasicConfiguration()==EL::StatusCode::FAILURE) {
    Error("configure()", "failed in LoadBasicConfiguration()");
    return EL::StatusCode::FAILURE;
  }
  
  Info("configure()", "Ends");
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ZJetBalanceMiniTreeAnaBase :: setupJob (EL::Job& job)
{
  Info("setupJob()", "called");
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  
  Info("setupJob()", "ends");
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTreeAnaBase :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  Info("histInitialize()", "called");
  this->InitializeCutflows();
  // histogram record will be done at the end of initialize()
  
  Info("histInitialize()", "ends");
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTreeAnaBase :: fileExecute ()
{
  Info("fileExecute()", "called");
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  Info("fileExecute()", "end");
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTreeAnaBase :: changeInput (bool firstFile)
{
  Info("changedInput", "called"); 
  
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  
  const std::pair<TH1F*, TH1F*> cutflows = ReturnCutflowPointers();  
  // if first time through - need to set bin labels
  if ( m_h_cutflow->GetXaxis()->GetBinLabel(1) != cutflows.first->GetXaxis()->GetBinLabel(1) ) {
    std::cout << "Labels do not agree" << std::endl;
    for (int iBin=1, nBins=cutflows.first->GetNbinsX(); iBin<=nBins; iBin++) {
      if( std::string(cutflows.first->GetXaxis()->GetBinLabel(iBin)) == "" ) { continue; }
      m_h_cutflow->GetXaxis()->SetBinLabel( iBin, 
          cutflows.first->GetXaxis()->GetBinLabel(iBin) );
      m_h_cutflow_weighted->GetXaxis()->SetBinLabel( iBin, 
          cutflows.first->GetXaxis()->GetBinLabel(iBin) );
      m_h_cutflow_weighted_final->GetXaxis()->SetBinLabel( iBin, 
          cutflows.first->GetXaxis()->GetBinLabel(iBin) );
      std::cout << "\t" << cutflows.first->GetXaxis()->GetBinLabel(iBin)  << std::endl;
    }
  }
  
  for (int iBin=1, nBins=m_h_cutflow->GetNbinsX(); iBin<=nBins; iBin++) { // update
    m_h_cutflow->SetBinContent(iBin, m_h_cutflow->GetBinContent(iBin)+cutflows.first->GetBinContent(iBin)); 
    Info("changeInput()", "iBin=%2d/%-2d new cutflow entry=%.1f (%s)", iBin, nBins, m_h_cutflow->GetBinContent(iBin), m_h_cutflow->GetXaxis()->GetBinLabel(iBin));
  }
  for (int iBin=1, nBins=m_h_cutflow_weighted->GetNbinsX(); iBin<=nBins; iBin++) { // update
    m_h_cutflow_weighted->SetBinContent(iBin, m_h_cutflow_weighted->GetBinContent(iBin)+cutflows.second->GetBinContent(iBin)); 
    Info("changeInput()", "iBin=%2d/%-2d new cutflow entry=%.1f (%s) - weighted ", iBin, nBins, m_h_cutflow_weighted->GetBinContent(iBin), m_h_cutflow_weighted->GetXaxis()->GetBinLabel(iBin));
  }
  
  TTree *tree = wk()->tree();
  InitTree(tree);
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTreeAnaBase :: initialize ()
{
  Info("initialize()", "Calling initialize");

  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  m_isMC          =  false;
  m_eventCounter  =  0;
  mcChannelNumber = -1; // needs for isMC decision
  
  // configuration from ENV file
  if ( this->configure() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
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


EL::StatusCode ZJetBalanceMiniTreeAnaBase :: execute ()
{
  // Here you do everything that needs to be done on every single
  // event, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  this->LoadMiniTree();
  
  double weight_final=1.0;
  if (m_isMC) {
    double pileup_reweighting_factor = GetPileupReweightingFactor();
    weight_final = 
      mcEventWeight*weight_xs*pileup_reweighting_factor*m_additional_weight;
  }
  
  // use non-weighted to check that the correct number of events was processes
  // use weighted to see the effect of the weights calculated in xAH
  FillCutflowHistograms("Process NTuple", mcEventWeight, weight_final);
  
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ZJetBalanceMiniTreeAnaBase :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTreeAnaBase :: finalize ()
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
  std::cout << "Finialize!" << std::endl;
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTreeAnaBase :: histFinalize ()
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

void ZJetBalanceMiniTreeAnaBase::InitializeCutflows()
{
  const std::pair<TH1F*, TH1F*> cutflows = ReturnCutflowPointers();
  int nBinsCutflow = cutflows.first->GetNbinsX();
  int xminCutflow  = cutflows.first->GetXaxis()->GetXmin();
  int xmaxCutflow  = cutflows.first->GetXaxis()->GetXmax();
  
  m_h_cutflow = new TH1F(Form("%scutflow",m_histPrefix.c_str()), "", nBinsCutflow, xminCutflow, xmaxCutflow);
  m_h_cutflow->SetCanExtend( TH1::kXaxis );
  m_h_cutflow_weighted = new TH1F(Form("%scutflow_weighted",m_histPrefix.c_str()), "ONLY MC WEIGHT given by generator", nBinsCutflow, xminCutflow, xmaxCutflow);
  m_h_cutflow_weighted->SetCanExtend( TH1::kXaxis );
  m_h_cutflow_weighted_final = new TH1F(Form("%scutflow_weighted_final",m_histPrefix.c_str()), "", nBinsCutflow, xminCutflow, xmaxCutflow);
  m_h_cutflow_weighted_final->SetCanExtend( TH1::kXaxis );
  
  wk()->addOutput( m_h_cutflow );
  wk()->addOutput( m_h_cutflow_weighted );
  wk()->addOutput( m_h_cutflow_weighted_final );
}

EL::StatusCode ZJetBalanceMiniTreeAnaBase::InitializePileupReweightingTool()
{
  m_pileuptool = new CP::PileupReweightingTool("Pileup");
  
  std::vector<std::string> PRWFiles;
  std::vector<std::string> lumiCalcFiles;
    
  std::string tmp_lumiCalcFileNames = m_lumiCalcFileNames;
  std::string tmp_PRWFileNames = m_PRWFileNames;
    
  // Parse all comma seperated files
  while( tmp_PRWFileNames.size() > 0){
    std::string::size_type pos = tmp_PRWFileNames.find_first_of(',');
    if( pos == std::string::npos){
      pos = tmp_PRWFileNames.size();
      PRWFiles.push_back( gSystem->ExpandPathName( (tmp_PRWFileNames.substr(0, pos)).c_str() ) );
      tmp_PRWFileNames.erase(0, pos);
    }else{
      PRWFiles.push_back( gSystem->ExpandPathName( (tmp_PRWFileNames.substr(0, pos)).c_str() ) );
      tmp_PRWFileNames.erase(0, pos+1);
    }
  }
  while( tmp_lumiCalcFileNames.size() > 0){
    std::string::size_type pos = tmp_lumiCalcFileNames.find_first_of(',');
    if( pos == std::string::npos){
      pos = tmp_lumiCalcFileNames.size();
      lumiCalcFiles.push_back(tmp_lumiCalcFileNames.substr(0, pos));
      tmp_lumiCalcFileNames.erase(0, pos);
    }else{
      lumiCalcFiles.push_back(tmp_lumiCalcFileNames.substr(0, pos));
      tmp_lumiCalcFileNames.erase(0, pos+1);
    }
  }

  std::cout << "PileupReweighting Tool is adding Pileup files:" << std::endl;
  for( unsigned int i=0; i < PRWFiles.size(); ++i){
    std::cout << "    " << PRWFiles.at(i) << std::endl;
  }
  std::cout << "PileupReweighting Tool is adding Lumi Calc files:" << std::endl;
  for( unsigned int i=0; i < lumiCalcFiles.size(); ++i){
    std::cout << "    " << lumiCalcFiles.at(i) << std::endl;
  }
    
  RETURN_CHECK("ZJetBalance::initialize()", m_pileuptool->setProperty("ConfigFiles", PRWFiles), "");
  RETURN_CHECK("ZJetBalance::initialize()", m_pileuptool->setProperty("LumiCalcFiles", lumiCalcFiles), "");
  RETURN_CHECK("ZJetBalance::initialize()", m_pileuptool->initialize(), "");
  std::cout << "ZJetBalance::initialize() Initialized PileupReweightingTool " << m_pileuptool->name() << std::endl; 
  
  return EL::StatusCode::SUCCESS;
}

void ZJetBalanceMiniTreeAnaBase::LoadMiniTree()
{
  //===============================
  // load tree
  //===============================
  TTree* tree = wk()->tree();
  if(m_btagJets) { if(!jet_isBTag || !jet_SFBTag) { this->SetBTagAddresses(tree); } }
  const int jentry = wk()->treeEntry();
  tree->LoadTree (jentry);
  tree->GetEntry (jentry);
  
  m_isMC = (mcChannelNumber!=-1);  
}


double ZJetBalanceMiniTreeAnaBase::GetPileupReweightingFactor()
{
  // what actually done in getCombinedWeight( const xAOD::EventInfo& eventInfo )
  // PileupReweighting/Root/PileupReweightingTool.cxx
  return ((m_doPUreweighting) ? 
	  m_pileuptool->GetCombinedWeight(runNumber, mcChannelNumber, averageInteractionsPerCrossing) :
	  1.0);
}

std::pair<TH1F*, TH1F*> ZJetBalanceMiniTreeAnaBase::ReturnCutflowPointers()
{
  std::pair<TH1F*, TH1F*> rc(0, 0);
  
  TFile* inputFile = wk()->inputFile();
  TIter next(inputFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    std::string keyName = key->GetName();
    
    std::size_t found = keyName.find("cutflow");
    bool foundCutFlow = (found!=std::string::npos);
    
    found = keyName.find("weighted");
    bool foundWeighted = (found!=std::string::npos);
    
    if(foundCutFlow && !foundWeighted){
      rc.first  = ((TH1F*)key->ReadObj());
    } else if(foundCutFlow && foundWeighted) {
      rc.second = ((TH1F*)key->ReadObj());
    }
  }//over Keys
  
  // rebin
  int nValidBins=0;
  // if first time through - need to set bin labels
  for (int iBin=1, nBins=rc.first->GetNbinsX(); iBin<=nBins; iBin++) {
    if( std::string(rc.first->GetXaxis()->GetBinLabel(iBin)) == "" ) { continue; }
    nValidBins++;
  }
  
  // reset bins = in order to avoid bins without label set (needed for TH1::FindBin(label) function )
  Info("ReturnCutflowPointers()", "New Binning nBins=%d min=%f max=%f", nValidBins, 0.5, nValidBins+0.5);
  rc.first->SetBins(nValidBins, 0.5, nValidBins+0.5);
  rc.second->SetBins(nValidBins, 0.5, nValidBins+0.5);

  return rc;
}



TH1F* ZJetBalanceMiniTreeAnaBase::book(std::string name,
                             std::string xlabel, int xbins, double xlow, double xhigh)
{
  TH1F* tmp = new TH1F( (m_histPrefix+name).c_str(), name.c_str(), xbins, xlow, xhigh);
  SetLabel(tmp, xlabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TH2F* ZJetBalanceMiniTreeAnaBase::book(std::string name,
                             std::string xlabel, int xbins, double xlow, double xhigh,
                             std::string ylabel, int ybins, double ylow, double yhigh)
{
  TH2F* tmp = new TH2F( (m_histPrefix+name).c_str(), name.c_str(), xbins, xlow, xhigh, ybins, ylow, yhigh);
  SetLabel(tmp, xlabel, ylabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}


/////// Variable Binned Histograms ///////
TH1F* ZJetBalanceMiniTreeAnaBase::book(std::string name,
    std::string xlabel, int xbins, const Double_t* xbinArr)
{
  TH1F* tmp = new TH1F( (m_histPrefix+name).c_str(), name.c_str(), xbins, xbinArr);
  SetLabel(tmp, xlabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TH2F* ZJetBalanceMiniTreeAnaBase::book(std::string name,
    std::string xlabel, int xbins, const Double_t* xbinArr,
    std::string ylabel, int ybins, double ylow, double yhigh)
{
  TH2F* tmp = new TH2F( (m_histPrefix+name).c_str(), name.c_str(), xbins, xbinArr, ybins, ylow, yhigh);
  SetLabel(tmp, xlabel, ylabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TH2F* ZJetBalanceMiniTreeAnaBase::book(std::string name,
    std::string xlabel, int xbins, double xlow, double xhigh,
    std::string ylabel, int ybins, const Double_t* ybinArr)
{
  TH2F* tmp = new TH2F( (m_histPrefix+name).c_str(), name.c_str(), xbins, xlow, xhigh, ybins, ybinArr);
  SetLabel(tmp, xlabel, ylabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

TH2F* ZJetBalanceMiniTreeAnaBase::book(std::string name,
    std::string xlabel, int xbins, const Double_t* xbinArr,
    std::string ylabel, int ybins, const Double_t* ybinArr)
{
  TH2F* tmp = new TH2F( (m_histPrefix+name).c_str(), name.c_str(), xbins, xbinArr, ybins, ybinArr);
  SetLabel(tmp, xlabel, ylabel);
  this->Sumw2(tmp);
  this->record(tmp);
  return tmp;
}

void ZJetBalanceMiniTreeAnaBase::Sumw2(TH1* hist, bool flag /*=true*/) {
  hist->Sumw2(flag);
}

void ZJetBalanceMiniTreeAnaBase::SetLabel(TH1* hist, std::string xlabel)
{
  hist->GetXaxis()->SetTitle(xlabel.c_str());
}

void ZJetBalanceMiniTreeAnaBase::SetLabel(TH1* hist, std::string xlabel, std::string ylabel)
{
  hist->GetYaxis()->SetTitle(ylabel.c_str());
  this->SetLabel(hist, xlabel);
}

void ZJetBalanceMiniTreeAnaBase::record(TH1* hist) {
  m_allHists.push_back( hist );
}

void ZJetBalanceMiniTreeAnaBase::record(EL::Worker* wk) {
  for( auto hist : m_allHists ){
    Info("record()", "histogram record : %s", hist->GetName());
    wk->addOutput(hist);
  }
}

void ZJetBalanceMiniTreeAnaBase::FillCutflowHistograms(const std::string& label, const double& xAHWeight, const double& weightFinal)
{
  // create first for both
  m_h_cutflow->GetXaxis()->FindBin(label.c_str());
  m_h_cutflow_weighted->GetXaxis()->FindBin(label.c_str());
  m_h_cutflow_weighted_final->GetXaxis()->FindBin(label.c_str());
  
  const int bin = m_h_cutflow->GetXaxis()->FindBin(label.c_str());
  const double xx = m_h_cutflow->GetXaxis()->GetBinCenter(bin);
  
  m_h_cutflow->Fill( xx, 1 );
  m_h_cutflow_weighted->Fill( xx, xAHWeight );
  m_h_cutflow_weighted_final->Fill( xx, weightFinal );
}

EL::StatusCode ZJetBalanceMiniTreeAnaBase::LoadBasicConfiguration() 
{
  Info("LoadBasicConfiguration()", "called");
  m_configName = gSystem->ExpandPathName( m_configName.c_str() );
  Info("LoadBasicConfiguration()", "Configuing ZJetBalanceMiniTreeAnaBase Interface. User configuration read from : %s \n", m_configName.c_str());
  TEnv* config = new TEnv(m_configName.c_str());
  if( !config ) {
    Error("configure()", "Failed to read config file!");
    Error("configure()", "config name : %s",m_configName.c_str());
    return EL::StatusCode::FAILURE;
  }
  
  //
  // Read Input from .config file
  //
  m_debug                    = config->GetValue("Debug" ,      false );
  m_doPUreweighting          = config->GetValue("DoPileupReweighting", false);
  m_lumiCalcFileNames        = config->GetValue("LumiCalcFiles", "");
  m_PRWFileNames             = config->GetValue("PRWFiles", "");
  m_additional_weight        = config->GetValue("AdditionalWeight", 1.0);
  m_btagJets                 = config->GetValue("BTagJets", false);
  m_btagOP                   = config->GetValue("BTagOP", "70"); // 85, 77, 70, 60
  m_histPrefix               = config->GetValue("HistPrefix" ,      "" );
  
  if( m_btagJets ) {
    if( m_btagOP != "85" && m_btagOP != "77" && m_btagOP != "70" && m_btagOP != "60" ) {
      std::cout << "Invalid b-tag operating point " << m_btagOP << std::endl;
      return EL::StatusCode::FAILURE;
    }
    std::cout << "BTagging using the " << m_btagOP << "% operating point" << std::endl;
  }
  
  std::cout << "The basic configuration has been done in ZJetBalanceMiniTreeAnaBase :: LoadBasicConfiguration()" << std::endl;
  std::cout << "   - m_debug             " << ((m_debug) ? "True" : "False") << std::endl;
  std::cout << "   - m_doPUreweighting   " << ((m_doPUreweighting) ? "True" : "False") << std::endl;
  std::cout << "   - m_lumiCalcFileNames " << m_lumiCalcFileNames << std::endl;
  std::cout << "   - m_PRWFileNames      " << m_PRWFileNames << std::endl;
  std::cout << "   - m_additional_weight " << m_additional_weight << std::endl;
  std::cout << "   - m_btagJets          " << m_btagJets << std::endl;
  std::cout << "   - m_btagOP            " << m_btagOP << std::endl;
  std::cout << "   - m_histPrefix        " << m_histPrefix << std::endl;
  
  Info("LoadBasicConfiguration()", "ZJetBalanceMiniTreeAnaBase Interface succesfully configured! \n");
  return EL::StatusCode::SUCCESS;
}

void ZJetBalanceMiniTreeAnaBase :: SetBTagAddresses(TTree* tree)
{
  if(!m_btagJets) { return; }
  std::string btagAddress("jet_MV2c20_is"+m_btagOP);
  std::string btagSFAddress("jet_MV2c20_SF"+m_btagOP);
  tree->SetBranchAddress(btagAddress.c_str(),   &jet_isBTag, &b_jet_isBTag);
  tree->SetBranchAddress(btagSFAddress.c_str(), &jet_SFBTag, &b_jet_SFBTag);
  if( !(jet_isBTag && jet_SFBTag) ) {
    std::cout << "CANNOT Set B-Tagging Branches : " << std::endl;
    std::cout << "\t" << btagAddress.c_str() << std::endl;
    std::cout << "\t" << btagSFAddress.c_str() << std::endl;
  }
  // FIXME - induce an exit
}

void ZJetBalanceMiniTreeAnaBase :: InitTree(TTree* tree)
{
  // Set object pointer
  weight_muon_trig = 0;
  jet_E = 0;
  jet_pt = 0;
  jet_phi = 0;
  jet_eta = 0;
  jet_rapidity = 0;
  jet_HECFrac = 0;
  jet_EMFrac = 0;
  jet_CentroidR = 0;
  jet_FracSamplingMax = 0;
  jet_FracSamplingMaxIndex = 0;
  jet_LowEtConstituentsFrac = 0;
  jet_GhostMuonSegmentCount = 0;
  jet_Width = 0;
  jet_NumTrkPt1000PV = 0;
  jet_SumPtTrkPt1000PV = 0;
  jet_TrackWidthPt1000PV = 0;
  jet_NumTrkPt500PV = 0;
  jet_SumPtTrkPt500PV = 0;
  jet_TrackWidthPt500PV = 0;
  jet_JVFPV = 0;
  jet_Jvt = 0;
  jet_JvtJvfcorr = 0;
  jet_JvtRpt = 0;
  jet_SV0 = 0;
  jet_SV1 = 0;
  jet_IP3D = 0;
  jet_SV1IP3D = 0;
  jet_MV1 = 0;
  jet_MV2c00 = 0;
  jet_MV2c20 = 0;
  jet_MV2c20_is85 = 0;
  jet_MV2c20_SF85 = 0;
  jet_MV2c20_is77 = 0;
  jet_MV2c20_SF77 = 0;
  jet_MV2c20_is70 = 0;
  jet_MV2c20_SF70 = 0;
  jet_MV2c20_is60 = 0;
  jet_MV2c20_SF60 = 0;
  jet_isBTag = 0;
  jet_SFBTag = 0;
  jet_GhostArea = 0;
  jet_ActiveArea = 0;
  jet_VoronoiArea = 0;
  jet_ActiveArea4vec_pt = 0;
  jet_ActiveArea4vec_eta = 0;
  jet_ActiveArea4vec_phi = 0;
  jet_ActiveArea4vec_m = 0;
  jet_ConeTruthLabelID = 0;
  jet_TruthCount = 0;
  jet_TruthLabelDeltaR_B = 0;
  jet_TruthLabelDeltaR_C = 0;
  jet_TruthLabelDeltaR_T = 0;
  jet_PartonTruthLabelID = 0;
  jet_GhostTruthAssociationFraction = 0;
  jet_truth_E = 0;
  jet_truth_pt = 0;
  jet_truth_phi = 0;
  jet_truth_eta = 0;
  jet_constitScaleEta = 0;
  jet_emScaleEta = 0;
  jet_mucorrected_pt = 0;
  jet_mucorrected_eta = 0;
  jet_mucorrected_phi = 0;
  jet_mucorrected_m = 0;
  muon_pt = 0;
  muon_eta = 0;
  muon_phi = 0;
  muon_m = 0;
  muon_effSF = 0;
   
  tree->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
  tree->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
  tree->SetBranchAddress("mcEventNumber", &mcEventNumber, &b_mcEventNumber);
  tree->SetBranchAddress("mcChannelNumber", &mcChannelNumber, &b_mcChannelNumber);
  tree->SetBranchAddress("mcEventWeight", &mcEventWeight, &b_mcEventWeight);
  tree->SetBranchAddress("weight_pileup", &weight_pileup, &b_weight_pileup);
  tree->SetBranchAddress("NPV", &NPV, &b_NPV);
  tree->SetBranchAddress("actualInteractionsPerCrossing", &actualInteractionsPerCrossing, &b_actualInteractionsPerCrossing);
  tree->SetBranchAddress("averageInteractionsPerCrossing", &averageInteractionsPerCrossing, &b_averageInteractionsPerCrossing);
  tree->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
  tree->SetBranchAddress("rhoEM", &rhoEM, &b_rhoEM);
  tree->SetBranchAddress("pdgId1", &pdgId1, &b_pdgId1);
  tree->SetBranchAddress("pdgId2", &pdgId2, &b_pdgId2);
  tree->SetBranchAddress("pdfId1", &pdfId1, &b_pdfId1);
  tree->SetBranchAddress("pdfId2", &pdfId2, &b_pdfId2);
  tree->SetBranchAddress("x1", &x1, &b_x1);
  tree->SetBranchAddress("x2", &x2, &b_x2);
  tree->SetBranchAddress("xf1", &xf1, &b_xf1);
  tree->SetBranchAddress("xf2", &xf2, &b_xf2);
  tree->SetBranchAddress("ZpT", &ZpT, &b_ZpT);
  tree->SetBranchAddress("Zeta", &Zeta, &b_Zeta);
  tree->SetBranchAddress("Zphi", &Zphi, &b_Zphi);
  tree->SetBranchAddress("ZM", &ZM, &b_ZM);
  tree->SetBranchAddress("dPhiZJet1", &dPhiZJet1, &b_dPhiZJet1);
  tree->SetBranchAddress("dEtaZJet1", &dEtaZJet1, &b_dEtaZJet1);
  tree->SetBranchAddress("pTRef1", &pTRef1, &b_pTRef1);
  tree->SetBranchAddress("dPhiZJet2", &dPhiZJet2, &b_dPhiZJet2);
  tree->SetBranchAddress("dEtaZJet2", &dEtaZJet2, &b_dEtaZJet2);
  tree->SetBranchAddress("pTRef2", &pTRef2, &b_pTRef2);
  tree->SetBranchAddress("jetDPhi", &jetDPhi, &b_jetDPhi);
  tree->SetBranchAddress("jetDEta", &jetDEta, &b_jetDEta);
  tree->SetBranchAddress("jetPtRatio", &jetPtRatio, &b_jetPtRatio);
  tree->SetBranchAddress("weight", &weight, &b_weight);
  tree->SetBranchAddress("weight_xs", &weight_xs, &b_weight_xs);
  tree->SetBranchAddress("weight_prescale", &weight_prescale, &b_weight_prescale);
  tree->SetBranchAddress("weight_muon_trig", &weight_muon_trig, &b_weight_muon_trig);
  tree->SetBranchAddress("njets", &njets, &b_njets);
  tree->SetBranchAddress("jet_E", &jet_E, &b_jet_E);
  tree->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
  tree->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
  tree->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
  tree->SetBranchAddress("jet_rapidity", &jet_rapidity, &b_jet_rapidity);
  tree->SetBranchAddress("jet_HECFrac", &jet_HECFrac, &b_jet_HECFrac);
  tree->SetBranchAddress("jet_EMFrac", &jet_EMFrac, &b_jet_EMFrac);
  tree->SetBranchAddress("jet_CentroidR", &jet_CentroidR, &b_jet_CentroidR);
  tree->SetBranchAddress("jet_FracSamplingMax", &jet_FracSamplingMax, &b_jet_FracSamplingMax);
  tree->SetBranchAddress("jet_FracSamplingMaxIndex", &jet_FracSamplingMaxIndex, &b_jet_FracSamplingMaxIndex);
  tree->SetBranchAddress("jet_LowEtConstituentsFrac", &jet_LowEtConstituentsFrac, &b_jet_LowEtConstituentsFrac);
  tree->SetBranchAddress("jet_GhostMuonSegmentCount", &jet_GhostMuonSegmentCount, &b_jet_GhostMuonSegmentCount);
  tree->SetBranchAddress("jet_Width", &jet_Width, &b_jet_Width);
  tree->SetBranchAddress("jet_NumTrkPt1000PV", &jet_NumTrkPt1000PV, &b_jet_NumTrkPt1000PV);
  tree->SetBranchAddress("jet_SumPtTrkPt1000PV", &jet_SumPtTrkPt1000PV, &b_jet_SumPtTrkPt1000PV);
  tree->SetBranchAddress("jet_TrackWidthPt1000PV", &jet_TrackWidthPt1000PV, &b_jet_TrackWidthPt1000PV);
  tree->SetBranchAddress("jet_NumTrkPt500PV", &jet_NumTrkPt500PV, &b_jet_NumTrkPt500PV);
  tree->SetBranchAddress("jet_SumPtTrkPt500PV", &jet_SumPtTrkPt500PV, &b_jet_SumPtTrkPt500PV);
  tree->SetBranchAddress("jet_TrackWidthPt500PV", &jet_TrackWidthPt500PV, &b_jet_TrackWidthPt500PV);
  tree->SetBranchAddress("jet_JVFPV", &jet_JVFPV, &b_jet_JVFPV);
  tree->SetBranchAddress("jet_Jvt", &jet_Jvt, &b_jet_Jvt);
  tree->SetBranchAddress("jet_JvtJvfcorr", &jet_JvtJvfcorr, &b_jet_JvtJvfcorr);
  tree->SetBranchAddress("jet_JvtRpt", &jet_JvtRpt, &b_jet_JvtRpt);
  tree->SetBranchAddress("jet_SV0", &jet_SV0, &b_jet_SV0);
  tree->SetBranchAddress("jet_SV1", &jet_SV1, &b_jet_SV1);
  tree->SetBranchAddress("jet_IP3D", &jet_IP3D, &b_jet_IP3D);
  tree->SetBranchAddress("jet_SV1IP3D", &jet_SV1IP3D, &b_jet_SV1IP3D);
  tree->SetBranchAddress("jet_MV1", &jet_MV1, &b_jet_MV1);
  tree->SetBranchAddress("jet_MV2c00", &jet_MV2c00, &b_jet_MV2c00);
  tree->SetBranchAddress("jet_MV2c20", &jet_MV2c20, &b_jet_MV2c20);
  tree->SetBranchAddress("jet_MV2c20_is85", &jet_MV2c20_is85, &b_jet_MV2c20_is85);
  tree->SetBranchAddress("jet_MV2c20_SF85", &jet_MV2c20_SF85, &b_jet_MV2c20_SF85);
  tree->SetBranchAddress("jet_MV2c20_is77", &jet_MV2c20_is77, &b_jet_MV2c20_is77);
  tree->SetBranchAddress("jet_MV2c20_SF77", &jet_MV2c20_SF77, &b_jet_MV2c20_SF77);
  tree->SetBranchAddress("jet_MV2c20_is70", &jet_MV2c20_is70, &b_jet_MV2c20_is70);
  tree->SetBranchAddress("jet_MV2c20_SF70", &jet_MV2c20_SF70, &b_jet_MV2c20_SF70);
  tree->SetBranchAddress("jet_MV2c20_is60", &jet_MV2c20_is60, &b_jet_MV2c20_is60);
  tree->SetBranchAddress("jet_MV2c20_SF60", &jet_MV2c20_SF60, &b_jet_MV2c20_SF60);
  if(m_btagJets && !m_btagOP.empty() ) { this->SetBTagAddresses(tree); }
  tree->SetBranchAddress("jet_GhostArea", &jet_GhostArea, &b_jet_GhostArea);
  tree->SetBranchAddress("jet_ActiveArea", &jet_ActiveArea, &b_jet_ActiveArea);
  tree->SetBranchAddress("jet_VoronoiArea", &jet_VoronoiArea, &b_jet_VoronoiArea);
  tree->SetBranchAddress("jet_ActiveArea4vec_pt", &jet_ActiveArea4vec_pt, &b_jet_ActiveArea4vec_pt);
  tree->SetBranchAddress("jet_ActiveArea4vec_eta", &jet_ActiveArea4vec_eta, &b_jet_ActiveArea4vec_eta);
  tree->SetBranchAddress("jet_ActiveArea4vec_phi", &jet_ActiveArea4vec_phi, &b_jet_ActiveArea4vec_phi);
  tree->SetBranchAddress("jet_ActiveArea4vec_m", &jet_ActiveArea4vec_m, &b_jet_ActiveArea4vec_m);
  tree->SetBranchAddress("jet_ConeTruthLabelID", &jet_ConeTruthLabelID, &b_jet_ConeTruthLabelID);
  tree->SetBranchAddress("jet_TruthCount", &jet_TruthCount, &b_jet_TruthCount);
  tree->SetBranchAddress("jet_TruthLabelDeltaR_B", &jet_TruthLabelDeltaR_B, &b_jet_TruthLabelDeltaR_B);
  tree->SetBranchAddress("jet_TruthLabelDeltaR_C", &jet_TruthLabelDeltaR_C, &b_jet_TruthLabelDeltaR_C);
  tree->SetBranchAddress("jet_TruthLabelDeltaR_T", &jet_TruthLabelDeltaR_T, &b_jet_TruthLabelDeltaR_T);
  tree->SetBranchAddress("jet_PartonTruthLabelID", &jet_PartonTruthLabelID, &b_jet_PartonTruthLabelID);
  tree->SetBranchAddress("jet_GhostTruthAssociationFraction", &jet_GhostTruthAssociationFraction, &b_jet_GhostTruthAssociationFraction);
  tree->SetBranchAddress("jet_truth_E", &jet_truth_E, &b_jet_truth_E);
  tree->SetBranchAddress("jet_truth_pt", &jet_truth_pt, &b_jet_truth_pt);
  tree->SetBranchAddress("jet_truth_phi", &jet_truth_phi, &b_jet_truth_phi);
  tree->SetBranchAddress("jet_truth_eta", &jet_truth_eta, &b_jet_truth_eta);
  tree->SetBranchAddress("jet_constitScaleEta", &jet_constitScaleEta, &b_jet_constitScaleEta);
  tree->SetBranchAddress("jet_emScaleEta", &jet_emScaleEta, &b_jet_emScaleEta);
  tree->SetBranchAddress("jet_mucorrected_pt", &jet_mucorrected_pt, &b_jet_mucorrected_pt);
  tree->SetBranchAddress("jet_mucorrected_eta", &jet_mucorrected_eta, &b_jet_mucorrected_eta);
  tree->SetBranchAddress("jet_mucorrected_phi", &jet_mucorrected_phi, &b_jet_mucorrected_phi);
  tree->SetBranchAddress("jet_mucorrected_m", &jet_mucorrected_m, &b_jet_mucorrected_m);
  tree->SetBranchAddress("nmuon", &nmuon, &b_nmuon);
  tree->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
  tree->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
  tree->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
  tree->SetBranchAddress("muon_m", &muon_m, &b_muon_m);
  tree->SetBranchAddress("muon_effSF", &muon_effSF, &b_muon_effSF);
}
