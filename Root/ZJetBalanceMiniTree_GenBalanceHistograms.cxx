#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include <AthContainers/ConstDataVector.h>

#include <xAODTracking/VertexContainer.h>
#include <xAODJet/JetContainer.h>
#include <xAODEventInfo/EventInfo.h>
#include <ZJetBalance/ZJetBalanceMiniTree_GenBalanceHistograms.h>
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
ClassImp(ZJetBalanceMiniTree_GenBalanceHistograms)

ZJetBalanceMiniTree_GenBalanceHistograms :: ZJetBalanceMiniTree_GenBalanceHistograms ()
{
  Info("ZJetBalanceMiniTree_GenBalanceHistograms()", "Calling constructor");
}

EL::StatusCode  ZJetBalanceMiniTree_GenBalanceHistograms :: configure ()
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
  Info("configure()", "Configuing ZJetBalanceMiniTree_GenBalanceHistograms Interface. User configuration read from : %s \n", m_configName.c_str());
  TEnv* config = new TEnv(m_configName.c_str());
  if( !config ) {
    Error("configure()", "Failed to read config file!");
    Error("configure()", "config name : %s",m_configName.c_str());
    return EL::StatusCode::FAILURE;
  }
  
  TString pT_binning_str     = config->GetValue("pT_binning", "20,25,30,35,45,60,80,110,160,210,260,310,500");
  TString eta_binning_str    = config->GetValue("eta_binning", "-4.5,-3.2,-2.5,-1.0,0,1.0,2.5,3.2,4.5");
  m_nBinsXForResponseHist    = config->GetValue("nBinsXForResponseHist", 50.);
  m_maxXForResponseHist      = config->GetValue("maxXForResponseHist", 5.0);
  m_minXForResponseHist      = config->GetValue("minXForResponseHist", 0.0);
  m_cutDPhiZJet              = config->GetValue("cutDPhiZJet", 2.8);
  m_ZMassWindow              = config->GetValue("ZMassWindow", 15);
  m_fillLeptonBefore         = config->GetValue("FillLeptonBefore", false);
  m_debug                    = config->GetValue("Debug", false);  

  // create object before configuration
  m_pT_binning = new Double_t[BUFSIZ];
  m_eta_binning = new Double_t[BUFSIZ];
  
  
  DecodeBinning(pT_binning_str, m_pT_binning, m_n_pT_binning);
  Info("configure()", "DecodeBinning() gives for pT  : nBins=%d, first=%.1f last=%.1f", m_n_pT_binning, m_pT_binning[0], m_pT_binning[m_n_pT_binning]);
  DecodeBinning(eta_binning_str, m_eta_binning, m_n_eta_binning);
  Info("configure()", "DecodeBinning() gives for eta : nBins=%d, first=%.1f last=%.1f", m_n_eta_binning, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  
  config->Print();
  Info("configure()", "ZJetBalanceMiniTree_GenBalanceHistograms Interface succesfully configured! \n");
  
  Info("configure()", "Ends");
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ZJetBalanceMiniTree_GenBalanceHistograms :: setupJob (EL::Job& job)
{
  Info("setupJob()", "called");
  
  Info("setupJob()", "ends");
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTree_GenBalanceHistograms :: histInitialize ()
{
  Info("histInitialize()", "called");
  
  // configuration from ENV file (note : histInistialize called before initialize)
  if ( this->configure() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  // list of the histograms
  TH1::SetDefaultSumw2();
  
  this->InitializeCutflows();
  
  // add your own histgrams
  m_h_RunNumber = book("RunNumber", "Run Number", 1000, 271000.5, 272000.5);
  m_h_muon1_pT = book("muon1_pT", "#mu_{1} p_{T} [GeV]", 120, 0, 120);
  m_h_muon2_pT = book("muon2_pT", "#mu_{2} p_{T} [GeV]", 120, 0, 120);
  m_h_muon1_eta = book("muon1_eta", "#mu_{1} #eta", 60, -3.0, 3.0);
  m_h_muon2_eta = book("muon2_eta", "#mu_{2} #eta", 60, -3.0, 3.0);
  m_h_muon1_phi = book("muon1_phi", "#mu_{1} #phi", 64, -TMath::Pi(), TMath::Pi());
  m_h_muon2_phi = book("muon2_phi", "#mu_{2} #phi", 64, -TMath::Pi(), TMath::Pi());

  m_h_electron1_pT = book("electron1_pT", "e_{1} p_{T} [GeV]", 120, 0, 120);
  m_h_electron2_pT = book("electron2_pT", "e_{2} p_{T} [GeV]", 120, 0, 120);
  m_h_electron1_eta = book("electron1_eta", "e_{1} #eta", 60, -3.0, 3.0);
  m_h_electron2_eta = book("electron2_eta", "e_{2} #eta", 60, -3.0, 3.0);
  m_h_electron1_phi = book("electron1_phi", "e_{1} #phi", 64, -TMath::Pi(), TMath::Pi());
  m_h_electron2_phi = book("electron2_phi", "e_{2} #phi", 64, -TMath::Pi(), TMath::Pi());

  // plots of Z itself
  m_h_ZpT    = book("ZpT",   "Z p_{T} [GeV]",  120,  0,    240);
  m_h_ZpTRef = book("ZpTRef","Z p_{T}^{ref} [GeV]",  120,  0,   240);
  m_h_Zeta   = book("Zeta",  "Z #eta",         60,  -3.0, 3.0);
  m_h_Zphi   = book("Zphi",  "Z #phi",         64, -TMath::Pi(), TMath::Pi());
  m_h_ZM     = book("ZM",    "m_{Z} [GeV]",    60, 60, 120);

  m_h_ZpT_beforecut   = book("ZpT_beforecut",   "Z p_{T} [GeV]",  120,  0,    240);
  m_h_Zeta_beforecut  = book("Zeta_beforecut",  "Z #eta",         60,  -3.0, 3.0);
  m_h_Zphi_beforecut  = book("Zphi_beforecut",  "Z #phi",         64, -TMath::Pi(), TMath::Pi());
  m_h_ZM_beforecut    = book("ZM_beforecut",    "m_{Z} [GeV]",    60, 60, 120);

  // Jet Plots
  // before cuts
  m_h_nJets_beforecut = book("nJets_beforecut", "N Jets Before Cuts", 10, -0.5, 9.5);
  m_h_jet_eta_beforecut = book("jet_eta_beforecut", "",   64, -3.2, 3.2);
  m_h_jet_pt_beforecut  = book("jet_pt_beforecut", "",    30,    0, 300);
  m_h_jet_eta_beforecut_b = book("jet_eta_beforecut_b", "",   64, -3.2, 3.2);
  m_h_jet_pt_beforecut_b  = book("jet_pt_beforecut_b", "",    30,    0, 300);
  m_h_jet_eta_beforecut_c = book("jet_eta_beforecut_c", "",   64, -3.2, 3.2);
  m_h_jet_pt_beforecut_c  = book("jet_pt_beforecut_c", "",    30,    0, 300);
  m_h_jet_eta_beforecut_l = book("jet_eta_beforecut_l", "",   64, -3.2, 3.2);
  m_h_jet_pt_beforecut_l  = book("jet_pt_beforecut_l", "",    30,    0, 300);
  
  // leading Jet Plots
  // before cuts  
  m_h_1st_jet_eta_beforecut = book("1st_jet_eta_beforecut", "",   64, -3.2, 3.2);
  m_h_1st_jet_pt_beforecut  = book("1st_jet_pt_beforecut", "",    30,    0, 300);
  m_h_1st_jet_eta_beforecut_b = book("1st_jet_eta_beforecut_b", "",   64, -3.2, 3.2);
  m_h_1st_jet_pt_beforecut_b  = book("1st_jet_pt_beforecut_b", "",    30,    0, 300);
  m_h_1st_jet_eta_beforecut_c = book("1st_jet_eta_beforecut_c", "",   64, -3.2, 3.2);
  m_h_1st_jet_pt_beforecut_c  = book("1st_jet_pt_beforecut_c", "",    30,    0, 300);
  m_h_1st_jet_eta_beforecut_l = book("1st_jet_eta_beforecut_l", "",   64, -3.2, 3.2);
  m_h_1st_jet_pt_beforecut_l  = book("1st_jet_pt_beforecut_l", "",    30,    0, 300);

  // after cuts
  m_h_nJets = book("nJets",           "N Jets", 10, -0.5, 9.5);
  m_h_jet_phi   = book("jet_phi", "jet_{1} #phi",   64, -TMath::Pi(), TMath::Pi());
  m_h_jet_phi_b = book("jet_phi_b", "jet_{1} #phi MC truth b",   64, -TMath::Pi(), TMath::Pi());
  m_h_jet_phi_c = book("jet_phi_c", "jet_{1} #phi MC truth c",   64, -TMath::Pi(), TMath::Pi());
  m_h_jet_phi_l = book("jet_phi_l", "jet_{1} #phi MC truth l",   64, -TMath::Pi(), TMath::Pi());
  m_h_jet_eta = book("jet_eta", "jet^{1} #eta", 50, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  m_h_jet_eta_b = book("jet_eta_b", "jet^{1} #eta MC truth b", 50, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  m_h_jet_eta_c = book("jet_eta_c", "jet^{1} #eta MC truth c", 50, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  m_h_jet_eta_l = book("jet_eta_l", "jet^{1} p_{T} MC truth l", 50, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  m_h_jet_pt  = book("jet_pt",  "jet^{1} p_{T}", 50, 0, 300.);
  m_h_jet_pt_b  = book("jet_pt_b",  "jet^{1} p_{T} MC truth b", 50, 0, 300.);
  m_h_jet_pt_c  = book("jet_pt_c",  "jet^{1} p_{T} MC truth c", 50, 0, 300.);
  m_h_jet_pt_l  = book("jet_pt_l",  "jet^{1} #eta MC truth l", 50, 0, 300.);
  m_h_SumPtTrkPt500PV_beforecut = book("track_sum_pt_500_MeV_beforecut","", 50, 0, 300.);
  m_h_SumPtTrkPt500PV = book("track_sum_pt_500_MeV","", 50, 0, 300.);

  m_h_jet_pt_bin = new TH2F("jet_pt_bin", "", 100, 0, 500, m_n_pT_binning, -0.5, -0.5+m_n_pT_binning); // validation purpose
  wk()->addOutput( m_h_jet_pt_bin );
  
  m_h_pt_binning_info = book("pt_binning_info", "", m_n_pT_binning, m_pT_binning);
  m_h_eta_binning_info = book("eta_binning_info", "", m_n_eta_binning, m_eta_binning);
  
  // Z-Jet relationship
  m_h_Z_jet_dPhi = book("Z_jet_dPhi", "", 100, -TMath::Pi(), TMath::Pi());
  m_h_Z_jet_dEta = book("Z_jet_dEta", "", 64, -3.2, 3.2);
  // balance hist instanted elsewhere as uses custom/configurable binning
  // weights
  m_h_prwfactor  = book("prwfactor", "Pile-up Weight", 100, 0, 3.0);
  m_h_muonTrigFactor =  book("muonTrigFactor", "Muon Trigger Weight", 100, 0.0, 2.0);
  m_h_muon1EffFactor =  book("muon1EffFactor", "muon_{1} Efficiency SF", 100, 0.0, 2.0);
  m_h_muon2EffFactor =  book("muon2EffFactor", "muon_{2} Efficiency SF", 100, 0.0, 2.0);    

  m_h_electronTrigFactor =  book("electronTrigFactor", "Electron Trigger Weight", 100, 0.0, 2.0);
  m_h_electron1EffFactor =  book("electron1EffFactor", "electron_{1} Efficiency SF (Reco*PID)", 100, 0.0, 2.0);
  m_h_electron2EffFactor =  book("electron2EffFactor", "electron_{2} Efficiency SF (Reco*PID)", 100, 0.0, 2.0);
  
  m_h_averageInteractionsPerCrossing = book("averageInteractionsPerCrossing", "", 50, 0, 50.);
  
  // before cuts muon plots
  if( m_fillLeptonBefore ) {
    m_h_muon1_pT_beforecut = book("muon1_pT_beforecut", "#mu_{1} p_{T} [GeV]", 120, 0, 120);
    m_h_muon2_pT_beforecut = book("muon2_pT_beforecut", "#mu_{2} p_{T} [GeV]", 120, 0, 120);
    m_h_muon1_eta_beforecut = book("muon1_eta_beforecut", "#mu_{1} #eta", 60, -3.0, 3.0);
    m_h_muon2_eta_beforecut = book("muon2_eta_beforecut", "#mu_{2} #eta", 60, -3.0, 3.0);
    m_h_muon1_phi_beforecut = book("muon1_phi_beforecut", "#mu_{1} #phi", 64, -TMath::Pi(), TMath::Pi());
    m_h_muon2_phi_beforecut = book("muon2_phi_beforecut", "#mu_{2} #phi", 64, -TMath::Pi(), TMath::Pi());

    m_h_electron1_pT_beforecut = book("electron1_pT_beforecut", "e_{1} p_{T} [GeV]", 120, 0, 120);
    m_h_electron2_pT_beforecut = book("electron2_pT_beforecut", "e_{2} p_{T} [GeV]", 120, 0, 120);
    m_h_electron1_eta_beforecut = book("electron1_eta_beforecut", "e_{1} #eta", 60, -3.0, 3.0);
    m_h_electron2_eta_beforecut = book("electron2_eta_beforecut", "e_{2} #eta", 60, -3.0, 3.0);
    m_h_electron1_phi_beforecut = book("electron1_phi_beforecut", "e_{1} #phi", 64, -TMath::Pi(), TMath::Pi());
    m_h_electron2_phi_beforecut = book("electron2_phi_beforecut", "e_{2} #phi", 64, -TMath::Pi(), TMath::Pi());

  }
  
  // additinoal histograms according to configuration parameters
  for (int iPtBin=1; iPtBin<m_n_pT_binning+1; iPtBin++) {
    std::vector<TH1F*> tmp_hist_container;
    std::vector<TH1F*> tmp_hist_container_b;
    std::vector<TH1F*> tmp_hist_container_c;
    std::vector<TH1F*> tmp_hist_container_l;
    for (int iEtaBin=1; iEtaBin<m_n_eta_binning+1; iEtaBin++) {
      Info("Initialize()", "%s", Form("DB_RefEtaBin%d_PtBin%d", iEtaBin, iPtBin));
      TH1F* h = new TH1F(Form("DB_RefEtaBin%d_PtBin%d", iEtaBin, iPtBin), 
			 Form("%.1f<p_{T}<%.1f, %.1f<#eta<%.1f", 
			      m_pT_binning[iPtBin-1], m_pT_binning[iPtBin],
			      m_eta_binning[iEtaBin-1], m_eta_binning[iEtaBin]
			      ), 
			 m_nBinsXForResponseHist, m_minXForResponseHist, m_maxXForResponseHist);
      TH1F* h_b = new TH1F(Form("DB_RefEtaBin%d_PtBin%d_b", iEtaBin, iPtBin), 
			   Form("%.1f<p_{T}<%.1f, %.1f<#eta<%.1f", 
				m_pT_binning[iPtBin-1], m_pT_binning[iPtBin],
				m_eta_binning[iEtaBin-1], m_eta_binning[iEtaBin]
				), 
			   m_nBinsXForResponseHist, m_minXForResponseHist, m_maxXForResponseHist);
      TH1F* h_c = new TH1F(Form("DB_RefEtaBin%d_PtBin%d_c", iEtaBin, iPtBin), 
			   Form("%.1f<p_{T}<%.1f, %.1f<#eta<%.1f", 
				m_pT_binning[iPtBin-1], m_pT_binning[iPtBin],
				m_eta_binning[iEtaBin-1], m_eta_binning[iEtaBin]
			      ), 
			   m_nBinsXForResponseHist, m_minXForResponseHist, m_maxXForResponseHist);
      TH1F* h_l = new TH1F(Form("DB_RefEtaBin%d_PtBin%d_l", iEtaBin, iPtBin), 
			   Form("%.1f<p_{T}<%.1f, %.1f<#eta<%.1f", 
				m_pT_binning[iPtBin-1], m_pT_binning[iPtBin],
				m_eta_binning[iEtaBin-1], m_eta_binning[iEtaBin]
				), 
			   m_nBinsXForResponseHist, m_minXForResponseHist, m_maxXForResponseHist);
      tmp_hist_container.push_back(h);
      tmp_hist_container_b.push_back(h_b);
      tmp_hist_container_c.push_back(h_c);
      tmp_hist_container_l.push_back(h_l);
      wk()->addOutput( h );
      wk()->addOutput( h_b );
      wk()->addOutput( h_c );
      wk()->addOutput( h_l );
    }
    
    m_balance_hists.push_back(tmp_hist_container);
    m_balance_hists_b.push_back(tmp_hist_container_b);
    m_balance_hists_c.push_back(tmp_hist_container_c);
    m_balance_hists_l.push_back(tmp_hist_container_l);
  }  
  
  // pass all histograms instanted with the book function are passed
  // to the worker through 
  this->record( wk() );

  Info("histInitialize()", "ends");
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTree_GenBalanceHistograms :: fileExecute ()
{
  // Once per file, determine whether the file is a muon sample or not
  m_isMuonSample = true;
  IsMuonSample();
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ZJetBalanceMiniTree_GenBalanceHistograms :: initialize ()
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
  
  Info("initialize()", "Succesfully initialized! \n");
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ZJetBalanceMiniTree_GenBalanceHistograms :: execute ()
{
  this->LoadMiniTree();

  double weight_final=1.0;
  if (m_isMC) {
    if(m_debug) { 
      Info("execute()","Line %d: Getting PRW factor", __LINE__); 
    }

    double pileup_reweighting_factor = GetPileupReweightingFactor();
    if(m_debug) { 
      Info("execute()","Line %d: PRW factor successfully obtained!", __LINE__); 
    }
    weight_final = 
      mcEventWeight*weight_xs*pileup_reweighting_factor*m_additional_weight;
    
    if (m_isMuonSample){
      // trigger weight: 0 is nominal, rest are systematics +,- 1 sigma for each
      m_h_muonTrigFactor->Fill( weight_muon_trig->at(0) );
      if( weight_muon_trig->at(0)!= -999.0 )
        weight_final *= weight_muon_trig->at(0);
      // muon efficiency scale factors: 0 is nominal, rest are systematics +,- 1 sigma for each
      m_h_muon1EffFactor->Fill( muon_RecoEff_SF->at(0)[0] * muon_IsoEff_SF_Gradient->at(0)[0] );
      m_h_muon2EffFactor->Fill( muon_RecoEff_SF->at(1)[0] * muon_IsoEff_SF_Gradient->at(1)[0]);
      if( muon_RecoEff_SF->at(0)[0]!= -999.0 && muon_RecoEff_SF->at(1)[0]!= -999.0 )
        weight_final *= muon_RecoEff_SF->at(0)[0] * muon_RecoEff_SF->at(1)[0];
      if( muon_IsoEff_SF_Gradient->at(0)[0] != -999.0 && muon_IsoEff_SF_Gradient->at(1)[0] != -999.0)
        weight_final *= muon_IsoEff_SF_Gradient->at(0)[0] * muon_IsoEff_SF_Gradient->at(1)[0];
    } else {
      
      // trigger weight: 0 is nominal, rest are systematics +,- 1 sigma for each
      m_h_electronTrigFactor->Fill( weight_electron_trig->at(0) );
      if(weight_electron_trig->at(0) != -999.0 ) 
        weight_final *= weight_electron_trig->at(0);
      // muon efficiency scale factors: 0 is nominal, rest are systematics +,- 1 sigma for each
      
      m_h_electron1EffFactor->Fill( el_PIDEff_SF_LHMedium->at(0)[0]*el_RecoEff_SF->at(0)[0] );
      m_h_electron2EffFactor->Fill( el_PIDEff_SF_LHMedium->at(1)[0]*el_RecoEff_SF->at(1)[0] );
      
      if( el_PIDEff_SF_LHMedium->at(0)[0]!= -999.0  && el_PIDEff_SF_LHMedium->at(1)[0]!= -999.0)
        weight_final *= el_PIDEff_SF_LHMedium->at(0)[0] * el_PIDEff_SF_LHMedium->at(1)[0];
      if( el_RecoEff_SF->at(0)[0]!= -999.0 && el_RecoEff_SF->at(1)[0]!= -999.0 )
        weight_final *= el_RecoEff_SF->at(0)[0] * el_RecoEff_SF->at(1)[0] ;
      
    }
  }
  
  FillCutflowHistograms("Process NTuple", mcEventWeight, weight_final);
  m_h_RunNumber->Fill(runNumber, weight_final);
  m_h_averageInteractionsPerCrossing->Fill(averageInteractionsPerCrossing, weight_final); // for validation


  if(m_debug) Info("execute()", "Processing Event @ RunNumber=%10d, EventNumber=%d", runNumber, eventNumber);
  ++m_eventCounter;
  if (m_eventCounter%10000==0) {
    Info("execute()", "%10d th event is been processed.", m_eventCounter);
  }
  
  //
  if(m_debug) Info("execute()","Line %d: Accessing jet kinematic variables.",__LINE__); 
  const float& lead_jet_pt      = jet_pt->at(0);
  const float& lead_jet_eta     = jet_eta->at(0);
  const float& lead_jet_constScale_eta = jet_constitScaleEta->at(0);
  const float& lead_jet_phi     = jet_phi->at(0);
  const int    lead_jet_truthLabel = (m_isMC && (jet_ConeTruthLabelID)) ? jet_ConeTruthLabelID->at(0) : -1;
  const int    Z_pt_ref_bin  = GetPtBin(pTRef1);
  //const int    lead_jet_eta_bin = GetEtaBin(lead_jet_constScale_eta);
  const int    lead_jet_eta_bin = GetEtaBin(lead_jet_eta);

  if(m_debug) Info("execute()","Line %d: Accessing sum track pt.",__LINE__); 
  const float& sumPtTrkPt500PV = jet_SumPtTrkPt500PV->at(0);
  
  if (lead_jet_eta_bin==-1) {return EL::StatusCode::SUCCESS;} // out of eta range (defined as binning)
  FillCutflowHistograms("jet_{1} #eta", mcEventWeight, weight_final);
  
  // selection criteria need to be applied
  m_h_ZM_beforecut->Fill(ZM, weight_final);
  
  if (TMath::Abs(ZM-91)>m_ZMassWindow)      { return EL::StatusCode::SUCCESS; }
  FillCutflowHistograms("m_{Z} Window", mcEventWeight, weight_final);
  
  // for valiadtion
  int nJetsBeforeCut = 0;
  
  m_h_1st_jet_eta_beforecut->Fill(jet_eta->at(0), weight_final);
  m_h_1st_jet_pt_beforecut->Fill(jet_pt->at(0), weight_final);
  if (m_isMC) {    
    FillFlavorHistograms(m_h_1st_jet_pt_beforecut_b, m_h_1st_jet_pt_beforecut_c, m_h_1st_jet_pt_beforecut_l, 
			 lead_jet_truthLabel, jet_pt->at(0), weight_final);
    FillFlavorHistograms(m_h_1st_jet_eta_beforecut_b, m_h_1st_jet_eta_beforecut_c, m_h_1st_jet_eta_beforecut_l,
			 lead_jet_truthLabel, jet_eta->at(0), weight_final);
  }
  
  for (int iJet=0, nJets=jet_pt->size(); iJet<nJets; iJet++) {
    const float& pt  = jet_pt->at(iJet);
    const float& eta = jet_eta->at(iJet);
    const int& truthLabel = (m_isMC && (jet_ConeTruthLabelID)) ? jet_ConeTruthLabelID->at(iJet) : -1;
    
    if ( pt < 30. ) continue;
    nJetsBeforeCut++;
    
    m_h_jet_eta_beforecut->Fill(eta, weight_final);
    m_h_jet_pt_beforecut->Fill(pt, weight_final);
    FillFlavorHistograms(m_h_jet_pt_beforecut_b, m_h_jet_pt_beforecut_c, m_h_jet_pt_beforecut_l, 
			 truthLabel, pt, weight_final);
    FillFlavorHistograms(m_h_jet_eta_beforecut_b, m_h_jet_eta_beforecut_c, m_h_jet_eta_beforecut_l,
			 truthLabel, eta, weight_final);
  }
  m_h_nJets_beforecut->Fill(nJetsBeforeCut, weight_final);
  m_h_ZpT_beforecut->Fill(ZpT, weight_final);
  m_h_Zeta_beforecut->Fill(Zeta, weight_final);
  m_h_Zphi_beforecut->Fill(Zphi, weight_final);
  m_h_ZM_beforecut->Fill(ZM, weight_final);
  m_h_SumPtTrkPt500PV_beforecut->Fill(sumPtTrkPt500PV, weight_final);

  // muon before cut
  if( m_fillLeptonBefore && m_isMuonSample ) {
    m_h_muon1_pT_beforecut->Fill ( muon_pt ->at(0), weight_final );
    m_h_muon1_eta_beforecut->Fill( muon_eta->at(0), weight_final );
    m_h_muon1_phi_beforecut->Fill( muon_phi->at(0), weight_final );
    m_h_muon2_pT_beforecut->Fill ( muon_pt ->at(1), weight_final );
    m_h_muon2_eta_beforecut->Fill( muon_eta->at(1), weight_final );
    m_h_muon2_phi_beforecut->Fill( muon_phi->at(1), weight_final );
  } else if( m_fillLeptonBefore && (!m_isMuonSample) ){
    m_h_electron1_pT_beforecut->Fill ( el_pt ->at(0), weight_final );
    m_h_electron1_eta_beforecut->Fill( el_eta->at(0), weight_final );
    m_h_electron1_phi_beforecut->Fill( el_phi->at(0), weight_final );
    m_h_electron2_pT_beforecut->Fill ( el_pt ->at(1), weight_final );
    m_h_electron2_eta_beforecut->Fill( el_eta->at(1), weight_final );
    m_h_electron2_phi_beforecut->Fill( el_phi->at(1), weight_final );
  }
  
  if (TMath::Abs(dPhiZJet1)<m_cutDPhiZJet)  { return EL::StatusCode::SUCCESS; }
  FillCutflowHistograms("#Delta#phi(Z,jet)", mcEventWeight, weight_final);
  
  
  //Info("execute()", "weight=%.f m_isMC=%s ptbin=%d etabin=%d", weight_final, (m_isMC ? "True" : "False"), Z_pt_ref_bin, lead_jet_eta_bin);

  //Info("execute()", "lead_jet_eta=%.1f (%d) lead_jet_pt=%.1f (%d)",
  //lead_jet_eta, lead_jet_eta_bin, lead_jet_pt, Z_pt_ref_bin);
  
  if (jet_pt->size()>1) { if (jet_pt->at(1)>pTRef1*0.2) {return EL::StatusCode::SUCCESS;} } // event with second jet is vetoed
  FillCutflowHistograms("jet_{2} p_{T}", mcEventWeight, weight_final);

  if( m_btagJets ) {
    // b-tag if asked for
    if( !jet_isBTag->at(0) ) { return EL::StatusCode::SUCCESS; }
    FillCutflowHistograms("jet_{1} b-tagged", mcEventWeight, weight_final);
    // apply b-tagging weight
    if (m_isMC) weight_final *= jet_SFBTag->at(0)[0]; // other weights are for systematics
    //Info("execute()", "jet_SFBTag->at(0)[0]=%f \n", jet_SFBTag->at(0)[0]);
    
    FillCutflowHistograms("jet_{1} b-tag SF", mcEventWeight, weight_final);
  } else {
    if( jet_isBTag->at(0) ) { 
      // only in order to store the cutflow 
      FillCutflowHistograms("jet_{1} b-tagged", mcEventWeight, weight_final);
      FillCutflowHistograms("jet_{1} b-tag SF", mcEventWeight, ((m_isMC) ? weight_final*jet_SFBTag->at(0)[0] : weight_final));
    }
  }

  // lepton plots
  if( m_isMuonSample ) {
    m_h_muon1_pT->Fill ( muon_pt ->at(0), weight_final );
    m_h_muon1_eta->Fill( muon_eta->at(0), weight_final );
    m_h_muon1_phi->Fill( muon_phi->at(0), weight_final );
    m_h_muon2_pT->Fill ( muon_pt ->at(1), weight_final );
    m_h_muon2_eta->Fill( muon_eta->at(1), weight_final );
    m_h_muon2_phi->Fill( muon_phi->at(1), weight_final );
  } else{
    m_h_electron1_pT->Fill ( el_pt ->at(0), weight_final );
    m_h_electron1_eta->Fill( el_eta->at(0), weight_final );
    m_h_electron1_phi->Fill( el_phi->at(0), weight_final );
    m_h_electron2_pT->Fill ( el_pt ->at(1), weight_final );
    m_h_electron2_eta->Fill( el_eta->at(1), weight_final );
    m_h_electron2_phi->Fill( el_phi->at(1), weight_final );
  }
  
  // plots of Z itself
  m_h_ZpT->Fill(ZpT, weight_final);
  m_h_ZpTRef->Fill(pTRef1, weight_final);
  m_h_Zeta->Fill(Zeta, weight_final);
  m_h_Zphi->Fill(Zphi, weight_final);
  m_h_ZM->Fill(ZM, weight_final);
  m_h_SumPtTrkPt500PV->Fill(sumPtTrkPt500PV, weight_final);
  // jet plots after cuts
  m_h_nJets->Fill(njets, weight_final);
  m_h_jet_eta->Fill(lead_jet_eta, weight_final);
  m_h_jet_pt->Fill(lead_jet_pt, weight_final);
  m_h_jet_phi->Fill(lead_jet_phi, weight_final);
  FillFlavorHistograms(m_h_jet_pt_b, m_h_jet_pt_c, m_h_jet_pt_l, 
		       lead_jet_truthLabel, lead_jet_pt, weight_final);
  FillFlavorHistograms(m_h_jet_eta_b, m_h_jet_eta_c, m_h_jet_eta_l, 
		       lead_jet_truthLabel, lead_jet_eta, weight_final);
  FillFlavorHistograms(m_h_jet_phi_b, m_h_jet_phi_c, m_h_jet_phi_l, 
		       lead_jet_truthLabel, lead_jet_phi, weight_final);
  // Z-Jet relationship
  m_h_Z_jet_dPhi->Fill(dPhiZJet1, weight_final);
  m_h_Z_jet_dEta->Fill(dEtaZJet1, weight_final);
  
  if (Z_pt_ref_bin>=0) { // valid Z_pT bin
    m_h_jet_pt_bin->Fill(lead_jet_pt, Z_pt_ref_bin, weight_final);
    (m_balance_hists[Z_pt_ref_bin])[lead_jet_eta_bin]->Fill(lead_jet_pt/pTRef1, weight_final);
    FillFlavorHistograms((m_balance_hists_b[Z_pt_ref_bin])[lead_jet_eta_bin],
			 (m_balance_hists_c[Z_pt_ref_bin])[lead_jet_eta_bin],
			 (m_balance_hists_l[Z_pt_ref_bin])[lead_jet_eta_bin],
			 lead_jet_truthLabel, lead_jet_pt/pTRef1, weight_final);    
  }
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ZJetBalanceMiniTree_GenBalanceHistograms :: postExecute ()
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTree_GenBalanceHistograms :: finalize ()
{
  std::cout << "Finialize!" << std::endl;
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ZJetBalanceMiniTree_GenBalanceHistograms :: histFinalize ()
{
  return EL::StatusCode::SUCCESS;
}

void ZJetBalanceMiniTree_GenBalanceHistograms::DecodeBinning(TString binning_str, Double_t* binning_array, Int_t& n_binning)
{
  TString tok;
  Ssiz_t  from = 0;
  binning_str.ReplaceAll(" ", "");
  int iArrayIndex=0;
  while (binning_str.Tokenize(tok, from, ",")) {  
    const double x = strtod(tok.Data(), NULL);
    binning_array[iArrayIndex]=x;
    iArrayIndex++;
  }
  n_binning=iArrayIndex-1;
  
  return;
}

int ZJetBalanceMiniTree_GenBalanceHistograms::GetPtBin(const double& _pt)
{
  int rc=-1;
  for (int iBin=0; iBin<m_n_pT_binning; iBin++) {
    if (m_pT_binning[iBin]<_pt && _pt<m_pT_binning[iBin+1]) {
      rc=iBin;
      break;
    }
  }
  return rc;
}

int ZJetBalanceMiniTree_GenBalanceHistograms::GetEtaBin(const double& _eta)
{
  int rc=-1;
  for (int iBin=0; iBin<m_n_eta_binning; iBin++) {
    if (m_eta_binning[iBin]<_eta && _eta<m_eta_binning[iBin+1]) {
      rc = iBin;
      break;
    }
  }
  return rc;
}

void ZJetBalanceMiniTree_GenBalanceHistograms::FillFlavorHistograms(TH1F* h_b, TH1F* h_c, TH1F* h_l, const int& truthLabel, const float& value, const float& weight)
{
  switch (TMath::Abs(truthLabel)) {
  case 5:
    h_b->Fill(value, weight);
    break;
  case 4:
    h_c->Fill(value, weight);
    break;
  default:
    h_l->Fill(value, weight);
    break;
  }
}

void ZJetBalanceMiniTree_GenBalanceHistograms::IsMuonSample()
{
  Info("IsMuonSample()", "Checking to see if sample contains muons or electrons");
  TFile* inputFile = wk()->inputFile();
  m_isMuonSample = ((TTree*)inputFile->Get("outTree"))->GetBranch("muon_pt");
  if( true ){
    std::string info = m_isMuonSample?"Processing muon sample":"Processing electron sample";
    Info("IsMuonSample()", info.c_str());
  }
}
