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

#include <stdlib.h>

using namespace std;

// this is needed to distribute the algorithm to the workers
ClassImp(ProcessZJetBalanceMiniTree)

ProcessZJetBalanceMiniTree :: ProcessZJetBalanceMiniTree ()
{
  Info("ProcessZJetBalanceMiniTree()", "Calling constructor");
}


EL::StatusCode  ProcessZJetBalanceMiniTree :: configure ()
{
  Info("configure()", "called");
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
  m_debug                    = config->GetValue("Debug" ,      false );
  TString pT_binning_str     = config->GetValue("pT_binning", "20,25,30,35,45,60,80,110,160,210,260,310,500");
  TString eta_binning_str    = config->GetValue("eta_binning", "-4.5,-3.2,-2.5,-1.0,0,1.0,2.5,3.2,4.5");
  m_nBinsXForResponseHist    = config->GetValue("nBinsXForResponseHist", 50.);
  m_maxXForResponseHist      = config->GetValue("maxXForResponseHist", 5.0);
  m_minXForResponseHist      = config->GetValue("minXForResponseHist", 0.0);
  m_doPUreweighting          = config->GetValue("DoPileupReweighting", false);
  m_lumiCalcFileNames        = config->GetValue("LumiCalcFiles", "");
  m_PRWFileNames             = config->GetValue("PRWFiles", "");
  m_cutDPhiZJet              = config->GetValue("cutDPhiZJet", 2.8);
  m_ZMassWindow              = config->GetValue("ZMassWindow", 15);
  m_MV2c20threshold          = config->GetValue("MV2c20threshold", -0.5911); 
  
  DecodeBinning(pT_binning_str, m_pT_binning, m_n_pT_binning);
  Info("configure()", "DecodeBinning() gives for pT  : nBins=%d, first=%.1f last=%.1f", m_n_pT_binning, m_pT_binning[0], m_pT_binning[m_n_pT_binning]);
  DecodeBinning(eta_binning_str, m_eta_binning, m_n_eta_binning);
  Info("configure()", "DecodeBinning() gives for eta : nBins=%d, first=%.1f last=%.1f", m_n_eta_binning, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  
  config->Print();
  Info("configure()", "ProcessZJetBalanceMiniTree Interface succesfully configured! \n");
  
  Info("configure()", "Ends");
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ProcessZJetBalanceMiniTree :: setupJob (EL::Job& job)
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



EL::StatusCode ProcessZJetBalanceMiniTree :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  Info("histInitialize()", "called");
  
  // list of the histograms
  TH1::SetDefaultSumw2();
  m_h_RunNumber = new TH1D("h_RunNumber", "", 1000, 271000.5, 272000.5);
  m_h_ZpT   = new TH1D("h_ZpT", "", 100, 0, 100); // validation purpose
  m_h_ZM    = new TH1D("h_ZM",  "", 100, 71, 110); // validation purpose
  m_h_nJets = new TH1D("h_nJets", "", 10, -0.5, 9.5); // validation purpose
  m_h_Z_jet_dPhi = new TH1D("h_Z_jet_dPhi", "", 100, -TMath::Pi(), TMath::Pi());
  m_h_prwfactor  = new TH1D("h_prwfactor", "", 100, 0, 3.0);
  m_h_njets_beforecut   = new TH1D("h_njets_beforecut", "", 7, -0.5, 6.5);
  m_h_jet_eta_beforecut = new TH1D("h_jet_eta_beforecut", "", 32, -3.2, 3.2);
  m_h_jet_pt_beforecut  = new TH1D("h_jet_pt_beforecut", "", 30, 0, 300);  
  m_h_nbjets_beforecut   = new TH1D("h_nbjets_beforecut", "", 7, -0.5, 6.5);
  m_h_bjet_eta_beforecut = new TH1D("h_bjet_eta_beforecut", "", 32, -3.2, 3.2);
  m_h_bjet_pt_beforecut  = new TH1D("h_bjet_pt_beforecut", "", 30, 0, 300);
  m_h_bjet_eta_beforecut_b = new TH1D("h_bjet_eta_beforecut_b", "", 32, -3.2, 3.2);
  m_h_bjet_pt_beforecut_b  = new TH1D("h_bjet_pt_beforecut_b", "", 30, 0, 300);
  m_h_bjet_eta_beforecut_c = new TH1D("h_bjet_eta_beforecut_c", "", 32, -3.2, 3.2);
  m_h_bjet_pt_beforecut_c  = new TH1D("h_bjet_pt_beforecut_c", "", 30, 0, 300);
  m_h_bjet_eta_beforecut_l = new TH1D("h_bjet_eta_beforecut_l", "", 32, -3.2, 3.2);
  m_h_bjet_pt_beforecut_l  = new TH1D("h_bjet_pt_beforecut_l", "", 30, 0, 300);
  
  wk()->addOutput( m_h_RunNumber );
  wk()->addOutput( m_h_ZpT );
  wk()->addOutput( m_h_ZM );
  wk()->addOutput( m_h_nJets );
  wk()->addOutput( m_h_Z_jet_dPhi );
  wk()->addOutput( m_h_prwfactor );
  wk()->addOutput( m_h_njets_beforecut );
  wk()->addOutput( m_h_jet_eta_beforecut );
  wk()->addOutput( m_h_jet_pt_beforecut );
  wk()->addOutput( m_h_nbjets_beforecut );
  wk()->addOutput( m_h_bjet_eta_beforecut );
  wk()->addOutput( m_h_bjet_pt_beforecut );
  wk()->addOutput( m_h_bjet_eta_beforecut_b );
  wk()->addOutput( m_h_bjet_pt_beforecut_b );
  wk()->addOutput( m_h_bjet_eta_beforecut_c );
  wk()->addOutput( m_h_bjet_pt_beforecut_c );
  wk()->addOutput( m_h_bjet_eta_beforecut_l );
  wk()->addOutput( m_h_bjet_pt_beforecut_l );
  
  const std::pair<TH1F*, TH1F*> cutflows = ReturnCutflowPointers();
  int nBinsCutflow = cutflows.first->GetNbinsX();
  int xminCutflow  = cutflows.first->GetXaxis()->GetXmin();
  int xmaxCutflow  = cutflows.first->GetXaxis()->GetXmax();
  
  m_h_cutflow = new TH1F("cutflow", "", nBinsCutflow, xminCutflow, xmaxCutflow); //!
  m_h_cutflow_weighted = new TH1F("cutflow_weighted", "", nBinsCutflow, xminCutflow, xmaxCutflow); //!
  
  wk()->addOutput( m_h_cutflow );
  wk()->addOutput( m_h_cutflow_weighted );

  Info("histInitialize()", "ends");
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessZJetBalanceMiniTree :: fileExecute ()
{
  Info("fileExecute()", "called");
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  Info("fileExecute()", "end");
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessZJetBalanceMiniTree :: changeInput (bool firstFile)
{
  Info("changedInput", "called"); 
  
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  
  const std::pair<TH1F*, TH1F*> cutflows = ReturnCutflowPointers();  
  for (int iBin=1, nBins=cutflows.first->GetNbinsX(); iBin<=nBins; iBin++) { // update
    m_h_cutflow->SetBinContent(iBin, m_h_cutflow->GetBinContent(iBin)+cutflows.first->GetBinContent(iBin)); 
  }
  for (int iBin=1, nBins=cutflows.second->GetNbinsX(); iBin<=nBins; iBin++) { // update
    m_h_cutflow_weighted->SetBinContent(iBin, m_h_cutflow_weighted->GetBinContent(iBin)+cutflows.second->GetBinContent(iBin)); 
  }
  
  TTree *tree = wk()->tree();
  InitTree(tree);
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessZJetBalanceMiniTree :: initialize ()
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
  m_eventCounter = 0;
  
  // create object before configuration
  m_pT_binning = new Double_t[BUFSIZ];
  m_eta_binning = new Double_t[BUFSIZ];
  
  // configuration from ENV file
  if ( this->configure() == EL::StatusCode::FAILURE ) {
    Error("initialize()", "Failed to properly configure. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  // additinoal histograms according to configuration parameters
  for (int iPtBin=1; iPtBin<m_n_pT_binning+1; iPtBin++) {
    std::vector<TH1D*> tmp_hist_container;
    for (int iEtaBin=1; iEtaBin<m_n_eta_binning+1; iEtaBin++) {
      Info("Initialize()", "%s", Form("DB_RefEtaBin%d_PtBin%d", iEtaBin, iPtBin));
      TH1D* h = new TH1D(Form("DB_RefEtaBin%d_PtBin%d", iEtaBin, iPtBin), 
			 Form("%.1f<p_{T}<%.1f, %.1f<#eta<%.1f", 
			      m_pT_binning[iPtBin-1], m_pT_binning[iPtBin],
			      m_eta_binning[iEtaBin-1], m_eta_binning[iEtaBin]
			      ), 
			 m_nBinsXForResponseHist, m_minXForResponseHist, m_maxXForResponseHist);
      tmp_hist_container.push_back(h);
      wk()->addOutput( h );
    }
    m_balance_hists.push_back(tmp_hist_container);
  }
  m_h_jet_pt_bin = new TH2D("h_jet_pt_bin", "", 100, 0, 500, m_n_pT_binning, -0.5, -0.5+m_n_pT_binning); // validation purpose
  wk()->addOutput( m_h_jet_pt_bin );
  m_h_pt_binning_info = new TH1D("h_pt_binning_info", "", m_n_pT_binning, m_pT_binning);
  wk()->addOutput( m_h_pt_binning_info );
  m_h_eta_binning_info = new TH1D("h_eta_binning_info", "", m_n_eta_binning, m_eta_binning);
  wk()->addOutput( m_h_eta_binning_info );
  
  m_h_jet_eta = new TH1D("h_jet_eta", "", 25, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  m_h_jet_pt  = new TH1D("h_jet_pt",  "", 50, 0, 300.);
  m_h_bjet_eta = new TH1D("h_bjet_eta", "", 10, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  m_h_bjet_pt  = new TH1D("h_bjet_pt", "", 25, 0, 300.);

  m_h_bjet_eta_b = new TH1D("h_bjet_eta_b", "", 10, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  m_h_bjet_pt_b  = new TH1D("h_bjet_pt_b", "", 25, 0, 300.);
  m_h_bjet_eta_c = new TH1D("h_bjet_eta_c", "", 10, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  m_h_bjet_pt_c  = new TH1D("h_bjet_pt_c", "", 25, 0, 300.);
  m_h_bjet_eta_l = new TH1D("h_bjet_eta_l", "", 10, m_eta_binning[0], m_eta_binning[m_n_eta_binning]);
  m_h_bjet_pt_l  = new TH1D("h_bjet_pt_l", "", 25, 0, 300.);
  
  m_h_averageInteractionsPerCrossing = new TH1D("h_averageInteractionsPerCrossing", "", 50, 0, 50.);

  
  wk()->addOutput( m_h_jet_eta );
  wk()->addOutput( m_h_jet_pt );
  wk()->addOutput( m_h_averageInteractionsPerCrossing );
  wk()->addOutput( m_h_bjet_eta );
  wk()->addOutput( m_h_bjet_pt );
  
  wk()->addOutput( m_h_bjet_eta_b );
  wk()->addOutput( m_h_bjet_pt_b );
  wk()->addOutput( m_h_bjet_eta_c );
  wk()->addOutput( m_h_bjet_pt_c );
  wk()->addOutput( m_h_bjet_eta_l );
  wk()->addOutput( m_h_bjet_pt_l );

  // Pileup RW Tool //
  if ( m_doPUreweighting ) {
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
        PRWFiles.push_back(tmp_PRWFileNames.substr(0, pos));
        tmp_PRWFileNames.erase(0, pos);
      }else{
        PRWFiles.push_back(tmp_PRWFileNames.substr(0, pos));
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
    
    RETURN_CHECK("BasicEventSelection::initialize()", m_pileuptool->setProperty("ConfigFiles", PRWFiles), "");
    RETURN_CHECK("BasicEventSelection::initialize()", m_pileuptool->setProperty("LumiCalcFiles", lumiCalcFiles), "");
    RETURN_CHECK("BasicEventSelection::initialize()", m_pileuptool->initialize(), "");
  }  
  
  
  Info("initialize()", "Succesfully initialized! \n");
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ProcessZJetBalanceMiniTree :: execute ()
{
  // Here you do everything that needs to be done on every single
  // event, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  
  //===============================
  // load tree
  //===============================
  TTree* tree = wk()->tree();
  const int jentry = wk()->treeEntry();
  tree->LoadTree (jentry);
  tree->GetEntry (jentry);
    
  bool isMC = (mcChannelNumber!=-1);
  
  double weight_final=1.0;
  if (isMC) {
    double pileup_reweighting_factor = GetPileupReweightingFactor();
    m_h_prwfactor->Fill(pileup_reweighting_factor);
    weight_final = mcEventWeight*weight_xs*pileup_reweighting_factor;
    
    // Info("execute()", "mcEventWeight=%.4e weight_xs=%.4e pileup_factor=%.1e weight_final=%.1e",
    // 	 mcEventWeight, weight_xs, pileup_reweighting_factor, weight_final);
  }
  
  if(m_debug) Info("execute()", "Processing Event @ RunNumber=%10d, EventNumber=%d", runNumber, eventNumber);
  ++m_eventCounter;
  if (m_eventCounter%10000==0) {
    Info("execute()", "%10d th event is been processed.", m_eventCounter);
  }
  
  m_h_RunNumber->Fill(runNumber, weight_final);
  m_h_averageInteractionsPerCrossing->Fill(averageInteractionsPerCrossing, weight_final); // for validation
  
  // for valiadtion
  int nJetsBeforeCut = 0;
  int nBJetsBeforeCut = 0;
  for (int iJet=0, nJets=jet_pt->size(); iJet<nJets; iJet++) {
    const float& pt  = jet_pt->at(iJet);
    const float& eta = jet_eta->at(iJet);
    
    if ( pt < 30. ) continue;
    nJetsBeforeCut++;
    
    m_h_jet_eta_beforecut->Fill(eta, weight_final);
    m_h_jet_pt_beforecut->Fill(pt, weight_final);

    if ( jet_MV2c20->at(iJet) > m_MV2c20threshold ) {
      nBJetsBeforeCut++;
      m_h_bjet_eta_beforecut->Fill(eta, weight_final);
      m_h_bjet_pt_beforecut->Fill(pt, weight_final);
      
      double bsf = 1.0;
      const int& ConeTruthLabelID = jet_ConeTruthLabelID->at(iJet);
      
      if (TMath::Abs(ConeTruthLabelID)==5) {
	m_h_bjet_eta_beforecut_b->Fill(eta, weight_final*bsf);
	m_h_bjet_pt_beforecut_b->Fill(pt, weight_final*bsf);
      } else if (TMath::Abs(ConeTruthLabelID)==4) {
	m_h_bjet_eta_beforecut_c->Fill(eta, weight_final*bsf);
	m_h_bjet_pt_beforecut_c->Fill(pt, weight_final*bsf);
      } else {
	m_h_bjet_eta_beforecut_l->Fill(eta, weight_final*bsf);
	m_h_bjet_pt_beforecut_l->Fill(pt, weight_final*bsf);
      }
    }
    
  }
  m_h_njets_beforecut->Fill(nJetsBeforeCut, weight_final);
  m_h_nbjets_beforecut->Fill(nBJetsBeforeCut, weight_final);
  
  // selection criteria need to be applied
  if (TMath::Abs(ZM-91)>m_ZMassWindow) {return EL::StatusCode::SUCCESS;}
  if (TMath::Abs(dPhiZJet1)<m_cutDPhiZJet) {return EL::StatusCode::SUCCESS;}
  
  // 
  const float& lead_jet_pt      = jet_pt->at(0);
  const float& lead_jet_eta     = jet_eta->at(0);
  const int    lead_jet_pt_bin  = GetPtBin(pTRef1);
  const int    lead_jet_eta_bin = GetEtaBin(lead_jet_eta);
  
  //Info("execute()", "lead_jet_eta=%.1f (%d) lead_jet_pt=%.1f (%d)",
  //lead_jet_eta, lead_jet_eta_bin, lead_jet_pt, lead_jet_pt_bin);
  
  if (lead_jet_eta_bin==-1) {return EL::StatusCode::SUCCESS;} // out of eta range (defined as binning)
  
  if (jet_pt->size()>1) { if (jet_pt->at(1)>pTRef1*0.2) {return EL::StatusCode::SUCCESS;} } // event with second jet is vetoed
  
  m_h_ZpT->Fill(ZpT, weight_final);
  m_h_nJets->Fill(njets, weight_final);  // for validation
  m_h_ZM->Fill(ZM, weight_final);  // for validation
  m_h_jet_pt_bin->Fill(lead_jet_pt, lead_jet_pt_bin, weight_final);  // for validation
  m_h_Z_jet_dPhi->Fill(dPhiZJet1, weight_final); // for validation
  m_h_jet_eta->Fill(lead_jet_eta, weight_final);
  m_h_jet_pt->Fill(lead_jet_pt, weight_final);
  (m_balance_hists[lead_jet_pt_bin])[lead_jet_eta_bin]->Fill(lead_jet_pt/pTRef1, weight_final);
  
  if (jet_MV2c20->at(0)>m_MV2c20threshold) {
    m_h_bjet_eta->Fill(lead_jet_eta, weight_final);
    m_h_bjet_pt->Fill(lead_jet_pt, weight_final);
    
    double bsf = 1.0;
    const int& label = jet_ConeTruthLabelID->at(0);
    
    if (TMath::Abs(label)==5) {
      m_h_bjet_eta_b->Fill(lead_jet_eta, weight_final*bsf);
      m_h_bjet_pt_b->Fill(lead_jet_pt, weight_final*bsf);
    } else if (TMath::Abs(label)==4) {
      m_h_bjet_eta_c->Fill(lead_jet_eta, weight_final*bsf);
      m_h_bjet_pt_c->Fill(lead_jet_pt, weight_final*bsf);
    } else {
      m_h_bjet_eta_l->Fill(lead_jet_eta, weight_final*bsf);
      m_h_bjet_pt_l->Fill(lead_jet_pt, weight_final*bsf);
    }
  }
  
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

void ProcessZJetBalanceMiniTree::DecodeBinning(TString binning_str, Double_t* binning_array, Int_t& n_binning)
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

int ProcessZJetBalanceMiniTree::GetPtBin(const double& _pt)
{
  int rc=0;
  for (; rc<m_n_pT_binning-1; rc++) {
    if (_pt<m_pT_binning[rc+1]) {break;}
  }
  return rc;
}

int ProcessZJetBalanceMiniTree::GetEtaBin(const double& _eta)
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

double ProcessZJetBalanceMiniTree::GetPileupReweightingFactor()
{
  // what actually done in getCombinedWeight( const xAOD::EventInfo& eventInfo )
  // PileupReweighting/Root/PileupReweightingTool.cxx
  return ((m_doPUreweighting) ? 
	  m_pileuptool->GetCombinedWeight(runNumber, mcChannelNumber, averageInteractionsPerCrossing) :
	  1.0);
}

std::pair<TH1F*, TH1F*> ProcessZJetBalanceMiniTree::ReturnCutflowPointers()
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
  
  return rc;
}


void ProcessZJetBalanceMiniTree :: InitTree(TTree* tree)
{
  // Set object pointer
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
  tree->SetBranchAddress("pTRef1", &pTRef1, &b_pTRef1);
  tree->SetBranchAddress("dPhiZJet2", &dPhiZJet2, &b_dPhiZJet2);
  tree->SetBranchAddress("pTRef2", &pTRef2, &b_pTRef2);
  tree->SetBranchAddress("jetDPhi", &jetDPhi, &b_jetDPhi);
  tree->SetBranchAddress("jetDEta", &jetDEta, &b_jetDEta);
  tree->SetBranchAddress("jetPtRatio", &jetPtRatio, &b_jetPtRatio);
  tree->SetBranchAddress("weight", &weight, &b_weight);
  tree->SetBranchAddress("weight_xs", &weight_xs, &b_weight_xs);
  tree->SetBranchAddress("weight_prescale", &weight_prescale, &b_weight_prescale);
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
}
