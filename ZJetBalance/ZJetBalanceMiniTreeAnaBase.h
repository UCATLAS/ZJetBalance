#ifndef ZJetBalanceMiniTreeAnaBase_H
#define ZJetBalanceMiniTreeAnaBase_H

#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>
#include <EventLoop/Worker.h>

//algorithm wrapper
#include <xAODAnaHelpers/Algorithm.h>
#include "PileupReweighting/PileupReweightingTool.h"


#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TBranch.h>

#include <map>
#include <string>
#include <vector>

class ZJetBalanceMiniTreeAnaBase : public xAH::Algorithm
{
 protected:
  // float cutValue;
  int m_eventCounter;  //!
  std::string m_histPrefix;  //!
  bool m_isMC;
  
  bool m_debug; //! set verbose mode
  bool m_doPUreweighting; //! configurable parameter
  std::string m_lumiCalcFileNames; //! configurable parameter
  std::string m_PRWFileNames; //! configurable parameter
  bool m_btagJets;      //! configurable parameter
  double m_additional_weight;   //! configurable parameter
  std::string m_btagOP; //! configurable parameter
  
  // tools 
  CP::PileupReweightingTool* m_pileuptool; //!
  
  // histograms
  TH1F* m_h_cutflow; //!
  TH1F* m_h_cutflow_weighted; //!
  TH1F* m_h_cutflow_weighted_final; //!
  
 public: 
  // this is a standard constructor
  ZJetBalanceMiniTreeAnaBase();
 
  virtual EL::StatusCode configure (); 
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();
  
 public:
  // this is needed to distribute the algorithm to the workers
  ClassDef(ZJetBalanceMiniTreeAnaBase, 1);
  
 protected:
  // common functions to be used
  std::vector< TH1* > m_allHists; //!

  // histogram functions
  TH1F* book(std::string title,
      std::string xlabel, int xbins, double xlow, double xhigh);

  TH2F* book(std::string title,
      std::string xlabel, int xbins, double xlow, double xhigh,
      std::string xyabel, int ybins, double ylow, double yhigh);
  //// Variable Binned Histograms ////
  TH1F* book(std::string title,
      std::string xlabel, int xbins, const Double_t* xbinsArr);

  TH2F* book(std::string title,
      std::string xlabel, int xbins, const Double_t* xbinsArr,
      std::string ylabel, int ybins, double ylow, double yhigh);
  TH2F* book(std::string title,
      std::string xyabel, int xbins, double xlow, double xhigh,
      std::string ylabel, int ybins, const Double_t* ybinsArr);
  TH2F* book(std::string title,
      std::string xyabel, int xbins, const Double_t* xbinsArr,
      std::string ylabel, int ybins, const Double_t* ybinsArr);
  
  // Record all histograms from m_allHists to the worker
  void record(EL::Worker* wk);
  
  // Turn on Sumw2 for the histogram
  void Sumw2(TH1* hist, bool flag=true);

  // Push the new histogram to m_allHists
  void record(TH1* hist);

  // Set the xlabel
  void SetLabel(TH1* hist, std::string xlabel);
  
  // Set the xlabel, ylabel
  void SetLabel(TH1* hist, std::string xlabel, std::string ylabel);
  
  // Return Cutflow pointer
  std::pair<TH1F*, TH1F*> ReturnCutflowPointers();
  
  // initialize cutflow (please call it in histInitialize())
  void InitializeCutflows();
  
  // fill three different types of histograms
  void FillCutflowHistograms(const std::string& label, const double& xAHWeight, const double& weightFinal);
  
  // initialize pileup reweighting tool (please call it in initialize())
  EL::StatusCode InitializePileupReweightingTool();
  
  // get pileup reweight handled with 
  double GetPileupReweightingFactor(); 
  
  // load mini tree please call this in initialize 
  void LoadMiniTree();
  
  // load basic configuration
  EL::StatusCode LoadBasicConfiguration();
  
 protected:
  // commmon functin related to TTree access
  void            SetBTagAddresses(TTree* tree);  
  void            InitTree(TTree* tree);
  // copied from MakeClass function and add //! for all the variables
  // Declaration of leaf types

  // Declaration of leaf types
  Int_t           runNumber; //!
  Int_t           eventNumber; //!
  Int_t           mcEventNumber; //!
  Int_t           mcChannelNumber; //!
  Float_t         mcEventWeight; //!
  Float_t         weight_pileup; //!
  Int_t           NPV; //!
  Float_t         actualInteractionsPerCrossing; //!
  Float_t         averageInteractionsPerCrossing; //!
  Int_t           lumiBlock; //!
  Double_t        rhoEM; //!
  Int_t           pdgId1; //!
  Int_t           pdgId2; //!
  Int_t           pdfId1; //!
  Int_t           pdfId2; //!
  Float_t         x1; //!
  Float_t         x2; //!
  Float_t         xf1; //!
  Float_t         xf2; //!
  Float_t         ZpT; //!
  Float_t         Zeta; //!
  Float_t         Zphi; //!
  Float_t         ZM; //!
  Float_t         dPhiZJet1; //!
  Float_t         dEtaZJet1; //!
  Float_t         pTRef1; //!
  Float_t         dPhiZJet2; //!
  Float_t         dEtaZJet2; //!
  Float_t         pTRef2; //!
  Float_t         jetDPhi; //!
  Float_t         jetDEta; //!
  Float_t         jetPtRatio; //!
  Float_t         weight; //!
  Float_t         weight_xs; //!
  Float_t         weight_prescale; //!
  std::vector<double>  *weight_muon_trig; //!
  Int_t           njets; //!
  std::vector<float>   *jet_E; //!
  std::vector<float>   *jet_pt; //!
  std::vector<float>   *jet_phi; //!
  std::vector<float>   *jet_eta; //!
  std::vector<float>   *jet_rapidity; //!
  std::vector<float>   *jet_HECFrac; //!
  std::vector<float>   *jet_EMFrac; //!
  std::vector<float>   *jet_CentroidR; //!
  std::vector<float>   *jet_FracSamplingMax; //!
  std::vector<float>   *jet_FracSamplingMaxIndex; //!
  std::vector<float>   *jet_LowEtConstituentsFrac; //!
  std::vector<float>   *jet_GhostMuonSegmentCount; //!
  std::vector<float>   *jet_Width; //!
  std::vector<int>     *jet_NumTrkPt1000PV; //!
  std::vector<float>   *jet_SumPtTrkPt1000PV; //!
  std::vector<float>   *jet_TrackWidthPt1000PV; //!
  std::vector<int>     *jet_NumTrkPt500PV; //!
  std::vector<float>   *jet_SumPtTrkPt500PV; //!
  std::vector<float>   *jet_TrackWidthPt500PV; //!
  std::vector<float>   *jet_JVFPV; //!
  std::vector<float>   *jet_Jvt; //!
  std::vector<float>   *jet_JvtJvfcorr; //!
  std::vector<float>   *jet_JvtRpt; //!
  std::vector<float>   *jet_SV0; //!
  std::vector<float>   *jet_SV1; //!
  std::vector<float>   *jet_IP3D; //!
  std::vector<float>   *jet_SV1IP3D; //!
  std::vector<float>   *jet_MV1; //!
  std::vector<float>   *jet_MV2c00; //!
  std::vector<float>   *jet_MV2c20; //!
  std::vector<int>     *jet_MV2c20_is85; //!
  std::vector<std::vector<float> > *jet_MV2c20_SF85; //!
  std::vector<int>     *jet_MV2c20_is77; //!
  std::vector<std::vector<float> > *jet_MV2c20_SF77; //!
  std::vector<int>     *jet_MV2c20_is70; //!
  std::vector<std::vector<float> > *jet_MV2c20_SF70; //!
  std::vector<int>     *jet_MV2c20_is60; //!
  std::vector<std::vector<float> > *jet_MV2c20_SF60; //!
  std::vector<int>     *jet_isBTag; //!
  std::vector<std::vector<float> > *jet_SFBTag; //!
  std::vector<float>   *jet_GhostArea; //!
  std::vector<float>   *jet_ActiveArea; //!
  std::vector<float>   *jet_VoronoiArea; //!
  std::vector<float>   *jet_ActiveArea4vec_pt; //!
  std::vector<float>   *jet_ActiveArea4vec_eta; //!
  std::vector<float>   *jet_ActiveArea4vec_phi; //!
  std::vector<float>   *jet_ActiveArea4vec_m; //!
  std::vector<int>     *jet_ConeTruthLabelID; //!
  std::vector<int>     *jet_TruthCount; //!
  std::vector<float>   *jet_TruthLabelDeltaR_B; //!
  std::vector<float>   *jet_TruthLabelDeltaR_C; //!
  std::vector<float>   *jet_TruthLabelDeltaR_T; //!
  std::vector<int>     *jet_PartonTruthLabelID; //!
  std::vector<float>   *jet_GhostTruthAssociationFraction; //!
  std::vector<float>   *jet_truth_E; //!
  std::vector<float>   *jet_truth_pt; //!
  std::vector<float>   *jet_truth_phi; //!
  std::vector<float>   *jet_truth_eta; //!
  std::vector<float>   *jet_constitScaleEta; //!
  std::vector<float>   *jet_emScaleEta; //!
  std::vector<float>   *jet_mucorrected_pt; //!
  std::vector<float>   *jet_mucorrected_eta; //!
  std::vector<float>   *jet_mucorrected_phi; //!
  std::vector<float>   *jet_mucorrected_m; //!
  std::vector<float>   *muon_pt; //!
  std::vector<float>   *muon_eta; //!
  std::vector<float>   *muon_phi; //!
  std::vector<float>   *muon_m; //!
  std::vector< std::vector<double> >  *muon_effSF; //!
  Int_t           nmuon; //!

  // List of branches
  TBranch        *b_runNumber;   //!
  TBranch        *b_eventNumber;   //!
  TBranch        *b_mcEventNumber;   //!
  TBranch        *b_mcChannelNumber;   //!
  TBranch        *b_mcEventWeight;   //!
  TBranch        *b_weight_pileup;   //!
  TBranch        *b_NPV;   //!
  TBranch        *b_actualInteractionsPerCrossing;   //!
  TBranch        *b_averageInteractionsPerCrossing;   //!
  TBranch        *b_lumiBlock;   //!
  TBranch        *b_rhoEM;   //!
  TBranch        *b_pdgId1;   //!
  TBranch        *b_pdgId2;   //!
  TBranch        *b_pdfId1;   //!
  TBranch        *b_pdfId2;   //!
  TBranch        *b_x1;   //!
  TBranch        *b_x2;   //!
  TBranch        *b_xf1;   //!
  TBranch        *b_xf2;   //!
  TBranch        *b_ZpT;   //!
  TBranch        *b_Zeta;   //!
  TBranch        *b_Zphi;   //!
  TBranch        *b_ZM;   //!
  TBranch        *b_dPhiZJet1;   //!
  TBranch        *b_dEtaZJet1;   //!
  TBranch        *b_pTRef1;   //!
  TBranch        *b_dPhiZJet2;   //!
  TBranch        *b_dEtaZJet2;   //!
  TBranch        *b_pTRef2;   //!
  TBranch        *b_jetDPhi;   //!
  TBranch        *b_jetDEta;   //!
  TBranch        *b_jetPtRatio;   //!
  TBranch        *b_weight;   //!
  TBranch        *b_weight_xs;   //!
  TBranch        *b_weight_prescale;   //!
  TBranch        *b_weight_muon_trig;   //!
  TBranch        *b_njets;   //!
  TBranch        *b_jet_E;   //!
  TBranch        *b_jet_pt;   //!
  TBranch        *b_jet_phi;   //!
  TBranch        *b_jet_eta;   //!
  TBranch        *b_jet_rapidity;   //!
  TBranch        *b_jet_HECFrac;   //!
  TBranch        *b_jet_EMFrac;   //!
  TBranch        *b_jet_CentroidR;   //!
  TBranch        *b_jet_FracSamplingMax;   //!
  TBranch        *b_jet_FracSamplingMaxIndex;   //!
  TBranch        *b_jet_LowEtConstituentsFrac;   //!
  TBranch        *b_jet_GhostMuonSegmentCount;   //!
  TBranch        *b_jet_Width;   //!
  TBranch        *b_jet_NumTrkPt1000PV;   //!
  TBranch        *b_jet_SumPtTrkPt1000PV;   //!
  TBranch        *b_jet_TrackWidthPt1000PV;   //!
  TBranch        *b_jet_NumTrkPt500PV;   //!
  TBranch        *b_jet_SumPtTrkPt500PV;   //!
  TBranch        *b_jet_TrackWidthPt500PV;   //!
  TBranch        *b_jet_JVFPV;   //!
  TBranch        *b_jet_Jvt;   //!
  TBranch        *b_jet_JvtJvfcorr;   //!
  TBranch        *b_jet_JvtRpt;   //!
  TBranch        *b_jet_SV0;   //!
  TBranch        *b_jet_SV1;   //!
  TBranch        *b_jet_IP3D;   //!
  TBranch        *b_jet_SV1IP3D;   //!
  TBranch        *b_jet_MV1;   //!
  TBranch        *b_jet_MV2c00;   //!
  TBranch        *b_jet_MV2c20;   //!
  TBranch        *b_jet_MV2c20_is85;   //!
  TBranch        *b_jet_MV2c20_SF85;   //!
  TBranch        *b_jet_MV2c20_is77;   //!
  TBranch        *b_jet_MV2c20_SF77;   //!
  TBranch        *b_jet_MV2c20_is70;   //!
  TBranch        *b_jet_MV2c20_SF70;   //!
  TBranch        *b_jet_MV2c20_is60;   //!
  TBranch        *b_jet_MV2c20_SF60;   //!
  TBranch        *b_jet_isBTag;   //!
  TBranch        *b_jet_SFBTag;   //!
  TBranch        *b_jet_GhostArea;   //!
  TBranch        *b_jet_ActiveArea;   //!
  TBranch        *b_jet_VoronoiArea;   //!
  TBranch        *b_jet_ActiveArea4vec_pt;   //!
  TBranch        *b_jet_ActiveArea4vec_eta;   //!
  TBranch        *b_jet_ActiveArea4vec_phi;   //!
  TBranch        *b_jet_ActiveArea4vec_m;   //!
  TBranch        *b_jet_ConeTruthLabelID;   //!
  TBranch        *b_jet_TruthCount;   //!
  TBranch        *b_jet_TruthLabelDeltaR_B;   //!
  TBranch        *b_jet_TruthLabelDeltaR_C;   //!
  TBranch        *b_jet_TruthLabelDeltaR_T;   //!
  TBranch        *b_jet_PartonTruthLabelID;   //!
  TBranch        *b_jet_GhostTruthAssociationFraction;   //!
  TBranch        *b_jet_truth_E;   //!
  TBranch        *b_jet_truth_pt;   //!
  TBranch        *b_jet_truth_phi;   //!
  TBranch        *b_jet_truth_eta;   //!
  TBranch        *b_jet_constitScaleEta;   //!
  TBranch        *b_jet_emScaleEta;   //!
  TBranch        *b_jet_mucorrected_pt;   //!
  TBranch        *b_jet_mucorrected_eta;   //!
  TBranch        *b_jet_mucorrected_phi;   //!
  TBranch        *b_jet_mucorrected_m;   //!
  TBranch        *b_nmuon;   //!
  TBranch        *b_muon_pt;   //!
  TBranch        *b_muon_eta;   //!
  TBranch        *b_muon_phi;   //!
  TBranch        *b_muon_m;   //!
  TBranch        *b_muon_effSF;   //!
};  
#endif
