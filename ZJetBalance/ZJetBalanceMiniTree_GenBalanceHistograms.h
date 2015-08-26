#ifndef ZJetBalanceMiniTree_GenBalanceHistograms_H
#define ZJetBalanceMiniTree_GenBalanceHistograms_H

#include <ZJetBalance/ZJetBalanceMiniTreeAnaBase.h>

//algorithm wrapper
#include <TH1F.h>
#include <TH2F.h>

class ZJetBalanceMiniTree_GenBalanceHistograms : public ZJetBalanceMiniTreeAnaBase
{
 public:
  // this is a standard constructor
  ZJetBalanceMiniTree_GenBalanceHistograms();
  
  virtual EL::StatusCode configure (); 
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(ZJetBalanceMiniTree_GenBalanceHistograms, 1);
  
 private:
  int    m_nBinsXForResponseHist; //! configurable parameter
  double m_maxXForResponseHist; //! configurable parameter
  double m_minXForResponseHist; //! configurable parameter
  bool   m_doPUreweighting; //! configurable parameter
  double m_cutDPhiZJet; //! configurable parameter
  double m_ZMassWindow; //! configurable parameter
  bool m_fillLeptonBefore; //! configurable parameter
  bool m_isMuonSample; //! set automatically for each input file

  // declaration for your histograms
  // Muon histograms
  TH1F* m_h_muon1_pT_beforecut; //!
  TH1F* m_h_muon1_eta_beforecut; //!
  TH1F* m_h_muon1_phi_beforecut; //!
  TH1F* m_h_muon2_pT_beforecut; //!
  TH1F* m_h_muon2_eta_beforecut; //!
  TH1F* m_h_muon2_phi_beforecut; //!
  TH1F* m_h_muon1_pT; //!
  TH1F* m_h_muon1_eta; //!
  TH1F* m_h_muon1_phi; //!
  TH1F* m_h_muon2_pT; //!
  TH1F* m_h_muon2_eta; //!
  TH1F* m_h_muon2_phi; //!

  // electrons histograms
  TH1F* m_h_electron1_pT_beforecut; //!
  TH1F* m_h_electron1_eta_beforecut; //!
  TH1F* m_h_electron1_phi_beforecut; //!
  TH1F* m_h_electron2_pT_beforecut; //!
  TH1F* m_h_electron2_eta_beforecut; //!
  TH1F* m_h_electron2_phi_beforecut; //!
  TH1F* m_h_electron1_pT; //!
  TH1F* m_h_electron1_eta; //!
  TH1F* m_h_electron1_phi; //!
  TH1F* m_h_electron2_pT; //!
  TH1F* m_h_electron2_eta; //!
  TH1F* m_h_electron2_phi; //!

  // histograms
  TH1F* m_h_RunNumber; //!
  TH1F* m_h_ZpT; //!
  TH1F* m_h_Zeta; //!
  TH1F* m_h_Zphi; //!
  TH1F* m_h_ZM; //!
  TH1F* m_h_ZpTRef; //!
  TH1F* m_h_ZpT_beforecut; //!
  TH1F* m_h_Zeta_beforecut; //!
  TH1F* m_h_Zphi_beforecut; //!
  TH1F* m_h_ZM_beforecut; //!
  TH1F* m_h_Z_jet_dPhi; //!
  TH1F* m_h_Z_jet_dEta; //!
  TH1F* m_h_nJets; //!
  TH1F* m_h_jet_eta; //!
  TH1F* m_h_jet_pt; //!
  TH1F* m_h_jet_phi; //!
  TH1F* m_h_jet_eta_b; //!
  TH1F* m_h_jet_pt_b; //!
  TH1F* m_h_jet_phi_b; //!
  TH1F* m_h_jet_eta_c; //!
  TH1F* m_h_jet_pt_c; //!
  TH1F* m_h_jet_phi_c; //!
  TH1F* m_h_jet_eta_l; //!
  TH1F* m_h_jet_pt_l; //!
  TH1F* m_h_jet_phi_l; //!
  TH1F* m_h_averageInteractionsPerCrossing; //!
  TH2F* m_h_jet_pt_bin; //!
  TH1F* m_h_pt_binning_info; //!
  TH1F* m_h_eta_binning_info; //!
  TH1F* m_h_prwfactor; //!
  TH1F* m_h_muonTrigFactor; //!
  TH1F* m_h_muon1EffFactor; //!
  TH1F* m_h_muon2EffFactor; //!
  TH1F* m_h_electronTrigFactor; //!
  TH1F* m_h_electron1EffFactor; //!
  TH1F* m_h_electron2EffFactor; //!
  TH1F* m_h_nJets_beforecut; //!
  TH1F* m_h_jet_eta_beforecut; //!
  TH1F* m_h_jet_pt_beforecut; //!
  TH1F* m_h_jet_eta_beforecut_b; //!
  TH1F* m_h_jet_pt_beforecut_b; //!
  TH1F* m_h_jet_eta_beforecut_c; //!
  TH1F* m_h_jet_pt_beforecut_c; //!
  TH1F* m_h_jet_eta_beforecut_l; //!
  TH1F* m_h_jet_pt_beforecut_l; //!
  TH1F* m_h_1st_jet_eta_beforecut; //!
  TH1F* m_h_1st_jet_pt_beforecut; //!
  TH1F* m_h_1st_jet_eta_beforecut_b; //!
  TH1F* m_h_1st_jet_pt_beforecut_b; //!
  TH1F* m_h_1st_jet_eta_beforecut_c; //!
  TH1F* m_h_1st_jet_pt_beforecut_c; //!
  TH1F* m_h_1st_jet_eta_beforecut_l; //!
  TH1F* m_h_1st_jet_pt_beforecut_l; //!
  TH1F* m_h_SumPtTrkPt500PV_beforecut; //!
  TH1F* m_h_SumPtTrkPt500PV; //!
  TH1F* m_h_cutflow; //!
  TH1F* m_h_cutflow_weighted; //!
  TH1F* m_h_cutflow_weighted_final; //!
  
  
  Double_t*       m_pT_binning; //!
  Int_t           m_n_pT_binning; //!
  Double_t*       m_eta_binning; //!
  Int_t           m_n_eta_binning; //!

  static void DecodeBinning(TString binning_str, Double_t* binning_array, Int_t& n_binning);
  int GetPtBin(const double& _pt);
  int GetEtaBin(const double& _eta);
  inline void FillFlavorHistograms(TH1F* h_b, TH1F* h_c, TH1F* h_l, const int& truthLabel, const float& value, const float& weight);
  void IsMuonSample();

  std::pair<TH1F*, TH1F*> ReturnCutflowPointers();
  std::vector< std::vector<TH1F*> > m_balance_hists;
  std::vector< std::vector<TH1F*> > m_balance_hists_b;
  std::vector< std::vector<TH1F*> > m_balance_hists_c;
  std::vector< std::vector<TH1F*> > m_balance_hists_l;
};  
#endif
