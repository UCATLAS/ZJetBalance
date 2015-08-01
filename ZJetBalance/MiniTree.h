#ifndef AnalysisExample_MiniTree_H
#define AnalysisExample_MiniTree_H

#include "xAODAnaHelpers/HelpTreeBase.h"
#include "TTree.h"

class MiniTree : public HelpTreeBase
{

  private:

    float m_ZpT;
    float m_Zeta;
    float m_Zphi;
    float m_ZM;
    
    float m_dPhiZJet1;
    float m_dEtaZJet1;
    float m_pTRef1;
    float m_dPhiZJet2;
    float m_dEtaZJet2;
    float m_pTRef2;
    float m_jetDPhi;
    float m_jetDEta;
    float m_jetPtRatio;
    
    float m_weight;
    float m_weight_corr;
    float m_weight_xs;
    float m_weight_prescale;

    std::vector<float> m_jet_constitScaleEta;
    std::vector<float> m_jet_emScaleEta;
    std::vector<float> m_jet_mucorrected_pt;
    std::vector<float> m_jet_mucorrected_eta;
    std::vector<float> m_jet_mucorrected_phi;
    std::vector<float> m_jet_mucorrected_m;
    
  public:

    MiniTree(xAOD::TEvent * event, TTree* tree, TFile* file);
    ~MiniTree();

    void AddEventUser( const std::string detailStr = "" );
    void AddJetsUser( const std::string detailStr = "" );
    void FillEventUser( const xAOD::EventInfo* eventInfo );
    void FillJetsUser( const xAOD::Jet* jet );
    void ClearEventUser();
    void ClearJetsUser();

};
#endif
