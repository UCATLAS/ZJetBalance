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

    float m_weight;
    float m_weight_corr;
    float m_weight_xs;
    float m_weight_prescale;

    std::vector<float> m_jet_constitScaleEta;
    std::vector<float> m_jet_emScaleEta;


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
