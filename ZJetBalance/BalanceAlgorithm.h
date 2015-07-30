#ifndef ZJetBalance_BalanceAlgorithm_H
#define ZJetBalance_BalanceAlgorithm_H

#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>

//algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

#ifndef __CINT__
  #include "xAODJet/JetContainer.h"
  #include "xAODMuon/MuonContainer.h"
  #include "xAODEventInfo/EventInfo.h"
  #include "ZJetBalance/MiniTree.h"
#endif

// ROOT include(s):
#include "TH1D.h"

#include <sstream>

class BalanceAlgorithm : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  public:
    // float cutValue;
    int m_eventCounter;     //!

    TH1D* m_cutflowHist;    //!
    TH1D* m_cutflowHistW;   //!
    int m_cutflowFirst;     //!
    int m_iCutflow;         //!
    float m_mcEventWeight;  //!

  private:
    //configuration variables
    std::string m_inputJetContainerName;    //! input container name
    std::string m_inputJetAlgo;             //! input algo for when running systs
    std::string m_inputMuonContainerName;   //! input container name
    std::string m_inputMuonAlgo;            //! input algo for when running systs
    std::string m_inputMuonForMuonInJetCorrectionContainerName; //! input container name
    std::string m_inputMuonForMuonInJetCorrectionAlgo;          //! input algo for when running systs
    bool m_isMC;                      //! Is MC
    bool m_useCutFlow;                //! true will write out cutflow histograms
    bool m_writeTree;                 //! true will write out a TTree
    float m_leadingJetPtCut;          //! Leading jet Pt cut
    bool m_truthLevelOnly;            //! truthLevelOnly info
    std::string m_eventDetailStr;     //! event info add to tree
    std::string m_trigDetailStr;      //! trigger info add to tree
    std::string m_jetDetailStr;       //! jet info add to tree
    std::string m_jetDetailStrSyst;   //! jetsyst info add to tree
    std::string m_muonDetailStr;      //! muon info add to tree
    std::string m_muonDetailStrSyst;  //! muonsyst info add to tree

    float m_xs; //!
    float m_filtEff; //!
    int m_numAMIEvents; //!
    int m_mcChannelNumber; //!
    std::stringstream m_ss; //!

    const float GeV = 1000.;
    std::string m_treeStream;
    void passCut();
    TLorentzVector getFourMomentumOfMuonInJet(const xAOD::Muon* muon);
    EL::StatusCode muonInJetCorrection (const xAOD::JetContainer* signalJets);

#ifndef __CINT__
    EL::StatusCode getLumiWeights(const xAOD::EventInfo* eventInfo);
    std::map< std::string, MiniTree* > m_myTrees; //!
#endif // not __CINT__

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  // this is a standard constructor
  BalanceAlgorithm ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // these are the functions not inherited from Algorithm
  virtual EL::StatusCode configure ();
#ifndef __CINT__
  bool executeAnalysis( const xAOD::EventInfo* eventInfo,
      const xAOD::JetContainer* signalJets,
      const xAOD::MuonContainer* signalMuons,
      const xAOD::VertexContainer* vertices,
      bool count,
      std::string systName = "");
#endif // not __CINT__
  void AddTree( std::string );

  // this is needed to distribute the algorithm to the workers
  ClassDef(BalanceAlgorithm, 1);
};

#endif
