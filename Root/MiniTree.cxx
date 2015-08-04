#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEventInfo/EventInfo.h"

#include "ZJetBalance/MiniTree.h"

MiniTree :: MiniTree(xAOD::TEvent * event, TTree* tree, TFile* file) :
  HelpTreeBase(event, tree, file, 1e3)
{
  Info("MiniTree", "Creating output TTree");
}

MiniTree :: ~MiniTree()
{
}
//////////////////// Connect Defined variables to branches here /////////////////////////////
void MiniTree::AddEventUser(const std::string detailStr)
{
  // event variables
  m_tree->Branch("ZpT",   &m_ZpT,   "ZpT/F"  );
  m_tree->Branch("Zeta",  &m_Zeta,  "Zeta/F" );
  m_tree->Branch("Zphi",  &m_Zphi,  "Zphi/F" );
  m_tree->Branch("ZM",    &m_ZM,    "ZM/F"   );
  m_tree->Branch("dPhiZJet1",  &m_dPhiZJet1,  "dPhiZJet1/F");
  m_tree->Branch("dEtaZJet1",  &m_dEtaZJet1,  "dEtaZJet1/F");
  m_tree->Branch("pTRef1",     &m_pTRef1,     "pTRef1/F");
  m_tree->Branch("dPhiZJet2",  &m_dPhiZJet2,  "dPhiZJet2/F");
  m_tree->Branch("dEtaZJet2",  &m_dEtaZJet2,  "dEtaZJet2/F");
  m_tree->Branch("pTRef2",     &m_pTRef2,     "pTRef2/F");
  m_tree->Branch("jetDPhi",    &m_jetDPhi,    "jetDPhi/F");
  m_tree->Branch("jetDEta",    &m_jetDEta,    "jetDEta/F");
  m_tree->Branch("jetPtRatio", &m_jetPtRatio, "jetPtRatio/F");


  m_tree->Branch("weight", &m_weight, "weight/F");
  m_tree->Branch("weight_xs", &m_weight_xs, "weight_xs/F");
  m_tree->Branch("weight_prescale", &m_weight_prescale, "weight_prescale/F");
}

void MiniTree::AddJetsUser(const std::string detailStr)
{
  m_tree->Branch("jet_constitScaleEta", &m_jet_constitScaleEta);
  m_tree->Branch("jet_emScaleEta", &m_jet_emScaleEta);
  m_tree->Branch("jet_mucorrected_pt", &m_jet_mucorrected_pt);
  m_tree->Branch("jet_mucorrected_eta", &m_jet_mucorrected_eta);
  m_tree->Branch("jet_mucorrected_phi", &m_jet_mucorrected_phi);
  m_tree->Branch("jet_mucorrected_m",  &m_jet_mucorrected_m);
}

//////////////////// Clear any defined vectors here ////////////////////////////
void MiniTree::ClearEventUser() {
  m_ZpT      = -999;
  m_Zeta     = -999;
  m_Zphi     = -999;
  m_ZM       = -999;

  m_dPhiZJet1  = -999.;
  m_dEtaZJet1  = -999.;
  m_pTRef1     = -999.;
  m_dPhiZJet2  = -999.;
  m_dEtaZJet2  = -999.;
  m_pTRef2     = -999.;
  m_jetDPhi    = -999.;
  m_jetDEta    = -999.;
  m_jetPtRatio = -999.;
  
  m_weight_corr = -999;
  m_weight    = -999;
  m_weight_xs = -999;
  m_weight_prescale = -999;
}

void MiniTree::ClearJetsUser() {
  m_jet_constitScaleEta.clear();
  m_jet_emScaleEta.clear();
  m_jet_mucorrected_pt.clear();
  m_jet_mucorrected_eta.clear();
  m_jet_mucorrected_phi.clear();
  m_jet_mucorrected_m.clear();
}

/////////////////// Assign values to defined event variables here ////////////////////////
void MiniTree::FillEventUser( const xAOD::EventInfo* eventInfo ) {
  
  
  if( eventInfo->isAvailable< float >( "ZpT" ) )
    m_ZpT = eventInfo->auxdecor< float >( "ZpT" );
  if( eventInfo->isAvailable< float >( "Zeta" ) )
    m_Zeta = eventInfo->auxdecor< float >( "Zeta" );
  if( eventInfo->isAvailable< float >( "Zphi" ) )
    m_Zphi = eventInfo->auxdecor< float >( "Zphi" );
  if( eventInfo->isAvailable< float >( "ZM" ) )
    m_ZM = eventInfo->auxdecor< float >( "ZM" );
  if( eventInfo->isAvailable< float >( "dPhiZJet1" ) ) 
    m_dPhiZJet1 = eventInfo->auxdecor< float >( "dPhiZJet1" );
  if( eventInfo->isAvailable< float >( "dEtaZJet1" ) ) 
    m_dEtaZJet1 = eventInfo->auxdecor< float >( "dEtaZJet1" );
  if( eventInfo->isAvailable< float >( "pTRef1" ) ) 
    m_pTRef1    = eventInfo->auxdecor< float >( "pTRef1" );
  if( eventInfo->isAvailable< float >( "dPhiZJet2" ) ) 
    m_dPhiZJet2 = eventInfo->auxdecor< float >( "dPhiZJet2" );
  if( eventInfo->isAvailable< float >( "dEtaZJet2" ) ) 
    m_dEtaZJet2 = eventInfo->auxdecor< float >( "dEtaZJet2" );
  if( eventInfo->isAvailable< float >( "pTRef2" ) ) 
    m_pTRef2    = eventInfo->auxdecor< float >( "pTRef2" );
  if( eventInfo->isAvailable< float >( "jetDPhi" ) ) 
    m_jetDPhi   = eventInfo->auxdecor< float >( "jetDPhi" );
  if( eventInfo->isAvailable< float >( "jetDEta" ) ) 
    m_jetDEta   = eventInfo->auxdecor< float >( "jetDEta" );
  if( eventInfo->isAvailable< float >( "jetPtRatio" ) ) 
    m_jetPtRatio = eventInfo->auxdecor< float >( "jetPtRatio" );
  

  
  if( eventInfo->isAvailable< float >( "weight" ) )
    m_weight = eventInfo->auxdecor< float >( "weight" );
  if( eventInfo->isAvailable< float >( "weight_xs" ) )
    m_weight_xs = eventInfo->auxdecor< float >( "weight_xs" );
  if( eventInfo->isAvailable< float >( "weight_prescale" ) )
    m_weight_prescale = eventInfo->auxdecor< float >( "weight_prescale" );

}

/////////////////// Assign values to defined jet variables here //////////////////
void MiniTree::FillJetsUser( const xAOD::Jet* jet ) {
  if( jet->isAvailable< float >( "constitScaleEta" ) ) {
    m_jet_constitScaleEta.push_back( jet->auxdata< float >("constitScaleEta") );
  } else {
    m_jet_constitScaleEta.push_back( -999 );
  }
  if( jet->isAvailable< float >( "emScaleEta" ) ) {
    m_jet_emScaleEta.push_back( jet->auxdata< float >("emScaleEta") );
  } else {
    m_jet_emScaleEta.push_back( -999 );
  }
  
  m_jet_mucorrected_pt.push_back(jet->auxdata< float >("mucorrected_pt")/m_units);
  m_jet_mucorrected_eta.push_back(jet->auxdata< float >("mucorrected_eta"));
  m_jet_mucorrected_phi.push_back(jet->auxdata< float >("mucorrected_phi"));
  m_jet_mucorrected_m.push_back(jet->auxdata< float >("mucorrected_m")/m_units);
}

