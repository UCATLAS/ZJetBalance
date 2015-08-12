#include <ZJetBalance/DrawingHelperOk.h>

#include <string>

int main()
{
  ZJetBalance::DrawingHelperOk drawer("output");
  
  drawer.SetDataFileName("/afs/cern.ch/work/o/okumura/public/analysi2015/ZBjetBalance/runscripts/all_histData50ns.root");
  drawer.SetLuminosity(0.0783190);
  drawer.AddMC("/afs/cern.ch/work/o/okumura/public/analysi2015/ZBjetBalance/runscripts/all_hist361107.root",
	       "Z#rightarrow#mu#mu",
	       "Zmumu",
	       kYellow, 8, true);
  drawer.AddMC("/afs/cern.ch/work/o/okumura/public/analysi2015/ZBjetBalance/runscripts/all_hist361108.root",
  	       "Z#rightarrow#tau#tau",
	       "Ztautau",
  	       kYellow+3, 8, true);
  
  drawer.AddMC("/afs/cern.ch/work/o/okumura/public/analysi2015/ZBjetBalance/runscripts/all_hist410004.root",
	       "t#bar{t}",
	       "Tt",
	       kGreen, 8, true);
  
  const std::string period("C2-C4");
  drawer.MyDataMcComparisonTH1F("ZpT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{Z} [GeV]");
  drawer.MyDataMcComparisonTH1F("ZM", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "M_{Z} [GeV]");
  drawer.MyDataMcComparisonTH1F("Z_jet_dPhi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#Delta (jet, Z)");
  drawer.MyDataMcComparisonTH1F("muon1_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{muon1}");
  drawer.MyDataMcComparisonTH1F("muon1_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#phi^{muon1}");
  drawer.MyDataMcComparisonTH1F("muon1_pT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{muon1}");
  drawer.MyDataMcComparisonTH1F("muon2_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{muon2}");
  drawer.MyDataMcComparisonTH1F("muon2_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#phi^{muon2}");
  drawer.MyDataMcComparisonTH1F("muon2_pT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{muon2}");
  drawer.MyDataMcComparisonTH1F("jet_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{jet}");
  drawer.MyDataMcComparisonTH1F("jet_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#phi^{jet}");
  drawer.MyDataMcComparisonTH1F("jet_pt", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{jet}");
  drawer.MyDataMcComparisonTH1F("averageInteractionsPerCrossing", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#mu (preselection) ");
  drawer.MyDataMcComparisonTH1F("1st_jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "leading jet #eta (preselection)");
  drawer.MyDataMcComparisonTH1F("1st_jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "leading jet p_{T} (preselection)");
  drawer.MyDataMcComparisonTH1F("nJets_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "N_{jets} (preselection)");
  drawer.MyDataMcComparisonTH1F("jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{jet} (preselection)");
  drawer.MyDataMcComparisonTH1F("jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{jet} (preselection)");
  drawer.MyDataMcComparisonTH1F("ZpT_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{Z} [GeV]");
  drawer.MyDataMcComparisonTH1F("Zeta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{Z}");
  drawer.MyDataMcComparisonTH1F("Zphi_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#phi^{Z}");
  drawer.MyDataMcComparisonTH1F("ZM_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "M_{Z} [GeV]");
  
  
  drawer.DumpCutFlow("DiMuon");
}
