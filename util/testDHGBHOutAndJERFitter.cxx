#include <ZJetBalance/DH_GBHOut.h>

#include <string>

int main()
{
  const double luminosity=0.0783190;
  const std::string label="Internal";
  
  ZJetBalance::DH_GBHOut drawer("output");
  
  drawer.SetDataFileName("all_histData50ns.root");
  drawer.SetLuminosity(luminosity);
  drawer.AddMC("all_hist361107.root",
	       "Z#rightarrow#mu#mu",
	       "Zmumu",
	       kYellow, 8, true);
  drawer.AddMC("all_hist361108.root",
  	       "Z#rightarrow#tau#tau",
	       "Ztautau",
  	       kYellow+3, 8, true);
  
  drawer.AddMC("all_hist410004.root",
	       "t#bar{t}",
	       "Tt",
	       kGreen, 8, true);
  
  const std::string period("C2-C4");
  drawer.MyDataMcComparisonTH1F("ZpT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{Z} [GeV]", label);
  drawer.MyDataMcComparisonTH1F("ZM", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "M_{Z} [GeV]", label);
  drawer.MyDataMcComparisonTH1F("Z_jet_dPhi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#Delta (jet, Z)", label);
  drawer.MyDataMcComparisonTH1F("muon1_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{muon1}", label);
  drawer.MyDataMcComparisonTH1F("muon1_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#phi^{muon1}", label);
  drawer.MyDataMcComparisonTH1F("muon1_pT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{muon1}", label);
  drawer.MyDataMcComparisonTH1F("muon2_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{muon2}", label);
  drawer.MyDataMcComparisonTH1F("muon2_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#phi^{muon2}", label);
  drawer.MyDataMcComparisonTH1F("muon2_pT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{muon2}", label);
  drawer.MyDataMcComparisonTH1F("jet_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{jet}", label);
  drawer.MyDataMcComparisonTH1F("jet_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#phi^{jet}", label);
  drawer.MyDataMcComparisonTH1F("jet_pt", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{jet}", label);
  drawer.DrawFlavorComposition("jet_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{jet}", label);
  drawer.DrawFlavorComposition("jet_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#phi^{jet}", label);
  drawer.DrawFlavorComposition("jet_pt", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{jet}", label);
  drawer.MyDataMcComparisonTH1F("averageInteractionsPerCrossing", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#mu (preselection) ", label);
  drawer.MyDataMcComparisonTH1F("1st_jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "leading jet #eta (preselection)", label);
  drawer.MyDataMcComparisonTH1F("1st_jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "leading jet p_{T} (preselection)", label);
  drawer.MyDataMcComparisonTH1F("nJets_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "N_{jets} (preselection)", label);
  drawer.MyDataMcComparisonTH1F("jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{jet} (preselection)", label);
  drawer.MyDataMcComparisonTH1F("jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{jet} (preselection)", label);
  drawer.DrawFlavorComposition("jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{jet} (preselection)", label);
  drawer.DrawFlavorComposition("jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{jet} (preselection)", label);  
  drawer.MyDataMcComparisonTH1F("ZpT_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "p_{T}^{Z} [GeV]", label);
  drawer.MyDataMcComparisonTH1F("Zeta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#eta^{Z}", label);
  drawer.MyDataMcComparisonTH1F("Zphi_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "#phi^{Z}", label);
  drawer.MyDataMcComparisonTH1F("ZM_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer.GetLuminosity()), "M_{Z} [GeV]", label);
  bool DrawBalanceHistograms = true;
  if (DrawBalanceHistograms) {
    TFile* inputFile = ZJetBalance::DH_GBHOut::GetTFile("all_histData50ns.root");
    TIter next(inputFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
      TString keyName = key->GetName();
      if (not keyName.BeginsWith("DB_")) { continue; }
      if (keyName.EndsWith("_b") or keyName.EndsWith("_c") or keyName.EndsWith("_l")) { continue; }
      std::string title = ((TH1F*)key->ReadObj())->GetTitle();
      
      drawer.MyDataMcComparisonTH1F(keyName.Data(), Form("%s", title.c_str()), "p_{T}^{jet}/p_{T, Z}^{ref}", label, "H", false, -1, -1, true, 0, 5); // set xrange by hand
      drawer.DrawFlavorComposition(keyName.Data(), Form("%s", title.c_str()), "p_{T}^{jet}/p_{T, Z}^{ref}", label, "H", false, -1, -1, true, 0, 5); // set xrange by hand
    }
  }
  
  bool DrawZJetBalancePlotterOut = true;
  if (DrawZJetBalancePlotterOut) {
    ZJetBalance::DH_GBHOut drawer2("output2");
    
    drawer2.SetDataFileName("ZJetBalancePlotterOutallData.root");
    drawer2.SetLuminosity(luminosity);
    drawer2.AddMC("ZJetBalancePlotterOutallMC.root", 
		  "MC total", "MCtotal", kBlue-3, 8, true, false); // no normalization (i.e. efficiency etc)
    
    TFile* inputFile = ZJetBalance::DH_GBHOut::GetTFile("ZJetBalancePlotterOutallData.root");
    TIter next(inputFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
      std::string keyName = key->GetName();
      
      std::size_t found = keyName.find("mean");
      if (found==std::string::npos) { continue; }
      std::string xtitle = ((TH1F*)key->ReadObj())->GetXaxis()->GetTitle();
      std::string title = ((TH1F*)key->ReadObj())->GetTitle();
      drawer2.MyDataMcComparisonTH1F_GraphStyle(keyName, Form("%s", title.c_str()), xtitle, label, true, 0.5, 1.5); // set Y range by hand
      //drawer2.MyDataMcComparisonTH1F(keyName, Form("%s", title.c_str()), xtitle);      
    }//over Keys
  }
  
  // dump latex table input
  drawer.DumpCutFlow("DiMuon");
  
}
