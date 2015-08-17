#include <ZJetBalance/DH_GBHOut.h>

#include <string>

int main()
{
  const double luminosity_muon=0.0783190;
  const double luminosity_electron=0.0620279;
  
  const std::string label="Internal";
  
  // ========================================
  ZJetBalance::DH_GBHOut drawer_electron("output_electron");
  
  drawer_electron.SetDataFileName("all_el_histData50ns.root");
  drawer_electron.SetLuminosity(luminosity_electron);
  drawer_electron.AddMC("all_el_hist361106.root",
			"Z#rightarrow ee",
			"Zee",
			kYellow, 8, true);
  drawer_electron.AddMC("all_el_hist361108.root",
			"Z#rightarrow#tau#tau",
			"Ztautau",
			kYellow+3, 8, true);
  
  drawer_electron.AddMC("all_el_hist410004.root",
			"t#bar{t}",
			"Tt",
			kGreen, 8, true);
  
  const std::string period("C2-C4");
  drawer_electron.MyDataMcComparisonTH1F("ZpT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "p_{T}^{Z} [GeV]", label);
  drawer_electron.MyDataMcComparisonTH1F("ZM", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "M_{Z} [GeV]", label);
  drawer_electron.MyDataMcComparisonTH1F("Z_jet_dPhi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#Delta (jet, Z)", label);
  drawer_electron.MyDataMcComparisonTH1F("electron1_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#eta^{electron1}", label);
  drawer_electron.MyDataMcComparisonTH1F("electron1_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#phi^{electron1}", label);
  drawer_electron.MyDataMcComparisonTH1F("electron1_pT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "p_{T}^{electron1}", label);
  drawer_electron.MyDataMcComparisonTH1F("electron2_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#eta^{electron2}", label);
  drawer_electron.MyDataMcComparisonTH1F("electron2_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#phi^{electron2}", label);
  drawer_electron.MyDataMcComparisonTH1F("electron2_pT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "p_{T}^{electron2}", label);
  drawer_electron.MyDataMcComparisonTH1F("jet_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#eta^{jet}", label);
  drawer_electron.MyDataMcComparisonTH1F("jet_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#phi^{jet}", label);
  drawer_electron.MyDataMcComparisonTH1F("jet_pt", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "p_{T}^{jet}", label);
  drawer_electron.DrawFlavorComposition("jet_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#eta^{jet}", label);
  drawer_electron.DrawFlavorComposition("jet_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#phi^{jet}", label);
  drawer_electron.DrawFlavorComposition("jet_pt", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "p_{T}^{jet}", label);
  drawer_electron.MyDataMcComparisonTH1F("averageInteractionsPerCrossing", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#mu (preselection) ", label);
  drawer_electron.MyDataMcComparisonTH1F("1st_jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "leading jet #eta (preselection)", label);
  drawer_electron.MyDataMcComparisonTH1F("1st_jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "leading jet p_{T} (preselection)", label);
  drawer_electron.MyDataMcComparisonTH1F("nJets_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "N_{jets} (preselection)", label);
  drawer_electron.MyDataMcComparisonTH1F("jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#eta^{jet} (preselection)", label);
  drawer_electron.MyDataMcComparisonTH1F("jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "p_{T}^{jet} (preselection)", label);
  drawer_electron.DrawFlavorComposition("jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#eta^{jet} (preselection)", label);
  drawer_electron.DrawFlavorComposition("jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "p_{T}^{jet} (preselection)", label);  
  drawer_electron.MyDataMcComparisonTH1F("ZpT_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "p_{T}^{Z} [GeV]", label);
  drawer_electron.MyDataMcComparisonTH1F("Zeta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#eta^{Z}", label);
  drawer_electron.MyDataMcComparisonTH1F("Zphi_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "#phi^{Z}", label);
  drawer_electron.MyDataMcComparisonTH1F("ZM_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_electron.GetLuminosity()), "M_{Z} [GeV]", label);
  bool DrawBalanceHistograms = true;
  if (DrawBalanceHistograms) {
    TFile* inputFile = ZJetBalance::DH_GBHOut::GetTFile("all_el_histData50ns.root");
    TIter next(inputFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
      TString keyName = key->GetName();
      if (not keyName.BeginsWith("DB_")) { continue; }
      if (keyName.EndsWith("_b") or keyName.EndsWith("_c") or keyName.EndsWith("_l")) { continue; }
      std::string title = ((TH1F*)key->ReadObj())->GetTitle();
      
      drawer_electron.MyDataMcComparisonTH1F(keyName.Data(), Form("%s", title.c_str()), "p_{T}^{jet}/p_{T, Z}^{ref}", label, "H", false, -1, -1, true, 0, 2.5); // set xrange by hand
      drawer_electron.DrawFlavorComposition(keyName.Data(), Form("%s", title.c_str()), "p_{T}^{jet}/p_{T, Z}^{ref}", label, "H", false, -1, -1, true, 0, 2.5); // set xrange by hand
    }
  }
  
  bool DrawZJetBalancePlotterOut = true;
  if (DrawZJetBalancePlotterOut) {
    ZJetBalance::DH_GBHOut drawer2_electron("output2_electron");
    
    drawer2_electron.SetDataFileName("ZJetBalancePlotterOutallData_el.root");
    drawer2_electron.SetLuminosity(luminosity_electron);
    drawer2_electron.AddMC("ZJetBalancePlotterOutallMC_el.root", 
			   "MC total", "MCtotal", kBlue-3, 8, true, false); // no normalization (i.e. efficiency etc)
    
    TFile* inputFile = ZJetBalance::DH_GBHOut::GetTFile("ZJetBalancePlotterOutallData_el.root");
    TIter next(inputFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
      std::string keyName = key->GetName();
      
      std::size_t found = keyName.find("mean");
      if (found==std::string::npos) { continue; }
      std::string xtitle = ((TH1F*)key->ReadObj())->GetXaxis()->GetTitle();
      std::string title = ((TH1F*)key->ReadObj())->GetTitle();
      drawer2_electron.MyDataMcComparisonTH1F_GraphStyle(keyName, Form("%s", title.c_str()), xtitle, label, true, 0.5, 1.5); // set Y range by hand
      //drawer2_electron.MyDataMcComparisonTH1F(keyName, Form("%s", title.c_str()), xtitle);      
    }//over Keys
  }


  // ========================================
  ZJetBalance::DH_GBHOut drawer_muon("output_muon");
  
  drawer_muon.SetDataFileName("all_mu_histData50ns.root");
  drawer_muon.SetLuminosity(luminosity_muon);
  drawer_muon.AddMC("all_mu_hist361107.root",
		    "Z#rightarrow#mu#mu",
		    "Zmumu",
		    kYellow, 8, true);
  drawer_muon.AddMC("all_mu_hist361108.root",
		    "Z#rightarrow#tau#tau",
		    "Ztautau",
		    kYellow+3, 8, true);
  
  drawer_muon.AddMC("all_mu_hist410004.root",
		    "t#bar{t}",
		    "Tt",
		    kGreen, 8, true);
  
  drawer_muon.MyDataMcComparisonTH1F("ZpT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "p_{T}^{Z} [GeV]", label);
  drawer_muon.MyDataMcComparisonTH1F("ZM", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "M_{Z} [GeV]", label);
  drawer_muon.MyDataMcComparisonTH1F("Z_jet_dPhi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#Delta (jet, Z)", label);
  drawer_muon.MyDataMcComparisonTH1F("muon1_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#eta^{muon1}", label);
  drawer_muon.MyDataMcComparisonTH1F("muon1_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#phi^{muon1}", label);
  drawer_muon.MyDataMcComparisonTH1F("muon1_pT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "p_{T}^{muon1}", label);
  drawer_muon.MyDataMcComparisonTH1F("muon2_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#eta^{muon2}", label);
  drawer_muon.MyDataMcComparisonTH1F("muon2_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#phi^{muon2}", label);
  drawer_muon.MyDataMcComparisonTH1F("muon2_pT", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "p_{T}^{muon2}", label);
  drawer_muon.MyDataMcComparisonTH1F("jet_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#eta^{jet}", label);
  drawer_muon.MyDataMcComparisonTH1F("jet_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#phi^{jet}", label);
  drawer_muon.MyDataMcComparisonTH1F("jet_pt", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "p_{T}^{jet}", label);
  drawer_muon.DrawFlavorComposition("jet_eta", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#eta^{jet}", label);
  drawer_muon.DrawFlavorComposition("jet_phi", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#phi^{jet}", label);
  drawer_muon.DrawFlavorComposition("jet_pt", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "p_{T}^{jet}", label);
  drawer_muon.MyDataMcComparisonTH1F("averageInteractionsPerCrossing", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#mu (preselection) ", label);
  drawer_muon.MyDataMcComparisonTH1F("1st_jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "leading jet #eta (preselection)", label);
  drawer_muon.MyDataMcComparisonTH1F("1st_jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "leading jet p_{T} (preselection)", label);
  drawer_muon.MyDataMcComparisonTH1F("nJets_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "N_{jets} (preselection)", label);
  drawer_muon.MyDataMcComparisonTH1F("jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#eta^{jet} (preselection)", label);
  drawer_muon.MyDataMcComparisonTH1F("jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "p_{T}^{jet} (preselection)", label);
  drawer_muon.DrawFlavorComposition("jet_eta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#eta^{jet} (preselection)", label);
  drawer_muon.DrawFlavorComposition("jet_pt_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "p_{T}^{jet} (preselection)", label);  
  drawer_muon.MyDataMcComparisonTH1F("ZpT_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "p_{T}^{Z} [GeV]", label);
  drawer_muon.MyDataMcComparisonTH1F("Zeta_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#eta^{Z}", label);
  drawer_muon.MyDataMcComparisonTH1F("Zphi_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "#phi^{Z}", label);
  drawer_muon.MyDataMcComparisonTH1F("ZM_beforecut", Form("%s   (%.2f fb^{-1})", period.c_str(), drawer_muon.GetLuminosity()), "M_{Z} [GeV]", label);

  if (DrawBalanceHistograms) {
    TFile* inputFile = ZJetBalance::DH_GBHOut::GetTFile("all_mu_histData50ns.root");
    TIter next(inputFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
      TString keyName = key->GetName();
      if (not keyName.BeginsWith("DB_")) { continue; }
      if (keyName.EndsWith("_b") or keyName.EndsWith("_c") or keyName.EndsWith("_l")) { continue; }
      std::string title = ((TH1F*)key->ReadObj())->GetTitle();
      
      drawer_muon.MyDataMcComparisonTH1F(keyName.Data(), Form("%s", title.c_str()), "p_{T}^{jet}/p_{T, Z}^{ref}", label, "H", false, -1, -1, true, 0, 2.5); // set xrange by hand
      drawer_muon.DrawFlavorComposition(keyName.Data(), Form("%s", title.c_str()), "p_{T}^{jet}/p_{T, Z}^{ref}", label, "H", false, -1, -1, true, 0, 2.5); // set xrange by hand
    }
  }
  
  if (DrawZJetBalancePlotterOut) {
    ZJetBalance::DH_GBHOut drawer2_muon("output2_muon");
    
    drawer2_muon.SetDataFileName("ZJetBalancePlotterOutallData_mu.root");
    drawer2_muon.SetLuminosity(luminosity_muon);
    drawer2_muon.AddMC("ZJetBalancePlotterOutallMC_mu.root", 
		       "MC total", "MCtotal", kBlue-3, 8, true, false); // no normalization (i.e. efficiency etc)
    
    TFile* inputFile = ZJetBalance::DH_GBHOut::GetTFile("ZJetBalancePlotterOutallData_mu.root");
    TIter next(inputFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
      std::string keyName = key->GetName();
      
      std::size_t found = keyName.find("mean");
      if (found==std::string::npos) { continue; }
      std::string xtitle = ((TH1F*)key->ReadObj())->GetXaxis()->GetTitle();
      std::string title = ((TH1F*)key->ReadObj())->GetTitle();
      drawer2_muon.MyDataMcComparisonTH1F_GraphStyle(keyName, Form("%s", title.c_str()), xtitle, label, true, 0.5, 1.5); // set Y range by hand
      //drawer2_muon.MyDataMcComparisonTH1F(keyName, Form("%s", title.c_str()), xtitle);      
    }//over Keys
  }
  
  // dump latex table input
  drawer_electron.DumpCutFlow("DiElectron");
  drawer_muon.DumpCutFlow("DiMuon");
  
}
