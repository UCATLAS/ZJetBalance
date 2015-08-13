// to do list 
// re-write in a better way...

#include <TObject.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TROOT.h>
#include <TTree.h>
#include <TF1.h>
#include <TKey.h>

#include <string>
#include <iostream>

#include <string.h>
#include <stdlib.h>

// ======================================
TFile* GetTFile(const std::string& filename)
{
  TFile* rc = TFile::Open(filename.c_str());
  if (!rc) {
    std::cout << "ERROR> give file name = " << filename << " does not exist." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  return rc;
}

// ======================================
TObject* GetObject(TFile* f, const std::string& name)
{
  TObject* rc = f->Get(name.c_str());
  if (!rc) {
    std::cout << "ERROR> " << name << " does not exist." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  return rc;
}

// ======================================
void CustimizeCampus(TCanvas& c1)
{
  c1.SetRightMargin(0.25);
  c1.SetLeftMargin(0.15);
  c1.SetBottomMargin(0.15);
  c1.SetTopMargin(0.10);
  c1.SetGridy(1);
  c1.SetGridx(1);
}

// ======================================
void MyDrawTH1FForTruthFlavors(TCanvas& c1, 
			       const std::string& outputfile,
			       TFile* fMC,
			       TFile* fData,
			       const std::string& histoname,
			       const std::string& comment,
			       const std::string& xtitle,
			       const double& weight
			       )
{
  TH1F* hMC_b_orig = (TH1F*)GetObject(fMC, histoname+"_b");
  TH1F* hMC_c_orig = (TH1F*)GetObject(fMC, histoname+"_c");
  TH1F* hMC_l_orig = (TH1F*)GetObject(fMC, histoname+"_l");
  TH1F* hData_orig = (TH1F*)GetObject(fData, histoname);
  
  // copy in order to avoid scaling original one
  TH1F hMC_b, hMC_b2, hMC_c, hMC_l, hData;
  hMC_b_orig->Copy(hMC_b);
  hMC_c_orig->Copy(hMC_c);
  hMC_l_orig->Copy(hMC_l);
  hData_orig->Copy(hData);
  
  // Scale MC
  hData.SetMinimum(hData.GetMinimum()<0 ? hData.GetMinimum() : 0.);
  hMC_b.SetMinimum(hMC_b.GetMinimum()<0 ? hMC_b.GetMinimum() : 0.);
  hMC_c.SetMinimum(hMC_c.GetMinimum()<0 ? hMC_c.GetMinimum() : 0.);
  hMC_l.SetMinimum(hMC_l.GetMinimum()<0 ? hMC_l.GetMinimum() : 0.);
  
  hData.SetTitle("");
  hMC_b.SetTitle("");
  hMC_c.SetTitle("");
  hMC_l.SetTitle("");

  hMC_b.Scale(weight);
  hMC_c.Scale(weight);
  hMC_l.Scale(weight);
  
  const double entry_b = hMC_b.Integral(-1, -1);
  const double entry_c = hMC_c.Integral(-1, -1);
  const double entry_l = hMC_l.Integral(-1, -1);
  
  hMC_b = hMC_b + hMC_c + hMC_l;
  hMC_c = hMC_c + hMC_l;
  hMC_b.Copy(hMC_b2);
  
  hMC_b.SetFillColor(kYellow);
  hMC_c.SetFillColor(kOrange);
  hMC_l.SetFillColor(kBlue);
  
  hMC_b.SetLineColor(kYellow);
  hMC_c.SetLineColor(kOrange);
  hMC_l.SetLineColor(kBlue);
  
  hMC_b2.SetLineColor(kBlack);
  
  hData.SetMarkerStyle(8);
  hData.SetLineColor(kBlack);
  hData.SetMarkerColor(kBlack);  
  
  if (hMC_b.GetMaximum()<hData.GetMaximum()) {
    hData.GetXaxis()->SetTitle(xtitle.c_str());
    hData.Draw("PE");
    hMC_b.Draw("H SAME");
    hMC_c.Draw("H SAME");
    hMC_l.Draw("H SAME");
    hMC_b2.Draw("H SAME");
    hData.Draw("PE SAME");
  } else {
    hMC_b.GetXaxis()->SetTitle(xtitle.c_str());
    hMC_b.Draw("H");
    hMC_c.Draw("H SAME");
    hMC_l.Draw("H SAME");
    hMC_b2.Draw("H SAME");
    hData.Draw("PE SAME");
  }
  
  TLegend leg(0.76, 0.15, 0.98, 0.90);
  
  leg.AddEntry(&hData, Form("#splitline{#splitline{Data}{(%.0f)}}{#splitline{mean=%.1f}{RMS=%.1f}}", hData.Integral(-1, -1), hData.GetMean(), hData.GetRMS()), "PE");
  leg.AddEntry(&hMC_b, Form("#splitline{#splitline{MC - b}{(%.0f)}}{#splitline{mean=%.1f}{RMS=%.1f}}", entry_b, hMC_b.GetMean(), hMC_b.GetRMS()), "F");
  leg.AddEntry(&hMC_c, Form("#splitline{#splitline{MC - c}{(%.0f)}}{#splitline{mean=%.1f}{RMS=%.1f}}", entry_c, hMC_c.GetMean(), hMC_c.GetRMS()), "F");
  leg.AddEntry(&hMC_l, Form("#splitline{#splitline{MC - l}{(%.0f)}}{#splitline{mean=%.1f}{RMS=%.1f}}", entry_l, hMC_l.GetMean(), hMC_l.GetRMS()), "F");
  leg.SetTextSize(0.03);
  leg.Draw();
  
  TLatex myComment(0.20, 0.94, comment.c_str());
  myComment.SetNDC();
  myComment.Draw();
  
  c1.Print(Form("%s.pdf", outputfile.c_str()));  
}

// ======================================
void MyDrawTH1F(TCanvas& c1, 
		const std::string& outputfile,
		TFile* fMC,
		TFile* fData,
		const std::string& histoname,
		const std::string& comment,
		const std::string& xtitle,
		const double& weight
		)
{
  TH1F* hMC_orig   = (TH1F*)GetObject(fMC, histoname);
  TH1F* hData_orig = (TH1F*)GetObject(fData, histoname);
  
  // copy in order to avoid scaling original one
  TH1F hMC, hData;
  hMC_orig->Copy(hMC);
  hData_orig->Copy(hData);
    
  hMC.SetTitle("");
  hData.SetTitle("");
  
  // Scale MC
  hMC.Scale(weight);
  
  hMC.SetFillColor(kYellow);
  hData.SetMarkerStyle(8);
  hData.SetLineColor(kBlack);
  hMC.SetLineColor(kBlack);
  hData.SetMarkerColor(kBlack);
  
  if (hMC.GetMaximum()<hData.GetMaximum()) {
    hData.SetMinimum(hData.GetMinimum()<0 ? hData.GetMinimum() : 0.);
    hData.GetXaxis()->SetTitle(xtitle.c_str());
    hData.Draw("PE");
    hMC.Draw("H SAME");
    hData.Draw("PE SAME");
  } else {
    hMC.SetMinimum(hMC.GetMinimum()<0 ? hMC.GetMinimum() : 0.);
    hMC.GetXaxis()->SetTitle(xtitle.c_str());
    hMC.Draw("H");
    hData.Draw("PE SAME");
  }
  
  TLegend leg(0.76, 0.15, 0.98, 0.90);
  leg.AddEntry(&hData, Form("#splitline{#splitline{Data}{(%.0f)}}{#splitline{mean=%.1f}{RMS=%.1f}}", hData.Integral(-1, -1), hData.GetMean(), hData.GetRMS()), "PE");
  leg.AddEntry(&hMC,   Form("#splitline{#splitline{MC}{(%.0f)}}{#splitline{mean=%.1f}{RMS=%.1f}}", hMC.Integral(-1, -1), hMC.GetMean(), hMC.GetRMS()), "F");
  leg.Draw();
  
  TLatex myComment(0.20, 0.94, comment.c_str());
  myComment.SetNDC();
  myComment.Draw();
  
  c1.Print(Form("%s.pdf", outputfile.c_str()));
}

// ======================================
void MyDrawResponse(TCanvas& c1, 
		    const std::string& outputfile,
		    TFile* fMC,
		    TFile* fData,
		    std::string histoname,
		    std::string comment,
		    std::string xtitle
		    )
{
  TH1F* hMC_orig   = (TH1F*)GetObject(fMC, histoname);
  TH1F* hData_orig = (TH1F*)GetObject(fData, histoname);
  
  // copy in order to avoid scaling original one
  TH1F hMC, hData;
  hMC_orig->Copy(hMC);
  hData_orig->Copy(hData);
  
  hData.GetXaxis()->SetTitleSize(0.04);
  hData.SetMarkerStyle(8);
  hData.SetLineColor(kBlack);
  hData.SetMarkerColor(kBlack);
  
  hMC.GetXaxis()->SetTitleSize(0.04);
  hMC.SetMarkerStyle(8);
  hMC.SetLineColor(0);
  hMC.SetMarkerColor(kBlue);
  hMC.SetFillColor(kBlue-10);
  
  hData.GetXaxis()->SetTitle(xtitle.c_str());
  hMC.GetXaxis()->SetTitle(xtitle.c_str());
  
  if (hMC.GetMaximum()<hData.GetMaximum()) {
    hData.Draw("PE");
    hMC.Draw("PE2 SAME");
    hData.Draw("PE SAME");
  } else {
    hMC.Draw("PE2");
    hData.Draw("PE SAME");
  }
  
  TLegend leg(0.76, 0.75, 0.98, 0.90);
  leg.AddEntry(&hData, "Data", "PE");
  leg.AddEntry(&hMC, "MC", "PEF");
  leg.Draw();
  
  TLatex myComment(0.20, 0.94, comment.c_str());
  myComment.SetNDC();
  myComment.Draw();
  
  c1.Print(Form("%s.pdf", outputfile.c_str()));
}

// ======================================
void DrawJESResponseFitterOut(TCanvas& c1, 
			      const std::string& outputfile, 
			      const std::string& data_file, 
			      const std::string& mc_file)
{
  // JRF = Jet Response Fitter out
  TFile* fMCJRF   = GetTFile(mc_file);
  TFile* fDataJRF = GetTFile(data_file);
  
  TIter next(fMCJRF->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    std::string keyName = key->GetName();
    std::string title = ((TH1F*)GetObject(fMCJRF, keyName))->GetTitle();
    std::string xtitle = ((TH1F*)GetObject(fMCJRF, keyName))->GetXaxis()->GetTitle();
    MyDrawResponse(c1, outputfile, fMCJRF, fDataJRF, keyName, "", xtitle);
  }//over Keys
  
  fMCJRF->Close();
  fDataJRF->Close();
}

// ======================================
void DrawPtBalanceForEachEtaPtBin(TCanvas& c1, 
				  const std::string& outputfile, 
				  const std::string& data_file, 
				  const std::string& mc_file,
				  double luminosity, double mcLuminosityWeight)
{
  TFile* fMC   = GetTFile(mc_file);
  TFile* fData = GetTFile(data_file);
  
  TIter next(fMC->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    std::string keyName = key->GetName();
    std::size_t found = keyName.find("DB_Ref");
    bool foundBalacePlot = (found!=std::string::npos);
    
    if (!foundBalacePlot) continue;
    std::string title = ((TH1F*)GetObject(fMC, keyName))->GetTitle();
    MyDrawTH1F(c1, outputfile, fMC, fData, keyName, Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), 
	       "Jet p_{T} / p_{T}^{ref}", mcLuminosityWeight);
  }//over Keys
  
  fMC->Close();
  fData->Close();
}

// ======================================
void DrawValidationPlots(TCanvas& c1, 
			 const std::string& outputfile,
			 const std::string& data_file,
			 const std::string& mc_file,
			 const double& luminosity, 
			 const double& mcLuminosityWeight)
{
  TFile* fMC   = GetTFile(mc_file);
  TFile* fData = GetTFile(data_file);
  
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalanceZpT", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "p_{T}^{Z} [GeV]", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalanceZM", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "M_{Z} Z [GeV]", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalanceZ_jet_dPhi", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#Delta (jet, Z)", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancemuon1_eta", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#eta^{muon1}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancemuon1_phi", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#phi^{muon1}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancemuon1_pT", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "p_{T}^{muon1}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancemuon2_eta", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#eta^{muon2}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancemuon2_phi", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#phi^{muon2}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancemuon2_pT", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "p_{T}^{muon2}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalanceelectron1_eta", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#eta^{electron1}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalanceelectron1_phi", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#phi^{electron1}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalanceelectron1_pT", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "p_{T}^{electron1}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalanceelectron2_eta", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#eta^{electron2}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalanceelectron2_phi", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#phi^{electron2}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalanceelectron2_pT", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "p_{T}^{electron2}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancejet_eta", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#eta^{jet}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancejet_phi", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#phi^{jet}", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancejet_pt", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "p_{T}^{jet}", mcLuminosityWeight);
  
  MyDrawTH1FForTruthFlavors(c1, outputfile, fMC, fData, "ProcZJetBalancejet_eta", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#eta^{jet}", mcLuminosityWeight);
  MyDrawTH1FForTruthFlavors(c1, outputfile, fMC, fData, "ProcZJetBalancejet_phi", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#phi^{jet}", mcLuminosityWeight);
  MyDrawTH1FForTruthFlavors(c1, outputfile, fMC, fData, "ProcZJetBalancejet_pt", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "p_{T}^{jet}", mcLuminosityWeight);
  
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalanceaverageInteractionsPerCrossing", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "Average Interaction per Bunch Crossing (preselection) ", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancenJets_beforecut", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "N_{jets} (preselection)", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancejet_eta_beforecut", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#eta^{jet} (preselection)", mcLuminosityWeight);
  MyDrawTH1F(c1, outputfile, fMC, fData, "ProcZJetBalancejet_pt_beforecut", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "p_{T}^{jet} (preselection)", mcLuminosityWeight);
  
  MyDrawTH1FForTruthFlavors(c1, outputfile, fMC, fData, "ProcZJetBalancejet_eta_beforecut", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "#eta^{jet} (preselection)", mcLuminosityWeight);
  MyDrawTH1FForTruthFlavors(c1, outputfile, fMC, fData, "ProcZJetBalancejet_pt_beforecut", Form("#it{ATLAS} internal Luminosity=%.2f fb^{-1}", luminosity), "p_{T}^{jet} (preselection)", mcLuminosityWeight);
  
  
  fMC->Close();
  fData->Close();
}

void print_usage(const std::string commandname)
{
  printf("%s -d <data histogram file> -m <MC histogram file> -l <luminosity in fb-1> \n",
	 commandname.c_str());
  printf("options: \n");
  printf("-a <additional weight i.e. skiming efficiency in DxAOD (default 1.)>");
  printf("-s <draw function selection (default 1)>\n");
  printf("  1: call DrawValidationPlots \n");
  printf("  2: call DrawPtBalanceForEachEtaPtBin \n");
  printf("  3: call DrawJESResponseFitterOut \n");
  printf("-w if weighted number of events is used for normalization \n"); 
  printf("-o <output file (default drawZJetBalanceOutput)> \n");
}

// ======================================
int main(int argc, char* argv[])
{
  int result;
  int modeSelection = 1;
  float luminosity = -1000.;
  float additional_weight = 1.0;
  bool useWeightedForNormalization = false;
  std::string data_file(""), mc_file(""), outputfile("drawZJetBalanceOutput");
  
  while((result=getopt(argc,argv,"d:m:l:o:s:a:w"))!=-1){
    switch(result){  
    case 'd':
      data_file = optarg;
      break;
    case 'm':
      mc_file = optarg;
      break;
    case 'l':
      luminosity = strtof(optarg, NULL);
      break;
    case 'o':
      outputfile = optarg;
      break;
    case 'w':
      useWeightedForNormalization = true;
      break;
    case 'a':
      additional_weight = strtof(optarg, NULL);
      break;
    case 's':
      modeSelection = strtoull(optarg, NULL, 0);
      break;
    default:
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  
  // check input parameter
  if ((mc_file.empty()) or (data_file.empty()) or (luminosity<0)) {
    printf("ERROR : please give all the needed input parameters.\n");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);    
  }
  
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);  
  
  // for normalization calculation
  TFile* fMC   = GetTFile(mc_file);
  double MCTotal = useWeightedForNormalization ? 
    ((TH1F*)GetObject(fMC, "cutflow_weighted"))->GetBinContent(1) : 
    ((TH1F*)GetObject(fMC, "cutflow"))->GetBinContent(1);
  const double McLuminosityWeight = luminosity/MCTotal*additional_weight;
  
  TCanvas c1;
  CustimizeCampus(c1);
  c1.Print(Form("%s.pdf[", outputfile.c_str()));  
  
  switch (modeSelection) {
  case 1:
    DrawValidationPlots(c1, outputfile, data_file, mc_file, luminosity, McLuminosityWeight);
    break;
  case 2:
    DrawPtBalanceForEachEtaPtBin(c1, outputfile, data_file, mc_file, luminosity, McLuminosityWeight);
    break;
  case 3:
    DrawJESResponseFitterOut(c1, outputfile, data_file, mc_file);
    break;
  default:
    printf("ERROR : invalid function selection = %d \n", modeSelection);
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  c1.Print(Form("%s.pdf]", outputfile.c_str()));  
}
