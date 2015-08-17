// borrowing JES_ResponseFitter/util/FitZjetDB.cxx

#include <JES_ResponseFitter/JES_BalanceFitter.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TList.h>

using namespace std;

void GetBinningInformation(TH1F* h, Double_t* binEdges, Int_t& nbins)
{
  nbins = h->GetNbinsX();
  for (int iBin=1; iBin<=nbins+1; iBin++) {
    binEdges[iBin-1] = h->GetBinLowEdge(iBin);
  }
}

void DrawHisto(TH1F* h, TString ytit, double min, double max, TString xtit) {
  h->SetTitleSize(0);
  h->GetYaxis()->SetRangeUser(min,max); 
  h->SetXTitle(xtit); 
  h->SetYTitle(ytit);
  h->SetMarkerStyle(20); 
  h->SetMarkerSize(0.8); 
  h->SetLineColor(kBlack);
  h->GetXaxis()->SetMoreLogLabels(); 
  h->SetStats(0); 
  h->Draw();
}

int main(int argc, char* argv[]) {
  
  TString inFile;
  TString outTag("output");
  Bool_t  usePoisson=true;
  TList outputList;
  
  /////////// Retrieve job arguments //////////////////////////
  std::vector< std::string> options;
  for(int ii=1; ii < argc; ++ii){
    options.push_back( std::string(argv[ii]) );
  }
  
  if (argc > 1 && options.at(0).compare("-h") == 0) {
    std::cout << std::endl
	      << " job submission" << std::endl
	      << std::endl
	      << " Optional arguments:" << std::endl
	      << "  -h               Prints this menu" << std::endl
	      << "  -usePoisson      Use poisson note : without arguement (default gausian) " << std::endl
	      << "  -inFile          file" << std::endl
	      << "  -outTag          tag for output files" << std::endl
	      << std::endl;
    exit(EXIT_SUCCESS);
  }
  
  int iArg = 0;
  while(iArg < argc-1) {
    
    if (options.at(iArg).compare("-h") == 0) {
      // Ignore if not first argument
      ++iArg;
    } else if (options.at(iArg).compare("-usePoisson")==0) {
      usePoisson=true;
      ++iArg;
    } else if (options.at(iArg).compare("-outTag") == 0) {
      outTag=options.at(iArg+1);
      iArg += 2;
    } else if (options.at(iArg).compare("-inFile") == 0) {
      char tmpChar = options.at(iArg+1)[0];
      if (iArg+1 == argc || tmpChar == '-' ) {
	std::cout << " -inFile should be followed by a file or folder" << std::endl;
	return 1;
      } else {
	inFile = options.at(iArg+1);
	iArg += 2;
      }
    }else{
      std::cout << "Couldn't understand argument " << options.at(iArg) << std::endl;
      return 1;
    }
  }//while arguments
  
  TString jetDesc="Anti k_{t} #it{R} = 0.4, EM+JES";
  TString fitDesc= usePoisson ? "Modified Poisson" : "Gaussian fit";
  TString pdf= usePoisson ? 
    Form("Zjet_DB_Poisson_fits_%s.pdf", outTag.Data()) :
    Form("Zjet_DB_Gauss_fits_%s.pdf", outTag.Data());
  
  gErrorIgnoreLevel=2000; // removes Canvas print statements
  TFile *f = TFile::Open(inFile.Data());
  if (!f) {
    std::cout << "inFile=[" << inFile.Data() << "] is not valid file name." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // set cutflow histogram
  outputList.Add((TH1F*)f->Get("cutflow"));  
  outputList.Add((TH1F*)f->Get("cutflow_weighted"));  
  outputList.Add((TH1F*)f->Get("cutflow_weighted_final"));
  
  Double_t* ptbins  = new Double_t[BUFSIZ];
  Double_t* etabins = new Double_t[BUFSIZ];
  
  int       nptbins  = -1;
  int       netabins = -1;
  
  TH1F* hRef_pt = (TH1F*)f->Get("pt_binning_info");
  GetBinningInformation(hRef_pt, ptbins, nptbins);
  
  TH1F* hRef_eta = (TH1F*)f->Get("eta_binning_info");
  GetBinningInformation(hRef_eta, etabins, netabins);
  
  
  double NsigmaForFit = 1.8;
  JES_BalanceFitter *myFitter = new JES_BalanceFitter(NsigmaForFit);
  myFitter->SetPoisson();
  
  for (int iEtaBin=1; iEtaBin<netabins+1;++iEtaBin) {
    TH1F *h_mean  = new TH1F(Form("mean_pt_eta%d", iEtaBin), 
			     Form("%.1f<#eta<%.1f", etabins[iEtaBin-1], etabins[iEtaBin]),
			     nptbins, ptbins);
    TH1F *h_width = new TH1F(Form("width_pt_eta%d", iEtaBin),
			     Form("%.1f<#eta<%.1f", etabins[iEtaBin-1], etabins[iEtaBin]),
			     nptbins, ptbins);
    //TH1F *h_chi2  = new TH1F(Form("chi2_pt_eta%d", iEtaBin), "", nptbins, ptbins);
    h_mean->GetXaxis()->SetTitle("p_{T, Z}^{ref} [GeV]");
    h_width->GetXaxis()->SetTitle("p_{T, Z}^{ref} [GeV]");
    
    Info("main()", "iEtaBin=%d : histograms %s %s created", 
	 iEtaBin, h_mean->GetName(), h_width->GetName());
    outputList.Add(h_mean);
    outputList.Add(h_width);
    //outputList.Add(h_chi2);
  }
  
  for (int iPtBin=1; iPtBin<nptbins+1;++iPtBin) {
    TH1F *h_mean  = new TH1F(Form("mean_eta_pt%d", iPtBin),
			     Form("%.1f<p_{T}^{ref}<%.1f", ptbins[iPtBin-1], ptbins[iPtBin]),
			     netabins, etabins);
    TH1F *h_width = new TH1F(Form("width_eta_pt%d", iPtBin),
			     Form("%.1f<p_{T}^{ref}<%.1f", ptbins[iPtBin-1], ptbins[iPtBin]),
			     netabins, etabins);
    //TH1F *h_chi2  = new TH1F(Form("chi2_eta_pt%d", iPtBin), "", netabins, etabins);
    h_mean->GetXaxis()->SetTitle("jet #eta");
    h_width->GetXaxis()->SetTitle("jet #eta");
    
    Info("main()", "iPtBin=%d : histograms %s %s created", 
	 iPtBin, h_mean->GetName(), h_width->GetName());
    outputList.Add(h_mean);
    outputList.Add(h_width);
    //outputList.Add(h_chi2);
  }
  
  TH1F *h_fitQuants[5], *h_quantiles[5];
  for (int i=0;i<5;++i) {
    h_fitQuants[i] = new TH1F(Form("fitQuant%d",i), "", nptbins, ptbins);
    h_quantiles[i] = new TH1F(Form("quantile%d",i), "", nptbins, ptbins);
  }
  
  TCanvas *can = new TCanvas();
  can->SetMargin(0.12,0.04,0.12,0.04);
  
  can->Print(pdf+"[");
  for (int iEtaBin=1;iEtaBin<netabins+1;++iEtaBin) {
    for (int iPtBin=1;iPtBin<nptbins+1;++iPtBin) {
      double ptlow=ptbins[iPtBin]; //, ptup=ptbins[bin+1], pt=(ptlow+ptup)/2;
      double fitMin = 17.0/ptlow;
      TH1F *h = (TH1F*)f->Get(Form("DB_RefEtaBin%d_PtBin%d", iEtaBin, iPtBin));
      h->SetStats(0);
      myFitter->FitAndDraw(h,fitMin);
      can->Print(pdf);
      
      TH1F *h_mean_pt  = (TH1F*)gROOT->FindObject(Form("mean_pt_eta%d", iEtaBin));
      TH1F *h_width_pt = (TH1F*)gROOT->FindObject(Form("width_pt_eta%d", iEtaBin));
      //TH1F *h_chi2_pt  = (TH1F*)gROOT->FindObject(Form("chi2_pt_eta%d", iEtaBin));
      TH1F *h_mean_eta  = (TH1F*)gROOT->FindObject(Form("mean_eta_pt%d", iPtBin));
      TH1F *h_width_eta = (TH1F*)gROOT->FindObject(Form("width_eta_pt%d", iPtBin));
      //TH1F *h_chi2_eta  = (TH1F*)gROOT->FindObject(Form("chi2_eta_pt%d", iPtBin));
      
      h_mean_pt->SetBinContent(iPtBin,myFitter->GetMean());
      h_mean_pt->SetBinError(iPtBin,myFitter->GetMeanError());
      h_width_pt->SetBinContent(iPtBin,myFitter->GetSigma()/myFitter->GetMean());
      h_width_pt->SetBinError(iPtBin,myFitter->GetSigmaError()/myFitter->GetMean());
      //h_chi2_pt->SetBinContent(iPtBin,myFitter->GetChi2Ndof());
      
      h_mean_eta->SetBinContent(iEtaBin,myFitter->GetMean());
      h_mean_eta->SetBinError(iEtaBin,myFitter->GetMeanError());
      h_width_eta->SetBinContent(iEtaBin,myFitter->GetSigma()/myFitter->GetMean());
      h_width_eta->SetBinError(iEtaBin,myFitter->GetSigmaError()/myFitter->GetMean());
      //h_chi2_eta->SetBinContent(iEtaBin,myFitter->GetChi2Ndof());
      
      // Quantiles
      h_fitQuants[0]->SetBinContent(iPtBin,myFitter->GetMedian());
      h_fitQuants[1]->SetBinContent(iPtBin,myFitter->GetNeg2SigQuantile());
      h_fitQuants[2]->SetBinContent(iPtBin,myFitter->GetNeg1SigQuantile());
      h_fitQuants[3]->SetBinContent(iPtBin,myFitter->GetPos1SigQuantile());
      h_fitQuants[4]->SetBinContent(iPtBin,myFitter->GetPos2SigQuantile());
      
      h_quantiles[0]->SetBinContent(iPtBin,myFitter->GetHistoMedian());
      h_quantiles[1]->SetBinContent(iPtBin,myFitter->GetNeg2SigHistoQuantile());
      h_quantiles[2]->SetBinContent(iPtBin,myFitter->GetNeg1SigHistoQuantile());
      h_quantiles[3]->SetBinContent(iPtBin,myFitter->GetPos1SigHistoQuantile());
      h_quantiles[4]->SetBinContent(iPtBin,myFitter->GetPos2SigHistoQuantile());
      h_quantiles[0]->SetBinError(iPtBin,myFitter->GetHistoQuantileError());
      h_quantiles[1]->SetBinError(iPtBin,myFitter->GetHistoQuantileError());
      h_quantiles[2]->SetBinError(iPtBin,myFitter->GetHistoQuantileError());
      h_quantiles[3]->SetBinError(iPtBin,myFitter->GetHistoQuantileError());
      h_quantiles[4]->SetBinError(iPtBin,myFitter->GetHistoQuantileError());
    }
  }
  
  can->SetLogx(1);
  for (int iEtaBin=1;iEtaBin<netabins+1;++iEtaBin) {
    TH1F *h_mean_pt  = (TH1F*)gROOT->FindObject(Form("mean_pt_eta%d", iEtaBin));
    TH1F *h_width_pt = (TH1F*)gROOT->FindObject(Form("width_pt_eta%d", iEtaBin));
    //TH1F *h_chi2_pt  = (TH1F*)gROOT->FindObject(Form("chi2_pt_eta%d", iEtaBin));
    
    TLatex title(0.15, 0.2, Form("#it{ATLAS} internal %1.f<#eta<%1.f", 
				 etabins[iEtaBin-1], etabins[iEtaBin]
				 ));
    title.SetNDC();
    DrawHisto(h_mean_pt,"Mean of fit",0.8,1.2,"p_{T, Z}^{ref} [GeV]");
    myFitter->ResetTextCounters();
    myFitter->DrawTextLeft(jetDesc); 
    myFitter->DrawTextLeft(fitDesc);
    myFitter->DrawTextLeft(h_mean_pt->GetTitle());
    title.Draw();
    can->Print(pdf);
    
    DrawHisto(h_width_pt,"Width/Mean of fit",0.0,0.6,"p_{T, Z}^{ref} [GeV]");
    myFitter->ResetTextCounters();
    myFitter->DrawTextLeft(jetDesc);
    myFitter->DrawTextLeft(fitDesc);
    myFitter->DrawTextLeft(h_width_pt->GetTitle());
    title.Draw();
    can->Print(pdf);
      
    // DrawHisto(h_chi2_pt,"#it{#chi}^{2}/#it{n}_{dof} of fit",0,5);
    // myFitter->ResetTextCounters();
    // myFitter->DrawTextLeft(jetDesc); myFitter->DrawTextLeft(fitDesc);
    // title.Draw();
    // can->Print(pdf);
    
    DrawHisto(h_quantiles[0],"0#sigma, #pm1#sigma and #pm2#sigma quantiles",0,2,"p_{T, Z}^{ref} [GeV]");
    for (int i=1;i<5;++i)
      h_quantiles[i]->Draw("same");
    for (int i=0;i<5;++i) {
      h_fitQuants[i]->SetLineColor(kRed+1);
      h_fitQuants[i]->Draw("same l");
    }
    myFitter->ResetTextCounters();
    myFitter->DrawTextRight(jetDesc);
    myFitter->DrawTextRight(fitDesc,kRed+1);
    title.Draw();
    can->Print(pdf);
    
    // can->SetLogx(0);
    // gRandom->SetSeed(12345678);
    // TH1F *h = new TH1F("","",10000,-1,3);
    // h->Sumw2(); h->SetStats(0);
    // // J0
    // for (int i=0;i<100;++i) h->Fill(gRandom->Gaus(1,0.4),100); // J0
    // for (int i=0;i<100;++i) h->Fill(gRandom->Gaus(1,0.4),1); // J1
    // for (int i=0;i<100;++i) h->Fill(gRandom->Gaus(1,0.4),0.01); // J2
    // h->Draw();
    // TF1 *ff = new TF1("","gaus",-2,2);
    // h->Fit(ff,"RWL");
    // can->Print(pdf);
    
    // myFitter->FitAndDraw(h);
    // can->Print(pdf);
  }
  
  can->SetLogx(0);
  for (int iPtBin=1;iPtBin<nptbins+1;++iPtBin) {
    TH1F *h_mean_eta  = (TH1F*)gROOT->FindObject(Form("mean_eta_pt%d", iPtBin));
    TH1F *h_width_eta = (TH1F*)gROOT->FindObject(Form("width_eta_pt%d", iPtBin));
    //TH1F *h_chi2_eta  = (TH1F*)gROOT->FindObject(Form("chi2_eta_pt%d", iPtBin));
    
    TLatex title(0.15, 0.2, Form("#it{ATLAS} internal %1.f<p_{T}<%1.f", 
				 ptbins[iPtBin-1], ptbins[iPtBin]
				 ));
    title.SetNDC();
    
    DrawHisto(h_mean_eta,"Mean of fit",0.8,1.2,"jet #eta");
    myFitter->ResetTextCounters();
    myFitter->DrawTextLeft(jetDesc);
    myFitter->DrawTextLeft(fitDesc);
    myFitter->DrawTextLeft(h_mean_eta->GetTitle());
    title.Draw();
    can->Print(pdf);
    
    DrawHisto(h_width_eta,"Width/Mean of fit",0.0,0.6,"jet #eta");
    myFitter->ResetTextCounters();
    myFitter->DrawTextLeft(jetDesc);
    myFitter->DrawTextLeft(fitDesc);
    myFitter->DrawTextLeft(h_width_eta->GetTitle());
    title.Draw();
    can->Print(pdf);
    
    // DrawHisto(h_chi2_eta,"#it{#chi}^{2}/#it{n}_{dof} of fit",0,5);
    // myFitter->ResetTextCounters();
    // myFitter->DrawTextLeft(jetDesc); myFitter->DrawTextLeft(fitDesc);
    // title.Draw();
    // can->Print(pdf);
  }
  
  can->Print(pdf+"]");
  
  printf("\nProduced:\n  %s\n\n",pdf.Data());
  
  // write to the output
  TFile fileOutput(Form("ZJetBalancePlotterOut%s.root", outTag.Data()), "RECREATE");
  outputList.Print();
  
  TIter next(&outputList);
  TObject* object = 0;
  while ((object = next())) {
    object->Write();
  }//over Keys  
  
  fileOutput.Write();
}
