// borrowing JES_ResponseFitter/util/FitZjetDB.cxx

#include <JES_ResponseFitter/JES_BalanceFitter.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom.h>

using namespace std;

void GetBinningInformation(TH1F* h, Double_t* binEdges, Int_t& nbins)
{
  nbins = h->GetNbinsX();
  for (int iBin=1; iBin<=nbins+1; iBin++) {
    binEdges[iBin-1] = h->GetBinLowEdge(iBin);
  }
}

void DrawHisto(TH1F* h, TString ytit, double min, double max) {
  h->GetYaxis()->SetRangeUser(min,max); h->SetXTitle("#it{p}_{T}^{ref} [GeV]"); h->SetYTitle(ytit);
  h->SetMarkerStyle(20); h->SetMarkerSize(0.8); h->SetLineColor(kBlack);
  h->GetXaxis()->SetMoreLogLabels(); h->SetStats(0); h->Draw();
}

int main(int argc, char* argv[]) {
  
  TString inFile;

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
	      << "  -inFile          file" << std::endl
	      << std::endl;
    exit(EXIT_SUCCESS);
  }
  
  int iArg = 0;
  while(iArg < argc-1) {
    
    if (options.at(iArg).compare("-h") == 0) {
      // Ignore if not first argument
      ++iArg;
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
  
  
  
  TString arg, fitDesc="Gaussian fit", jetDesc="Anti k_{t} #it{R} = 0.4, EM+JES";
  TString pdf="Zjet_DB_Gauss_fits.pdf";
  
  gErrorIgnoreLevel=2000; // removes Canvas print statements
  TFile *f = TFile::Open(inFile.Data());
  if (!f) {
    std::cout << "inFile=[" << inFile.Data() << "] is not valid file name." << std::endl;
    exit(EXIT_FAILURE);
  }

  Double_t* ptbins  = new Double_t[BUFSIZ];
  int       nptbins = -1;
  
  TH1F* hRef = (TH1F*)f->Get("h_pt_binning_info");
  GetBinningInformation(hRef, ptbins, nptbins);
  

  double NsigmaForFit = 1.8;
  JES_BalanceFitter *myFitter = new JES_BalanceFitter(NsigmaForFit);
  
  TH1F *h_mean  = new TH1F("mean", "", nptbins, ptbins);
  TH1F *h_width = new TH1F("width","", nptbins, ptbins);
  TH1F *h_chi2  = new TH1F("chi2", "", nptbins, ptbins);
  
  TH1F *h_fitQuants[5], *h_quantiles[5];
  for (int i=0;i<5;++i) {
    h_fitQuants[i] = new TH1F(Form("fitQuant%d",i), "", nptbins, ptbins);
    h_quantiles[i] = new TH1F(Form("quantile%d",i), "", nptbins, ptbins);
  }
  
  TCanvas *can = new TCanvas();
  can->SetMargin(0.12,0.04,0.12,0.04);
  
  can->Print(pdf+"[");
  for (int bin=1;bin<nptbins+1;++bin) {
    double ptlow=ptbins[bin]; //, ptup=ptbins[bin+1], pt=(ptlow+ptup)/2;
    double fitMin = 17.0/ptlow;
    printf("INFO> %s \n", Form("DB_RefEtaBin_PtBin%d",bin));
    TH1F *h = (TH1F*)f->Get(Form("DB_RefEtaBin_PtBin%d",bin));
    h->SetStats(0);
    myFitter->FitAndDraw(h,fitMin);
    can->Print(pdf);

    h_mean->SetBinContent(bin,myFitter->GetMean());
    h_mean->SetBinError(bin,myFitter->GetMeanError());
    h_width->SetBinContent(bin,myFitter->GetSigma()/myFitter->GetMean());
    h_width->SetBinError(bin,myFitter->GetSigmaError()/myFitter->GetMean());
    h_chi2->SetBinContent(bin,myFitter->GetChi2Ndof());
    
    // Quantiles
    h_fitQuants[0]->SetBinContent(bin,myFitter->GetMedian());
    h_fitQuants[1]->SetBinContent(bin,myFitter->GetNeg2SigQuantile());
    h_fitQuants[2]->SetBinContent(bin,myFitter->GetNeg1SigQuantile());
    h_fitQuants[3]->SetBinContent(bin,myFitter->GetPos1SigQuantile());
    h_fitQuants[4]->SetBinContent(bin,myFitter->GetPos2SigQuantile());

    h_quantiles[0]->SetBinContent(bin,myFitter->GetHistoMedian());
    h_quantiles[1]->SetBinContent(bin,myFitter->GetNeg2SigHistoQuantile());
    h_quantiles[2]->SetBinContent(bin,myFitter->GetNeg1SigHistoQuantile());
    h_quantiles[3]->SetBinContent(bin,myFitter->GetPos1SigHistoQuantile());
    h_quantiles[4]->SetBinContent(bin,myFitter->GetPos2SigHistoQuantile());
    h_quantiles[0]->SetBinError(bin,myFitter->GetHistoQuantileError());
    h_quantiles[1]->SetBinError(bin,myFitter->GetHistoQuantileError());
    h_quantiles[2]->SetBinError(bin,myFitter->GetHistoQuantileError());
    h_quantiles[3]->SetBinError(bin,myFitter->GetHistoQuantileError());
    h_quantiles[4]->SetBinError(bin,myFitter->GetHistoQuantileError());
  }
  can->SetLogx();

  DrawHisto(h_mean,"Mean of fit",0.8,1.2);
  myFitter->ResetTextCounters();
  myFitter->DrawTextLeft(jetDesc); myFitter->DrawTextLeft(fitDesc);
  can->Print(pdf);

  DrawHisto(h_width,"Width/Mean of fit",0.0,0.6);
  myFitter->ResetTextCounters();
  myFitter->DrawTextLeft(jetDesc); myFitter->DrawTextLeft(fitDesc);
  can->Print(pdf);

  DrawHisto(h_chi2,"#it{#chi}^{2}/#it{n}_{dof} of fit",0,5);
  myFitter->ResetTextCounters();
  myFitter->DrawTextLeft(jetDesc); myFitter->DrawTextLeft(fitDesc);
  can->Print(pdf);

  DrawHisto(h_quantiles[0],"0#sigma, #pm1#sigma and #pm2#sigma quantiles",0,2);
  for (int i=1;i<5;++i)
    h_quantiles[i]->Draw("same");
  for (int i=0;i<5;++i) {
    h_fitQuants[i]->SetLineColor(kRed+1);
    h_fitQuants[i]->Draw("same l");
  }
  myFitter->ResetTextCounters();
  myFitter->DrawTextRight(jetDesc); myFitter->DrawTextRight(fitDesc,kRed+1);
  can->Print(pdf);

  can->SetLogx(0);
  gRandom->SetSeed(12345678);
  TH1F *h = new TH1F("","",10000,-1,3);
  h->Sumw2(); h->SetStats(0);
  // J0
  for (int i=0;i<100;++i) h->Fill(gRandom->Gaus(1,0.4),100); // J0
  for (int i=0;i<100;++i) h->Fill(gRandom->Gaus(1,0.4),1); // J1
  for (int i=0;i<100;++i) h->Fill(gRandom->Gaus(1,0.4),0.01); // J2
  h->Draw();
  TF1 *ff = new TF1("","gaus",-2,2);
  h->Fit(ff,"RWL");
  can->Print(pdf);

  myFitter->FitAndDraw(h);
  can->Print(pdf);
  can->Print(pdf+"]");

  printf("\nProduced:\n  %s\n\n",pdf.Data());
}
