#include <ZJetBalance/DrawingHelperOk.h>

// ======================================
ZJetBalance::DrawingHelperOk::DrawingHelperOk()
  : m_outputFileName("outpuf_default"),
    m_luminosity(-1000.),
    m_showStat(false)
{
  m_canvas = new TCanvas();
  m_canvas->Print(Form("%s.pdf[", m_outputFileName.c_str()));
  gStyle->SetOptStat(0);
  
  TH1::SetDefaultSumw2();
}

// ======================================
ZJetBalance::DrawingHelperOk::DrawingHelperOk(const std::string& outputFileName)
  : m_outputFileName(outputFileName),
    m_luminosity(-1000.)
{
  m_canvas = new TCanvas();
  m_canvas->Print(Form("%s.pdf[", m_outputFileName.c_str()));
  gStyle->SetOptStat(0);
  
  TH1::SetDefaultSumw2();
}

// ======================================
ZJetBalance::DrawingHelperOk::~DrawingHelperOk()
{
  m_canvas->Print(Form("%s.pdf]", m_outputFileName.c_str()));
}

// ======================================
TFile* 
ZJetBalance::DrawingHelperOk::GetTFile(const std::string& filename)
{
  TFile* rc = TFile::Open(filename.c_str());
  if (!rc) {
    Error("GetTFile()", "given file name = %s does not exist.", filename.c_str());
    exit(EXIT_FAILURE);
  }
  
  return rc;
}

// ======================================
void 
ZJetBalance::DrawingHelperOk::ATLASLabel(Double_t x,Double_t y,const char* text,Color_t color) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetTextSize(0.05); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);

  double delx = 0.11*696*gPad->GetWh()/(472*gPad->GetWw());

  l.DrawLatex(x,y,"ATLAS");
  if (text) {
    TLatex p; 
    p.SetTextSize(0.05); 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+delx,y,text);
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }
}

// ======================================
TObject* 
ZJetBalance::DrawingHelperOk::GetObject(TFile* f, const std::string& name)
{
  TObject* rc = f->Get(name.c_str());
  if (!rc) {
    Error("GetObject()", "given object name = %s does not exist in %s.", name.c_str(), f->GetName());
    exit(EXIT_FAILURE);
  }
  
  return rc;
}


// ======================================
void 
ZJetBalance::DrawingHelperOk::CustimizeCampusWithRightMargin()
{
  CustimizeCampus(0.25, 0.15, 0.15, 0.10, 1, 1);
}

// ======================================
void 
ZJetBalance::DrawingHelperOk::CustimizeCampusWithoutRightMargin()
{
  CustimizeCampus(0.10, 0.15, 0.15, 0.10, 1, 1);
}

// ======================================
void 
ZJetBalance::DrawingHelperOk::CustimizeCampus(const double& rightMargin,
					      const double& leftMargin,
					      const double& bottomMargin,
					      const double& topMargin,
					      const int& gridx,
					      const int& gridy)
{
  m_canvas->SetRightMargin(rightMargin);
  m_canvas->SetLeftMargin(leftMargin);
  m_canvas->SetBottomMargin(bottomMargin);
  m_canvas->SetTopMargin(topMargin);
  m_canvas->SetGridy(gridy);
  m_canvas->SetGridx(gridx);
}

// ======================================
void 
ZJetBalance::DrawingHelperOk::AddMC(const std::string& filename,
				    const std::string& title,
				    const std::string& symbol,
				    const int& color,
				    const int& marker_sylte,
				    const bool& useWeightedSum,
				    const bool& doNormalize)
{
  m_MC_fileNames.push_back(filename);
  m_MC_sampleTitles.push_back(title);
  m_MC_colors.push_back(color);
  m_MC_styles.push_back(marker_sylte);
  m_MC_symbols.push_back(symbol);
  
  // for normalization calculation
  TFile* fMC   = GetTFile(filename);
  double MCTotal = (!doNormalize) ?
    1. : ( useWeightedSum ? 
	   ((TH1F*)GetObject(fMC, "cutflow_weighted"))->GetBinContent(1) : 
	   ((TH1F*)GetObject(fMC, "cutflow"))->GetBinContent(1) );
  const double McLuminosityWeight = (!doNormalize) ?
    1. : ( (m_luminosity>0) ? m_luminosity/MCTotal : 1. );
  
  if (m_luminosity<0) {
    Info("AddMC()", "Normalization Factor is not set to %s", title.c_str());
  }
  
  m_MC_normalizationFactor.push_back(McLuminosityWeight);
  fMC->Close();
  
  Info("AddMC()", "==========================================");
  Info("AddMC()", "Adding MC samples");
  Info("AddMC()", "File             : %s", filename.c_str());
  Info("AddMC()", "Title            : %s", title.c_str());
  Info("AddMC()", "Symbol           : %s", symbol.c_str());
  Info("AddMC()", "NomFactor        : %.2e", McLuminosityWeight);
  Info("AddMC()", "Normalization    : %s", (doNormalize ? "YES" : "NO"));
  if (doNormalize) {
    Info("AddMC()", "Use Weighted Sum : %s", (useWeightedSum ? "YES" : "NO"));
    Info("AddMC()", "MCTotal          : %.2e", MCTotal);
  }
}

// ======================================
TH1F* 
ZJetBalance::DrawingHelperOk::PrepareTH1F(TFile* fileobj,
					  const std::string& histname,
					  const std::string& xtitle,
					  const int& color,
					  const int& markerStyle,
					  const bool& fillHistogram,
					  const double& normalization)
{
  // data preparation 
  TH1F* h = (TH1F*)GetObject(fileobj, histname);
  h->GetXaxis()->SetTitle(xtitle.c_str());
  h->SetTitle("");
  h->SetLineColor(color);
  if (fillHistogram) {    
    h->SetFillColor(color);
    h->SetMarkerStyle(0);
    h->SetMarkerSize(0);
  } else {
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerStyle);
  }
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetXaxis()->SetNdivisions(510);
  h->GetYaxis()->SetNdivisions(510);  
  h->Scale(normalization);
  
  return h;
}

// ======================================
std::vector<TFile*> 
ZJetBalance::DrawingHelperOk::OpenAndReturnMCFiles()
{
  // ======================================
  // this function leaks memory on purpose
  // please call CloseMCFiles at the end of use
  // ======================================
  std::vector<TFile*> rc;
  for (int iFile=0, nFiles=m_MC_fileNames.size(); iFile<nFiles; iFile++) {
    TFile* tmp = GetTFile(m_MC_fileNames.at(iFile));
    rc.push_back(tmp);
  }  
  
  return rc;
}

// ======================================
void 
ZJetBalance::DrawingHelperOk::CloseMCFiles(const std::vector<TFile*>& files)
{
  for (int iFile=0, nFiles=files.size(); iFile<nFiles; iFile++) {
    files.at(iFile)->Close();
  }
}

// ======================================
void 
ZJetBalance::DrawingHelperOk::RatioPlot(TH1F* hData,
					std::vector<TH1F>& mcHistStack,
					const std::vector<std::map<std::string, double> >& mcStats,
					const std::string& comment,
					const std::string& label,
					const double& ratio_plot_range_min,
					const double& ratio_plot_range_max,
					const std::string& mcDrawOption,
					const bool& setYRange,
					const double& yMinimum,
					const double& yMaximum,
					const bool& setXRange,
					const double& xMinimum,
					const double& xMaximum)
{
  // a function for backward compatibility
  RatioPlot(hData,
	    mcHistStack,
	    m_MC_sampleTitles,
	    mcStats,
	    comment,
	    label,
	    ratio_plot_range_min,
	    ratio_plot_range_max,
	    mcDrawOption,
	    setYRange,
	    yMinimum,
	    yMaximum,
	    setXRange,
	    xMinimum,
	    xMaximum);
}

// ======================================
void 
ZJetBalance::DrawingHelperOk::RatioPlot(TH1F* hData,
					std::vector<TH1F>& mcHistStack,
					const std::vector<std::string>& mcSampleTitles,
					const std::vector<std::map<std::string, double> >& mcStats,
					const std::string& comment,
					const std::string& label,
					const double& ratio_plot_range_min,
					const double& ratio_plot_range_max,
					const std::string& mcDrawOption,
					const bool& setYRange,
					const double& yMinimum,
					const double& yMaximum,
					const bool& setXRange,
					const double& xMinimum,
					const double& xMaximum
					)
{
  m_canvas->Clear();
  m_canvas->Update();
  CustimizeCampusWithRightMargin();
  
  TH1F* hMC = (&mcHistStack[0]);
  
  TString canvasname = m_canvas->GetName();
  m_canvas->Divide(1, 2);
  
  TH1D h_ratio_data;
  TH1D h_ratio_mc;
  hMC->Copy(h_ratio_mc);
  hData->Copy(h_ratio_data);
  
  h_ratio_mc.Clear();
  h_ratio_mc.SetName("ratio_data");
  h_ratio_data.Clear();
  h_ratio_data.SetName("ratio_data");
  
  h_ratio_mc.GetXaxis()->SetTitle(hMC->GetXaxis()->GetTitle());
  h_ratio_data.GetXaxis()->SetTitle(hMC->GetXaxis()->GetTitle());
  
  for (int iBin=1; iBin<=hMC->GetNbinsX(); iBin++) {
    const double mc_entry   = hMC->GetBinContent(iBin);
    const double data_entry = hData->GetBinContent(iBin);
    if (mc_entry==0) {continue;}
    
    const double mc_entry_uncert   = hMC->GetBinError(iBin);
    const double data_entry_uncert = hData->GetBinError(iBin);
    
    const double ratio_data       = data_entry/mc_entry;
    const double ratio_data_error = ratio_data*data_entry_uncert/data_entry;

    const double ratio_mc       = 1.0;
    const double ratio_mc_error = ratio_mc*mc_entry_uncert/mc_entry;
    
    h_ratio_data.SetBinContent(iBin, ratio_data);
    h_ratio_mc.SetBinContent(iBin, ratio_mc);

    h_ratio_data.SetBinError(iBin, ratio_data_error);
    h_ratio_mc.SetBinError(iBin, ratio_mc_error);
  }
  
  TPad* canvas_up = (TPad*) m_canvas->GetListOfPrimitives()->FindObject(canvasname+"_1");
  TPad* canvas_dw = (TPad*) m_canvas->GetListOfPrimitives()->FindObject(canvasname+"_2");
  
  // define the size
  double up_height     = 0.75; // please tune so that the upper figures size will meet your requirement
  double dw_correction = 1.28; // please tune so that the smaller canvas size will work in your environment
  //double font_size_dw  = 0.1; // please tune the font size parameter for bottom figure
  
  double dw_height    = (1. - up_height) * dw_correction;
  
  // set pad size
  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
  canvas_dw->SetPad(0., 0, 1., dw_height);
  canvas_up->SetFrameFillColor(0);
  canvas_up->SetLeftMargin(0.15);
  canvas_dw->SetLeftMargin(0.15);
  canvas_dw->SetRightMargin(0.25);
  canvas_up->SetRightMargin(0.25);
  canvas_up->SetGridx(1);
  canvas_dw->SetGridx(1);
  canvas_up->SetGridy(1);
  canvas_dw->SetGridy(1);

  canvas_up->SetFillColor(0);
  canvas_dw->SetFillColor(0);
  canvas_dw->SetBottomMargin(0.5);
  canvas_dw->SetFrameFillColor(0);
  
  // set top margin 0 for bottom figure
  canvas_dw->SetTopMargin(0);
  
  // draw top figure
  canvas_up->cd();
  
  if (setYRange) {
    hData->SetMaximum(yMaximum);
    hData->SetMinimum(yMinimum);
    hMC->SetMaximum(yMaximum);
    hMC->SetMinimum(yMinimum);
  } else if (mcDrawOption=="H") { // nominal histograms (auto range)
    hData->SetMinimum(hData->GetMinimum()<0 ? hData->GetMinimum() : 0.);
    hMC->SetMinimum(hMC->GetMinimum()<0 ? hMC->GetMinimum() : 0.);
  }
  
  if (setXRange) {
    hData->GetXaxis()->SetRangeUser(xMinimum, xMaximum);
    hMC->GetXaxis()->SetRangeUser(xMinimum, xMaximum);
  }
  
  if (hMC->GetMaximum()<hData->GetMaximum()) {
    hData->Draw("PE");
  } else {
    hMC->Draw(Form("%s", mcDrawOption.c_str()));
  }
  
  hMC->Draw(Form("%s SAME", mcDrawOption.c_str()));
  for (int iMC=0, nMCs=mcHistStack.size(); iMC<nMCs; iMC++) {
    mcHistStack.at(iMC).Draw(Form("%s SAME", mcDrawOption.c_str()));
  }
  hData->Draw("PE SAME");
  
  TLegend leg(0.76, 0.15, 0.98, 0.90);
  leg.SetLineStyle(0);
  leg.SetLineColor(0);
  leg.SetFillStyle(0);
  
  std::string dataLegend = ReturnLegend("Data", hData->Integral(-1, -1), hData->GetMean(), hData->GetRMS(), mcDrawOption, true);
  leg.AddEntry(hData, dataLegend.c_str(), "PL");
  for (int iMC=0, nMCs=mcHistStack.size(); iMC<nMCs; iMC++) {
    std::string mcLegend = ReturnLegend(mcSampleTitles.at(iMC),mcStats.at(iMC), mcDrawOption, false);
    leg.AddEntry( &(mcHistStack[iMC]), mcLegend.c_str(), "F");
  }
      
  ATLASLabel(0.20, 0.94, label.c_str(), kBlack);
  
  TLatex myComment(0.45, 0.94, comment.c_str());
  myComment.SetNDC();
  myComment.Draw();
  
  canvas_dw->cd();
  canvas_dw->SetGridy();
  const std::string xtitle = hData->GetXaxis()->GetTitle();
  TLine l1((setXRange ? xMinimum : h_ratio_data.GetXaxis()->GetXmin()), 1.0, 
	   (setXRange ? xMaximum : h_ratio_data.GetXaxis()->GetXmax()), 1.0);
  // font size
  h_ratio_mc.SetMaximum(ratio_plot_range_max);
  h_ratio_mc.SetMinimum(ratio_plot_range_min);
  h_ratio_mc.SetFillColor(kYellow+3);
  h_ratio_mc.GetXaxis()->SetLabelSize(0.15);
  h_ratio_mc.GetXaxis()->SetTitleSize(0.15);
  h_ratio_mc.GetYaxis()->SetLabelSize(0.1);
  h_ratio_mc.GetYaxis()->SetTitleSize(0.15);
  h_ratio_mc.GetYaxis()->SetTitleOffset(0.3);
  h_ratio_mc.GetXaxis()->SetTitle(xtitle.c_str());
  h_ratio_mc.GetXaxis()->SetNdivisions(505);
  h_ratio_mc.GetYaxis()->SetNdivisions(505);
  h_ratio_mc.GetYaxis()->SetTitle("Ratio");
  h_ratio_mc.GetXaxis()->SetTitleOffset(1.2);
  
  h_ratio_data.SetLineColor(kBlack);
  h_ratio_data.SetMarkerColor(kBlack);
  h_ratio_data.SetMarkerStyle(8);
  
  
  h_ratio_mc.Draw("E2");
  h_ratio_data.Draw("PE SAME");
  l1.Draw();
  
  canvas_up->cd();  
  leg.Draw();
  
  m_canvas->Print(Form("%s.pdf", m_outputFileName.c_str()));
  
  m_canvas->Clear();
  m_canvas->Update();
  CustimizeCampusWithRightMargin();
}

void 
ZJetBalance::DrawingHelperOk::MyDataMcComparisonTH1F_GraphStyle(const std::string& histname,
								const std::string& comment,
								const std::string& xtitle,
								const std::string& label,
								const bool& setYRange,
								const double& yMinimum,
								const double& yMaximum,
								const bool& setXRange,
								const double& xMinimum,
								const double& xMaximum,
								const double& ratio_plot_range_min,
								const double& ratio_plot_range_max)
{
  MyDataMcComparisonTH1F(histname,
			 comment,
			 xtitle,
			 label,
			 "E2",
			 setYRange,
			 yMinimum,
			 yMaximum,
			 setXRange,
			 xMinimum,
			 xMaximum,
			 ratio_plot_range_min,
			 ratio_plot_range_max);
}
  
// ======================================
void 
ZJetBalance::DrawingHelperOk::MyDataMcComparisonTH1F(const std::string& histname,
						     const std::string& comment,
						     const std::string& xtitle,
						     const std::string& label,
						     const std::string& mcDrawOption,
						     const bool& setYRange,
						     const double& yMinimum,
						     const double& yMaximum,
						     const bool& setXRange,
						     const double& xMinimum,
						     const double& xMaximum,
						     const double& ratio_plot_range_min,
						     const double& ratio_plot_range_max)
{    
  // data preparation 
  TFile* fData = GetTFile(m_Data_fileName);
  
  TH1F* hData = PrepareTH1F(fData,
			    histname,
			    xtitle,
			    kBlack,
			    8,
			    false, // not fillHistogram
			    1.0);
  
  std::vector<TFile*> mcFiles = OpenAndReturnMCFiles();
  std::vector<TH1F*> mcHists;
  std::vector<std::map<std::string, double> > mcStats;
  
  if (m_MC_fileNames.size()==0) {
    Error("MyDataMcComparison()", "no MC files registered. no draw for %s", 
	  histname.c_str());
    return;
  }
  
  for (int iMC=0, nMCs=m_MC_fileNames.size(); iMC<nMCs; iMC++)  {
    TH1F* tmp = PrepareTH1F(mcFiles.at(iMC),
			    histname,
			    xtitle,
			    m_MC_colors.at(iMC),
			    m_MC_styles.at(iMC),
			    true, // fillHistogram
			    m_MC_normalizationFactor.at(iMC));
    mcHists.push_back(tmp);
    std::map<std::string, double> stats = ReturnStatsMap(tmp);
    mcStats.push_back(stats);
  }
  
  std::vector<TH1F> mcHistStack(mcHists.size());
  for (int iFile=0, nFiles=mcHists.size(); iFile<nFiles; iFile++) {
    mcHists.at(iFile)->Copy(mcHistStack[iFile]);
    for (int kFile=iFile+1; kFile<nFiles; kFile++) {
      mcHistStack[iFile] = mcHistStack[iFile] + (*mcHists.at(kFile));
    }
  }
  
  TH1F* hMC = (& (mcHistStack[0]) );
  if (setYRange) {
    hData->SetMaximum(yMaximum);
    hData->SetMinimum(yMinimum);
    hMC->SetMaximum(yMaximum);
    hMC->SetMinimum(yMinimum);
  } else if (mcDrawOption=="H") { // nominal histograms (auto range)
    hData->SetMinimum(hData->GetMinimum()<0 ? hData->GetMinimum() : 0.);
    hMC->SetMinimum(hMC->GetMinimum()<0 ? hMC->GetMinimum() : 0.);
  } 
  
  if (setXRange) {
    hData->GetXaxis()->SetRangeUser(xMinimum, xMaximum);
    hMC->GetXaxis()->SetRangeUser(xMinimum, xMaximum);
  }
  
  if (hMC->GetMaximum()<hData->GetMaximum()) {
    hData->GetXaxis()->SetTitle(xtitle.c_str());
    hData->Draw("PE");
  } else {
    hMC->GetXaxis()->SetTitle(xtitle.c_str());
    hMC->Draw(Form("%s", mcDrawOption.c_str()));
  }
  
  hMC->Draw(Form("%s SAME", mcDrawOption.c_str()));
  for (int iMC=0, nMCs=mcHistStack.size(); iMC<nMCs; iMC++) {
    mcHistStack.at(iMC).Draw(Form("%s SAME", mcDrawOption.c_str()));
  }
  hData->Draw("PE SAME");
  
  TLegend leg(0.76, 0.15, 0.98, 0.90);
  leg.SetLineStyle(0);
  leg.SetLineColor(0);
  leg.SetFillStyle(0);
  
  std::string dataLegend = ReturnLegend("Data", hData->Integral(-1, -1), hData->GetMean(), hData->GetRMS(), mcDrawOption, true);
  leg.AddEntry(hData, dataLegend.c_str(), "PL");
  for (int iMC=0, nMCs=mcHistStack.size(); iMC<nMCs; iMC++) {
    std::string mcLegend = ReturnLegend(m_MC_sampleTitles.at(iMC),mcStats.at(iMC), mcDrawOption, false);
    leg.AddEntry( &(mcHistStack[iMC]), mcLegend.c_str(), "F");
  }
  
  leg.Draw();
  
  ATLASLabel(0.20, 0.94, label.c_str(), kBlack);

  TLatex myComment(0.45, 0.94, comment.c_str());
  myComment.SetNDC();
  myComment.Draw();
  
  CustimizeCampusWithRightMargin();
  m_canvas->Print(Form("%s.pdf", m_outputFileName.c_str()));

  
  // ratio plot
  RatioPlot(hData, mcHistStack, mcStats, comment, label, 
	    ratio_plot_range_min, ratio_plot_range_max, 
	    mcDrawOption, setYRange, yMinimum, yMaximum, setXRange, xMinimum, xMaximum);
  
  CloseMCFiles(mcFiles);
}

// ======================================
void
ZJetBalance::DrawingHelperOk::DumpCutFlow(const std::string& commandPrefix, // e.g. "WithBTag_"
					  const std::string& cutflowHistogramName)
{
  static const std::string CutStageIDs[] = 
    {"AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK", "AL", "AM", "AN", "AO", "AP", "AQ", "AR", "AS", "AT", "AU", "AV", "AW", "AX", "AY", "AZ",
     "BA", "BB", "BC", "BD", "BE", "BF", "BG", "BH", "BI", "BJ", "BK", "BL", "BM", "BN", "BO", "BP", "BQ", "BR", "BS", "BT", "BU", "BV", "BW", "BX", "BY", "BZ"};
  
  printf("\\newcommand{\\%sLuminosityFb}{%.1f} \n", commandPrefix.c_str(), m_luminosity);
  printf("\\newcommand{\\%sLuminosityPb}{%.1f} \n", commandPrefix.c_str(), m_luminosity*1000);
  
  std::vector<double> cutflow_data;
  
  TFile* fData = GetTFile(m_Data_fileName);
  TH1F* h_cutflow_data = (TH1F*)GetObject(fData, cutflowHistogramName);
  for (int ii=0; ii<h_cutflow_data->GetNbinsX(); ii++) {
    const int iBin = ii+1;
    cutflow_data.push_back(h_cutflow_data->GetBinContent(iBin));
  }
  
  std::map< std::string, std::vector<double> > cutflow_mc;
  std::map<std::string, double> cutflow_total_mc;
  
  std::vector<TFile*> mcFiles = OpenAndReturnMCFiles();
  for (int iMC=0, nMCs=mcFiles.size(); iMC<nMCs; iMC++)  {
    TFile* fMC = mcFiles.at(iMC);
    const double& weight = m_MC_normalizationFactor.at(iMC);
    const std::string& symbol = m_MC_symbols.at(iMC);
    
    TH1F* h_cutflow_mc = (TH1F*)GetObject(fMC, cutflowHistogramName);
    
    cutflow_mc[symbol]     = std::vector<double>();
    
    for (int ii=0; ii<h_cutflow_mc->GetNbinsX(); ii++) {
      const int iBin = ii+1;
      const std::string& cutId = CutStageIDs[iBin];
      //Info("DumpCutFlow()", "%10s weight=%.2e h_cutflow_mc_weighted->GetBinContent(%d)=%.2e",
      //      symbol.c_str(), weight, iBin, h_cutflow_mc_weighted->GetBinContent(iBin));
      cutflow_mc[symbol].push_back(weight*h_cutflow_mc->GetBinContent(iBin));
      if (cutflow_total_mc.find(cutId)==cutflow_total_mc.end()) {
	cutflow_total_mc[cutId] = weight*h_cutflow_mc->GetBinContent(iBin);
      } else {
	cutflow_total_mc[cutId] += weight*h_cutflow_mc->GetBinContent(iBin);
      }
    }
  }
  
  for (int ii=0, n=cutflow_data.size(); ii<n; ii++) {
    const int iBin = ii + 1;
    const std::string& cutId = CutStageIDs[iBin];
    printf("\\newcommand{\\%sCF%sData}{%.0f} \n", commandPrefix.c_str(), cutId.c_str(), cutflow_data.at(ii));
    if (ii!=0) {
      if (cutflow_data.at(ii-1)>0.001) {
	printf("\\newcommand{\\%sCF%sDataRelEffPerC}{%.0f} \n", commandPrefix.c_str(), cutId.c_str(), cutflow_data.at(ii)/cutflow_data.at(ii-1)*100.);
	printf("\\newcommand{\\%sCF%sDataRelEff}{%.2f} \n", commandPrefix.c_str(), cutId.c_str(), cutflow_data.at(ii)/cutflow_data.at(ii-1));
      }
      else {
	printf("\\newcommand{\\%sCF%sDataRelEffPerC}{-} \n", commandPrefix.c_str(), cutId.c_str());
	printf("\\newcommand{\\%sCF%sDataRelEff}{-} \n", commandPrefix.c_str(), cutId.c_str());
      }
    } else {
      printf("\\newcommand{\\%sCF%sDataRelEffPerC}{-} \n", commandPrefix.c_str(), cutId.c_str());
      printf("\\newcommand{\\%sCF%sDataRelEff}{-} \n", commandPrefix.c_str(), cutId.c_str());
    }
  }
  
  for (int iMC=0, nMCs=mcFiles.size(); iMC<nMCs; iMC++)  {
    const std::string& symbol = m_MC_symbols.at(iMC);
    for (int ii=0, n=cutflow_mc[symbol].size(); ii<n; ii++) {
      const int iBin = ii + 1;
      const std::string& cutId = CutStageIDs[iBin];
      printf("\\newcommand{\\%sCF%s%s}{%.1f} \n", commandPrefix.c_str(), cutId.c_str(), symbol.c_str(), cutflow_mc[symbol].at(ii));
      
      if (ii!=0) {
	if (cutflow_mc[symbol].at(ii-1)>0.001) {
	  printf("\\newcommand{\\%sCF%s%sRelEffPerC}{%.0f} \n", commandPrefix.c_str(), cutId.c_str(), symbol.c_str(), cutflow_mc[symbol].at(ii)/cutflow_mc[symbol].at(ii-1)*100.);
	  printf("\\newcommand{\\%sCF%s%sRelEff}{%.2f} \n", commandPrefix.c_str(), cutId.c_str(), symbol.c_str(), cutflow_mc[symbol].at(ii)/cutflow_mc[symbol].at(ii-1));
	}
	else {
	  printf("\\newcommand{\\%sCF%s%sRelEffPerC}{-} \n", commandPrefix.c_str(), cutId.c_str(), symbol.c_str());
	  printf("\\newcommand{\\%sCF%s%sRelEff}{-} \n", commandPrefix.c_str(), cutId.c_str(), symbol.c_str());
	}
      } else {
	printf("\\newcommand{\\%sCF%s%sRelEffPerC}{-} \n", commandPrefix.c_str(), cutId.c_str(), symbol.c_str());
	printf("\\newcommand{\\%sCF%s%sRelEff}{-} \n", commandPrefix.c_str(), cutId.c_str(), symbol.c_str());
      }     
    }
  }
  
  // MC total
  const std::string& symbol0 = m_MC_symbols.at(0);
  for (int ii=0, n=cutflow_mc[symbol0].size(); ii<n; ii++) {
    const int iBin = ii + 1;
    const std::string& cutId     = CutStageIDs[iBin];
    const std::string& cutIdPrev = CutStageIDs[iBin-1];
    printf("\\newcommand{\\%sCF%s%s}{%.1f} \n", commandPrefix.c_str(), cutId.c_str(), "MCTotal", cutflow_total_mc[cutId]);    
    
    if (ii!=0) {
      if (cutflow_total_mc[cutIdPrev]>0.001) {
	printf("\\newcommand{\\%sCF%s%sRelEffPerC}{%.0f} \n", commandPrefix.c_str(), cutId.c_str(), "MCTotal", cutflow_total_mc[cutId]/cutflow_total_mc[cutIdPrev]*100.);
	printf("\\newcommand{\\%sCF%s%sRelEff}{%.2f} \n", commandPrefix.c_str(), cutId.c_str(), "MCTotal", cutflow_total_mc[cutId]/cutflow_total_mc[cutIdPrev]);
      }
      else {
	printf("\\newcommand{\\%sCF%s%sRelEffPerC}{-} \n", commandPrefix.c_str(), cutId.c_str(), "MCTotal");
	printf("\\newcommand{\\%sCF%s%sRelEff}{-} \n", commandPrefix.c_str(), cutId.c_str(), "MCTotal");
	}
    } else {
      printf("\\newcommand{\\%sCF%s%sRelEffPerC}{-} \n", commandPrefix.c_str(), cutId.c_str(), "MCTotal");
      printf("\\newcommand{\\%sCF%s%sRelEff}{-} \n", commandPrefix.c_str(), cutId.c_str(), "MCTotal");
    }     
  }
  
  fData->Close();
  CloseMCFiles(mcFiles);
}
  
// ======================================
std::map<std::string, double> 
ZJetBalance::DrawingHelperOk::ReturnStatsMap(TH1* h)
{
  std::map<std::string, double> rc;
  rc["integral"] = h->Integral(-1, -1);
  rc["mean"] = h->GetMean();
  rc["rms"] = h->GetRMS();
  return rc;
}

// ======================================
std::string 
ZJetBalance::DrawingHelperOk::ReturnLegend(const std::string& title,
					   const std::map<std::string, double> stats,
					   const std::string& drawOption,
					   bool isData)
{
  return ReturnLegend(title,
		      stats.at("integral"),
		      stats.at("mean"),
		      stats.at("rms"),
		      drawOption,
		      isData);
}

// ======================================
std::string 
ZJetBalance::DrawingHelperOk::ReturnLegend(const std::string& title,
					   const double& entry,
					   const double& mean,
					   const double& rms,
					   const std::string& drawOption,
					   bool isData)
{
  std::string rc;
  
  if (drawOption!="H") {
    rc = title;
  } else {
    if (!m_showStat) {
      rc = (isData ? 
	    Form("%s (%.0f)", title.c_str(), entry) : 
	    Form("%s (%.1f)", title.c_str(), entry));
    } else {
      rc = (isData ? 
	    Form("#splitline{%s (%.0f)}{ave=%.1f rms=%.1f}", title.c_str(), entry, mean, rms) :
	    Form("#splitline{%s (%.1f)}{ave=%.1f rms=%.1f}", title.c_str(), entry, mean, rms));
	    //Form("#splitline{%s (%.0f)}{#splitline{mean=%.1f}{RMS=%.1f}}", title.c_str(), entry, mean, rms) :
	    //Form("#splitline{%s (%.1f)}{#splitline{mean=%.1f}{RMS=%.1f}}", title.c_str(), entry, mean, rms));
    }
  }
  
  return rc;
}
