#include <ZJetBalance/DH_GBHOut.h>


ZJetBalance::DH_GBHOut::DH_GBHOut(const std::string& outputFileName) : DrawingHelperOk(outputFileName), m_color_b(kBlue-3), m_color_c(kGreen), m_color_l(kYellow+2)
{
}


void 
ZJetBalance::DH_GBHOut::DrawFlavorComposition(const std::string& histname,
					      const std::string& comment,
					      const std::string& xtitle,
					      const std::string& label,
					      const std::string& mcDrawOption,
					      const bool& setYRange,
					      const double& yMinimum,
					      const double& yMaximum,
					      const bool& setXRange,
					      const double& xMinimum,
					      const double& xMaximum)
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
  std::vector<TH1F*>  mcHists(3); // for three flavor
  std::vector<double> mcEntries(3); // for three flavor
  std::vector<std::string> mcSampleTitle(3); // for three flavor
  
  // initializer
  mcEntries[0]=0.; mcEntries[1]=0.; mcEntries[2]=0.;
  mcSampleTitle[0]="b"; mcSampleTitle[1]="c"; mcSampleTitle[2]="others";
  mcHists[0] = new TH1F(); mcHists[1] = new TH1F(); mcHists[2] = new TH1F(); 
  
  for (int iMC=0, nMCs=m_MC_fileNames.size(); iMC<nMCs; iMC++)  {
    TH1F* tmp_b = PrepareTH1F(mcFiles.at(iMC),
			      histname+"_b",
			      xtitle,
			      m_MC_colors.at(iMC),
			      m_MC_styles.at(iMC),
			      true, // fillHistogram
			      m_MC_normalizationFactor.at(iMC));
    TH1F* tmp_c = PrepareTH1F(mcFiles.at(iMC),
			      histname+"_c",
			      xtitle,
			      m_MC_colors.at(iMC),
			      m_MC_styles.at(iMC),
			      true, // fillHistogram
			      m_MC_normalizationFactor.at(iMC));
    TH1F* tmp_l = PrepareTH1F(mcFiles.at(iMC),
			      histname+"_l",
			      xtitle,
			      m_MC_colors.at(iMC),
			      m_MC_styles.at(iMC),
			      true, // fillHistogram
			      m_MC_normalizationFactor.at(iMC));
    if (iMC==0) { // copy
      tmp_b->Copy(*(mcHists[0]));
      tmp_c->Copy(*(mcHists[1]));
      tmp_l->Copy(*(mcHists[2]));
      mcHists[0]->SetName("__tmp_b__");
      mcHists[1]->SetName("__tmp_c__");
      mcHists[2]->SetName("__tmp_l__");
    } else {
      mcHists[0]->Add(tmp_b);
      mcHists[1]->Add(tmp_c);
      mcHists[2]->Add(tmp_l);
    }
    
    (mcEntries[0]) += tmp_b->Integral(-1, -1);
    (mcEntries[1]) += tmp_c->Integral(-1, -1);
    (mcEntries[2]) += tmp_l->Integral(-1, -1);
  }
  
  
  std::vector<TH1F> mcHistStack(mcHists.size());
  for (int iFlavor=0, nFlavors=mcHists.size(); iFlavor<nFlavors; iFlavor++) {
    mcHists.at(iFlavor)->Copy(mcHistStack[iFlavor]);
    for (int kFlavor=iFlavor+1; kFlavor<nFlavors; kFlavor++) {
      mcHistStack[iFlavor] = mcHistStack[iFlavor] + (*mcHists.at(kFlavor));
    }
  }
  
  SetTH1FColors(&(mcHistStack[0]), m_color_b, true);
  SetTH1FColors(&(mcHistStack[1]), m_color_c, true);
  SetTH1FColors(&(mcHistStack[2]), m_color_l, true);
  

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
  
  leg.AddEntry(hData, 
	       Form("Data (%.0f)", 
		    hData->Integral(-1, -1)),
	       "PL");
  for (int iMC=0, nMCs=mcHistStack.size(); iMC<nMCs; iMC++) {
    leg.AddEntry( &(mcHistStack[iMC]),   
  		  Form("%s (%.0f)", 
  		       mcSampleTitle.at(iMC).c_str(),
  		       mcEntries.at(iMC)),
  		  "F");
  }
  
  leg.Draw();
  
  ATLASLabel(0.20, 0.94, label.c_str(), kBlack);
  
  TLatex myComment(0.45, 0.94, comment.c_str());
  myComment.SetNDC();
  myComment.Draw();
  
  CustimizeCampusWithRightMargin();
  m_canvas->Print(Form("%s.pdf", m_outputFileName.c_str()));
  
  // ratio plot
  RatioPlot(hData, mcHistStack, mcSampleTitle, mcEntries, comment, label, 0.5, mcDrawOption,
	    setYRange, yMinimum, yMaximum, setXRange, xMinimum, xMaximum);
  
  
  fData->Close();
  delete mcHists[0];
  delete mcHists[1];
  delete mcHists[2];
  CloseMCFiles(mcFiles);

  
  return;
}

void 
ZJetBalance::DH_GBHOut::SetTH1FColors(TH1F* h, 
				      const int& color, 
				      const bool& fillHistogram)
{
  h->SetLineColor(color);
  if (fillHistogram) {
    h->SetFillColor(color);
  } else {
    h->SetMarkerColor(color);
    h->SetMarkerStyle(8);
  }
}
