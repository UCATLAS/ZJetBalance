#ifndef DrawingHelperOk_H
#define DrawingHelperOk_H

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
#include <TLine.h>

#include <vector>
#include <string>
#include <iostream>

#include <string.h>
#include <stdlib.h>

namespace ZJetBalance {
  
  class DrawingHelperOk
  {
  public:
    DrawingHelperOk();
    DrawingHelperOk(const std::string& outputFileName);
    ~DrawingHelperOk();
    static TFile* GetTFile(const std::string& filename);
    static TObject* GetObject(TFile* f, const std::string& name);
    static void ATLASLabel(Double_t x, Double_t y, const char* text, Color_t color);
    
    void AddMC(const std::string& filename,
	       const std::string& title,
	       const std::string& symbol,
	       const int& color,
	       const int& marker_sylte,
	       const bool& useWeightedSum,
	       const bool& doNormalize = true /* could be false in case we do not need normalize the distribution, for e.g. efficiency */
	       );
    void SetShowStat(bool showStat) { m_showStat=showStat; }
    void SetDataFileName(const std::string& name) {m_Data_fileName=name;} 
    void SetLuminosity(const double& luminosity) {m_luminosity=luminosity;}
    double GetLuminosity() {return m_luminosity;}
    void MyDataMcComparisonTH1F(const std::string& histname,
				const std::string& comment,
				const std::string& xtitle,
				const std::string& label="Internal",
				const std::string& mcDrawOption="H",
				const bool& setYRange=false,
				const double& yMinimum=-1,
				const double& yMaximum=-1,
				const bool& setXRange=false,
				const double& xMinimum=-1,
				const double& xMaximum=-1,
				const double& ratio_plot_range_min=0.5,
				const double& ratio_plot_range_max=1.5);
    void MyDataMcComparisonTH1F_GraphStyle(const std::string& histname,
					   const std::string& comment,
					   const std::string& xtitle,
					   const std::string& label="Internal",
					   const bool& setYRange=false,
					   const double& yMinimum=-1,
					   const double& yMaximum=-1,
					   const bool& setXRange=false,
					   const double& xMinimum=-1,
					   const double& xMaximum=-1,
					   const double& ratio_plot_range_min=0.5,
					   const double& ratio_plot_range_max=1.5);
    void RatioPlot(TH1F* hData,
		   std::vector<TH1F>& mcHistStack,
		   const std::vector<std::string>& mcSampleTitles,
		   const std::vector<std::map<std::string, double> >& mcStats,
		   const std::string& comment,
		   const std::string& label,
		   const double& ratio_plot_range_min,
		   const double& ratio_plot_range_max,
		   const std::string& mcDrawOption="H",
		   const bool& setYRange=false,
		   const double& yMinimum=-1,
		   const double& yMaximum=-1,
		   const bool& setXRange=false,
		   const double& xMinimum=-1,
		   const double& xMaximum=-1);
    void RatioPlot(TH1F* hData,
		   std::vector<TH1F>& mcHistStack,
		   const std::vector<std::map<std::string, double> >& mcStats,
		   const std::string& comment,
		   const std::string& label,
		   const double& ratio_plot_range_min,
		   const double& ratio_plot_range_max,
		   const std::string& mcDrawOption="H",
		   const bool& setYRange=false,
		   const double& yMinimum=-1,
		   const double& yMaximum=-1,
		   const bool& setXRange=false,
		   const double& xMinimum=-1,
		   const double& xMaximum=-1);
    void DumpCutFlow(const std::string& commandPrefix="", // e.g. "WithBTag_"
		     const std::string& cutflowHistogramName="cutflow_weighted_final");
    std::string ReturnLegend(const std::string& title,
			     const double& entry,
			     const double& mean,
			     const double& rms,
			     const std::string& drawOption,
			     bool isData);
    std::string ReturnLegend(const std::string& title,
			     const std::map<std::string, double> stats,
			     const std::string& drawOption,
			     bool isData);
      
  protected:
    // const memeber
    const std::string        m_outputFileName;
    
    std::vector<std::string> m_MC_fileNames;
    std::vector<std::string> m_MC_sampleTitles;
    std::vector<int>         m_MC_colors;
    std::vector<int>         m_MC_styles;
    std::vector<double>      m_MC_normalizationFactor;
    std::vector<bool>        m_MC_userSumOfWeight;
    std::vector<std::string> m_MC_symbols;
    std::string              m_Data_fileName;
    double                   m_luminosity;
    bool                     m_showStat;
    
    TCanvas* m_canvas;
    void CustimizeCampusWithRightMargin();
    void CustimizeCampusWithoutRightMargin();
    void CustimizeCampus(const double& rightMargin,
			 const double& leftMargin,
			 const double& bottomMargin,
			 const double& topMargin,
			 const int& gridx,
			 const int& gridy);
    std::vector<TFile*> OpenAndReturnMCFiles();
    void CloseMCFiles(const std::vector<TFile*>& files);
    TH1F* PrepareTH1F(TFile* fileobj,
		      const std::string& histname,
		      const std::string& xtitle,
		      const int& color,
		      const int& markerStyle,
		      const bool& fillHistogram,
		      const double& normalization = 1.0);
    std::map<std::string, double> ReturnStatsMap(TH1* h);
  };

}



#endif 
