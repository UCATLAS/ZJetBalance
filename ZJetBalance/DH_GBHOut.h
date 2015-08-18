// Super class for Drawing Helper for ZJetBalanceMiniTree_GenBalanceHistograms.h output

#ifndef DH_GBHOut_H 
#define DH_GBHOut_H 

#include <ZJetBalance/DrawingHelperOk.h>

namespace ZJetBalance {
  class DH_GBHOut : public DrawingHelperOk {
  public:
    DH_GBHOut(const std::string& outputFileName);
    void DrawFlavorComposition(const std::string& histname,
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
    void SetFlavorColors(int color_b, int color_c, int color_l) { 
      m_color_b=color_b;
      m_color_c=color_c;
      m_color_l=color_l;
    }
    
  private:
    int m_color_b;
    int m_color_c;
    int m_color_l;
    void SetTH1FColors(TH1F* h, const int& color, const bool& fillHistogram);
  };
}

#endif 
