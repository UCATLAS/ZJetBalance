function RunZJetBalanceMiniTreeAnaWithCustomizedConfig {
    if [ $# -ne 9 ]; then
	echo "give configuration for 'AdditionalWeight', 'LumiCalcFiles', 'PRWFiles' '-inFile' 'pT_binning' 'eta_binning' 'BTagJets' 'BTagOP''-submitDir'"
    else 
	templaceConfigFile=${ROOTCOREBIN}/data/ZJetBalance/ZJetBalanceMiniTree_GenBalanceHistograms.config
	echo "using ${templaceConfigFile} as template"
	tmpConfFileName=__${RANDOM}_tmp.config
	grep -v 'AdditionalWeight' ${templaceConfigFile} | grep -v 'HistPrefix' | grep -v 'LumiCalcFiles' | grep -v 'PRWFiles' | grep -v 'pT_binning' | grep -v 'eta_binning' | grep -v 'BTagJets' | grep -v 'BTagOP' > ${tmpConfFileName}
	echo "AdditionalWeight ${1}" >> ${tmpConfFileName}
	echo "LumiCalcFiles ${2}" >> ${tmpConfFileName}
	echo "PRWFiles ${3}" >> ${tmpConfFileName}
	echo "pT_binning ${4}" >> ${tmpConfFileName}
	echo "eta_binning ${5}" >> ${tmpConfFileName}
	echo "BTagJets ${6}" >> ${tmpConfFileName}
	echo "BTagOP ${7}" >> ${tmpConfFileName}
	
	cat ${tmpConfFileName}
	
	runZJetBalanceMiniTreeAna -inFile ${8} -submitDir ${9} -algorithm 2 -algoConfig ${tmpConfFileName}
	\rm -r ${tmpConfFileName}
    fi
}

  
  


# Z to ee run
# xAOD  : 19890410
# DxAOD : 6380269
# additional weight =  0.320356
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 0.320356 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_HLT_e24_lhmedium_L1EM18VH_276262-276954.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361106.361107.361108.410000.r6282.root' '30,200' '-2.5,2.5' 'True' 'Fix70' '/afs/cern.ch/work/b/btuan/public/output/25ns/electron/mcZee/all.root' 'el_hist361106'
rm -f all_el_hist361106.root
hadd all_el_hist361106.root el_hist361106/hist-*.root

# Z to tau tau run
# xAOD  : 19152982
# DxAOD : 36342
# additional weight = 0.001893
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 0.001893 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_HLT_e24_lhmedium_L1EM18VH_276262-276954.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361106.361107.361108.410000.r6282.root' '30,200' '-2.5,2.5' 'True' 'Fix70' '/afs/cern.ch/work/b/btuan/public/output/25ns/electron/mcZtautau/all.root' 'el_hist361108'
rm -f all_el_hist361108.root
hadd all_el_hist361108.root el_hist361108/hist-*.root

# ttbar run
# xAOD  : 13392604
# DxAOD : 732290 
# additional weight = 0.036690
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 0.036690 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_HLT_e24_lhmedium_L1EM18VH_276262-276954.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361106.361107.361108.410000.r6282.root' '30,200' '-2.5,2.5' 'True' 'Fix70' '/afs/cern.ch/work/b/btuan/public/output/25ns/electron/mcttbar/all.root' 'el_hist410000' 
rm -f all_el_hist410000.root
hadd all_el_hist410000.root el_hist410000/hist-*.root

# Data
# with dummy input for PRW
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 1.0 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_HLT_e24_lhmedium_L1EM18VH_276262-276954.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361106.361107.361108.410000.r6282.root' '30,200' '-2.5,2.5' 'True' 'Fix70' '/afs/cern.ch/work/b/btuan/public/output/25ns/electron/data/all.root' 'el_histData25ns'
rm -f all_el_histData25ns.root
hadd all_el_histData25ns.root el_histData25ns/hist-*.root


# Run Jet Resoponse Fitter
rm -f all_el_MC.root 
hadd all_el_MC.root all_el_hist361106.root all_el_hist361108.root all_el_hist410000.root

ZJetBalancePlotter -usePoisson -inFile all_el_MC.root -outTag allMC_el
ZJetBalancePlotter -usePoisson -inFile all_el_histData25ns.root -outTag allData_el



############################################################
# Z to mumu run
# xAOD  : 19962997
# DxAOD : 8737695 
# additional weight = 0.438137 
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 0.438137 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_HLT_mu20_iloose_L1MU15_276262-276954.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361106.361107.361108.410000.r6282.root' '30,200' '-2.5,2.5' 'True' 'Fix70' '/afs/cern.ch/work/b/btuan/public/output/25ns/muon/mcZmumu/all.root' 'mu_hist361107'
rm -f all_mu_hist361107.root
hadd all_mu_hist361107.root mu_hist361107/hist-*.root

# Z to tau tau run
# xAOD  : 19152982
# DxAOD : 36342 
# additional weight = 0.001893 
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 0.001893 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_HLT_mu20_iloose_L1MU15_276262-276954.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361106.361107.361108.410000.r6282.root'  '30,200' '-2.5,2.5' 'True' 'Fix70' '/afs/cern.ch/work/b/btuan/public/output/25ns/muon/mcZtautau/all.root' 'mu_hist361108'
rm -f all_mu_hist361108.root
hadd all_mu_hist361108.root mu_hist361108/hist-*.root

# ttbar run
# xAOD  : 13392604
# DxAOD : 732290 
# additional weight = 0.036690
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 0.036690 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_HLT_mu20_iloose_L1MU15_276262-276954.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361106.361107.361108.410000.r6282.root' '30,200' '-2.5,2.5' 'True' 'Fix70' '/afs/cern.ch/work/b/btuan/public/output/25ns/muon/mcttbar/all.root' 'mu_hist410000'
rm -f all_mu_hist410000.root
hadd all_mu_hist410000.root mu_hist410000/hist-*.root

# Data
# with dummy input for PRW
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 1.0 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_HLT_mu20_iloose_L1MU15_276262-276954.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361106.361107.361108.410000.r6282.root' '30,200' '-2.5,2.5' 'True' 'Fix70' '/afs/cern.ch/work/b/btuan/public/output/25ns/muon/data/all.root' 'mu_histData25ns'
rm -f all_mu_histData25ns.root
hadd all_mu_histData25ns.root mu_histData25ns/hist-*.root


# Run Jet Resoponse Fitter
rm -f all_mu_MC.root 
hadd all_mu_MC.root all_mu_hist361107.root all_mu_hist361108.root all_mu_hist410000.root

ZJetBalancePlotter -usePoisson -inFile all_mu_MC.root -outTag allMC_mu
ZJetBalancePlotter -usePoisson -inFile all_mu_histData25ns.root -outTag allData_mu


# Run Plotter
testDHGBHOutAndJERFitter | tee constants.tex

#pdflatex table.tex
