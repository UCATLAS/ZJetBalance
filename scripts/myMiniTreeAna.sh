function RunZJetBalanceMiniTreeAnaWithCustomizedConfig {
    if [ $# -ne 5 ]; then
	echo "give configuration for 'AdditionalWeight', 'LumiCalcFiles', 'PRWFiles' '-inFile' '-submitDir'"
    else 
	templaceConfigFile=${ROOTCOREBIN}/data/ZJetBalance/ZJetBalanceMiniTree_GenBalanceHistograms.config
	echo "using ${templaceConfigFile} as template"
	tmpConfFileName=__${RANDOM}_tmp.config
	grep -v 'AdditionalWeight' ${templaceConfigFile} | grep -v 'HistPrefix' | grep -v 'LumiCalcFiles' | grep -v 'PRWFiles' > ${tmpConfFileName}
	echo "AdditionalWeight ${1}" >> ${tmpConfFileName}
	echo "LumiCalcFiles ${2}" >> ${tmpConfFileName}
	echo "PRWFiles ${3}" >> ${tmpConfFileName}
	cat ${tmpConfFileName}
	
	runZJetBalanceMiniTreeAna -inFile ${4} -submitDir ${5} -algorithm 2 -algoConfig ${tmpConfFileName}
	\rm -r ${tmpConfFileName}
    fi
}

# Z to mumu run
# xAOD  : 19962997
# DxAOD : 9881474
# additional weight = 0.4949895
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 0.4949895 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_None_271298-271595_RUN2-UPD4-04.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361107.e3601_a766_a767_r6264.v00_METADATA.root' '/eos/atlas/user/o/okumura/data/ZBJetBalanceStudy_v20150810/mc15_13TeV_50ns/mc15_13TeV_50ns_361107' 'hist361107'
rm -f all_hist361107.root
hadd all_hist361107.root hist361107/hist-*.root

# Z to tau tau run
# xAOD  : 19152982
# DxAOD : 438544
# additional weight = 0.0228969
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 0.0228969 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_None_271298-271595_RUN2-UPD4-04.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361108.e3601_a766_a767_r6264.v00_METADATA.root' '/eos/atlas/user/o/okumura/data/ZBJetBalanceStudy_v20150810/mc15_13TeV_50ns/mc15_13TeV_50ns_361108' 'hist361108'
rm -f all_hist361108.root
hadd all_hist361108.root hist361108/hist-*.root

# ttbar run
# xAOD  : 13392604
# DxAOD : 3714174
# additional weight = 0.2773302
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 0.2773302 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_None_271298-271595_RUN2-UPD4-04.root' '$ROOTCOREBIN/data/ZJetBalance/prw.410004.e3601_a766_a767_r6264.v00_METADATA.root' '/eos/atlas/user/o/okumura/data/ZBJetBalanceStudy_v20150810/mc15_13TeV_50ns/mc15_13TeV_50ns_410004' 'hist410004'
rm -f all_hist410004.root
hadd all_hist410004.root hist410004/hist-*.root

# Data
# with dummy input for PRW
RunZJetBalanceMiniTreeAnaWithCustomizedConfig 1.0 '$ROOTCOREBIN/data/ZJetBalance/ilumicalc_histograms_None_271298-271595_RUN2-UPD4-04.root' '$ROOTCOREBIN/data/ZJetBalance/prw.361107.e3601_a766_a767_r6264.v00_METADATA.root' '/eos/atlas/user/o/okumura/data/ZBJetBalanceStudy_v20150810/data15_13TeV_50ns' 'histData50ns'
rm -f all_histData50ns.root
hadd all_histData50ns.root histData50ns/hist-*.root
