# ZJetBalance
Run II studies for Z+Jet Balancing including b-tagging

```
setupATLAS
rcSetup Base,2.3.23
git clone https://github.com/UCATLAS/xAODAnaHelpers
git clone https://github.com/UCATLAS/ZJetBalance
python xAODAnaHelpers/scripts/checkoutASGtags.py 2.3.21
rc checkout_pkg atlasperf/CombPerf/JetETMiss/JetCalibrationTools/JES_ResponseFitter/trunk JES_ResponseFitter
patch -p0 -i ZJetBalance/patches/JESResponse.patch
rc find_packages
rc compile
```

There are four executables ( all in ZJetBalance/util ):
  runZBJetBalance - Process xAOD with basic object/event selection after object calibration and writes out small analysis NTuple
  runProcessZJetBalanceMiniTree - last of object/event selection and fills histograms
  ZJetBalancePlotter - run the JES_ResponseFitter/ and evaluate the JES response function
  drawZJetBlanace - run final plotter

Then to run an example xAOD to make an NTuple :
```
runZBJetBalance -inFile /afs/cern.ch/user/g/gfacini/public/mc15_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.DAOD_SUSY5.e3601_s2576_s2132_r6633_r6264_p2370_tid05768578_00/DAOD_SUSY5.05768578._000001.pool.root.1
```
The output tree is in submitDir/data-tree/.


Then, the next step is to fill histograms from the output NTuple
```
runProcessZJetBalanceMiniTree -inFile /eos/atlas/user/o/okumura/data/ZBJetBalanceStudy_v20150807/mc15_13TeV_50ns/mc15_13TeV_50ns_361107 -submitDir submitDirZMC

runProcessZJetBalanceMiniTree -inFile /eos/atlas/user/o/okumura/data/ZBJetBalanceStudy_v20150807/data15_13TeV_50ns -submitDir submitDirData

# do hadd
rm -f ZmumuAll.root
hadd ZmumuAll.root submitDirZMC/hist-user*.root

rm -f DataAll.root
hadd DataAll.root submitDirData/hist-user*.root
```

Then the output histogram will be fed into the final JES_Response_Fitter tool as following.

```
ZJetBalancePlotter -inFile ZmumuAll.root -outTag ZmumuAll

ZJetBalancePlotter -inFile DataAll.root -outTag DataAll
```

Then you will get the following four files.

```
Zjet_DB_Gauss_fits_ZmumuAll.pdf
ZJetBalancePlotterOutZmumuAll.root
Zjet_DB_Gauss_fits_DataAll.pdf
ZJetBalancePlotterOutDataAll.root
```

```
# (1) draw set of validation plots
# note : dxAOD skimming efficiency (Zmumu) = 9881474/19962997 = 4.949895e-01
drawZJetBlanace -d DataAll.root -m ZmumuAll.root -l 0.0783190 -w -a 4.949895e-01 -s 1

# (2) draw set of Jet Response Fitter
drawZJetBlanace -d ZJetBalancePlotterOutDataAll.root -m ZJetBalancePlotterOutZmumuAll.root -l 0.0783190 -w -a 4.949895e-01 -s 3
```

See TWiki for more information
