# ZJetBalance
Run II studies for Z+Jet Balancing including b-tagging

```
setupATLAS
rcSetup Base,2.3.21
git clone https://github.com/UCATLAS/xAODAnaHelpers
git clone https://github.com/UCATLAS/ZJetBalance
python xAODAnaHelpers/scripts/checkoutASGtags.py 2.3.21
rc checkout_pkg atlasperf/CombPerf/JetETMiss/JetCalibrationTools/JES_ResponseFitter/trunk JES_ResponseFitter
patch -p0 -i ZJetBalance/patches/JESResponse.patch
rc find_packages
rc compile
```

There are three executables ( all in ZJetBalance/util ):
  runZBJetBalance - Process xAOD with basic object/event selection after object calibration and writes out small analysis NTuple
  runProcessZJetBalanceMiniTree - last of object/event selection and fills histograms
  ZJetBalancePlotter - 

Then to run an example xAOD to make an NTuple :
```
runZBJetBalance -inFile /afs/cern.ch/user/g/gfacini/public/mc15_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.DAOD_SUSY5.e3601_s2576_s2132_r6633_r6264_p2370_tid05768578_00/DAOD_SUSY5.05768578._000001.pool.root.1
```
The output tree is in submitDir/data-tree/.


Then, the next step is to fill histograms from the output NTuple
```
runProcessZJetBalanceMiniTree  -inFile submitDir/data-tree/mc15_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.DAOD_SUSY5.e3601_s2576_s2132_r6633_r6264_p2370_tid05768578_00.root -submitDir submitDir2
```

Then the output histogram will be fed into the final JES_Response_Fitter tool as following.
You may use submitDir/hist-mc15_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.DAOD_SUSY5.e3601_s2576_s2132_r6633_r6264_p2370_tid05768578_00.root but use another example with enough statsitics. 
```
ZJetBalancePlotter -inFile ${ROOTCOREBIN}/data/ZJetBalance/PlotterExampleInput.root
```
Then you will see Zjet_DB_Gauss_fits.pdf as output.


See TWiki for more information
