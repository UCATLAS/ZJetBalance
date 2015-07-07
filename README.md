# ZJetBalance
Run II studies for Z+Jet Balancing including b-tagging

```
setupATLAS
rcSetup Base,2.3.14
git clone https://github.com/UCATLAS/xAODAnaHelpers
git clone https://github.com/UCATLAS/ZJetBalance
python xAODAnaHelpers/scripts/checkoutASGtags.py 2.3.14
rc find_packages
rc compile
```

Then to run an example file try:
```
runZJetBalance -inFile /afs/cern.ch/user/g/gfacini/public/mc15_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.DAOD_SUSY5.e3601_s2576_s2132_r6633_r6264_p2370_tid05768578_00/DAOD_SUSY5.05768578._000001.pool.root.1
```
The output tree is in submitDir/data-tree/.

See TWiki for more information
