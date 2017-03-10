# Setup

Install with:

    git clone git@github.com:matz-e/cms-hcal-debug.git Debug/HcalDebug
    scram b -j 8

# Running

Run with:

    runTheMatrix.py -w upgrade -l 10039
    cmsRun Debug/HcalDebug/test/cmp_legacy.py

With customizations:

    cmsDriver.py something \
      --conditions auto:phase1_2017_realistic \
      -s RAW2DIGI,L1 -n 10 --geometry DB:Extended --era Run2_2017 \
      --customise Debug/HcalDebug/customize.analyze_raw_tp \
      --customise Debug/HcalDebug/customize.analyze_reemul_tp \
      --filein file:step2.root
