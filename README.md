# Setup

Install with:

    git clone git@github.com:cms-hcal-trigger/cms-hcal-debug.git Debug/HcalDebug
    scram b -j 8

# Examples

## With Workflows from `runTheMatrix.py`

Run with:

    runTheMatrix.py -w upgrade -l 10039
    cmsRun Debug/HcalDebug/test/cmp_legacy.py

Change to the output directory and then analyze the second step:

    cmsDriver.py analyze \
      --conditions auto:phase1_2017_realistic \
      -s RAW2DIGI,DIGI --geometry DB:Extended --era Run2_2017 \
      --customise Debug/HcalDebug/customize.analyze_raw_tp \
      --customise Debug/HcalDebug/customize.analyze_reemul_tp \
      --filein file:step2.root \
      -n 10

## Datasets from DAS

Use as input to the `cmsDriver.py` command:

    cmsDriver.py analyze \
      --conditions auto:phase1_2017_realistic \
      -s RAW2DIGI,DIGI --geometry DB:Extended --era Run2_2017 \
      --customise Debug/HcalDebug/customize.analyze_raw_tp \
      --customise Debug/HcalDebug/customize.analyze_reemul_tp \
      --filein das:/RelValTTbarLepton_13/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v15-v1/GEN-SIM-DIGI-RAW \
      -n 1000

## From Data, Using L1T Digis

Using a run with HF FG bit mis-matches between L1T inputs (HCAL RAW does
not include FG bits) and re-emulation:

    cmsDriver.py analyze \
      --data --conditions auto:run2_data \
      -s RAW2DIGI --geometry DB:Extended --era Run2_2017 \
      --customise Debug/HcalDebug/customize.analyze_l1t_tp \
      --customise Debug/HcalDebug/customize.analyze_raw_tp \
      --customise Debug/HcalDebug/customize.analyze_reemul_tp \
      --customise Debug/HcalDebug/customize.compare_l1t_reemul_tp \
      --customise Debug/HcalDebug/customize.use_data_reemul_tp \
      --filein /store/data/Run2017C/HcalNZS/RAW/v1/000/299/844/00000/AE36B18A-5271-E711-A223-02163E013895.root,/store/data/Run2017C/HcalNZS/RAW/v1/000/299/844/00000/46B78BA1-5271-E711-8820-02163E01A60E.root \
      -n -1

## From Data, Using L1T Digis and comparing with RecHits

As before, but using files to contain primary and secondary inputs, and
adding TriggerPrimitive to RecHit comparisons:

    cmsDriver.py analyze \
      --data --conditions 92X_dataRun2_Prompt_v8 \
      -s RAW2DIGI --geometry DB:Extended --era Run2_2017 \
      --no_output \
      --customise Debug/HcalDebug/customize.analyze_l1t_tp \
      --customise Debug/HcalDebug/customize.analyze_raw_tp \
      --customise Debug/HcalDebug/customize.analyze_reemul_tp \
      --customise Debug/HcalDebug/customize.compare_l1t_reemul_tp \
      --customise Debug/HcalDebug/customize.compare_raw_reco_sev9 \
      --customise Debug/HcalDebug/customize.compare_raw_reco_sev9999 \
      --customise Debug/HcalDebug/customize.use_data_reemul_tp \
      --filein=$(<~/JetHTRECO.txt) \
      --secondfilein=$(<~/JetHT.txt) \
      -n 50000
