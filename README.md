# Setup

Install with:

    git clone git@github.com:matz-e/cms-hcal-debug.git Debug/HcalDebug
    scram b -j 8

# Running

Run with:

    runTheMatrix.py -w upgrade -l 10039
    cmsRun Debug/HcalDebug/test/cmp_legacy.py
