#!/bin/sh

run_ttbar() {
   dir=$1
   conditions=$2
   era=$3
   geometry=$4
   events=50

   [ -n "$geometry" ] && geometry="--geometry $geometry"

   [ ! -d $dir ] && mkdir $dir

   set -e
   set -x

   cd $dir

   cmsDriver.py TTbar_13TeV_TuneCUETP8M1_cfi  --conditions auto:$conditions -n $events --era $era --eventcontent FEVTDEBUG --relval 9000,50 -s GEN,SIM --datatier GEN-SIM --beamspot Realistic50ns13TeVCollision --fileout=step1.root $geometry &> step1.log
   cmsDriver.py step2  --conditions auto:$conditions -s DIGI:pdigi_valid,L1,DIGI2RAW,HLT:@fake --datatier GEN-SIM-DIGI-RAW -n $events $geometry --era $era --eventcontent FEVTDEBUGHLT --fileout=step2.root --filein=file:step1.root &> step2.log
   cmsDriver.py step3  --conditions auto:$conditions -n $events --era $era --eventcontent RECOSIM,MINIAODSIM,DQM --runUnscheduled  -s RAW2DIGI,L1Reco,RECO,EI,PAT,VALIDATION:@standardValidation+@miniAODValidation,DQM:@standardDQM+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --fileout=step3.root --filein=file:step2.root $geometry &> step3.log
}

run_ttbar 2017_HCAL phase1_2017_realistic Run2_2017 Extended2017new &
run_ttbar 2016_plain run2_mc Run2_2016 "" &

wait
