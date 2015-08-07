#!/usr/bin/env python
# vim: foldmethod=marker foldlevel=0

import FWCore.ParameterSet.Config as cms

process = cms.Process('test')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')

# from Configuration.AlCa.GlobalTag import GlobalTag as alcaGlobalTag
# from Configuration.AlCa.autoCond import conditions_L1_Run2012D as l1cond_raw
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'W19_500_62E2::All', '')

process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.Digi_cff')

# tag = "GR_R_53_V21::All" if data else "START53_V20::All"
# process.GlobalTag = alcaGlobalTag(process.GlobalTag,
        # globaltag=tag, conditions=l1cond)

process.chainplotter = cms.EDAnalyzer("HcalCompareLegacyChains",
        TriggerPrimitives = cms.InputTag('simHcalTriggerPrimitiveDigis'),
        RecHits = cms.InputTag('hbhereco'),
        DataFrames = cms.VInputTag(
                cms.InputTag("simHcalUnsuppressedDigis","HBHEDigiCollection","DIGI2RAW"),
                cms.InputTag("simHcalUnsuppressedDigis","HFDigiCollection","DIGI2RAW")
        )
)

process.p = cms.Path(process.chainplotter) # for plots

process.schedule = cms.Schedule(process.p)

process.TFileService = cms.Service("TFileService",
        closeFileFast = cms.untracked.bool(True),
        fileName = cms.string('legacy.root'))

process.source = cms.Source('PoolSource',
        fileNames =
        cms.untracked.vstring('file:10039_TTbar_14TeV+TTbar_Tauola_14TeV_2017_GenSimFull+DigiFull_2017+RecoFull_2017+HARVESTFull_2017/step3.root'),
        secondaryFileNames =
        cms.untracked.vstring('file:10039_TTbar_14TeV+TTbar_Tauola_14TeV_2017_GenSimFull+DigiFull_2017+RecoFull_2017+HARVESTFull_2017/step2.root'))
