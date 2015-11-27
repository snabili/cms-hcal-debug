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

process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.Digi_cff')

# tag = "GR_R_53_V21::All" if data else "START53_V20::All"
# process.GlobalTag = alcaGlobalTag(process.GlobalTag,
        # globaltag=tag, conditions=l1cond)

from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2019 

# call to customisation function cust_2019 imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
# process = cust_2019(process)

process.load("CalibCalorimetry/HcalPlugins/Hcal_Conditions_forGlobalTag_cff")
process.es_hardcode.toGet = cms.untracked.vstring(
            'GainWidths',
            'MCParams',
            'RecoParams',
            'RespCorrs',
            'QIEData',
            'Gains',
            'Pedestals',
            'PedestalWidths',
            'ChannelQuality',
            'ZSThresholds',
            'TimeCorrs',
            'LUTCorrs',
            'LutMetadata',
            'L1TriggerObjects',
            'PFCorrs',
            'ElectronicsMap',
            'CholeskyMatrices',
            'CovarianceMatrices'
            )

# process.es_hardcode.hcalTopologyConstants.mode=cms.string('HcalTopologyMode::SLHC')
# process.es_hardcode.hcalTopologyConstants.maxDepthHB=cms.int32(3)
# process.es_hardcode.hcalTopologyConstants.maxDepthHB=cms.int32(3)
# process.es_hardcode.hcalTopologyConstants.maxDepthHE=cms.int32(5)
# process.es_hardcode.HcalReLabel.RelabelHits=cms.untracked.bool(True)
# Special Upgrade trick (if absent - regular case assumed)
process.es_hardcode.GainWidthsForTrigPrims = cms.bool(True)
process.es_hardcode.HEreCalibCutoff = cms.double(100.) #for aging

# process.hcalTopologyIdeal.hcalTopologyConstants.mode=cms.string('HcalTopologyMode::SLHC')
# process.hcalTopologyIdeal.hcalTopologyConstants.maxDepthHB=cms.int32(3)
# process.hcalTopologyIdeal.hcalTopologyConstants.maxDepthHE=cms.int32(5)

if hasattr(process,'HcalTPGCoderULUT'):
    # print "MODIFY"
    process.HcalTPGCoderULUT.upgrade = cms.bool(True)
    process.CaloTPGTranscoder.upgrade = cms.bool(True)
    # process.HcalTPGCoderULUT.hcalTopologyConstants.mode=cms.string('HcalTopologyMode::SLHC')
    # process.HcalTPGCoderULUT.hcalTopologyConstants.maxDepthHB=cms.int32(3)
    # process.HcalTPGCoderULUT.hcalTopologyConstants.maxDepthHE=cms.int32(5)

process.chainplotter = cms.EDAnalyzer("HcalDebug",
        TriggerPrimitives = cms.InputTag('simHcalTriggerPrimitiveDigis'),
        RecHits = cms.InputTag('hbheUpgradeReco'),
        DataFrames = cms.VInputTag(
                cms.InputTag("simHcalUnsuppressedDigis","HBHEUpgradeDigiCollection","DIGI2RAW"),
                cms.InputTag("simHcalUnsuppressedDigis","HFUpgradeDigiCollection","DIGI2RAW")
        )
)

process.p = cms.Path(process.chainplotter) # for plots

process.schedule = cms.Schedule(process.p)

process.TFileService = cms.Service("TFileService",
        closeFileFast = cms.untracked.bool(True),
        fileName = cms.string('upgrade.root'))

process.source = cms.Source('PoolSource',
        fileNames =
        cms.untracked.vstring('file:10239_TTbar_14TeV+TTbar_Tauola_14TeV_2019_GenSimFull+DigiFull_2019+RecoFull_2019+HARVESTFull_2019/step3.root'),
        secondaryFileNames =
        cms.untracked.vstring('file:10239_TTbar_14TeV+TTbar_Tauola_14TeV_2019_GenSimFull+DigiFull_2019+RecoFull_2019+HARVESTFull_2019/step2.root'))
