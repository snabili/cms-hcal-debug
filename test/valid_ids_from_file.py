import FWCore.ParameterSet.Config as cms

process = cms.Process("HFCALIB")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '76X_mcRun2_HeavyIon_v11'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/0C2BBF49-3387-E511-ADD5-3417EBE6459A.root')
)

process.load("Geometry.HcalCommonData.testPhase0GeometryXML_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("Configuration.Geometry.GeometryReco_cff")

# process.HcalTrigTowerGeometryESProducer.useFullGranularityHF = cms.bool( True )

process.es_ascii = cms.ESSource("HcalTextCalibrations",
    input = cms.VPSet(
        cms.PSet(
            object = cms.string('LutMetadata'),
            # full path: /afs/cern.ch/user/a/akhukhun/public/HF1x1TPs/LutMetadata_1x1.txt
            file = cms.FileInPath('LutMetadata_1x1.txt')
        )
    )
)
process.es_prefer_es_ascii = cms.ESPrefer("HcalTextCalibrations", "es_ascii")

process.validateIds = cms.EDAnalyzer("HcalValidDetIds")
process.p = cms.Path(process.validateIds)
