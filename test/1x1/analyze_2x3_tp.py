import FWCore.ParameterSet.Config as cms

process = cms.Process("HFCALIB")

## Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = 'PRE_SHI72_V7'
process.GlobalTag.globaltag = '76X_mcRun2_HeavyIon_v11'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/0C2BBF49-3387-E511-ADD5-3417EBE6459A.root')
)

process.out = cms.OutputModule( "PoolOutputModule",
        fileName = cms.untracked.string("output.root"),
        outputCommands = cms.untracked.vstring( 'keep *' )
)
process.end = cms.EndPath(process.out)

process.load('L1Trigger.RegionalCaloTrigger.rctDigis_cfi')
process.rctDigis.hcalDigis = cms.VInputTag(cms.InputTag("simHcalTriggerPrimitiveDigis"))

process.load("Geometry.HcalCommonData.testPhase0GeometryXML_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("Configuration.Geometry.GeometryReco_cff")

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('simHcalUnsuppressedDigis'), cms.InputTag('simHcalUnsuppressedDigis') )
# process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('simHcalDigis'), cms.InputTag('simHcalDigis') )
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)

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

process.TFileService = cms.Service("TFileService",
        closeFileFast = cms.untracked.bool(True),
        fileName = cms.string('analyze23.root'))

process.analyze23 = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"))
process.analyzeOld23 = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HLT"))

process.p = cms.Path(process.simHcalTriggerPrimitiveDigis * process.analyze23 * process.analyzeOld23)

# print process.dumpPython()
