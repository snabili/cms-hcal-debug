import FWCore.ParameterSet.Config as cms

process = cms.Process("HFCALIB")

## Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['run2_data']
# process.GlobalTag.globaltag = "80X_dataRun2_HLT_Validation_HcalL1TriggerObjects_38T_week24_2016"
print process.GlobalTag.globaltag

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
# process.skipEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

lst = []
process.source = cms.Source("PoolSource",
        skipBadFiles = cms.untracked.bool(True),
        fileNames = cms.untracked.vstring(lst),
        firstRun = cms.untracked.uint32(275783)
)

process.source.fileNames.extend([
    '/store/data/Run2016C/HcalNZS/RAW/v2/000/275/783/00000/A491D9C5-433B-E611-BA33-02163E011A20.root',
    '/store/data/Run2016C/HcalNZS/RAW/v2/000/275/783/00000/BC0B523A-743B-E611-BFDD-02163E012443.root',
])

process.load("Configuration.Geometry.GeometryExtended2016Reco_cff")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")

# from SLHCUpgradeSimulations.Configuration.HCalCustoms import customise_Hcal2016
# customise_Hcal2016(process)

process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('hcalDigis'), cms.InputTag('hcalDigis') )
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)

process.TFileService = cms.Service("TFileService",
        closeFileFast = cms.untracked.bool(True),
        fileName = cms.string('analyze_data.root'))

process.analyze = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"))
process.analyzeRaw = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("hcalDigis", "" , "HFCALIB"))
process.compare = cms.EDAnalyzer("CompareTP",
        swapIphi = cms.bool(False),
        triggerPrimitives = cms.InputTag("hcalDigis", "" , "HFCALIB"),
        emulTriggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"))

# process.hcalDigis.InputLabel = cms.InputTag("rawDataRepacker")

process.p = cms.Path(
        process.hcalDigis
        * process.simHcalTriggerPrimitiveDigis
        * process.analyze
        * process.analyzeRaw
        * process.compare
)
