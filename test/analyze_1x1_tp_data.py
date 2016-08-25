import FWCore.ParameterSet.Config as cms

process = cms.Process("HFCALIB")

## Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['run2_data']
process.GlobalTag.globaltag = "80X_dataRun2_HLT_v6"
print process.GlobalTag.globaltag

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
# process.skipEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

lst = []
process.source = cms.Source("PoolSource",
        skipBadFiles = cms.untracked.bool(True),
        fileNames = cms.untracked.vstring(lst),
        # firstRun = cms.untracked.uint32(272818)
        # fileNames = cms.untracked.vstring('/store/data/Run2015E/HIEWQExo/RAW/v1/000/262/219/00000/1E4169BC-4891-E511-8D99-02163E0146CF.root')
)

process.source.fileNames.extend([

    '/store/data/Run2016F/JetHT/RAW/v1/000/278/808/00001/FEE2A608-3462-E611-A66E-02163E01370F.root'

])

# process.out = cms.OutputModule( "PoolOutputModule",
#         fileName = cms.untracked.string("output_data.root"),
#         outputCommands = cms.untracked.vstring( 'keep *' )
# )
# process.end = cms.EndPath(process.out)

process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
# process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('hcalDigis'), cms.InputTag('hcalDigis') )
# process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('simHcalDigis'), cms.InputTag('simHcalDigis') )
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)

process.load("Configuration.Geometry.GeometryExtended2016Reco_cff")

process.load("EventFilter.L1TRawToDigi.caloStage2Digis_cfi")
process.load("EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi")

# process.raw2digi_step = cms.Path(process.hcalDigis)
# from SLHCUpgradeSimulations.Configuration.HCalCustoms import customise_HcalPhase1
# customise_HcalPhase1(process)

# process.es_pool = cms.ESSource("PoolDBESSource",
#      process.CondDBSetup,
#      timetype = cms.string('runnumber'),
#      toGet = cms.VPSet(
#          # cms.PSet(record = cms.string("HcalLutMetadataRcd"),
#          #     tag = cms.string("HcalLutMetadata_HFTP_1x1")
#          #     ),
#          cms.PSet(record = cms.string("HcalElectronicsMapRcd"),
#              tag = cms.string("HcalElectronicsMap_HFTP_1x1")
#              )
#          ),
#      connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
#      authenticationMethod = cms.untracked.uint32(0)
#      )
# # process.es_hardcode.toGet.remove("LutMetadata")
# process.es_hardcode.toGet.remove("ElectronicsMap")
# process.es_prefer_es_pool = cms.ESPrefer("PoolDBESSource", "es_pool")

process.TFileService = cms.Service("TFileService",
        closeFileFast = cms.untracked.bool(True),
        fileName = cms.string('analyze_data.root'))

process.analyze = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"))
process.analyzeL1T = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("l1tCaloLayer1Digis", "" , "HFCALIB"))
process.analyzeRaw = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("hcalDigis", "" , "HFCALIB"))
process.compare = cms.EDAnalyzer("CompareTP",
        swapIphi = cms.bool(False),
        triggerPrimitives = cms.InputTag("hcalDigis", "" , "HFCALIB"),
        emulTriggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"))
process.compareL1T = cms.EDAnalyzer("CompareTP",
        swapIphi = cms.bool(False),
        triggerPrimitives = cms.InputTag("l1tCaloLayer1Digis", "" , "HFCALIB"),
        emulTriggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"))
process.analyzeCT = cms.EDAnalyzer("AnalyzeCT",
        caloTowers = cms.InputTag("caloStage2Digis", "CaloTower"))

# process.hcalDigis.InputLabel = cms.InputTag("rawDataRepacker")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
        process.hcalDigis
        * process.l1tCaloLayer1Digis
        * process.caloStage2Digis
        * process.simHcalTriggerPrimitiveDigis
        # * process.dump
        * process.analyze
        * process.analyzeL1T
        * process.analyzeRaw
        * process.analyzeCT
        * process.compare
        * process.compareL1T
)

# print process.dumpPython()
