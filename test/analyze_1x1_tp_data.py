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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
# process.skipEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

lst = []
process.source = cms.Source("PoolSource",
        skipBadFiles = cms.untracked.bool(True),
        fileNames = cms.untracked.vstring(lst),
        # firstRun = cms.untracked.uint32(272818)
        # fileNames = cms.untracked.vstring('/store/data/Run2015E/HIEWQExo/RAW/v1/000/262/219/00000/1E4169BC-4891-E511-8D99-02163E0146CF.root')
)

process.source.fileNames.extend([

    # '/store/data/Run2016B/HcalNZS/RAW/v1/000/272/554/00000/08ECD81D-6112-E611-9E5B-02163E0145F3.root'
    # '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/10E13AC8-7418-E611-818E-02163E01425B.root'
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/10E13AC8-7418-E611-818E-02163E01425B.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/24189AC8-7418-E611-9009-02163E014243.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/24F836BF-7418-E611-8736-02163E014475.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/3A8E5DBF-7418-E611-A74D-02163E01420A.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/3CF28DBF-7418-E611-B347-02163E0128C4.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/422D55D0-7418-E611-8521-02163E0142BE.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/74EE8EBC-7418-E611-8441-02163E014336.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/BEBBA81C-7518-E611-94D1-02163E014105.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/CA46F5ED-7418-E611-AFEA-02163E01438B.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/CCDBB4E2-7418-E611-BDB6-02163E0136C7.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/EADA80CB-7F18-E611-9F1A-02163E0142E8.root',
    '/store/data/Run2016B/HcalNZS/RAW/v2/000/273/158/00000/F27360D3-7418-E611-AB7C-02163E013653.root',

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
process.analyzeRaw = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("hcalDigis", "" , "HFCALIB"))
process.compare = cms.EDAnalyzer("CompareTP",
        swapIphi = cms.bool(False),
        triggerPrimitives = cms.InputTag("hcalDigis", "" , "HFCALIB"),
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
        * process.analyzeRaw
        * process.analyzeCT
        * process.compare
)

# print process.dumpPython()
