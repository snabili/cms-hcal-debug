import FWCore.ParameterSet.Config as cms

process = cms.Process("HFCALIB")

## Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['run2_data']

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
        skipBadFiles = cms.untracked.bool(True),
        fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring(),
        firstRun = cms.untracked.uint32(260627)
        # fileNames = cms.untracked.vstring('/store/data/Run2015E/HIEWQExo/RAW/v1/000/262/219/00000/1E4169BC-4891-E511-8D99-02163E0146CF.root')
)

process.source.fileNames.extend([
    '/store/data/Run2015D/HcalNZS/RECO/PromptReco-v4/000/260/627/00000/842C5A38-8884-E511-B203-02163E011CD0.root'
])
process.source.secondaryFileNames.extend([
    '/store/data/Run2015D/HcalNZS/RAW/v1/000/260/627/00000/22600A0D-0C82-E511-932F-02163E014276.root'
])

# process.out = cms.OutputModule( "PoolOutputModule",
#         fileName = cms.untracked.string("output_data.root"),
#         outputCommands = cms.untracked.vstring( 'keep *' )
# )
# process.end = cms.EndPath(process.out)

process.load('L1Trigger.RegionalCaloTrigger.rctDigis_cfi')
process.rctDigis.hcalDigis = cms.VInputTag(cms.InputTag("simHcalTriggerPrimitiveDigis"))

process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
# process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('hcalDigis'), cms.InputTag('hcalDigis') )
# process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('simHcalDigis'), cms.InputTag('simHcalDigis') )
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)

process.load("Configuration.Geometry.GeometryExtended2016Reco_cff")

# process.raw2digi_step = cms.Path(process.hcalDigis)
# from SLHCUpgradeSimulations.Configuration.HCalCustoms import customise_HcalPhase1
# customise_HcalPhase1(process)

process.es_pool = cms.ESSource("PoolDBESSource",
     process.CondDBSetup,
     timetype = cms.string('runnumber'),
     toGet = cms.VPSet(
         cms.PSet(record = cms.string("HcalLutMetadataRcd"),
             tag = cms.string("HcalLutMetadata_HFTP_1x1")
             ),
         cms.PSet(record = cms.string("HcalElectronicsMapRcd"),
             tag = cms.string("HcalElectronicsMap_HFTP_1x1")
             )
         ),
     connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
     authenticationMethod = cms.untracked.uint32(0)
     )
# process.es_hardcode.toGet.remove("LutMetadata")
# process.es_hardcode.toGet.remove("ElectronicsMap")
process.es_prefer_es_pool = cms.ESPrefer("PoolDBESSource", "es_pool")

process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

process.TFileService = cms.Service("TFileService",
        closeFileFast = cms.untracked.bool(True),
        fileName = cms.string('analyze_data.root'))

process.analyze = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"))
process.analyzeRaw = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("hcalDigis", "" , "HFCALIB"))
process.compare = cms.EDAnalyzer("CompareTP",
        triggerPrimitives = cms.InputTag("hcalDigis", "" , "HFCALIB"),
        emulTriggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"),
        swapIphi = cms.bool(False))

process.chainplotter = cms.EDAnalyzer("HcalCompareLegacyChains",
        triggerPrimitives = cms.InputTag('hcalDigis', '', 'HFCALIB'),
        recHits = cms.VInputTag('hbhereco', 'hfreco'),
        dataFrames = cms.VInputTag(
                cms.InputTag("hcalDigis", "", "HFCALIB"),
                cms.InputTag("hcalDigis", "", "HFCALIB")
        ),
        swapIphi = cms.bool(False)
)

process.compareSwapped = process.compare.clone()
process.compareSwapped.swapIphi = cms.bool(True)
process.chainplotterSwapped = process.chainplotter.clone()
process.chainplotterSwapped.swapIphi = cms.bool(True)

# process.hcalDigis.InputLabel = cms.InputTag("rawDataRepacker")

process.p = cms.Path(
        process.hcalDigis
        * process.simHcalTriggerPrimitiveDigis 
        * process.analyze
        * process.analyzeRaw
        * process.compare
        * process.compareSwapped
        * process.chainplotter
        * process.chainplotterSwapped
)

# print process.dumpPython()
