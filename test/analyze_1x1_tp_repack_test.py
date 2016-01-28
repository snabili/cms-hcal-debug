import FWCore.ParameterSet.Config as cms

process = cms.Process("HFCALIB")

## Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_v0'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(

                # '/store/relval/CMSSW_8_0_0_pre4/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_asymptotic_v13-v1/00000/0CA0E567-BFA5-E511-9E3B-0025905B85F6.root',
                'file:0CA0E567-BFA5-E511-9E3B-0025905B85F6.root',

        )
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

# process.es_ascii = cms.ESSource("HcalTextCalibrations",
#     input = cms.VPSet(
#         cms.PSet(
#             object = cms.string('LutMetadata'),
#             # full path: /afs/cern.ch/user/a/akhukhun/public/HF1x1TPs/LutMetadata_1x1.txt
#             file = cms.FileInPath('LutMetadata_1x1.txt')
#         )
#     )
# )
# process.es_prefer_es_ascii = cms.ESPrefer("HcalTextCalibrations", "es_ascii")

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
process.es_prefer_es_pool = cms.ESPrefer( "PoolDBESSource", "es_pool" )

# process.load('CalibCalorimetry.HcalPlugins.Hcal_Conditions_forGlobalTag_cff')
# process.es_hardcode.toGet.append("LutMetadata")

process.TFileService = cms.Service("TFileService",
        closeFileFast = cms.untracked.bool(True),
        fileName = cms.string('analyze.root'))

process.load("EventFilter.HcalRawToDigi.HcalDigiToRaw_cfi")
process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")

process.hcalRawData.TRIG = cms.untracked.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB")
process.hcalDigis.InputLabel = cms.InputTag("hcalRawData")

process.analyzeOutput = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("hcalDigis", "" , "HFCALIB"))
process.analyzeInput = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"))

process.p = cms.Path(
        process.simHcalTriggerPrimitiveDigis
        * process.hcalRawData
        * process.hcalDigis
        * process.analyzeInput
        * process.analyzeOutput
)

# print process.dumpPython()
