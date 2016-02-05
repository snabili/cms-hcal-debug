import FWCore.ParameterSet.Config as cms

process = cms.Process("HFCALIB")

## Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_v0'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(

                '/store/relval/CMSSW_7_6_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_asymptotic_v11-v1/00000/4431031E-9E7F-E511-9F42-0025905938A4.root',
                '/store/relval/CMSSW_7_6_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_asymptotic_v11-v1/00000/4C462F65-9F7F-E511-972A-0026189438A9.root',
                '/store/relval/CMSSW_7_6_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_asymptotic_v11-v1/00000/703E7EAB-9D7F-E511-B886-003048FFCBFC.root',
                '/store/relval/CMSSW_7_6_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_asymptotic_v11-v1/00000/8AF07AAB-9D7F-E511-B8B4-003048FFCBFC.root',
                '/store/relval/CMSSW_7_6_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_asymptotic_v11-v1/00000/962BEF7C-9D7F-E511-A2BB-0025905B85AA.root',
                '/store/relval/CMSSW_7_6_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_asymptotic_v11-v1/00000/C409A519-9E7F-E511-BD4C-0025905B8590.root',
                '/store/relval/CMSSW_7_6_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_asymptotic_v11-v1/00000/E8D41D6A-9F7F-E511-A10A-003048FFD740.root',
                '/store/relval/CMSSW_7_6_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_asymptotic_v11-v1/00000/EE048767-9E7F-E511-B1AA-0025905B8606.root',

        )
)

# process.out = cms.OutputModule( "PoolOutputModule",
#         fileName = cms.untracked.string("output.root"),
#         outputCommands = cms.untracked.vstring( 'keep *' )
# )
# process.end = cms.EndPath(process.out)

process.load('L1Trigger.RegionalCaloTrigger.rctDigis_cfi')
process.rctDigis.hcalDigis = cms.VInputTag(cms.InputTag("simHcalTriggerPrimitiveDigis"))

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('simHcalUnsuppressedDigis'), cms.InputTag('simHcalUnsuppressedDigis') )
# process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('simHcalDigis'), cms.InputTag('simHcalDigis') )
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)

process.load("Configuration.Geometry.GeometryExtended2016Reco_cff")
process.XMLIdealGeometryESSource.geomXMLFiles.remove('Geometry/HcalCommonData/data/Phase0/hcalRecNumbering.xml')
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/HcalCommonData/data/Phase0/hcalRecNumberingRun2.xml')

process.es_pool = cms.ESSource("PoolDBESSource",
     process.CondDBSetup,
     timetype = cms.string('runnumber'),
     toGet = cms.VPSet(
         cms.PSet(record = cms.string("HcalLutMetadataRcd"),
             tag = cms.string("HcalLutMetadata_HFTP_1x1")
             )
         ),
     connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
     authenticationMethod = cms.untracked.uint32(0)
     )
process.es_prefer_es_pool = cms.ESPrefer( "PoolDBESSource", "es_pool" )

process.TFileService = cms.Service("TFileService",
        closeFileFast = cms.untracked.bool(True),
        fileName = cms.string('analyze.root'))

process.analyze = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"))
process.analyzeOld = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HLT"))

process.p = cms.Path(process.simHcalTriggerPrimitiveDigis * process.analyze * process.analyzeOld)

# print process.dumpPython()
