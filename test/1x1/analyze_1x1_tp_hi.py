import FWCore.ParameterSet.Config as cms

process = cms.Process("HFCALIB")

## Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '80X_mcRun2_HeavyIon_v0'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
                '/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/0C2BBF49-3387-E511-ADD5-3417EBE6459A.root',
                '/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/200C7831-3987-E511-944F-00266CF20468.root',
                '/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/3849AE0A-3487-E511-90EF-00266CFAE228.root',
                '/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/40550625-3087-E511-A4C4-00266CFAE764.root',
                '/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/500B67F4-3487-E511-B068-00266CF2506C.root',
                '/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/5267F176-3887-E511-927F-00266CFADEC0.root',
                '/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/5EF74DD0-3187-E511-B2BF-00266CF9BEF8.root',
                '/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/6C09D335-4087-E511-AFE5-00266CF9B9F0.root',
                '/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/6EBBD119-3187-E511-ACB3-00266CF9B828.root',
                '/store/relval/CMSSW_8_0_0_pre1/RelValHydjetQ_MinBias_5020GeV/GEN-SIM-DIGI-RAW-HLTDEBUG/76X_mcRun2_HeavyIon_v11_resub-v1/00000/86C86721-3787-E511-A445-00A0D1EE271C.root',
        )
)

# process.out = cms.OutputModule( "PoolOutputModule",
#         fileName = cms.untracked.string("output_hi.root"),
#         outputCommands = cms.untracked.vstring( 'keep *' )
# )
# process.end = cms.EndPath(process.out)

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

process.load('CalibCalorimetry.HcalPlugins.Hcal_Conditions_forGlobalTag_cff')
process.es_hardcode.toGet.append("LutMetadata")

process.TFileService = cms.Service("TFileService",
        closeFileFast = cms.untracked.bool(True),
        fileName = cms.string('analyze_hi.root'))

process.analyze = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HFCALIB"))
process.analyzeOld = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , "HLT"))

process.p = cms.Path(process.simHcalTriggerPrimitiveDigis * process.analyze * process.analyzeOld)

# print process.dumpPython()
