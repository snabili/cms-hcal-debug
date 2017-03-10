import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.autoCond import autoCond
from CondCore.CondDB.CondDB_cfi import CondDB
from SLHCUpgradeSimulations.Configuration.HCalCustoms import customise_HcalPhase1

process = cms.Process("HFCALIB")

# Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = autoCond['run2_mc']

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1000))

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(

        '/store/mc/RunIISpring16DR80/SinglePiMinus_E1to1000_Eta5p2_13TeV_FlatRandom/GEN-SIM-RAW/NoPURAW_NZS_withHLT_80X_mcRun2_asymptotic_v14-v1/40000/0222BED7-F25B-E611-9CA3-1CC1DE18CFF0.root',

    )
)

# process.out = cms.OutputModule( "PoolOutputModule",
#         fileName = cms.untracked.string("output.root"),
#         outputCommands = cms.untracked.vstring( 'keep *' )
# )
# process.end = cms.EndPath(process.out)

process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag(cms.InputTag('simHcalUnsuppressedDigis'), cms.InputTag('simHcalUnsuppressedDigis'))
# process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('simHcalDigis'), cms.InputTag('simHcalDigis') )
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)

process.load("Configuration.Geometry.GeometryExtended2016Reco_cff")

process.CondDBSetup = CondDB.clone()
delattr(process.CondDBSetup, 'connect')

# process.es_pool = cms.ESSource("PoolDBESSource",
#      process.CondDBSetup,
#      timetype = cms.string('runnumber'),
#      toGet = cms.VPSet(
#          cms.PSet(record = cms.string("HcalLutMetadataRcd"),
#              tag = cms.string("HcalLutMetadata_HFTP_1x1")
#              ),
#          cms.PSet(record = cms.string("HcalElectronicsMapRcd"),
#              tag = cms.string("HcalElectronicsMap_HFTP_1x1")
#              )
#          ),
#      connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
#      authenticationMethod = cms.untracked.uint32(0)
#      )
# process.es_prefer_es_pool = cms.ESPrefer( "PoolDBESSource", "es_pool" )

# customise_HcalPhase1(process)

process.TFileService = cms.Service("TFileService",
                                   closeFileFast=cms.untracked.bool(True),
                                   fileName=cms.string('analyze.root'))

process.analyze = cms.EDAnalyzer("AnalyzeTP",
                                 triggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis", "", "HFCALIB"))
process.analyzeOld = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis", "", "HLT"))
process.analyzeRaw = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("hcalDigis", "", ""))
process.compare = cms.EDAnalyzer("CompareTP",
                                 triggerPrimitives=cms.InputTag(
                                     "simHcalTriggerPrimitiveDigis", "", "HLT"),
                                 emulTriggerPrimitives=cms.InputTag(
                                     "simHcalTriggerPrimitiveDigis", "", "HFCALIB"),
                                 swapIphi=cms.bool(False))

process.p = cms.Path(process.hcalDigis *
                     process.simHcalTriggerPrimitiveDigis *
                     process.analyze *
                     process.analyzeOld *
                     process.analyzeRaw *
                     process.compare)

# print process.dumpPython()
