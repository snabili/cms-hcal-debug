import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PLOT', eras.Run2_2017)

# Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
# process.GlobalTag.globaltag = '80X_dataRun2_Express_v12'

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(10000))

process.source = cms.Source(
    "HcalTBSource",
    fileNames=cms.untracked.vstring(
        '/store/group/dpg_hcal/comm_hcal/USC/run283104/USC_283104.root'
    )
)

process.load("CondCore.CondDB.CondDB_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")

# process.es_pool = cms.ESSource(
#     "PoolDBESSource",
#     process.CondDBSetup,
#     timetype=cms.string('runnumber'),
#     toGet=cms.VPSet(
#         cms.PSet(record=cms.string(
#             "HcalL1TriggerObjectsRcd"),
#             tag=cms.string("HcalL1TriggerObjects_Physics2016v5B38T")
#         )
#     ),
#     connect=cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
#     authenticationMethod=cms.untracked.uint32(0)
# )
# process.es_prefer_es_pool = cms.ESPrefer("PoolDBESSource", "es_pool")

process.es_ascii = cms.ESSource(
    'HcalTextCalibrations',
    input=cms.VPSet(
        cms.PSet(
            object=cms.string('ElectronicsMap'),
            file=cms.FileInPath('Debug/HcalDebug/test/version_G_emap_all_ngHF2016.txt')
        )
    )
)
process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')


# process.load('CalibCalorimetry.HcalPlugins.Hcal_Conditions_forGlobalTag_cff')
# process.load('Configuration.Geometry.GeometryExtended2017newReco_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')

# process.es_hardcode.toGet.append("HcalTPParametersRcd")

process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag(
    cms.InputTag('hcalDigis'), cms.InputTag('hcalDigis'))
process.simHcalTriggerPrimitiveDigis.inputUpgradeLabel = cms.VInputTag(
    cms.InputTag('hcalDigis'), cms.InputTag('hcalDigis'))
process.simHcalTriggerPrimitiveDigis.parameters = cms.untracked.PSet(
    ADCThresholdHF=cms.uint32(1024))
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)
process.simHcalTriggerPrimitiveDigis.upgradeHF = cms.bool(True)

process.TFileService = cms.Service("TFileService",
                                   closeFileFast=cms.untracked.bool(True),
                                   fileName=cms.string('analyze.root'))

process.hcalDigis.InputLabel = cms.InputTag("source")
process.analyzeRAW = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("hcalDigis", "", ""))
process.analyzeSIM = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis", "", ""))
process.compare = cms.EDAnalyzer("CompareTP",
                                 triggerPrimitives=cms.InputTag("hcalDigis"),
                                 emulTriggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis"),
                                 swapIphi=cms.bool(False))

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.hcalDigis *
    # process.dump *
    process.simHcalTriggerPrimitiveDigis *
    process.analyzeRAW *
    process.analyzeSIM *
    process.compare)

# print process.dumpPython()
