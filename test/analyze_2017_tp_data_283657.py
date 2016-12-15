import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PLOT', eras.Run2_2017)

# Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
process.GlobalTag.globaltag = '81X_upgrade2017_HCALdev_v1'

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(52))

process.source = cms.Source(
    "HcalTBSource",
    fileNames=cms.untracked.vstring(
        '/store/group/dpg_hcal/comm_hcal/USC/run283657/USC_283657.root'
    ),
    skipEvents=cms.untracked.uint32(51)
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

process.es_pool = cms.ESSource("PoolDBESSource",
                               #    process.CondDB,
                               timetype=cms.string("runnumber"),
                               connect=cms.string("sqlite:HcalCond_data_test_2016.db"),
                               authenticationMethod=cms.untracked.uint32(0),
                               toGet=cms.VPSet(
                                   cms.PSet(
                                       record=cms.string("HcalRecoParamsRcd"),
                                       tag=cms.string("HcalRecoParams")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalMCParamsRcd"),
                                       tag=cms.string("HcalMCParams")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalTimeCorrsRcd"),
                                       tag=cms.string("HcalTimeCorrs")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalQIETypesRcd"),
                                       tag=cms.string("HcalQIETypes")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalLongRecoParamsRcd"),
                                       tag=cms.string("HcalLongRecoParams")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalPFCorrsRcd"),
                                       tag=cms.string("HcalPFCorrs")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalFlagHFDigiTimeParamsRcd"),
                                       tag=cms.string("HcalFlagHFDigiTimeParams")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalChannelQualityRcd"),
                                       tag=cms.string("HcalChannelQuality")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalPedestalWidthsRcd"),
                                       tag=cms.string("HcalPedestalWidths")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalPedestalsRcd"),
                                       tag=cms.string("HcalPedestals")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalQIEDataRcd"),
                                       tag=cms.string("HcalQIEData")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalRespCorrsRcd"),
                                       tag=cms.string("HcalRespCorrs")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalGainsRcd"),
                                       tag=cms.string("HcalGains")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalLutMetadataRcd"),
                                       tag=cms.string("HcalLutMetadata")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalL1TriggerObjectsRcd"),
                                       tag=cms.string("HcalL1TriggerObjects")
                                   ),
                                   cms.PSet(
                                       record=cms.string("HcalZSThresholdsRcd"),
                                       tag=cms.string("HcalZSThresholds")
                                   ),
                                   # cms.PSet(
                                   #          record = cms.string("HcalElectronicsMapRcd"),
                                   #          tag = cms.string("HcalElectronicsMap")
                                   #          ),
                                   #                      cms.PSet(
                                   #                               record = cms.string("HcalGainWidthsRcd"),
                                   #                               tag = cms.string("HcalGainWidths_mc_test_2016")
                                   #                               ),
                                   cms.PSet(
                                       record=cms.string("HcalLUTCorrsRcd"),
                                       tag=cms.string("HcalLUTCorrs")
                                   )
                               ))
process.es_prefer_es_pool = cms.ESPrefer("PoolDBESSource", "es_pool")


# process.load('CalibCalorimetry.HcalPlugins.Hcal_Conditions_forGlobalTag_cff')
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2016dev_cff')
# process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
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
    TDCMaskHF=cms.uint64(0x8003FFFFFFFFFFFF))
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)
process.simHcalTriggerPrimitiveDigis.upgradeHF = cms.bool(True)

process.TFileService = cms.Service("TFileService",
                                   closeFileFast=cms.untracked.bool(True),
                                   fileName=cms.string('analyze_283657.root'))

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
