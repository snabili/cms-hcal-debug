import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.Utilities import convertToUnscheduled

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
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(10000))

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(

        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/A6ED3F67-3877-E611-AAAA-02163E0146FA.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/C2233C55-3877-E611-8494-02163E012213.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/401D9061-3877-E611-AA18-02163E01391F.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/BE02DB35-3877-E611-9CB3-02163E011F15.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/54F2571F-3877-E611-8B6C-FA163EC64E27.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/28D9326B-3877-E611-A9B8-02163E01458B.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/987EA939-3877-E611-A6C3-02163E01374D.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/C298E838-3877-E611-AB1F-02163E011885.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/14714179-3877-E611-9083-02163E0143F9.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/2CD0F07A-3877-E611-BC34-02163E01466F.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/6C76B06B-3877-E611-92F1-02163E01382C.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/4CE07138-3877-E611-94F1-02163E0119F4.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/80F66662-3877-E611-8F71-02163E014181.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/EAC82B3E-3877-E611-B97C-02163E0139D5.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/CCAA918C-3877-E611-AEE4-02163E0142D9.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/D21B5467-3877-E611-9F07-02163E01286F.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/DE5C3D45-3877-E611-817D-02163E014340.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/B0F2FD76-3877-E611-864B-02163E0136FC.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/06C15E74-3877-E611-B3CC-02163E0133E7.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/F66EB75A-3877-E611-91BB-02163E0138F9.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/FE84BC6E-3877-E611-ADA2-02163E0146B8.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/1EE10F6D-3877-E611-86FD-02163E0138C1.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/F03A178F-3877-E611-AB7E-02163E0138F6.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/F8A06465-3877-E611-AE96-02163E0142BF.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/66815774-3877-E611-9AE3-02163E013577.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/D6DCBF9A-3877-E611-98CA-02163E012B20.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/6C5D0567-3877-E611-8254-02163E0145F9.root',
        '/store/data/Run2016G/JetHT/RAW/v1/000/280/385/00000/28F19F77-3877-E611-893E-02163E013979.root',

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
process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")
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
                                   fileName=cms.string('analyze_jetht.root'))

# process.hcalDigis.InputLabel = cms.InputTag("source")
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

convertToUnscheduled(process)

with open('dump.txt', 'w') as fd:
    fd.write(process.dumpPython())
