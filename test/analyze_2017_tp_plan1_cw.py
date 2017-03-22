import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PLOT', eras.Run2_2017)

# Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
# print("Using GlobalTag {}".format(process.GlobalTag.globaltag.value()))
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '90X_upgrade2017_realistic_v20'

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('file:step3.root'),
                            secondaryFileNames=cms.untracked.vstring('file:step2.root'))

# process.out = cms.OutputModule( "PoolOutputModule",
#         fileName = cms.untracked.string("output.root"),
#         outputCommands = cms.untracked.vstring( 'keep *' )
# )
# process.end = cms.EndPath(process.out)

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag(
    cms.InputTag('simHcalUnsuppressedDigis'), cms.InputTag('simHcalUnsuppressedDigis'))
# process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('simHcalDigis'), cms.InputTag('simHcalDigis') )
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")

process.TFileService = cms.Service("TFileService",
                                   closeFileFast=cms.untracked.bool(True),
                                   fileName=cms.string('analyze.root'))

process.analyze = cms.EDAnalyzer("AnalyzeTP",
                                 triggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis", "", ""))
process.analyzeRaw = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("hcalDigis", "", ""))
process.chainplotter = cms.EDAnalyzer("HcalCompareLegacyChains",
                                      triggerPrimitives=cms.InputTag('simHcalTriggerPrimitiveDigis', '', ''),
                                      recHits=cms.VInputTag('hbhereco', 'hfreco'),
                                      dataFrames=cms.VInputTag(cms.InputTag("hcalDigis", "", ""), cms.InputTag("hcalDigis", "", "")),
                                      swapIphi=cms.bool(False)
                                      )

process.p = cms.Path(process.hcalDigis * process.analyze * process.analyzeRaw * process.chainplotter)

# print process.dumpPython()
process.source = cms.Source('PoolSource',
        fileNames =
        cms.untracked.vstring(
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_1.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_10.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_11.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_12.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_13.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_14.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_15.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_16.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_17.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_18.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_19.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_2.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_20.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_21.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_22.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_23.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_24.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_25.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_26.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_27.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_28.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_29.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_3.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_30.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_31.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_32.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_33.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_34.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_35.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_36.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_37.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_38.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_39.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_4.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_40.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_41.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_42.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_43.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_44.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_45.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_46.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_47.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_48.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_49.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_5.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_50.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_51.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_52.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_53.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_54.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_55.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_56.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_57.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_58.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_59.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_6.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_60.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_61.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_62.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_63.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_64.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_65.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_66.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_67.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_68.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_69.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_7.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_70.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_71.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_72.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_73.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_74.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_75.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_76.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_77.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_78.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_79.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_8.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_80.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_81.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_82.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_83.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_84.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_85.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_86.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_87.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_88.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_89.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_9.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step3_90.root'
),
        secondaryFileNames =
        cms.untracked.vstring(
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_1.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_10.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_11.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_12.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_13.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_14.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_15.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_16.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_17.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_18.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_19.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_2.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_20.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_21.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_22.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_23.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_24.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_25.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_26.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_27.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_28.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_29.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_3.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_30.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_31.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_32.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_33.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_34.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_35.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_36.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_37.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_38.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_39.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_4.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_40.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_41.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_42.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_43.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_44.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_45.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_46.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_47.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_48.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_49.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_5.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_50.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_51.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_52.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_53.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_54.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_55.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_56.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_57.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_58.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_59.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_6.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_60.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_61.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_62.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_63.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_64.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_65.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_66.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_67.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_68.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_69.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_7.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_70.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_71.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_72.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_73.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_74.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_75.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_76.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_77.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_78.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_79.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_8.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_80.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_81.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_82.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_83.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_84.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_85.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_86.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_87.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_88.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_89.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_9.root',
'/store/user/cawest/RelValTTbar/CMSSW_9_0_0_pre6-90X_upgrade2017_realistic_v20/step2_90.root'
))
