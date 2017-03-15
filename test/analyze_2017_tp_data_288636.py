import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PLOT', eras.Run2_2016)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '90X_dataRun2_Prompt_v2', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(100))

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('/store/data/Commissioning2017/HcalNZS/RAW-RECO/LogError-PromptReco-v1/000/288/636/00000/5A5C920D-D301-E711-A2FF-02163E01A3F5.root'))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')

process.TFileService = cms.Service("TFileService",
                                   closeFileFast=cms.untracked.bool(True),
                                   fileName=cms.string('analyze_288636.root'))

process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)
process.HcalTPGCoderULUT.read_XML_LUTs = cms.bool(True)
process.HcalTPGCoderULUT.inputLUTs = cms.FileInPath('Debug/HcalDebug/data/testPhysics2016v5ctry02newV2.xml')

process.simHcalTriggerPrimitiveDigis.upgradeHF = cms.bool(False)
process.simHcalTriggerPrimitiveDigis.upgradeHE = cms.bool(True)
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag("hcalDigis", "hcalDigis")
process.simHcalTriggerPrimitiveDigis.inputUpgradeLabel = cms.VInputTag("hcalDigis", "hcalDigis")

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

# print process.dumpPython()
