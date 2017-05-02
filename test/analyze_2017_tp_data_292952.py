import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PLOT', eras.Run2_2016)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '90X_dataRun2_HLT_Candidate_HCAL_payloads_tests_2017_03_24', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("HcalTBSource",
                            fileNames=cms.untracked.vstring('/store/group/dpg_hcal/comm_hcal/USC/run292952/USC_292952.root'))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')

process.hcalDigis.InputLabel = cms.InputTag("source")

process.TFileService = cms.Service("TFileService",
                                   closeFileFast=cms.untracked.bool(True),
                                   fileName=cms.string('analyze_292952.root'))

process.emulTP = process.simHcalTriggerPrimitiveDigis.clone()
process.emulTP.upgradeHF = cms.bool(True)
process.emulTP.upgradeHE = cms.bool(True)
process.emulTP.inputLabel = cms.VInputTag("hcalDigis", "hcalDigis")
process.emulTP.inputUpgradeLabel = cms.VInputTag("hcalDigis", "hcalDigis")

# process.hcalDigis.InputLabel = cms.InputTag("source")
process.analyzeRAW = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("hcalDigis", "", ""))
process.analyzeSIM = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("emulTP", "", ""))
process.compare = cms.EDAnalyzer("CompareTP",
                                 triggerPrimitives=cms.InputTag("hcalDigis"),
                                 emulTriggerPrimitives=cms.InputTag("emulTP"),
                                 swapIphi=cms.bool(False))

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.hcalDigis *
    # process.dump *
    process.emulTP *
    process.analyzeRAW *
    process.analyzeSIM *
    process.compare)

# print process.dumpPython()
