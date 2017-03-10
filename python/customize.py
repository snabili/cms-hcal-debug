import FWCore.ParameterSet.Config as cms


def add_fileservice(process):
    process.TFileService = cms.Service("TFileService",
                                       closeFileFast=cms.untracked.bool(True),
                                       fileName=cms.string('analyze.root'))


def analyze_raw_tp(process):
    add_fileservice(process)
    process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
    process.analyzeRaw = cms.EDAnalyzer("AnalyzeTP",
                                        triggerPrimitives=cms.InputTag("hcalDigis", "", ""))
    process.analyzeRawP = cms.Path(process.analyzeRaw)
    process.schedule.append(process.analyzeRawP)
    return process


def analyze_reemul_tp(process):
    add_fileservice(process)
    process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
    process.analyzeReemul = cms.EDAnalyzer("AnalyzeTP",
                                           triggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis", "", process.name_()))
    process.analyzeReemulP = cms.Path(process.analyzeReemul)
    process.schedule.append(process.analyzeReemulP)
    return process


def compare_raw_reemul_tp(process):
    add_fileservice(process)
    process.compare = cms.EDAnalyzer("CompareTP",
                                     swapIphi=cms.bool(False),
                                     triggerPrimitives=cms.InputTag("hcalDigis", "", process.name_()),
                                     emulTriggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis", "", process.name_()))
    process.compareP = cms.Path(process.compare)
    process.schedule.append(process.compareP)
    return process
