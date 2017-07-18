import FWCore.ParameterSet.Config as cms


def add_fileservice(process):
    process.TFileService = cms.Service("TFileService",
                                       closeFileFast=cms.untracked.bool(True),
                                       fileName=cms.string('analyze.root'))


def add_path(process):
    if not hasattr(process, 'tpCheck'):
        process.tpCheck = cms.Path()
        process.schedule.append(process.tpCheck)


def analyze_raw_tp(process):
    add_fileservice(process)
    add_path(process)
    process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
    process.analyzeRaw = cms.EDAnalyzer("AnalyzeTP",
                                        triggerPrimitives=cms.InputTag("hcalDigis", "", ""))
    process.tpCheck *= process.analyzeRaw
    return process


def analyze_reemul_tp(process):
    add_fileservice(process)
    add_path(process)
    process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
    process.analyzeReemul = cms.EDAnalyzer("AnalyzeTP",
                                           triggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis", "", process.name_()))
    process.tpCheck *= process.analyzeReemul
    return process


def compare_raw_reemul_tp(process):
    add_fileservice(process)
    add_path(process)
    process.compare = cms.EDAnalyzer("CompareTP",
                                     swapIphi=cms.bool(False),
                                     triggerPrimitives=cms.InputTag("hcalDigis", "", process.name_()),
                                     emulTriggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis", "", process.name_()))
    process.tpCheck *= process.compare
    return process
