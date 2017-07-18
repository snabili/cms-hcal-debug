import FWCore.ParameterSet.Config as cms


def add_fileservice(process):
    process.TFileService = cms.Service("TFileService",
                                       closeFileFast=cms.untracked.bool(True),
                                       fileName=cms.string('analyze.root'))


def add_path(process):
    if not hasattr(process, 'tpCheck'):
        process.tpCheck = cms.Path()
        process.schedule.append(process.tpCheck)


def analyze_tp(process, name, tag):
    add_fileservice(process)
    add_path(process)
    process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
    setattr(process, name, cms.EDAnalyzer("AnalyzeTP",
                                          triggerPrimitives=cms.InputTag(tag, "", "")))
    process.tpCheck *= getattr(process, name)
    return process


def analyze_raw_tp(process):
    return analyze_tp(process, 'analyzeRaw', 'hcalDigis')


def analyze_reemul_tp(process):
    return analyze_tp(process, 'analyzeReemul', 'simHcalTriggerPrimitiveDigis')


def compare_tp(process, name, tag1, tag2):
    add_fileservice(process)
    add_path(process)
    setattr(process, name, cms.EDAnalyzer("CompareTP",
                                          swapIphi=cms.bool(False),
                                          triggerPrimitives=cms.InputTag(tag1, "", process.name_()),
                                          emulTriggerPrimitives=cms.InputTag(tag2, "", process.name_())))
    process.tpCheck *= getattr(process, name)
    return process


def compare_raw_reemul_tp(process):
    return compare_tp(process, 'compare', 'hcalDigis', 'simHcalTriggerPrimitiveDigis')
