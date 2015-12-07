# Just look at stored TP in a file, don't run the algorithm at all.  Used
# to debug the PR adding 1x1 TP.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.parseArguments()

process = cms.Process("ANALYZE")

## Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(options.inputFiles)
)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Digi_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_mc', '')

process.TFileService = cms.Service("TFileService",
        closeFileFast = cms.untracked.bool(True),
        fileName = cms.string(options.outputFile))

process.analyze = cms.EDAnalyzer("AnalyzeTP",
        triggerPrimitives = cms.InputTag("simHcalTriggerPrimitiveDigis", "" , ""),
        useDetIdForUncompression = cms.untracked.bool(False),
        threshold = cms.untracked.double(0.00001))

process.p = cms.Path(process.analyze)
