import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PLOT', eras.Run2_2016)

# Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",
                            skipBadFiles=cms.untracked.bool(True),
                            fileNames=cms.untracked.vstring(),
                            secondaryFileNames=cms.untracked.vstring(),
                            firstRun=cms.untracked.uint32(260627)
                            )

process.source.fileNames.extend([
    '/store/data/Run2016H/JetHT/RECO/PromptReco-v2/000/282/735/00000/F44821A2-6890-E611-9FC9-02163E0139C8.root'
])
process.source.secondaryFileNames.extend([
    '/store/data/Run2016H/JetHT/RAW/v1/000/282/735/00000/F47698CE-5C8E-E611-B643-02163E011CA5.root'
])

process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)

process.load("Configuration.Geometry.GeometryExtended2016Reco_cff")

process.TFileService = cms.Service("TFileService",
                                   closeFileFast=cms.untracked.bool(True),
                                   fileName=cms.string('analyze_data.root'))

process.analyzeRaw = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("hcalDigis", "", "PLOT"))

process.chainplotter = cms.EDAnalyzer("HcalCompareLegacyChains",
                                      triggerPrimitives=cms.InputTag('hcalDigis', '', 'PLOT'),
                                      recHits=cms.VInputTag('hbhereco', 'hfreco'),
                                      dataFrames=cms.VInputTag(
                                          cms.InputTag("hcalDigis", "", "PLOT"),
                                          cms.InputTag("hcalDigis", "", "PLOT")
                                      ),
                                      swapIphi=cms.bool(False)
                                      )

process.p = cms.Path(
    process.hcalDigis
    * process.analyzeRaw
    * process.chainplotter
)

# print process.dumpPython()
