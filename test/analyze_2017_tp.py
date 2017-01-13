import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PLOT', eras.Run2_2017)

# Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1000))

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

process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')

# customise_Hcal2017Full(process)
# process.simHcalTriggerPrimitiveDigis.upgradeHE = cms.bool(True)
# process.simHcalTriggerPrimitiveDigis.upgradeHF = cms.bool(True)

process.TFileService = cms.Service("TFileService",
                                   closeFileFast=cms.untracked.bool(True),
                                   fileName=cms.string('analyze.root'))

process.analyze = cms.EDAnalyzer("AnalyzeTP",
                                 triggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis", "", ""))
process.chainplotter = cms.EDAnalyzer("HcalCompareLegacyChains",
                                      triggerPrimitives=cms.InputTag(
                                          'simHcalTriggerPrimitiveDigis', '', ''),
                                      recHits=cms.VInputTag('hbheprereco', 'hfreco'),
                                      dataFrames=cms.VInputTag(
                                          cms.InputTag("hcalDigis", "", ""),
                                          cms.InputTag("hcalDigis", "", "")
                                      ),
                                      swapIphi=cms.bool(False)
                                      )

process.p = cms.Path(process.analyze * process.chainplotter)

# print process.dumpPython()
