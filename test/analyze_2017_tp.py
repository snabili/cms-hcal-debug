import FWCore.ParameterSet.Config as cms

process = cms.Process("HFCALIB")

# Import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['run2_mc']

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('file:step2.root'))

# process.out = cms.OutputModule( "PoolOutputModule",
#         fileName = cms.untracked.string("output.root"),
#         outputCommands = cms.untracked.vstring( 'keep *' )
# )
# process.end = cms.EndPath(process.out)

process.load('L1Trigger.RegionalCaloTrigger.rctDigis_cfi')
process.rctDigis.hcalDigis = cms.VInputTag(cms.InputTag("simHcalTriggerPrimitiveDigis"))

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag(cms.InputTag('simHcalUnsuppressedDigis'), cms.InputTag('simHcalUnsuppressedDigis'))
# process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag( cms.InputTag('simHcalDigis'), cms.InputTag('simHcalDigis') )
process.simHcalTriggerPrimitiveDigis.FrontEndFormatError = cms.bool(False)

process.load("Configuration.Geometry.GeometryExtended2016Reco_cff")
from CondCore.CondDB.CondDB_cfi import CondDB

process.CondDBSetup = CondDB.clone()
delattr(process.CondDBSetup, 'connect')

# process.es_pool = cms.ESSource("PoolDBESSource",
#      process.CondDBSetup,
#      timetype = cms.string('runnumber'),
#      toGet = cms.VPSet(
#          cms.PSet(record = cms.string("HcalLutMetadataRcd"),
#              tag = cms.string("HcalLutMetadata_HFTP_1x1")
#              ),
#          cms.PSet(record = cms.string("HcalElectronicsMapRcd"),
#              tag = cms.string("HcalElectronicsMap_HFTP_1x1")
#              )
#          ),
#      connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
#      authenticationMethod = cms.untracked.uint32(0)
#      )
# process.es_prefer_es_pool = cms.ESPrefer( "PoolDBESSource", "es_pool" )

from SLHCUpgradeSimulations.Configuration.HCalCustoms import customise_Hcal2017Full
customise_Hcal2017Full(process)

process.TFileService = cms.Service("TFileService",
                                   closeFileFast=cms.untracked.bool(True),
                                   fileName=cms.string('analyze.root'))

process.analyze = cms.EDAnalyzer("AnalyzeTP",
                                 triggerPrimitives=cms.InputTag("simHcalTriggerPrimitiveDigis", "", ""))

process.p = cms.Path(process.analyze)

# print process.dumpPython()
