# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: analyze --data --conditions 101X_dataRun2_Prompt_v11 -s RAW2DIGI,L1Reco,RECO --geometry DB:Extended --era Run2_2018 --eventcontent FEVTDEBUG --customise Debug/HcalDebug/customize.compare_raw_reemul_tp --customise Debug/HcalDebug/customize.use_data_reemul_tp --customise_commands cms.untracked.vstring('keep *,keep FEDRawData_*_*_*') --filein=filelist:filelist_324077_ZeroBias.txt --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from EventFilter.L1TRawToDigi.caloStage2Digis_cfi import caloStage2Digis

process = cms.Process('RECO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:/data/users/snabili/CMSSW_10_3_2/src/0BF9DD6F-5FEA-F44B-AF2B-BF13A22CE1BE.root'),
    #fileNames = cms.untracked.vstring('file:/data/users/snabili/CMSSW_10_3_2/src/FOLDER/EF6F8572-12C6-A145-B6B4-4E850B119962.root'),
    fileNames = cms.untracked.vstring('file:/data/users/snabili/CMSSW_10_3_2/src/CD216871-A33B-3D49-B03A-81307557099B.root'),#JetHt samples
    #secondaryFileNames = cms.untracked.vstring(),
    #eventsToProcess = cms.untracked.VEventRange('324077:915920650-324077:max'),
    #eventsToSkip = cms.untracked.VEventRange('324077:916005548, 324077:916219900'), 
    #eventsToProcess = cms.untracked.VEventRange('324077:916544971-324077:max'),
    #eventsToProcess = cms.untracked.VEventRange('324077:917811743-324077:max'),
    #eventsToProcess = cms.untracked.VEventRange('324077:918423993-324077:max'),
    #eventsToProcess = cms.untracked.VEventRange('324077:920279702-324077:max'),
    #eventsToProcess = cms.untracked.VEventRange('324077:564574218-324077:max'),
    eventsToProcess = cms.untracked.VEventRange('325170:8173616-325170:max'),# JetHt samples: 8037535, 8138097, 8173615
    noEventSort = cms.untracked.bool(False)
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('analyze nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('test_file.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v11', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from Debug.HcalDebug.customize
from Debug.HcalDebug.customize import compare_raw_reemul_tp,use_data_reemul_tp,compare_raw_reco_sev9 

#call to customisation function compare_raw_reemul_tp imported from Debug.HcalDebug.customize
process = compare_raw_reemul_tp(process)

#call to customisation function use_data_reemul_tp imported from Debug.HcalDebug.customize
process = use_data_reemul_tp(process)

#call to customisation function compare_reemul_reco_sev9 imported from Debug.HcalDebug.customize
process = compare_raw_reco_sev9(process)

# End of customisation functions

# Customisation from command line
process.TFileService.fileName=cms.string('test.root')

cms.untracked.vstring('keep *,keep FEDRawData_*_*_*')
#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
