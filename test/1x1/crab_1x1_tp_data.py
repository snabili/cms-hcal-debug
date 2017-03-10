from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Debug/HcalDebug/test/analyze_1x1_tp_data.py'
config.section_("Data")
config.Data.inputDataset = '/HcalNZS/Run2016B-v2/RAW'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 10000
config.Data.runRange = '273158'
config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
