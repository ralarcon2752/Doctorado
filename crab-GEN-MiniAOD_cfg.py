from CRABClient.UserUtilities import config
import datetime
import time

config = config()

ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d-%H-%M')

channel = 'B0sToPhiMuMu'
year = '2018'
step = 'PrivateMC-'+year
nEvents = 1000
NJOBS = 1
myrun = 'step0-GS-B0sToPhiMuMu_cfg.py'
myname = step+'-'+channel

config.General.requestName = step+'-'+channel+'-'+st
config.General.transferOutputs = True
config.General.transferLogs = False
config.General.workArea = 'crab_'+step+'-'+channel

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = myrun
config.JobType.inputFiles = ['step1-DR-B0sToPhiMuMu_cfg.py',
                             'step2-DR-B0sToPhiMuMu_cfg.py',
                             'step3-MiniAOD-B0sToPhiMuMu_cfg.py',
                             'step4-NanoAOD-B0sToPhiMuMu_cfg.py'
                             #'Skimming_MC_cfg.py'
                             ]
config.JobType.disableAutomaticOutputCollection = True
config.JobType.eventsPerLumi = 10000
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 3300
config.JobType.scriptExe = 'MCcrabJobScript.sh'
config.JobType.scriptArgs = ['CHANNEL_DECAY='+channel,'YEAR='+year]
config.JobType.outputFiles = ['step0-GS-B0sToPhiMuMu.root ', 'step4-NanoAOD-B0sToPhiMuMu.root']
config.Data.outputPrimaryDataset = myname
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = nEvents
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ralarcon/MCTest'
config.Data.publication = False

config.Site.storageSite = 'T3_CH_CERNBOX'