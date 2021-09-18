#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_10_2_15_patch1/src ] ; then
  echo release CMSSW_10_2_15_patch1 already exists
else
  scram p CMSSW CMSSW_10_2_15_patch1
fi
cd CMSSW_10_2_15_patch1/src
eval `scram runtime -sh`

pyfile="py8_B0sToPhiMuMu_EvtGen_18GS_13TeV_cfi.py"
##PREGUNTAR POR EL PATH AL PYFILE
curl -s --insecure https://raw.githubusercontent.com/CesarMH18/MCProduction/master/$pyfile --retry 2 --create-dirs -o Configuration/GenProduction/python/$pyfile

scram b
cd ../../

cmsDriver.py Configuration/GenProduction/python/$pyfile --python_filename step0-GS-BsToPhiMuMu_cfg.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --fileout file:step0-GS-BsToPhiMuMu.root --conditions 102X_upgrade2018_realistic_v11 --beamspot Realistic25ns13TeVEarly2018Collision --step GEN,SIM --geometry DB:Extended --era Run2_2018 --no_exec --mc -n $EVENTS || exit $? ;#--fileout file:step0-GS-BsToPhiMuMu.root --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 102X_upgrade2018_realistic_v11 --beamspot Realistic25ns13TeVEarly2018Collision --step GEN,SIM --nThreads 1 --geometry DB:Extended --era Run2_2018 --customise Configuration/DataProcessing/Utils.addMonitoring --python_filename step0-GS-BsToPhiMuMu_cfg.py --no_exec -n 1000;
sed -i "20 a from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper \nrandSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)\nrandSvc.populate()" step0-GS-BsToPhiMuMu_cfg.py


export SCRAM_ARCH=slc6_amd64_gcc700

source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_10_2_13/src ] ; then
  echo release CMSSW_10_2_13 already exists
else
  scram p CMSSW CMSSW_10_2_13
fi
cd CMSSW_10_2_13/src
eval `scram runtime -sh`

scram b
cd ../..


cmsDriver.py step1 --filein file:step0-GS-BsToPhiMuMu.root --fileout file:step1-DR-BsToPhiMuMu.root --python_filename step1-DR-BsToPhiMuMu_cfg.py --eventcontent FEVTDEBUGHLT --pileup "AVE_25_BX_25ns,{'N': 20}" --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-DIGI-RAW --pileup_input "dbs:/MinBias_TuneCP5_13TeV-pythia8/RunIIFall18GS-102X_upgrade2018_realistic_v9-v1/GEN-SIM" --conditions 102X_upgrade2018_realistic_v15 --step DIGI,L1,DIGI2RAW,HLT:@relval2018 --geometry DB:Extended --era Run2_2018 --no_exec --mc -n $EVENTS || exit $? ;#--mc --eventcontent PREMIXRAW --datatier GEN-SIM-RAW --conditions 102X_upgrade2018_realistic_v15 --step DIGI,DATAMIX,L1,DIGI2RAW,HLT:@relval2018 --procModifiers premix_stage2 --nThreads 1 --geometry DB:Extended --datamix PreMix --era Run2_2018 --python_filename step1-DR-BsToPhiMuMu_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n -1;
sed -i "20 a from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper\nrandSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)\nrandSvc.populate() " step1-DR-BsToPhiMuMu_cfg.py

cmsDriver.py step2 --filein file:step1-DR-BsToPhiMuMu.root --fileout file:step2-DR-BsToPhiMuMu.root --python_filename step2-DR-BsToPhiMuMu_cfg.py --eventcontent AODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM --conditions 102X_upgrade2018_realistic_v15 --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI --geometry DB:Extended --era Run2_2018 --runUnscheduled --no_exec --mc -n $EVENTS || exit $? ;#--mc --eventcontent AODSIM --runUnscheduled --datatier AODSIM --conditions 102X_upgrade2018_realistic_v15 --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI --procModifiers premix_stage2 --nThreads 1 --era Run2_2018 --python_filename step2-DR-ups2s2ups1spipi_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n -1 ;
sed -i "20 a from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper\nrandSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)\nrandSvc.populate()" step2-DR-BsToPhiMuMu_cfg.py


source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_10_2_14/src ] ; then
  echo release CMSSW_10_2_14 already exists
else
  scram p CMSSW CMSSW_10_2_14
fi
cd CMSSW_10_2_14/src
eval `scram runtime -sh`

scram b
cd ../..

cmsDriver.py step3 --filein file:step2-DR-BsToPhiMuMu.root --fileout file:step3-MiniAOD-BsToPhiMuMu.root --python_filename step3-MiniAOD-BsToPhiMuMu_cfg.py --eventcontent MINIAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM --conditions 102X_upgrade2018_realistic_v15 --step PAT --geometry DB:Extended --era Run2_2018,bParking --runUnscheduled --no_exec --mc -n $EVENTS || exit $? ; #--mc --eventcontent MINIAODSIM --runUnscheduled --datatier MINIAODSIM --conditions 102X_upgrade2018_realistic_v15 --step PAT --nThreads 1 --geometry DB:Extended --era Run2_2018 --python_filename step3-MiniAOD-ups2s2ups1spipi_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n -1;
sed -i "20 a from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper\nrandSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)\nrandSvc.populate() " step3-MiniAOD-BsToPhiMuMu_cfg.py


export SCRAM_ARCH=slc6_amd64_gcc700

source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_10_2_11/src ] ; then
  echo release CMSSW_10_2_11 already exists
else
  scram p CMSSW CMSSW_10_2_11
fi
cd CMSSW_10_2_11/src
eval `scram runtime -sh`

scram b
cd ../..

cmsDriver.py step4 --filein file:step3-NanoAOD-BsToPhiMuMu.root --fileout file:step4-NanoAOD-BsToPhiMuMu.root --python_filename step4-NanoAOD-BsToPhiMuMu_cfg.py --eventcontent NANOEDMAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v16 --customise_commands 'process.particleLevelSequence.remove(process.genParticles2HepMCHiggsVtx);process.particleLevelSequence.remove(process.rivetProducerHTXS);process.particleLevelTables.remove(process.HTXSCategoryTable)' --step NANO --era Run2_2018,run2_nanoAOD_102Xv1 --no_exec --mc -n $EVENTS || exit $? ;



