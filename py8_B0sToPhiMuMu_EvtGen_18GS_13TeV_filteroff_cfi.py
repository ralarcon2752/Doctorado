import FWCore.ParameterSet.Config as cms
from GeneratorInterface.EvtGenInterface.EvtGenSetting_cff import *
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            convertPythiaCodes = cms.untracked.bool(False),
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            list_forced_decays = cms.vstring('MyB0s', 
                'Myanti-B0s'),
            operates_on_particles = cms.vint32(),
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2014.pdl'),
            user_decay_embedded = cms.vstring(
"""
Alias      MyB0s        B_s0
Alias      Myanti-B0s   anti-B_s0
ChargeConj MyB0s        Myanti-B0s
Alias      MyPhi        phi
ChargeConj MyPhi        MyPhi
#
Decay MyB0s
  1.000        MyPhi     mu+     mu-               PHOTOS PHSP;
Enddecay
Decay Myanti-B0s
  1.000        MyPhi     mu+     mu-            PHOTOS PHSP;
Enddecay
#
Decay MyPhi
  1.000        K+        K-                    PHOTOS PHSP;
Enddecay
End
"""
            )
        ),
        parameterSets = cms.vstring('EvtGen130')
    ),
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring('pythia8CommonSettings', 
            'pythia8CP5Settings', 
            'processParameters'),
        processParameters = cms.vstring('SoftQCD:nonDiffractive = on', 
            'PTFilter:filter = on', 
            'PTFilter:quarkToFilter = 5', 
            'PTFilter:scaleToFilter = 1.0'),
        pythia8CP5Settings = cms.vstring('Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:ecmPow=0.03344', 
            'PDF:pSet=20', 
            'MultipartonInteractions:bProfile=2', 
            'MultipartonInteractions:pT0Ref=1.41', 
            'MultipartonInteractions:coreRadius=0.7634', 
            'MultipartonInteractions:coreFraction=0.63', 
            'ColourReconnection:range=5.176', 
            'SigmaTotal:zeroAXB=off', 
            'SpaceShower:alphaSorder=2', 
            'SpaceShower:alphaSvalue=0.118', 
            'SigmaProcess:alphaSvalue=0.118', 
            'SigmaProcess:alphaSorder=2', 
            'MultipartonInteractions:alphaSvalue=0.118', 
            'MultipartonInteractions:alphaSorder=2', 
            'TimeShower:alphaSorder=2', 
            'TimeShower:alphaSvalue=0.118'),
        pythia8CommonSettings = cms.vstring('Tune:preferLHAPDF = 2', 
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on')
    ),
    comEnergy = cms.double(13000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(0)
)

###### Filters ##########

bfilter = cms.EDFilter("PythiaFilter",
    MaxEta = cms.untracked.double(9999999),
    MinEta = cms.untracked.double(-9999999),
    ParticleID = cms.untracked.int32(531)
)

decayfilter = cms.EDFilter("PythiaDauVFilter",
    DaughterIDs = cms.untracked.vint32(-13, 13, 333),
    MaxEta = cms.untracked.vdouble(9999999,9999999,9999999),
    MinEta = cms.untracked.vdouble(-9999999,-9999999,-9999999),
    MinPt = cms.untracked.vdouble(0.0,0.0,0.0),
    NumberDaughters = cms.untracked.int32(3),
    ParticleID = cms.untracked.int32(531),
    verbose = cms.untracked.int32(1)
)

phifilter = cms.EDFilter("PythiaDauVFilter",
    DaughterIDs = cms.untracked.vint32(321, -321),
    MaxEta = cms.untracked.vdouble(9999999,9999999),
    MinEta = cms.untracked.vdouble(-9999999,-9999999),
    MinPt = cms.untracked.vdouble(0.0,0.0),
    MotherID = cms.untracked.int32(531),
    NumberDaughters = cms.untracked.int32(2),
    ParticleID = cms.untracked.int32(333),
    verbose = cms.untracked.int32(1)
)

mu3filter = cms.EDFilter("MCMultiParticleFilter",
            src = cms.untracked.InputTag("generator", "unsmeared"),   
            Status = cms.vint32(1),
            ParticleID = cms.vint32(13),
            PtMin = cms.vdouble(0.0),
            NumRequired = cms.int32(3),
            EtaMax = cms.vdouble(9999999),
            AcceptMore = cms.bool(True)
            )

#mufilter = cms.EDFilter("PythiaFilter",  # bachelor muon with kinematic cuts.
#    MaxEta = cms.untracked.double(2.5),
#    MinEta = cms.untracked.double(-2.5),
#    MinPt = cms.untracked.double(5.),
#    ParticleID = cms.untracked.int32(13),
#)

#probefilter = cms.EDFilter(  
 #   "PythiaProbeFilter",
  #  verbose         = cms.untracked.int32(1),
   # NumberOfSisters = cms.untracked.int32(2),
    #ParticleID      = cms.untracked.int32(13),    
#    MomID           = cms.untracked.int32(531),
 #   SisterIDs     = cms.untracked.vint32(333,-13),
  #  MinPt           = cms.untracked.double(5.),
   # MinEta          = cms.untracked.double(-2.5),
    #MaxEta          = cms.untracked.double(2.5)
    #)



ProductionFilterSequence = cms.Sequence(generator*bfilter*decayfilter*phifilter*mu3filter)
