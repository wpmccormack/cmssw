import FWCore.ParameterSet.Config as cms

from CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi import tpScales
from Configuration.Eras.Modifier_run2_HE_2017_cff import run2_HE_2017
from Configuration.Eras.Modifier_run2_HF_2017_cff import run2_HF_2017
from Configuration.Eras.Modifier_run3_HB_cff import run3_HB

LSParameter =cms.untracked.PSet(
HcalFeatureHFEMBit= cms.bool(False),
Min_Long_Energy= cms.double(10),#makes a cut based on energy deposited in short vrs long
    Min_Short_Energy= cms.double(10),
    Long_vrs_Short_Slope= cms.double(100.2),
    Long_Short_Offset= cms.double(10.1))


printTheDetIDs = cms.EDProducer("ProducerTest",
    #pf_src = cms.InputTag( "hbheprereco" ),
    #pf_src = cms.InputTag( "hltParticleFlowRecHitHBHE" ),
    pf_src = cms.InputTag( "hltHbhereco" ),
    clus_src = cms.InputTag("hltParticleFlowClusterHCAL")
)

printTheDetIDsReco = cms.EDProducer("ProducerTest",
    pf_src = cms.InputTag( "hbhereco" ),
    clus_src = cms.InputTag("particleFlowClusterHCAL")
)


printTheDetIDsWithRadialCutoff = cms.EDProducer("ProducerTestWithRadialCutoff",
    #pf_src = cms.InputTag( "hbheprereco" ),
    #pf_src = cms.InputTag( "hltParticleFlowRecHitHBHE" ),
    pf_src = cms.InputTag( "hltHbhereco" ),
    clus_src = cms.InputTag("hltParticleFlowClusterHCAL")
)

printTheDetIDsRecoWithRadialCutoff = cms.EDProducer("ProducerTestWithRadialCutoff",
    pf_src = cms.InputTag( "hbhereco" ),
    clus_src = cms.InputTag("particleFlowClusterHCAL")
)

printTheDetIDsRecoWithRadialCutoffECAL = cms.EDProducer("ProducerTestWithRadialCutoff_ECAL",
    pf_src = cms.InputTag( "hbhereco" ),
    clus_src = cms.InputTag("particleFlowClusterECAL")
)

printTheDetIDsWithRadialCutoffEXPANDEDINFO = cms.EDProducer("ProducerTestWithRadialCutoffEXPANDEDINFO",
    #pf_src = cms.InputTag( "hbheprereco" ),
    #pf_src = cms.InputTag( "hltParticleFlowRecHitHBHE" ),
    pf_src = cms.InputTag( "hltHbhereco" ),
    clus_src = cms.InputTag("hltParticleFlowClusterHCAL")
)

printTheDetIDsRecoWithRadialCutoffEXPANDEDINFO = cms.EDProducer("ProducerTestWithRadialCutoffEXPANDEDINFO",
    pf_src = cms.InputTag( "hbhereco" ),
    clus_src = cms.InputTag("particleFlowClusterHCAL")
)
