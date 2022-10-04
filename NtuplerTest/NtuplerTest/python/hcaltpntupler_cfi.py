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


simHcalTriggerPrimitiveNtuple = cms.EDProducer("HcalTrigPrimDigiNtupler",
    clus_src = cms.InputTag("printTheDetIDs","HBHERecHitTruth")
)

simHcalTriggerPrimitiveNtupleReco = cms.EDProducer("HcalTrigPrimDigiNtupler",
    clus_src = cms.InputTag("printTheDetIDsReco","HBHERecHitTruth")
)

simHcalTriggerPrimitiveNtupleWithRadialCutoff = cms.EDProducer("HcalTrigPrimDigiNtupler",
    clus_src = cms.InputTag("printTheDetIDsWithRadialCutoff","HBHERecHitTruth")
)

simHcalTriggerPrimitiveNtupleRecoWithRadialCutoff = cms.EDProducer("HcalTrigPrimDigiNtupler",
    clus_src = cms.InputTag("printTheDetIDsRecoWithRadialCutoff","HBHERecHitTruth")
)

#,
"""                                             
  GenInfo = cms.untracked.PSet(
    isActive            = ( cms.untracked.bool(True) ),
    edmGenEventInfoName     = cms.untracked.string('generator'),
    edmGenParticlesName     = cms.untracked.string('prunedGenParticles'),
    edmGenPackParticlesName = cms.untracked.string('packedGenParticles'),
    fillAllGen          = cms.untracked.bool(True),
    fillLHEWeights      = cms.untracked.bool(False)
  ),

  GenJet  = cms.untracked.PSet(
    isActive            = ( cms.untracked.bool(True)),
    isActiveFatJet      = ( cms.untracked.bool(True)),
    #isActive            = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    #isActiveFatJet      = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmGenParticlesName = cms.untracked.string('prunedGenParticles'),
    genJetName          = cms.untracked.string('slimmedGenJets'),
    genFatJetName       = cms.untracked.string('slimmedGenJetsAK8'),
    fillAllGen          = cms.untracked.bool(True)
  ),

  Info = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    edmPFCandName        = cms.untracked.string('packedPFCandidates'),
    edmPileupInfoName    = cms.untracked.string('slimmedAddPileupInfo'),
    #edmPileupInfoName    = cms.untracked.string('addPileupInfo'),
    edmBeamspotName      = cms.untracked.string('offlineBeamSpot'),
    edmECALPFClusterName = cms.untracked.string('particleFlowClusterECAL'),
    edmMETName           = cms.untracked.string('slimmedMETs'),
    edmPFMETName         = cms.untracked.InputTag('slimmedMETs'),
    #edmPFMETName         = cms.untracked.InputTag('slimmedMETsV2','','MakingBacon'),
    edmMVAMETName        = cms.untracked.string(''),
    #edmPuppETName        = cms.untracked.string('slimmedMETsPuppi'),
    edmPuppETName        = cms.untracked.InputTag('slimmedMETsPuppi','',''), 
    #edmAlpacaMETName     = cms.untracked.string(alpacaMet),
    #edmPupAlpacaMETName  = cms.untracked.string(alpacaPuppiMet),
    edmRhoForIsoName     = cms.untracked.string('fixedGridRhoFastjetAll'),
    edmRhoForJetEnergy   = cms.untracked.string('fixedGridRhoFastjetAll'),
    doFillMETFilters     = cms.untracked.bool(True),
    doFillMET            = cms.untracked.bool(True)
  ),

    peakFilter = cms.bool(True),
    weights = cms.vdouble(1.0, 0.0), ##hardware algo        
    #weights = cms.vdouble(1.0, 1.0), ##hardware algo        
    latency = cms.int32(1),
    FG_threshold = cms.uint32(12), ## threshold for setting fine grain bit
    FG_HF_threshold = cms.uint32(17), ## threshold for setting fine grain bit
    ZS_threshold = cms.uint32(1),  ## threshold for setting fine grain bit
    numberOfSamples = cms.int32(4),
    numberOfPresamples = cms.int32(2),
    numberOfSamplesHF = cms.int32(2),
    numberOfPresamplesHF = cms.int32(1),
    MinSignalThreshold = cms.uint32(0), # For HF PMT veto
    PMTNoiseThreshold = cms.uint32(0),  # For HF PMT veto
    LSConfig=LSParameter,
    upgradeHF = cms.bool(True),
    edmPileupInfoName = cms.untracked.string("addPileupInfo"),                
    upgradeHB = cms.bool(False),
    upgradeHE = cms.bool(True),
    inputLabel = cms.VInputTag(cms.InputTag('simHcalUnsuppressedDigis'),
                               cms.InputTag('simHcalUnsuppressedDigis')),
    inputUpgradeLabel = cms.VInputTag(
        cms.InputTag('simHcalUnsuppressedDigis:HBHEQIE11DigiCollection'),
        cms.InputTag('simHcalUnsuppressedDigis:HFQIE10DigiCollection')),
    #edmRecHitName = cms.untracked.InputTag("hbheprereco"),
    edmRecHitName = cms.untracked.InputTag("HBHERecHitTruth"),
    edmRecHitNameEE = cms.untracked.InputTag("hbheprereco"),
    edmRecHitNameEB = cms.untracked.InputTag("hbheprereco"),
    #edmRecHitName = cms.InputTag("hbheprerecoM2"),
    #edmRecHitName = cms.InputTag("hbheprerecoM3"),
    #edmRecHitName = cms.InputTag("hbheprerecoMahi"),
    edmChanInfoName = cms.untracked.InputTag("hbheprereco"),                                           
    InputTagFEDRaw = cms.InputTag("rawDataCollector"),
    RunZS = cms.bool(False),
    FrontEndFormatError = cms.bool(False), # Front End Format Error, for real data only
    PeakFinderAlgorithm = cms.int32(2),
    tpScales = tpScales,
"""
#)
