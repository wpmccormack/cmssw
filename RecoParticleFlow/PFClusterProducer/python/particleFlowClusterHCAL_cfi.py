import FWCore.ParameterSet.Config as cms

_thresholdsHB = cms.vdouble(0.8, 0.8, 0.8, 0.8)
_thresholdsHE = cms.vdouble(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)
_thresholdsHBphase1 = cms.vdouble(0.1, 0.2, 0.3, 0.3)
_thresholdsHEphase1 = cms.vdouble(0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)

particleFlowClusterHCAL = cms.EDProducer( "PFClusterSonicProducer",
    pf_src = cms.InputTag( "particleFlowRecHitHBHE" ),
    Client = cms.PSet(
        timeout = cms.untracked.uint32(300),
        mode = cms.string("Async"),
        modelName = cms.string("spvcnn"),
        modelConfigPath = cms.FileInPath("RecoParticleFlow/PFClusterProducer/data/models/spvcnn/config.pbtxt"),
        modelVersion = cms.string("1"),
        verbose = cms.untracked.bool(False),
        allowedTries = cms.untracked.uint32(0),
        useSharedMemory = cms.untracked.bool(True),
        compression = cms.untracked.string(""),
    ),
    allCellsPositionCalc = cms.PSet(
        algoName = cms.string("Basic2DGenericPFlowPositionCalc"),
        minFractionInCalc = cms.double(1e-9),    
        posCalcNCrystals = cms.int32(-1),
        logWeightDenominatorByDetector = cms.VPSet(
            cms.PSet( detector = cms.string("HCAL_BARREL1"),
                      depths = cms.vint32(1, 2, 3, 4),
                      logWeightDenominator = _thresholdsHB,
                  ),
            cms.PSet( detector = cms.string("HCAL_ENDCAP"),
                      depths = cms.vint32(1, 2, 3, 4, 5, 6, 7),
                      logWeightDenominator = _thresholdsHE,
                  )
        ),
        minAllowedNormalization = cms.double(1e-9)
    )
)

# offline 2018 -- uncollapsed
from Configuration.Eras.Modifier_run2_HE_2018_cff import run2_HE_2018
from Configuration.ProcessModifiers.run2_HECollapse_2018_cff import run2_HECollapse_2018
(run2_HE_2018 & ~run2_HECollapse_2018).toModify(particleFlowClusterHCAL,
    allCellsPositionCalc = dict(logWeightDenominatorByDetector = {1 : dict(logWeightDenominator = _thresholdsHEphase1) } ),
)

# offline 2021
from Configuration.Eras.Modifier_run3_HB_cff import run3_HB
run3_HB.toModify(particleFlowClusterHCAL,
    allCellsPositionCalc = dict(logWeightDenominatorByDetector = {0 : dict(logWeightDenominator = _thresholdsHBphase1) } ),
)

# HCALonly WF
particleFlowClusterHCALOnly = particleFlowClusterHCAL.clone(
    pf_src = "particleFlowRecHitHBHEOnly"
)
