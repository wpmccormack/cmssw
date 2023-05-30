import FWCore.ParameterSet.Config as cms
from RecoParticleFlow.PFClusterProducer.particleFlowCaloResolution_cfi import _timeResolutionHCALMaxSample

_thresholdsHB = cms.vdouble(0.8, 0.8, 0.8, 0.8)
_thresholdsHE = cms.vdouble(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)
_thresholdsHBphase1 = cms.vdouble(0.1, 0.2, 0.3, 0.3)
_thresholdsHEphase1 = cms.vdouble(0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)

_thresholdsHB = cms.vdouble(0.8, 0.8, 0.8, 0.8)
_thresholdsHE = cms.vdouble(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)
_thresholdsHBphase1 = cms.vdouble(0.1, 0.2, 0.3, 0.3)
_thresholdsHEphase1 = cms.vdouble(0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
_seedingThresholdsHB = cms.vdouble(1.0, 1.0, 1.0, 1.0)
_seedingThresholdsHE = cms.vdouble(1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1)
_seedingThresholdsHBphase1 = cms.vdouble(0.125, 0.25, 0.35, 0.35)
_seedingThresholdsHEphase1 = cms.vdouble(0.1375, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275)


particleFlowClusterHCAL = cms.EDProducer( "PFTruthClusterProducer2",
    recHitsSource = cms.InputTag( "particleFlowRecHitHBHE" ),
    recHitCleaners = cms.VPSet(
    ),
    seedCleaners = cms.VPSet(
    ),
    positionReCalc = cms.PSet(  ),
    energyCorrector = cms.PSet(  ),
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
      ),
    seedFinder = cms.PSet(
        algoName = cms.string("LocalMaximumSeedFinder"),
        thresholdsByDetector = cms.VPSet(
              cms.PSet( detector = cms.string("HCAL_BARREL1"),
                        depths = cms.vint32(1, 2, 3, 4),
                        seedingThreshold = _seedingThresholdsHB,
                        seedingThresholdPt = cms.vdouble(0.0, 0.0, 0.0, 0.0)
                        ),
              cms.PSet( detector = cms.string("HCAL_ENDCAP"),
                        depths = cms.vint32(1, 2, 3, 4, 5, 6, 7),
                        seedingThreshold = _seedingThresholdsHE,
                        seedingThresholdPt = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                        )
              ),
        nNeighbours = cms.int32(4)
    ),
    initialClusteringStep = cms.PSet(
        algoName = cms.string("Basic2DGenericTopoClusterizer"),    
        thresholdsByDetector = cms.VPSet(
        cms.PSet( detector = cms.string("HCAL_BARREL1"),
                  depths = cms.vint32(1, 2, 3, 4),
                  gatheringThreshold = _thresholdsHB,
                  gatheringThresholdPt = cms.vdouble(0.0, 0.0, 0.0, 0.0)
                  ),
        cms.PSet( detector = cms.string("HCAL_ENDCAP"),
                  depths = cms.vint32(1, 2, 3, 4, 5, 6, 7),
                  gatheringThreshold = _thresholdsHE,
                  gatheringThresholdPt = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                  )
        ),
        useCornerCells = cms.bool(True)
    ),
    pfClusterBuilder = cms.PSet(
           algoName = cms.string("Basic2DGenericPFlowClusterizer"),
           #pf clustering parameters
           minFractionToKeep = cms.double(1e-7),
           positionCalc = cms.PSet(
                 algoName = cms.string("Basic2DGenericPFlowPositionCalc"),
                 minFractionInCalc = cms.double(1e-9),    
                 posCalcNCrystals = cms.int32(5),
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
           ),
           allCellsPositionCalc =cms.PSet(
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
           ),
           

           timeSigmaEB = cms.double(10.),
           timeSigmaEE = cms.double(10.),
           maxNSigmaTime = cms.double(10.),
           minChi2Prob = cms.double(0.),
           clusterTimeResFromSeed = cms.bool(False),
           timeResolutionCalcBarrel = _timeResolutionHCALMaxSample,
           timeResolutionCalcEndcap = _timeResolutionHCALMaxSample,
           showerSigma = cms.double(10.0),
           stoppingTolerance = cms.double(1e-8),
           maxIterations = cms.uint32(50),
           excludeOtherSeeds = cms.bool(True),
           minFracTot = cms.double(1e-20), ## numerical stabilization
           recHitEnergyNorms = cms.VPSet(
            cms.PSet( detector = cms.string("HCAL_BARREL1"),
                      depths = cms.vint32(1, 2, 3, 4),
                      recHitEnergyNorm = _thresholdsHB,
                      ),
            cms.PSet( detector = cms.string("HCAL_ENDCAP"),
                      depths = cms.vint32(1, 2, 3, 4, 5, 6, 7),
                      recHitEnergyNorm = _thresholdsHE,
                      )
            )
    ),
    pfClusterBuilderDepth =cms.PSet(
       algoName = cms.string("PFMultiDepthClusterizer"),
       nSigmaEta = cms.double(2.),
       nSigmaPhi = cms.double(2.),
       #pf clustering parameters
       minFractionToKeep = cms.double(1e-7),
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
    recHitsSource = "particleFlowRecHitHBHEOnly"
)
