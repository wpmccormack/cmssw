import FWCore.ParameterSet.Config as cms

hltJetMatcherAK4 = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
    src         = cms.InputTag("hltAK4PFJets"),      # RECO jets (any View<Jet> is ok)
    matched     = cms.InputTag("ak4GenJets"),        # GEN jets  (must be GenJetCollection)
    mcPdgId     = cms.vint32(),                      # n/a
    mcStatus    = cms.vint32(),                      # n/a
    checkCharge = cms.bool(False),                   # n/a
    maxDeltaR   = cms.double(0.4),                   # Minimum deltaR for the match
    #maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first
)

hltJetMatcherAK8 = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
    src         = cms.InputTag("hltAK8PFJets"),      # RECO jets (any View<Jet> is ok)
    matched     = cms.InputTag("ak8GenJets"),        # GEN jets  (must be GenJetCollection)
    mcPdgId     = cms.vint32(),                      # n/a
    mcStatus    = cms.vint32(),                      # n/a
    checkCharge = cms.bool(False),                   # n/a
    maxDeltaR   = cms.double(0.4),                   # Minimum deltaR for the match
    #maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first
)
