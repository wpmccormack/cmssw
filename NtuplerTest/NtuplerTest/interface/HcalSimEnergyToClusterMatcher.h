#ifndef HcalSimEnergyToClusterMatcher_h
#define HcalSimEnergyToClusterMatcher_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimCalorimetry/HcalTrigPrimAlgos/interface/HcalTriggerPrimitiveAlgo.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterBuilderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/RecHitTopologicalCleanerBase.h"
#include "RecoParticleFlow/PFProducer/interface/PFAlgo.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHitTruth.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitDefs.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHitTruth.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimCalorimetry/CaloSimAlgos/interface/CaloHitResponse.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"


class HcalSimEnergyToClusterMatcher : public edm::stream::EDProducer<> {
public:

  explicit HcalSimEnergyToClusterMatcher(const edm::ParameterSet& ps);
  ~HcalSimEnergyToClusterMatcher() override { endJob(); }//delete fTrackArr; delete fVertexArr; }// delete fFillerRH; delete fEvtArr; delete fFillerEventInfo; delete fFillerGenInfo; delete fGenParArr;}
  //~HcalSimEnergyToClusterMatcher() override { endJob(); delete fPFParArr; delete fFillerPF; delete fRHParArr; delete fFillerRH;}
  /**Produces the EDM products,*/
  void produce(edm::Event& e, const edm::EventSetup& c) override;
  void endJob();
  
private:

  edm::EDGetTokenT<std::vector<HBHERecHitTruth>> RHTToken;

  edm::EDGetTokenT<edm::PCaloHitContainer> fSHitToken;
  edm::EDGetTokenT<edm::SimTrackContainer> SimTrackToken;
  edm::EDGetTokenT<edm::SimVertexContainer> SimVertexToken;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> GenParticleToken;

  edm::EDGetTokenT<reco::PFClusterCollection> _clustersLabel;

  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit> > _rechitsLabel;
  edm::ESGetToken<HcalDDDRecConstants, HcalRecNumberingRecord> hdrcToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;

  //edm::EDGetTokenT<reco::PFClusterCollection> _clustersLabel;

  const HcalDDDRecConstants *fRecNumber;
  HcalSimParameterMap fSimParameterMap;
  CaloHitResponse* fResponse;

  double fMinLongEnergy, fMinShortEnergy, fLongShortSlope, fLongShortOffset;
  bool fUpgrade;
  bool fLegacy;
  bool fHFEMB;
  edm::ParameterSet fLongShortCut;
  //const HcalDDDRecConstants *fRecNumber;
  const CaloGeometry* fGeometry;
  TFile *fFile;
  TTree *fTree;


  const static int max_simPart = 20000;
  UInt_t value_simPart_n;
  int value_simPart_ID[max_simPart];
  float value_simPart_energy[max_simPart];
  float value_simPart_eta[max_simPart];
  int value_simPart_PID[max_simPart];
  float value_simPart_clusterEnergy[max_simPart];



};

#endif
