#ifndef HcalTrigPrimDigiNtuplerEXPANDEDINFO_h
#define HcalTrigPrimDigiNtuplerEXPANDEDINFO_h

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


class HcalTrigPrimDigiNtuplerEXPANDEDINFO : public edm::stream::EDProducer<> {
public:

  explicit HcalTrigPrimDigiNtuplerEXPANDEDINFO(const edm::ParameterSet& ps);
  ~HcalTrigPrimDigiNtuplerEXPANDEDINFO() override { endJob(); }//delete fTrackArr; delete fVertexArr; }// delete fFillerRH; delete fEvtArr; delete fFillerEventInfo; delete fFillerGenInfo; delete fGenParArr;}
  //~HcalTrigPrimDigiNtuplerEXPANDEDINFO() override { endJob(); delete fPFParArr; delete fFillerPF; delete fRHParArr; delete fFillerRH;}
  /**Produces the EDM products,*/
  void produce(edm::Event& e, const edm::EventSetup& c) override;
  void endJob();
  
private:
  //FillerRH  *fFillerRH;
  //FillerGenParticle *fFillerGenParticle;
  //baconhep::FillerGenJets *fFillerGenJets;
  //baconhep::FillerGenInfo *fFillerGenInfo;
  //baconhep::FillerEventInfo *fFillerEventInfo;
  //FillerSimTrack *fFillerSimTrack; 
  //FillerSimVertex *fFillerSimVertex; 

  /// input tags for HCAL digis
  //std::vector<edm::InputTag> fInputLabel;
  //std::vector<edm::InputTag> fInputUpgradeLabel;
  // this seems a strange way of doing things
  //edm::EDGetTokenT<edm::PCaloHitContainer> fSHitToken;
  //edm::EDGetTokenT<edm::PCaloHitContainer> fSHitTokenEE;
  //edm::EDGetTokenT<edm::PCaloHitContainer> fSHitTokenEB;
  //edm::EDGetTokenT<edm::SimTrackContainer> fSimTrackToken;
  //edm::EDGetTokenT<edm::SimVertexContainer> fSimVertexToken;

  edm::EDGetTokenT<std::vector<HBHERecHitTruth>> RHTToken;

  edm::EDGetTokenT<edm::PCaloHitContainer> fSHitToken;
  edm::EDGetTokenT<edm::SimTrackContainer> SimTrackToken;
  edm::EDGetTokenT<edm::SimVertexContainer> SimVertexToken;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> GenParticleToken;

  edm::EDGetTokenT<reco::PFClusterCollection> _clustersLabel;

  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit> > _rechitsLabel;
  edm::ESGetToken<HcalDDDRecConstants, HcalRecNumberingRecord> hdrcToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;

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

  const static int max_hit = 20000;
  UInt_t value_hit_n;
  float value_hit_energy[max_hit];
  float value_hit_x[max_hit];
  float value_hit_y[max_hit];
  float value_hit_z[max_hit];
  float value_hit_eta[max_hit];
  float value_hit_phi[max_hit];
  int value_hit_depth[max_hit];
  bool value_hit_isHB[max_hit];
  bool value_hit_cleanHLT[max_hit];
  bool value_hit_cleanReco[max_hit];
  float value_hit_time[max_hit];
  float value_hit_timeFalling[max_hit];
  int value_hit_PFcluster0Id[max_hit];
  float value_hit_PFcluster0frac[max_hit];
  int value_hit_PFcluster1Id[max_hit];
  float value_hit_PFcluster1frac[max_hit];
  int value_hit_PFcluster2Id[max_hit];
  float value_hit_PFcluster2frac[max_hit];
  int value_hit_parentPart[max_hit];
  float value_simHit_linkedEnergy[max_hit];
  //int value_simHit_depth[max_hit];
  //float value_simHit_time[max_hit];


  const static int max_simHit = 20000;
  UInt_t value_simHit_n;
  float value_simHit_energy[max_simHit];
  float value_simHit_x[max_simHit];
  float value_simHit_y[max_simHit];
  float value_simHit_z[max_simHit];
  float value_simHit_eta[max_simHit];
  float value_simHit_phi[max_simHit];
  int value_simHit_depth[max_simHit];

  const static int max_simHit_TRACEBACK = 20000;
  UInt_t value_simHit_TRACEBACK_n;
  int value_simHit_TRACEBACK_trackID[max_simHit_TRACEBACK];
  float value_simHit_TRACEBACK_energy[max_simHit_TRACEBACK];
  float value_simHit_TRACEBACK_x[max_simHit_TRACEBACK];
  float value_simHit_TRACEBACK_y[max_simHit_TRACEBACK];
  float value_simHit_TRACEBACK_z[max_simHit_TRACEBACK];
  float value_simHit_TRACEBACK_eta[max_simHit_TRACEBACK];
  float value_simHit_TRACEBACK_phi[max_simHit_TRACEBACK];
  int value_simHit_TRACEBACK_depth[max_simHit_TRACEBACK];

  const static int max_simTrack = 20000;
  UInt_t value_simTrack_n;
  float value_simTrack_energy[max_simTrack];
  float value_simTrack_x[max_simTrack];
  float value_simTrack_y[max_simTrack];
  float value_simTrack_z[max_simTrack];
  float value_simTrack_eta[max_simTrack];
  float value_simTrack_phi[max_simTrack];
  int value_simTrack_depth[max_simTrack];

  const static int max_simTrack_TRACEBACK = 20000;
  UInt_t value_simTrack_TRACEBACK_n;
  int value_simTrack_TRACEBACK_ID[max_simTrack_TRACEBACK];
  float value_simTrack_TRACEBACK_genEnergy[max_simTrack_TRACEBACK];
  int value_simTrack_TRACEBACK_genPDGID[max_simTrack_TRACEBACK];
  float value_simTrack_TRACEBACK_simmedEnergy[max_simTrack_TRACEBACK];
  float value_simTrack_TRACEBACK_x[max_simTrack_TRACEBACK];
  float value_simTrack_TRACEBACK_y[max_simTrack_TRACEBACK];
  float value_simTrack_TRACEBACK_z[max_simTrack_TRACEBACK];
  float value_simTrack_TRACEBACK_genEta[max_simTrack_TRACEBACK];
  float value_simTrack_TRACEBACK_genPhi[max_simTrack_TRACEBACK];
  int value_simTrack_TRACEBACK_depth[max_simTrack_TRACEBACK];
  int value_simTrack_TRACEBACK_hasHit_withDepth[max_simTrack_TRACEBACK][7];
  float value_simTrack_TRACEBACK_energy_withDepth[max_simTrack_TRACEBACK][7];
  float value_simTrack_TRACEBACK_avgPhi_withDepth[max_simTrack_TRACEBACK][7];
  float value_simTrack_TRACEBACK_avgEta_withDepth[max_simTrack_TRACEBACK][7];
  float value_simTrack_TRACEBACK_widthPhi_withDepth[max_simTrack_TRACEBACK][7];
  float value_simTrack_TRACEBACK_widthEta_withDepth[max_simTrack_TRACEBACK][7];
  float value_simTrack_TRACEBACK_EWavgPhi_withDepth[max_simTrack_TRACEBACK][7];
  float value_simTrack_TRACEBACK_EWavgEta_withDepth[max_simTrack_TRACEBACK][7];
  float value_simTrack_TRACEBACK_EWwidthPhi_withDepth[max_simTrack_TRACEBACK][7];
  float value_simTrack_TRACEBACK_EWwidthEta_withDepth[max_simTrack_TRACEBACK][7];




};

#endif

