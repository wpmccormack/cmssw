#ifndef HcalTrigPrimDigiNtupler_h
#define HcalTrigPrimDigiNtupler_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimCalorimetry/HcalTrigPrimAlgos/interface/HcalTriggerPrimitiveAlgo.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

//#include "CaloTrigNN/CaloNtupler/interface/FillerRH.hh"
//#include "CaloTrigNN/CaloNtupler/interface/FillerSimTrack.hh"
//#include "CaloTrigNN/CaloNtupler/interface/FillerSimVertex.hh"
//#include "CaloTrigNN/CaloNtupler/interface/FillerGenJets.hh"
//#include "CaloTrigNN/CaloNtupler/interface/FillerGenInfo.hh"
//#include "CaloTrigNN/CaloNtupler/interface/FillerEventInfo.hh"
////#include "CaloTrigNN/CaloNtupler/interface/FillerGenParticle.hh"
//#include "CaloTrigNN/DataFormats/interface/TRHPart.hh"
//#include "CaloTrigNN/DataFormats/interface/TGenJet.hh"
//#include "CaloTrigNN/DataFormats/interface/TEventInfo.hh"
//#include "CaloTrigNN/DataFormats/interface/TGenEventInfo.hh"
//#include "CaloTrigNN/DataFormats/interface/TGenParticle.hh"
//#include "CaloTrigNN/DataFormats/interface/TSimTrack.hh"
//#include "CaloTrigNN/DataFormats/interface/TSimVertex.hh"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "RecoParticleFlow/PFProducer/interface/PFAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
//#include "RecoParticleFlow/PFProducer/plugins/PFProducer.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHitTruth.h"


class HcalTrigPrimDigiNtupler : public edm::stream::EDProducer<> {
public:

  explicit HcalTrigPrimDigiNtupler(const edm::ParameterSet& ps);
  ~HcalTrigPrimDigiNtupler() override { endJob(); }//delete fTrackArr; delete fVertexArr; }// delete fFillerRH; delete fEvtArr; delete fFillerEventInfo; delete fFillerGenInfo; delete fGenParArr;}
  //~HcalTrigPrimDigiNtupler() override { endJob(); delete fPFParArr; delete fFillerPF; delete fRHParArr; delete fFillerRH;}
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

  double fMinLongEnergy, fMinShortEnergy, fLongShortSlope, fLongShortOffset;
  bool fUpgrade;
  bool fLegacy;
  bool fHFEMB;
  edm::ParameterSet fLongShortCut;
  const HcalDDDRecConstants *fRecNumber;
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
  float value_simHit_energy[max_hit];
  int value_simHit_depth[max_hit];
  float value_simHit_time[max_hit];


  //TClonesArray *fEnergyArr;

  //TClonesArray *fRHParArr;
  //TClonesArray *fGenJetArr;
  //TClonesArray *fGenFatJetArr;
  //TClonesArray *fGenParticleArr;
  //TClonesArray *fECALPFClusterArr;
  //TClonesArray *fGenEventInfoArr;
  //TClonesArray *fEventInfoArr;
  //baconhep::TEventInfo *fEvtInfo;
  //baconhep::TGenEventInfo *fGenEvtInfo;
  //TClonesArray *fTrackArr;
  //TClonesArray *fVertexArr;
};

#endif

