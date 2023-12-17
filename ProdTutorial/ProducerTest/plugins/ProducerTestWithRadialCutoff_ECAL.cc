// -*- C++ -*-
//
// Package:    ProdTutorial/ProducerTestWithRadialCutoff_ECAL
// Class:      ProducerTestWithRadialCutoff_ECAL
//
/**\class ProducerTestWithRadialCutoff_ECAL ProducerTestWithRadialCutoff_ECAL.cc ProdTutorial/ProducerTest/plugins/ProducerTestWithRadialCutoff_ECAL.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  William McCormack
//         Created:  Wed, 11 May 2022 21:10:34 GMT
//
//

// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "DataFormats/DetId/interface/DetId.h"

#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterBuilderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/RecHitTopologicalCleanerBase.h"

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

#include "DataFormats/EcalDetId/interface/EcalScDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHitTruth.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"

#include "SimCalorimetry/CaloSimAlgos/interface/CaloHitResponse.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalSimParameterMap.h"

//
// class declaration
//

class ProducerTestWithRadialCutoff_ECAL : public edm::stream::EDProducer<> {
public:
  explicit ProducerTestWithRadialCutoff_ECAL(const edm::ParameterSet&);
  ~ProducerTestWithRadialCutoff_ECAL() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  const HcalDDDRecConstants *fRecNumber;
  //const CaloGeometry* fGeometry;
  edm::EDGetTokenT<edm::PCaloHitContainer> fSHitTokenHCAL;
  edm::EDGetTokenT<edm::PCaloHitContainer> fSHitToken;
  edm::EDGetTokenT<edm::PCaloHitContainer> fSHitTokenEE;
  edm::EDGetTokenT<edm::PCaloHitContainer> fSHitTokenES;
  edm::EDGetTokenT<edm::SimTrackContainer> SimTrackToken;
  edm::EDGetTokenT<edm::SimVertexContainer> SimVertexToken;
  //edm::EDGetTokenT<reco::PFRecHitCollection> _rechitsLabel;
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit> > _rechitsLabel;
  edm::ESGetToken<HcalDDDRecConstants, HcalRecNumberingRecord> hdrcToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;

  //std::string fPFClustHName;
  //edm::EDGetTokenT<reco::PFClusterCollection>   fTokPFClustHName;
  edm::EDGetTokenT<reco::PFClusterCollection> _clustersLabel;

  typedef std::vector<HBHERecHitTruth> HBHERecHitTruthCollection;

  HcalSimParameterMap fSimParameterMap;
  //EcalSimParameterMap fSimParameterMap;
  CaloHitResponse* fResponse;

  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ProducerTestWithRadialCutoff_ECAL::ProducerTestWithRadialCutoff_ECAL(const edm::ParameterSet& iConfig){
  //register your products
  /* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
  */

  produces<HBHERecHitTruthCollection>("HBHERecHitTruth");

  //now do what ever other initialization is needed
  fSHitTokenHCAL = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));
  fSHitToken = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","EcalHitsEB"));
  fSHitTokenEE = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","EcalHitsEE"));
  fSHitTokenES = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","EcalHitsES"));
  SimTrackToken = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits",""));
  SimVertexToken = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits",""));
  //_rechitsLabel = consumes<reco::PFRecHitCollection>(iConfig.getParameter<edm::InputTag>("pf_src"));
  _rechitsLabel = consumes<edm::SortedCollection<HBHERecHit> >(iConfig.getParameter<edm::InputTag>("pf_src"));
  //hdrcToken_ = iC.esConsumes<edm::Transition::BeginRun>();
  hdrcToken_ = esConsumes<HcalDDDRecConstants, HcalRecNumberingRecord>();
  //hdrcToken_ = consumes<edm::Transition::BeginRun>();
  //hdrcToken_ =

  //fPFClustHName = iConfig.getUntrackedParameter<std::string>("edmName","particleFlowClusterHCAL");
  //fTokPFClustHName    = consumes<reco::PFClusterCollection>    (fPFClustHName);
  _clustersLabel = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("clus_src"));

  geomToken_ = esConsumes();

}

ProducerTestWithRadialCutoff_ECAL::~ProducerTestWithRadialCutoff_ECAL() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty

  //fSHitToken = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));

}

//
// member functions
//

// ------------ method called to produce the data  ------------
void ProducerTestWithRadialCutoff_ECAL::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  /* This is an event example
  //Read 'ExampleData' from the Event
  ExampleData const& in = iEvent.get(inToken_);

  //Use the ExampleData to create an ExampleData2 which 
  // is put into the Event
  iEvent.put(std::make_unique<ExampleData2>(in));
  */

  /* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  SetupData& setup = iSetup.getData(setupToken_);
  */

  edm::ESHandle<CaloGeometry> geoHandle = iSetup.getHandle(geomToken_);
  const CaloSubdetectorGeometry* hcalBarrelGeo = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  const CaloSubdetectorGeometry* hcalEndcapGeo = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalEndcap);

  const HcalDDDRecConstants* pHRNDC;
  //const edm::ESGetToken<HcalDDDRecConstants, HcalRecNumberingRecord> hdrcToken_;
  pHRNDC = &iSetup.getData(hdrcToken_);
  fRecNumber = &*pHRNDC;
  //const HcalDDDRecConstants fRecNumber = iSetup.getData(hdrcToken_);


  edm::Handle<edm::SimTrackContainer> hSimTracks;      // create handle
  iEvent.getByToken(SimTrackToken, hSimTracks);   // SimTracks

  edm::Handle<edm::SimVertexContainer> hSimVertices;      // create handle
  iEvent.getByToken(SimVertexToken, hSimVertices);   // SimVertices

  edm::Handle<edm::PCaloHitContainer> hSimHitsHCAL;      // create handle
  iEvent.getByToken(fSHitTokenHCAL, hSimHitsHCAL);   // SimHits

  edm::Handle<edm::PCaloHitContainer> hSimHits;      // create handle
  iEvent.getByToken(fSHitToken, hSimHits);   // SimHits

  edm::Handle<edm::PCaloHitContainer> hSimHitsEE;      // create handle
  iEvent.getByToken(fSHitTokenEE, hSimHitsEE);   // SimHits

  edm::Handle<edm::PCaloHitContainer> hSimHitsES;      // create handle
  iEvent.getByToken(fSHitTokenES, hSimHitsES);   // SimHits


  edm::PCaloHitContainer lSimHitsHCAL  = *hSimHitsHCAL;

  edm::PCaloHitContainer lSimHits  = *hSimHits;
  edm::PCaloHitContainer lSimHitsEE  = *hSimHitsEE;
  edm::PCaloHitContainer lSimHitsES  = *hSimHitsES;

  std::map<unsigned int, int> detId_map; //maps keys = detectorID and values = Sim Track geantID
  std::map<unsigned int, int> simhitIndex_map; //maps keys = detectorID and values = simHit index that already corresponds to that det ID
  std::map<unsigned int, float> simHit_energy_map; //maps keys = detectorID and values = Sim Track energy
  std::map<unsigned int, int> simHit_depth_map; //maps keys = detectorID and values = Sim Track depth
  std::map<unsigned int, float> simHit_time_map; //maps keys = detectorID and values = Sim Track time

  std::map<unsigned int, int> detId_map_2; //maps keys = detectorID and values = Sim Track geantID
  std::map<unsigned int, int> simhitIndex_map_2; //maps keys = detectorID and values = simHit index that already corresponds to that det ID
  std::map<unsigned int, float> simHit_energy_map_2; //maps keys = detectorID and values = Sim Track energy
  std::map<unsigned int, int> simHit_depth_map_2; //maps keys = detectorID and values = Sim Track depth
  std::map<unsigned int, float> simHit_time_map_2; //maps keys = detectorID and values = Sim Track time

  std::map<unsigned int, int> detId_map_3; //maps keys = detectorID and values = Sim Track geantID
  std::map<unsigned int, int> simhitIndex_map_3; //maps keys = detectorID and values = simHit index that already corresponds to that det ID
  std::map<unsigned int, float> simHit_energy_map_3; //maps keys = detectorID and values = Sim Track energy
  std::map<unsigned int, int> simHit_depth_map_3; //maps keys = detectorID and values = Sim Track depth
  std::map<unsigned int, float> simHit_time_map_3; //maps keys = detectorID and values = Sim Track time

  int largestTrackID = 0;

  fResponse = new CaloHitResponse(NULL, (CaloShapes*)NULL);
  fResponse->setGeometry(&*geoHandle);

  float totalSimEnergy = 0.;

  std::cout<<"loop over simhits"<<std::endl;
  for (int j=0; j < (int) lSimHits.size(); j++) {
    double samplingFactor = 1;
    //std::cout<<"pre relabel (lSimHits)[j].id() = "<<(lSimHits)[j].id()<<std::endl;
    HcalDetId simId = HcalHitRelabeller::relabel((lSimHits)[j].id(), fRecNumber);
    //std::cout<<"about to test subdet "<<simId.subdet()<<std::endl;
    //samplingFactor = fSimParameterMap.simParameters().samplingFactor(simId);
    //std::cout<<"samplingFactor? "<<samplingFactor<<std::endl;
    /*
    if(simId.subdet() == HcalBarrel) {
      samplingFactor = fSimParameterMap.hbParameters().samplingFactor(simId);
    } else if (simId.subdet() == HcalEndcap) {
      samplingFactor = fSimParameterMap.heParameters().samplingFactor(simId);
    }
    */
    //std::cout<<"j = "<<j<<" simID "<<simId.rawId()<<" geantTrackId() = "<<(lSimHits)[j].geantTrackId()<<" time = "<<(lSimHits)[j].time()<<" energy = "<<(lSimHits)[j].energy()<<" detid = "<<(lSimHits)[j].id()<<" depth = "<<(lSimHits)[j].depth()<<std::endl;
    totalSimEnergy += (lSimHits)[j].energy();
    if(!detId_map.count(simId.rawId())){ //new detID
      //std::cout<<"in first if"<<std::endl;
      detId_map[simId.rawId()] = (lSimHits)[j].geantTrackId(); //associate simId to highest energy geantTrackId

      if((lSimHits)[j].geantTrackId() > largestTrackID) largestTrackID = (lSimHits)[j].geantTrackId();

      simhitIndex_map[simId.rawId()] = j;
      //std::cout<<"(lSimHits)[j].geantTrackId() = "<<(lSimHits)[j].geantTrackId()<<std::endl;
      simHit_energy_map[simId.rawId()] = samplingFactor*((lSimHits)[j].energy());
      //std::cout<<"samplingFactor = "<<samplingFactor<<" samplingFactor*((lSimHits)[j].energy()) = "<<samplingFactor*((lSimHits)[j].energy())<<std::endl;
      HcalDetId simIdForDepth = (lSimHits)[j].id();
      //simHit_depth_map[simId.rawId()] = ((lSimHits)[j].id()).depth();
      simHit_depth_map[simId.rawId()] = simIdForDepth.depth();
      //std::cout<<"simIdForDepth.depth() = "<<simIdForDepth.depth()<<std::endl;
      //std::cout<<"(lSimHits)[j].time() = "<<(lSimHits)[j].time()<<" fResponse->timeOfFlight((lSimHits)[j].id()) = "<<fResponse->timeOfFlight((lSimHits)[j].id())<<std::endl;
      //DetId detId((lSimHits)[j].id());
      //std::cout<<"(lSimHits)[j].time() = "<<(lSimHits)[j].time()<<" fResponse->timeOfFlight((lSimHits)[j].id()) = "<<fResponse->timeOfFlight(detId)<<std::endl;
      //simHit_time_map[simId.rawId()] = (lSimHits)[j].time() - fResponse->timeOfFlight((lSimHits)[j].id());
      //simHit_time_map[simId.rawId()] = (lSimHits)[j].time() - fResponse->timeOfFlight(detId);
    } 
    else if((lSimHits)[j].energy() > (lSimHits)[simhitIndex_map[simId.rawId()]].energy() ){
      //else if( (lSimHits)[j].energy() > (lSimHits)[detId_map[simId.rawId()]].energy() ){ //higher energy deposit than previous max
      //std::cout<<"in elseif"<<std::endl;
      detId_map[simId.rawId()] = (lSimHits)[j].geantTrackId();

      if((lSimHits)[j].geantTrackId() > largestTrackID) largestTrackID = (lSimHits)[j].geantTrackId();

      simhitIndex_map[simId.rawId()] = j;
      simHit_energy_map[simId.rawId()] += samplingFactor*((lSimHits)[j].energy());
      HcalDetId simIdForDepth = (lSimHits)[j].id();
      simHit_depth_map[simId.rawId()] = simIdForDepth.depth();
      //simHit_time_map[simId.rawId()] = (lSimHits)[j].time() - fResponse->timeOfFlight((lSimHits)[j].id());
      //DetId detId((lSimHits)[j].id());
      //simHit_time_map[simId.rawId()] = (lSimHits)[j].time() - fResponse->timeOfFlight(detId);
    } else{ //lower energy deposit than max energy deposit
      //std::cout<<"in else"<<std::endl;
      simHit_energy_map[simId.rawId()] += samplingFactor*((lSimHits)[j].energy());
    }
  }

  std::cout<<"     totalSimEnergy EB = "<<totalSimEnergy<<std::endl;

  float totalSimEnergyEE = 0.;
  for (int j=0; j < (int) lSimHitsEE.size(); j++) {
    //double samplingFactor = 1;
    //std::cout<<"pre relabel (lSimHits)[j].id() = "<<(lSimHits)[j].id()<<std::endl;
    //HcalDetId simId = HcalHitRelabeller::relabel((lSimHits)[j].id(), fRecNumber);
    totalSimEnergyEE += (lSimHitsEE)[j].energy();
  }
  std::cout<<"     totalSimEnergy EE = "<<totalSimEnergyEE<<std::endl;

  float totalSimEnergyES = 0.;
  for (int j=0; j < (int) lSimHitsES.size(); j++) {
    //double samplingFactor = 1;
    //std::cout<<"pre relabel (lSimHits)[j].id() = "<<(lSimHits)[j].id()<<std::endl;
    //HcalDetId simId = HcalHitRelabeller::relabel((lSimHits)[j].id(), fRecNumber);
    totalSimEnergyES += (lSimHitsES)[j].energy();
  }



  std::cout<<"     totalSimEnergy ES = "<<totalSimEnergyES<<std::endl;


  float totalSimEnergyHB = 0.;
  float totalSimEnergyHB0 = 0.;
  float totalSimEnergyHB1 = 0.;
  float totalSimEnergyHB2 = 0.;
  float totalSimEnergyHB3 = 0.;
  float totalSimEnergyHB4 = 0.;
  float totalSimEnergyHB5 = 0.;
  float totalSimEnergyHB6 = 0.;
  float totalSimEnergyHB7 = 0.;
  float totalSimEnergyHB8 = 0.;
  float totalSimEnergyHB9 = 0.;
  float totalSimEnergyHB10 = 0.;
  float totalSimEnergyHB11 = 0.;
  float totalSimEnergyHB12 = 0.;
  float totalSimEnergyHB13 = 0.;
  float totalSimEnergyHB14 = 0.;
  float totalSimEnergyHB15 = 0.;
  float totalSimEnergyHB16 = 0.;
  float totalSimEnergyHE = 0.;
  for (int j=0; j < (int) lSimHitsHCAL.size(); j++) {
    double samplingFactor = 1;
    //std::cout<<"pre relabel (lSimHits)[j].id() = "<<(lSimHits)[j].id()<<std::endl;
    HcalDetId simId = HcalHitRelabeller::relabel((lSimHitsHCAL)[j].id(), fRecNumber);
    //std::cout<<"about to test subdet "<<simId.subdet()<<std::endl;
    if(simId.subdet() == HcalBarrel) {
      samplingFactor = fSimParameterMap.hbParameters().samplingFactor(simId);
      totalSimEnergyHB += samplingFactor*(lSimHitsHCAL)[j].energy();
      HcalDetId simIdForDepth = (lSimHitsHCAL)[j].id();
      std::cout<<"j = "<<j<<" simID "<<simId.rawId()<<" geantTrackId() = "<<(lSimHitsHCAL)[j].geantTrackId()<<" time = "<<(lSimHitsHCAL)[j].time()<<" energy = "<<(lSimHitsHCAL)[j].energy()<<" detid = "<<(lSimHitsHCAL)[j].id()<<" depth = "<<(lSimHitsHCAL)[j].depth()<<"     sampled e = "<<samplingFactor*(lSimHitsHCAL)[j].energy()<<"    other depth? "<<simIdForDepth.depth()<<std::endl;
      if(simIdForDepth.depth() == 0) totalSimEnergyHB0 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 1) totalSimEnergyHB1 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 2) totalSimEnergyHB2 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 3) totalSimEnergyHB3 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 4) totalSimEnergyHB4 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 5) totalSimEnergyHB5 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 6) totalSimEnergyHB6 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 7) totalSimEnergyHB7 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 8) totalSimEnergyHB8 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 9) totalSimEnergyHB9 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 10) totalSimEnergyHB10 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 11) totalSimEnergyHB11 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 12) totalSimEnergyHB12 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 13) totalSimEnergyHB13 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 14) totalSimEnergyHB14 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 15) totalSimEnergyHB15 += samplingFactor*(lSimHitsHCAL)[j].energy();
      if(simIdForDepth.depth() == 16) totalSimEnergyHB16 += samplingFactor*(lSimHitsHCAL)[j].energy();
    } else if (simId.subdet() == HcalEndcap) {
      samplingFactor = fSimParameterMap.heParameters().samplingFactor(simId);
      totalSimEnergyHE += samplingFactor*(lSimHitsHCAL)[j].energy();
    }
  }
  std::cout<<"     totalSimEnergy HB = "<<totalSimEnergyHB<<std::endl;
  std::cout<<"            totalSimEnergy HB0 = "<<totalSimEnergyHB0<<std::endl;
  std::cout<<"            totalSimEnergy HB1 = "<<totalSimEnergyHB1<<std::endl;
  std::cout<<"            totalSimEnergy HB2 = "<<totalSimEnergyHB2<<std::endl;
  std::cout<<"            totalSimEnergy HB3 = "<<totalSimEnergyHB3<<std::endl;
  std::cout<<"            totalSimEnergy HB4 = "<<totalSimEnergyHB4<<std::endl;
  std::cout<<"            totalSimEnergy HB5 = "<<totalSimEnergyHB5<<std::endl;
  std::cout<<"            totalSimEnergy HB6 = "<<totalSimEnergyHB6<<std::endl;
  std::cout<<"            totalSimEnergy HB7 = "<<totalSimEnergyHB7<<std::endl;
  std::cout<<"            totalSimEnergy HB8 = "<<totalSimEnergyHB8<<std::endl;
  std::cout<<"            totalSimEnergy HB9 = "<<totalSimEnergyHB9<<std::endl;
  std::cout<<"            totalSimEnergy HB10 = "<<totalSimEnergyHB10<<std::endl;
  std::cout<<"            totalSimEnergy HB11 = "<<totalSimEnergyHB11<<std::endl;
  std::cout<<"            totalSimEnergy HB12 = "<<totalSimEnergyHB12<<std::endl;
  std::cout<<"            totalSimEnergy HB13 = "<<totalSimEnergyHB13<<std::endl;
  std::cout<<"            totalSimEnergy HB14 = "<<totalSimEnergyHB14<<std::endl;
  std::cout<<"            totalSimEnergy HB15 = "<<totalSimEnergyHB15<<std::endl;
  std::cout<<"            totalSimEnergy HB16 = "<<totalSimEnergyHB16<<std::endl;
  std::cout<<"     totalSimEnergy HE = "<<totalSimEnergyHE<<std::endl;



  std::map<unsigned int, unsigned int> track_map; //map keys = Sim Track geantID and values = Sim Track index
  //int st=0;
  std::cout<<"loop over simtracks"<<std::endl;
  for (unsigned int i = 0; i < hSimTracks->size(); ++i) {
    std::cout<<"simTrack "<<i<<" geantID = "<<hSimTracks->at(i).genpartIndex()<<" vertID = "<<hSimTracks->at(i).vertIndex()<<" trackID = "<<hSimTracks->at(i).trackId()<<" energy = "<<hSimTracks->at(i).momentum()<<std::endl;
    track_map[hSimTracks->at(i).trackId()] = i;
    //i++;
  }

  //int sv=0;
  for (unsigned i = 0; i < hSimVertices->size(); ++i) {
    std::cout<<"simVertex "<<i<<" vertexID = "<<hSimVertices->at(i).vertexId()<<" parentID = "<<hSimVertices->at(i).parentIndex()<<" position = "<<hSimVertices->at(i).position().X()<<" "<<hSimVertices->at(i).position().Y()<<" "<<hSimVertices->at(i).position().Z()<<std::endl;
    //sv++;
  }


  //edm::Handle<reco::PFRecHitCollection> rechits;
  edm::Handle<edm::SortedCollection<HBHERecHit> > rechits;
  iEvent.getByToken(_rechitsLabel, rechits);
  //std::cout<<"SONIC rechits->size() = "<<rechits->size()<<std::endl;
  //const int sizeRH = rechits->size()<<std::endl;
  int sizeRH = 0;
  for (const auto& erh : *rechits) {
    if(erh.energy() < 0.8){
      continue;
    }
    sizeRH++;
  }


  //edm::Handle<reco::PFClusterCollection> hPFCHCAL;
  //iEvent.getByToken(fTokPFClustHName,hPFCHCAL);

  edm::Handle<reco::PFClusterCollection> hPFCHCAL;
  iEvent.getByToken(_clustersLabel, hPFCHCAL);
  assert(hPFCHCAL.isValid());
  //int pId = 0;
  std::cout<<"loop over hpfchcal"<<std::endl;
  std::map<HcalDetId,std::vector<std::pair<int, float>>> RH_PFClusterAndFrac;
  for(reco::PFClusterCollection::const_iterator itC = hPFCHCAL->begin(); itC!=hPFCHCAL->end(); itC++){
    const std::vector<reco::PFRecHitFraction>& hitsandfractions = itC->recHitFractions();
    //std::cout <<"cluster " << std::distance(hPFCHCAL->begin(),itC) << std::endl;
    for(long unsigned int it=0; it<hitsandfractions.size(); it++){
      const auto & RH = *hitsandfractions[it].recHitRef();
      HcalDetId pDetId(RH.detId());
      if (pDetId.det() != DetId::Hcal) continue;
      int cluster_id = std::distance(hPFCHCAL->begin(), itC);
      float hit_frac = hitsandfractions[it].fraction();
      RH_PFClusterAndFrac[pDetId].push_back(std::make_pair(cluster_id, hit_frac));
      //std::cout << pDetId.ieta() <<"/" << pDetId.iphi() << "/"<< RH.energy() << "/" << hitsandfractions[it].fraction() << std::endl;
    }
  }



  //std::vector<float> thresholdsHLT_HB = {0.0, 0.8, 0.8, 0.8, 0.8};
  //std::vector<float> thresholdsHLT_HE = {0.0, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
  std::vector<float> thresholdsReco_HB = {0.0, 0.125, 0.25, 0.35, 0.35}; //actually, I'm not sure what's going on here with these cuts.  need thorough investigation
  std::vector<float> thresholdsReco_HE = {0.0, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};


  //std::auto_ptr<HBHERecHitTruthCollection> hbheRHT( new HBHERecHitTruthCollection );
  auto hbheRHT = std::make_unique<HBHERecHitTruthCollection>();
  //hbheRHT->reserve( sizeRH );

  int h=0;
  //Loop for HBHE hits
  std::cout<<"loop over rechits"<<std::endl;
  //for (unsigned i = 0; i < nhits; ++i) {

  
  std::map<int, int> splitting_map; //keys are assigned parent IDs, values are the parent particle that the rechit was originally linked to
  std::map<int, double> parenteta; //
  std::map<int, double> parentphi; //
  std::cout<<"The largest geant ID was "<<largestTrackID<<std::endl;


  for (const auto& erh : *rechits) {
  //for (auto erh : *rechits) {

    bool cleanHLT = false;
    bool cleanReco = false;

    HcalDetId pId = erh.idFront();

    const HcalDetId detid = erh.idFront();
    HcalSubdetector esd = (HcalSubdetector)detid.subdetId();

    std::shared_ptr<const CaloCellGeometry> thisCell = nullptr;
    PFLayer::Layer layer = PFLayer::HCAL_BARREL1;
    switch (esd) {
    case HcalBarrel:
      thisCell = hcalBarrelGeo->getGeometry(detid);
      layer = PFLayer::HCAL_BARREL1;
      if(erh.energy() < 0.8){
	cleanHLT = false;
      } else{
	cleanHLT = true;
      }
      if(erh.energy() < thresholdsReco_HB.at(detid.depth())){
	//cleanReco = false;
	cleanReco = true;
      } else{
        cleanReco = true;
      }
      break;

    case HcalEndcap:
      thisCell = hcalEndcapGeo->getGeometry(detid);
      layer = PFLayer::HCAL_ENDCAP;
      if(erh.energy() < 0.8){
	cleanHLT = false;
      } else{
	cleanHLT = true;
      }
      if(erh.energy() < thresholdsReco_HE.at(detid.depth())){
	cleanReco = false;
      } else{
        cleanReco = true;
      }
      break;
    default:
      break;
    } 

    // find rechit geometry
    if (!thisCell) {
      edm::LogError("PFHBHERecHitCreator") << "warning detid " << std::hex << detid.rawId() << std::dec << " "
					   << detid << " not found in geometry" << std::endl;
      continue;
    }


    //CMS PF cluster matching with fractions
    std::vector<std::pair<int,float>> PFClusterMatches = RH_PFClusterAndFrac[pId];
    std::sort(PFClusterMatches.begin(), PFClusterMatches.end(), [](auto &left, auto &right){return left.second > right.second;});    
    int PFcluster0Id = -1;
    float PFcluster0frac = -1.;
    int PFcluster1Id = -1;
    float PFcluster1frac = -1.;
    int PFcluster2Id = -1;
    float PFcluster2frac = -1.;
    if (PFClusterMatches.size() > 0){
      //std::cout <<"clust "<<PFClusterMatches[0].first<< "   1/" << PFClusterMatches[0].second << std::endl;
      PFcluster0Id   = PFClusterMatches[0].first;
      PFcluster0frac = PFClusterMatches[0].second;
    }
    if (PFClusterMatches.size() > 1){
      //std::cout <<"clust "<<PFClusterMatches[1].first<< "   2/" << PFClusterMatches[1].second << std::endl;
      PFcluster1Id   = PFClusterMatches[1].first;
      PFcluster1frac = PFClusterMatches[1].second;
    } 
    if (PFClusterMatches.size() > 2){
      //std::cout <<"clust "<<PFClusterMatches[2].first<<  "   3/" << PFClusterMatches[2].second << std::endl;
      PFcluster2Id   = PFClusterMatches[2].first;
      PFcluster2frac = PFClusterMatches[2].second;
    }


    //std::cout<<"recHit energy = "<<erh.energy()<<" det id = "<<erh.detid().rawId()<<std::endl;
    if(!detId_map.count(erh.detid().rawId())){ //no simHit for recHit
      
      HBHERecHitTruth rht(erh, -1, thisCell->getPosition().basicVector().x(), thisCell->getPosition().basicVector().y(), thisCell->getPosition().basicVector().z(), thisCell->repPos().eta(), thisCell->repPos().phi(), detid.depth(), layer, cleanHLT, cleanReco, PFcluster0Id, PFcluster0frac, PFcluster1Id, PFcluster1frac, PFcluster2Id, PFcluster2frac);
      hbheRHT->push_back(rht);
      //std::cout<<"simHit-less recHit?"<<std::endl;

      h++;
      continue;
    }

    //std::cout<<"simVertex "<<i<<" vertexID = "<<hSimVertices->at(i).vertexId()<<" parentID = "<<hSimVertices->at(i).parentIndex()<<" position = "<<hSimVertices->at(i).position().X()<<" "<<hSimVertices->at(i).position().Y()<<" "<<hSimVertices->at(i).position().Z()<<std::endl;

    switch (esd) {
    case HcalBarrel:
      std::cout<<"Here is a BARREL recHit."<<std::endl;
      break;

    case HcalEndcap:
      std::cout<<"Here is an ENDCAP recHit."<<std::endl;
      break;
    default:
      break;
    }

    std::cout<<"recHit "<<h<<" energy "<<erh.energy()<<" time "<<erh.time()<<" simHit energy = "<<simHit_energy_map[erh.detid().rawId()]<<" simHit depth = "<<simHit_depth_map[erh.detid().rawId()]<<" simHit time = "<<simHit_time_map[erh.detid().rawId()]<<std::endl; //for HBHE rec hits
    std::cout<<thisCell->getPosition().basicVector().x()<<" "<<thisCell->getPosition().basicVector().y()<<" "<<thisCell->getPosition().basicVector().y()<<" ("<<thisCell->repPos().eta()<<", "<<thisCell->repPos().phi()<<")"<<std::endl;
    std::cout<<"NOW LOOKING FOR PARENT"<<std::endl;

    unsigned int steps = 0;
    //unsigned int vertValue = hSimTracks->at(track_map[detId_map[rechits->at(i).detId()]]).vertIndex();
    //unsigned int vertValue = hSimTracks->at(track_map[detId_map[rechits->at(i).detid().rawId()]]).vertIndex(); //for hbhe rec hits
    //std::cout<<"erh.detid().rawId() = "<<erh.detid().rawId()<<" detId_map[erh.detid().rawId()] = "<<detId_map[erh.detid().rawId()]<<" track_map[detId_map[erh.detid().rawId()]] = "<<track_map[detId_map[erh.detid().rawId()]]<<" hSimTracks->at(track_map[detId_map[erh.detid().rawId()]]).vertIndex() = "<<hSimTracks->at(track_map[detId_map[erh.detid().rawId()]]).vertIndex()<<std::endl;
    unsigned int vertValue = hSimTracks->at(track_map[detId_map[erh.detid().rawId()]]).vertIndex(); //for hbhe rec hits
    std::cout<<"    just a test of the trackersurfacemomentum stuff...  what is hSimTracks->at(track_map[detId_map[erh.detid().rawId()]]).trackerSurfaceMomentum()? "<<hSimTracks->at(track_map[detId_map[erh.detid().rawId()]]).trackerSurfaceMomentum()<<std::endl;
    //int parentPart = hSimVertices->at(vertValue).parentIndex();
    int parentPart = detId_map[erh.detid().rawId()];
    std::cout<<"first check verteval = "<<vertValue<<" parentPart = "<<parentPart<<" detId_map[erh.detid().rawId()] = "<<detId_map[erh.detid().rawId()]<<std::endl;
    std::cout<<"starting vertex at "<<hSimVertices->at(vertValue).position().X()<<" "<<hSimVertices->at(vertValue).position().Y()<<" "<<hSimVertices->at(vertValue).position().Z()<<std::endl;
    //if(vertValue == 0){
    if(parentPart == -1){
      parentPart = detId_map[erh.detid().rawId()];
    }
    //unsigned int gID = hSimTracks->at(track_map[detId_map[rechits->at(i).detId()]]).genpartIndex();
    //int particleOrigin = 0;
    //std::cout<<"verteValue = "<<vertValue<<" parentPart = "<<parentPart<<std::endl;
    //while(steps < 100 && vertValue > 0){
    while(steps < 100 && parentPart >= 0){
      vertValue = hSimTracks->at(track_map[parentPart]).vertIndex();
      bool ecp = (esd == HcalEndcap);
      bool brl = (esd == HcalBarrel);
      std::cout<<"new vertvalue = "<<vertValue<<".  Position is "<<hSimVertices->at(vertValue).position().X()<<", "<<hSimVertices->at(vertValue).position().Y()<<", "<<hSimVertices->at(vertValue).position().Z()<<".  parentIndex = "<<hSimVertices->at(vertValue).parentIndex()<<" endcap? "<<ecp<<" barrel? "<<brl<<std::endl;
      //if(vertValue == 0){
      if(hSimVertices->at(vertValue).parentIndex() < 0){
	break;
      }
      if(esd == HcalBarrel && sqrt(hSimVertices->at(vertValue).position().X()*hSimVertices->at(vertValue).position().X() + hSimVertices->at(vertValue).position().Y()*hSimVertices->at(vertValue).position().Y()) < 180 ){
      //if(esd == HcalBarrel && sqrt(hSimVertices->at(vertValue).position().X()*hSimVertices->at(vertValue).position().X() + hSimVertices->at(vertValue).position().Y()*hSimVertices->at(vertValue).position().Y()) < 120 ){ //ecal barrel
	break;
      }
      if(esd == HcalEndcap && abs(hSimVertices->at(vertValue).position().Z()) < 390 ){
	break;
      }
      std::cout<<"parent vertex at "<<hSimVertices->at(vertValue).position().X()<<" "<<hSimVertices->at(vertValue).position().Y()<<" "<<hSimVertices->at(vertValue).position().Z()<<std::endl;
      parentPart = hSimVertices->at(vertValue).parentIndex();
      std::cout<<"verteValue = "<<vertValue<<" parentPart = "<<parentPart<<std::endl;
      steps++;
    }
    
    std::cout<<"recHit "<<h<<" energy "<<erh.energy()<<" time "<<erh.time()<<" simHit energy = "<<simHit_energy_map[erh.detid().rawId()]<<" simHit depth = "<<simHit_depth_map[erh.detid().rawId()]<<" simHit time = "<<simHit_time_map[erh.detid().rawId()]<<" parent = "<<parentPart<<std::endl; //for HBHE rec hits
    if(parentPart>0){
    //  double deleta = (hSimTracks->at(track_map[parentPart]).momentum().Eta() - thisCell->repPos().eta());
    //  double delphi, ogphi, newphi;
    //  newphi = thisCell->repPos().phi();
    //  ogphi = hSimTracks->at(track_map[parentPart]).momentum().Phi();
    //  if(newphi < 0){
    //newphi = newphi + 6.28319;
    //  }
    //  if(ogphi < 0){
    //ogphi = ogphi + 6.28319;
    //  }
    //  std::cout<<"delta eta = "<<deleta<<" and delta phi = "<<delphi<<std::endl;
    //  double deltaR = sqrt( deleta*deleta + delphi*delphi );
    //  std::cout<<"RH with eta = "<<thisCell->repPos().eta()<<" and phi = "<<thisCell->repPos().phi()<<" matched to truth track with eta = "<<hSimTracks->at(track_map[parentPart]).momentum().Eta()<<" phi = "<<hSimTracks->at(track_map[parentPart]).momentum().Phi()<<" so deltaR = "<<deltaR<<std::endl;
    //}
    //erh.setParentPart(parentPart);

      bool checkingRH = true;
      while(parenteta.count(parentPart) && checkingRH){

	double deleta = (parenteta[parentPart] - thisCell->repPos().eta());
        double delphi, ogphi, newphi;
	newphi = thisCell->repPos().phi();
	ogphi = parentphi[parentPart];
	if(newphi < 0){
	  newphi = newphi + 6.28319;
	}
	if(ogphi < 0){
	  ogphi = ogphi + 6.28319;
	}
	delphi = abs(newphi - ogphi);
	if(delphi > 3.14159){
	  delphi = 3.14159 - std::fmod(delphi, 3.14159);
	}
	std::cout<<"delta eta = "<<deleta<<" and delta phi = "<<delphi<<std::endl;
	double deltaR = sqrt( deleta*deleta + delphi*delphi );
	if(deltaR > 1.){
	  if(splitting_map.count(parentPart)){
	    parentPart = splitting_map[parentPart];
	  } else{
	    largestTrackID++;
	    splitting_map[parentPart] = largestTrackID;
	    parentPart = largestTrackID;
	  }
	} else{
	  checkingRH = false;
	}

      }
  
      if(!parenteta.count(parentPart)){
	parenteta[parentPart] = thisCell->repPos().eta();
	parentphi[parentPart] = thisCell->repPos().phi();
      }
      
    }

    std::cout<<thisCell->getPosition().basicVector().x()<<" "<<thisCell->getPosition().basicVector().y()<<" "<<thisCell->getPosition().basicVector().y()<<" ("<<thisCell->repPos().eta()<<", "<<thisCell->repPos().phi()<<")"<<std::endl;


    
    HBHERecHitTruth rht(erh, parentPart, thisCell->getPosition().basicVector().x(), thisCell->getPosition().basicVector().y(), thisCell->getPosition().basicVector().z(), thisCell->repPos().eta(), thisCell->repPos().phi(), detid.depth(), layer, cleanHLT, cleanReco, PFcluster0Id, PFcluster0frac, PFcluster1Id, PFcluster1frac, PFcluster2Id, PFcluster2frac, simHit_energy_map[erh.detid().rawId()], simHit_depth_map[erh.detid().rawId()],simHit_time_map[erh.detid().rawId()]);
    //std::cout<<rht.parentPart()<<std::endl;

    hbheRHT->push_back(rht);

    //std::cout<<"recHit "<<i<<" position = ("<<rechits->at(i).position().x()<<", "<<rechits->at(i).position().y()<<", "<<rechits->at(i).position().z()<<") energy = "<<rechits->at(i).energy()<<" (eta, phi) = ("<<rechits->at(i).positionREP().eta()<<", "<<rechits->at(i).positionREP().phi()<<") time = "<<rechits->at(i).time()<<", id = "<<rechits->at(i).detId()<<" check detID map "<<detId_map[rechits->at(i).detId()]<<" and track index map "<<track_map[detId_map[rechits->at(i).detId()]]<<" parentPart = "<<parentPart<<std::endl;
    //rechits->at(i).setTruthLabel(parentPart);
    //std::cout<<"recHit "<<h<<" position = ("<<rechits->at(i).position().x()<<", "<<rechits->at(i).position().y()<<", "<<rechits->at(i).position().z()<<") energy = "rechits->at(i).energy()<<" (eta, phi) = ("<<rechits->at(i).positionREP().eta()<<", "<<rechits->at(i).positionREP().phi()<<") time = "<<rechits->at(i).time()<<", id = "<<rechits->at(i).id().det()<<" "<<rechits->at(i).id().subdetId()<<" raw id = "<<rechits->at(i).id().rawId()<<std::endl;
    //pfdata.push_back(rechits->at(i).position().x());
    //pfdata.push_back(rechits->at(i).position().y());
    //pfdata.push_back(rechits->at(i).position().z());
    //pfdata.push_back(rechits->at(i).energy());
    //pfdata.push_back(rechits->at(i).positionREP().eta());
    //pfdata.push_back(rechits->at(i).positionREP().phi());
    //pfdata.push_back(rechits->at(i).time());
    h++;
  }

  /*
    auto nhits = rechits->size();
  for (unsigned i = 0; i < nhits; ++i) {

    unsigned int steps = 0;
    unsigned int vertValue = hSimTracks->at(track_map[detId_map[rechits->at(i).detId()]]).vertIndex();
    unsigned int parentPart = hSimVertices->at(vertValue).parentIndex();
    //unsigned int gID = hSimTracks->at(track_map[detId_map[rechits->at(i).detId()]]).genpartIndex();
    //int particleOrigin = 0;
    //std::cout<<"verteValue = "<<vertValue<<" parentPart = "<<parentPart<<std::endl;
    while(steps < 100 && vertValue > 0){
      vertValue = hSimTracks->at(track_map[parentPart]).vertIndex();
      if(vertValue == 0){
	break;
      }
      parentPart = hSimVertices->at(vertValue).parentIndex();
      //std::cout<<"verteValue = "<<vertValue<<" parentPart = "<<parentPart<<std::endl;
      steps++;
    }
    

    std::cout<<"recHit "<<i<<" position = ("<<rechits->at(i).position().x()<<", "<<rechits->at(i).position().y()<<", "<<rechits->at(i).position().z()<<") energy = "<<rechits->at(i).energy()<<" (eta, phi) = ("<<rechits->at(i).positionREP().eta()<<", "<<rechits->at(i).positionREP().phi()<<") time = "<<rechits->at(i).time()<<", id = "<<rechits->at(i).detId()<<" check detID map "<<detId_map[rechits->at(i).detId()]<<" and track index map "<<track_map[detId_map[rechits->at(i).detId()]]<<" parentPart = "<<parentPart<<std::endl;
    //rechits->at(i).setTruthLabel(parentPart);
    
    //std::cout<<"recHit "<<h<<" position = ("<<rechits->at(i).position().x()<<", "<<rechits->at(i).position().y()<<", "<<rechits->at(i).position().z()<<") energy = "rechits->at(i).energy()<<" (eta, phi) = ("<<rechits->at(i).positionREP().eta()<<", "<<rechits->at(i).positionREP().phi()<<") time = "<<rechits->at(i).time()<<", id = "<<rechits->at(i).id().det()<<" "<<rechits->at(i).id().subdetId()<<" raw id = "<<rechits->at(i).id().rawId()<<std::endl;
    //pfdata.push_back(rechits->at(i).position().x());
    //pfdata.push_back(rechits->at(i).position().y());
    //pfdata.push_back(rechits->at(i).position().z());
    //pfdata.push_back(rechits->at(i).energy());
    //pfdata.push_back(rechits->at(i).positionREP().eta());
    //pfdata.push_back(rechits->at(i).positionREP().phi());
    //pfdata.push_back(rechits->at(i).time());
    //h++;
  }
*/


  //iEvent.put( hbheRHT, "HBHERecHitTruth" );
  iEvent.put( std::move(hbheRHT), "HBHERecHitTruth" );
  //iEvent.put( std::move(hbheRHT));

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void ProducerTestWithRadialCutoff_ECAL::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void ProducerTestWithRadialCutoff_ECAL::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
ProducerTestWithRadialCutoff_ECAL::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
ProducerTestWithRadialCutoff_ECAL::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ProducerTestWithRadialCutoff_ECAL::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ProducerTestWithRadialCutoff_ECAL::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ProducerTestWithRadialCutoff_ECAL::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
  desc.add<edm::InputTag>("pf_src");
}

//define this as a plug-in
DEFINE_FWK_MODULE(ProducerTestWithRadialCutoff_ECAL);
