// -*- C++ -*-
//
// Package:    ProdTutorial/ProducerTestWithRadialCutoffEXPANDEDINFO
// Class:      ProducerTestWithRadialCutoffEXPANDEDINFO
//
/**\class ProducerTestWithRadialCutoffEXPANDEDINFO ProducerTestWithRadialCutoffEXPANDEDINFO.cc ProdTutorial/ProducerTest/plugins/ProducerTestWithRadialCutoffEXPANDEDINFO.cc

 Description: Create rechit object with expanded information: HBHERecHitTruth, which has truth-labelling information for rechits.  HBHERecHitTruth object contains extra information used for training SPVCNN and potential variations on this algorithm

 Implementation:
     Truth-matching scheme still under development
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
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHitTruth.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"

#include "SimCalorimetry/CaloSimAlgos/interface/CaloHitResponse.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"

//
// class declaration
//

class ProducerTestWithRadialCutoffEXPANDEDINFO : public edm::stream::EDProducer<> {
public:
  explicit ProducerTestWithRadialCutoffEXPANDEDINFO(const edm::ParameterSet&);
  ~ProducerTestWithRadialCutoffEXPANDEDINFO() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  const HcalDDDRecConstants *fRecNumber;
  edm::EDGetTokenT<edm::PCaloHitContainer> fSHitToken;
  edm::EDGetTokenT<edm::SimTrackContainer> SimTrackToken;
  edm::EDGetTokenT<edm::SimVertexContainer> SimVertexToken;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> GenParticleToken;
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit> > _rechitsLabel;
  edm::ESGetToken<HcalDDDRecConstants, HcalRecNumberingRecord> hdrcToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;

  edm::EDGetTokenT<reco::PFClusterCollection> _clustersLabel;

  typedef std::vector<HBHERecHitTruth> HBHERecHitTruthCollection;

  HcalSimParameterMap fSimParameterMap;
  CaloHitResponse* fResponse;

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
ProducerTestWithRadialCutoffEXPANDEDINFO::ProducerTestWithRadialCutoffEXPANDEDINFO(const edm::ParameterSet& iConfig){
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
  fSHitToken = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));
  SimTrackToken = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits",""));
  SimVertexToken = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits",""));
  GenParticleToken = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles",""));
  _rechitsLabel = consumes<edm::SortedCollection<HBHERecHit> >(iConfig.getParameter<edm::InputTag>("pf_src"));
  hdrcToken_ = esConsumes<HcalDDDRecConstants, HcalRecNumberingRecord>();

  //use commeted out lines if you want to get some exotic cluster collection
  //fPFClustHName = iConfig.getUntrackedParameter<std::string>("edmName","particleFlowClusterHCAL");
  //fTokPFClustHName    = consumes<reco::PFClusterCollection>    (fPFClustHName);
  _clustersLabel = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("clus_src"));

  geomToken_ = esConsumes();

}

ProducerTestWithRadialCutoffEXPANDEDINFO::~ProducerTestWithRadialCutoffEXPANDEDINFO() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //

}

//
// member functions
//

// ------------ method called to produce the data  ------------
void ProducerTestWithRadialCutoffEXPANDEDINFO::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
  pHRNDC = &iSetup.getData(hdrcToken_);
  fRecNumber = &*pHRNDC;


  edm::Handle<edm::SimTrackContainer> hSimTracks;      // create handle
  iEvent.getByToken(SimTrackToken, hSimTracks);   // SimTracks

  edm::Handle<edm::SimVertexContainer> hSimVertices;      // create handle
  iEvent.getByToken(SimVertexToken, hSimVertices);   // SimVertices

  edm::Handle<edm::PCaloHitContainer> hSimHits;      // create handle
  iEvent.getByToken(fSHitToken, hSimHits);   // SimHits

  edm::Handle<std::vector<reco::GenParticle>> hGenParticles;      // create handle
  iEvent.getByToken(GenParticleToken, hGenParticles);   // SimHits

  /*
  //Print out generated particles for reference
  for (unsigned i = 0; i < hGenParticles->size(); ++i) {
    std::cout<<"GenParticle "<<i<<" pdgid "<<hGenParticles->at(i).pdgId()<<" energy "<<hGenParticles->at(i).energy()<<" (px "<<hGenParticles->at(i).px()<<" py "<<hGenParticles->at(i).py()<<" pz "<<hGenParticles->at(i).pz()<<") eta "<<hGenParticles->at(i).eta()<<", phi "<<hGenParticles->at(i).phi()<<" (vx "<<hGenParticles->at(i).vx()<<" vy "<<hGenParticles->at(i).vy()<<" vz "<<hGenParticles->at(i).vz()<<")"<<" numberOfDaughters() "<<hGenParticles->at(i).numberOfDaughters()<<" numberOfMothers() "<<hGenParticles->at(i).numberOfMothers()<<" isPromptFinalState() "<<hGenParticles->at(i).isPromptFinalState()<<" isPromptDecayed() "<<hGenParticles->at(i).isPromptDecayed()<<" isHardProcess() "<<hGenParticles->at(i).isHardProcess()<<" fromHardProcessFinalState() "<<hGenParticles->at(i).fromHardProcessFinalState()<<" fromHardProcessDecayed() "<<hGenParticles->at(i).fromHardProcessDecayed()<<" isLastCopy() "<<hGenParticles->at(i).isLastCopy()<<std::endl;
  }
  */


  edm::PCaloHitContainer lSimHits  = *hSimHits;

  std::map<unsigned int, std::vector<int>> detId_map; //maps keys = detectorID and values = Sim Track geantID
  std::map<unsigned int, std::vector<int>> simhitIndex_map; //maps keys = detectorID and values = simHit index that already corresponds to that det ID
  std::map<unsigned int, std::vector<float>> simHit_energy_map; //maps keys = detectorID and values = Sim Track energy
  std::map<unsigned int, std::vector<int>> simHit_depth_map; //maps keys = detectorID and values = Sim Track depth
  std::map<unsigned int, std::vector<float>> simHit_time_map; //maps keys = detectorID and values = Sim Track time
  std::map<int, float> trackID_energy_map; //map keys = track ID and values = total simhit energy for that particle

  int largestTrackID = 0;

  fResponse = new CaloHitResponse(NULL, (CaloShapes*)NULL);
  fResponse->setGeometry(&*geoHandle);

  for (int j=0; j < (int) lSimHits.size(); j++) { //loop over simulated hits
    double samplingFactor = 1;
    HcalDetId simId = HcalHitRelabeller::relabel((lSimHits)[j].id(), fRecNumber);
    //determine calorimeter sampling fraction based on HCAL component
    if(simId.subdet() == HcalBarrel) {
      samplingFactor = fSimParameterMap.hbParameters().samplingFactor(simId);
    } else if (simId.subdet() == HcalEndcap) {
      samplingFactor = fSimParameterMap.heParameters().samplingFactor(simId);
    }

    HcalDetId simIdForDepth = (lSimHits)[j].id(); //HcalDetId stores HCAL depth information in the correct format
    
    trackID_energy_map[(lSimHits)[j].geantTrackId()] += samplingFactor * (lSimHits)[j].energy(); //increment energy deposited by sim particle in question (index by GEANT ID)

    detId_map[simId.rawId()].push_back((lSimHits)[j].geantTrackId());
    simhitIndex_map[simId.rawId()].push_back(j);
    simHit_energy_map[simId.rawId()].push_back(samplingFactor*((lSimHits)[j].energy()));
    simHit_depth_map[simId.rawId()].push_back(simIdForDepth.depth());
    
    if((lSimHits)[j].geantTrackId() > largestTrackID) largestTrackID = (lSimHits)[j].geantTrackId();

  }

  //hash map to get tracks based on their trackID, rather than having to loop over track container to find correct track
  std::map<unsigned int, unsigned int> track_map; //map keys = Sim Track geantID and values = Sim Track index
  for (unsigned int i = 0; i < hSimTracks->size(); ++i) {
    //retaining comment as reference for potentially relevant track content
    /*std::cout<<"simTrack "<<i<<" geantID = "<<hSimTracks->at(i).genpartIndex()<<" vertID = "<<hSimTracks->at(i).vertIndex()<<" trackID = "<<hSimTracks->at(i).trackId()<<" type = "<<hSimTracks->at(i).type()<<" energy = "<<hSimTracks->at(i).momentum().E()<<" momentum = ("<<hSimTracks->at(i).momentum().X()<<", "<<hSimTracks->at(i).momentum().Y()<<", "<<hSimTracks->at(i).momentum().Z()<<") and genpartIndex = "<<hSimTracks->at(i).genpartIndex()<<" tracker surface position = "<<hSimTracks->at(i).trackerSurfacePosition().X()<<", "<<hSimTracks->at(i).trackerSurfacePosition().Y()<<", "<<hSimTracks->at(i).trackerSurfacePosition().Z()<<" tracker surface momentum = "<<hSimTracks->at(i).trackerSurfaceMomentum().X()<<", "<<hSimTracks->at(i).trackerSurfaceMomentum().Y()<<", "<<hSimTracks->at(i).trackerSurfaceMomentum().Z()<<" crossed boundary? "<<hSimTracks->at(i).crossedBoundary()<<" getIDAtBoundary = "<<hSimTracks->at(i).getIDAtBoundary()<<" momentum at boundary = "<<hSimTracks->at(i).getMomentumAtBoundary().X()<<", "<<hSimTracks->at(i).getMomentumAtBoundary().Y()<<", "<<hSimTracks->at(i).getMomentumAtBoundary().Z()<<std::endl;*/
    track_map[hSimTracks->at(i).trackId()] = i;
  }

  /*
  //Can be useful to look at simulated vertices for reference
  //int sv=0;
  for (unsigned i = 0; i < hSimVertices->size(); ++i) {
    std::cout<<"simVertex "<<i<<" vertexID = "<<hSimVertices->at(i).vertexId()<<" parentID = "<<hSimVertices->at(i).parentIndex()<<" position = "<<hSimVertices->at(i).position().X()<<" "<<hSimVertices->at(i).position().Y()<<" "<<hSimVertices->at(i).position().Z()<<std::endl;
    //sv++;
  }
  */


  edm::Handle<edm::SortedCollection<HBHERecHit> > rechits;
  iEvent.getByToken(_rechitsLabel, rechits);
  int sizeRH = 0;
  for (const auto& erh : *rechits) {
    if(erh.energy() < 0.8){
      continue;
    }
    sizeRH++;
  }


  edm::Handle<reco::PFClusterCollection> hPFCHCAL;
  iEvent.getByToken(_clustersLabel, hPFCHCAL);
  assert(hPFCHCAL.isValid());

  //use simhitIndex_map to find simhits in each cluster
  std::map<HcalDetId,std::vector<std::pair<int, float>>> RH_PFClusterAndFrac;
  std::map<int, std::vector<std::vector<float>>> simID_toClustEnEtaPhi;
  for(reco::PFClusterCollection::const_iterator itC = hPFCHCAL->begin(); itC!=hPFCHCAL->end(); itC++){ //loop over clusters in event

    std::map<int, float> cluster_geantid_energy_map; //key = simulated particle geant id; value = energy contributed to the cluster
    
    const std::vector<reco::PFRecHitFraction>& hitsandfractions = itC->recHitFractions();
    for(long unsigned int it=0; it<hitsandfractions.size(); it++){ //loop over rechits in cluster
      const auto & RH = *hitsandfractions[it].recHitRef();

      HcalDetId pDetId(RH.detId());
      if (pDetId.det() != DetId::Hcal) continue;
      if(simHit_energy_map.count(pDetId.rawId())){
	
	for( unsigned int i = 0; i < simHit_energy_map[pDetId.rawId()].size(); i++ ){ //Loop over simulated hits in rechit cell.  Length of simHit_energy_map[pDetId.rawId()] vector will be the same as e.g. simhitIndex_map[pDetId.rawId()], which gives us the index of the simhit in question from the simhit collection
	  std::map<int, float> rh_geantid_energy_map;

	  if( !cluster_geantid_energy_map.count((lSimHits)[simhitIndex_map[pDetId.rawId()][i]].geantTrackId()) ){
	    cluster_geantid_energy_map[(lSimHits)[simhitIndex_map[pDetId.rawId()][i]].geantTrackId()] = hitsandfractions[it].fraction()*simHit_energy_map[pDetId.rawId()][i];
	  }
	  else{
	    cluster_geantid_energy_map[(lSimHits)[simhitIndex_map[pDetId.rawId()][i]].geantTrackId()] += hitsandfractions[it].fraction()*simHit_energy_map[pDetId.rawId()][i];
	  }

	  if( !rh_geantid_energy_map.count((lSimHits)[simhitIndex_map[pDetId.rawId()][i]].geantTrackId()) ){
	    rh_geantid_energy_map[(lSimHits)[simhitIndex_map[pDetId.rawId()][i]].geantTrackId()] = hitsandfractions[it].fraction()*simHit_energy_map[pDetId.rawId()][i];
	  }
	  else{
	    rh_geantid_energy_map[(lSimHits)[simhitIndex_map[pDetId.rawId()][i]].geantTrackId()] += hitsandfractions[it].fraction()*simHit_energy_map[pDetId.rawId()][i];
	  }

	}
      }



      int cluster_id = std::distance(hPFCHCAL->begin(), itC);
      float hit_frac = hitsandfractions[it].fraction();
      RH_PFClusterAndFrac[pDetId].push_back(std::make_pair(cluster_id, hit_frac));
    }
    
    std::map<int,float>::iterator r;
    for (r = cluster_geantid_energy_map.begin(); r!= cluster_geantid_energy_map.end(); r++){
      
      std::vector<std::vector<float>> dummy;
      std::vector<float> clusenetaphi;
      clusenetaphi.push_back(std::distance(hPFCHCAL->begin(), itC));
      clusenetaphi.push_back(itC->energy());
      clusenetaphi.push_back(itC->positionREP().eta());
      clusenetaphi.push_back(itC->positionREP().phi());
      clusenetaphi.push_back((*r).second);
      
      if( !simID_toClustEnEtaPhi.count((*r).first) ){
	dummy.push_back(clusenetaphi);
	simID_toClustEnEtaPhi[(*r).first] = dummy;
      }
      else{
	simID_toClustEnEtaPhi[(*r).first].push_back(clusenetaphi);
      }
    }
    
  } //end loop over clusters
  
  /*
  //For each simulated particles, you can look at which clusters it contributed to and how much energy it contributed
  std::cout<<"LOOP OVER SIM IDS"<<std::endl;
  std::map<int,std::vector<std::vector<float>>>::iterator r;
  for (r = simID_toClustEnEtaPhi.begin(); r!= simID_toClustEnEtaPhi.end(); r++){
    std::cout<<(*r).first<<std::endl;
    for(unsigned int cl = 0; cl<(*r).second.size(); cl++){
      std::cout<<"    clust "<<(*r).second.at(cl).at(0)<<" en "<<(*r).second.at(cl).at(1)<<" eta "<<(*r).second.at(cl).at(2)<<" phi "<<(*r).second.at(cl).at(3)<<" sim en in clust = "<<(*r).second.at(cl).at(4)<<" which is "<<(*r).second.at(cl).at(4)/trackID_energy_map[(*r).first]<<"\n";
    }
  }
  */
  


  //if you want to use HLT-level thresholds only
  //std::vector<float> thresholdsHLT_HB = {0.0, 0.8, 0.8, 0.8, 0.8};
  //std::vector<float> thresholdsHLT_HE = {0.0, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};

  //if you want to use offline-level thresholds
  std::vector<float> thresholdsReco_HB = {0.0, 0.125, 0.25, 0.35, 0.35}; //actually, I'm not sure what's going on here with these cuts.  need thorough investigation
  std::vector<float> thresholdsReco_HE = {0.0, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};

  auto hbheRHT = std::make_unique<HBHERecHitTruthCollection>();

  int h=0;

  
  std::map<int, int> splitting_map; //keys are assigned parent IDs, values are the parent particle that the rechit was originally linked to
  std::map<int, double> parenteta; //keep track of energy weighted average eta of rechits that have a given truth label
  std::map<int, double> parentphi; //keep track of energy weighted average phi of rechits that have a given truth label

  for (const auto& erh : *rechits) { //loop over rechits.  This is where we will create the rechittruth objects
    
    //boolean for whether or not rechit passes energy thresholds
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

    /*
    //for each rechit, you can use the below to print out the simhits that contribute to this rechit and their energies
    if(simHit_energy_map.count(pId)){
      for( unsigned int i = 0; i < simHit_energy_map[pId].size(); i++ ){
	std::cout<<"   contribution from "<<(lSimHits)[simhitIndex_map[pId][i]].geantTrackId()<<" simHitEnergyMap: "<<simHit_energy_map[pId][i]<<std::endl;
      }
    }
    */




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
      PFcluster0Id   = PFClusterMatches[0].first;
      PFcluster0frac = PFClusterMatches[0].second;
    }
    if (PFClusterMatches.size() > 1){
      PFcluster1Id   = PFClusterMatches[1].first;
      PFcluster1frac = PFClusterMatches[1].second;
    } 
    if (PFClusterMatches.size() > 2){
      PFcluster2Id   = PFClusterMatches[2].first;
      PFcluster2frac = PFClusterMatches[2].second;
    }

    if(!detId_map.count(erh.detid().rawId())){ //no simHit for recHit.  Treat this as a noise hit.  Set truth label to -1 and create HBHERecHitTruth object
      
      HBHERecHitTruth rht(erh, -1, thisCell->getPosition().basicVector().x(), thisCell->getPosition().basicVector().y(), thisCell->getPosition().basicVector().z(), thisCell->repPos().eta(), thisCell->repPos().phi(), detid.depth(), layer, cleanHLT, cleanReco, PFcluster0Id, PFcluster0frac, PFcluster1Id, PFcluster1frac, PFcluster2Id, PFcluster2frac);
      hbheRHT->push_back(rht);

      h++;
      continue; //move on to the next rechit
    }

    //AT LEAST ONE SIM PARTICLE CONTRIBUTED TO THIS CLUSTER
    std::vector<int> parentPart; //parentPart will keep track of candidate truth labels

    for( unsigned int i = 0; i < detId_map[erh.detid().rawId()].size(); i++ ){
      unsigned int steps = 0;

      unsigned int vertValue = hSimTracks->at(track_map[detId_map[erh.detid().rawId()][i]]).vertIndex(); //for hbhe rec hits
      int tmpparentPart = detId_map[erh.detid().rawId()][i];
      if(tmpparentPart == -1){
	tmpparentPart = detId_map[erh.detid().rawId()][i];
      }
      while(steps < 100 && tmpparentPart >= 0){ //trace back potentially 100 steps in the parent-daughter particle shower.  Shower is never that deep
	vertValue = hSimTracks->at(track_map[tmpparentPart]).vertIndex();
	bool ecp = (esd == HcalEndcap);
	bool brl = (esd == HcalBarrel);
	if(hSimVertices->at(vertValue).parentIndex() < 0){
	  break;
	}
	//only trace back to the inner radius of the HCAL barrel or endcap's face
	if(esd == HcalBarrel && sqrt(hSimVertices->at(vertValue).position().X()*hSimVertices->at(vertValue).position().X() + hSimVertices->at(vertValue).position().Y()*hSimVertices->at(vertValue).position().Y()) < 180 ){
	  break;
	}
	if(esd == HcalEndcap && abs(hSimVertices->at(vertValue).position().Z()) < 390 ){
	  break;
	}
	tmpparentPart = hSimVertices->at(vertValue).parentIndex();
	steps++;
      }
      parentPart.push_back(tmpparentPart);
    }
    
    
    //Check if this rechit is within deltaR of 1 of existing average position for rechits with this truth-label, or within delR of other break-offs that had that truth label
    if(parentPart.size()>0){

      bool checkingRH = true;
      while(parenteta.count(parentPart.at(0)) && checkingRH){

	//quick delR calculation
	double deleta = (parenteta[parentPart.at(0)] - thisCell->repPos().eta());
        double delphi, ogphi, newphi;
	newphi = thisCell->repPos().phi();
	ogphi = parentphi[parentPart.at(0)];
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
	double deltaR = sqrt( deleta*deleta + delphi*delphi );

	if(deltaR > 1.){
	  if(splitting_map.count(parentPart.at(0))){
	    parentPart.at(0) = splitting_map[parentPart.at(0)];
	  } else{
	    largestTrackID++;
	    splitting_map[parentPart.at(0)] = largestTrackID; //make sure that the "breakoff" cluster has a unique label
	    parentPart.at(0) = largestTrackID;
	  }
	} else{
	  checkingRH = false;
	}

      }
  
      if(!parenteta.count(parentPart.at(0))){
	parenteta[parentPart.at(0)] = thisCell->repPos().eta();
	parentphi[parentPart.at(0)] = thisCell->repPos().phi();
      }
      
    }


    float relevant_energy = simHit_energy_map[erh.detid().rawId()][0];
    int relevant_depth = simHit_depth_map[erh.detid().rawId()][0];

    //Finally, create the labelled HBHERecHitTruth for this rechit :)
    HBHERecHitTruth rht(erh, parentPart.at(0), thisCell->getPosition().basicVector().x(), thisCell->getPosition().basicVector().y(), thisCell->getPosition().basicVector().z(), thisCell->repPos().eta(), thisCell->repPos().phi(), detid.depth(), layer, cleanHLT, cleanReco, PFcluster0Id, PFcluster0frac, PFcluster1Id, PFcluster1frac, PFcluster2Id, PFcluster2frac, relevant_energy, relevant_depth, 0.);

    hbheRHT->push_back(rht);

    h++;
  } // end loop over rec hits in the event

  iEvent.put( std::move(hbheRHT), "HBHERecHitTruth" ); //Add collection of HBHERecHitTruth objects to the event so that we can make an nTuple for training models later on

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void ProducerTestWithRadialCutoffEXPANDEDINFO::beginStream(edm::StreamID) {
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void ProducerTestWithRadialCutoffEXPANDEDINFO::endStream() {
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ProducerTestWithRadialCutoffEXPANDEDINFO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
  desc.add<edm::InputTag>("pf_src");
}

//define this as a plug-in
DEFINE_FWK_MODULE(ProducerTestWithRadialCutoffEXPANDEDINFO);
