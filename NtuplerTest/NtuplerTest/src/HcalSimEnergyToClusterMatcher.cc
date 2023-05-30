#include "NtuplerTest/NtuplerTest/interface/HcalSimEnergyToClusterMatcher.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalDigi/interface/HFDataFrame.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHitTruth.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/HcalTPGCompressor.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "CondFormats/HcalObjects/interface/HcalLutMetadata.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoParticleFlow/PFProducer/interface/PFAlgo.h"

#include <algorithm>

using namespace std;
using namespace edm;
using namespace boost; 

HcalSimEnergyToClusterMatcher::HcalSimEnergyToClusterMatcher(const edm::ParameterSet& iPS)  
  {
    //fSHitToken = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));
    //fSHitTokenEE = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","EcalHitsEE"));
    //fSHitTokenEB = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","EcalHitsEB"));
    //fSimTrackToken = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits",""));
    //fSimVertexToken = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits",""));
    
    //RHTToken = consumes<std::vector<HBHERecHitTruth>>(edm::InputTag("printTheDetIDs","HBHERecHitTruth"));
    RHTToken = consumes<std::vector<HBHERecHitTruth>>(iPS.getParameter<edm::InputTag>("clus_src"));
    
    fSHitToken = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));
    SimTrackToken = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits",""));
    SimVertexToken = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits",""));
    GenParticleToken = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles",""));

    _clustersLabel = consumes<reco::PFClusterCollection>(iPS.getParameter<edm::InputTag>("regular_clus_src"));

    hdrcToken_ = esConsumes<HcalDDDRecConstants, HcalRecNumberingRecord>();

    geomToken_ = esConsumes();
    

    fFile              = new TFile("Output_SimToClusterMatcher.root","RECREATE");

    fTree        = new TTree("Events","Events");

    fTree->Branch("nSimPart", &value_simPart_n, "nSimPart/i");
    fTree->Branch("simPart_energy", value_simPart_energy, "simPart_energy[nSimPart]/F");
    fTree->Branch("simPart_eta", value_simPart_eta, "simPart_eta[nSimPart]/F");
    fTree->Branch("simPart_PID", value_simPart_PID, "simPart_PID[nSimPart]/I");
    fTree->Branch("simPart_clusterEnergy", value_simPart_clusterEnergy, "simPart_clusterEnergy[nSimPart]/F");


}
void HcalSimEnergyToClusterMatcher::endJob() {
  fFile->cd();
  fTree->Write();
  fFile->Write();
  fFile->Close();
}
void HcalSimEnergyToClusterMatcher::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::cout<<"IN Sim-to-cluster matching NTUPLER PRODUCE"<<std::endl;

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




  edm::PCaloHitContainer lSimHits  = *hSimHits;

  std::map<unsigned int, std::vector<int>> detId_map; //maps keys = detectorID and values = Sim Track geantID
  std::map<unsigned int, std::vector<int>> simhitIndex_map; //maps keys = detectorID and values = simHit index that already corresponds to that det ID
  std::map<unsigned int, std::vector<float>> simHit_energy_map; //maps keys = detectorID and values = Sim Track energy
  std::map<unsigned int, std::vector<int>> simHit_depth_map; //maps keys = detectorID and values = Sim Track depth
  std::map<unsigned int, std::vector<float>> simHit_time_map; //maps keys = detectorID and values = Sim Track time
  std::map<int, float> trackID_energy_map; //map keys = track ID and values = total simhit energy for that particle
  int largestTrackID = 0;

  struct simEnergyStruct{
    int ID;
    float energy;
  };
  std::map<int, simEnergyStruct> simID_to_energy_map; //maps keys = Sim Track geantID and values = total simmed energy for this geantID


  fResponse = new CaloHitResponse(NULL, (CaloShapes*)NULL);
  fResponse->setGeometry(&*geoHandle);

  for (int j=0; j < (int) lSimHits.size(); j++) {

    HcalDetId simId = HcalHitRelabeller::relabel((lSimHits)[j].id(), fRecNumber);

    HcalSubdetector esd = (HcalSubdetector)simId.subdetId();
    std::shared_ptr<const CaloCellGeometry> thisCell = nullptr;
    switch (esd) {
    case HcalBarrel:
      thisCell = hcalBarrelGeo->getGeometry(simId);
      break;
    case HcalEndcap:
      thisCell = hcalEndcapGeo->getGeometry(simId);
      break;
    default:
      break;
    }
    if(!thisCell){
      continue;
    }


    double samplingFactor = 1;
    //HcalDetId simId = HcalHitRelabeller::relabel((lSimHits)[j].id(), fRecNumber);
    if(simId.subdet() == HcalBarrel) {
      samplingFactor = fSimParameterMap.hbParameters().samplingFactor(simId);
    } else if (simId.subdet() == HcalEndcap) {
      samplingFactor = fSimParameterMap.heParameters().samplingFactor(simId);
    }
    HcalDetId simIdForDepth = (lSimHits)[j].id();
    
    trackID_energy_map[(lSimHits)[j].geantTrackId()] += samplingFactor * (lSimHits)[j].energy();


    if(!simID_to_energy_map.count((lSimHits)[j].geantTrackId())){
      simEnergyStruct newstruct;
      newstruct.ID = (lSimHits)[j].geantTrackId();
      newstruct.energy = samplingFactor*((lSimHits)[j].energy());
      simID_to_energy_map[(lSimHits)[j].geantTrackId()] = newstruct;
    }
    else{
      simID_to_energy_map[(lSimHits)[j].geantTrackId()].energy += samplingFactor*((lSimHits)[j].energy());
    }


    if(!detId_map.count(simId.rawId())){ //new detID
      detId_map[simId.rawId()].push_back((lSimHits)[j].geantTrackId());
      simhitIndex_map[simId.rawId()].push_back(j);
      simHit_energy_map[simId.rawId()].push_back(samplingFactor*((lSimHits)[j].energy()));
      simHit_depth_map[simId.rawId()].push_back(simIdForDepth.depth());

      if((lSimHits)[j].geantTrackId() > largestTrackID) largestTrackID = (lSimHits)[j].geantTrackId();

    }
    else{
      detId_map[simId.rawId()].push_back((lSimHits)[j].geantTrackId());
      simhitIndex_map[simId.rawId()].push_back(j);
      simHit_energy_map[simId.rawId()].push_back(samplingFactor*((lSimHits)[j].energy()));
      simHit_depth_map[simId.rawId()].push_back(simIdForDepth.depth());

      if((lSimHits)[j].geantTrackId() > largestTrackID) largestTrackID = (lSimHits)[j].geantTrackId();
    }
  }





  std::map<unsigned int, unsigned int> track_map; //map keys = Sim Track geantID and values = Sim Track index
  for (unsigned int i = 0; i < hSimTracks->size(); ++i) {
    track_map[hSimTracks->at(i).trackId()] = i;
  }

  

  /*edm::Handle<edm::SortedCollection<HBHERecHit> > rechits;
  iEvent.getByToken(_rechitsLabel, rechits);
  int sizeRH = 0;
  for (const auto& erh : *rechits) {
    if(erh.energy() < 0.8){
      continue;
    }
    sizeRH++;
    }*/

  edm::Handle<reco::PFClusterCollection> hPFCHCAL;
  iEvent.getByToken(_clustersLabel, hPFCHCAL);
  assert(hPFCHCAL.isValid());
  std::map<HcalDetId,std::vector<std::pair<int, float>>> RH_PFClusterAndFrac;
  std::map<int, std::vector<std::vector<float>>> simID_toClustEnEtaPhi;
  for(reco::PFClusterCollection::const_iterator itC = hPFCHCAL->begin(); itC!=hPFCHCAL->end(); itC++){

    std::cout<<"NEW CLUSTER"<<std::endl;
    std::cout<<"cluster energy = "<<itC->energy()<<", cluster depth = "<<itC->depth()<<", cluster eta, phi = "<<itC->positionREP().eta()<<", "<<itC->positionREP().phi()<<" has nrechits = "<<itC->recHitFractions().size()<<std::endl;
    std::map<int, float> cluster_geantid_energy_map;
    

    const std::vector<reco::PFRecHitFraction>& hitsandfractions = itC->recHitFractions();
    for(long unsigned int it=0; it<hitsandfractions.size(); it++){
      const auto & RH = *hitsandfractions[it].recHitRef();

      HcalDetId pDetId(RH.detId());
      if (pDetId.det() != DetId::Hcal) continue;
      if(simHit_energy_map.count(pDetId.rawId())){
	std::cout<<"rechit "<<it<<" has energy "<<RH.energy()<<" but frac "<<hitsandfractions[it].fraction()<<". there are "<<simHit_energy_map[pDetId.rawId()].size()<<" simhits for this rechit at depth "<<RH.depth()<<" ("<<RH.positionREP().eta()<<", "<<RH.positionREP().phi()<<")"<<std::endl;
	
	for( unsigned int i = 0; i < simHit_energy_map[pDetId.rawId()].size(); i++ ){
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
      std::cout<<"    "<<(*r).first<<"    "<<(*r).second<<" which is "<<(*r).second/itC->energy()<<"\n";

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
    
  }

  std::vector<simEnergyStruct> simEnergyStructVector;
  std::map<int,simEnergyStruct>::iterator r;
  for (r = simID_to_energy_map.begin(); r!= simID_to_energy_map.end(); r++){
    simEnergyStructVector.push_back((*r).second);
  }

  std::sort(simEnergyStructVector.begin(), simEnergyStructVector.end(), [](const simEnergyStruct & a, const simEnergyStruct & b) -> bool{ return a.energy > b.energy; });


  std::vector<bool> clusterMask(hPFCHCAL->size(), true);

  value_simPart_n = 0;

  std::cout<<"LOOP OVER SORTED SIM IDS"<<std::endl;
  for(unsigned int i = 0; i < simEnergyStructVector.size(); i++){
    std::cout<<simEnergyStructVector.at(i).ID<<" "<<simEnergyStructVector.at(i).energy<<std::endl;
    std::vector<std::vector<float>> dummy = simID_toClustEnEtaPhi[simEnergyStructVector.at(i).ID];
    float bestMatchFrac = 0.;
    int bestMatchClusterIndex = -1;
    int bestMatchClusterInternalIndex = -1;
    for(unsigned int cl = 0; cl<dummy.size(); cl++){
      //std::cout<<"    clust "<<dummy.at(cl).at(0)<<" en "<<dummy.at(cl).at(1)<<" eta "<<dummy.at(cl).at(2)<<" phi "<<dummy.at(cl).at(3)<<" sim en in clust = "<<dummy.at(cl).at(4)<<" which is "<<dummy.at(cl).at(4)/trackID_energy_map[(*r).first]<<"\n";
      std::cout<<"    clust "<<dummy.at(cl).at(0)<<" en "<<dummy.at(cl).at(1)<<" eta "<<dummy.at(cl).at(2)<<" phi "<<dummy.at(cl).at(3)<<" sim en in clust = "<<dummy.at(cl).at(4)<<" which is "<<dummy.at(cl).at(4)/simEnergyStructVector.at(i).energy<<"\n";
      //if(dummy.at(cl).at(4)/trackID_energy_map[(*r).first] > bestMatchFrac && clusterMask[dummy.at(cl).at(0)]){
      if(dummy.at(cl).at(4)/simEnergyStructVector.at(i).energy > bestMatchFrac && clusterMask[dummy.at(cl).at(0)]){
	//bestMatchFrac = dummy.at(cl).at(4)/trackID_energy_map[(*r).first];
	bestMatchFrac = dummy.at(cl).at(4)/simEnergyStructVector.at(i).energy;
	bestMatchClusterIndex = dummy.at(cl).at(0);
	bestMatchClusterInternalIndex = cl;
      }
    }

    

    value_simPart_energy[value_simPart_n] = simEnergyStructVector.at(i).energy;
    value_simPart_eta[value_simPart_n] = hSimTracks->at(track_map[simEnergyStructVector.at(i).ID]).momentum().Eta();
    value_simPart_PID[value_simPart_n] = hSimTracks->at(track_map[simEnergyStructVector.at(i).ID]).type();
    if(bestMatchClusterInternalIndex >= 0){
      clusterMask[bestMatchClusterIndex] = 0;
      value_simPart_clusterEnergy[value_simPart_n] = dummy.at(bestMatchClusterInternalIndex).at(1);
    } else{
      value_simPart_clusterEnergy[value_simPart_n] = 0.;
    }
    value_simPart_ID[value_simPart_n] = simEnergyStructVector.at(i).ID;
    value_simPart_n++;
    

  }

  /*
  std::cout<<"LOOP OVER SIM IDS"<<std::endl;
  std::map<int,std::vector<std::vector<float>>>::iterator r;
  for (r = simID_toClustEnEtaPhi.begin(); r!= simID_toClustEnEtaPhi.end(); r++){
    std::cout<<(*r).first<<std::endl;
    for(unsigned int cl = 0; cl<(*r).second.size(); cl++){
      std::cout<<"    clust "<<(*r).second.at(cl).at(0)<<" en "<<(*r).second.at(cl).at(1)<<" eta "<<(*r).second.at(cl).at(2)<<" phi "<<(*r).second.at(cl).at(3)<<" sim en in clust = "<<(*r).second.at(cl).at(4)<<" which is "<<(*r).second.at(cl).at(4)/trackID_energy_map[(*r).first]<<"\n";
    }
  }
  */










  fTree->Fill();
  // Step C.1: Run FE Format Error / ZS for real data.
  //iEvent.put(std::move(result));
}

//DEFINE_FWK_MODULE(HcalSimEnergyToClusterMatcher);
