#include "NtuplerTest/NtuplerTest/interface/HcalTrigPrimDigiNtuplerEXPANDEDINFO.h"
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

HcalTrigPrimDigiNtuplerEXPANDEDINFO::HcalTrigPrimDigiNtuplerEXPANDEDINFO(const edm::ParameterSet& iPS)  
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

    hdrcToken_ = esConsumes<HcalDDDRecConstants, HcalRecNumberingRecord>();

    geomToken_ = esConsumes();
    

  fFile              = new TFile("Output_old.root","RECREATE");
  //fECALPFClusterArr  = new TClonesArray("baconhep::TECALPFCluster",200);
  //fEnergyArr = new TClonesArray("int",20000);
  //fRHParArr          = new TClonesArray("baconhep::TRHPart",20000);
  //fGenJetArr         = new TClonesArray("baconhep::TGenJet",20);
  //fGenFatJetArr      = new TClonesArray("baconhep::TGenJet",20);
  //fGenEventInfoArr   = new TClonesArray("baconhep::TGenEventInfo",20);
  //fEvtInfo         = new baconhep::TEventInfo();
  //fGenEvtInfo        = new baconhep::TGenEventInfo();
  //fGenParticleArr    = new TClonesArray("baconhep::TGenParticle",500);
  //fTrackArr    = new TClonesArray("baconhep::TSimTrack",40000);
  //fVertexArr    = new TClonesArray("baconhep::TSimVertex",40000);
 
  /*fFillerRH = new FillerRH(iPS,consumesCollector());
  edm::ParameterSet cfg(iPS.getUntrackedParameter<edm::ParameterSet>("GenJet"));
  fFillerGenJets = new baconhep::FillerGenJets(cfg,consumesCollector());
  edm::ParameterSet cfg1(iPS.getUntrackedParameter<edm::ParameterSet>("GenInfo"));
  fFillerGenInfo = new baconhep::FillerGenInfo(cfg1,consumesCollector());
  edm::ParameterSet cfg2(iPS.getUntrackedParameter<edm::ParameterSet>("Info"));
  fFillerEventInfo = new baconhep::FillerEventInfo(cfg2, 0, consumesCollector());
  //fFillerGenParticle = new baconhep::FillerGenParticle(iPS,consumesCollector());
  //fFillerSimTrack = new FillerSimTrack(iPS,consumesCollector());
  //fFillerSimVertex = new FillerSimVertex(iPS,consumesCollector());
  */
  fTree        = new TTree("Events","Events");

  fTree->Branch("nHit", &value_hit_n, "nHit/i");
  fTree->Branch("Hit_energy", value_hit_energy, "Hit_energy[nHit]/F");
  fTree->Branch("Hit_x", value_hit_x, "Hit_x[nHit]/F");
  fTree->Branch("Hit_y", value_hit_y, "Hit_y[nHit]/F");
  fTree->Branch("Hit_z", value_hit_z, "Hit_z[nHit]/F");
  fTree->Branch("Hit_eta", value_hit_eta, "Hit_eta[nHit]/F");
  fTree->Branch("Hit_phi", value_hit_phi, "Hit_phi[nHit]/F");
  fTree->Branch("Hit_depth", value_hit_depth, "Hit_depth[nHit]/I");
  fTree->Branch("Hit_isHB", value_hit_isHB, "Hit_isHB[nHit]/b");
  fTree->Branch("Hit_cleanHLT", value_hit_cleanHLT, "Hit_cleanHLT[nHit]/b");
  fTree->Branch("Hit_cleanReco", value_hit_cleanReco, "Hit_cleanReco[nHit]/b");
  fTree->Branch("Hit_time", value_hit_time, "Hit_time[nHit]/F");
  fTree->Branch("Hit_timeFalling", value_hit_timeFalling, "Hit_timeFalling[nHit]/F");
  fTree->Branch("Hit_PFcluster0Id", value_hit_PFcluster0Id, "Hit_PFcluster0Id[nHit]/I");
  fTree->Branch("Hit_PFcluster0frac", value_hit_PFcluster0frac, "Hit_PFcluster0frac[nHit]/F");
  fTree->Branch("Hit_PFcluster1Id", value_hit_PFcluster1Id, "Hit_PFcluster1Id[nHit]/I");
  fTree->Branch("Hit_PFcluster1frac", value_hit_PFcluster1frac, "Hit_PFcluster1frac[nHit]/F");
  fTree->Branch("Hit_PFcluster2Id", value_hit_PFcluster2Id, "Hit_PFcluster2Id[nHit]/I");
  fTree->Branch("Hit_PFcluster2frac", value_hit_PFcluster2frac, "Hit_PFcluster2frac[nHit]/F");
  fTree->Branch("Hit_parentPart", value_hit_parentPart, "Hit_parentPart[nHit]/I");
  fTree->Branch("SimHit_linkedEnergy", value_simHit_linkedEnergy, "SimHit_linkedEnergy[nHit]/F");
  //fTree->Branch("SimHit_depth", value_simHit_depth, "SimHit_depth[nHit]/I");
  //fTree->Branch("SimHit_time", value_simHit_time, "SimHit_time[nHit]/F");

  fTree->Branch("nSimHit", &value_simHit_n, "nSimHit/i");
  fTree->Branch("simHit_energy", value_simHit_energy, "simHit_energy[nSimHit]/F");
  fTree->Branch("simHit_x", value_simHit_x, "simHit_x[nSimHit]/F");
  fTree->Branch("simHit_y", value_simHit_y, "simHit_y[nSimHit]/F");
  fTree->Branch("simHit_z", value_simHit_z, "simHit_z[nSimHit]/F");
  fTree->Branch("simHit_eta", value_simHit_eta, "simHit_eta[nSimHit]/F");
  fTree->Branch("simHit_phi", value_simHit_phi, "simHit_phi[nSimHit]/F");
  fTree->Branch("simHit_depth", value_simHit_depth, "simHit_depth[nSimHit]/I");

  fTree->Branch("nSimHit_TRACEBACK", &value_simHit_TRACEBACK_n, "nSimHit_TRACEBACK/i");
  fTree->Branch("simHit_TRACEBACK_trackID", value_simHit_TRACEBACK_trackID, "simHit_TRACEBACK_trackID[nSimHit_TRACEBACK]/I");
  fTree->Branch("simHit_TRACEBACK_energy", value_simHit_TRACEBACK_energy, "simHit_TRACEBACK_energy[nSimHit_TRACEBACK]/F");
  fTree->Branch("simHit_TRACEBACK_x", value_simHit_TRACEBACK_x, "simHit_TRACEBACK_x[nSimHit_TRACEBACK]/F");
  fTree->Branch("simHit_TRACEBACK_y", value_simHit_TRACEBACK_y, "simHit_TRACEBACK_y[nSimHit_TRACEBACK]/F");
  fTree->Branch("simHit_TRACEBACK_z", value_simHit_TRACEBACK_z, "simHit_TRACEBACK_z[nSimHit_TRACEBACK]/F");
  fTree->Branch("simHit_TRACEBACK_eta", value_simHit_TRACEBACK_eta, "simHit_TRACEBACK_eta[nSimHit_TRACEBACK]/F");
  fTree->Branch("simHit_TRACEBACK_phi", value_simHit_TRACEBACK_phi, "simHit_TRACEBACK_phi[nSimHit_TRACEBACK]/F");
  fTree->Branch("simHit_TRACEBACK_depth", value_simHit_TRACEBACK_depth, "simHit_TRACEBACK_depth[nSimHit_TRACEBACK]/I");

  fTree->Branch("nSimTrack", &value_simTrack_n, "nSimTrack/i");
  fTree->Branch("simTrack_energy", value_simTrack_energy, "simTrack_energy[nSimTrack]/F");
  fTree->Branch("simTrack_x", value_simTrack_x, "simTrack_x[nSimTrack]/F");
  fTree->Branch("simTrack_y", value_simTrack_y, "simTrack_y[nSimTrack]/F");
  fTree->Branch("simTrack_z", value_simTrack_z, "simTrack_z[nSimTrack]/F");
  fTree->Branch("simTrack_eta", value_simTrack_eta, "simTrack_eta[nSimTrack]/F");
  fTree->Branch("simTrack_phi", value_simTrack_phi, "simTrack_phi[nSimTrack]/F");
  fTree->Branch("simTrack_depth", value_simTrack_depth, "simTrack_depth[nSimTrack]/I");

  fTree->Branch("nSimTrack_TRACEBACK", &value_simTrack_TRACEBACK_n, "nSimTrack_TRACEBACK/i");
  fTree->Branch("simTrack_TRACEBACK_ID", value_simTrack_TRACEBACK_ID, "simTrack_TRACEBACK_ID[nSimTrack_TRACEBACK]/I");
  fTree->Branch("simTrack_TRACEBACK_genEnergy", value_simTrack_TRACEBACK_genEnergy, "simTrack_TRACEBACK_genEnergy[nSimTrack_TRACEBACK]/F");
  fTree->Branch("simTrack_TRACEBACK_genPDGID", value_simTrack_TRACEBACK_genPDGID, "simTrack_TRACEBACK_genPDGID[nSimTrack_TRACEBACK]/I");
  fTree->Branch("simTrack_TRACEBACK_simmedEnergy", value_simTrack_TRACEBACK_simmedEnergy, "simTrack_TRACEBACK_simmedEnergy[nSimTrack_TRACEBACK]/F");
  fTree->Branch("simTrack_TRACEBACK_x", value_simTrack_TRACEBACK_x, "simTrack_TRACEBACK_x[nSimTrack_TRACEBACK]/F");
  fTree->Branch("simTrack_TRACEBACK_y", value_simTrack_TRACEBACK_y, "simTrack_TRACEBACK_y[nSimTrack_TRACEBACK]/F");
  fTree->Branch("simTrack_TRACEBACK_z", value_simTrack_TRACEBACK_z, "simTrack_TRACEBACK_z[nSimTrack_TRACEBACK]/F");
  fTree->Branch("simTrack_TRACEBACK_genEta", value_simTrack_TRACEBACK_genEta, "simTrack_TRACEBACK_genEta[nSimTrack_TRACEBACK]/F");
  fTree->Branch("simTrack_TRACEBACK_genPhi", value_simTrack_TRACEBACK_genPhi, "simTrack_TRACEBACK_genPhi[nSimTrack_TRACEBACK]/F");
  fTree->Branch("simTrack_TRACEBACK_depth", value_simTrack_TRACEBACK_depth, "simTrack_TRACEBACK_depth[nSimTrack_TRACEBACK]/I");
  fTree->Branch("simTrack_TRACEBACK_hasHit_withDepth", value_simTrack_TRACEBACK_hasHit_withDepth, "simTrack_TRACEBACK_hasHit_withDepth[nSimTrack_TRACEBACK][7]/I");
  fTree->Branch("simTrack_TRACEBACK_energy_withDepth", value_simTrack_TRACEBACK_energy_withDepth, "simTrack_TRACEBACK_energy_withDepth[nSimTrack_TRACEBACK][7]/F");
  fTree->Branch("simTrack_TRACEBACK_avgPhi_withDepth", value_simTrack_TRACEBACK_avgPhi_withDepth, "simTrack_TRACEBACK_avgPhi_withDepth[nSimTrack_TRACEBACK][7]/F");
  fTree->Branch("simTrack_TRACEBACK_avgEta_withDepth", value_simTrack_TRACEBACK_avgEta_withDepth, "simTrack_TRACEBACK_avgEta_withDepth[nSimTrack_TRACEBACK][7]/F");
  fTree->Branch("simTrack_TRACEBACK_widthPhi_withDepth", value_simTrack_TRACEBACK_widthPhi_withDepth, "simTrack_TRACEBACK_widthPhi_withDepth[nSimTrack_TRACEBACK][7]/F");
  fTree->Branch("simTrack_TRACEBACK_widthEta_withDepth", value_simTrack_TRACEBACK_widthEta_withDepth, "simTrack_TRACEBACK_widthEta_withDepth[nSimTrack_TRACEBACK][7]/F");
  fTree->Branch("simTrack_TRACEBACK_EWavgPhi_withDepth", value_simTrack_TRACEBACK_EWavgPhi_withDepth, "simTrack_TRACEBACK_EWavgPhi_withDepth[nSimTrack_TRACEBACK][7]/F");
  fTree->Branch("simTrack_TRACEBACK_EWavgEta_withDepth", value_simTrack_TRACEBACK_EWavgEta_withDepth, "simTrack_TRACEBACK_EWavgEta_withDepth[nSimTrack_TRACEBACK][7]/F");
  fTree->Branch("simTrack_TRACEBACK_EWwidthPhi_withDepth", value_simTrack_TRACEBACK_EWwidthPhi_withDepth, "simTrack_TRACEBACK_EWwidthPhi_withDepth[nSimTrack_TRACEBACK][7]/F");
  fTree->Branch("simTrack_TRACEBACK_EWwidthEta_withDepth", value_simTrack_TRACEBACK_EWwidthEta_withDepth, "simTrack_TRACEBACK_EWwidthEta_withDepth[nSimTrack_TRACEBACK][7]/F");

}
void HcalTrigPrimDigiNtuplerEXPANDEDINFO::endJob() {
  fFile->cd();
  fTree->Write();
  fFile->Write();
  fFile->Close();
}
void HcalTrigPrimDigiNtuplerEXPANDEDINFO::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::cout<<"IN NTUPLER PRODUCE"<<std::endl;

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


  int largestTrackID = 0;

  fResponse = new CaloHitResponse(NULL, (CaloShapes*)NULL);
  fResponse->setGeometry(&*geoHandle);

  //std::map<unsigned int, int> detId_map; //maps keys = detectorID and values = Sim Track geantID
  //std::map<unsigned int, int> simhitIndex_map; //maps keys = detectorID and values = simHit index that already corresponds to that det ID
  //std::map<unsigned int, float> simHit_energy_map; //maps keys = detectorID and values = Sim Track energy
  //std::map<unsigned int, int> simHit_depth_map; //maps keys = detectorID and values = Sim Track depth
  //std::map<unsigned int, float> simHit_time_map; //maps keys = detectorID and values = Sim Track time

  std::map<unsigned int, unsigned int> track_map; //map keys = Sim Track geantID and values = Sim Track index
  for (unsigned int i = 0; i < hSimTracks->size(); ++i) {
    track_map[hSimTracks->at(i).trackId()] = i;
    //if(trackID_energy_map.count(hSimTracks->at(i).trackId())){
    //  std::cout<<"total simhit energy for this trackid = "<<trackID_energy_map[hSimTracks->at(i).trackId()]<<" full simtrack energy was "<<hSimTracks->at(i).momentum().E()<<std::endl;
    //}
  }

  std::map<std::pair<int, HcalDetId>, float> simhit_energy_map; //key is pair of geant track id and HcalDetId.  Value is energy (have to sum it up)
  std::map<std::pair<int, HcalDetId>, float> simhit_energy_map_TRACEBACK; //key is pair of geant track id and HcalDetId.  Value is energy (have to sum it up)
  std::vector<int> uniqueTrackIDs;
  std::vector<int> uniqueTrackIDs_TRACEBACK;
  std::map<int, float> trackID_energy_map; //map keys = track ID and values = total simhit energy for that particle
  std::map<int, float> trackID_energy_map_TRACEBACK; //map keys = track ID and values = total simhit energy for that particle
  std::map<int, int> full_tracebacks; //keys = trackIDs, values = tracedback trackIDs

  value_simHit_TRACEBACK_n = 0;

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

    /*
      std::shared_ptr<const CaloCellGeometry> thisCell = nullptr;
    std::shared_ptr<const CaloCellGeometry> thisCell = nullptr;
    PFLayer::Layer layer = PFLayer::HCAL_BARREL1;
    switch (esd) {
    case HcalBarrel:
      thisCell = hcalBarrelGeo->getGeometry(detid);
      layer = PFLayer::HCAL_BARREL1;
      break;

    case HcalEndcap:
      thisCell = hcalEndcapGeo->getGeometry(detid);
      layer = PFLayer::HCAL_ENDCAP;
      break;
    default:
      break;
    }

    // find rechit geometry
    if (!thisCell) {
      edm::LogError("HcalTrigPrimDigiNtuplerEXPANDEDINFO") << "warning detid " << std::hex << detid.rawId() << std::dec << " "
                                           << detid << " not found in geometry" << std::endl;
      continue;
    }
    */

    if(!full_tracebacks.count((lSimHits)[j].geantTrackId())){
      unsigned int steps = 0;
      unsigned int vertValue = hSimTracks->at(track_map[(lSimHits)[j].geantTrackId()]).vertIndex();
      int tmpparentPart = (lSimHits)[j].geantTrackId();
      while(steps < 100 && tmpparentPart > 0){
	vertValue = hSimTracks->at(track_map[tmpparentPart]).vertIndex();
	bool ecp = (esd == HcalEndcap);
	bool brl = (esd == HcalBarrel);
	if(hSimVertices->at(vertValue).parentIndex() <= 0){
	  break;
	}
	if(esd == HcalBarrel && sqrt(hSimVertices->at(vertValue).position().X()*hSimVertices->at(vertValue).position().X() + hSimVertices->at(vertValue).position().Y()*hSimVertices->at(vertValue).position().Y()) < 180 ){
	  break;
	}
	if(esd == HcalEndcap && abs(hSimVertices->at(vertValue).position().Z()) < 390 ){
	  break;
	}
	tmpparentPart = hSimVertices->at(vertValue).parentIndex();
	std::cout<<"simhit for trackid "<<(lSimHits)[j].geantTrackId()<<" tracedback to parent "<<tmpparentPart<<std::endl;
	steps++;
      }
      full_tracebacks[(lSimHits)[j].geantTrackId()] = tmpparentPart;
    }

    std::cout<<"j = "<<j<<" simID "<<simId.rawId()<<" geantTrackId() = "<<(lSimHits)[j].geantTrackId()<<" time = "<<(lSimHits)[j].time()<<" energy = "<<samplingFactor * (lSimHits)[j].energy()<<" detid = "<<(lSimHits)[j].id()<<" depth = "<<simId.depth()<<std::endl;

    trackID_energy_map[(lSimHits)[j].geantTrackId()] += samplingFactor * (lSimHits)[j].energy();
    trackID_energy_map_TRACEBACK[full_tracebacks[(lSimHits)[j].geantTrackId()]] += samplingFactor * (lSimHits)[j].energy();

    simhit_energy_map[ std::make_pair((lSimHits)[j].geantTrackId(), simId) ] += samplingFactor * (lSimHits)[j].energy();
    simhit_energy_map_TRACEBACK[ std::make_pair(full_tracebacks[(lSimHits)[j].geantTrackId()], simId) ] += samplingFactor * (lSimHits)[j].energy();

    if(!(std::count(uniqueTrackIDs.begin(), uniqueTrackIDs.end(), (lSimHits)[j].geantTrackId()))){
      uniqueTrackIDs.push_back((lSimHits)[j].geantTrackId());
    }
    if(!(std::count(uniqueTrackIDs_TRACEBACK.begin(), uniqueTrackIDs_TRACEBACK.end(), full_tracebacks[(lSimHits)[j].geantTrackId()]))){
      uniqueTrackIDs_TRACEBACK.push_back(full_tracebacks[(lSimHits)[j].geantTrackId()]);
    }
    /*
    value_simHit_TRACEBACK_trackID[j] = full_tracebacks[(lSimHits)[j].geantTrackId()];
    value_simHit_TRACEBACK_energy[j] = samplingFactor * (lSimHits)[j].energy();
    value_simHit_TRACEBACK_x[j] = thisCell->getPosition().basicVector().x();
    value_simHit_TRACEBACK_y[j] = thisCell->getPosition().basicVector().y();
    value_simHit_TRACEBACK_z[j] = thisCell->getPosition().basicVector().z();
    value_simHit_TRACEBACK_eta[j] = thisCell->repPos().eta();
    value_simHit_TRACEBACK_phi[j] = thisCell->repPos().phi();
    value_simHit_TRACEBACK_depth[j] = simId.depth();
    */
    value_simHit_TRACEBACK_trackID[value_simHit_TRACEBACK_n] = full_tracebacks[(lSimHits)[j].geantTrackId()];
    value_simHit_TRACEBACK_energy[value_simHit_TRACEBACK_n] = samplingFactor * (lSimHits)[j].energy();
    value_simHit_TRACEBACK_x[value_simHit_TRACEBACK_n] = thisCell->getPosition().basicVector().x();
    value_simHit_TRACEBACK_y[value_simHit_TRACEBACK_n] = thisCell->getPosition().basicVector().y();
    value_simHit_TRACEBACK_z[value_simHit_TRACEBACK_n] = thisCell->getPosition().basicVector().z();
    value_simHit_TRACEBACK_eta[value_simHit_TRACEBACK_n] = thisCell->repPos().eta();
    value_simHit_TRACEBACK_phi[value_simHit_TRACEBACK_n] = thisCell->repPos().phi();
    value_simHit_TRACEBACK_depth[value_simHit_TRACEBACK_n] = simId.depth();

    value_simHit_TRACEBACK_n++;
    std::cout<<"value_simHit_TRACEBACK_n = "<<value_simHit_TRACEBACK_n<<std::endl;

  }

  std::cout<<"There were "<<uniqueTrackIDs.size()<<" unique trackIDs with simhits"<<std::endl;
  std::cout<<"There were "<<uniqueTrackIDs_TRACEBACK.size()<<" unique trackIDs with simhits after tracing back"<<std::endl;
  std::cout<<"There were "<<simhit_energy_map_TRACEBACK.size()<<" unique trackID+cell combos"<<std::endl;
  
  value_simTrack_n = uniqueTrackIDs.size();
  value_simTrack_TRACEBACK_n = uniqueTrackIDs_TRACEBACK.size();

  /*
  for (auto const& [key, val] : symbolTable)
    {
      std::cout << key        // string (key)
		<< ':'  
		<< val        // string's value
		<< std::endl;
    }
  */

  int totalskipped = 0;
  for (auto const& it : simhit_energy_map_TRACEBACK){
    HcalSubdetector esd = (HcalSubdetector)it.first.second.subdetId();
    std::shared_ptr<const CaloCellGeometry> thisCell = nullptr;
    switch (esd) {
    case HcalBarrel:
      thisCell = hcalBarrelGeo->getGeometry(it.first.second);
      break;
    case HcalEndcap:
      thisCell = hcalEndcapGeo->getGeometry(it.first.second);
      break;
    default:
      break;
    }
    if(!thisCell){
      totalskipped++;
    }
  }
  std::cout<<"total skipped = "<<totalskipped<<std::endl;
  
  
  std::map<int, std::vector<float>> trackID_average_phi_map_TRACEBACK; //map keys = track ID and values = vector of average phi of simhits for trackid (vector index is depth)
  std::map<int, std::vector<float>> trackID_average_eta_map_TRACEBACK; //map keys = track ID and values = vector of average phi of simhits for trackid (vector index is depth)
  std::map<int, std::vector<float>> trackID_width_phi_map_TRACEBACK; //map keys = track ID and values = vector of average phi of simhits for trackid (vector index is depth)
  std::map<int, std::vector<float>> trackID_width_eta_map_TRACEBACK; //map keys = track ID and values = vector of average phi of simhits for trackid (vector index is depth)
  std::map<int, std::vector<float>> trackID_EWaverage_phi_map_TRACEBACK; //map keys = track ID and values = vector of average phi of simhits for trackid (vector index is depth)
  std::map<int, std::vector<float>> trackID_EWaverage_eta_map_TRACEBACK; //map keys = track ID and values = vector of average phi of simhits for trackid (vector index is depth)
  std::map<int, std::vector<float>> trackID_EWwidth_phi_map_TRACEBACK; //map keys = track ID and values = vector of average phi of simhits for trackid (vector index is depth)
  std::map<int, std::vector<float>> trackID_EWwidth_eta_map_TRACEBACK; //map keys = track ID and values = vector of average phi of simhits for trackid (vector index is depth)
  std::map<int, std::vector<bool>> trackID_hasHit_map_TRACEBACK; //map keys = track ID and values = vector of average phi of simhits for trackid (vector index is depth)
  std::map<int, std::vector<float>> trackID_energyDepth_map_TRACEBACK; //map keys = track ID and values = vector of average phi of simhits for trackid (vector index is depth)
  for(unsigned int t=0; t<uniqueTrackIDs_TRACEBACK.size();t++){
    for(int d = 1; d < 8; d++){
      
      bool atLeast1Hit = false;
      float avgX = 0.;
      float avgY = 0.;
      float avgZ = 0.;
      float totalN = 0.;
      float avgEX = 0.;
      float avgEY = 0.;
      float avgEZ = 0.;
      float totalE = 0.;
      std::vector<float> energies;
      std::vector<float> normedXs;
      std::vector<float> normedYs;
      std::vector<float> normedZs;
      std::vector<float> etas;
      std::vector<float> phis;
      
      
      for (auto const& it : simhit_energy_map_TRACEBACK){
	if(it.first.first == uniqueTrackIDs_TRACEBACK.at(t) && it.first.second.depth() == d){
	  
	  atLeast1Hit = true;
	  
	  HcalSubdetector esd = (HcalSubdetector)it.first.second.subdetId();
	  std::shared_ptr<const CaloCellGeometry> thisCell = nullptr;
	  switch (esd) {
	  case HcalBarrel:
	    thisCell = hcalBarrelGeo->getGeometry(it.first.second);
	    break;
	    
	  case HcalEndcap:
	    thisCell = hcalEndcapGeo->getGeometry(it.first.second);
	    break;
	  default:
	    break;
	  }
	  
	  // find rechit geometry
	  if (!thisCell) {
	    continue;
	  }
	  
	  float normedx = thisCell->getPosition().basicVector().x()/thisCell->getPosition().basicVector().mag();
	  float normedy = thisCell->getPosition().basicVector().y()/thisCell->getPosition().basicVector().mag();
	  float normedz = thisCell->getPosition().basicVector().z()/thisCell->getPosition().basicVector().mag();
	  float yoverx = normedy/normedx;
	  float phi1 = atan2( normedy, normedx );
	  float cylindricalR = sqrt(normedx*normedx + normedy*normedy);
	  float cylindricalTheta = atan2(cylindricalR, normedz);
	  float computeEta = -1.*log(tan(cylindricalTheta/2.));
	  
	  std::cout<<"simhit for trackid "<<it.first.first<<" in depth "<<d<<" pos = ("<<thisCell->getPosition().basicVector().x()<<", "<<thisCell->getPosition().basicVector().y()<<", "<<thisCell->getPosition().basicVector().z()<<") or (eta, phi) = ("<<thisCell->repPos().eta()<<", "<<thisCell->repPos().phi()<<") has energy = "<<it.second<<" normed vector? ("<<thisCell->getPosition().basicVector().x()/thisCell->getPosition().basicVector().mag()<<", "<<thisCell->getPosition().basicVector().y()/thisCell->getPosition().basicVector().mag()<<", "<<thisCell->getPosition().basicVector().z()/thisCell->getPosition().basicVector().mag()<<") and test derived phi = "<<phi1<<" and calculated eta = "<<computeEta<<std::endl;
	  std::cout<<"to be clear atan( yoverx ) = "<<atan2( normedy, normedx )<<std::endl;
	  
	  avgX += normedx;
	  avgY += normedy;
	  avgZ += normedz;
	  totalN += 1;
	  avgEX += it.second*normedx;
	  avgEY += it.second*normedy;
	  avgEZ += it.second*normedz;
	  totalE += it.second;
	  energies.push_back(it.second);
	  normedXs.push_back(normedx);
	  normedYs.push_back(normedy);
	  normedZs.push_back(normedz);
	  etas.push_back(computeEta);
	  phis.push_back(phi1);
	  
	}
      }
      
      if(atLeast1Hit){

	float normedAvgX = avgX/totalN;
	float normedAvgY = avgY/totalN;
	float normedAvgZ = avgZ/totalN;
	float avgPhi = atan2( normedAvgY, normedAvgX );
	float avgR = sqrt(normedAvgX*normedAvgX + normedAvgY*normedAvgY);
	float avgTheta = atan2(avgR, normedAvgZ);
	float avgEta = -1.*log(tan(avgTheta/2.));
	std::cout<<"for trackid "<<uniqueTrackIDs_TRACEBACK.at(t)<<" and depth "<<d<<" average x = "<<avgX/totalE<<" and average y = "<<avgY/totalE<<" so avg Phi = "<<avgPhi<<" and average eta = "<<avgEta<<std::endl;
	float etawidth = 0.;
	float phiwidth = 0.;
	for(unsigned int v = 0; v < energies.size(); v++){
	  etawidth += (etas.at(v) - avgEta)*(etas.at(v) - avgEta);
	  float phidiff = abs(phis.at(v) - avgPhi);
	  if(phidiff > M_PI){
	    phidiff = 2.*M_PI - phidiff;
	  }
	  phiwidth += phidiff*phidiff;
	}
	etawidth = etawidth/totalN;
	etawidth = sqrt(etawidth);
	phiwidth = phiwidth/totalN;
	phiwidth = sqrt(phiwidth);

	float normedAvgEX = avgEX/totalE;
	float normedAvgEY = avgEY/totalE;
	float normedAvgEZ = avgEZ/totalE;
	float avgEPhi = atan2( normedAvgEY, normedAvgEX );
	float avgER = sqrt(normedAvgEX*normedAvgEX + normedAvgEY*normedAvgEY);
	float avgETheta = atan2(avgER, normedAvgEZ);
	float avgEEta = -1.*log(tan(avgETheta/2.));
	std::cout<<"for trackid "<<uniqueTrackIDs_TRACEBACK.at(t)<<" and depth "<<d<<" average EW x = "<<avgEX/totalE<<" and average EW y = "<<avgEY/totalE<<" so avg EW Phi = "<<avgEPhi<<" and average EW eta = "<<avgEEta<<std::endl;
	float etaEwidth = 0.;
	float phiEwidth = 0.;
	for(unsigned int v = 0; v < energies.size(); v++){
	  etaEwidth += energies.at(v)*(etas.at(v) - avgEEta)*(etas.at(v) - avgEEta);
	  float phidiff = abs(phis.at(v) - avgEPhi);
	  if(phidiff > M_PI){
	    phidiff = 2.*M_PI - phidiff;
	  }
	  phiEwidth += energies.at(v)*phidiff*phidiff;
	}
	etaEwidth = etaEwidth/totalE;
	etaEwidth = sqrt(etaEwidth);
	phiEwidth = phiEwidth/totalE;
	phiEwidth = sqrt(phiEwidth);

	std::cout<<"also etawidth = "<<etawidth<<" and phiwidth = "<<phiwidth<<std::endl;
	trackID_hasHit_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(atLeast1Hit);
	trackID_average_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(avgPhi);
	trackID_average_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(avgEta);
	trackID_EWaverage_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(avgEPhi);
	trackID_EWaverage_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(avgEEta);
	trackID_width_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(phiwidth);
	trackID_width_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(etawidth);
	trackID_EWwidth_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(phiEwidth);
	trackID_EWwidth_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(etaEwidth);
	trackID_energyDepth_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(totalE);
      }
      else{      
	trackID_hasHit_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(0);
	trackID_average_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(-99);
	trackID_average_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(-99);
	trackID_EWaverage_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(-99);
	trackID_EWaverage_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(-99);
	trackID_width_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(-99);
	trackID_width_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(-99);
	trackID_EWwidth_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(-99);
	trackID_EWwidth_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(-99);
	trackID_energyDepth_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].push_back(-99);
      }
    }
  }




  for(unsigned int t=0; t<uniqueTrackIDs_TRACEBACK.size();t++){
    /*
    value_simTrack_TRACEBACK_simmedEnergy.push_back(trackID_energy_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)]);
    value_simTrack_TRACEBACK_genEnergy.push_back( hSimTracks->at(track_map[uniqueTrackIDs_TRACEBACK.at(t)]).momentum().E() );
    value_simTrack_TRACEBACK_genPDGID.push_back( hSimTracks->at(track_map[uniqueTrackIDs_TRACEBACK.at(t)]).type() );
    value_simTrack_TRACEBACK_genEta.push_back( hGenParticles->at(hSimTracks->at(track_map[uniqueTrackIDs_TRACEBACK.at(t)]).genpartIndex() - 1).eta() );
    value_simTrack_TRACEBACK_genPhi.push_back( hGenParticles->at(hSimTracks->at(track_map[uniqueTrackIDs_TRACEBACK.at(t)]).genpartIndex() - 1).phi() );
    */
    //std::cout<<"doing track id "<<uniqueTrackIDs_TRACEBACK.at(t)<<" genparticle index? "<<hSimTracks->at(track_map[uniqueTrackIDs_TRACEBACK.at(t)]).genpartIndex() - 1<<std::endl;
    value_simTrack_TRACEBACK_ID[t] = uniqueTrackIDs_TRACEBACK.at(t);
    value_simTrack_TRACEBACK_simmedEnergy[t] = trackID_energy_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)];
    value_simTrack_TRACEBACK_genEnergy[t] = hSimTracks->at(track_map[uniqueTrackIDs_TRACEBACK.at(t)]).momentum().E();
    value_simTrack_TRACEBACK_genPDGID[t] = hSimTracks->at(track_map[uniqueTrackIDs_TRACEBACK.at(t)]).type();
    int genparticleindex = hSimTracks->at(track_map[uniqueTrackIDs_TRACEBACK.at(t)]).genpartIndex() - 1;
    if(genparticleindex >= 0){
      value_simTrack_TRACEBACK_genEta[t] = hGenParticles->at(hSimTracks->at(track_map[uniqueTrackIDs_TRACEBACK.at(t)]).genpartIndex() - 1).eta();
      value_simTrack_TRACEBACK_genPhi[t] = hGenParticles->at(hSimTracks->at(track_map[uniqueTrackIDs_TRACEBACK.at(t)]).genpartIndex() - 1).phi();
    }
    else{
      value_simTrack_TRACEBACK_genEta[t] = -99;
      value_simTrack_TRACEBACK_genPhi[t] = -99;
    }
    for(unsigned int d = 0; d<7; d++){
      value_simTrack_TRACEBACK_hasHit_withDepth[t][d] = trackID_hasHit_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].at(d);
      value_simTrack_TRACEBACK_energy_withDepth[t][d] = trackID_energyDepth_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].at(d);
      value_simTrack_TRACEBACK_avgPhi_withDepth[t][d] = trackID_average_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].at(d);
      value_simTrack_TRACEBACK_avgEta_withDepth[t][d] = trackID_average_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].at(d);
      value_simTrack_TRACEBACK_widthPhi_withDepth[t][d] = trackID_width_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].at(d);
      value_simTrack_TRACEBACK_widthEta_withDepth[t][d] = trackID_width_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].at(d);
      value_simTrack_TRACEBACK_EWavgPhi_withDepth[t][d] = trackID_EWaverage_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].at(d);
      value_simTrack_TRACEBACK_EWavgEta_withDepth[t][d] = trackID_EWaverage_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].at(d);
      value_simTrack_TRACEBACK_EWwidthPhi_withDepth[t][d] = trackID_EWwidth_phi_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].at(d);
      value_simTrack_TRACEBACK_EWwidthEta_withDepth[t][d] = trackID_EWwidth_eta_map_TRACEBACK[uniqueTrackIDs_TRACEBACK.at(t)].at(d);
    } 
    
  }



  edm::Handle<std::vector<HBHERecHitTruth>> hRHT;
  iEvent.getByToken(RHTToken, hRHT);

  value_hit_n = 0;

  for(unsigned int r = 0; r < hRHT->size(); r++){
    value_hit_energy[r] = hRHT->at(r).energy();
    value_hit_x[r] = hRHT->at(r).x();
    value_hit_y[r] = hRHT->at(r).y();
    value_hit_z[r] = hRHT->at(r).z();
    value_hit_eta[r] = hRHT->at(r).eta();
    value_hit_phi[r] = hRHT->at(r).phi();
    value_hit_depth[r] = hRHT->at(r).depth();
    value_hit_isHB[r] = hRHT->at(r).isHB();
    value_hit_cleanHLT[r] = hRHT->at(r).isCleanHLT();
    value_hit_cleanReco[r] = hRHT->at(r).isCleanReco();
    value_hit_time[r] = hRHT->at(r).time();
    value_hit_timeFalling[r] = hRHT->at(r).timeFalling();
    value_hit_PFcluster0Id[r] = hRHT->at(r).PFcluster0Id();
    value_hit_PFcluster0frac[r] = hRHT->at(r).PFcluster0frac();
    value_hit_PFcluster1Id[r] = hRHT->at(r).PFcluster1Id();
    value_hit_PFcluster1frac[r] = hRHT->at(r).PFcluster1frac();
    value_hit_PFcluster2Id[r] = hRHT->at(r).PFcluster2Id();
    value_hit_PFcluster2frac[r] = hRHT->at(r).PFcluster2frac();
    value_hit_parentPart[r] = hRHT->at(r).parentPart();
    value_simHit_linkedEnergy[r] = hRHT->at(r).simHitEnergy();
    //value_simHit_depth[r] = hRHT->at(r).simHitDepth();
    //value_simHit_time[r] = hRHT->at(r).simHitTime();
    value_hit_n++;
  }


  fTree->Fill();
  // Step C.1: Run FE Format Error / ZS for real data.
  //iEvent.put(std::move(result));
}

//DEFINE_FWK_MODULE(HcalTrigPrimDigiNtuplerEXPANDEDINFO);
