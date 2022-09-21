#include "NtuplerTest/NtuplerTest/interface/HcalTrigPrimDigiNtupler.h"
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

HcalTrigPrimDigiNtupler::HcalTrigPrimDigiNtupler(const edm::ParameterSet& iPS)  
  {
    //fSHitToken = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));
  //fSHitTokenEE = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","EcalHitsEE"));
  //fSHitTokenEB = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","EcalHitsEB"));
  //fSimTrackToken = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits",""));
  //fSimVertexToken = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits",""));

  //RHTToken = consumes<std::vector<HBHERecHitTruth>>(edm::InputTag("printTheDetIDs","HBHERecHitTruth"));
    RHTToken = consumes<std::vector<HBHERecHitTruth>>(iPS.getParameter<edm::InputTag>("clus_src"));

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
  fTree->Branch("SimHit_energy", value_simHit_energy, "SimHit_energy[nHit]/F");
  fTree->Branch("SimHit_depth", value_simHit_depth, "SimHit_depth[nHit]/I");
  fTree->Branch("SimHit_time", value_simHit_time, "SimHit_time[nHit]/F");

  //fTree->Branch("RecHit", &fRHParArr); 
  //fTree->Branch("Energy", &fEnergyArr); 
  //fTree->Branch("GenParticle", &fGenParticleArr); 
  //fTree->Branch("GenEventInfo", &fGenEvtInfo); 
  //fTree->Branch("EventInfo", &fEvtInfo); 
  //fTree->Branch("GenJet", &fGenJetArr);
  //fTree->Branch("GenFatJet", &fGenFatJetArr); 
  //fTree->Branch("ECALPFCluster", &fECALPFClusterArr); 
  //fTree->Branch("SimTrack", &fTrackArr); 
  //fTree->Branch("SimVertex", &fVertexArr); 
  //fTree->Branch("PFDepth"   , &fPFParArr);
  //fTree->Branch("GenParticle", &fGenParArr); 
}
void HcalTrigPrimDigiNtupler::endJob() {
  fFile->cd();
  fTree->Write();
  fFile->Write();
  fFile->Close();
}
void HcalTrigPrimDigiNtupler::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //fRHParArr->Clear();

  //TClonesArray &rArray = *fRHParArr;

  edm::Handle<std::vector<HBHERecHitTruth>> hRHT;
  iEvent.getByToken(RHTToken, hRHT);
  //recHitHCAL = hRHT.product();

  //std::cout<<"new event "<<hRHT->size()<<std::endl;

  //baconhep::TRHPart *pRH = (baconhep::TRHPart*) rArray[index];

  value_hit_n = 0;

  for(unsigned int r = 0; r < hRHT->size(); r++){

    //assert(rArray.GetEntries() < rArray.GetSize());
    //const int index = rArray.GetEntries();
    //new(rArray[index]) baconhep::TRHPart();
    //baconhep::TRHPart *pRH = (baconhep::TRHPart*) rArray[index];

    //assert(fRHParArr->GetEntries() < fRHParArr->GetSize());
    //const int index = fRHParArr->GetEntries();
    //std::cout<<"index "<<index<<std::endl;
    //new(fRHParArr[index]) baconhep::TRHPart();
    //baconhep::TRHPart *pRH = (baconhep::TRHPart*) fRHParArr[index];

    //std::cout<<"rechit "<<r<<" energy "<<hRHT->at(r).energy()<<" parent part = "<<hRHT->at(r).parentPart()<<" pos = "<<hRHT->at(r).x()<<", "<<hRHT->at(r).y()<<", "<<hRHT->at(r).z()<<std::endl;
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
    value_simHit_energy[r] = hRHT->at(r).simHitEnergy();
    value_simHit_depth[r] = hRHT->at(r).simHitDepth();
    value_simHit_time[r] = hRHT->at(r).simHitTime();
    //baconhep::TRHPart pRH;
    //pRH->timeFalling   = hRHT->at(r).timeFalling();
    //pRH->energy = hRHT->at(r).energy();
    //pRH->eta    = hRHT->at(r).eta();
    //pRH->phi    = hRHT->at(r).phi();
    //pRH->time   = hRHT->at(r).time();
    //pRH->x      = hRHT->at(r).x();
    //pRH->y      = hRHT->at(r).y();
    //pRH->z      = hRHT->at(r).z();
    //pRH->time      = hRHT->at(r).time();
    //pRH->trackId = hRHT->at(r).parentPart();
    //fRHParArr[r] = pRH;
    //fRHParArr->push_back(pRH);
    //pRH.depth  = itRH->id().depth();
    //pRH.genE   = genE;

    value_hit_n++;
  }

  //for(HBHERecHitCollection::const_iterator itRH = recHitHCAL->begin(); itRH != recHitHCAL->end(); itRH++) {

  /*
  edm::ESHandle<HcalDDDRecConstants> pHRNDC;
  iSetup.get<HcalRecNumberingRecord>().get( pHRNDC );
  fRecNumber= &*pHRNDC;
  
 
  edm::ESHandle<HcalTrigTowerGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  //Get the Gen Info
  edm::Handle<edm::PCaloHitContainer> hSimHits;      // create handle
  iEvent.getByToken(fSHitToken, hSimHits);   // SimHits

  edm::ESHandle<CaloGeometry> geometry;
  iEvent.getByToken(fSHitToken, hSimHits);   // SimHits
  iSetup.get<HcalRecNumberingRecord>().get( pHRNDC );
  iSetup.get<CaloGeometryRecord>().get( geometry );
  
  edm::PCaloHitContainer lSimHits  = *hSimHits;
  fGeometry = &*geometry;
  fRecNumber= &*pHRNDC;

  edm::Handle<edm::SimTrackContainer> hSimTracks;
  iEvent.getByToken(fSimTrackToken, hSimTracks);
  edm::SimTrackContainer lSimTracks = *hSimTracks;
  
  edm::Handle<edm::SimVertexContainer> hSimVertex;
  iEvent.getByToken(fSimVertexToken, hSimVertex);
  edm::SimVertexContainer lSimVertex = *hSimVertex;

  edm::Handle<edm::PCaloHitContainer> hSimHitsEE;      // create handle
  iEvent.getByToken(fSHitTokenEE, hSimHitsEE);   // SimHits
  edm::PCaloHitContainer lSimHitsEE  = *hSimHitsEE;

  edm::Handle<edm::PCaloHitContainer> hSimHitsEB;      // create handle
  iEvent.getByToken(fSHitTokenEB, hSimHitsEB);   // SimHits
  edm::PCaloHitContainer lSimHitsEB  = *hSimHitsEB;
  */
  // fRHParArr->Clear();
  //fGenParticleArr->Clear();
  //fGenEvtInfo->Clear();
  //fEvtInfo->Clear();
  //fGenJetArr->Clear();
  //fGenFatJetArr->Clear();
  //fECALPFClusterArr->Clear();
  //fTrackArr->Clear();
  //fVertexArr->Clear();
  //const reco::Vertex *pv = 0;
  //fFillerGenInfo->fill(fGenEvtInfo,fGenParticleArr,0,iEvent,0);
  //fFillerGenJets->fill(fGenJetArr,fGenFatJetArr,iEvent);
  //fFillerEventInfo->fill(fEvtInfo,fECALPFClusterArr,iEvent,*pv,0,0);
  //fFillerRH->fill(fRHParArr,iEvent,iSetup,lSimHits,lSimHitsEE,lSimHitsEB,fRecNumber);
  //fFillerSimTrack->fill(fTrackArr, iEvent,iSetup, lSimTracks);
  //fFillerSimVertex->fill(fVertexArr, iEvent,iSetup, lSimVertex);

  fTree->Fill();
  // Step C.1: Run FE Format Error / ZS for real data.
  //iEvent.put(std::move(result));
}

//DEFINE_FWK_MODULE(HcalTrigPrimDigiNtupler);
