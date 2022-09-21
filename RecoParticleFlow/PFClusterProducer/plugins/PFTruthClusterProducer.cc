#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterBuilderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/RecHitTopologicalCleanerBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/SeedFinderBase.h"


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

#include "DataFormats/HcalRecHit/interface/HBHERecHitTruth.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"

#include "SimCalorimetry/CaloSimAlgos/interface/CaloHitResponse.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"

#include <memory>

class PFTruthClusterProducer : public edm::stream::EDProducer<> {
  typedef RecHitTopologicalCleanerBase RHCB;
  typedef InitialClusteringStepBase ICSB;
  typedef PFClusterBuilderBase PFCBB;
  typedef PFCPositionCalculatorBase PosCalc;

public:
  PFTruthClusterProducer(const edm::ParameterSet&);
  ~PFTruthClusterProducer() override = default;

  void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  // inputs
  edm::EDGetTokenT<reco::PFRecHitCollection> _rechitsLabel;
  // options
  const bool _prodInitClusters;
  // the actual algorithm
  std::vector<std::unique_ptr<RecHitTopologicalCleanerBase> > _cleaners;
  std::vector<std::unique_ptr<RecHitTopologicalCleanerBase> > _seedcleaners;
  std::unique_ptr<SeedFinderBase> _seedFinder;
  std::unique_ptr<InitialClusteringStepBase> _initialClustering;
  std::unique_ptr<PFClusterBuilderBase> _pfClusterBuilder;
  std::unique_ptr<PFCPositionCalculatorBase> _positionReCalc;
  std::unique_ptr<PFClusterEnergyCorrectorBase> _energyCorrector;

  std::unique_ptr<PFCPositionCalculatorBase> _allCellsPosCalc;

  const HcalDDDRecConstants *fRecNumber;
  edm::EDGetTokenT<edm::PCaloHitContainer> fSHitToken;
  edm::EDGetTokenT<edm::SimTrackContainer> SimTrackToken;
  edm::EDGetTokenT<edm::SimVertexContainer> SimVertexToken;
  //edm::EDGetTokenT<edm::SortedCollection<HBHERecHit> > _rechitsLabel;
  edm::ESGetToken<HcalDDDRecConstants, HcalRecNumberingRecord> hdrcToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
  edm::EDGetTokenT<reco::PFClusterCollection> _clustersLabel;

  HcalSimParameterMap fSimParameterMap;
  CaloHitResponse* fResponse;

protected:
  reco::PFRecHitRef makeRefhit(const edm::Handle<reco::PFRecHitCollection>& h, const unsigned i) const {
    return reco::PFRecHitRef(h, i);
  }
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFTruthClusterProducer);

#ifdef PFLOW_DEBUG
#define LOGVERB(x) edm::LogVerbatim(x)
#define LOGWARN(x) edm::LogWarning(x)
#define LOGERR(x) edm::LogError(x)
#define LOGDRESSED(x) edm::LogInfo(x)
#else
#define LOGVERB(x) LogTrace(x)
#define LOGWARN(x) edm::LogWarning(x)
#define LOGERR(x) edm::LogError(x)
#define LOGDRESSED(x) LogDebug(x)
#endif

PFTruthClusterProducer::PFTruthClusterProducer(const edm::ParameterSet& conf)
    : _prodInitClusters(conf.getUntrackedParameter<bool>("prodInitialClusters", false)) {
  _rechitsLabel = consumes<reco::PFRecHitCollection>(conf.getParameter<edm::InputTag>("recHitsSource"));
  edm::ConsumesCollector cc = consumesCollector();

  const edm::ParameterSet& acConf = conf.getParameterSet("allCellsPositionCalc");
  const std::string& algoac = acConf.getParameter<std::string>("algoName");
  _allCellsPosCalc = PFCPositionCalculatorFactory::get()->create(algoac, acConf, cc);

  //std::cout<<"In PFTruthClusterProducer"<<std::endl;

  //setup rechit cleaners
  const edm::VParameterSet& cleanerConfs = conf.getParameterSetVector("recHitCleaners");
  for (const auto& conf : cleanerConfs) {
    const std::string& cleanerName = conf.getParameter<std::string>("algoName");
    _cleaners.emplace_back(RecHitTopologicalCleanerFactory::get()->create(cleanerName, conf, cc));
  }

  if (conf.exists("seedCleaners")) {
    const edm::VParameterSet& seedcleanerConfs = conf.getParameterSetVector("seedCleaners");

    for (const auto& conf : seedcleanerConfs) {
      const std::string& seedcleanerName = conf.getParameter<std::string>("algoName");
      _seedcleaners.emplace_back(RecHitTopologicalCleanerFactory::get()->create(seedcleanerName, conf, cc));
    }
  }


  //setup (possible) recalcuation of positions
  const edm::ParameterSet& pConf = conf.getParameterSet("positionReCalc");
  if (!pConf.empty()) {
    const std::string& pName = pConf.getParameter<std::string>("algoName");
    _positionReCalc = PFCPositionCalculatorFactory::get()->create(pName, pConf, cc);
  }
  // see if new need to apply corrections, setup if there.
  const edm::ParameterSet& cConf = conf.getParameterSet("energyCorrector");
  if (!cConf.empty()) {
    const std::string& cName = cConf.getParameter<std::string>("algoName");
    _energyCorrector = PFClusterEnergyCorrectorFactory::get()->create(cName, cConf);
  }

  //now do what ever other initialization is needed
  fSHitToken = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));
  SimTrackToken = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits",""));
  SimVertexToken = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits",""));
  //_rechitsLabel = consumes<edm::SortedCollection<HBHERecHit> >(conf.getParameter<edm::InputTag>("pf_src"));
  hdrcToken_ = esConsumes<HcalDDDRecConstants, HcalRecNumberingRecord>();
  //_clustersLabel = consumes<reco::PFClusterCollection>(conf.getParameter<edm::InputTag>("clus_src"));
  geomToken_ = esConsumes();

  if (_prodInitClusters) {
    produces<reco::PFClusterCollection>("initialClusters");
  }
  produces<reco::PFClusterCollection>();
}

void PFTruthClusterProducer::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& es) {
  for (const auto& cleaner : _cleaners)
    cleaner->update(es);
  for (const auto& cleaner : _seedcleaners)
    cleaner->update(es);
}

void PFTruthClusterProducer::produce(edm::Event& e, const edm::EventSetup& es) {
  using namespace edm;

  edm::Handle<reco::PFRecHitCollection> rechits;
  e.getByToken(_rechitsLabel, rechits);

  std::vector<bool> mask(rechits->size(), true);
  for (const auto& cleaner : _cleaners) {
    cleaner->clean(rechits, mask);
  }

  auto nhits = rechits->size();
  int nhits2_ = 0;
  for (unsigned i = 0; i < nhits; ++i) {
    if (!mask[i]){ 
      continue;  // cannot seed masked objects
    }
    nhits2_++;
  }
  
  // no seeding on these hits
  std::vector<bool> seedmask = mask;
  for (const auto& cleaner : _seedcleaners) {
    cleaner->clean(rechits, seedmask);
  }

  edm::ESHandle<CaloGeometry> geoHandle = es.getHandle(geomToken_);
  const CaloSubdetectorGeometry* hcalBarrelGeo = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  const CaloSubdetectorGeometry* hcalEndcapGeo = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalEndcap);

  const HcalDDDRecConstants* pHRNDC;
  pHRNDC = &es.getData(hdrcToken_);
  fRecNumber = &*pHRNDC;

  edm::Handle<edm::SimTrackContainer> hSimTracks;      // create handle
  e.getByToken(SimTrackToken, hSimTracks);   // SimTracks

  edm::Handle<edm::SimVertexContainer> hSimVertices;      // create handle
  e.getByToken(SimVertexToken, hSimVertices);   // SimVertices

  edm::Handle<edm::PCaloHitContainer> hSimHits;      // create handle
  e.getByToken(fSHitToken, hSimHits);   // SimHits


  edm::PCaloHitContainer lSimHits  = *hSimHits;

  std::map<unsigned int, int> detId_map; //maps keys = detectorID and values = Sim Track geantID
  std::map<unsigned int, float> simHit_energy_map; //maps keys = detectorID and values = Sim Track energy
  std::map<unsigned int, int> simHit_depth_map; //maps keys = detectorID and values = Sim Track depth
  std::map<unsigned int, float> simHit_time_map; //maps keys = detectorID and values = Sim Track time

  fResponse = new CaloHitResponse(NULL, (CaloShapes*)NULL);
  fResponse->setGeometry(&*geoHandle);

  for (int j=0; j < (int) lSimHits.size(); j++) {
    double samplingFactor = 1;
    HcalDetId simId = HcalHitRelabeller::relabel((lSimHits)[j].id(), fRecNumber);
    if(simId.subdet() == HcalBarrel) {
      samplingFactor = fSimParameterMap.hbParameters().samplingFactor(simId);
    } else if (simId.subdet() == HcalEndcap) {
      samplingFactor = fSimParameterMap.heParameters().samplingFactor(simId);
    }
    if(!detId_map.count(simId.rawId())){ //new detID
      detId_map[simId.rawId()] = (lSimHits)[j].geantTrackId(); //associate simId to highest energy geantTrackId
      simHit_energy_map[simId.rawId()] = samplingFactor*((lSimHits)[j].energy());
      HcalDetId simIdForDepth = (lSimHits)[j].id();
      simHit_depth_map[simId.rawId()] = simIdForDepth.depth();
    } else if( (lSimHits)[j].energy() > (lSimHits)[detId_map[simId.rawId()]].energy() ){ //higher energy deposit than previous max
      detId_map[simId.rawId()] = (lSimHits)[j].geantTrackId();
      simHit_energy_map[simId.rawId()] += samplingFactor*((lSimHits)[j].energy());
      HcalDetId simIdForDepth = (lSimHits)[j].id();
      simHit_depth_map[simId.rawId()] = simIdForDepth.depth();
    } else{ //lower energy deposit than max energy deposit
      simHit_energy_map[simId.rawId()] += samplingFactor*((lSimHits)[j].energy());
    }
  }


  std::map<unsigned int, unsigned int> track_map; //map keys = Sim Track geantID and values = Sim Track index
  for (unsigned int i = 0; i < hSimTracks->size(); ++i) {
    track_map[hSimTracks->at(i).trackId()] = i;
  }




  int dummyParent = -1; //a dummy parent particle for noise hits
  std::vector<int> parentIndices;

  for (const auto& erh : *rechits) {

    if(!detId_map.count(erh.detId())){ //no simHit for recHit
      parentIndices.push_back(dummyParent);
      dummyParent = dummyParent-1; //each noise hit can have a unique ID; do I need to do this for the SPVCNN version?!
      continue;
    }

    unsigned int steps = 0;
    unsigned int vertValue = hSimTracks->at(track_map[detId_map[erh.detId()]]).vertIndex(); //for hbhe rec hits
    int parentPart = hSimVertices->at(vertValue).parentIndex();
    if(parentPart == -1){
      parentPart = detId_map[erh.detId()];
    }
    while(steps < 100 && parentPart >= 0){
      vertValue = hSimTracks->at(track_map[parentPart]).vertIndex();
      if(hSimVertices->at(vertValue).parentIndex() < 0){
	break;
      }
      parentPart = hSimVertices->at(vertValue).parentIndex();
      steps++;
    }
    
    parentIndices.push_back(parentPart);
    

  }




  std::map<int, reco::PFCluster> cluster_map; //map keys = cluster index from SPVCNN and values = clusters that may be added to

  auto pfClusters = std::make_unique<reco::PFClusterCollection>();
  for (int i = 0; i < nhits2_; ++i) { //this probably isn't correct due to masking?...
    reco::PFCluster current;
    auto ref = makeRefhit(rechits, i); //reco::PFRecHitRef(h, i);
    current.addRecHitFraction(reco::PFRecHitFraction(ref, 1));

    const auto rh_energy = ref->energy();
    current.setSeed(ref->detId());
    current.setEnergy(rh_energy);
    current.setTime(ref->time());
    current.setLayer(ref->layer());
    current.setPosition(math::XYZPoint(ref->position().x(), ref->position().y(), ref->position().z()));
    current.calculatePositionREP();
    current.setDepth(ref->depth());

    if( !cluster_map.count( parentIndices.at(i) ) ){ //new cluster index
      cluster_map[parentIndices.at(i)] = current;
    } else{ //code for combining clusters taken from PFMultiDepthClusterizer
      double e1 = 0.0;
      double e2 = 0.0;
      reco::PFCluster main = cluster_map[parentIndices.at(i)];
      for (const auto& fraction : main.recHitFractions()){
	if (fraction.recHitRef()->detId() == main.seed()) {
	  e1 = fraction.recHitRef()->energy();
	}
      }
      for (const auto& fraction : current.recHitFractions()) {
	main.addRecHitFraction(fraction);
	if (fraction.recHitRef()->detId() == current.seed()) {
	  e2 = fraction.recHitRef()->energy();
	}
      }
      if (e2 > e1){
	main.setSeed(current.seed());
      }
      cluster_map[parentIndices.at(i)] = main;
    }

  }

  for (auto const& cl : cluster_map){
    _allCellsPosCalc->calculateAndSetPosition(cluster_map[cl.first]);
    if(cl.second.energy() > 0){
      pfClusters->push_back(cl.second);
    }
  }


  if (_positionReCalc) {
    _positionReCalc->calculateAndSetPositions(*pfClusters);
  }

  if (_energyCorrector) {
    _energyCorrector->correctEnergies(*pfClusters);
  }

  //if (_prodInitClusters)
  //e.put(std::move(initialClusters), "initialClusters");
  e.put(std::move(pfClusters));

}
