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

//#include "DataFormats/HcalRecHit/interface/HBHERecHitTruth.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"

#include "SimCalorimetry/CaloSimAlgos/interface/CaloHitResponse.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"

#include <memory>

class PFTruthClusterProducer2 : public edm::stream::EDProducer<> {
  typedef RecHitTopologicalCleanerBase RHCB;
  typedef InitialClusteringStepBase ICSB;
  typedef PFClusterBuilderBase PFCBB;
  typedef PFCPositionCalculatorBase PosCalc;

public:
  PFTruthClusterProducer2(const edm::ParameterSet&);
  ~PFTruthClusterProducer2() override = default;

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
  std::unique_ptr<PFClusterBuilderBase> _pfClusterBuilderDepth;
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

  bool buildTopoCluster(const edm::Handle<reco::PFRecHitCollection>&,
                        const std::vector<bool>&,  // masked rechits
                        unsigned int,              //present rechit
                        std::vector<bool>&,        // hit usage state
                        reco::PFCluster&,         // the topocluster
			float,         // the simmed energy of the simParticle
			int,         // the hit depth
			int,         // the alpha rec hit
			std::vector<int>);         // allowed rec hits

  const bool _useCornerCells = true;
  void buildTopoCluster2(const edm::Handle<reco::PFRecHitCollection>&,
                        const std::vector<bool>&,  // masked rechits
                        unsigned int,              //present rechit
                        std::vector<bool>&,        // hit usage state
                        reco::PFCluster&);         // the topocluster

  float calculateDR(float alphaEta, float newEta, float alphaPhi, float newPhi){
    double deleta = (alphaEta - newEta);
    double delphi, ogphi, newphi;
    newphi = newPhi;
    ogphi = alphaPhi;
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
    return sqrt( deleta*deleta + delphi*delphi );
  }


protected:
  reco::PFRecHitRef makeRefhit(const edm::Handle<reco::PFRecHitCollection>& h, const unsigned i) const {
    return reco::PFRecHitRef(h, i);
  }
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFTruthClusterProducer2);

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

PFTruthClusterProducer2::PFTruthClusterProducer2(const edm::ParameterSet& conf)
    : _prodInitClusters(conf.getUntrackedParameter<bool>("prodInitialClusters", false)) {
  _rechitsLabel = consumes<reco::PFRecHitCollection>(conf.getParameter<edm::InputTag>("recHitsSource"));
  edm::ConsumesCollector cc = consumesCollector();

  const edm::ParameterSet& acConf = conf.getParameterSet("allCellsPositionCalc");
  const std::string& algoac = acConf.getParameter<std::string>("algoName");
  _allCellsPosCalc = PFCPositionCalculatorFactory::get()->create(algoac, acConf, cc);

  //std::cout<<"In PFTruthClusterProducer2"<<std::endl;

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

  const edm::ParameterSet& sfConf = conf.getParameterSet("seedFinder");
  const std::string& sfName = sfConf.getParameter<std::string>("algoName");
  _seedFinder = SeedFinderFactory::get()->create(sfName, sfConf);

  //setup topo cluster builder
  const edm::ParameterSet& initConf = conf.getParameterSet("initialClusteringStep");
  const std::string& initName = initConf.getParameter<std::string>("algoName");
  _initialClustering = InitialClusteringStepFactory::get()->create(initName, initConf, cc);

  const edm::ParameterSet& pfcConf = conf.getParameterSet("pfClusterBuilder");
  if (!pfcConf.empty()) {
    const std::string& pfcName = pfcConf.getParameter<std::string>("algoName");
    _pfClusterBuilder = PFClusterBuilderFactory::get()->create(pfcName, pfcConf, cc);
  }

  const edm::ParameterSet& pfcConfDepth = conf.getParameterSet("pfClusterBuilderDepth");
  if (!pfcConfDepth.empty()) {
    const std::string& pfcNameDepth = pfcConfDepth.getParameter<std::string>("algoName");
    _pfClusterBuilderDepth = PFClusterBuilderFactory::get()->create(pfcNameDepth, pfcConfDepth, cc);
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

void PFTruthClusterProducer2::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& es) {
  _initialClustering->update(es);
  if (_pfClusterBuilder)
    _pfClusterBuilder->update(es);
  if (_pfClusterBuilderDepth)
    _pfClusterBuilderDepth->update(es);
  for (const auto& cleaner : _cleaners)
    cleaner->update(es);
  for (const auto& cleaner : _seedcleaners)
    cleaner->update(es);
}

void PFTruthClusterProducer2::produce(edm::Event& e, const edm::EventSetup& es) {
  using namespace edm;

  _initialClustering->reset();
  if (_pfClusterBuilder)
    _pfClusterBuilder->reset();
  if (_pfClusterBuilderDepth)
    _pfClusterBuilderDepth->reset();

  edm::Handle<reco::PFRecHitCollection> rechits;
  e.getByToken(_rechitsLabel, rechits);

  _initialClustering->updateEvent(e);

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
















  std::vector<bool> used(rechits->size(), false);







  /*


  std::vector<bool> maskHB(rechits->size(), true);
  for (const auto& cleaner : _cleaners) {
    cleaner->clean(rechits, maskHB);
  }

  for (unsigned int i = 0; i < rechits->size(); ++i) {
    //if (used[i] || rechits->at(i).energy() < 4.5){ //0.8 seems to have no effect relative to default.  I've tried 1.6.  100 should exclude basically everything
    if (used[i]){// || rechits->at(i).idFront()){
      maskHB[i] = 0;
    }
    const HcalDetId detid = rechits->at(i).detId();//rechits->at(i).idFront();
    //rechits->at(i).detId()
    HcalSubdetector esd = (HcalSubdetector)detid.subdetId();
    if(esd == HcalEndcap){
      maskHB[i] = 0;
    }
    if(esd == HcalBarrel){
      used[i] = 1;
    }
  }

  std::vector<bool> seedmaskHB = maskHB;
  for (const auto& cleaner : _seedcleaners) {
    cleaner->clean(rechits, seedmaskHB);
  }

  std::vector<bool> seedableHB(rechits->size(), false);
  _seedFinder->findSeeds(rechits, seedmaskHB, seedableHB);

  auto initialClustersHB = std::make_unique<reco::PFClusterCollection>();
  _initialClustering->buildClusters(rechits, maskHB, seedableHB, *initialClustersHB);
  LOGVERB("PFClusterProducer::produce()") << *_initialClustering;

  auto pfClustersHB = std::make_unique<reco::PFClusterCollection>();
  pfClustersHB = std::make_unique<reco::PFClusterCollection>();
  if (_pfClusterBuilder) {  // if we've defined a re-clustering step execute it
    _pfClusterBuilder->buildClusters(*initialClustersHB, seedableHB, *pfClustersHB);
    LOGVERB("PFClusterProducer::produce()") << *_pfClusterBuilder;
  } else {
    pfClustersHB->insert(pfClustersHB->end(), initialClustersHB->begin(), initialClustersHB->end());
  }

  std::cout<<"finished splitting initial HB clusters. "<<pfClustersHB->size()<<std::endl;

  std::vector<bool> seedable2HB;
  auto pfClusters3HB = std::make_unique<reco::PFClusterCollection>();
  _pfClusterBuilderDepth->buildClusters(*pfClustersHB, seedable2HB, *pfClusters3HB);
  //to insert generic clustered leftovers:

  std::cout<<"finished building nontruth seeded HB clusters. "<<pfClusters3HB->size()<<std::endl; 

  //pfClusters->insert(pfClusters->end(), pfClusters3->begin(), pfClusters3->end());


  */





















  //auto const& hits = ;
  //std::vector<bool> used(rechits->size(), false);
  std::vector<unsigned int> seeds;

  // get the seeds and sort them descending in energy
  seeds.reserve(rechits->size());
  for (unsigned int i = 0; i < rechits->size(); ++i) {
    //if (!rechitMask[i] || !seedable[i] || used[i])
    if (used[i])
      continue;
    seeds.emplace_back(i);
  }
  // maxHeap would be better
  std::sort( seeds.begin(), seeds.end(), [&](unsigned int i, unsigned int j) { return rechits->at(i).energy() > rechits->at(j).energy(); });

  //std::cout<<"PRINT ORDERED SEEDS ("<<seeds.size()<<")"<<std::endl;
  //for(unsigned int j = 0; j < seeds.size(); j++){
  //  std::cout<<" "<<seeds.at(j);
  //}
  //std::cout<<"\n";

  int sum_of_elems1 = 0;
  for(std::vector<bool>::iterator it = used.begin(); it != used.end(); ++it){
    if(*it){
      sum_of_elems1 += 1;
    }
  }
  //std::cout<<"number of unused rechits before starting: "<<sum_of_elems1<<std::endl;


  std::map<unsigned int, int> detID_toRecHitIndex; //maps keys = detectorID and values = rechit index, which should be iterable as seeds
  for (unsigned int i = 0; i < rechits->size(); ++i) {
    detID_toRecHitIndex[rechits->at(i).detId()] = i;
  }




  edm::PCaloHitContainer lSimHits  = *hSimHits;

  std::map<unsigned int, int> detId_map; //maps keys = detectorID and values = Sim Track geantID
  std::map<unsigned int, int> simhitIndex_map; //maps keys = detectorID and values = simHit index that already corresponds to that det ID
  std::map<unsigned int, float> simHit_energy_map; //maps keys = detectorID and values = Sim Track energy
  std::map<unsigned int, int> simHit_depth_map; //maps keys = detectorID and values = Sim Track depth
  std::map<unsigned int, float> simHit_time_map; //maps keys = detectorID and values = Sim Track time

  struct simEnergyStruct{
    int ID;
    float energy;
  };
  std::map<int, simEnergyStruct> simID_to_energy_map; //maps keys = Sim Track geantID and values = total simmed energy for this geantID
  std::map<int, std::vector<int>> simID_to_rechitindex_map; //maps keys = Sim Track geantID and values = index of rechits containing to this simID


  std::map<unsigned int, std::map<int, float>> detID_to_simID_to_energy_map; //maps keys = detectorID and values =  map with keys = Sim Track geantID and values = simmed enegy for this geantID (for the cell in question)


  fResponse = new CaloHitResponse(NULL, (CaloShapes*)NULL);
  fResponse->setGeometry(&*geoHandle);

  for (int j=0; j < (int) lSimHits.size(); j++) {
    std::cout<<"simhit j = "<<j<<" (lSimHits)[j].geantTrackId() = "<<(lSimHits)[j].geantTrackId()<<std::endl;
    double samplingFactor = 1;
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

    if(simId.subdet() == HcalBarrel) {
      samplingFactor = fSimParameterMap.hbParameters().samplingFactor(simId);
    } else if (simId.subdet() == HcalEndcap) {
      samplingFactor = fSimParameterMap.heParameters().samplingFactor(simId);
    }


    if(!simID_to_energy_map.count((lSimHits)[j].geantTrackId())){
      simEnergyStruct newstruct;
      //newstruct.ID = simId.rawId();
      newstruct.ID = (lSimHits)[j].geantTrackId();
      newstruct.energy = samplingFactor*((lSimHits)[j].energy());
      simID_to_energy_map[(lSimHits)[j].geantTrackId()] = newstruct;
    }
    else{
      simID_to_energy_map[(lSimHits)[j].geantTrackId()].energy += samplingFactor*((lSimHits)[j].energy());
    }


    //simID_to_rechitindex_map
    if(!simID_to_rechitindex_map.count((lSimHits)[j].geantTrackId())){
      std::vector<int> dummy;
      dummy.push_back(detID_toRecHitIndex[simId.rawId()]);
      simID_to_rechitindex_map[(lSimHits)[j].geantTrackId()] = dummy;
    }
    else{
      if( !( std::count(simID_to_rechitindex_map[(lSimHits)[j].geantTrackId()].begin(), simID_to_rechitindex_map[(lSimHits)[j].geantTrackId()].end(), detID_toRecHitIndex[simId.rawId()]) ) ){
	simID_to_rechitindex_map[(lSimHits)[j].geantTrackId()].push_back(detID_toRecHitIndex[simId.rawId()]);
      }
    }
    


    if(!detID_to_simID_to_energy_map.count(simId.rawId())){ //new detID
      std::map<int, float> dummy;
      dummy[(lSimHits)[j].geantTrackId()] = samplingFactor*((lSimHits)[j].energy());
      detID_to_simID_to_energy_map[simId.rawId()] = dummy;
    }
    else if(!detID_to_simID_to_energy_map[simId.rawId()].count((lSimHits)[j].geantTrackId())){
      detID_to_simID_to_energy_map[simId.rawId()][(lSimHits)[j].geantTrackId()] = samplingFactor*((lSimHits)[j].energy());
    }
    else{
      detID_to_simID_to_energy_map[simId.rawId()][(lSimHits)[j].geantTrackId()] += samplingFactor*((lSimHits)[j].energy());
    }

    if(!detId_map.count(simId.rawId())){ //new detID
      std::cout<<"new simhit in this detID, which is "<<simId.rawId()<<" with energy "<<((lSimHits)[j].energy())<<" and samplingFactor = "<<samplingFactor<<std::endl;
      detId_map[simId.rawId()] = (lSimHits)[j].geantTrackId(); //associate simId to highest energy geantTrackId
      simhitIndex_map[simId.rawId()] = j;
      simHit_energy_map[simId.rawId()] = samplingFactor*((lSimHits)[j].energy());
      HcalDetId simIdForDepth = (lSimHits)[j].id();
      simHit_depth_map[simId.rawId()] = simIdForDepth.depth();
    } 
    else if((lSimHits)[j].energy() > (lSimHits)[simhitIndex_map[simId.rawId()]].energy() ){
      //else if( (lSimHits)[j].energy() > (lSimHits)[detId_map[simId.rawId()]].energy() ){ //higher energy deposit than previous max
      //std::cout<<"Replacing a simhit in detId_map which had geantid = "<<detId_map[simId.rawId()]<<". new track id = "<<(lSimHits)[j].geantTrackId()<<". Old energy was "<<(lSimHits)[detId_map[simId.rawId()]].energy()<<" new energy is "<<(lSimHits)[j].energy()<<" well with samplingFactor = "<<samplingFactor<<std::endl;
      std::cout<<"Replacing a simhit in detId_map which had geantid = "<<detId_map[simId.rawId()]<<". new track id = "<<(lSimHits)[j].geantTrackId()<<". Old energy was "<<(lSimHits)[simhitIndex_map[simId.rawId()]].energy()<<" new energy is "<<(lSimHits)[j].energy()<<" well with samplingFactor = "<<samplingFactor<<std::endl;
      detId_map[simId.rawId()] = (lSimHits)[j].geantTrackId();
      simhitIndex_map[simId.rawId()] = j;
      simHit_energy_map[simId.rawId()] += samplingFactor*((lSimHits)[j].energy());
      HcalDetId simIdForDepth = (lSimHits)[j].id();
      simHit_depth_map[simId.rawId()] = simIdForDepth.depth();
    } else{ //lower energy deposit than max energy deposit
      std::cout<<"simhit not enough energy to take over hit, so the track id is still "<<detId_map[simId.rawId()]<<".  the new track id is "<<(lSimHits)[j].geantTrackId()<<", which has energy "<<(lSimHits)[j].energy()<<" and samplingFactor = "<<samplingFactor<<" by the way the comparison was between original energy: "<<(lSimHits)[simhitIndex_map[simId.rawId()]].energy()<<" and new energy "<<(lSimHits)[j].energy()<<std::endl;
      simHit_energy_map[simId.rawId()] += samplingFactor*((lSimHits)[j].energy());
    }
  }

  std::vector<simEnergyStruct> simEnergyStructVector;
  std::map<int,simEnergyStruct>::iterator r;
  for (r = simID_to_energy_map.begin(); r!= simID_to_energy_map.end(); r++){
    simEnergyStructVector.push_back((*r).second);
  }
  //std::sort(simEnergyStructVector.begin(), simEnergyStructVector.end(), [&](unsigned int i, unsigned int j) { return hits[i].energy() > hits[j].energy(); });
  std::sort(simEnergyStructVector.begin(), simEnergyStructVector.end(), [](const simEnergyStruct & a, const simEnergyStruct & b) -> bool{ return a.energy > b.energy; });




  std::map<unsigned int, unsigned int> track_map; //map keys = Sim Track geantID and values = Sim Track index
  for (unsigned int i = 0; i < hSimTracks->size(); ++i) {
    std::cout<<"simTrack "<<i<<" type = "<<hSimTracks->at(i).type()<<" pt = "<<hSimTracks->at(i).momentum().Pt()<<" eta = "<<hSimTracks->at(i).momentum().Eta()<<" phi = "<<hSimTracks->at(i).momentum().Phi()<<" E = "<<hSimTracks->at(i).momentum().E()<<" geantID = "<<hSimTracks->at(i).genpartIndex()<<" vertID = "<<hSimTracks->at(i).vertIndex()<<" trackID = "<<hSimTracks->at(i).trackId()<<std::endl; //WPM added late
    track_map[hSimTracks->at(i).trackId()] = i;
  }

  //for (unsigned i = 0; i < hSimVertices->size(); ++i) {
  //  std::cout<<"simVertex "<<i<<" vertexID = "<<hSimVertices->at(i).vertexId()<<" processtype = "<<hSimVertices->at(i).processType()<<" parentID = "<<hSimVertices->at(i).parentIndex()<<" position = "<<hSimVertices->at(i).position().X()<<" "<<hSimVertices->at(i).position().Y()<<" "<<hSimVertices->at(i).position().Z()<<std::endl;
  //} //WPM added late

  std::cout<<"PRINT SIMIDS AND SIMMEDENERGIES"<<std::endl;
  for(unsigned int i = 0; i < simEnergyStructVector.size(); i++){
    std::cout<<simEnergyStructVector.at(i).ID<<" "<<simEnergyStructVector.at(i).energy<<std::endl;

    std::vector<int> dummy;

    //  std::sort( seeds.begin(), seeds.end(), [&](unsigned int i, unsigned int j) { return rechits->at(i).energy() > rechits->at(j).energy(); });
    //std::sort( simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].begin(), simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].end(), [&](unsigned int i, unsigned int j) { return rechits->at(i).energy() > rechits->at(j).energy(); } );
    std::sort( simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].begin(), simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].end(), [&](unsigned int j, unsigned int k) { return detID_to_simID_to_energy_map[rechits->at(j).detId()][simEnergyStructVector.at(i).ID] > detID_to_simID_to_energy_map[rechits->at(k).detId()][simEnergyStructVector.at(i).ID]; } );

    float simeng = 0.;
    float recoeng = 0.;
    float simexclude = 0.;
    float recoexclude = 0.;
    for(unsigned int j = 0; j < simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].size(); j++){
      //std::cout<<" "<<simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)<<" ("<<rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).energy()<<")";
      
      std::cout<<" "<<simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)<<" ("<<rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).energy()<<", "<<detID_to_simID_to_energy_map[rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).detId()][simEnergyStructVector.at(i).ID]<<")";
      
      if(detID_to_simID_to_energy_map[rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).detId()][simEnergyStructVector.at(i).ID] > 0.01*simEnergyStructVector.at(i).energy){
	simeng += detID_to_simID_to_energy_map[rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).detId()][simEnergyStructVector.at(i).ID];
	recoeng += rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).energy();
      }

      //if( detID_to_simID_to_energy_map[rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).detId()][simEnergyStructVector.at(i).ID] / rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).energy() < 0.33){
	//simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].erase(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].begin()+j, simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].begin()+j+1);
	///std::cout<<" THIS IS A BAD BOY ";
	//continue;
      //}

      //if(detID_to_simID_to_energy_map[rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).detId()][simEnergyStructVector.at(i).ID] < 0.05*simEnergyStructVector.at(i).energy){

      if(detID_to_simID_to_energy_map[rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).detId()][simEnergyStructVector.at(i).ID] < 0.01*simEnergyStructVector.at(i).energy){
	std::cout<<" START EXCLUSION";
	for(unsigned int j2 = j; j2 < simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].size(); j2++){
	  std::cout<<" "<<simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j2)<<" ("<<rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j2)).energy()<<", "<<detID_to_simID_to_energy_map[rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j2)).detId()][simEnergyStructVector.at(i).ID]<<")";
	  simexclude += detID_to_simID_to_energy_map[rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j2)).detId()][simEnergyStructVector.at(i).ID];
	  recoexclude += rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j2)).energy();
	}
	simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].erase(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].begin()+j, simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].end());
	break;
      }

      dummy.push_back(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j));

    }
    std::cout<<"\n";
    std::cout<<"                  simeng in core = "<<simeng<<", which is "<<simeng/simEnergyStructVector.at(i).energy<<" of truth.  Recoeng in core = "<<recoeng<<" which is "<<recoeng/simEnergyStructVector.at(i).energy<<" of truth."<<std::endl;
    std::cout<<"                  simexclude from core = "<<simexclude<<", which is "<<simexclude/simEnergyStructVector.at(i).energy<<" of truth.  Recoexclude from core = "<<recoexclude<<" which is "<<recoexclude/simEnergyStructVector.at(i).energy<<" of truth."<<std::endl;

    simID_to_rechitindex_map[simEnergyStructVector.at(i).ID] = dummy;

  }
  std::cout<<"PRINT DETIDS AND SIMIDS AND SIMMEDENERGIES"<<std::endl;
  std::map<unsigned int, std::map<int, float>>::iterator r1;
  float totalSimEnergy = 0.;
  float totalSimEnergy1 = 0.;
  float totalSimEnergy2 = 0.;
  float totalSimEnergy3 = 0.;
  float totalSimEnergy4 = 0.;
  float totalSimEnergy5 = 0.;
  float totalSimEnergy6 = 0.;
  float totalSimEnergy7 = 0.;
  for (r1 = detID_to_simID_to_energy_map.begin(); r1!= detID_to_simID_to_energy_map.end(); r1++){
    std::map<int, float>::iterator r2;
    for (r2 = (*r1).second.begin(); r2!= (*r1).second.end(); r2++){
      std::cout<<(*r1).first<<" "<<(*r2).first<<" "<<(*r2).second<<"      layer? "<<rechits->at(detID_toRecHitIndex[(*r1).first]).depth()<<std::endl;
      totalSimEnergy += (*r2).second;
      if(rechits->at(detID_toRecHitIndex[(*r1).first]).depth() == 1) totalSimEnergy1 += (*r2).second;
      if(rechits->at(detID_toRecHitIndex[(*r1).first]).depth() == 2) totalSimEnergy2 += (*r2).second;
      if(rechits->at(detID_toRecHitIndex[(*r1).first]).depth() == 3) totalSimEnergy3 += (*r2).second;
      if(rechits->at(detID_toRecHitIndex[(*r1).first]).depth() == 4) totalSimEnergy4 += (*r2).second;
      if(rechits->at(detID_toRecHitIndex[(*r1).first]).depth() == 5) totalSimEnergy5 += (*r2).second;
      if(rechits->at(detID_toRecHitIndex[(*r1).first]).depth() == 6) totalSimEnergy6 += (*r2).second;
      if(rechits->at(detID_toRecHitIndex[(*r1).first]).depth() == 7) totalSimEnergy7 += (*r2).second;
    }
    //for(unsigned int j = 0; j < simID_to_rechitindex_map[(*r1).first].size(); j++){
    //for(unsigned int j = 0; j < detID_toRecHitIndex[(*r1).first].size(); j++){
    //  std::cout<<"    layer? "<<rechits->at(detID_toRecHitIndex[(*r1).first].at(j)).depth()<<std::endl;
    //}
  }

  std::cout<<"    in truthclustering.  Hcal totalsimenergy = "<<totalSimEnergy<<std::endl;
  std::cout<<"        totalsimenergy 1 = "<<totalSimEnergy1<<std::endl;
  std::cout<<"        totalsimenergy 2 = "<<totalSimEnergy2<<std::endl;
  std::cout<<"        totalsimenergy 3 = "<<totalSimEnergy3<<std::endl;
  std::cout<<"        totalsimenergy 4 = "<<totalSimEnergy4<<std::endl;
  std::cout<<"        totalsimenergy 5 = "<<totalSimEnergy5<<std::endl;
  std::cout<<"        totalsimenergy 6 = "<<totalSimEnergy6<<std::endl;
  std::cout<<"        totalsimenergy 7 = "<<totalSimEnergy7<<std::endl;


  reco::PFCluster temp;
  auto pfClusters = std::make_unique<reco::PFClusterCollection>();
  for(unsigned int i = 0; i < simEnergyStructVector.size(); i++){
    std::cout<<"starting cluster for "<<simEnergyStructVector.at(i).ID<<std::endl;
    temp.reset();
    for(unsigned int j = 0; j < simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].size(); j++){
      int seedInQuestion = simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j);
      if(used[seedInQuestion]){
	continue;
      }
      std::cout<<"     starting seed is "<<seedInQuestion<<std::endl;
      //auto ref = makeRefhit(rechits, seedInQuestion); //reco::PFRecHitRef(h, i);
      int alphaDepth = rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).depth();
      float alphaPhi = rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).positionREP().phi();
      float alphaEta = rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(j)).positionREP().eta();
      buildTopoCluster(rechits, mask, seedInQuestion, used, temp, simEnergyStructVector.at(i).energy, alphaDepth, seedInQuestion, simID_to_rechitindex_map[simEnergyStructVector.at(i).ID]);

      for(int k = 1; k < 8; k++){
	bool extended = false;
	for(unsigned int l = 0; l < simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].size(); l++){
	  int seedInQuestion2 = simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(l);
	  if(used[seedInQuestion2] || rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(l)).depth() != alphaDepth+k){
	    continue;
	  }
	  double deleta = (alphaEta - rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(l)).positionREP().eta());
	  double delphi, ogphi, newphi;
	  newphi = rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(l)).positionREP().phi();
	  ogphi = alphaPhi;
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
	  //if( deltaR > 0.2 ){
	  if( deltaR > 0.4 ){
	  //if( deltaR > 0.1 ){
	    continue;
	  }
	  extended = buildTopoCluster(rechits, mask, seedInQuestion2, used, temp, simEnergyStructVector.at(i).energy, alphaDepth, seedInQuestion, simID_to_rechitindex_map[simEnergyStructVector.at(i).ID]);
	  if(extended){
	    break;
	  }
	}
	if(!extended){
	  break;
	}
      }

      for(int k = 1; k < 8; k++){
	bool extended = false;
	for(unsigned int l = 0; l < simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].size(); l++){
	  int seedInQuestion2 = simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(l);
	  if(used[seedInQuestion2] || rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(l)).depth() != alphaDepth-k){
	    continue;
	  }
	  double deleta = (alphaEta - rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(l)).positionREP().eta());
	  double delphi, ogphi, newphi;
	  newphi = rechits->at(simID_to_rechitindex_map[simEnergyStructVector.at(i).ID].at(l)).positionREP().phi();
	  ogphi = alphaPhi;
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
	  //if( deltaR > 0.2 ){
	  if( deltaR > 0.4 ){
	  //if( deltaR > 0.1 ){
	    continue;
	  }
	  extended = buildTopoCluster(rechits, mask, seedInQuestion2, used, temp, simEnergyStructVector.at(i).energy, alphaDepth, seedInQuestion, simID_to_rechitindex_map[simEnergyStructVector.at(i).ID]);
	  if(extended){
	    break;
	  }
	}
	if(!extended){
	  break;
	}
      }

      if (!temp.recHitFractions().empty()){
	pfClusters->push_back(temp);
      }

      break;
    }
  }

  std::cout<<"finished building truth seeded clusters."<<std::endl;



  //pfClusters->insert(pfClusters->end(), pfClusters3HB->begin(), pfClusters3HB->end());




  _initialClustering->reset();
  if (_pfClusterBuilder)
    _pfClusterBuilder->reset();
  if (_pfClusterBuilderDepth)
    _pfClusterBuilderDepth->reset();

  _initialClustering->updateEvent(e);


  //edm::Handle<reco::PFRecHitCollection> rechits2;
  //auto rechits2 = std::make_unique<reco::PFRecHitCollection>();
  //for (unsigned int i = 0; i < rechits->size(); ++i) {
  //  //if (!rechitMask[i] || !seedable[i] || used[i])
  //  if (used[i])
  //    continue;
  //  //rechits2.emplace_back(rechits->at(i));
  //  //rechits2->push_back(rechits->at(i));
  //  rechits2->emplace_back(rechits->at(i));
  //}

  std::vector<bool> mask2(rechits->size(), true);
  for (const auto& cleaner : _cleaners) {
    cleaner->clean(rechits, mask2);
  }

  for (unsigned int i = 0; i < rechits->size(); ++i) {
    //if (!rechitMask[i] || !seedable[i] || used[i])
    //if ( used[i] ){ //I used this for my first attempt
    if (used[i] || rechits->at(i).energy() < 4.5){ //0.8 seems to have no effect relative to default.  I've tried 1.6.  100 should exclude basically everything
      mask2[i] = 0;
    }
  }

  //std::cout<<"what is mask2?"<<std::endl;
  //for(unsigned int i = 0; i < mask2.size(); i++){
  //  std::cout<<" "<<mask2.at(i);
  //}
  //std::cout<<"\n";
  
  // no seeding on these hits
  std::vector<bool> seedmask = mask2;
  for (const auto& cleaner : _seedcleaners) {
    cleaner->clean(rechits, seedmask);
  }

  //std::cout<<"what is seedmask?"<<std::endl;
  //for(unsigned int i = 0; i < seedmask.size(); i++){
  //  std::cout<<" "<<seedmask.at(i);
  //}
  //std::cout<<"\n";

  int sum_of_elems = 0;
  for(std::vector<bool>::iterator it = used.begin(); it != used.end(); ++it){
    if(*it){
      sum_of_elems += 1;
    }
  }
  std::cout<<"number of unused rechits: "<<sum_of_elems<<std::endl;

  //std::vector<bool> seedable = used;
  std::vector<bool> seedable(rechits->size(), false);
  _seedFinder->findSeeds(rechits, seedmask, seedable);

  auto initialClusters = std::make_unique<reco::PFClusterCollection>();
  _initialClustering->buildClusters(rechits, mask2, seedable, *initialClusters);
  LOGVERB("PFClusterProducer::produce()") << *_initialClustering;

  //auto initialClusters = std::make_unique<reco::PFClusterCollection>();
  //for (auto seed : seeds) {
  //  if (!mask[seed] || used[seed])
  //    continue;
  //  temp.reset();
  //  buildTopoCluster2(rechits, mask, seed, used, temp);
  //  if (!temp.recHitFractions().empty())
  //    initialClusters->push_back(temp);
  //}

  std::cout<<"finished building initial clusters. "<<initialClusters->size()<<std::endl;

  //std::vector<bool> seedable(rechits->size(), true);
  //std::vector<bool> seedable(rechits->size(), true);
  auto pfClusters2 = std::make_unique<reco::PFClusterCollection>();
  pfClusters2 = std::make_unique<reco::PFClusterCollection>();
  if (_pfClusterBuilder) {  // if we've defined a re-clustering step execute it
    //std::cout<<"i passed _pfClusterBuilder"<<std::endl;
    //std::cout<<"initialclusters size = "<<initialClusters->size()<<std::endl;
    _pfClusterBuilder->buildClusters(*initialClusters, seedable, *pfClusters2);
    LOGVERB("PFClusterProducer::produce()") << *_pfClusterBuilder;
  } else {
    pfClusters2->insert(pfClusters2->end(), initialClusters->begin(), initialClusters->end());
  }

  std::cout<<"finished splitting initial clusters. "<<pfClusters2->size()<<std::endl;

  std::vector<bool> seedable2;
  auto pfClusters3 = std::make_unique<reco::PFClusterCollection>();
  _pfClusterBuilderDepth->buildClusters(*pfClusters2, seedable2, *pfClusters3);
  //to insert generic clustered leftovers:

  std::cout<<"finished building nontruth seeded clusters. "<<pfClusters3->size()<<std::endl; 

 pfClusters->insert(pfClusters->end(), pfClusters3->begin(), pfClusters3->end());

  if (_positionReCalc) {
    _positionReCalc->calculateAndSetPositions(*pfClusters);
  }

  std::cout<<"finished position recalc."<<std::endl;

  if (_energyCorrector) {
    _energyCorrector->correctEnergies(*pfClusters);
  }

  //if (_prodInitClusters)
  //e.put(std::move(initialClusters), "initialClusters");
  std::cout<<"TRUTH clusters size = "<<pfClusters->size()<<std::endl;
  e.put(std::move(pfClusters));

}



bool PFTruthClusterProducer2::buildTopoCluster(const edm::Handle<reco::PFRecHitCollection>& input,
                                                     const std::vector<bool>& rechitMask,
                                                     unsigned int kcell,
                                                     std::vector<bool>& used,
                                                     reco::PFCluster& topocluster,
						     float simmedEnergy,
						     int alphaDepth,
						     int alpha,
						     std::vector<int> allowedRHs) {
  if(!( std::count(allowedRHs.begin(), allowedRHs.end(), kcell) )){
    return false;
  }
  bool foundSomething = false;
  auto const& cell = (*input)[kcell];
  int cell_layer = (int)cell.layer();
  if (cell_layer == PFLayer::HCAL_BARREL2 && std::abs(cell.positionREP().eta()) > 0.34) {
    cell_layer *= 100;
  }

  /*
  auto const& thresholds = _thresholds.find(cell_layer)->second;
  double thresholdE = 0.;
  double thresholdPT2 = 0.;

  for (unsigned int j = 0; j < (std::get<1>(thresholds)).size(); ++j) {
    int depth = std::get<0>(thresholds)[j];

    if ((cell_layer == PFLayer::HCAL_BARREL1 && cell.depth() == depth) ||
        (cell_layer == PFLayer::HCAL_ENDCAP && cell.depth() == depth) ||
        (cell_layer != PFLayer::HCAL_BARREL1 && cell_layer != PFLayer::HCAL_ENDCAP)) {
      thresholdE = std::get<1>(thresholds)[j];
      thresholdPT2 = std::get<2>(thresholds)[j];
    }
  }



  if (cell.energy() < thresholdE || cell.pt2() < thresholdPT2) {
    LOGDRESSED("GenericTopoCluster::buildTopoCluster()")
        << "RecHit " << cell.detId() << " with enegy " << cell.energy() << " GeV was rejected!." << std::endl;
    return false;
  }
  */

  auto k = kcell;
  used[k] = true;
  auto ref = makeRefhit(input, k);
  bool newcluster = topocluster.recHitFractions().empty();
  topocluster.addRecHitFraction(reco::PFRecHitFraction(ref, 1.0));

  if(newcluster){
    const auto rh_energy = ref->energy();
    topocluster.setSeed(ref->detId());
    topocluster.setEnergy(rh_energy);
    topocluster.setTime(ref->time());
    topocluster.setLayer(ref->layer());
    topocluster.setPosition(math::XYZPoint(ref->position().x(), ref->position().y(), ref->position().z()));
    topocluster.calculatePositionREP();
    topocluster.setDepth(ref->depth());
  } else{
    _allCellsPosCalc->calculateAndSetPosition(topocluster);
  }

  foundSomething = true;

  auto const& neighbours = (_useCornerCells ? cell.neighbours8() : cell.neighbours4());

  for (auto nb : neighbours) {
    //if (used[nb] || !rechitMask[nb]) {
    //if (used[nb] || calculateDR(input->at(alpha).positionREP().eta(), input->at(nb).positionREP().eta(), input->at(alpha).positionREP().phi(), input->at(nb).positionREP().phi()) > 0.3) {
    //if (used[nb] || calculateDR(input->at(alpha).positionREP().eta(), input->at(nb).positionREP().eta(), input->at(alpha).positionREP().phi(), input->at(nb).positionREP().phi()) > 0.2) {
    //if (used[nb] || calculateDR(input->at(alpha).positionREP().eta(), input->at(nb).positionREP().eta(), input->at(alpha).positionREP().phi(), input->at(nb).positionREP().phi()) > 0.4) {
    if (used[nb] || calculateDR(input->at(alpha).positionREP().eta(), input->at(nb).positionREP().eta(), input->at(alpha).positionREP().phi(), input->at(nb).positionREP().phi()) > 0.5) {
      LOGDRESSED("GenericTopoCluster::buildTopoCluster()")
          << "  RecHit " << cell.detId() << "\'s"
          << " neighbor RecHit " << input->at(nb).detId() << " with enegy " << input->at(nb).energy()
          << " GeV was rejected!"
          << " Reasons : " << used[nb] << " (used) " << !rechitMask[nb] << " (masked)." << std::endl;
      continue;
    }
    //buildTopoCluster(input, rechitMask, nb, used, topocluster);
    buildTopoCluster(input, rechitMask, nb, used, topocluster, simmedEnergy, alphaDepth, alpha, allowedRHs);
  }
  
  return foundSomething;
}


void PFTruthClusterProducer2::buildTopoCluster2(const edm::Handle<reco::PFRecHitCollection>& input,
                                                     const std::vector<bool>& rechitMask,
                                                     unsigned int kcell,
                                                     std::vector<bool>& used,
                                                     reco::PFCluster& topocluster) {
  auto const& cell = (*input)[kcell];
  int cell_layer = (int)cell.layer();
  if (cell_layer == PFLayer::HCAL_BARREL2 && std::abs(cell.positionREP().eta()) > 0.34) {
    cell_layer *= 100;
  }

  /*
  auto const& thresholds = _thresholds.find(cell_layer)->second;
  double thresholdE = 0.;
  double thresholdPT2 = 0.;

  for (unsigned int j = 0; j < (std::get<1>(thresholds)).size(); ++j) {
    int depth = std::get<0>(thresholds)[j];

    if ((cell_layer == PFLayer::HCAL_BARREL1 && cell.depth() == depth) ||
        (cell_layer == PFLayer::HCAL_ENDCAP && cell.depth() == depth) ||
        (cell_layer != PFLayer::HCAL_BARREL1 && cell_layer != PFLayer::HCAL_ENDCAP)) {
      thresholdE = std::get<1>(thresholds)[j];
      thresholdPT2 = std::get<2>(thresholds)[j];
    }
  }

  if (cell.energy() < thresholdE || cell.pt2() < thresholdPT2) {
    LOGDRESSED("GenericTopoCluster::buildTopoCluster()")
        << "RecHit " << cell.detId() << " with enegy " << cell.energy() << " GeV was rejected!." << std::endl;
    return;
  }
  */

  auto k = kcell;
  used[k] = true;
  auto ref = makeRefhit(input, k);
  bool newcluster = topocluster.recHitFractions().empty();
  topocluster.addRecHitFraction(reco::PFRecHitFraction(ref, 1.0));

  if(newcluster){
    const auto rh_energy = ref->energy();
    topocluster.setSeed(ref->detId());
    topocluster.setEnergy(rh_energy);
    topocluster.setTime(ref->time());
    topocluster.setLayer(ref->layer());
    topocluster.setPosition(math::XYZPoint(ref->position().x(), ref->position().y(), ref->position().z()));
    topocluster.calculatePositionREP();
    topocluster.setDepth(ref->depth());
  } else{
    _allCellsPosCalc->calculateAndSetPosition(topocluster);
  }
 
  auto const& neighbours = (_useCornerCells ? cell.neighbours8() : cell.neighbours4());

  for (auto nb : neighbours) {
    if (used[nb] || !rechitMask[nb]) {
      LOGDRESSED("GenericTopoCluster::buildTopoCluster()")
          << "  RecHit " << cell.detId() << "\'s"
          << " neighbor RecHit " << input->at(nb).detId() << " with enegy " << input->at(nb).energy()
          << " GeV was rejected!"
          << " Reasons : " << used[nb] << " (used) " << !rechitMask[nb] << " (masked)." << std::endl;
      continue;
    }
    buildTopoCluster2(input, rechitMask, nb, used, topocluster);
  }
}
