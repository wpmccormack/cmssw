#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterBuilderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/RecHitTopologicalCleanerBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/SeedFinderBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"

#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonData.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <memory>


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


class PFClusterSonicProducer : public TritonEDProducer<> {
  typedef RecHitTopologicalCleanerBase RHCB;
  typedef InitialClusteringStepBase ICSB;
  typedef PFClusterBuilderBase PFCBB;
  typedef PFCPositionCalculatorBase PosCalc;

public:
  PFClusterSonicProducer(const edm::ParameterSet&);
  ~PFClusterSonicProducer() override = default;

  void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
  void acquire(edm::Event const &iEvent, edm::EventSetup const &iSetup, Input &iInput) override;
  void produce(edm::Event &iEvent, edm::EventSetup const &iSetup, Output const &iOutput) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  //void update(const edm::EventSetup& es) override { _allCellsPosCalc->update(es); }
  //void produce(edm::Event&, const edm::EventSetup&) override;

private:
  // inputs
  edm::EDGetTokenT<reco::PFRecHitCollection> _rechitsLabel;
  //edm::EDGetTokenT<reco::PFClusterCollection> _rechitsLabel;
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

  int nhits2_;

protected:
  reco::PFRecHitRef makeRefhit(const edm::Handle<reco::PFRecHitCollection>& h, const unsigned i) const {
    return reco::PFRecHitRef(h, i);
  }
};

PFClusterSonicProducer::PFClusterSonicProducer(const edm::ParameterSet& conf)
  : TritonEDProducer<>(conf),
    _prodInitClusters(conf.getUntrackedParameter<bool>("prodInitialClusters", false)) {
  //_rechitsLabel = consumes<reco::PFRecHitCollection>(conf.getParameter<edm::InputTag>("recHitsSource"));
  _rechitsLabel = consumes<reco::PFRecHitCollection>(conf.getParameter<edm::InputTag>("pf_src"));
  //_rechitsLabel = consumes<reco::PFClusterCollection>(conf.getParameter<edm::InputTag>("pf_src"));
  //std::cout<<"in PFClusterSonicProducer::PFClusterSonicProducer, recHitsSource = "<<conf.getParameter<edm::InputTag>("pf_src")<<std::endl;
  

  edm::ConsumesCollector cc = consumesCollector();
  ////if (conf.exists("allCellsPositionCalc")) {
  const edm::ParameterSet& acConf = conf.getParameterSet("allCellsPositionCalc");
  const std::string& algoac = acConf.getParameter<std::string>("algoName");
  _allCellsPosCalc = PFCPositionCalculatorFactory::get()->create(algoac, acConf, cc);
  ////}

  //std::cout<<"HELLO I'M IN PFClusterSonicProducer::PFClusterSonicProducer"<<std::endl;
  /*
  //setup rechit cleaners
  const edm::VParameterSet& cleanerConfs = conf.getParameterSetVector("recHitCleaners");
  for (const auto& conf : cleanerConfs) {
    const std::string& cleanerName = conf.getParameter<std::string>("algoName");
    std::cout<<"SONIC cleanerName = "<<cleanerName<<std::endl;
    _cleaners.emplace_back(RecHitTopologicalCleanerFactory::get()->create(cleanerName, conf, cc));
  }

  if (conf.exists("seedCleaners")) {
    const edm::VParameterSet& seedcleanerConfs = conf.getParameterSetVector("seedCleaners");

    for (const auto& conf : seedcleanerConfs) {
      const std::string& seedcleanerName = conf.getParameter<std::string>("algoName");
      std::cout<<"SONIC seedcleanerName = "<<seedcleanerName<<std::endl;
      _seedcleaners.emplace_back(RecHitTopologicalCleanerFactory::get()->create(seedcleanerName, conf, cc));
    }
  }

  // setup seed finding
  const edm::ParameterSet& sfConf = conf.getParameterSet("seedFinder");
  const std::string& sfName = sfConf.getParameter<std::string>("algoName");
  std::cout<<"SONIC what is sfName? "<<sfName<<std::endl;
  _seedFinder = SeedFinderFactory::get()->create(sfName, sfConf);
  //setup topo cluster builder
  const edm::ParameterSet& initConf = conf.getParameterSet("initialClusteringStep");
  const std::string& initName = initConf.getParameter<std::string>("algoName");
  std::cout<<"SONIC what is initName? "<<initName<<std::endl;
  _initialClustering = InitialClusteringStepFactory::get()->create(initName, initConf, cc);
  //setup pf cluster builder if requested
  const edm::ParameterSet& pfcConf = conf.getParameterSet("pfClusterBuilder");
  if (!pfcConf.empty()) {
    const std::string& pfcName = pfcConf.getParameter<std::string>("algoName");
    //std::cout<<"SONIC what is pfcName? "<<pfcName<<" how about (*pfcName)? "<<(*pfcName)<<std::endl;
    std::cout<<"SONIC what is pfcName? "<<pfcName<<std::endl;
    _pfClusterBuilder = PFClusterBuilderFactory::get()->create(pfcName, pfcConf, cc);
  }
  //setup (possible) recalcuation of positions
  const edm::ParameterSet& pConf = conf.getParameterSet("positionReCalc");
  if (!pConf.empty()) {
    const std::string& pName = pConf.getParameter<std::string>("algoName");
    std::cout<<"SONIC what is pName? "<<pName<<std::endl;
    _positionReCalc = PFCPositionCalculatorFactory::get()->create(pName, pConf, cc);
  }
  // see if new need to apply corrections, setup if there.
  const edm::ParameterSet& cConf = conf.getParameterSet("energyCorrector");
  if (!cConf.empty()) {
    const std::string& cName = cConf.getParameter<std::string>("algoName");
    std::cout<<"SONIC what is cName? "<<cName<<std::endl;
    _energyCorrector = PFClusterEnergyCorrectorFactory::get()->create(cName, cConf);
  }
*/
  if (_prodInitClusters) {
    produces<reco::PFClusterCollection>("initialClusters");
  }
  produces<reco::PFClusterCollection>();
}

void PFClusterSonicProducer::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& es) {
  //std::cout<<"SONIC HELLO I'M IN PFClusterSonicProducer::beginLuminosityBlock"<<std::endl;
  //_initialClustering->update(es);
  //if (_pfClusterBuilder)
  //_pfClusterBuilder->update(es);
  //if (_positionReCalc)
  //_positionReCalc->update(es);
  //for (const auto& cleaner : _cleaners)
  //cleaner->update(es);
  //for (const auto& cleaner : _seedcleaners)
  //cleaner->update(es);
}

//void PFClusterSonicProducer::produce(edm::Event& e, const edm::EventSetup& es) {
void PFClusterSonicProducer::acquire(edm::Event const &iEvent, edm::EventSetup const &iSetup, Input &iInput) {
  client_->setBatchSize(1);
  //std::cout<<"SONIC HELLO I'M IN PFClusterSonicProducer::acquire"<<std::endl;
  //_initialClustering->reset();
  //if (_pfClusterBuilder)
  //_pfClusterBuilder->reset();

  edm::Handle<reco::PFRecHitCollection> rechits;
  //edm::Handle<reco::PFClusterCollection> rechits;
  iEvent.getByToken(_rechitsLabel, rechits);

  //_initialClustering->updateEvent(iEvent);

  //std::cout<<"SONIC rechits->size() = "<<rechits->size()<<std::endl;
  std::vector<bool> mask(rechits->size(), true);
  //for (const auto& cleaner : _cleaners) {
  //cleaner->clean(rechits, mask);
  //}

  // no seeding on these hits
  std::vector<bool> seedmask = mask;
  //for (const auto& cleaner : _seedcleaners) {
  //cleaner->clean(rechits, seedmask);
  //}

  auto nhits = rechits->size();
  nhits2_ = 0;
  for (unsigned i = 0; i < nhits; ++i) {
    if (!mask[i]){
      //std::cout<<"SONIC failed mask?"<<std::endl;
      continue;  // cannot seed masked objects
    }
    nhits2_++;
  }
  //std::cout<<"SONIC in PFClusterSonicProducer::acquire, what is nhits2_? = "<<nhits2_<<std::endl;
  if(nhits2_ < 1) return;
  auto &input = iInput.at("INPUT0");
  input.setShape(0, nhits2_);
  auto tdata = input.allocate<float>();
  auto& pfdata = (*tdata)[0];

  int h=0;
  for (unsigned i = 0; i < nhits; ++i) {
    if (!mask[i]){
      //std::cout<<"SONIC failed mask?"<<std::endl;
      continue;  // cannot seed masked objects
    }
    //auto& pfdata = (*tdata)[h];
    //std::cout<<"input rechit energy "<<rechits->at(i).energy()<<std::endl;
    pfdata.push_back(rechits->at(i).position().x());
    pfdata.push_back(rechits->at(i).position().y());
    pfdata.push_back(rechits->at(i).position().z());
    pfdata.push_back(rechits->at(i).energy());
    pfdata.push_back(rechits->at(i).positionREP().eta());
    pfdata.push_back(rechits->at(i).positionREP().phi());
    pfdata.push_back(rechits->at(i).time());
    //pfdata.push_back( sqrt(rechits->at(i).position().x()*rechits->at(i).position().x() + rechits->at(i).position().y()*rechits->at(i).position().y() ) );
    //pfdata.push_back( 2. * atan( exp(-1. * rechits->at(i).positionREP().eta() ) ) );
    pfdata.push_back( 2. * atan( exp(-1. * rechits->at(i).positionREP().eta() ) ) );
    pfdata.push_back( sqrt(rechits->at(i).position().x()*rechits->at(i).position().x() + rechits->at(i).position().y()*rechits->at(i).position().y() ) );
    pfdata.push_back(rechits->at(i).depth());
    //pfdata.push_back(rechits->at(i).time());
    h++;
  }

  //std::cout<<"about to send data to server"<<std::endl;

  input.toServer(tdata);

  //std::cout<<"sent data to server"<<std::endl;
  /*
  std::vector<bool> seedable(rechits->size(), false);
  _seedFinder->findSeeds(rechits, seedmask, seedable);

  auto initialClusters = std::make_unique<reco::PFClusterCollection>();
  _initialClustering->buildClusters(rechits, mask, seedable, *initialClusters);
  LOGVERB("PFClusterSonicProducer::produce()") << *_initialClustering;

  auto pfClusters = std::make_unique<reco::PFClusterCollection>();
  pfClusters = std::make_unique<reco::PFClusterCollection>();
  std::cout<<"SONIC first check of pfClusters.size() = "<<pfClusters->size()<<std::endl;
  if (_pfClusterBuilder) {  // if we've defined a re-clustering step execute it
    _pfClusterBuilder->buildClusters(*initialClusters, seedable, *pfClusters);
    LOGVERB("PFClusterSonicProducer::produce()") << *_pfClusterBuilder;
  } else {
    pfClusters->insert(pfClusters->end(), initialClusters->begin(), initialClusters->end());
  }
  std::cout<<"SONIC final check of pfClusters.size() = "<<pfClusters->size()<<std::endl;

  if (_positionReCalc) {
    _positionReCalc->calculateAndSetPositions(*pfClusters);
  }

  if (_energyCorrector) {
    _energyCorrector->correctEnergies(*pfClusters);
  }

  if (_prodInitClusters)
    e.put(std::move(initialClusters), "initialClusters");
  e.put(std::move(pfClusters));
  */
}


void PFClusterSonicProducer::produce(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) {
  //std::cout<<"in PFClusterSonicProducer::produce.  nhits2_ = "<<nhits2_<<std::endl;
  edm::Handle<reco::PFRecHitCollection> rechits;
  iEvent.getByToken(_rechitsLabel, rechits);

  if(nhits2_ < 1) return;
  const auto& output1 = iOutput.begin()->second;
  const auto& outputs = output1.fromServer<float>();

  //std::map<int, reco::PFCluster> cluster_map; //map keys = cluster index from SPVCNN and values = clusters that may be added to
  std::map<int, std::vector<reco::PFCluster>> cluster_map; //map keys = cluster index from SPVCNN and values = clusters that may be added to

  auto pfClusters = std::make_unique<reco::PFClusterCollection>();
  for (int i = 0; i < nhits2_; ++i) {
    //std::cout<<"SONIC outputs[0]["<<i<<"] = "<<outputs[0][i]<<" while x = "<<rechits->at(i).position().x()<<" y = "<<rechits->at(i).position().y()<<" z = "<<rechits->at(i).position().z()<<" energy = "<<rechits->at(i).energy()<<" eta = "<<rechits->at(i).positionREP().eta()<<" phi = "<<rechits->at(i).positionREP().phi()<<" and time = "<<rechits->at(i).time()<<std::endl;

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

    //std::cout<<"cluster got seed = "<<ref->detId()<<" energy = "<<rh_energy<<" time = "<<ref->time()<<" layer = "<<ref->layer()<<" position "<<ref->position().x()<<", "<<ref->position().y()<<", "<<ref->position().z()<<" depth = "<<ref->depth()<<std::endl;

    if( !cluster_map.count( outputs[0][i] ) ){ //new cluster index
      cluster_map[outputs[0][i]].push_back(current);
      //std::cout<<"NEW cluster "<<outputs[0][i]<<std::endl;
    } else{ //code for combining clusters taken from PFMultiDepthClusterizer
      //std::cout<<"EXISTING cluster "<<outputs[0][i]<<std::endl;


      bool newCluster = true;
      unsigned int theClusterIndex = 0;
      for(unsigned int cs = 0; cs < cluster_map[outputs[0][i]].size(); cs++){
        double deleta = (cluster_map[outputs[0][i]].at(cs).eta() - current.eta());
        double delphi, ogphi, newphi;
        newphi = current.phi();
        ogphi = cluster_map[outputs[0][i]].at(cs).phi();
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
	//std::cout<<"delta eta = "<<deleta<<" and delta phi = "<<delphi<<std::endl;
        double deltaR = sqrt( deleta*deleta + delphi*delphi );
	//std::cout<<"delta eta = "<<deleta<<" and delta phi = "<<delphi<<" so deltaR = "<<deltaR<<std::endl;
        if(deltaR < 1.){
          newCluster=false;
          theClusterIndex = cs;
          break;
        }
      }

      if(newCluster){
        cluster_map[outputs[0][i]].push_back(current);
	//std::cout<<"ANOTHER NEW cluster "<<outputs[0][i]<<std::endl;
      }

      else{
        double e1 = 0.0;
        double e2 = 0.0;
	reco::PFCluster main = cluster_map[outputs[0][i]].at(theClusterIndex);
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
        _allCellsPosCalc->calculateAndSetPosition(main);
	cluster_map[outputs[0][i]].at(theClusterIndex) = main;
      }

    }

    /*
      double e1 = 0.0;
      double e2 = 0.0;
      reco::PFCluster main = cluster_map[outputs[0][i]];
      for (const auto& fraction : main.recHitFractions()){
	if (fraction.recHitRef()->detId() == main.seed()) {
	  e1 = fraction.recHitRef()->energy();
	}
      }
      //std::cout<<"cluster "<<outputs[0][i]<<" seed energy = "<<e1<<std::endl;

      for (const auto& fraction : current.recHitFractions()) {
	main.addRecHitFraction(fraction);
	if (fraction.recHitRef()->detId() == current.seed()) {
	  e2 = fraction.recHitRef()->energy();
	}
      }
      //std::cout<<"the current e2 = "<<e2<<std::endl;
      if (e2 > e1){
	main.setSeed(current.seed());
      }
      //std::cout<<"main energy = "<<main.energy()<<std::endl;
      cluster_map[outputs[0][i]] = main;
      */
  }

    //pfClusters->push_back(current);
    //pfClusters->push_back(&rechits->at(i));
  //}
  
  //for (it = cluster_map.begin(); it != cluster_map.end(); it++){
  //for (auto const& cl : cluster_map){
  //std::cout<<"soniccluster energy = "<<cl.second.energy()<<std::endl;
  //_allCellsPosCalc->calculateAndSetPosition(cluster_map[cl.first]);
  //  //std::cout<<"did i recalibrate the energy? energy = "<<cl.second.energy()<<std::endl;
  //  //std::cout<<"adding a cluster at layer "<<cl.second.layer()<<" with energy "<<cl.second.energy()<<std::endl;
  //  if(cl.second.energy() > 0){
  //    pfClusters->push_back(cl.second);
  //  }
  //}
  for (auto const& cl : cluster_map){
    for(unsigned int cs2 = 0; cs2 < cluster_map[cl.first].size(); cs2++){
      _allCellsPosCalc->calculateAndSetPosition(cluster_map[cl.first].at(cs2));
      if(cl.second.at(cs2).energy() > 0){
	pfClusters->push_back(cl.second.at(cs2));
      }
    }
  }

  //std::cout<<"SONIC sonicclusters size = "<<pfClusters->size()<<std::endl;
  iEvent.put(std::move(pfClusters));

  // outputs are px and py
  //float px = outputs[0][0] * norm_;
  //float py = outputs[0][1] * norm_;

  // subtract the lepton pt contribution
  //px -= px_leptons_;
  //py -= py_leptons_;

  //if (debug_) {
  //  std::cout << "MET from DeepMET Sonic Producer is MET_x " << px << " and MET_y " << py << std::endl;
  //}

  //auto pf_mets = std::make_unique<pat::METCollection>();
  //const reco::Candidate::LorentzVector p4(px, py, 0., std::hypot(px, py));
  //pf_mets->emplace_back(reco::MET(p4, {}));
  //iEvent.put(std::move(pf_mets));
}


void PFClusterSonicProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  TritonClient::fillPSetDescription(desc);
  desc.add<edm::InputTag>("pf_src");
  
  //desc.add<edm::ParameterSet>("allCellsPositionCalc");
  //desc.add<edm::ParameterSet>("allCellsPositionCalc");

  edm::ParameterSetDescription reClusterizer;
  reClusterizer.add<double>("minAllowedNormalization", 1.0E-9);
  reClusterizer.add<int>("posCalcNCrystals", -1);
  reClusterizer.add<std::string>("algoName", "Basic2DGenericPFlowPositionCalc");
  reClusterizer.add<double>("minFractionInCalc", 1.0E-9);
  /*
  std::vector<edm::ParameterSetDescription> logWeightDenominatorByDetector;
  edm::ParameterSetDescription det1;
  det1.add<std::vector<int>>("depths", {1, 2, 3, 4});
  det1.add<std::string>("detector", "HCAL_BARREL1");
  det1.add<std::vector<double>>("logWeightDenominator", {0.1, 0.2, 0.3, 0.3});
  logWeightDenominatorByDetector.push_back(det1);
  edm::ParameterSetDescription det2;
  det2.add<std::vector<int>>("depths", {1, 2, 3, 4, 5, 6, 7});
  det2.add<std::string>("detector", "HCAL_ENDCAP");
  det2.add<std::vector<double>>("logWeightDenominator", {0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2});
  logWeightDenominatorByDetector.push_back(det2);

  reClusterizer.add<std::vector<edm::ParameterSetDescription>>("logWeightDenominatorByDetector", logWeightDenominatorByDetector);
  */
  /*
  std::vector<edm::ParameterSetDescription> logWeightDenominatorByDetector;
  edm::ParameterSetDescription det1;
  //det1.add<std::vector<int>>("depths", {1, 2, 3, 4});
  det1.add<std::string>("detector", "HCAL_BARREL1");
  //det1.add<std::vector<double>>("logWeightDenominator", {0.1, 0.2, 0.3, 0.3});
  logWeightDenominatorByDetector.push_back(det1);
  edm::ParameterSetDescription det2;
  //det2.add<std::vector<int>>("depths", {1, 2, 3, 4, 5, 6, 7});
  det2.add<std::string>("detector", "HCAL_ENDCAP");
  //det2.add<std::vector<double>>("logWeightDenominator", {0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2});
  logWeightDenominatorByDetector.push_back(det2);

  //reClusterizer.add<std::vector<edm::ParameterSetDescription>>("logWeightDenominatorByDetector", logWeightDenominatorByDetector);
  reClusterizer.addVPSet("logWeightDenominatorByDetector", logWeightDenominatorByDetector);
  */
 

  edm::ParameterSetDescription logWeightDenominatorByDetector;
  logWeightDenominatorByDetector.add<std::vector<int>>("depths", {1, 2, 3, 4});
  logWeightDenominatorByDetector.add<std::string>("detector", "HCAL_BARREL1");
  logWeightDenominatorByDetector.add<std::vector<double>>("logWeightDenominator", {0.1, 0.2, 0.3, 0.3});
  std::vector<edm::ParameterSet> vDefaults;
  edm::ParameterSet vDefaults0;
  vDefaults0.addParameter<std::vector<int>>("depths", {1, 2, 3, 4});
  vDefaults0.addParameter<std::string>("detector", "HCAL_BARREL1");
  vDefaults0.addParameter<std::vector<double>>("logWeightDenominator", {0.1, 0.2, 0.3, 0.3});
  vDefaults.push_back(vDefaults0);
  edm::ParameterSet vDefaults1;
  vDefaults1.addParameter<std::vector<int>>("depths", {1, 2, 3, 4, 5, 6, 7});
  vDefaults1.addParameter<std::string>("detector", "HCAL_ENDCAP");
  vDefaults1.addParameter<std::vector<double>>("logWeightDenominator", {0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2});
  vDefaults.push_back(vDefaults1);

  reClusterizer.addVPSet("logWeightDenominatorByDetector", logWeightDenominatorByDetector, vDefaults);
  desc.add<edm::ParameterSetDescription>("allCellsPositionCalc", reClusterizer);
  /*
  edm::ParameterSetDescription desc;
  edm::ParameterSetDescription validator;
  validator.add<int>("x", 7);

  std::vector<edm::ParameterSet> vDefaults;
  edm::ParameterSet vDefaults0;
  vDefaults.push_back(vDefaults0);
  edm::ParameterSet vDefaults1;
  vDefaults1.addParameter<int>("x", 100);
  vDefaults.push_back(vDefaults1);

  desc.addVPSet("nameForVPSet", validator, vDefaults);

  desc.add<edm::ParameterSetDescription>("allCellsPositionCalc", reClusterizer);
  */
  /*
    allCellsPositionCalc = cms.PSet(
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        logWeightDenominatorByDetector = cms.VPSet(
          cms.PSet(  depths = cms.vint32( 1, 2, 3, 4 ),
            detector = cms.string( "HCAL_BARREL1" ),
            logWeightDenominator = cms.vdouble( 0.1, 0.2, 0.3, 0.3 )
          ),
          cms.PSet(  depths = cms.vint32( 1, 2, 3, 4, 5, 6, 7 ),
            detector = cms.string( "HCAL_ENDCAP" ),
            logWeightDenominator = cms.vdouble( 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 )
          )
        )
      )
   */



  //desc.add<bool>("ignore_leptons", false);
  //desc.add<double>("norm_factor", 50.);
  //desc.add<unsigned int>("max_n_pf", 4500);
  //desc.addOptionalUntracked<bool>("debugMode", false);
  descriptions.add("PFClusterSonicProducer", desc);
}

DEFINE_FWK_MODULE(PFClusterSonicProducer);
