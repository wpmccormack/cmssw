// -*- C++ -*-
//
// Package:    ProdTutorial/ProducerTest
// Class:      HLTJetMatcher
//
/**\class HLTJetMatcher HLTJetMatcher.cc ProdTutorial/ProducerTest/plugins/HLTJetMatcher.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  William McCormack
//         Created:  Tues, 14 June 2022 21:10:34 GMT
//
//


// system include files
#include <algorithm>
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

#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/transform.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

//
// class declaration
//

class HLTJetMatcher : public edm::stream::EDProducer<> {
public:
explicit HLTJetMatcher(const edm::ParameterSet&);
~HLTJetMatcher() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  edm::EDGetTokenT<edm::View<reco::PFJet>> jetsToken_;
  bool addGenJetMatch_;
  bool embedGenJetMatch_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetToken_;

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
HLTJetMatcher::HLTJetMatcher(const edm::ParameterSet& iConfig){
//register your products
/* Examples
   produces<ExampleData2>();
   
   //if do put with a label
   produces<ExampleData2>("label");
   
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/

jetsToken_ = consumes<edm::View<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jetSource"));

addGenJetMatch_ = iConfig.getParameter<bool>("addGenJetMatch");
embedGenJetMatch_ = iConfig.getParameter<bool>("embedGenJetMatch");
if (addGenJetMatch_){
genJetToken_ = consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("genJetMatch"));
}

produces<reco::GenJetCollection>("genJets");
}



HLTJetMatcher::~HLTJetMatcher() {
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
void HLTJetMatcher::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

// Get the vector of jets
edm::Handle<edm::View<reco::PFJet>> jets;
iEvent.getByToken(jetsToken_, jets);

 edm::Handle<edm::Association<reco::GenJetCollection>> genJetMatch;
 if (addGenJetMatch_){
   iEvent.getByToken(genJetToken_, genJetMatch);
 }

auto genJetsOut = std::make_unique<reco::GenJetCollection>();

edm::RefProd<reco::GenJetCollection> h_genJetsOut = iEvent.getRefBeforePut<reco::GenJetCollection>("genJets");
 
 for(unsigned int j = 0; j < jets->size(); j++){
   //for (edm::View<reco::PFJet>::const_iterator itJet = jets->begin(); itJet != jets->end(); itJet++) {
   // construct the Jet from the ref -> save ref to original object
   //unsigned int idx = itJet - jets->begin();
   //edm::RefToBase<reco::PFJet> jetRef = jets->refAt(idx);
   edm::RefToBase<reco::PFJet> jetRef = jets->refAt(j);
   //edm::Ptr<reco::Jet> jetPtr = jets->ptrAt(idx);
   reco::PFJet test = jets->at(j);
   //reco::PFJet ajet(jetRef);
   //Jet ajet(jetRef);
   
   // store the match to the GenJets
   if (addGenJetMatch_) {
     reco::GenJetRef genjet = (*genJetMatch)[jetRef];
     if (genjet.isNonnull() && genjet.isAvailable()) {
       genJetsOut->push_back(*genjet);
       // set the "forward" ref to the thinned collection
       edm::Ref<reco::GenJetCollection> genForwardRef(h_genJetsOut, genJetsOut->size() - 1);
       // set the "backward" ref to the original collection
       const edm::Ref<reco::GenJetCollection> &genBackRef(genjet);
       // make the FwdPtr
       edm::FwdRef<reco::GenJetCollection> genjetFwdRef(genForwardRef, genBackRef);
       //ajet.setGenJetRef(genjetFwdRef);
       //jetRef->setGenJetRef(genjetFwdRef);
       //test.setGenJetRef(genjetFwdRef);
     }  // leave empty if no match found
   }
   
 }

 iEvent.put(std::move(genJetsOut), "genJets");

}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void HLTJetMatcher::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void HLTJetMatcher::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
HLTJetMatcher::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
HLTJetMatcher::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HLTJetMatcher::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HLTJetMatcher::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HLTJetMatcher::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  desc.add<edm::InputTag>("jetSource", edm::InputTag("no default"))->setComment("input collection");
  desc.add<bool>("addGenJetMatch", true)->setComment("add MC matching");
  desc.add<bool>("embedGenJetMatch", false)->setComment("embed MC matched MC information");
  desc.add<edm::InputTag>("genJetMatch", edm::InputTag())->setComment("input with MC match information");

}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTJetMatcher);
