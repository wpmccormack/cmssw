#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/BTauReco/interface/DeepBoostedJetTagInfo.h"

#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"

#include "HeterogeneousCore/SonicTriton/interface/TritonData.h"

#include "RecoBTag/FeatureTools/interface/deep_helpers.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <nlohmann/json.hpp>

using namespace btagbtvdeep;

class ParticleNetSonicJetTagsProducer : public TritonEDProducer<> {
public:
  explicit ParticleNetSonicJetTagsProducer(const edm::ParameterSet &);
  ~ParticleNetSonicJetTagsProducer() override;

  void acquire(edm::Event const &iEvent, edm::EventSetup const &iSetup, Input &iInput) override;
  void produce(edm::Event &iEvent, edm::EventSetup const &iSetup, Output const &iOutput) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &);

private:
  typedef std::vector<reco::DeepBoostedJetTagInfo> TagInfoCollection;
  typedef reco::JetTagCollection JetTagCollection;

  void make_inputs(const reco::DeepBoostedJetTagInfo &taginfo);

  const edm::EDGetTokenT<TagInfoCollection> src_;
  std::vector<std::string> flav_names_;             // names of the output scores
  std::vector<std::string> input_names_;            // names of each input group - the ordering is important!
  std::vector<std::vector<int64_t>> input_shapes_;  // shapes of each input group (-1 for dynamic axis)
  std::vector<unsigned> input_sizes_;               // total length of each input vector
  std::unordered_map<std::string, PreprocessParams> prep_info_map_;  // preprocessing info for each input group
  unsigned batchSize_;
  bool debug_ = false;
  int testint_;
};

ParticleNetSonicJetTagsProducer::ParticleNetSonicJetTagsProducer(const edm::ParameterSet &iConfig)
  : TritonEDProducer<>(iConfig, "ParticleNetSonicJetTagsProducer"),
    src_(consumes<TagInfoCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    flav_names_(iConfig.getParameter<std::vector<std::string>>("flav_names")),
    batchSize_(iConfig.getParameter<unsigned>("batchSize")),
    debug_(iConfig.getUntrackedParameter<bool>("debugMode", false)) {
  // load preprocessing info
  auto json_path = iConfig.getParameter<std::string>("preprocess_json");
  if (!json_path.empty()) {
    // use preprocessing json file if available
    std::ifstream ifs(edm::FileInPath(json_path).fullPath());
    nlohmann::json js = nlohmann::json::parse(ifs);
    js.at("input_names").get_to(input_names_);
    for (const auto &group_name : input_names_) {
      const auto &group_pset = js.at(group_name);
      auto &prep_params = prep_info_map_[group_name];
      group_pset.at("var_names").get_to(prep_params.var_names);
      if (group_pset.contains("var_length")) {
        prep_params.min_length = group_pset.at("var_length");
        prep_params.max_length = prep_params.min_length;
      } else {
        prep_params.min_length = group_pset.at("min_length");
        prep_params.max_length = group_pset.at("max_length");
        input_shapes_.push_back({1, (int64_t)prep_params.var_names.size(), -1});
      }
      const auto &var_info_pset = group_pset.at("var_infos");
      for (const auto &var_name : prep_params.var_names) {
        const auto &var_pset = var_info_pset.at(var_name);
        double median = var_pset.at("median");
        double norm_factor = var_pset.at("norm_factor");
        double replace_inf_value = var_pset.at("replace_inf_value");
        double lower_bound = var_pset.at("lower_bound");
        double upper_bound = var_pset.at("upper_bound");
        double pad = var_pset.contains("pad") ? double(var_pset.at("pad")) : 0;
        prep_params.var_info_map[var_name] =
            PreprocessParams::VarInfo(median, norm_factor, replace_inf_value, lower_bound, upper_bound, pad);
      }
    }
  } else {
    // otherwise use the PSet in the python config file
    const auto &prep_pset = iConfig.getParameterSet("preprocessParams");
    input_names_ = prep_pset.getParameter<std::vector<std::string>>("input_names");
    for (const auto &group_name : input_names_) {
      const auto &group_pset = prep_pset.getParameterSet(group_name);
      auto &prep_params = prep_info_map_[group_name];
      prep_params.var_names = group_pset.getParameter<std::vector<std::string>>("var_names");
      prep_params.min_length = group_pset.getParameter<unsigned>("var_length");
      prep_params.max_length = prep_params.min_length;
      const auto &var_info_pset = group_pset.getParameterSet("var_infos");
      for (const auto &var_name : prep_params.var_names) {
        const auto &var_pset = var_info_pset.getParameterSet(var_name);
        double median = var_pset.getParameter<double>("median");
        double norm_factor = var_pset.getParameter<double>("norm_factor");
        double replace_inf_value = var_pset.getParameter<double>("replace_inf_value");
        double lower_bound = var_pset.getParameter<double>("lower_bound");
        double upper_bound = var_pset.getParameter<double>("upper_bound");
        prep_params.var_info_map[var_name] =
            PreprocessParams::VarInfo(median, norm_factor, replace_inf_value, lower_bound, upper_bound, 0);
      }
    }
  }
  if (debug_) {
    for (unsigned i = 0; i < input_names_.size(); ++i) {
      const auto &group_name = input_names_.at(i);
      if (!input_shapes_.empty()) {
	std::cout << group_name << "\nshapes: ";
        for (const auto &x : input_shapes_.at(i)) {
	  std::cout << x << ", ";
        }
      }
      std::cout << "\nvariables: ";
      for (const auto &x : prep_info_map_.at(group_name).var_names) {
	std::cout << x << ", ";
      }
      std::cout << "\n";
    }
    std::cout << "flav_names: ";
    for (const auto &flav_name : flav_names_) {
      std::cout << flav_name << ", ";
    }
    std::cout << "\n";
  }

  // get output names from flav_names
  for (const auto &flav_name : flav_names_) {
    produces<JetTagCollection>(flav_name);
  }
  //preprocessInfoLoader(&iConfig);
}

ParticleNetSonicJetTagsProducer::~ParticleNetSonicJetTagsProducer() {}

void ParticleNetSonicJetTagsProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // pfDeepBoostedJetTags
  edm::ParameterSetDescription desc;
  TritonClient::fillPSetDescription(desc);
  desc.add<unsigned>("batchSize", 1);
  desc.add<edm::InputTag>("src", edm::InputTag("pfDeepBoostedJetTagInfos"));
  desc.add<std::string>("preprocess_json", "");
  // `preprocessParams` is deprecated -- use the preprocessing json file instead
  edm::ParameterSetDescription preprocessParams;
  preprocessParams.setAllowAnything();
  preprocessParams.setComment("`preprocessParams` is deprecated, please use `preprocess_json` instead.");
  desc.addOptional<edm::ParameterSetDescription>("preprocessParams", preprocessParams);
  desc.add<std::vector<std::string>>("flav_names",
                                     std::vector<std::string>{
				       "probTbcq",
                                         "probTbqq",
                                         "probTbc",
                                         "probTbq",
                                         "probWcq",
                                         "probWqq",
                                         "probZbb",
                                         "probZcc",
                                         "probZqq",
                                         "probHbb",
                                         "probHcc",
                                         "probHqqqq",
                                         "probQCDbb",
                                         "probQCDcc",
                                         "probQCDb",
                                         "probQCDc",
                                         "probQCDothers",
					 });
  desc.addOptionalUntracked<bool>("debugMode", false);

  descriptions.addWithDefaultLabel(desc);
}

void ParticleNetSonicJetTagsProducer::acquire(edm::Event const &iEvent, edm::EventSetup const &iSetup, Input &iInput) {
  edm::Handle<TagInfoCollection> tag_infos;
  iEvent.getByToken(src_, tag_infos);
  client_->setBatchSize(tag_infos->size());
  if (!tag_infos->empty()) {
    for (unsigned igroup = 0; igroup < input_names_.size(); ++igroup) {
      const auto &group_name = input_names_[igroup];
      auto &input = iInput.at(group_name);
      auto tdata = input.allocate<float>(true);
      for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n){
	const auto &taginfo = (*tag_infos)[jet_n];
	auto &vdata = (*tdata)[jet_n];
	const auto &prep_params = prep_info_map_.at(group_name);
	unsigned curr_pos = 0;
	// transform/pad
	for (unsigned i = 0; i < prep_params.var_names.size(); ++i) {
	  const auto &varname = prep_params.var_names[i];
	  const auto &raw_value = taginfo.features().get(varname);
	  const auto &info = prep_params.info(varname);
	  int insize = center_norm_pad(raw_value,
				     info.center,
				     info.norm_factor,
				     prep_params.min_length,
				     prep_params.max_length,
				     &vdata,
				     curr_pos,
				     info.pad,
				     info.replace_inf_value,
				     info.lower_bound,
				     info.upper_bound);
	  curr_pos += insize;
	  if (i == 0 && (!input_shapes_.empty())) {
	    input_shapes_[igroup][2] = insize;
	  }

	  if (debug_) {
	    std::cout << " -- var=" << varname << ", center=" << info.center << ", scale=" << info.norm_factor
		      << ", replace=" << info.replace_inf_value << ", pad=" << info.pad << std::endl;
	    for(unsigned i = curr_pos - insize; i < curr_pos; i++){
	      std::cout << vdata[i] << ",";
	    }
	    std::cout << std::endl;
	  }

	}
      }
      input.toServer(tdata);
    }
  }
}

void ParticleNetSonicJetTagsProducer::produce(edm::Event &iEvent,
                                              const edm::EventSetup &iSetup,
                                              Output const &iOutput) {
  edm::Handle<TagInfoCollection> tag_infos;
  iEvent.getByToken(src_, tag_infos);

  // initialize output collection
  std::vector<std::unique_ptr<JetTagCollection>> output_tags;
  if (!tag_infos->empty()) {
    auto jet_ref = tag_infos->begin()->jet();
    auto ref2prod = edm::makeRefToBaseProdFrom(jet_ref, iEvent);
    for (std::size_t i = 0; i < flav_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>(ref2prod));
    }
  } else {
    for (std::size_t i = 0; i < flav_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>());
    }
  }

  if(tag_infos->size() > 0){
    const auto &output1 = iOutput.begin()->second;
    const auto &outputs_from_server = output1.fromServer<float>();
    
    for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
      const auto &taginfo = (*tag_infos)[jet_n];
      
      if (!taginfo.features().empty()) {
	
	const auto &jet_ref = tag_infos->at(jet_n).jet();
	for (std::size_t flav_n = 0; flav_n < flav_names_.size(); flav_n++) {
	  (*(output_tags[flav_n]))[jet_ref] = outputs_from_server[jet_n][flav_n];
	}
      }
      else {
	const auto &jet_ref = tag_infos->at(jet_n).jet();
	for (std::size_t flav_n = 0; flav_n < flav_names_.size(); flav_n++) {
	  (*(output_tags[flav_n]))[jet_ref] = 0.;
	}
      }
    }
    
  }
  
  if (debug_) {
    std::cout << "=== " << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":" << iEvent.id().event()
              << " ===" << std::endl;
    for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
      const auto &jet_ref = tag_infos->at(jet_n).jet();
      std::cout << " - Jet #" << jet_n << ", pt=" << jet_ref->pt() << ", eta=" << jet_ref->eta()
                << ", phi=" << jet_ref->phi() << std::endl;
      for (std::size_t flav_n = 0; flav_n < flav_names_.size(); ++flav_n) {
	std::cout << "    " << flav_names_.at(flav_n) << " = " << (*(output_tags.at(flav_n)))[jet_ref] << std::endl;
      }
    }
  }

  // put into the event
  for (std::size_t flav_n = 0; flav_n < flav_names_.size(); ++flav_n) {
    iEvent.put(std::move(output_tags[flav_n]), flav_names_[flav_n]);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(ParticleNetSonicJetTagsProducer);
