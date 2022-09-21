#ifndef DATAFORMATS_HCALRECHIT_HBHERECHITTRUTH_H
#define DATAFORMATS_HCALRECHIT_HBHERECHITTRUTH_H 1

#include <vector>

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"

#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"

/** \class HBHERecHitTruth
 *  
 * \author J. Mans - Minnesota
 */
class HBHERecHitTruth : public CaloRecHit {
public:
  typedef HcalDetId key_type;

  constexpr HBHERecHitTruth()
      : CaloRecHit(),
        timeFalling_(0),
        chiSquared_(-1),
        rawEnergy_(-1.0e21),
        auxEnergy_(-1.0e21),
        auxHBHE_(0),
        auxPhase1_(0),
        auxTDC_(0),
        parentPart_(0),
        posX_(0),
        posY_(0),
        posZ_(0),
        eta_(0),
        phi_(0),
        isHB_(false),
        depth_(0),
        isCleanHLT_(false),
        isCleanReco_(false),
        layer_(PFLayer::NONE),
        PFcluster0Id_(-1), PFcluster0frac_(-1), PFcluster1Id_(-1), PFcluster1frac_(-1), PFcluster2Id_(-1), PFcluster2frac_(-1),
        simHitEnergy_(0), simHitDepth_(0),simHitTime_(0) {}

  constexpr HBHERecHitTruth(const HBHERecHit& rh, uint32_t p = 0, float x = 0, float y = 0, float z = 0, float eta = 0 , float phi = 0, int depth = 0, PFLayer::Layer layer = PFLayer::NONE, bool cleanHLT = false, bool cleanReco = false, int id0 = -1, float fr0 = -1, int id1 = -1, float fr1 = -1, int id2 = -1, float fr2 = -1, float she = 0, int shd = -1, float sht = -99)
    : CaloRecHit(rh.id(), rh.energy(), rh.time()),
        timeFalling_(rh.timeFalling()),
        chiSquared_(rh.chi2()),
        rawEnergy_(rh.eraw()),
        auxEnergy_(rh.eaux()),
        auxHBHE_(rh.auxHBHE()),
        auxPhase1_(rh.auxPhase1()),
        auxTDC_(rh.auxTDC()),
        parentPart_(p),
        posX_(x),
        posY_(y),
        posZ_(z),
        eta_(eta),
        phi_(phi),
        isHB_(layer == PFLayer::HCAL_BARREL1 ? true : false),
        depth_(depth),
        isCleanHLT_(cleanHLT),
        isCleanReco_(cleanReco),
        layer_(layer),
        PFcluster0Id_(id0), PFcluster0frac_(fr0), PFcluster1Id_(id1), PFcluster1frac_(fr1), PFcluster2Id_(id2), PFcluster2frac_(fr2),
        simHitEnergy_(she), simHitDepth_(shd),simHitTime_(sht) {}

  /// get the hit falling time
  constexpr inline float timeFalling() const { return timeFalling_; }
  constexpr inline void setTimeFalling(float timeFalling) { timeFalling_ = timeFalling; }
  /// get the id
  constexpr inline HcalDetId id() const { return HcalDetId(detid()); }

  constexpr inline void setChiSquared(const float chi2) { chiSquared_ = chi2; }
  constexpr inline float chi2() const { return chiSquared_; }

  constexpr inline void setRawEnergy(const float en) { rawEnergy_ = en; }
  constexpr inline float eraw() const { return rawEnergy_; }

  constexpr inline void setAuxEnergy(const float en) { auxEnergy_ = en; }
  constexpr inline float eaux() const { return auxEnergy_; }

  constexpr inline void setAuxHBHE(const uint32_t aux) { auxHBHE_ = aux; }
  constexpr inline uint32_t auxHBHE() const { return auxHBHE_; }

  constexpr inline void setAuxPhase1(const uint32_t aux) { auxPhase1_ = aux; }
  constexpr inline uint32_t auxPhase1() const { return auxPhase1_; }

  constexpr inline void setAuxTDC(const uint32_t aux) { auxTDC_ = aux; }
  constexpr inline uint32_t auxTDC() const { return auxTDC_; }

  constexpr inline void setParentPart(const uint32_t parent) { parentPart_ = parent; }
  constexpr inline uint32_t parentPart() const { return parentPart_; }

  constexpr inline void setPosition(const float x, const float y, const float z, const float eta, const float phi) {
    posX_ = x;
    posY_ = y;
    posZ_ = z;
    eta_ = eta;
    phi_ = phi;
  }
  constexpr inline float x() const { return posX_; }
  constexpr inline float y() const { return posY_; }
  constexpr inline float z() const { return posZ_; }
  constexpr inline float eta() const { return eta_; }
  constexpr inline float phi() const { return phi_; }
  PFLayer::Layer layer() const { return layer_; }

  constexpr inline bool isHB() const { return isHB_; }
  constexpr inline bool isCleanHLT() const { return isCleanHLT_; }
  constexpr inline bool isCleanReco() const { return isCleanReco_; }
  constexpr inline int depth() const { return depth_; }

  constexpr inline void setClusters(const int id0, const float fr0, const int id1, const float fr1, const int id2, const float fr2) {
    PFcluster0Id_ = id0;
    PFcluster0frac_ = fr0;
    PFcluster1Id_ = id1;
    PFcluster1frac_ = fr1;
    PFcluster2Id_ = id2;
    PFcluster2frac_ = fr2;
  }
  constexpr inline int PFcluster0Id() const { return PFcluster0Id_; }
  constexpr inline float PFcluster0frac() const { return PFcluster0frac_; }
  constexpr inline int PFcluster1Id() const { return PFcluster1Id_; }
  constexpr inline float PFcluster1frac() const { return PFcluster1frac_; }
  constexpr inline int PFcluster2Id() const { return PFcluster2Id_; }
  constexpr inline float PFcluster2frac() const { return PFcluster2frac_; }

  constexpr inline void setSimHitEnergy(const float she) { simHitEnergy_ = she; }
  constexpr inline float simHitEnergy() const { return simHitEnergy_; }
  constexpr inline void setSimHitDepth(const float shd) { simHitDepth_ = shd; }
  constexpr inline int simHitDepth() const { return simHitDepth_; }
  constexpr inline void setSimHitTime(const float sht) { simHitTime_ = sht; }
  constexpr inline float simHitTime() const { return simHitTime_; }

  //void setParentPart(const uint32_t parent) { parentPart_ = parent; }
  //uint32_t parentPart() const { return parentPart_; }

  // The following method returns "true" for "Plan 1" merged rechits
  bool isMerged() const;

  // The following method fills the vector with the ids of the
  // rechits that have been merged to construct the "Plan 1" rechit.
  // For normal (i.e., not merged) rechits the vector will be cleared.
  void getMergedIds(std::vector<HcalDetId>* ids) const;

  // Returns the DetId of the front Id if it is a merged RecHit in "Plan 1"
  HcalDetId idFront() const;

private:
  float timeFalling_;
  float chiSquared_;
  float rawEnergy_;
  float auxEnergy_;
  uint32_t auxHBHE_;
  uint32_t auxPhase1_;
  uint32_t auxTDC_;
  uint32_t parentPart_;
  float posX_;
  float posY_;
  float posZ_;
  float eta_;
  float phi_;
  bool isHB_;
  int depth_;
  bool isCleanHLT_;
  bool isCleanReco_;
  PFLayer::Layer layer_ = PFLayer::NONE;
  int PFcluster0Id_;
  float PFcluster0frac_;
  int PFcluster1Id_;
  float PFcluster1frac_;
  int PFcluster2Id_;
  float PFcluster2frac_;
  float simHitEnergy_;
  int simHitDepth_;
  float simHitTime_;
};

std::ostream& operator<<(std::ostream& s, const HBHERecHitTruth& hit);

#endif
