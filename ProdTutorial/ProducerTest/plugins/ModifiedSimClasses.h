#ifndef ModifiedSimClasses_H
#define ModifiedSimClasses_H


#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

class SimVertex_v2 : public SimVertex {
 public:
  typedef SimTrack ST;
  
 SimTrack_v2(const SimTrack& t, int it) : itrk(it){}
  
  int itrk;
  
}

//#include <iosfwd>
//#std::ostream& operator<<(std::ostream& o, const SimTrack& t);

#endif
