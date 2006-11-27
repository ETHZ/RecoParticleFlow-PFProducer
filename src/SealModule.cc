#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "RecoParticleFlow/PFProducer/interface/PFProducer.h"
#include "RecoParticleFlow/PFProducer/interface/TauHadronDecayFilter.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(PFProducer);
DEFINE_ANOTHER_FWK_MODULE(TauHadronDecayFilter);