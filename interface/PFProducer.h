#ifndef RecoParticleFlow_PFProducer_h_
#define RecoParticleFlow_PFProducer_h_

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

#include "RecoTracker/TrackProducer/interface/TrackProducerAlgorithm.h"


/**\class PFProducer 
\brief Producer for particle flow tracks, particles and reconstructed particles 

This producer makes use of PFAlgo, the particle flow algorithm.

\author Colin Bernet, Renaud Bruneliere
\date   July 2006
*/


class PFProducer : public edm::EDProducer {
 public:
  explicit PFProducer(const edm::ParameterSet&);
  ~PFProducer();
  
  virtual void produce(edm::Event&, const edm::EventSetup&);

 private:
  // ----------member data ---------------------------
  std::string tcCollection_;
  std::string pfRecTrackCollection_;
  edm::ParameterSet conf_;
  TrackProducerAlgorithm trackAlgo_;
};

#endif
