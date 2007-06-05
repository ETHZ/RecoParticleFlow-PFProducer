#ifndef RecoParticleFlow_PFProducer_TauBenchmarkAnalyzer
#define RecoParticleFlow_PFProducer_TauBenchmarkAnalyzer


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowReco/interface/PFSimParticleFwd.h"
#include "DataFormats/JetReco/interface/BasicJetfwd.h"
#include "DataFormats/JetReco/interface/CaloJetfwd.h"
#include "DataFormats/JetReco/interface/PFJetfwd.h"
// #include <TH1F.h>

class TH1F;
class TFile;

class TauBenchmarkAnalyzer : public edm::EDAnalyzer {
 public:
  explicit TauBenchmarkAnalyzer(const edm::ParameterSet&);
  ~TauBenchmarkAnalyzer();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  //   virtual bool remove_leptonic_decay();
  //    virtual void makeJets();
  // ----------member data ---------------------------
  unsigned int minTracks_;
  edm::Handle<reco::PFSimParticleCollection> trueParticles_;
  edm::Handle<reco::CaloJetCollection> caloJets_;
  edm::Handle<reco::PFJetCollection> pfJets_;
  
  bool testflag_;
  int testcounter_;
  int total_;
  int VERBOSE;
  int verbosity_;
  TH1F* h_deltaETvisible_CMSSW_MCEHT_;
  TH1F* h_deltaETvisible_CMSSW_MCPF_;

  /// output root file
  TFile* file_;
};

#endif
