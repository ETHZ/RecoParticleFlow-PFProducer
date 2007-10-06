#ifndef RecoParticleFlow_PFProducer_TauBenchmarkAnalyzer
#define RecoParticleFlow_PFProducer_TauBenchmarkAnalyzer


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowReco/interface/PFSimParticleFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

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
  bool enablegenjets_;
  bool enablepfsimparticles_;
  std::string jetcollection_;
  edm::Handle<reco::PFSimParticleCollection> trueParticles_;
  edm::Handle<reco::CaloJetCollection> caloJets_;
  edm::Handle<reco::PFJetCollection> pfJets_;
  edm::Handle<reco::GenJetCollection> genJets_;
  edm::Handle<edm::HepMCProduct> hepMC_; 
  bool testflag_;
  int testcounter_;
  int total_;
  int VERBOSE;
  int verbosity_;
  TH1F* h_deltaETvisible_MCPFSIM_EHT_;
  TH1F* h_deltaETvisible_MCPFSIM_PF_;
  TH1F* h_deltaETvisible_MCGENJET_EHT_;
  TH1F* h_deltaETvisible_MCGENJET_PF_; 
  TH1F* h_deltaETvisible_MCHEPMC_PF_;	                      
  TH1F* h_deltaETvisible_MCHEPMC_EHT_;
  TH1F* h_deltaETvisible_GENJET_PFSIM_;                   
  TH1F* h_deltaETvisible_GENJET_HEPMC_;
  // multiplicities
  TH1F* h_genJetMatching_MCEHTnum_;
  TH1F* h_genJetMatching_MCPFnum_; 
  TH1F* h_genJetMatching_MCnum_;

  /// output root file
  TFile* file_;
};

#endif
