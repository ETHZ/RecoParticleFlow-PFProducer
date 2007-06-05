// system include files
#include <memory>

// user include files

#include "RecoParticleFlow/PFProducer/test/TauBenchmarkAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrajectoryPoint.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrajectoryPointFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/LorentzVectorFwd.h"
#include "DataFormats/Math/interface/Vector3D.h"


#include "RecoParticleFlow/PFRootEvent/interface/EventColin.h" 

#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TLorentzVector.h>

// #include "RecoParticleFlow/PFRootEvent/interface/IO.h"
// #include "RecoParticleFlow/PFRootEvent/interface/Utils.h" 


using namespace edm;
using namespace std;
using namespace math;

//define this as a plug-in
DEFINE_FWK_MODULE(TauBenchmarkAnalyzer);

TauBenchmarkAnalyzer::TauBenchmarkAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  testcounter_=0;
  total_=0;    

  string outputRootFileName 
    = iConfig.getUntrackedParameter< string >("outputRootFileName", 
					      "tauBenchmark.root");

  LogDebug("TauBenchmarkAnalyzer")
    <<"opening output root file "<<outputRootFileName<<endl;
  
  file_ = TFile::Open(outputRootFileName.c_str(), "RECREATE");
 
  if(file_->IsZombie() ) {
    string err = "output file ";
    err += outputRootFileName;
    err += " can't be opened";
    throw cms::Exception("OutputFileOpen",err);
  }

  h_deltaETvisible_CMSSW_MCEHT_
    = new TH1F("h_deltaETvisible_MCEHT","Jet Et difference CaloTowers-MC"
	       ,500,-50,50);
  
  h_deltaETvisible_CMSSW_MCPF_
    = new TH1F("h_deltaETvisible_MCPF" ,
	       "Jet Et difference ParticleFlow-MC"
	       ,500,-50,50);

}


TauBenchmarkAnalyzer::~TauBenchmarkAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void TauBenchmarkAnalyzer::beginJob(const edm::EventSetup&)
{
  std::cout<<"Starting MyAnalyzer"<<std::endl;
  verbosity_=VERBOSE=1;
}


void TauBenchmarkAnalyzer::endJob() {
  // std::cout<<"Testcounter: "<<testcounter_<<std::endl;
  // std::cout<<"Total: "<<total_<<std::endl; 
  //   TCanvas c1("c1","delta et",200,200,800,600);
  //   h_deltaETvisible_CMSSW_MCPF_->Draw("");
  //   h_deltaETvisible_CMSSW_MCEHT_->Draw("Same");
  //   c1.Print("test.jpg");
//   h_deltaETvisible_CMSSW_MCEHT_->Write();
//   h_deltaETvisible_CMSSW_MCPF_->Write();
  file_->Write();
  file_->Close();
}





void TauBenchmarkAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  cout<<endl<<endl<<endl;
  //this is the event loop
  double total_et=0.0;
  //  double total_et_calo=0.0;
  trueParticles_.clear();
  caloJets_.clear();
  pfJets_.clear();
  //recoPFSimParticles_particleFlow__USER.obj
  
  iEvent.getByLabel("particleFlowSimParticle", trueParticles_);
  iEvent.getByLabel("Fastjet10CaloJets", caloJets_);
  iEvent.getByLabel("Fastjet10PFJets", pfJets_);

  //========================================================================
  //REMOVE EVENTS LEPTONIC DECAY OF THE TAUS================================
  //========================================================================

  
  vector<reco::PFSimParticle> vectPART;
  for ( unsigned i=0;  i < (*trueParticles_).size(); i++) {
    const reco::PFSimParticle& ptc = (*trueParticles_)[i];
    vectPART.push_back(ptc);
  }//loop 
  for ( unsigned i=0;  i < (*trueParticles_).size(); i++) {
    const reco::PFSimParticle& ptc = (*trueParticles_)[i];
    const std::vector<int>& ptcdaughters = ptc.daughterIds();
    if (abs(ptc.pdgCode()) == 15) {
      for ( unsigned int dapt=0; dapt < ptcdaughters.size(); ++dapt) {
	unsigned int pdgdaugter = abs(vectPART[ptcdaughters[dapt]].pdgCode());
	cout << "DAUGHTER=" << pdgdaugter << endl;
	if (pdgdaugter == 11 || pdgdaugter == 13) { 
	  if (verbosity_ == VERBOSE) cout << "REMOVED" << endl; 
	  return;
	  
	}//electron or muons?
      }//loop daughter
    }//tau
  }//loop particle
  



  
  // ====================================================================
  //MAKING TRUE PARTICLE JETS============================================
  // ====================================================================
  
  TLorentzVector partTOTMC;
  //ctor<reco::PFSimParticle> vectPART;
  for ( unsigned i=0;  i < (*trueParticles_).size(); i++) {
    const reco::PFSimParticle& ptc = (*trueParticles_)[i];
    vectPART.push_back(ptc);
  }//loop
  
  for ( unsigned i=0;  i < (*trueParticles_).size(); i++) 
    {
      const reco::PFSimParticle& ptc = (*trueParticles_)[i];
      const std::vector<int>& ptcdaughters = ptc.daughterIds();
      
      if (abs(ptc.pdgCode()) == 15) {
	for ( unsigned int dapt=0; dapt < ptcdaughters.size(); ++dapt) {
	  
	  const reco::PFTrajectoryPoint& tpatvtx 
	    = vectPART[ptcdaughters[dapt]].trajectoryPoint(0);
	  TLorentzVector partMC;
	  partMC.SetPxPyPzE(tpatvtx.momentum().Px(),
			    tpatvtx.momentum().Py(),
			    tpatvtx.momentum().Pz(),
			    tpatvtx.momentum().E());
	  partTOTMC += partMC;
	  int pdgcode = vectPART[ptcdaughters[dapt]].pdgCode();
	  cout << "pdgcode:"<<pdgcode << endl;
	  //cout << tpatvtx << endl;
	  //cout << partMC.Px() << " " << partMC.Py() << " " 
	  //    << partMC.Pz() << " " << partMC.E()
	  //     << " PT=" 
	  //     << sqrt(partMC.Px()*partMC.Px()+partMC.Py()*partMC.Py()) 
	  //     << endl;
	}//loop daughter
      }//tau?
    }//loop particles
  total_et=partTOTMC.Et();
  
  
  
  cout<<"total et: "<<total_et<<endl;
  


  // ==================================================================
  //CALO TOWER JETS (ECAL+HCAL Towers)=================================
  // ==================================================================

  double JetEHTETmax=0.0;
  double jetcalo_et=0.0;
  for ( unsigned int i = 0; i < (*caloJets_).size(); ++i)
    {
      jetcalo_et=(*caloJets_)[i].et();
      //if (jetcalo_et >= JetEHTETmax) 
      JetEHTETmax += jetcalo_et;
    }//loop calo towers
  
  cout<<"total et calo: "<<JetEHTETmax<<endl;
  
  
  // ==================================================================
  //PF Jets ===========================================================
  // ==================================================================
  
  double pfjet_et=0.0;
  double JetPFETmax=0.0;
  for ( unsigned int i = 0; i < (*pfJets_).size(); ++i)
    {
      pfjet_et=(*pfJets_)[i].et();
      JetPFETmax += pfjet_et;
    }//loop calo towers
  cout<<"total et pfjet: "<<JetPFETmax<<endl;
  
  // ==================================================================
  // Status output ====================================================
  // ==================================================================
  
  
  cout<<"delta et normal: "<<JetEHTETmax-total_et<<endl;
  cout<<"delta et pflow : "<<JetPFETmax-total_et<<endl; 
  
  // ==================================================================
  //Filling Histogramms ===============================================
  // ==================================================================

  h_deltaETvisible_CMSSW_MCPF_->Fill(JetPFETmax-total_et);
  h_deltaETvisible_CMSSW_MCEHT_->Fill(JetEHTETmax-total_et); 
  
}











