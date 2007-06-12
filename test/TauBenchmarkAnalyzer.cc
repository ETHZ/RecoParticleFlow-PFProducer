// system include files
#include <memory>

// user include files

#include "RecoParticleFlow/PFProducer/test/TauBenchmarkAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"

#include "DataFormats/ParticleFlowReco/interface/PFTrajectoryPoint.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrajectoryPointFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/HepMCCandidate.h"


#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"
#include "RecoParticleFlow/PFClusterAlgo/interface/PFClusterAlgo.h"
#include "RecoParticleFlow/PFAlgo/interface/PFBlock.h"
#include "RecoParticleFlow/PFAlgo/interface/PFBlockElement.h"
#include "RecoParticleFlow/PFAlgo/interface/PFGeometry.h"
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
	enablegenjets_=iConfig.getUntrackedParameter<bool>("enablegenjets",false);
        enablepfsimparticles_=iConfig.getUntrackedParameter<bool>("enablepfsimparticles",false);
	jetcollection_=iConfig.getUntrackedParameter< string >("jetcollection", 
	"Fastjet10");
	LogDebug("TauBenchmarkAnalyzer")
	<<"opening output root file "<<outputRootFileName<<endl;
	
	file_ = TFile::Open(outputRootFileName.c_str(), "RECREATE");
	
	if(file_->IsZombie() ) {
		string err = "output file ";
		err += outputRootFileName;
		err += " can't be opened";
		throw cms::Exception("OutputFileOpen",err);
	}
	// Comparion RecJets with PFSIMParticle 
	h_deltaETvisible_MCPFSIM_EHT_
	= new TH1F("h_deltaETvisible_MCPFSIM_EHT","Jet Et difference CaloTowers-PFSIM"
	,500,-50,50);	
	h_deltaETvisible_MCPFSIM_PF_
	= new TH1F("h_deltaETvisible_MCPFSIM_PF_" ,"Jet Et difference ParticleFlow-PFSIM"
	,500,-50,50);
		
	//Comparion RecJets with HepmMC 
	h_deltaETvisible_MCHEPMC_PF_
	= new TH1F("h_deltaETvisible_MCHEPMC_PF_" , "Jet Et difference PF -HepMC Truth"
	,500,-50,50);
	h_deltaETvisible_MCHEPMC_EHT_
	= new TH1F("h_deltaETvisible_MCHEPMC_EHT_" , "Jet Et difference EHT - HepMC Truth"
	,500,-50,50);
	
	// Comparison RecJets with Gen Jets 
	h_deltaETvisible_MCGENJET_EHT_
        = new TH1F("h_deltaETvisible_MCGENJET_EHT_","Jet Et difference CaloTowers-GenJets"
        ,500,-50,50);
	h_deltaETvisible_MCGENJET_PF_
        = new TH1F("h_deltaETvisible_MCGENJET_PF_","Jet Et difference PFJets-GenJets"
        ,500,-50,50);

	// Comparison HepMC/PFSIMParticle  with GenJets 
	h_deltaETvisible_GENJET_PFSIM_
	= new TH1F("h_deltaETvisible_MCGENJET_PFSIM_" , "Jet Et difference GenJet - PFSimP"
	,500,-50,50);
       	h_deltaETvisible_GENJET_HEPMC_
	= new TH1F("h_deltaETvisible_MCGENJET_HEPMC" , "Jet Et difference GenJet - HepMC"
	,500,-50,50);

	// multiplicities
	h_genJetMatching_MCEHTnum_
        = new TH1F("h_genJetMatching_MCEHTnum_","Jet num difference CaloTowers-GenJets"
        ,21,-10,10);
	h_genJetMatching_MCPFnum_
        = new TH1F("h_genJetMatching_MCPFnum_","Jet num difference PFJets-GenJets"
        ,21,-10,10);
	h_genJetMatching_MCnum_
        = new TH1F("h_genJetMatching_MCnum_","Jet num difference Trueparticle-GenJets"
        ,21,-10,10); 
}


TauBenchmarkAnalyzer::~TauBenchmarkAnalyzer()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}


void TauBenchmarkAnalyzer::beginJob(const edm::EventSetup&)
{
        LogInfo("Info")<<"Starting MyAnalyzer"<<endl;
	
}


void TauBenchmarkAnalyzer::endJob() {
	
	file_->Write();
	file_->Close();
}


void TauBenchmarkAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	LogInfo("Info")<<endl<<endl<<endl;
	//this is the event loop
	//double total_et=0.0;
	//double total_et_calo=0.0;
	trueParticles_.clear();
	caloJets_.clear();
	pfJets_.clear();
	genJets_.clear();
	if(enablepfsimparticles_==true)
	  iEvent.getByLabel("particleFlowSimParticle", trueParticles_);
	iEvent.getByLabel(jetcollection_+"CaloJets", caloJets_);
	iEvent.getByLabel(jetcollection_+"PFJets", pfJets_);
	if(enablegenjets_==true)
	  {
	    iEvent.getByLabel(jetcollection_+"GenJets", genJets_);
	    iEvent.getByLabel("VtxSmeared", hepMC_);
	  }
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
				LogInfo("Info")<< "DAUGHTER=" << pdgdaugter << endl;
				if (pdgdaugter == 11 || pdgdaugter == 13) { 
					return;
					
				}//electron or muons?
			}//loop daughter
		}//tau
	}//loop particle
	// ====================================================================
	//MAKING TRUE PARTICLES FROM PFSIMPARTICLES ===========================
	// ====================================================================
	
	int jettrue_num=0;
	double trueparticle_et=0.0;
	if(enablepfsimparticles_==true)
	  {
	    TLorentzVector partTOTMC;
	    for ( unsigned i=0;  i < (*trueParticles_).size(); i++) {
	      const reco::PFSimParticle& ptc = (*trueParticles_)[i];
	      vectPART.push_back(ptc);
	    }//loop
	    
	    for ( unsigned i=0;  i < (*trueParticles_).size(); i++) 
	      {  
	        LogInfo("Info")<< "pdgcode:"<<  (*trueParticles_)[i].pdgCode()<< endl;
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
		    LogInfo("Info")<< "pdgcode:"<<pdgcode << endl;
		  }//loop daughter
		  jettrue_num++;	
		}//tau?
	      }//loop particles
	    trueparticle_et=partTOTMC.Et();
	
	    LogInfo("Info")<<"total et PFSIMPARTICLE "<<trueparticle_et<<endl;
	  }
	
	// ====================================================================
	
	//MAKING TRUE PARTICLES FROM HEPMC-PARTICLES ==========================
	// ====================================================================
	double JetGENETmax=0.0;
	double jetgen_et=0.0;
	int jetgen_num=0;
	double neutrino_et=0.0;
	double hepmc_et =0.0; 
	if(enablegenjets_==true)
	  {
	    const HepMC::GenEvent * myGenEvent = hepMC_->GetEvent();
	  
	    for (HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();   p != myGenEvent->particles_end(); ++p )
	      {
		if (abs((*p)->pdg_id()) == 15) {
		  // retrieve decay vertex
		  const HepMC::GenVertex * decayVtx = (*p)->end_vertex();
		  if ( decayVtx != 0 ) {
		    HepMC::GenVertex::particles_in_const_iterator child_mcpartItr  = decayVtx->particles_out_const_begin();
		    HepMC::GenVertex::particles_in_const_iterator child_mcpartItrE = decayVtx->particles_out_const_end();
		    
		    for (; child_mcpartItr != child_mcpartItrE; ++child_mcpartItr) {
		      
		      HepMC::GenParticle * child = (*child_mcpartItr);
		      if (std::abs(child->pdg_id())!=12 && std::abs(child->pdg_id())!=14 && std::abs(child->pdg_id())!=16 ){	
			hepmc_et+=child->momentum().et();
			LogInfo("Info")<< "pdgcode:"<< child->momentum().et()
				       << endl;}
		    }//loop daughter
		  }//  decayVtx
		}//if tau//
	      }
	    LogInfo("Info")<<"HepMC truth MC without neutrino " << hepmc_et<< endl; 
	    
	    
	    // ====================================================================
	    //MAKING TRUE GEN JETS============================================
	    // ====================================================================
	    
	 
	    //const HepMC::GenEvent * myGenEvent = hepMC_->GetEvent();
	    std::vector<const HepMC::GenParticle *> nu;
	    for (HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();   p != myGenEvent->particles_end(); ++p )
	      {
		LogInfo("Info")<<" HEPMC PDG "<< (*p)->pdg_id() <<endl;
		if ( (std::abs((*p)->pdg_id())==12 || std::abs((*p)->pdg_id())==14 || std::abs((*p)->pdg_id())==16 ) && (*p)->status()==1 )
		  nu.push_back(*p);
		//if(mu.size()>1) break;
	      }
	    LogInfo("Info")<<"Number of Neutrinos in Jets: "<<nu.size()<<endl;
	    for(unsigned int i=0; i<nu.size();i++)
	      {
		neutrino_et+=nu[i]->momentum().et();
	      }
	    LogInfo("Info")<<"sum of neutrino et: "<<neutrino_et<<endl;
	    LogInfo("Info")<<"genjet.size: "<<(*genJets_).size()<<endl;
	    for ( unsigned int i = 0; i < (*genJets_).size(); i++)
	      {
		//Taking the most energetic jet
		jetgen_et=(*genJets_)[i].et();
		if(JetGENETmax<jetgen_et)
		  JetGENETmax = jetgen_et;
		
		std::vector <const reco::GenParticleCandidate*> jetConstituents;		
		jetConstituents=(*genJets_)[i].getConstituents();
		LogInfo("Info")<<"Number of Jet Constituents:"<<jetConstituents.size()<<endl;
		LogInfo("Info")<<"Number of Jet Daughters: "<<(*genJets_)[i].numberOfDaughters()<<endl;
	      }//loop calo towers
	    //jetgen_num=(*genJets_).size();
	    LogInfo("Info")<<"total et jetgen_et: "<<JetGENETmax<<endl;
	    //cout<<"number of particle gen: "<<jetgen_num<<endl;
	  }//enablegenjets_
	// ==================================================================
	//CALO TOWER JETS (ECAL+HCAL Towers)=================================
	// ==================================================================
	
	double JetEHTETmax=0.0;
	double jetcalo_et=0.0;
	int jetcalo_num=0;
	
	for ( unsigned int i = 0; i < (*caloJets_).size(); ++i)
	{
		jetcalo_et=(*caloJets_)[i].et();
		//if (jetcalo_et >= JetEHTETmax) 
		JetEHTETmax += jetcalo_et;
	}//loop calo towers
	jetcalo_num=(*caloJets_).size();
	LogInfo("Info")<<"total et calo: "<<JetEHTETmax<<endl;
	LogInfo("Info")<<"number of particle calo: "<<jetcalo_num<<endl;	
	
	// ==================================================================
	//PF Jets ===========================================================
	// ==================================================================
	
	double pfjet_et=0.0;
	double JetPFETmax=0.0;
	int jetpf_num=0;
	for ( unsigned int i = 0; i < (*pfJets_).size(); ++i)
	{
		pfjet_et=(*pfJets_)[i].et();
		JetPFETmax += pfjet_et;
	}//loop calo towers
	jetpf_num=(*pfJets_).size();
	LogInfo("Info")<<"total et pfjet: "<<JetPFETmax<<endl;
	LogInfo("Info")<<"number of particle pf: "<<jetpf_num<<endl;	
	
	// ==================================================================
	// Status output ====================================================
	// ==================================================================
	
	LogInfo("Info")<<"delta et 1_1:  "<<JetEHTETmax-JetGENETmax <<endl;
	//  cout<<"delta et 1_2 : "<<JetEHTETmax-trueparticle_et<<endl; 
	LogInfo("Info")<<"delta et 2_1 : "<<JetPFETmax-JetGENETmax+neutrino_et<<endl; 
	LogInfo("Info")<<"delta et 1_2:  "<<JetEHTETmax-JetGENETmax<<endl;
	LogInfo("Info")<<"delta et 2_2 : "<<JetPFETmax-JetGENETmax<<endl; 
	//cout<<"delta et 2_2 : "<<JetPFETmax-trueparticle_et<<endl; 
	//cout<<"delta num 1:  "<<jetcalo_num-jetgen_num<<endl;
	//cout<<"delta num 2 : "<<jetpf_num-jetgen_num<<endl;
	
	
	// ==================================================================
	//Filling Histogramms ===============================================
	// ==================================================================
	// PFSIM Comp
	h_deltaETvisible_MCPFSIM_PF_->Fill(JetPFETmax-trueparticle_et);
	h_deltaETvisible_MCPFSIM_EHT_->Fill(JetEHTETmax-trueparticle_et);
	//HepMC Comp
	h_deltaETvisible_MCHEPMC_PF_->Fill(JetPFETmax-hepmc_et);
	h_deltaETvisible_MCHEPMC_EHT_->Fill(JetEHTETmax-hepmc_et); 
	//GenJet Comp
	h_deltaETvisible_MCGENJET_PF_->Fill(JetPFETmax-JetGENETmax);
	h_deltaETvisible_MCGENJET_EHT_->Fill(JetEHTETmax-JetGENETmax); 		
	// GenJet / True Particle  (later parton) 
	h_deltaETvisible_GENJET_PFSIM_->Fill(JetGENETmax-trueparticle_et);
	h_deltaETvisible_GENJET_HEPMC_->Fill(JetGENETmax-hepmc_et);
	//multi
	h_genJetMatching_MCEHTnum_->Fill(jetcalo_num-jetgen_num);
	h_genJetMatching_MCPFnum_->Fill(jetpf_num-jetgen_num);
	h_genJetMatching_MCnum_->Fill(jettrue_num-jetgen_num);
	
}












