#include "RecoParticleFlow/PFProducer/plugins/PFBlockProducer.h"

#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFNuclearInteraction.h"
#include "DataFormats/ParticleFlowReco/interface/PFConversionFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFConversion.h"
#include "DataFormats/ParticleFlowReco/interface/PFV0Fwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFV0.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <set>

using namespace std;
using namespace edm;

PFBlockProducer::PFBlockProducer(const edm::ParameterSet& iConfig) {
  

  // use configuration file to setup input/output collection names
  inputTagRecTracks_ 
    = iConfig.getParameter<InputTag>("RecTracks");

  inputTagGsfRecTracks_ 
    = iConfig.getParameter<InputTag>("GsfRecTracks");

  inputTagRecMuons_ 
    = iConfig.getParameter<InputTag>("RecMuons");

  inputTagPFNuclear_ 
    = iConfig.getParameter<InputTag>("PFNuclear");

  inputTagPFConversions_ 
    = iConfig.getParameter<InputTag>("PFConversions");

  inputTagPFV0_ 
    = iConfig.getParameter<InputTag>("PFV0");

  inputTagPFClustersECAL_ 
    = iConfig.getParameter<InputTag>("PFClustersECAL");

  inputTagPFClustersHCAL_ 
    = iConfig.getParameter<InputTag>("PFClustersHCAL");

  inputTagPFClustersHFEM_ 
    = iConfig.getParameter<InputTag>("PFClustersHFEM");

  inputTagPFClustersHFHAD_ 
    = iConfig.getParameter<InputTag>("PFClustersHFHAD");

  inputTagPFClustersPS_ 
    = iConfig.getParameter<InputTag>("PFClustersPS");



  verbose_ = 
    iConfig.getUntrackedParameter<bool>("verbose",false);

  bool debug_ = 
    iConfig.getUntrackedParameter<bool>("debug",false);

  usePFatHLT_ = iConfig.getParameter<bool>("usePFatHLT");

  useNuclear_ = iConfig.getParameter<bool>("useNuclear");

  useConversions_ = iConfig.getParameter<bool>("useConversions");

  useV0_ = iConfig.getParameter<bool>("useV0");

  produces<reco::PFBlockCollection>();
  

  
  // particle flow parameters  -----------------------------------

  string map_ECAL_eta 
    = iConfig.getParameter<string>("pf_resolution_map_ECAL_eta");  
  string map_ECAL_phi 
    = iConfig.getParameter<string>("pf_resolution_map_ECAL_phi");  
  //   will be necessary when preshower is used:
  //   string map_ECALec_x 
  //     = iConfig.getParameter<string>("pf_resolution_map_ECALec_x");  
  //   string map_ECALec_y 
  //     = iConfig.getParameter<string>("pf_resolution_map_ECALec_y");  
  string map_HCAL_eta 
    = iConfig.getParameter<string>("pf_resolution_map_HCAL_eta");  
  string map_HCAL_phi 
    = iConfig.getParameter<string>("pf_resolution_map_HCAL_phi"); 
	
  std::vector<double> DPtovPtCut 
     = iConfig.getParameter<std::vector<double> >("pf_DPtoverPt_Cut");   

  std::vector<unsigned> NHitCut 
     = iConfig.getParameter<std::vector<unsigned> >("pf_NHit_Cut");   

  double chi2_ECAL_PS 
     = iConfig.getParameter<double>("pf_chi2_ECAL_PS");  
//   double chi2_HCAL_PS 
//     = iConfig.getParameter<double>("pf_chi2_HCAL_PS");  

  double chi2_ECAL_Track 
    = iConfig.getParameter<double>("pf_chi2_ECAL_Track");  
  double chi2_ECAL_GSF 
    = iConfig.getParameter<double>("pf_chi2_ECAL_GSF");  
  double chi2_HCAL_Track 
    = iConfig.getParameter<double>("pf_chi2_HCAL_Track");  
  double chi2_ECAL_HCAL 
    = iConfig.getParameter<double>("pf_chi2_ECAL_HCAL");  
  double chi2_PS_Track 
    = iConfig.getParameter<double>("pf_chi2_PS_Track");  
  double chi2_PSH_PSV 
    = iConfig.getParameter<double>("pf_chi2_PSH_PSV");  
  
  bool multiLink = 
    iConfig.getParameter<bool>("pf_multilink");
  
  //energyCalibration_ = new PFEnergyCalibration(iConfig);

  
  //   PFBlock::setEnergyResolution(energyResolution_);

  edm::FileInPath path_ECAL_eta( map_ECAL_eta.c_str() );
  edm::FileInPath path_ECAL_phi( map_ECAL_phi.c_str() );
  edm::FileInPath path_HCAL_eta( map_HCAL_eta.c_str() );
  edm::FileInPath path_HCAL_phi( map_HCAL_phi.c_str() );
   
  pfBlockAlgo_.setParameters( path_ECAL_eta.fullPath().c_str(),
			      path_ECAL_phi.fullPath().c_str(),
			      path_HCAL_eta.fullPath().c_str(),
			      path_HCAL_phi.fullPath().c_str(),
			      DPtovPtCut,
			      NHitCut,
			      chi2_ECAL_Track,
			      chi2_ECAL_GSF,
			      chi2_HCAL_Track,
			      chi2_ECAL_HCAL,
			      chi2_ECAL_PS,
			      chi2_PS_Track,
			      chi2_PSH_PSV,
			      multiLink );
  
  pfBlockAlgo_.setDebug(debug_);

//   energyCalibration_ = new PFEnergyCalibration();
//   double calibParamECAL_slope_ 
//     = iConfig.getParameter<double>("pf_ECAL_calib_p1");
//   double calibParamECAL_offset_ 
//     = iConfig.getParameter<double>("pf_ECAL_calib_p0");
  
//   energyCalibration_->setCalibrationParametersEm(calibParamECAL_slope_, calibParamECAL_offset_); 
  
//   //   PFBlock::setEnergyCalibration(energyCalibration_);
//   //energyResolution_ = new PFEnergyResolution(iConfig);

//   energyResolution_ = new PFEnergyResolution();


}



PFBlockProducer::~PFBlockProducer() { }



void 
PFBlockProducer::beginJob() { }

void 
PFBlockProducer::beginRun(edm::Run & run, 
			  const edm::EventSetup & es) { }


void 
PFBlockProducer::produce(Event& iEvent, 
			 const EventSetup& iSetup) {
  
  LogDebug("PFBlockProducer")<<"START event: "<<iEvent.id().event()
			     <<" in run "<<iEvent.id().run()<<endl;
  
  
  // get rectracks
  
  Handle< reco::PFRecTrackCollection > recTracks;
  
  // LogDebug("PFBlockProducer")<<"get reco tracks"<<endl;
  bool found = iEvent.getByLabel(inputTagRecTracks_, recTracks);
    
  if(!found )
    LogError("PFBlockProducer")<<" cannot get rectracks: "
			       <<inputTagRecTracks_<<endl;



  // get GsfTracks 
  Handle< reco::GsfPFRecTrackCollection > GsfrecTracks;

  if(!usePFatHLT_) {
    found = iEvent.getByLabel(inputTagGsfRecTracks_,GsfrecTracks);

    if(!found )
      LogError("PFBlockProducer")<<" cannot get Gsfrectracks: "
				 << inputTagGsfRecTracks_ <<endl;
  }
  // get recmuons

  Handle< reco::MuonCollection > recMuons;

  // LogDebug("PFBlockProducer")<<"get reco muons"<<endl;
  if(!usePFatHLT_) {
    found = iEvent.getByLabel(inputTagRecMuons_, recMuons);
  
    //if(!found )
    //  LogError("PFBlockProducer")<<" cannot get recmuons: "
    //			 <<inputTagRecMuons_<<endl;

  // get PFNuclearInteractions
  }
  Handle< reco::PFNuclearInteractionCollection > pfNuclears;
  if( useNuclear_ ) {
    found = iEvent.getByLabel(inputTagPFNuclear_, pfNuclears);


    if(!found )
      LogError("PFBlockProducer")<<" cannot get PFNuclearInteractions : "
                               <<inputTagPFNuclear_<<endl;
  }



 
  // get conversions
  Handle< reco::PFConversionCollection > pfConversions;
  if( useConversions_ ) {
    found = iEvent.getByLabel(inputTagPFConversions_, pfConversions);
    
    if(!found )
      LogError("PFBlockProducer")<<" cannot get PFConversions : "
				 <<inputTagPFConversions_<<endl;
  }
  

  // get V0s
  Handle< reco::PFV0Collection > pfV0;
  if( useV0_ ) {
    found = iEvent.getByLabel(inputTagPFV0_, pfV0);
    
    if(!found )
      LogError("PFBlockProducer")<<" cannot get PFV0 : "
				 <<inputTagPFV0_<<endl;
  }


  
  // get ECAL, HCAL and PS clusters
  
  
  Handle< reco::PFClusterCollection > clustersECAL;
  found = iEvent.getByLabel(inputTagPFClustersECAL_, 
			    clustersECAL);      
  if(!found )
    LogError("PFBlockProducer")<<" cannot get ECAL clusters: "
			       <<inputTagPFClustersECAL_<<endl;
    
  
  Handle< reco::PFClusterCollection > clustersHCAL;
  found = iEvent.getByLabel(inputTagPFClustersHCAL_, 
			    clustersHCAL);      
  if(!found )
    LogError("PFBlockProducer")<<" cannot get HCAL clusters: "
			       <<inputTagPFClustersHCAL_<<endl;
    
  Handle< reco::PFClusterCollection > clustersHFEM;
  found = iEvent.getByLabel(inputTagPFClustersHFEM_, 
			    clustersHFEM);      
  if(!found )
    LogError("PFBlockProducer")<<" cannot get HFEM clusters: "
			       <<inputTagPFClustersHFEM_<<endl;
    
  Handle< reco::PFClusterCollection > clustersHFHAD;
  found = iEvent.getByLabel(inputTagPFClustersHFHAD_, 
			    clustersHFHAD);      
  if(!found )
    LogError("PFBlockProducer")<<" cannot get HFHAD clusters: "
			       <<inputTagPFClustersHFHAD_<<endl;
    

  Handle< reco::PFClusterCollection > clustersPS;
  found = iEvent.getByLabel(inputTagPFClustersPS_, 
			    clustersPS);      
  if(!found )
    LogError("PFBlockProducer")<<" cannot get PS clusters: "
			       <<inputTagPFClustersPS_<<endl;
    
  if( usePFatHLT_  ) {
     pfBlockAlgo_.setInput( recTracks, 			   
			   clustersECAL,
			   clustersHCAL,
			   clustersHFEM,
			   clustersHFHAD,
			   clustersPS );
  } else { 
    pfBlockAlgo_.setInput( recTracks, 
			   GsfrecTracks,
			   recMuons, 
			   pfNuclears,
			   pfConversions,
			   pfV0,
			   clustersECAL,
			   clustersHCAL,
			   clustersHFEM,
			   clustersHFHAD,
			   clustersPS );
  }
  pfBlockAlgo_.findBlocks();
  
  if(verbose_) {
    ostringstream  str;
    str<<pfBlockAlgo_<<endl;
    LogInfo("PFBlockProducer") << str.str()<<endl;
  }    
  
  auto_ptr< reco::PFBlockCollection > 
    pOutputBlockCollection( pfBlockAlgo_.transferBlocks() ); 
  
  
  iEvent.put(pOutputBlockCollection);
  
  LogDebug("PFBlockProducer")<<"STOP event: "<<iEvent.id().event()
			     <<" in run "<<iEvent.id().run()<<endl;
}