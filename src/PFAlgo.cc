#include "RecoParticleFlow/PFProducer/interface/PFAlgo.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"  //PFMuons
#include "RecoParticleFlow/PFProducer/interface/PFElectronAlgo.h"  
#include "RecoParticleFlow/PFProducer/interface/PFConversionAlgo.h"  

#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFClusterCalibration.h" 

#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// #include "FWCore/Framework/interface/OrphanHandle.h"
#include "DataFormats/Common/interface/OrphanHandle.h"
// #include "DataFormats/Common/interface/ProductID.h"
#include "DataFormats/Provenance/interface/ProductID.h"


#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "Math/PxPyPzM4D.h"
#include "Math/LorentzVector.h"
#include "Math/DisplacementVector3D.h"
#include "Math/SMatrix.h"
#include "TDecompChol.h"
#include "TMVA/Config.h" // needed for TMVA::gConfig() workaround below

#include "boost/graph/adjacency_matrix.hpp" 
#include "boost/graph/graph_utility.hpp" 

using namespace std;
using namespace reco;
using namespace boost;

float PFAlgo::maxMvaCut_ = PFCandidate::bigMva_;


typedef std::list< reco::PFBlockRef >::iterator IBR;



PFAlgo::PFAlgo()
  : pfCandidates_( new PFCandidateCollection),
    nSigmaECAL_(0), 
    nSigmaHCAL_(1),
    PSCut_(999.),
    mergedPhotonsMVA_( 0 ),
    mvaCut_(PFAlgo::maxMvaCut_),
    algo_(1),
    debug_(false),
    pfele_(0)
{}

PFAlgo::~PFAlgo() {
  if (usePFElectrons_) delete pfele_;
  if (usePFConversions_) delete pfConversion_;
}


void 
PFAlgo::setParameters(double nSigmaECAL,
                      double nSigmaHCAL, 
                      const shared_ptr<PFEnergyCalibration>& calibration,
                      const shared_ptr<pftools::PFClusterCalibration>& clusterCalibration,
		      unsigned int newCalib,
                      double PSCut,
                      double mvaCut,
                      const  char* mvaWeightFile ) {

  nSigmaECAL_ = nSigmaECAL;
  nSigmaHCAL_ = nSigmaHCAL;

  calibration_ = calibration;
  clusterCalibration_ = clusterCalibration;
  newCalib_ = newCalib;
  // std::cout << "Cluster calibration parameters : " << *clusterCalibration_ << std::endl;

  PSCut_=PSCut;
  mvaCut_ = mvaCut;
  
  if( mergedPhotonsMVA_ ) delete mergedPhotonsMVA_;
  mergedPhotonsMVA_ = new TMVA::Reader();
  TMVA::gConfig().SetSilent(kTRUE); // workaround for ROOT 5.18/00 bug
  mergedPhotonsMVA_->AddVariable("eECAL/pTrack", &eECALOverpTrack_);
  mergedPhotonsMVA_->AddVariable("chi2ECAL", &chi2ECAL_);
  mergedPhotonsMVA_->AddVariable("ptTrack", &ptTrack_);
  
  if( mvaCut < PFAlgo::maxMvaCut_ ) {
    FILE * file = fopen(mvaWeightFile, "r");
    if (file) {
      fclose(file);
      // file is readable, book MVA
      mergedPhotonsMVA_->BookMVA( "MVA", mvaWeightFile );
    }
    else {
      // weight file is not readable
      string err = "PFAlgo: cannot open weight file '";
      err += mvaWeightFile;
      err += "'";
      throw invalid_argument( err );
    }
  }
}

//PFElectrons: a new method added to set the parameters for electron reconstruction. 
void 
PFAlgo::setPFEleParameters(double chi2EcalGSF,
			   double chi2EcalBrem,
			   double chi2HcalGSF,
			   double chi2HcalBrem,
			   double chi2PsGSF,
			   double chi2PsBrem,
			   double mvaEleCut,
			   string mvaWeightFileEleID,
			   bool usePFElectrons) {
  setchi2Values_.push_back(chi2EcalGSF);
  setchi2Values_.push_back(chi2EcalBrem);
  setchi2Values_.push_back(chi2HcalGSF);
  setchi2Values_.push_back(chi2HcalBrem);
  setchi2Values_.push_back(chi2PsGSF);
  setchi2Values_.push_back(chi2PsBrem);
  mvaEleCut_ = mvaEleCut;
  usePFElectrons_ = usePFElectrons;
  if(!usePFElectrons_) return;
  mvaWeightFileEleID_ = mvaWeightFileEleID;
  FILE * fileEleID = fopen(mvaWeightFileEleID_.c_str(), "r");
  if (fileEleID) {
    fclose(fileEleID);
  }
  else {
    string err = "PFAlgo: cannot open weight file '";
    err += mvaWeightFileEleID;
    err += "'";
    throw invalid_argument( err );
  }
  pfele_= new PFElectronAlgo(setchi2Values_,mvaEleCut_,mvaWeightFileEleID_);
}

void 
PFAlgo::setPFMuonAndFakeParameters(std::vector<double> muonHCAL,
				   std::vector<double> muonECAL,
				   double nSigmaTRACK,
				   double ptError,
				   std::vector<double> factors45) 
{
  muonHCAL_ = muonHCAL;
  muonECAL_ = muonECAL;
  nSigmaTRACK_ = nSigmaTRACK;
  ptError_ = ptError;
  factors45_ = factors45;
}
  

void PFAlgo::setPFConversionParameters(bool usePFConversions ) {

  usePFConversions_ = usePFConversions;
  pfConversion_ = new PFConversionAlgo();


}



void PFAlgo::reconstructParticles( const reco::PFBlockHandle& blockHandle ) {

  blockHandle_ = blockHandle;
  reconstructParticles( *blockHandle_ ); 
}



void PFAlgo::reconstructParticles( const reco::PFBlockCollection& blocks ) {

  // reset output collection
  if(pfCandidates_.get() )
    pfCandidates_->clear();
  else 
    pfCandidates_.reset( new reco::PFCandidateCollection );

  if(pfElectronCandidates_.get() )
    pfElectronCandidates_->clear();
  else
    pfElectronCandidates_.reset( new reco::PFCandidateCollection);
  

  if( debug_ ) {
    cout<<"*********************************************************"<<endl;
    cout<<"*****           Particle flow algorithm             *****"<<endl;
    cout<<"*********************************************************"<<endl;
  }

  // sort elements in three lists:
  std::list< reco::PFBlockRef > hcalBlockRefs;
  std::list< reco::PFBlockRef > ecalBlockRefs;
  std::list< reco::PFBlockRef > otherBlockRefs;
  
  for( unsigned i=0; i<blocks.size(); ++i ) {
    // reco::PFBlockRef blockref( blockh,i );
    reco::PFBlockRef blockref = createBlockRef( blocks, i);
    
    const reco::PFBlock& block = *blockref;
   const edm::OwnVector< reco::PFBlockElement >& 
      elements = block.elements();
        
    bool singleEcalOrHcal = false;
    if( elements.size() == 1 ){
      if( elements[0].type() == reco::PFBlockElement::ECAL ){
        ecalBlockRefs.push_back( blockref );
        singleEcalOrHcal = true;
      }
      if( elements[0].type() == reco::PFBlockElement::HCAL ){
        hcalBlockRefs.push_back( blockref );
        singleEcalOrHcal = true;
      }
    }      
    
    if(!singleEcalOrHcal) {
      otherBlockRefs.push_back( blockref );
    }
  }//loop blocks
  
  if( debug_ ){
    cout<<"# Ecal blocks: "<<ecalBlockRefs.size()
        <<", # Hcal blocks: "<<hcalBlockRefs.size()
        <<", # Other blocks: "<<otherBlockRefs.size()<<endl;
  }


  // loop on blocks that are not single ecal, 
  // and not single hcal.

  for( IBR io = otherBlockRefs.begin(); io!=otherBlockRefs.end(); ++io) {
    processBlock( *io, hcalBlockRefs, ecalBlockRefs );
  }

  std::list< reco::PFBlockRef > empty;

  // process remaining single hcal blocks
  for( IBR ih = hcalBlockRefs.begin(); ih!=hcalBlockRefs.end(); ++ih) {
    processBlock( *ih, empty, empty );
  }

  // process remaining single ecal blocks
  for( IBR ie = ecalBlockRefs.begin(); ie!=ecalBlockRefs.end(); ++ie) {
    processBlock( *ie, empty, empty );
  }
}


void PFAlgo::processBlock( const reco::PFBlockRef& blockref,
                           std::list<reco::PFBlockRef>& hcalBlockRefs, 
                           std::list<reco::PFBlockRef>& ecalBlockRefs ) { 
  
  assert(!blockref.isNull() );
  const reco::PFBlock& block = *blockref;

  typedef std::multimap<double, unsigned>::iterator IE;
  typedef std::multimap<double, std::pair<unsigned,double> >::iterator IS;
  typedef std::multimap<double, std::pair<unsigned,bool> >::iterator IT;

  if(debug_) {
    cout<<"#########################################################"<<endl;
    cout<<"#####           Process Block:                      #####"<<endl;
    cout<<"#########################################################"<<endl;
    cout<<block<<endl;
  }

  const edm::OwnVector< reco::PFBlockElement >& elements = block.elements();

  // make a copy of the link data, which will be edited.
  PFBlock::LinkData linkData =  block.linkData();
  
  // keep track of the elements which are still active.
  vector<bool>   active( elements.size(), true );


  // //PFElectrons:
  // usePFElectrons_ external configurable parameter to set the usage of pf electron
  if (usePFElectrons_) {
    if (pfele_->isElectronValidCandidate(blockref,active)){
      // if there is at least a valid candidate is get the vector of pfcandidates
      std::vector<reco::PFCandidate> PFElectCandidates_;
      PFElectCandidates_ = pfele_->getElectronCandidates();
      for ( std::vector<reco::PFCandidate>::iterator ec = PFElectCandidates_.begin();
	    ec != PFElectCandidates_.end(); ++ec )
	{
	  pfCandidates_->push_back(*ec);
	  // the pfalgo candidates vector is filled
	}
      // The vector active is automatically changed (it is passed by ref) in PFElectronAlgo
      // for all the electron candidate      
    }
    pfElectronCandidates_->insert(pfElectronCandidates_->end(), 
                                  pfele_->getAllElectronCandidates().begin(), 
                                  pfele_->getAllElectronCandidates().end()); 
    std::vector<reco::PFCandidate>::const_iterator ec=pfele_->getAllElectronCandidates().begin(); 
    /*
    for(;ec!=pfele_->getAllElectronCandidates().end();++ec) 
      { 
        unsigned nd=ec->numberOfDaughters(); 
        for(unsigned i=0;i<nd;++i) 
	  std::cout << *(PFCandidate*)ec->daughter(i) << std::endl; 
      } 
    */
  }

  
  // usePFConversions_ is used to switch ON/OFF the use of the PFConversionAlgo
  if (usePFConversions_) {
    if (pfConversion_->isConversionValidCandidate(blockref, active )){
      // if at least one conversion candidate is found ,it is fed to the final list of Pflow Candidates 
      std::vector<reco::PFCandidate> PFConversionCandidates = pfConversion_->conversionCandidates();
            
      for ( std::vector<reco::PFCandidate>::iterator iConv = PFConversionCandidates.begin(); 
	    iConv != PFConversionCandidates.end(); ++iConv ) {
	pfCandidates_->push_back(*iConv);
      }
      for(unsigned iEle=0; iEle<elements.size(); iEle++) {
	//cout << " PFAlgo::processBlock conversion element " << iEle << " activ e " << active[iEle] << endl;
      }
    }
  }



  if(debug_) 
    cout<<endl<<"--------------- loop 1 ------------------"<<endl;

  //COLINFEB16
  // In loop 1, we loop on all elements. 

  // The primary goal is to deal with tracks that are:
  // - not associated to an HCAL cluster
  // - not identified as an electron. 
  // Those tracks should be predominantly relatively low energy charged 
  // hadrons which are not detected in the ECAL. 

  // The secondary goal is to prepare for the next loops
  // - The ecal and hcal elements are sorted in separate vectors
  // which will be used as a base for the corresponding loops.
  // - For tracks which are connected to more than one HCAL cluster, 
  // the links between the track and the cluster are cut for all clusters 
  // but the closest one. 
  // - HF only blocks ( HFEM, HFHAD, HFEM+HFAD) are identified 

  // obsolete comments? 
  // loop1: 
  // - sort ecal and hcal elements in separate vectors
  // - for tracks:
  //       - lock closest ecal cluster 
  //       - cut link to farthest hcal cluster, if more than 1.

  // vectors to store indices to hcal and ecal elements 
  vector<unsigned> hcalIs;
  vector<unsigned> ecalIs;
  vector<unsigned> trackIs;
  vector<unsigned> ps1Is;
  vector<unsigned> ps2Is;

  vector<unsigned> hfEmIs;
  vector<unsigned> hfHadIs;
  
  for(unsigned iEle=0; iEle<elements.size(); iEle++) {

    PFBlockElement::Type type = elements[iEle].type();
        
    if(debug_ && type != PFBlockElement::BREM ) cout<<endl<<elements[iEle];   

    switch( type ) {
    case PFBlockElement::TRACK:
      if ( active[iEle] ) { 
	trackIs.push_back( iEle );
	if(debug_) cout<<"TRACK, stored index, continue"<<endl;
      }
      break;
    case PFBlockElement::ECAL: 
      if ( active[iEle]  ) { 
	ecalIs.push_back( iEle );
	if(debug_) cout<<"ECAL, stored index, continue"<<endl;
      }
      continue;
    case PFBlockElement::HCAL:
      if ( active[iEle] ) { 
	hcalIs.push_back( iEle );
	if(debug_) cout<<"HCAL, stored index, continue"<<endl;
      }
      continue;
    case PFBlockElement::HFEM:
      if ( active[iEle] ) { 
	hfEmIs.push_back( iEle );
	if(debug_) cout<<"HFEM, stored index, continue"<<endl;
      }
      continue;
    case PFBlockElement::HFHAD:
      if ( active[iEle] ) { 
	hfHadIs.push_back( iEle );
	if(debug_) cout<<"HFHAD, stored index, continue"<<endl;
      }
      continue;

      //COLINFEB16 the following is not used... removing
      
      //     case PFBlockElement::PS1:
      //        if ( active[iEle] ) { 
      // 	 ps1Is.push_back( iEle );
      // 	 if(debug_) cout<<"PS1, stored index, continue"<<endl;
      //        }
      //        continue;
      //     case PFBlockElement::PS2:
      //       if ( active[iEle] ) { 
      // 	ps2Is.push_back( iEle );
      // 	if(debug_) cout<<"PS2, stored index, continue"<<endl;
      //       }
      //       continue;
    default:
      continue;
    }

    // we're now dealing with a track
    unsigned iTrack = iEle;
    if(debug_) { 
      cout<<"TRACK"<<endl;
      if ( !active[iTrack] ) 
	cout << "Already used by electrons, muons, conversions" << endl;
    }
    
    // Track already used as electron, muon, conversion?
    // Added a check on the activated element
    if ( ! active[iTrack] ) continue;
  
    reco::TrackRef trackRef = elements[iTrack].trackRef();
    assert( !trackRef.isNull() ); 
                      

    // look for associated elements of all types
    //COLINFEB16 
    // all types of links are considered. 
    // the elements are sorted by increasing distance, not chi2
    std::multimap<double, unsigned> ecalElems;
    block.associatedElements( iTrack,  linkData,
                              ecalElems ,
                              reco::PFBlockElement::ECAL,
			      reco::PFBlock::LINKTEST_ALL );
    
    std::multimap<double, unsigned> hcalElems;
    block.associatedElements( iTrack,  linkData,
                              hcalElems,
                              reco::PFBlockElement::HCAL,
                              reco::PFBlock::LINKTEST_ALL ); 
    //MICHELE 
    //TEMPORARY SOLUTION FOR ELECTRON REJECTION IN PFTAU
    //COLINFEB16 
    // in case particle flow electrons are not reconstructed, 
    // the mva_e_pi of the charged hadron will be set to 1 
    // if a GSF element is associated to the current TRACK element
    // This information will be used in the electron rejection for tau ID.
    std::multimap<double,unsigned> gsfElems;
    if (usePFElectrons_ == false) {
      block.associatedElements( iTrack,  linkData,
				gsfElems ,
				reco::PFBlockElement::GSF );
    }
    //

    if(hcalElems.empty() && debug_) {
      cout<<"no hcal element connected to track "<<iTrack<<endl;
    }

    // will now loop on associated elements ... 
    // typedef std::multimap<double, unsigned>::iterator IE;
   
    bool hcalFound = false;

    if(debug_) 
      cout<<"now looping on elements associated to the track"<<endl;
    
    // ... first on associated ECAL elements
    // Check if there is still a free ECAL for this track
    for(IE ie = ecalElems.begin(); ie != ecalElems.end(); ++ie ) {
      
      unsigned index = ie->second;
      // Sanity checks and optional printout
      PFBlockElement::Type type = elements[index].type();
      if(debug_) {
	double  dist  = ie->first;
	//double chi2 = block.chi2(iTrack, index, linkData,
	//			 reco::PFBlock::LINKTEST_ALL);
        cout<<"\telement "<<elements[index]<<" linked with distance = "<< dist <<endl;
	if ( ! active[index] ) cout << "This ECAL is already used - skip it" << endl;      
      }
      assert( type == PFBlockElement::ECAL );

      // This ECAL is not free (taken by an electron?) - just skip it
      if ( ! active[index] ) continue;      

      // Flag ECAL clusters for which the corresponding track is not linked to an HCAL cluster

      //reco::PFClusterRef clusterRef = elements[index].clusterRef();
      //assert( !clusterRef.isNull() );
      if( !hcalElems.empty() && debug_)
	cout<<"\t\tat least one hcal element connected to the track."
	    <<" Sparing Ecal cluster for the hcal loop"<<endl; 

    } //loop print ecal elements

  
    // tracks which are not linked to an HCAL 
    // are reconstructed now. 

    if( hcalElems.empty() ) {

      if ( debug_ ) std::cout << "Now deals with tracks linked to no HCAL clusters" << std::endl; 
      // vector<unsigned> elementIndices;
      // elementIndices.push_back(iTrack); 

      // The track momentum
      double trackMomentum = elements[iTrack].trackRef()->p();      
      if ( debug_ ) std::cout << "Track momentum = " << trackMomentum << " GeV." << std::endl; 
      
      // Is it a "tight" muon ?
      bool thisIsAMuon = PFMuonAlgo::isMuon(elements[iTrack]);
      if ( thisIsAMuon ) trackMomentum = 0.;
      
      // Create a track candidate       
      unsigned tmpi = reconstructTrack( elements[iTrack] );
      active[iTrack] = false;
      
      double slopeEcal = 1.;
      bool connectedToEcal = false;
      unsigned iEcal = 99999;
      double calibEcal = 0.;
      double calibHcal = 0.;
      double totalEcal = thisIsAMuon ? -muonECAL_[0] : 0.;
      double totalHcal = 0.;
      std::multimap<double, unsigned> sortedTracks;

      // Loop over all ECAL linked clusters ordered by increasing distance.
      for(IE ie = ecalElems.begin(); ie != ecalElems.end(); ++ie ) {
	
	unsigned index = ie->second;
	PFBlockElement::Type type = elements[index].type();
	assert( type == PFBlockElement::ECAL );

	// Just skip clusters already taken
	if ( ! active[index] ) continue;

	// Just skip this cluster if it's closer to another track
	block.associatedElements( index,  linkData,
				  sortedTracks,
				  reco::PFBlockElement::TRACK,
				  reco::PFBlock::LINKTEST_ALL );
	if ( sortedTracks.begin()->second != iTrack ) continue;  
	
	// The corresponding PFCluster and energy
	reco::PFClusterRef clusterRef = elements[index].clusterRef();
	assert( !clusterRef.isNull() );

	if ( debug_ ) {
	  double dist = ie->first;
	  std::cout <<"Ecal cluster with raw energy = " << clusterRef->energy() 
		    <<" linked with distance = " << dist << std::endl;
	}
	/*
	double dist = ie->first;
	if ( !connectedToEcal && dist > 0.1 ) {
	  std::cout << "Warning - first ECAL cluster at a distance of " << dist << "from the closest track!" << std::endl;
	  cout<<"\telement "<<elements[index]<<" linked with distance = "<< dist <<endl;
	  cout<<"\telement "<<elements[iTrack]<<" linked with distance = "<< dist <<endl;
	}
	*/ 

	// Check the presence of preshower clusters in the vicinity
	// Preshower cluster closer to another ECAL cluster are ignored.
	vector<double> ps1Ene(1,static_cast<double>(0.));
	associatePSClusters(index, reco::PFBlockElement::PS1, block, elements, linkData, active, ps1Ene);
	vector<double> ps2Ene(1,static_cast<double>(0.));
	associatePSClusters(index, reco::PFBlockElement::PS2, block, elements, linkData, active, ps2Ene);
	
	// Get the energy calibrated (for photons)
	double ecalEnergy = calibration_->energyEm(*clusterRef,ps1Ene,ps2Ene);
	if ( debug_ )
	  std::cout << "Corrected ECAL(+PS) energy = " << ecalEnergy << std::endl;

	// Since the electrons were found beforehand, this track must be a hadron. Calibrate 
	// the energy under the hadron hypothesis.
	totalEcal += ecalEnergy;
	totalHcal = 0.;
	double previousCalibEcal = calibEcal;
	double previousSlopeEcal = slopeEcal;
	calibEcal = std::max(totalEcal,0.);
	calibHcal = 0.;
	if ( newCalib_ == 1 ) { 
	  // Warning ! This function changed the value of calibEcal and calibHcal
	  clusterCalibration_->
	    getCalibratedEnergyEmbedAInHcal(calibEcal, calibHcal,
					    clusterRef->positionREP().Eta(),
					    clusterRef->positionREP().Phi());
	  if ( totalEcal > 0. ) slopeEcal = calibEcal/totalEcal;
	} else if ( newCalib_ == 0 ) { 
	  // So called "Colin's calibration" - done with FAMOS in its early times
	  slopeEcal = calibration_->paramECALplusHCAL_slopeECAL();
	  calibEcal *= slopeEcal; 
	} else {
	  // Here enters the most recent calibration 
	  // The calibration with E+H or with E only is actually very similar.
	  // but it certainly deserves a special look (PJ, 24-Feb-2009)
	  calibration_->energyEmHad(trackMomentum,calibEcal,calibHcal,
	  			    clusterRef->positionREP().Eta(),
	  			    clusterRef->positionREP().Phi());
	  if ( totalEcal > 0.) slopeEcal = calibEcal/totalEcal;
	}

	if ( debug_ )
	  std::cout << "The total calibrated energy so far amounts to = " << calibEcal << std::endl;
	

	// Stop the loop when adding more ECAL clusters ruins the compatibility
	if ( connectedToEcal && calibEcal - trackMomentum >= 0. ) {
	// if ( connectedToEcal && calibEcal - trackMomentum >=
	//     nSigmaECAL_*neutralHadronEnergyResolution(trackMomentum,clusterRef->positionREP().Eta())  ) { 
	  calibEcal = previousCalibEcal;
	  slopeEcal = previousSlopeEcal;

	  // Turn this last cluster in a photon 
	  // (The PS clusters are already locked in "associatePSClusters")
	  active[index] = false;
	  unsigned tmpi = reconstructCluster( *clusterRef, ecalEnergy ); 
	  (*pfCandidates_)[tmpi].setEcalEnergy( ecalEnergy );
	  (*pfCandidates_)[tmpi].setHcalEnergy( 0 );
	  (*pfCandidates_)[tmpi].setPs1Energy( ps1Ene[0] );
	  (*pfCandidates_)[tmpi].setPs2Energy( ps2Ene[0] );
	  (*pfCandidates_)[tmpi].addElementInBlock( blockref, index );
	  break;
	}

	// Lock used clusters.
	connectedToEcal = true;
	iEcal = index;
	active[index] = false;

      } // Loop ecal elements

      

      (*pfCandidates_)[tmpi].setEcalEnergy( calibEcal );
      (*pfCandidates_)[tmpi].setHcalEnergy( 0 );
      (*pfCandidates_)[tmpi].setPs1Energy( 0 );
      (*pfCandidates_)[tmpi].setPs2Energy( 0 );
      (*pfCandidates_)[tmpi].addElementInBlock( blockref, iTrack ); 

      // Create a photon if the ecal energy is too large
      if( connectedToEcal ) {
	reco::PFClusterRef pivotalRef = elements[iEcal].clusterRef();
	(*pfCandidates_)[tmpi].addElementInBlock( blockref, iEcal ); 

	double neutralEnergy = calibEcal - trackMomentum;

	// Check if there are other tracks linked to that ECAL
	for(IE ie = sortedTracks.begin(); ie != sortedTracks.end() && neutralEnergy > 0; ++ie ) {
	  unsigned jTrack = ie->second;
	  
	  // The loop is on the other tracks !
	  if ( jTrack == iTrack ) continue;

	  // Check if the ECAL closest to this track is the current ECAL
	  // Otherwise ignore this track in the neutral energy determination
	  std::multimap<double, unsigned> sortedECAL;
	  block.associatedElements( jTrack,  linkData,
				    sortedECAL,
				    reco::PFBlockElement::ECAL,
				    reco::PFBlock::LINKTEST_ALL );
	  if ( sortedECAL.begin()->second != iEcal ) continue;
	  
	  // Check if this track is also linked to an HCAL 
	  // (in which case it goes to the next loop and is ignored here)
	  std::multimap<double, unsigned> sortedHCAL;
	  block.associatedElements( jTrack,  linkData,
				    sortedHCAL,
				    reco::PFBlockElement::HCAL,
				    reco::PFBlock::LINKTEST_ALL );
	  if ( sortedHCAL.size() ) continue;
	  
	  // Otherwise, subtract this track momentum from the ECAL energy 
	  reco::TrackRef trackRef = elements[jTrack].trackRef();
	  neutralEnergy -= trackRef->p();

	} // End other tracks
	
	// Add a photon is the energy excess is large enough
	if ( neutralEnergy > nSigmaECAL_*neutralHadronEnergyResolution(trackMomentum,pivotalRef->positionREP().Eta()) ) { 
	  neutralEnergy /= slopeEcal;
	  unsigned tmpj = reconstructCluster( *pivotalRef, neutralEnergy ); 
	  (*pfCandidates_)[tmpj].setEcalEnergy( neutralEnergy );
	  (*pfCandidates_)[tmpj].setHcalEnergy( 0. );
	  (*pfCandidates_)[tmpj].setPs1Energy( -1 );
	  (*pfCandidates_)[tmpj].setPs2Energy( -1 );
	  (*pfCandidates_)[tmpj].addElementInBlock(blockref, iEcal);
	  (*pfCandidates_)[tmpj].addElementInBlock( blockref, iTrack ); 
	} // End neutral energy
      } // End connected ECAL
      
      //MICHELE 
      //TEMPORARY SOLUTION FOR ELECTRON REJECTION IN PFTAU
      //COLINFEB16 
      // in case particle flow electrons are not reconstructed, 
      // the mva_e_pi of the charged hadron will be set to 1 
      // if a GSF element is associated to the current TRACK element
      if (usePFElectrons_ == false) 
	(*pfCandidates_)[tmpi].set_mva_e_pi(gsfElems.size()>0);
    }


    //COLINFEB16 
    // In case several HCAL elements are linked to this track, 
    // unlinking all of them except the closest. 
    // BUG?? currently forgetting to unlink the distance link!!
    for(IE ie = hcalElems.begin(); ie != hcalElems.end(); ++ie ) {
      
      // double   chi2  = ie->first;
      unsigned index = ie->second;
      
      PFBlockElement::Type type = elements[index].type();

      if(debug_) {
	double chi2 = block.chi2(iTrack,index,linkData,reco::PFBlock::LINKTEST_ALL);
        cout<<"\telement "<<elements[index]<<" linked with chi2 "<<chi2<<endl;
      }

      assert( type == PFBlockElement::HCAL );
      
      // all hcal clusters except the closest 
      // will be unlinked from the track
      if( !hcalFound ) { // closest hcal
        if(debug_) 
          cout<<"\t\tclosest hcal cluster, doing nothing"<<endl;
        
        hcalFound = true;
        
        // active[index] = false;
        // hcalUsed.push_back( index );
      }
      else { // other associated hcal
        // unlink from the track
        if(debug_) 
          cout<<"\t\tsecondary hcal cluster. unlinking"<<endl;
        //block.setLink( iTrack, index,  -1, -1, linkData );
        //Note ALex: a fonction should be implemented to
        //force link to -1 for all possible tests
        //For the moment (15/11/2007) only 2 possible
        //tests: CHI2 and RECHIT.

	//COLINFEB16 
	// There is no LINKTEST_DISTANCE
	// the true distance and the chi2 are computed for the 2 linktests,
	// which are LINKTEST_CHI2 and LINKTEST_RECHIT.
	// and those 2 values are the same, whatever the type of link...
	// we should review the PFBlock dataformat.
        block.setLink( iTrack, index,  -1, -1, linkData, 
                       PFBlock::LINKTEST_CHI2 );
        block.setLink( iTrack, index,  -1, -1, linkData,
                       PFBlock::LINKTEST_RECHIT );
      }
    } //loop hcal elements   
  } // end of loop 1 on elements iEle of any type



  // deal with HF. 
  if( !(hfEmIs.empty() &&  hfHadIs.empty() ) ) {
    // there is at least one HF element in this block. 
    // so all elements must be HF.
    assert( hfEmIs.size() + hfHadIs.size() == elements.size() );

    if( elements.size() == 1 ) {
      //Auguste: HAD-only calibration here
      reco::PFClusterRef clusterRef = elements[0].clusterRef();
      double energyHF = 0.;
      unsigned tmpi = 0;
      switch( clusterRef->layer() ) {
      case PFLayer::HF_EM:
	// do EM-only calibration here
	energyHF = clusterRef->energy();
	tmpi = reconstructCluster( *clusterRef, energyHF );     
	(*pfCandidates_)[tmpi].setEcalEnergy( energyHF );
	(*pfCandidates_)[tmpi].setHcalEnergy( 0.);
	(*pfCandidates_)[tmpi].setPs1Energy( -1 );
	(*pfCandidates_)[tmpi].setPs2Energy( -1 );
	(*pfCandidates_)[tmpi].addElementInBlock( blockref, hfEmIs[0] );
	//std::cout << "HF EM alone ! " << energyHF << std::endl;
	break;
      case PFLayer::HF_HAD:
	// do HAD-only calibration here
	energyHF = clusterRef->energy();
	tmpi = reconstructCluster( *clusterRef, energyHF );     
	(*pfCandidates_)[tmpi].setHcalEnergy( energyHF );
	(*pfCandidates_)[tmpi].setEcalEnergy( 0.);
	(*pfCandidates_)[tmpi].setPs1Energy( -1 );
	(*pfCandidates_)[tmpi].setPs2Energy( -1 );
	(*pfCandidates_)[tmpi].addElementInBlock( blockref, hfHadIs[0] );
	//std::cout << "HF Had alone ! " << energyHF << std::endl;
	break;
      default:
	assert(0);
      }
    }
    else if(  elements.size() == 2 ) {
      //Auguste: EM + HAD calibration here
      reco::PFClusterRef c0 = elements[0].clusterRef();
      reco::PFClusterRef c1 = elements[1].clusterRef();
      // 2 HF elements. Must be in each layer.
      reco::PFClusterRef cem = c1;
      reco::PFClusterRef chad = c0;
      
      if( c0->layer()== PFLayer::HF_EM ) {
	if ( c1->layer()== PFLayer::HF_HAD ) {
	  cem = c0; chad = c1;
	} else { 
	  assert(0);
	}
	// do EM+HAD calibration here
	double energyHfEm = cem->energy();
	double energyHfHad = chad->energy();
	unsigned tmpi = reconstructCluster( *chad, energyHfEm+energyHfHad );     
	(*pfCandidates_)[tmpi].setEcalEnergy( energyHfEm );
	(*pfCandidates_)[tmpi].setHcalEnergy( energyHfHad);
	(*pfCandidates_)[tmpi].setPs1Energy( -1 );
	(*pfCandidates_)[tmpi].setPs2Energy( -1 );
	(*pfCandidates_)[tmpi].addElementInBlock( blockref, hfEmIs[0] );
	(*pfCandidates_)[tmpi].addElementInBlock( blockref, hfHadIs[0] );
	//std::cout << "HF EM+HAD found ! " << energyHfEm << " " << energyHfHad << std::endl;

      }
      else {
	cerr<<"Warning: 2 elements, but not 1 HFEM and 1 HFHAD"<<endl;
	cerr<<block<<endl;
// 	assert( c1->layer()== PFLayer::HF_EM &&
// 		c0->layer()== PFLayer::HF_HAD );
      }      
    }
    else {
      // 1 HF element in the block, 
      // but number of elements not equal to 1 or 2 
      cerr<<"Warning: HF, but n elem different from 1 or 2"<<endl;
      cerr<<block<<endl;
//       assert(0);
//       cerr<<"not ready for navigation in the HF!"<<endl;
    }
  }



  if(debug_) { 
    cout<<endl;
    cout<<endl<<"--------------- loop hcal ---------------------"<<endl;
  }
  
  // track information:
  // trackRef
  // ecal energy
  // hcal energy
  // rescale 

  for(unsigned i=0; i<hcalIs.size(); i++) {
    
    unsigned iHcal= hcalIs[i];    
    PFBlockElement::Type type = elements[iHcal].type();

    assert( type == PFBlockElement::HCAL );

    if(debug_) cout<<endl<<elements[iHcal]<<endl;

    // vector<unsigned> elementIndices;
    // elementIndices.push_back(iHcal);

    // associated tracks
    std::multimap<double, unsigned> sortedTracks;
    block.associatedElements( iHcal,  linkData,
			      sortedTracks,
			      reco::PFBlockElement::TRACK,
			      reco::PFBlock::LINKTEST_ALL );
   
    std::map< unsigned, std::pair<double, unsigned> > associatedEcals;

    std::map< unsigned, std::pair<double, double> > associatedPSs;
 
    std::multimap<double, std::pair<unsigned,bool> > associatedTracks;
    
    // A temporary maps for ECAL satellite clusters
    std::multimap<double,std::pair<unsigned,double> > ecalSatellites;
    std::pair<unsigned,double> fakeSatellite = make_pair(iHcal,0.);
    ecalSatellites.insert( make_pair(-1., fakeSatellite) );

    PFClusterRef hclusterref = elements[iHcal].clusterRef();
    assert(!hclusterref.isNull() );


    //if there is no track attached to that HCAL, then do not
    //reconstruct an HCAL alone, check if it can be recovered 
    //first
    
    // if no associated tracks, create a neutral hadron
    //if(sortedTracks.empty() ) {

    if( sortedTracks.empty() ) {
      if(debug_) 
        cout<<"\tno associated tracks, keep for later"<<endl;
      continue;
    }

    // Lock the HCAL cluster
    active[iHcal] = false;


    // in the following, tracks are associated to this hcal cluster.
    // will look for an excess of energy in the calorimeters w/r to 
    // the charged energy, and turn this excess into a neutral hadron or 
    // a photon.

    if(debug_) cout<<"\t"<<sortedTracks.size()<<" associated tracks:"<<endl;

    double totalChargedMomentum = 0;
    double sumpError2 = 0.;
    double totalEcal = 0.;
    double totalHcal = hclusterref->energy();
    vector<double> hcalP;
    vector<double> hcalDP;
    vector<unsigned> tkIs;
    double maxDPovP = -9999.;
    
    vector< unsigned > chargedHadronsIndices;
    double mergedNeutralHadronEnergy = 0;
    double muonHCALEnergy = 0.;
    double muonECALEnergy = 0.;
    double muonHCALError = 0.;
    double muonECALError = 0.;
    unsigned nMuons = 0; 

    // Loop over all tracks associated to this HCAL cluster
    for(IE ie = sortedTracks.begin(); ie != sortedTracks.end(); ++ie ) {

      unsigned iTrack = ie->second;

      // Check the track has not already been used (e.g., in electrons, conversions...)
      if ( !active[iTrack] ) continue;
      // Sanity check 1
      PFBlockElement::Type type = elements[iTrack].type();
      assert( type == reco::PFBlockElement::TRACK );
      // Sanity check 2
      reco::TrackRef trackRef = elements[iTrack].trackRef();
      assert( !trackRef.isNull() ); 

      // Create a PF Candidate right away if the track is a tight muon
      reco::MuonRef muonRef = elements[iTrack].muonRef();
      bool thisIsAMuon = PFMuonAlgo::isMuon(elements[iTrack]);
      bool thisIsALooseMuon = PFMuonAlgo::isLooseMuon(elements[iTrack]);
      if ( thisIsAMuon ) { 
	if ( debug_ ) { 
	  std::cout << "\t\tThis track is identified as a muon - remove it from the stack" << std::endl;
	  std::cout << "\t\t" << elements[iTrack] << std::endl;
	}
	// Estimate of the energy deposit & resolution in the calorimeters
	muonHCALEnergy += muonHCAL_[0];
	muonHCALError += muonHCAL_[1]*muonHCAL_[1];
	muonECALEnergy += muonECAL_[0];
	muonECALError += muonECAL_[1]*muonECAL_[1];	
	// Create a muon.
	unsigned tmpi = reconstructTrack( elements[iTrack] );
	(*pfCandidates_)[tmpi].addElementInBlock( blockref, iTrack );
	(*pfCandidates_)[tmpi].addElementInBlock( blockref, iHcal );
	// Remove it from the stack
	active[iTrack] = false;
	// Go to next track
	continue;
      }
 
      // PJ - 01-Mar-09 : Do we still need this temporary solution ?
      //MICHELE 
      //TEMPORARY SOLUTION FOR ELECTRON REJECTION IN PFTAU

      /*						  
      if (usePFElectrons_ == false) {
	std::multimap<double,unsigned> gsfElems;
	block.associatedElements( iTrack,  linkData,
				  gsfElems ,
				  reco::PFBlockElement::GSF );
	currentChargedHadron.set_mva_e_pi(gsfElems.size()>0);
      }
      */
      //

      if(debug_) cout<<"\t\t"<<elements[iTrack]<<endl;

      // introduce  tracking errors 
      double trackMomentum = trackRef->p();
      totalChargedMomentum += trackMomentum; 

      // If the track is not a tight muon, but still resembles a muon
      // keep it for later for blocks with too large a charged energy
      if ( thisIsALooseMuon && !thisIsAMuon ) nMuons += 1;	

      // ... and keep anyway the pt error for possible fake rejection
      double DPt = trackRef->ptError();
      std::pair<unsigned,bool> tkmuon = make_pair(iTrack,thisIsALooseMuon);
      associatedTracks.insert( make_pair(-DPt, tkmuon) );     
 
      // Also keep the total track momentum error for comparison with the calo energy
      double Dp = trackRef->qoverpError()*trackMomentum*trackMomentum;
      sumpError2 += Dp*Dp;

      // look for ECAL elements associated to iTrack (associated to iHcal)
      std::multimap<double, unsigned> sortedEcals;
      block.associatedElements( iTrack,  linkData,
                                sortedEcals,
				reco::PFBlockElement::ECAL,
				reco::PFBlock::LINKTEST_ALL );
    
      if(debug_) cout<<"\t\t\tnumber of Ecal elements linked to this track: "
                     <<sortedEcals.size()<<endl;
      
      bool connectedToEcal = false; // Will become true if there is at least one ECAL cluster connected
      if( !sortedEcals.empty() ) 
        { // start case: at least one ecal element associated to iTrack
	  
	  // Loop over all connected ECAL clusters
	  for ( IE iec=sortedEcals.begin(); 
		iec!=sortedEcals.end(); ++iec ) {

	    unsigned iEcal = iec->second;
	    double dist = iec->first;
	    
	    // Ignore ECAL cluters already used
	    if( !active[iEcal] ) { 
	      if(debug_) cout<<"\t\tcluster locked"<<endl;          
	      continue;
	    }

	    // Sanity checks
	    PFBlockElement::Type type = elements[ iEcal ].type();
	    assert( type == PFBlockElement::ECAL );
	    PFClusterRef eclusterref = elements[iEcal].clusterRef();
	    assert(!eclusterref.isNull() );
	      
	    // Check if this ECAL is not closer to another track - ignore it in that case
	    std::multimap<double, unsigned> sortedTracksEcal;
	    block.associatedElements( iEcal,  linkData,
				      sortedTracksEcal,
				      reco::PFBlockElement::TRACK,
				      reco::PFBlock::LINKTEST_ALL );
	    unsigned jTrack = sortedTracksEcal.begin()->second;
	    if ( jTrack != iTrack ) continue; 

	    double chi2Ecal = block.chi2(jTrack,iEcal,linkData,
					 reco::PFBlock::LINKTEST_ALL);
	    // totalEcal += eclusterref->energy();
	    // float ecalEnergyUnCalibrated = eclusterref->energy();
	    //std::cout << "Ecal Uncalibrated " << ecalEnergyUnCalibrated << std::endl;
	    
	    // check the presence of preshower clusters in the vicinity
	    vector<double> ps1Ene(1,static_cast<double>(0.));
	    associatePSClusters(iEcal, reco::PFBlockElement::PS1, 
				block, elements, linkData, active, 
				ps1Ene);
	    vector<double> ps2Ene(1,static_cast<double>(0.));
	    associatePSClusters(iEcal, reco::PFBlockElement::PS2, 
				block, elements, linkData, active, 
				ps2Ene); 
	    std::pair<double,double> psEnes = make_pair(ps1Ene[0],
							ps2Ene[0]);
	    associatedPSs.insert(make_pair(iEcal,psEnes));
	      
	    // Calibrate the ECAL energy for photons
	    float ecalEnergyCalibrated = calibration_->energyEm(*eclusterref,ps1Ene,ps2Ene);

	    if ( !connectedToEcal ) { // This is the closest ECAL cluster - will add its energy later
	      
	      if(debug_) cout<<"\t\t\tclosest: "
			     <<elements[iEcal]<<endl;
	      
	      connectedToEcal = true;
	      // PJ 1st-April-09 : To be done somewhere !!! (Add to comment it, but it is needed)
	      // currentChargedHadron.addElementInBlock( blockref, iEcal );
       
	      std::pair<double, unsigned> associatedEcal 
		= make_pair( chi2Ecal, iEcal );
	      associatedEcals.insert( make_pair(iTrack, associatedEcal) );
	      std::pair<unsigned,double> satellite = make_pair(iEcal,ecalEnergyCalibrated);
	      ecalSatellites.insert( make_pair(-1., satellite) );

	    } else { // Keep satellite clusters for later
	      
	      std::pair<unsigned,double> satellite = make_pair(iEcal,ecalEnergyCalibrated);
	      ecalSatellites.insert( make_pair(dist, satellite) );
	      
	    }

	  } // End loop ecal associated to iTrack
	} // end case: at least one ecal element associated to iTrack
      
    } // end loop on tracks associated to hcal element iHcal
    
    // test compatibility between calo and tracker. //////////////

    double caloEnergy = 0.;
    double slopeEcal = 1.0;
    double calibEcal = 0.;
    double calibHcal = 0.;

    // Determine the expected calo resolution from the total charged momentum
    double Caloresolution = neutralHadronEnergyResolution( totalChargedMomentum, hclusterref->positionREP().Eta());    
    Caloresolution *= totalChargedMomentum;
    // Account for muons
    Caloresolution = std::sqrt(Caloresolution*Caloresolution + muonHCALError + muonECALError);
    totalEcal -= std::min(totalEcal,muonECALEnergy);
    totalHcal -= std::min(totalEcal,muonHCALEnergy);

    // Loop over all ECAL satellites, starting for the closest to the various tracks
    // and adding other satellites until saturation of the total track momentum
    // Note : for code simplicity, the first element of the loop is the HCAL cluster
    // with 0 energy in the ECAL
    for ( IS is = ecalSatellites.begin(); is != ecalSatellites.end(); ++is ) { 

      // Add the energy of this ECAL cluster
      double previousCalibEcal = calibEcal;
      double previousCalibHcal = calibHcal;
      double previousCaloEnergy = caloEnergy;
      double previousSlopeEcal = slopeEcal;
      //
      totalEcal += is->second.second;
      calibEcal = std::max(0.,totalEcal);
      calibHcal = std::max(0.,totalHcal);

      // Calibrate ECAL and HCAL energy under the hadron hypothesis.
      if ( newCalib_ == 1) {
	clusterCalibration_->
	  getCalibratedEnergyEmbedAInHcal(calibEcal, calibHcal,
					  hclusterref->positionREP().Eta(),
					  hclusterref->positionREP().Phi());
	caloEnergy = calibEcal+calibHcal;
	if ( totalEcal > 0.) slopeEcal = calibEcal/totalEcal;
      } else if ( newCalib_ == 0 ) { 
	if( totalEcal>0) { 
	  caloEnergy = calibration_->energyEmHad( totalEcal, totalHcal );
	  slopeEcal = calibration_->paramECALplusHCAL_slopeECAL();
	  calibEcal = totalEcal * slopeEcal;
	  calibHcal = caloEnergy-calibEcal;
	} else { 
	  caloEnergy = calibration_->energyHad( totalHcal );
	  calibEcal = totalEcal;
	  calibHcal = caloEnergy-calibEcal;
	}
      } else { 
	calibration_->energyEmHad(totalChargedMomentum,calibEcal,calibHcal,
				  hclusterref->positionREP().Eta(),
				  hclusterref->positionREP().Phi());
	caloEnergy = calibEcal+calibHcal;
	if ( totalEcal > 0.) slopeEcal = calibEcal/totalEcal;
      }
      
      // Continue looping until all closest clusters are exhausted and as long as
      // the calorimetric energy does not saturate the total momentum.
      if ( is->first < 0. || caloEnergy - totalChargedMomentum <= 0. ) { 
	if(debug_) cout<<"\t\t\tactive, adding "<<is->second.second
		       <<" to ECAL energy, and locking"<<endl;
	active[is->second.first] = false;
	continue;
      }

      // Otherwise, do not consider the last cluster examined and exit.
      // active[is->second.first] = true;
      totalEcal -= is->second.second;
      calibEcal = previousCalibEcal;
      calibHcal = previousCalibHcal;
      slopeEcal = previousSlopeEcal;
      caloEnergy = previousCaloEnergy;

      break;

    }

    // Sanity check !
    assert(caloEnergy>=0);


    // All other satellites are treated as photons.
    for ( IS is = ecalSatellites.begin(); is != ecalSatellites.end(); ++is ) { 
      
      // Ignore satellites already taken
      unsigned iEcal = is->second.first;
      if ( !active[iEcal] ) continue;

      // Sanity checks again (well not useful, this time!)
      PFBlockElement::Type type = elements[ iEcal ].type();
      assert( type == PFBlockElement::ECAL );
      PFClusterRef eclusterref = elements[iEcal].clusterRef();
      assert(!eclusterref.isNull() );

      // Lock the cluster
      active[iEcal] = false;
      
      // Create a photon
      unsigned tmpi = reconstructCluster( *eclusterref, is->second.second ); 
      (*pfCandidates_)[tmpi].setEcalEnergy( is->second.second );
      (*pfCandidates_)[tmpi].setHcalEnergy( 0 );
      (*pfCandidates_)[tmpi].setPs1Energy( associatedPSs[iEcal].first );
      (*pfCandidates_)[tmpi].setPs2Energy( associatedPSs[iEcal].second );
      (*pfCandidates_)[tmpi].addElementInBlock( blockref, iEcal );

    }

    // And now check for hadronic energy excess...

    //colin: resolution should be measured on the ecal+hcal case. 
    // however, the result will be close. 
    // double Caloresolution = neutralHadronEnergyResolution( caloEnergy );
    // Caloresolution *= caloEnergy;
    // PJ The resolution is on the expected charged calo energy !
    //double Caloresolution = neutralHadronEnergyResolution( totalChargedMomentum, hclusterref->positionREP().Eta());
    //Caloresolution *= totalChargedMomentum;
    // that of the charged particles linked to the cluster!
    double TotalError = sqrt(sumpError2 + Caloresolution*Caloresolution);

    /* */
    ////////////////////// TRACKER MUCH LARGER THAN CALO /////////////////////////
    if ( totalChargedMomentum - caloEnergy > nSigmaTRACK_*Caloresolution ) { 
      /* 
      cout<<"\tCompare Calo Energy to total charged momentum "<<endl;
      cout<<"\t\tsum p    = "<<totalChargedMomentum<<" +- "<<sqrt(sumpError2)<<endl;
      cout<<"\t\tsum ecal = "<<totalEcal<<endl;
      cout<<"\t\tsum hcal = "<<totalHcal<<endl;
      cout<<"\t\t => Calo Energy = "<<caloEnergy<<" +- "<<Caloresolution<<endl;
      cout<<"\t\t => Calo Energy- total charged momentum = "
         <<caloEnergy-totalChargedMomentum<<" +- "<<TotalError<<endl;
      cout<<endl;
      cout << "\t\tNumber/momentum of muons found in the block : " << nMuons << std::endl;
      */
      // First consider loose muons
      if ( nMuons > 0 ) { 
	for ( IT it = associatedTracks.begin(); it != associatedTracks.end(); ++it ) { 
	  unsigned iTrack = it->second.first;
	  // Only active tracks
	  if ( !active[iTrack] ) continue;
	  // Only muons
	  if ( !it->second.second ) continue;
	  // Consider only muons with pt>20 (maybe put this selection in isLooseMuon?)
	  double trackMomentum = elements[it->second.first].trackRef()->p();
	  double muonMomentum = elements[it->second.first].muonRef()->combinedMuon()->p();
	  double staMomentum = elements[it->second.first].muonRef()->standAloneMuon()->p();
	  double trackPt = elements[it->second.first].trackRef()->pt();
	  double muonPt = elements[it->second.first].muonRef()->combinedMuon()->pt();
	  double staPt = elements[it->second.first].muonRef()->standAloneMuon()->pt();
	  if ( staPt < 20. && muonPt < 20. && trackPt < 20. ) continue;	  
	  // Add the muon to the PF candidate list
	  unsigned tmpi = reconstructTrack( elements[iTrack] );
	  (*pfCandidates_)[tmpi].addElementInBlock( blockref, iTrack );
	  (*pfCandidates_)[tmpi].addElementInBlock( blockref, iHcal );
	  reco::PFCandidate::ParticleType particleType = reco::PFCandidate::mu;
	  (*pfCandidates_)[tmpi].setParticleType(particleType);
	  // Take the best pt measurement
	  double trackPtError = elements[it->second.first].trackRef()->ptError();
	  double muonPtError = elements[it->second.first].muonRef()->combinedMuon()->ptError();
	  double staPtError = elements[it->second.first].muonRef()->standAloneMuon()->ptError();
	  double globalCorr = 1;
	  if ( muonPt > 20. && muonPtError < trackPtError && muonPtError < staPtError ) globalCorr = muonMomentum/trackMomentum;
	  if ( staPt > 20. && staPtError < trackPtError && staPtError < muonPtError ) globalCorr = staMomentum/trackMomentum;
	  (*pfCandidates_)[tmpi].rescaleMomentum(globalCorr);
	  if (debug_) std::cout << "\tElement  " << elements[iTrack] << std::endl 
				<< "PFAlgo: particle type set to muon" << std::endl; 
	  // Remove it from the block
	  totalChargedMomentum -= trackMomentum;
	  // Update the calo energies
	  calibEcal -= std::min(calibEcal,muonECAL_[0]*calibEcal/totalEcal);
	  calibHcal -= std::min(calibHcal,muonHCAL_[0]*calibHcal/totalHcal);
	  totalEcal -= std::min(totalEcal,muonECAL_[0]);
	  totalHcal -= std::min(totalHcal,muonHCAL_[0]);
	  caloEnergy = calibEcal+calibHcal;
	  muonHCALEnergy += muonHCAL_[0];
	  muonHCALError += muonHCAL_[1]*muonHCAL_[1];
	  muonECALEnergy += muonECAL_[0];
	  muonECALError += muonECAL_[1]*muonECAL_[1];
	  active[iTrack] = false;
	  // Stop the loop whenever enough muons are removed
	  if ( totalChargedMomentum < caloEnergy ) break;	
	}
	// New calo resolution.
	Caloresolution = neutralHadronEnergyResolution( totalChargedMomentum, hclusterref->positionREP().Eta());    
	Caloresolution *= totalChargedMomentum;
	Caloresolution = std::sqrt(Caloresolution*Caloresolution + muonHCALError + muonECALError);
      }
    }

    // Second consider bad tracks (if still needed after muon removal)
    unsigned corrTrack = 10000000;
    double corrFact = 1.;
    if ( totalChargedMomentum - caloEnergy > nSigmaTRACK_*Caloresolution ) { 
      for ( IT it = associatedTracks.begin(); it != associatedTracks.end(); ++it ) { 
	unsigned iTrack = it->second.first;
	// Only active tracks
	if ( !active[iTrack] ) continue;
	// Consider only bad tracks
	if ( fabs(it->first) < ptError_ ) break;
	// What would become the block charged momentum if this track were removed
	double wouldBeTotalChargedMomentum = 
	  totalChargedMomentum - elements[it->second.first].trackRef()->p();
	// Reject worst tracks, as long as the total charged momentum 
	// is larger than the calo energy
	if ( wouldBeTotalChargedMomentum > caloEnergy ) { 
	  active[iTrack] = false;
	  totalChargedMomentum = wouldBeTotalChargedMomentum;
	  if ( debug_ )
	    std::cout << "\tElement  " << elements[iTrack] 
		      << " rejected (DPt = " << -it->first << " GeV/c)" << std::endl;
	// Just rescale the nth worst track momentum to equalize the calo energy
	} else {	
	  corrTrack = iTrack;
	  corrFact = (caloEnergy - wouldBeTotalChargedMomentum)/elements[it->second.first].trackRef()->p();
	  totalChargedMomentum -= elements[it->second.first].trackRef()->p()*(1.-corrFact);
	  if ( debug_ ) 
	    std::cout << "\tElement  " << elements[iTrack] 
		      << " rescaled by " << corrFact << std::endl;
	  break;
	}
      }
    }


    // New determination of the calo and track resolution avec track deletion/rescaling.
    Caloresolution = neutralHadronEnergyResolution( totalChargedMomentum, hclusterref->positionREP().Eta());    
    Caloresolution *= totalChargedMomentum;
    Caloresolution = std::sqrt(Caloresolution*Caloresolution + muonHCALError + muonECALError);

    // Make PF candidates with the remaining tracks in the block
    sumpError2 = 0.;
    for ( IT it = associatedTracks.begin(); it != associatedTracks.end(); ++it ) { 
      unsigned iTrack = it->second.first;
      if ( !active[iTrack] ) continue;
      reco::TrackRef trackRef = elements[iTrack].trackRef();
      double trackMomentum = trackRef->p();
      double Dp = trackRef->qoverpError()*trackMomentum*trackMomentum;
      unsigned tmpi = reconstructTrack( elements[iTrack] );
      (*pfCandidates_)[tmpi].addElementInBlock( blockref, iTrack );
      (*pfCandidates_)[tmpi].addElementInBlock( blockref, iHcal );
      if ( iTrack == corrTrack ) { 
	(*pfCandidates_)[tmpi].rescaleMomentum(corrFact);
	trackMomentum *= corrFact;
      }
      chargedHadronsIndices.push_back( tmpi );
      active[iTrack] = false;
      hcalP.push_back(trackMomentum);
      hcalDP.push_back(Dp);
      if (Dp/trackMomentum > maxDPovP) maxDPovP = Dp/trackMomentum;
      sumpError2 += Dp*Dp;
    }

    // The total uncertainty of the difference Calo-Track
    TotalError = sqrt(sumpError2 + Caloresolution*Caloresolution);

    if ( debug_ ) { 
      cout<<"\tCompare Calo Energy to total charged momentum "<<endl;
      cout<<"\t\tsum p    = "<<totalChargedMomentum<<" +- "<<sqrt(sumpError2)<<endl;
      cout<<"\t\tsum ecal = "<<totalEcal<<endl;
      cout<<"\t\tsum hcal = "<<totalHcal<<endl;
      cout<<"\t\t => Calo Energy = "<<caloEnergy<<" +- "<<Caloresolution<<endl;
      cout<<"\t\t => Calo Energy- total charged momentum = "
	  <<caloEnergy-totalChargedMomentum<<" +- "<<TotalError<<endl;
      cout<<endl;
    }

    /* */

    /////////////// TRACKER AND CALO COMPATIBLE  ////////////////
    double nsigma = nSigmaHCAL_;
    if ( abs(totalChargedMomentum-caloEnergy)<nsigma*TotalError ) {

      // deposited caloEnergy compatible with total charged momentum
      // if tracking errors are large take weighted average
      
      if(debug_) {
        cout<<"\t\tcase 1: COMPATIBLE "
            <<"|Calo Energy- total charged momentum| = "
            <<abs(caloEnergy-totalChargedMomentum)
            <<" < "<<nsigma<<" x "<<TotalError<<endl;
        if (maxDPovP < 0.1 )
          cout<<"\t\t\tmax DP/P = "<< maxDPovP 
              <<" less than 0.1: do nothing "<<endl;
        else 
          cout<<"\t\t\tmax DP/P = "<< maxDPovP 
              <<" >  0.1: take weighted averages "<<endl;
      }
      
      // if max DP/P < 10%  do nothing
      if (maxDPovP > 0.1) {
      
        // for each track associated to hcal
        //      int nrows = tkIs.size();
        int nrows = chargedHadronsIndices.size();
        TMatrixTSym<double> a (nrows);
        TVectorD b (nrows);
        TVectorD check(nrows);
        double sigma2E = Caloresolution*Caloresolution;
        for(int i=0; i<nrows; i++) {
          double sigma2i = hcalDP[i]*hcalDP[i];
          if (debug_){
            cout<<"\t\t\ttrack associated to hcal "<<i 
                <<" P = "<<hcalP[i]<<" +- "
                <<hcalDP[i]<<endl;
          }
          a(i,i) = 1./sigma2i + 1./sigma2E;
          b(i) = hcalP[i]/sigma2i + caloEnergy/sigma2E;
          for(int j=0; j<nrows; j++) {
            if (i==j) continue;
            a(i,j) = 1./sigma2E;
          } // end loop on j
        } // end loop on i
        
        // solve ax = b
        //if (debug_){// debug1
        //cout<<" matrix a(i,j) "<<endl;
        //a.Print();
        //cout<<" vector b(i) "<<endl;
        //b.Print();
        //} // debug1
        TDecompChol decomp(a);
        bool ok = false;
        TVectorD x = decomp.Solve(b, ok);
        // for each track create a PFCandidate track
        // with a momentum rescaled to weighted average
        if (ok) {
          for (int i=0; i<nrows; i++){
            //      unsigned iTrack = trackInfos[i].index;
            unsigned ich = chargedHadronsIndices[i];
            double rescaleFactor =  x(i)/hcalP[i];
            (*pfCandidates_)[ich].rescaleMomentum( rescaleFactor );

            if(debug_){
              cout<<"\t\t\told p "<<hcalP[i]
                  <<" new p "<<x(i) 
                  <<" rescale "<<rescaleFactor<<endl;
            }
          }
        }
        else {
          cerr<<"not ok!"<<endl;
          assert(0);
        }
      }
    }

    /////////////// NEUTRAL DETECTION  ////////////////
    else if( caloEnergy > totalChargedMomentum ) {
      
      //case 2: caloEnergy > totalChargedMomentum + nsigma*TotalError
      //there is an excess of energy in the calos
      //create a neutral hadron or a photon
        
      double eNeutralHadron = caloEnergy - totalChargedMomentum;
      double ePhoton = (caloEnergy - totalChargedMomentum) / slopeEcal;
        
      if(debug_) {
        if(!sortedTracks.empty() ){
          cout<<"\t\tcase 2: NEUTRAL DETECTION "
              <<caloEnergy<<" > "<<nsigma<<"x"<<TotalError
              <<" + "<<totalChargedMomentum<<endl;
          cout<<"\t\tneutral activity detected: "<<endl
              <<"\t\t\t           photon = "<<ePhoton<<endl
              <<"\t\t\tor neutral hadron = "<<eNeutralHadron<<endl;
            
          cout<<"\t\tphoton or hadron ?"<<endl;}
          
        if(sortedTracks.empty() )
          cout<<"\t\tno track -> hadron "<<endl;
        else 
          cout<<"\t\t"<<sortedTracks.size()
              <<"tracks -> loop and compute mva"<<endl;
      }
        
      // IE maxMvaIe = sortedTracks.end(); 
      reco::PFClusterRef maxMvaEcalRef;
      reco::PFClusterRef maxPSClusterRef;
      double maxMva = -PFCandidate::bigMva_;
      bool PSsignature = false;
      unsigned maxMvaiEcal= 9999;       
      unsigned maxPSiEcal= 9999;  
      // for each track associated to hcal: iterator IE ie :
        
      for(IE ie = sortedTracks.begin(); ie != sortedTracks.end(); ++ie ) {
          
        unsigned iTrack = ie->second;
          
        PFBlockElement::Type type = elements[iTrack].type();
        assert( type == reco::PFBlockElement::TRACK );
          
        reco::TrackRef trackRef = elements[iTrack].trackRef();
        assert( !trackRef.isNull() ); 
          
        std::map< unsigned, std::pair<double, unsigned> >::iterator 
          iae = associatedEcals.find(iTrack);
          
        if( iae == associatedEcals.end() ) continue;
          
        double chi2ECAL = iae->second.first;
        unsigned iEcal = iae->second.second;
          
        PFBlockElement::Type typeEcal = elements[iEcal].type();
        assert(typeEcal == reco::PFBlockElement::ECAL);
          
        reco::PFClusterRef clusterRef = elements[iEcal].clusterRef();
        assert( !clusterRef.isNull() );
          
        chi2ECAL_ = static_cast<float>( chi2ECAL );
        ptTrack_ = static_cast<float>( trackRef->pt() );
          
        float pTrack = static_cast<float>( trackRef->p() );
        float eECAL = static_cast<float>( clusterRef->energy() );
        eECALOverpTrack_ = eECAL / pTrack;
          
        double value = -PFAlgo::maxMvaCut_;
        if(mvaCut_ < PFAlgo::maxMvaCut_ ) {
          value = mergedPhotonsMVA_->EvaluateMVA( "MVA" );
          if( value>maxMva ) {
            maxMva = value;
            maxMvaEcalRef = clusterRef;
            maxMvaiEcal = iEcal;
          }
        }         
          
        if(debug_) {
          cout<<"\t\t"<<elements[iTrack]<<endl;
          cout<<"\t\t\tassociated ECAL "<<elements[iEcal]<<endl;
          cout<<"\t\t\tmva evaluation : "<<endl;
          cout<<"\t\t\t\teECALOverpTrack = "<<eECALOverpTrack_<<endl;
          cout<<"\t\t\t\tptTrack         = "<<ptTrack_<<endl;
          cout<<"\t\t\t\tchi2ECAL        = "<<chi2ECAL_<<endl;
          cout<<"\t\t\t ==>        mva = "<<value<<endl;
            
        }

        // look for PS1/PS2 elements associated to iEcal
        double PS1energy = associatedPSs[iEcal].first;
        double PS2energy = associatedPSs[iEcal].second;
	bool PS1found = PS1energy > PSCut_;
	bool PS2found = PS2energy > PSCut_;
	
        // evaluate PS signature
        if (PSsignature) continue;
        //PSsignature = PS1found&&PS2found;
        PSsignature = PS1found||PS2found;
        maxPSClusterRef = clusterRef;
        maxPSiEcal = iEcal;
        if (debug_) {
          cout<<" PS signature "<<PSsignature 
              <<" PS1 energy "<<PS1energy
              <<" PS2 energy "<<PS2energy<<endl;
        }
	
      }  // end loop on tracks associated to hcal: iterator IE ie
    
      std::vector<reco::PFClusterRef> pivotalClusterRef;
      std::vector<unsigned> iPivotal;
      std::vector<double> particleEnergy, ecalEnergy, hcalEnergy;
      // vector<unsigned> elementIndices;
      
      // double particleEnergy = -1;
      // float  ecalEnergy = 0;
      // float  hcalEnergy = 0;
      
      if( PSsignature) { //  overwrites MVA decision
        if(debug_) {
          cout<<" PS signature is true"
              <<", reconstruct a photon"<<endl;
        }
      
	if ( ePhoton < totalEcal || eNeutralHadron-calibEcal < 1E-10 ) { 
	  
	  particleEnergy.push_back(ePhoton);
	  ecalEnergy.push_back(ePhoton);
	  hcalEnergy.push_back(0.);
	  pivotalClusterRef.push_back(maxPSClusterRef);
	  iPivotal.push_back(maxPSiEcal);
	  
	} else { 
	  
	  // Assign the ECAL energy to the photon
	  pivotalClusterRef.push_back(maxMvaEcalRef);
	  particleEnergy.push_back(totalEcal);
	  ecalEnergy.push_back(totalEcal);
	  hcalEnergy.push_back(0.);
	  iPivotal.push_back(maxPSiEcal);
	  
	  // Assign the remaining excess to a neutral hadron
	  pivotalClusterRef.push_back(hclusterref);
	  particleEnergy.push_back(eNeutralHadron-calibEcal);
	  ecalEnergy.push_back(0.);
	  hcalEnergy.push_back(eNeutralHadron-calibEcal);
	  iPivotal.push_back(iHcal);
	  
	  mergedNeutralHadronEnergy = eNeutralHadron-calibEcal;
	  
	}
	
      } // end case: PS signature is true
      else { // start case PS signature is false
        if( !maxMvaEcalRef.isNull() && 
            maxMva > mvaCut_) { //  start case: maxMva > mvaCut_
	  
	  if ( ePhoton < totalEcal || eNeutralHadron-calibEcal < 1E-10 ) { // PJ !!! waiting for a better MVA ;-)
	    if(debug_) {
	      cout<<"\t\tmaxMva = "<<maxMva<<">"<< mvaCut_
		  <<", reconstruct a photon"<<endl;
	      
	    }
	    
	    // PJ BUG !
	    //if( ePhoton>totalEcal ) 
	    //  ePhoton = totalEcal;
	    
            pivotalClusterRef.push_back(maxMvaEcalRef);
	    particleEnergy.push_back(ePhoton);
	    ecalEnergy.push_back(ePhoton);
	    hcalEnergy.push_back(0.);
	    iPivotal.push_back(maxMvaiEcal);
	    
	  } else { 
	    
	    // Assign the ECAL energy to the photon
  	    pivotalClusterRef.push_back(maxMvaEcalRef);
	    particleEnergy.push_back(totalEcal);
	    ecalEnergy.push_back(totalEcal);
	    hcalEnergy.push_back(0.);
	    iPivotal.push_back(maxMvaiEcal);
	    
	    // Assign the remaining excess to a neutral hadron
            pivotalClusterRef.push_back(hclusterref);
	    particleEnergy.push_back(eNeutralHadron-calibEcal);
	    ecalEnergy.push_back(0.);
	    hcalEnergy.push_back(eNeutralHadron-calibEcal);
	    iPivotal.push_back(iHcal);
	    
	    mergedNeutralHadronEnergy = eNeutralHadron-calibEcal;
	    
	  }
	  
	}   // end case: maxMva > mvaCut_
        else {   //  start case : maxMva < mvaCut_
          if(debug_) {
            cout<<"\t\treconstruct a hadron"<<endl;
          }
          
	  // if ( totalEcal > 0. ) { 
	  if ( ePhoton < totalEcal || eNeutralHadron-calibEcal <= 1E-10 ) { // PJ !!! waiting for a better MVA ;-)
	    // Assign the ECAL energy to the photon
	    pivotalClusterRef.push_back(maxMvaEcalRef);
	    particleEnergy.push_back(ePhoton);
	    ecalEnergy.push_back(ePhoton);
	    hcalEnergy.push_back(0.);
	    iPivotal.push_back(maxMvaiEcal);
	    
	  } else { 
	    
	    // Assign the ECAL energy to the photon
	    if ( totalEcal > 0. ) { 
	      pivotalClusterRef.push_back(maxMvaEcalRef);
	      particleEnergy.push_back(totalEcal);
	      ecalEnergy.push_back(totalEcal);
	      hcalEnergy.push_back(0.);
	      iPivotal.push_back(maxMvaiEcal);
	    }
	    
	    // Assign the remaining excess to a neutral hadron
	    pivotalClusterRef.push_back(hclusterref);
	    particleEnergy.push_back(eNeutralHadron-calibEcal);
	    ecalEnergy.push_back(0.);
	    hcalEnergy.push_back(eNeutralHadron-calibEcal);
	    iPivotal.push_back(iHcal);
	    
	    /*
	      pivotalClusterRef.push_back(hclusterref);
	      particleEnergy.push_back(eNeutralHadron);
	      ecalEnergy.push_back(0.);
	      hcalEnergy.push_back(eNeutralHadron);
	      iPivotal.push_back(iHcal);
	    */
	    
	    // keep track globally of the merged neutral hadron
	    mergedNeutralHadronEnergy = eNeutralHadron-calibEcal;
	  } 
        }   // end case : maxMva < mvaCut_
      }   // end case PS signature is false 
    
      
      // reconstructing a merged neutral
      // the type of PFCandidate is known from the 
      // reference to the pivotal cluster. 
      
      for ( unsigned iPivot=0; iPivot<iPivotal.size(); ++iPivot ) { 
	
	if ( particleEnergy[iPivot] < 0. ) 
	  std::cout << "ALARM = Negative energy ! " 
		    << particleEnergy[iPivot] << std::endl;

	unsigned tmpi = reconstructCluster( *pivotalClusterRef[iPivot], 
					    particleEnergy[iPivot] ); 
      
	(*pfCandidates_)[tmpi].setEcalEnergy( ecalEnergy[iPivot] );
	(*pfCandidates_)[tmpi].setHcalEnergy( hcalEnergy[iPivot] );
	(*pfCandidates_)[tmpi].setPs1Energy( -1 );
	(*pfCandidates_)[tmpi].setPs2Energy( -1 );
	(*pfCandidates_)[tmpi].set_mva_nothing_gamma( maxMva );
	//       (*pfCandidates_)[tmpi].addElement(&elements[iPivotal]);
	(*pfCandidates_)[tmpi].addElementInBlock(blockref, iPivotal[iPivot]);
	// add PS elements to be done

      }

    } // excess of energy
  

    // will now share the hcal energy between the various charged hadron 
    // candidates, taking into account the potential neutral hadrons
    
    double totalHcalEnergyCalibrated = calibHcal;
    /*
    if(totalEcal>0) {
      // removing ecal energy from abc calibration
      totalHcalEnergyCalibrated  -= calibEcal;
      // totalHcalEnergyCalibrated -= calibration_->paramECALplusHCAL_slopeECAL() * totalEcal;
    }
    */
    // else caloEnergy = hcal only calibrated energy -> ok.

    // remove the energy of the potential neutral hadron
    totalHcalEnergyCalibrated -= mergedNeutralHadronEnergy;

    // share between the charged hadrons

    //COLIN can compute this before
    // not exactly equal to sum p, this is sum E
    double chargedHadronsTotalEnergy = 0;
    for( unsigned ich=0; ich<chargedHadronsIndices.size(); ++ich ) {
      unsigned index = chargedHadronsIndices[ich];
      reco::PFCandidate& chargedHadron = (*pfCandidates_)[index];
      chargedHadronsTotalEnergy += chargedHadron.energy();
    }

    for( unsigned ich=0; ich<chargedHadronsIndices.size(); ++ich ) {
      unsigned index = chargedHadronsIndices[ich];
      reco::PFCandidate& chargedHadron = (*pfCandidates_)[index];
      float fraction = chargedHadron.energy()/chargedHadronsTotalEnergy;
      chargedHadron.setHcalEnergy( fraction * totalHcalEnergyCalibrated );          
    }

  } // end loop on hcal element iHcal= hcalIs[i] 


  // Processing the remaining HCAL clusters
  if(debug_) { 
    cout<<endl;
    if(debug_) 
      cout<<endl
          <<"---- loop remaining hcal ------- "<<endl;
  }
  

  //COLINFEB16 
  // now dealing with the HCAL elements that are not linked to any track 
  for(unsigned ihcluster=0; ihcluster<hcalIs.size(); ihcluster++) {
    unsigned iHcal = hcalIs[ihcluster];
    
    if(debug_)
      cout<<endl<<elements[iHcal]<<" ";
    
    
    if( !active[iHcal] ) {
      if(debug_) 
        cout<<"not active"<<endl;         
      continue;
    }
    
    // Find the ECAL elements linked to it
    std::multimap<double, unsigned> ecalElems;
    block.associatedElements( iHcal,  linkData,
                              ecalElems ,
                              reco::PFBlockElement::ECAL,
			      reco::PFBlock::LINKTEST_ALL );

    // Loop on these ECAL elements
    float totalEcal = 0.;
    float ecalMax = 0.;
    unsigned jEcal = 0;
    reco::PFClusterRef eClusterRef;
    for(IE ie = ecalElems.begin(); ie != ecalElems.end(); ++ie ) {
      
      unsigned iEcal = ie->second;
      double dist = ie->first;
      double chi2 = block.chi2(iHcal,iEcal,linkData,reco::PFBlock::LINKTEST_ALL);
      PFBlockElement::Type type = elements[iEcal].type();
      assert( type == PFBlockElement::ECAL );
      
      // Check if already used
      if( !active[iEcal] ) continue;

      // std::cout << "chi2, dist : " << chi2 << ", " << dist << std::endl;
      // Check the chi2
      //COLINFEB16 
      // these are parameters and should not be hardcoded. 
      // PJ FEB17 : 
      // 1. Link by chi**2 are on the verge of disappearing
      // 2. The cut on the chi2 is anyway set to 0 in the particleFlow.opt.
      // But I agree that no parameter must be hardcoded
      if ( chi2 > 50. ) continue;

      // Check the distance (one HCALPlusECAL tower, roughly)
      if ( dist > 0.15 ) continue;

      //COLINFEB16 
      // what could be done is to
      // - link by rechit.
      // - take in the neutral hadron all the ECAL clusters 
      // which are within the same CaloTower, according to the distance, 
      // except the ones which are closer to another HCAL cluster.
      // - all the other ECAL linked to this HCAL are photons. 
      // 
      // about the closest HCAL cluster. 
      // it could maybe be easier to loop on the ECAL clusters first
      // to cut the links to all HCAL clusters except the closest, as is 
      // done in the first track loop. But maybe not!
      // or add an helper function to the PFAlgo class to ask 
      // if a given element is the closest of a given type to another one?
    
      // Check if not closer from another free HCAL 
      std::multimap<double, unsigned> hcalElems;
      block.associatedElements( iEcal,  linkData,
				hcalElems ,
				reco::PFBlockElement::HCAL,
				reco::PFBlock::LINKTEST_ALL );

      bool isClosest = true;
      for(IE ih = hcalElems.begin(); ih != hcalElems.end(); ++ih ) {

	unsigned jHcal = ih->second;
	double distH = ih->first;
	
	if ( !active[jHcal] ) continue;

	if ( distH < dist ) { 
	  isClosest = false;
	  break;
	}

      }

      if (!isClosest) continue;


      if(debug_) {
        cout<<"\telement "<<elements[iEcal]<<" linked with chi2/dist "<<chi2<<" / " << dist<<endl;
	cout<<"Added to HCAL cluster to form a neutral hadron"<<endl;
      }
      
      reco::PFClusterRef eclusterRef = elements[iEcal].clusterRef();
      assert( !eclusterRef.isNull() );

      // Check the presence of ps clusters in the vicinity
      vector<double> ps1Ene(1,static_cast<double>(0.));
      associatePSClusters(iEcal, reco::PFBlockElement::PS1, block, elements, linkData, active, ps1Ene);
      vector<double> ps2Ene(1,static_cast<double>(0.));
      associatePSClusters(iEcal, reco::PFBlockElement::PS2, block, elements, linkData, active, ps2Ene);
      double ecalEnergy = calibration_->energyEm(*eclusterRef,ps1Ene,ps2Ene);
      
      //std::cout << "EcalEnergy, ps1, ps2 = " << ecalEnergy 
      //          << ", " << ps1Ene[0] << ", " << ps2Ene[0] << std::endl;
      totalEcal += ecalEnergy;
      if ( ecalEnergy > ecalMax ) { 
	ecalMax = ecalEnergy;
	eClusterRef = eclusterRef;
	jEcal = iEcal;
      }
	
      active[iEcal] = false;

    }// End loop ECAL

    
    PFClusterRef hclusterRef 
      = elements[iHcal].clusterRef();
    assert( !hclusterRef.isNull() );
    
    // HCAL energy
    double totalHcal =  hclusterRef->energy();
    // std::cout << "Hcal Energy,eta : " << totalHcal 
    //          << ", " << hclusterRef->positionREP().Eta()
    //          << std::endl;
    // Calibration
    double caloEnergy = totalHcal;
    double slopeEcal = 1.0;
    double calibEcal = totalEcal > 0. ? totalEcal : 0.;
    double calibHcal = std::max(0.,totalHcal);
    if ( newCalib_ == 1 ) {
      clusterCalibration_->
	getCalibratedEnergyEmbedAInHcal(calibEcal, calibHcal,
					hclusterRef->positionREP().Eta(),
					hclusterRef->positionREP().Phi());
      if ( calibEcal == 0. ) calibEcal = totalEcal;
      if ( calibHcal == 0. ) calibHcal = totalHcal;
      caloEnergy = calibEcal+calibHcal;
      if ( totalEcal > 0. ) slopeEcal = calibEcal/totalEcal;
    } else if ( newCalib_ == 0 ) { 
      if( totalEcal>0) { 
	caloEnergy = calibration_->energyEmHad( totalEcal, totalHcal );
	slopeEcal = calibration_->paramECALplusHCAL_slopeECAL();
	calibEcal = totalEcal * slopeEcal;
	calibHcal = caloEnergy - calibEcal;
      } else if  ((  hclusterRef->layer() != PFLayer::HF_HAD ) && (hclusterRef->layer() != PFLayer::HF_EM)){
	caloEnergy = calibration_->energyHad( totalHcal );
	calibEcal = totalEcal;
	calibHcal = caloEnergy;
      } else { 
	caloEnergy = totalHcal/0.7;
	calibEcal = totalEcal;
	calibHcal = caloEnergy;
      }
    } else {
      if  (  hclusterRef->layer() == PFLayer::HF_HAD  ||
	     hclusterRef->layer() == PFLayer::HF_EM ) { 
	caloEnergy = totalHcal/0.7;
	calibEcal = totalEcal;
	calibHcal = caloEnergy;
      } else { 
	calibration_->energyEmHad(-1.,calibEcal,calibHcal,
				  hclusterRef->positionREP().Eta(),
				  hclusterRef->positionREP().Phi());
	caloEnergy = calibEcal+calibHcal;
      }
    }

    // std::cout << "CalibEcal,HCal = " << calibEcal << ", " << calibHcal << std::endl;
    // std::cout << "-------------------------------------------------------------------" << std::endl;
    // ECAL energy : calibration

    // double particleEnergy = caloEnergy;
    // double particleEnergy = totalEcal + calibHcal;
    // particleEnergy /= (1.-0.724/sqrt(particleEnergy)-0.0226/particleEnergy);

    unsigned tmpi = reconstructCluster( *hclusterRef, 
                                        calibEcal+calibHcal ); 

    
    (*pfCandidates_)[tmpi].setEcalEnergy( 0. );
    (*pfCandidates_)[tmpi].setHcalEnergy( calibHcal +calibEcal);
    (*pfCandidates_)[tmpi].setPs1Energy( -1 );
    (*pfCandidates_)[tmpi].setPs2Energy( -1 );
    (*pfCandidates_)[tmpi].addElementInBlock( blockref, iHcal );
      
    /*
    if ( totalEcal > 0. ) { 
      unsigned tmpi = reconstructCluster( *eClusterRef, 
					  totalEcal ); 
      
      (*pfCandidates_)[tmpi].setEcalEnergy( totalEcal );
      (*pfCandidates_)[tmpi].setHcalEnergy( 0. );
      (*pfCandidates_)[tmpi].setPs1Energy( -1 );
      (*pfCandidates_)[tmpi].setPs2Energy( -1 );
      (*pfCandidates_)[tmpi].addElementInBlock( blockref, jEcal );
    }
    */
      
  }//loop hcal elements




  if(debug_) { 
    cout<<endl;
    if(debug_) cout<<endl<<"---- loop ecal------- "<<endl;
  }
  
  // for each ecal element iEcal = ecalIs[i] in turn:

  for(unsigned i=0; i<ecalIs.size(); i++) {
    unsigned iEcal = ecalIs[i];
    
    if(debug_) 
      cout<<endl<<elements[iEcal]<<" ";
    
    if( ! active[iEcal] ) {
      if(debug_) 
        cout<<"not active"<<endl;         
      continue;
    }
    
    PFBlockElement::Type type = elements[ iEcal ].type();
    assert( type == PFBlockElement::ECAL );

    PFClusterRef clusterref = elements[iEcal].clusterRef();
    assert(!clusterref.isNull() );

    active[iEcal] = false;


    // Check the presence of ps clusters in the vicinity
    vector<double> ps1Ene(1,static_cast<double>(0.));
    associatePSClusters(iEcal, reco::PFBlockElement::PS1, block, elements, linkData, active, ps1Ene);
    vector<double> ps2Ene(1,static_cast<double>(0.));
    associatePSClusters(iEcal, reco::PFBlockElement::PS2, block, elements, linkData, active, ps2Ene);
    float ecalEnergy = calibration_->energyEm(*clusterref,ps1Ene,ps2Ene);
    // float ecalEnergy = calibration_->energyEm( clusterref->energy() );
    double particleEnergy = ecalEnergy;
    
    unsigned tmpi = reconstructCluster( *clusterref, 
                                        particleEnergy );
 
    (*pfCandidates_)[tmpi].setEcalEnergy( ecalEnergy );
    (*pfCandidates_)[tmpi].setHcalEnergy( 0 );
    (*pfCandidates_)[tmpi].setPs1Energy( -1 );
    (*pfCandidates_)[tmpi].setPs2Energy( -1 );
    (*pfCandidates_)[tmpi].addElementInBlock( blockref, iEcal );
    

  }  // end loop on ecal elements iEcal = ecalIs[i]
 
}  // end processBlock



unsigned PFAlgo::reconstructTrack( const reco::PFBlockElement& elt ) {

  const reco::PFBlockElementTrack* eltTrack 
    = dynamic_cast<const reco::PFBlockElementTrack*>(&elt);

  reco::TrackRef trackRef = eltTrack->trackRef();
  const reco::Track& track = *trackRef;

  reco::MuonRef muonRef = eltTrack->muonRef();

  int charge = track.charge()>0 ? 1 : -1;

  // Assign the pion mass to all charged particles
  double px = track.px();
  double py = track.py();
  double pz = track.pz();
  double energy = sqrt(track.p()*track.p() + 0.13957*0.13957);

  // Except if it is a muon, of course !
  bool thisIsAMuon = PFMuonAlgo::isMuon(elt);
  if ( thisIsAMuon ) { 
    reco::TrackRef combinedMu = muonRef->combinedMuon();
    px = combinedMu->px();
    py = combinedMu->py();
    pz = combinedMu->pz();
    energy = sqrt(combinedMu->p()*combinedMu->p() + 0.1057*0.1057);
  } 

  // Create a PF Candidate
  math::XYZTLorentzVector momentum(px,py,pz,energy);
  reco::PFCandidate::ParticleType particleType 
    = reco::PFCandidate::h;

  // Add it to the stack
  pfCandidates_->push_back( PFCandidate( charge, 
                                         momentum,
                                         particleType ) );
  
  pfCandidates_->back().setVertex( track.vertex() );
  pfCandidates_->back().setTrackRef( trackRef );
  pfCandidates_->back().setPositionAtECALEntrance( eltTrack->positionAtECALEntrance());


  // setting the muon ref if there is
  if (muonRef.isNonnull()) {
    pfCandidates_->back().setMuonRef( muonRef );

    // setting the muon particle type if it is a global muon
    if ( thisIsAMuon ) {
      particleType = reco::PFCandidate::mu;
      pfCandidates_->back().setParticleType( particleType );
      if (debug_) cout << "PFAlgo: particle type set to muon" << endl; 
    }
  }

  // nuclear
  if( particleType != reco::PFCandidate::mu )
    if( eltTrack->trackType(reco::PFBlockElement::T_FROM_NUCL)) {
      pfCandidates_->back().setFlag( reco::PFCandidate::T_FROM_NUCLINT, true);
      pfCandidates_->back().setNuclearRef( eltTrack->nuclearRef() );
    }
    else if( eltTrack->trackType(reco::PFBlockElement::T_TO_NUCL)) {
      pfCandidates_->back().setFlag( reco::PFCandidate::T_TO_NUCLINT, true);
      pfCandidates_->back().setNuclearRef( eltTrack->nuclearRef() );
    }

  // conversion...
  
  if(debug_) 
    cout<<"** candidate: "<<pfCandidates_->back()<<endl; 
  
  // returns index to the newly created PFCandidate
  return pfCandidates_->size()-1;
}





unsigned PFAlgo::reconstructCluster(const reco::PFCluster& cluster,
                                    double particleEnergy) {
  
  reco::PFCandidate::ParticleType particleType = reco::PFCandidate::X;

  // need to convert the math::XYZPoint data member of the PFCluster class=
  // to a displacement vector: 
  ROOT::Math::DisplacementVector3D< ROOT::Math::Cartesian3D<double> , ROOT::Math::DefaultCoordinateSystemTag > clusterPos( cluster.position().X(), cluster.position().Y(),cluster.position().Z() );
  

  clusterPos = clusterPos.Unit();
  clusterPos *= particleEnergy;

  // clusterPos is now a vector along the cluster direction, 
  // with a magnitude equal to the cluster energy.
  
  double mass = 0;
  switch( cluster.layer() ) {
  case PFLayer::ECAL_BARREL:
  case PFLayer::ECAL_ENDCAP:
    particleType = PFCandidate::gamma;
    break;
  case PFLayer::HCAL_BARREL1:
  case PFLayer::HCAL_ENDCAP:
    particleType = PFCandidate::h0;
    break;
  case PFLayer::HF_HAD:
    particleType = PFCandidate::h_HF;
    break;
  case PFLayer::HF_EM:
    particleType = PFCandidate::egamma_HF;
    break;
  default:
    assert(0);
  }


  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > 
    momentum( clusterPos.X(), clusterPos.Y(), clusterPos.Z(), mass); 

  int charge = 0;

  // mathcore is a piece of #$%
  math::XYZTLorentzVector  tmp;

  // implicit constructor not allowed
  tmp = momentum;

  pfCandidates_->push_back( PFCandidate( charge, 
                                         tmp, 
                                         particleType ) );

  if(debug_) 
    cout<<"** candidate: "<<pfCandidates_->back()<<endl; 

  // returns index to the newly created PFCandidate
  return pfCandidates_->size()-1;

}




double PFAlgo::neutralHadronEnergyResolution( double clusterEnergyHCAL, double eta ) 
  const {


  double resol = 0.;
  if ( newCalib_ == 1 ) 
    resol =   1.40/sqrt(clusterEnergyHCAL) +5.00/clusterEnergyHCAL;
  else if ( newCalib_ == 0 ) 
    resol =   1.50/sqrt(clusterEnergyHCAL) +3.00/clusterEnergyHCAL;
  else 
    resol =   fabs(eta) < 1.48 ? 
      //min(0.25,sqrt (1.02*1.02/clusterEnergyHCAL + 0.065*0.065)):
      //min(0.30,sqrt (1.35*1.35/clusterEnergyHCAL + 0.018*0.018));
      sqrt (1.02*1.02/clusterEnergyHCAL + 0.065*0.065):
      sqrt (1.35*1.35/clusterEnergyHCAL + 0.018*0.018);

  return resol;
}


ostream& operator<<(ostream& out, const PFAlgo& algo) {
  if(!out ) return out;

  out<<"====== Particle Flow Algorithm ======= ";
  out<<endl;
  out<<"nSigmaECAL_     "<<algo.nSigmaECAL_<<endl;
  out<<"nSigmaHCAL_     "<<algo.nSigmaHCAL_<<endl;
  out<<"mvaCut_         "<<algo.mvaCut_<<endl;
  out<<"PSCut_          "<<algo.PSCut_<<endl;
  out<<endl;
  out<<(*(algo.calibration_))<<endl;
  out<<endl;
  out<<"reconstructed particles: "<<endl;
   
  const std::auto_ptr< reco::PFCandidateCollection >& 
    candidates = algo.pfCandidates(); 

  if(!candidates.get() ) {
    out<<"candidates already transfered"<<endl;
    return out;
  }
  for(PFCandidateConstIterator ic=algo.pfCandidates_->begin(); 
      ic != algo.pfCandidates_->end(); ic++) {
    out<<(*ic)<<endl;
  }

  return out;
}


bool PFAlgo::isSatelliteCluster( const reco::PFRecTrack& track, 
                                 const PFCluster& cluster ){
  //This is to check whether a HCAL cluster could be 
  //a satellite cluster of a given track (charged hadron)
  //Definitions:
  //Track-Hcal: An Hcal cluster can be associated to a 
  //            track if it is found in a region around (dr<0.17 ~ 3x3 matrix)
  //            the position of the track extrapolated to the Hcal.

  bool link = false;

  if( cluster.layer() == PFLayer::HCAL_BARREL1 ||
      cluster.layer() == PFLayer::HCAL_ENDCAP ) { //Hcal case
    
    const reco::PFTrajectoryPoint& atHCAL 
      = track.extrapolatedPoint( reco::PFTrajectoryPoint::HCALEntrance );
              
    if( atHCAL.isValid() ){ //valid extrapolation?
      double tracketa = atHCAL.position().Eta();
      double trackphi = atHCAL.position().Phi();
      double hcaleta  = cluster.positionREP().Eta();
      double hcalphi  = cluster.positionREP().Phi();
                  
      //distance track-cluster
      double deta = hcaleta - tracketa;
      double dphi = acos(cos(hcalphi - trackphi));
      double dr   = sqrt(deta*deta + dphi*dphi);

      if( debug_ ){
        cout<<"\t\t\tSatellite Test " 
            <<tracketa<<" "<<trackphi<<" "
            <<hcaleta<<" "<<hcalphi<<" dr="<<dr 
            <<endl;
      }
      
      //looking if cluster is in the  
      //region around the track. 
      //Alex: Will have to adjust this cut?
      if( dr < 0.17 ) link = true;
    }//extrapolation
    
  }//HCAL
  
  return link; 
}

reco::PFBlockRef 
PFAlgo::createBlockRef( const reco::PFBlockCollection& blocks, 
			unsigned bi ) {

  if( blockHandle_.isValid() ) {
    return reco::PFBlockRef(  blockHandle_, bi );
  } 
  else {
    return reco::PFBlockRef(  &blocks, bi );
  }
}

void 
PFAlgo::associatePSClusters(unsigned iEcal,
			    reco::PFBlockElement::Type psElementType,
			    const reco::PFBlock& block,
			    const edm::OwnVector< reco::PFBlockElement >& elements, 
			    const reco::PFBlock::LinkData& linkData,
			    std::vector<bool>& active, 
			    std::vector<double>& psEne) {

  // Find all PS clusters with type psElement associated to ECAL cluster iEcal, 
  // within all PFBlockElement "elements" of a given PFBlock "block"
  // psElement can be reco::PFBlockElement::PS1 or reco::PFBlockElement::PS2
  // Returns a vector of PS cluster energies, and updates the "active" vector.

  // Find all PS clusters linked to the iEcal cluster
  std::multimap<double, unsigned> sortedPS;
  typedef std::multimap<double, unsigned>::iterator IE;
  block.associatedElements( iEcal,  linkData,
			    sortedPS, psElementType,
			    reco::PFBlock::LINKTEST_ALL );

  // Loop over these PS clusters
  double totalPS = 0.;
  for ( IE ips=sortedPS.begin(); ips!=sortedPS.end(); ++ips ) {

    // CLuster index and distance to iEcal
    unsigned iPS = ips->second;
    // double distPS = ips->first;

    // Ignore clusters already in use
    if (!active[iPS]) continue;

    // Check that this cluster is not closer to another ECAL cluster
    std::multimap<double, unsigned> sortedECAL;
    block.associatedElements( iPS,  linkData,
			      sortedECAL,
			      reco::PFBlockElement::ECAL,
			      reco::PFBlock::LINKTEST_ALL );
    unsigned jEcal = sortedECAL.begin()->second;
    if ( jEcal != iEcal ) continue; 
    
    // Update PS energy
    PFBlockElement::Type pstype = elements[ iPS ].type();
    assert( pstype == psElementType );
    PFClusterRef psclusterref = elements[iPS].clusterRef();
    assert(!psclusterref.isNull() );
    totalPS += psclusterref->energy(); 
    psEne[0] += psclusterref->energy();
    active[iPS] = false;
  }
	    

}