#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerTrackEcalAOD.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "TMath.h"


KDTreeLinkerTrackEcalAOD::KDTreeLinkerTrackEcalAOD()
  : KDTreeLinkerBase()
{}

KDTreeLinkerTrackEcalAOD::~KDTreeLinkerTrackEcalAOD()
{
  clear();
}

void
KDTreeLinkerTrackEcalAOD::insertTargetElt(reco::PFBlockElement	*track)
{
  targetSet_.insert(track);
}


void
KDTreeLinkerTrackEcalAOD::insertFieldClusterElt(reco::PFBlockElement	*ecalCluster)
{
  reco::PFClusterRef clusterref = ecalCluster->clusterRef();

  // This test is more or less done in PFBlockAlgo.h. In others cases, it should be switch on.
  //   if (!((clusterref->layer() == PFLayer::ECAL_ENDCAP) ||
  // 	(clusterref->layer() == PFLayer::ECAL_BARREL)))
  //     return;

  const std::vector< std::pair<DetId, float> >& hitsAndFracs = clusterref->hitsAndFractions();

  // We create a list of ecalCluster
  fieldClusterSet_.insert(ecalCluster);

  for(size_t iId = 0; iId < hitsAndFracs.size(); ++iId) {
    
	 const DetId* di = &(hitsAndFracs[iId].first);
    double fract = hitsAndFracs[iId].second;

//    if ((di.isNull()) || (fract < 1E-4))
    if (fract < 1E-4)  // add check for ID validity (AA)
	 continue;
      
    const DetId& detid = *di;
      
    // We save the links of detids to EcalClusters
    detid2ClusterLinks_[&detid].insert(ecalCluster);
    
    // We create a liste of detids
    detidsSet_.insert(&detid);
  } 

}

void 
KDTreeLinkerTrackEcalAOD::buildTree()
{
  // List of pseudo-rechits that will be used to create the KDTree
  std::vector<KDTreeNodeInfoAOD> eltList;

  // Filling of this list
  for(DetIdSet::const_iterator it = detidsSet_.begin(); 
      it != detidsSet_.end(); it++) {
    
	 // here we get the position of the cell, and convert to the formats in PF for now (AA)
    const CaloCellGeometry* cellGeometry = caloGeometry_->getGeometry(**it);
	 // IMPORTANT: see in the PFRecHitProducer how exactly is the position defined  <---------------- !!!!!
	 // this is temporary and uses the face of the cel (AA)
	 const GlobalPoint gpPosition = cellGeometry->getPosition();
	 const reco::PFRecHit::REPPoint posrep(gpPosition.perp(), gpPosition.eta(), gpPosition.phi());


    KDTreeNodeInfoAOD di1 (*it, posrep.Eta(), posrep.Phi());
    eltList.push_back(di1);
    
    // Here we solve the problem of phi circular set by duplicating some rechits
    // too close to -Pi (or to Pi) and adding (substracting) to them 2 * Pi.
    if (di1.dim2 > (M_PI - getPhiOffset())) {
      double phi = di1.dim2 - 2 * M_PI;
      KDTreeNodeInfoAOD di2(*it, posrep.Eta(), phi); 
      eltList.push_back(di2);
    }

    if (di1.dim2 < (M_PI * -1.0 + getPhiOffset())) {
      double phi = di1.dim2 + 2 * M_PI;
      KDTreeNodeInfoAOD di3(*it, posrep.Eta(), phi); 
      eltList.push_back(di3);
    }
  }

  // Here we define the upper/lower bounds of the 2D space (eta/phi).
  double phimin = -1.0 * M_PI - getPhiOffset();
  double phimax = M_PI + getPhiOffset();

  // etamin-etamax, phimin-phimax
  KDTreeBox region(-3.0, 3.0, phimin, phimax);

  // We may now build the KDTree
  tree_.build(eltList, region);
}

void
KDTreeLinkerTrackEcalAOD::searchLinks()
{
  // Most of the code has been taken from LinkByRecHit.cc

  // We iterate over the tracks.
  for(BlockEltSet::iterator it = targetSet_.begin(); 
      it != targetSet_.end(); it++) {
	
    reco::PFRecTrackRef trackref = (*it)->trackRefPF();

    // We set the multilinks flag of the track to true. It will allow us to 
    // use in an optimized way our algo results in the recursive linking algo.
    (*it)->setIsValidMultilinks(true);

    const reco::PFTrajectoryPoint& atECAL = 
      trackref->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax);

    // The track didn't reach ecal
    if( ! atECAL.isValid() ) continue;
    
    const reco::PFTrajectoryPoint& atVertex = 
      trackref->extrapolatedPoint( reco::PFTrajectoryPoint::ClosestApproach );
    
    double trackPt = sqrt(atVertex.momentum().Vect().Perp2());
    double tracketa = atECAL.positionREP().Eta();
    double trackphi = atECAL.positionREP().Phi();
    double trackx = atECAL.position().X();
    double tracky = atECAL.position().Y();
    double trackz = atECAL.position().Z();
    
    // Estimate the maximal envelope in phi/eta that will be used to find rechit candidates.
    // Same envelope for cap et barrel rechits.
    double range = getCristalPhiEtaMaxSize() * (2.0 + 1.0 / std::min(1., trackPt / 2.)); 

    // We search for all candidate detIds, ie all detIds contained in the maximal size envelope.
    std::vector<KDTreeNodeInfoAOD> detIds;
    KDTreeBox trackBox(tracketa-range, tracketa+range, trackphi-range, trackphi+range);
    tree_.search(trackBox, detIds);
    
    // Here we check all detid candidates using the non-approximated method.
    for(std::vector<KDTreeNodeInfoAOD>::const_iterator dId = detIds.begin(); 
	 dId != detIds.end(); ++dId) {
           
      // get the position and corners - check for depth corrections (AA) <---------------------------------- !!!!
		const CaloCellGeometry* cellGeometry = caloGeometry_->getGeometry(*(dId->ptr));
		const GlobalPoint gpPosition = cellGeometry->getPosition();
		CaloCellGeometry::CornersVec cellGeoCorners = cellGeometry->getCorners();
		// convert to types used in the PF framework
      const math::XYZPoint posxyz(gpPosition.x(), gpPosition.y(), gpPosition.z()); 
      const reco::PFRecHit::REPPoint posrep(gpPosition.perp(), gpPosition.eta(), gpPosition.phi());
		std::vector< math::XYZPoint > cornersxyz;
	   std::vector<reco::PFRecHit::REPPoint> corners;
	   for (uint i=0; i<4; ++i) {
			cornersxyz.push_back(math::XYZPoint(cellGeoCorners[3-i].x(),
															cellGeoCorners[3-i].y(),
															cellGeoCorners[3-i].z()
															));
			corners.push_back(reco::PFRecHit::REPPoint(cellGeoCorners[3-i].perp(),
						   										  cellGeoCorners[3-i].eta(),
																	  cellGeoCorners[3-i].phi()
																	  ));
		} // end of corner re-ordering according to the convention in PFRecHit

      if(corners.size() != 4) continue;
      
      double cellSizeEta = fabs(corners[0].Eta() - corners[2].Eta());
      double cellSizePhi = fabs(corners[0].Phi() - corners[2].Phi());
      if ( cellSizePhi > M_PI ) cellSizePhi = 2.*M_PI - cellSizePhi;
      
      double deta = fabs(posrep.Eta() - tracketa);
      double dphi = fabs(posrep.Phi() - trackphi);
      if ( dphi > M_PI ) dphi = 2.*M_PI - dphi;
      
      // Find all clusters associated to a given detid
      DetId2BlockEltMap::iterator ret = detid2ClusterLinks_.find(dId->ptr);
      
      for(BlockEltSet::const_iterator clusterIt = ret->second.begin(); 
	  clusterIt != ret->second.end(); clusterIt++) {
	
	reco::PFClusterRef clusterref = (*clusterIt)->clusterRef();
	double clusterz = clusterref->position().Z();
	int fracsNbr = clusterref->hitsAndFractions().size();

	if (clusterref->layer() == PFLayer::ECAL_BARREL){ // BARREL
	  // Check if the track is in the barrel
	  if (fabs(trackz) > 300.) continue;

	  double _cellSizeEta = cellSizeEta * (2.00 + 1.0 / (fracsNbr * std::min(1.,trackPt/2.)));
	  double _cellSizePhi = cellSizePhi * (2.00 + 1.0 / (fracsNbr * std::min(1.,trackPt/2.)));
	  
	  // Check if the track and the cluster are linked
	  if(deta < (_cellSizeEta / 2.) && dphi < (_cellSizePhi / 2.))
	    target2ClusterLinks_[*it].insert(*clusterIt);

	  
	} else { // ENDCAP

	  // Check if the track is in the cap
	  if (fabs(trackz) < 300.) continue;
	  if (trackz * clusterz < 0.) continue;
	  
	  double x[5];
	  double y[5];
	  for ( unsigned jc=0; jc<4; ++jc ) {
	    math::XYZPoint cornerposxyz = cornersxyz[jc];
	    x[jc] = cornerposxyz.X() + (cornerposxyz.X()-posxyz.X())
	      * (1.00+0.50/fracsNbr /std::min(1.,trackPt/2.));
	    y[jc] = cornerposxyz.Y() + (cornerposxyz.Y()-posxyz.Y())
	      * (1.00+0.50/fracsNbr /std::min(1.,trackPt/2.));
	  }
	  
	  x[4] = x[0];
	  y[4] = y[0];
	  
	  bool isinside = TMath::IsInside(trackx,
					  tracky,
					  5,x,y);
	  
	  // Check if the track and the cluster are linked
	  if( isinside )
	    target2ClusterLinks_[*it].insert(*clusterIt);
	}
      }
    }
  }
}

void
KDTreeLinkerTrackEcalAOD::updatePFBlockEltWithLinks()
{
  //TODO YG : Check if cluster positionREP() is valid ?

  // Here we save in each track the list of phi/eta values of linked clusters.
  for (BlockElt2BlockEltMap::iterator it = target2ClusterLinks_.begin();
       it != target2ClusterLinks_.end(); ++it) {
    reco::PFMultiLinksTC multitracks(true);

    for (BlockEltSet::iterator jt = it->second.begin();
	 jt != it->second.end(); ++jt) {

      double clusterPhi = (*jt)->clusterRef()->positionREP().Phi();
      double clusterEta = (*jt)->clusterRef()->positionREP().Eta();

      multitracks.linkedClusters.push_back(std::make_pair(clusterPhi, clusterEta));
    }

    it->first->setMultilinks(multitracks);
  }
}

void
KDTreeLinkerTrackEcalAOD::clear()
{
  targetSet_.clear();
  fieldClusterSet_.clear();

  detidsSet_.clear();

  detid2ClusterLinks_.clear();
  target2ClusterLinks_.clear();

  tree_.clear();
}
