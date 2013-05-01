#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerTrackHcalAOD.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "TMath.h"

KDTreeLinkerTrackHcalAOD::KDTreeLinkerTrackHcalAOD()
  : KDTreeLinkerBase()
{
  setCristalPhiEtaMaxSize(0.2);
  setPhiOffset(0.32);
}

KDTreeLinkerTrackHcalAOD::~KDTreeLinkerTrackHcalAOD()
{
  clear();
}

void
KDTreeLinkerTrackHcalAOD::insertTargetElt(reco::PFBlockElement	*track)
{
  targetSet_.insert(track);
}


void
KDTreeLinkerTrackHcalAOD::insertFieldClusterElt(reco::PFBlockElement	*hcalCluster)
{
  reco::PFClusterRef clusterref = hcalCluster->clusterRef();

  // This test is more or less done in PFBlockAlgo.h. In others cases, it should be switch on.
  //   if (!((clusterref->layer() == PFLayer::HCAL_ENDCAP) ||
  // 	(clusterref->layer() == PFLayer::HCAL_BARREL1)))
  //     return;

  const std::vector< std::pair<DetId, float> >& hitsAndFracs = clusterref->hitsAndFractions();

  // We create a list of hcalCluster
  fieldClusterSet_.insert(hcalCluster);

  for(size_t iId = 0; iId < hitsAndFracs.size(); ++iId) {
    
	 const DetId* di = &(hitsAndFracs[iId].first);
    double fract = hitsAndFracs[iId].second;

//    if ((di.isNull()) || (fract < 1E-4))
    if (fract < 1E-4)  // add check for ID validity (AA)
	 continue;
      
    const DetId& detid = *di;
      
    // We save the links of detids to HcalClusters
    detid2ClusterLinks_[&detid].insert(hcalCluster);
    
    // We create a liste of detids
    detidsSet_.insert(&detid);
  }

}

void 
KDTreeLinkerTrackHcalAOD::buildTree()
{
  // List of pseudo-rechits that will be used to create the KDTree
  std::vector<KDTreeNodeInfoAOD> eltList;

  // Filling of this list
  for(DetIdSet::const_iterator it = detidsSet_.begin(); 
      it != detidsSet_.end(); it++) {
    
	 // here we get the position of the cell, and convert to the formats in PF for now (AA)
    const CaloCellGeometry* cellGeometry = caloGeometry_->getGeometry(**it);
	 // IMPORTANT: see in the PFRecHitProducer how exactly is the position defined  <---------------- !!!!!
	 // this is temporary and uses the facoe of the cel (AA)
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
KDTreeLinkerTrackHcalAOD::searchLinks()
{
  // Most of the code has been taken from LinkByRecHit.cc

  // We iterate over the tracks.
  for(BlockEltSet::iterator it = targetSet_.begin(); 
      it != targetSet_.end(); it++) {
	
    reco::PFRecTrackRef trackref = (*it)->trackRefPF();

    const reco::PFTrajectoryPoint& atHCAL = 
      trackref->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance);
    const reco::PFTrajectoryPoint& atHCALExit = 
      trackref->extrapolatedPoint(reco::PFTrajectoryPoint::HCALExit);
    
    // The track didn't reach hcal
    if( ! atHCAL.isValid()) continue;
    
    double dHEta = atHCALExit.positionREP().Eta() - atHCAL.positionREP().Eta();
    double dHPhi = atHCALExit.positionREP().Phi() - atHCAL.positionREP().Phi(); 
    if ( dHPhi > M_PI ) dHPhi = dHPhi - 2. * M_PI;
    else if ( dHPhi < -M_PI ) dHPhi = dHPhi + 2. * M_PI; 
    
    double tracketa = atHCAL.positionREP().Eta() + 0.1 * dHEta;
    double trackphi = atHCAL.positionREP().Phi() + 0.1 * dHPhi;
    
    if (trackphi > M_PI) trackphi -= 2 * M_PI;
    else if (trackphi < -M_PI) trackphi += 2 * M_PI;

    // Estimate the maximal envelope in phi/eta that will be used to find rechit candidates.
    // Same envelope for cap et barrel rechits.
    double inflation = 1.;
    double rangeEta = (getCristalPhiEtaMaxSize() * (1.5 + 0.5) + 0.2 * fabs(dHEta)) * inflation; 
    double rangePhi = (getCristalPhiEtaMaxSize() * (1.5 + 0.5) + 0.2 * fabs(dHPhi)) * inflation; 

    // We search for all candidate detIds, ie all detIds contained in the maximal size envelope.
    std::vector<KDTreeNodeInfoAOD> detIds;
    KDTreeBox trackBox(tracketa - rangeEta, tracketa + rangeEta, 
		       trackphi - rangePhi, trackphi + rangePhi);
    tree_.search(trackBox, detIds);
    
    // Here we check all detid candidates using the non-approximated method.
    for(std::vector<KDTreeNodeInfoAOD>::const_iterator dId = detIds.begin(); 
	dId != detIds.end(); ++dId) {
      
      // get the position corners - check for depth corrections (AA) <---------------------------------- !!!!
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
	int fracsNbr = clusterref->hitsAndFractions().size();
	
	double _cellSizeEta = cellSizeEta * (1.5 + 0.5 / fracsNbr) + 0.2 * fabs(dHEta);
	double _cellSizePhi = cellSizePhi * (1.5 + 0.5 / fracsNbr) + 0.2 * fabs(dHPhi);
	
	// Check if the track and the cluster are linked
	if(deta < (_cellSizeEta / 2.) && dphi < (_cellSizePhi / 2.))
	  cluster2TargetLinks_[*clusterIt].insert(*it);
      }
    }
  }
}

void
KDTreeLinkerTrackHcalAOD::updatePFBlockEltWithLinks()
{
  //TODO YG : Check if cluster positionREP() is valid ?

  // Here we save in each HCAL cluster the list of phi/eta values of linked clusters.
  for (BlockElt2BlockEltMap::iterator it = cluster2TargetLinks_.begin();
       it != cluster2TargetLinks_.end(); ++it) {
    reco::PFMultiLinksTC multitracks(true);

    for (BlockEltSet::iterator jt = it->second.begin();
	 jt != it->second.end(); ++jt) {

      reco::PFRecTrackRef trackref = (*jt)->trackRefPF();
      const reco::PFTrajectoryPoint& atHCAL = 
	trackref->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance);
      double tracketa = atHCAL.positionREP().Eta();
      double trackphi = atHCAL.positionREP().Phi();
      
      multitracks.linkedClusters.push_back(std::make_pair(trackphi, tracketa));
    }

    it->first->setMultilinks(multitracks);
  }

  // We set the multilinks flag of the track to true. It will allow us to 
  // use in an optimized way our algo results in the recursive linking algo.
  for (BlockEltSet::iterator it = fieldClusterSet_.begin();
       it != fieldClusterSet_.end(); ++it)
    (*it)->setIsValidMultilinks(true);

}

void
KDTreeLinkerTrackHcalAOD::clear()
{
  targetSet_.clear();
  fieldClusterSet_.clear();

  detidsSet_.clear();

  detid2ClusterLinks_.clear();
  cluster2TargetLinks_.clear();

  tree_.clear();
}
