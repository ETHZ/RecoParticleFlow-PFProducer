#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerPSEcalAOD.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "TMath.h"

KDTreeLinkerPSEcalAOD::KDTreeLinkerPSEcalAOD()
  : KDTreeLinkerBase(),
    resPSpitch_ (0.19),
    resPSlength_ (6.1),
    ps1ToEcal_ (1.072),
    ps2ToEcal_ (1.057)
{}

KDTreeLinkerPSEcalAOD::~KDTreeLinkerPSEcalAOD()
{
  clear();
}


void
KDTreeLinkerPSEcalAOD::insertTargetElt(reco::PFBlockElement	*psCluster)
{
  reco::PFClusterRef clusterref = psCluster->clusterRef();

  // This test is more or less done in PFBlockAlgo.h. In others cases, it should be switch on.
  //   if (!((clusterref->layer() == PFLayer::PS1) || 
  // 	(clusterref->layer() == PFLayer::PS2)))
  //     return;

  targetSet_.insert(psCluster);
}


void
KDTreeLinkerPSEcalAOD::insertFieldClusterElt(reco::PFBlockElement	*ecalCluster)
{
  reco::PFClusterRef clusterref = ecalCluster->clusterRef();

  if (clusterref->layer() != PFLayer::ECAL_ENDCAP)
    return;

  const std::vector< std::pair<DetId, float> >& hitsAndFracs = clusterref->hitsAndFractions();

  // We create a list of cluster
  fieldClusterSet_.insert(ecalCluster);

  double clusterz = clusterref->position().Z();
  DetIdSet& detidsSet = (clusterz < 0) ? detidsNegSet_ : detidsPosSet_;

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
    detidsSet.insert(&detid);
  } 

}

void 
KDTreeLinkerPSEcalAOD::buildTree()
{
  buildTree(detidsNegSet_, treeNeg_);
  buildTree(detidsPosSet_, treePos_);
}

void 
KDTreeLinkerPSEcalAOD::buildTree(const DetIdSet	&detidsSet,
			      KDTreeLinkerAlgoAOD	&tree)
{
  // List of pseudo-rechits that will be used to create the KDTree
  std::vector<KDTreeNodeInfoAOD> eltList;

  // Filling of this eltList
  for(DetIdSet::const_iterator it = detidsSet.begin(); 
      it != detidsSet.end(); it++) {
   
	 // here we get the position of the cell, and convert to the formats in PF for now (AA)
    const CaloCellGeometry* cellGeometry = caloGeometry_->getGeometry(**it);
	 // IMPORTANT: see in the PFRecHitProducer how exactly is the position defined  <---------------- !!!!!
	 // this is temporary and uses the face of the cel (AA)
	 const GlobalPoint gpPosition = cellGeometry->getPosition();
	 const math::XYZPoint posxyz(gpPosition.x(), gpPosition.y(), gpPosition.z());

    KDTreeNodeInfoAOD diInfo (*it, posxyz.X(), posxyz.Y());
    eltList.push_back(diInfo);
  }

  // xmin-xmax, ymain-ymax
  KDTreeBox region(-150., 150., -150., 150.);

  // We may now build the KDTree
  tree.build(eltList, region);
}

void
KDTreeLinkerPSEcalAOD::searchLinks()
{
  // Must of the code has been taken from LinkByRecHit.cc

  // We iterate over the PS clusters.
  for(BlockEltSet::iterator it = targetSet_.begin(); 
      it != targetSet_.end(); it++) {

    (*it)->setIsValidMultilinks(true);
	
    reco::PFClusterRef clusterPSRef = (*it)->clusterRef();
    const reco::PFCluster& clusterPS = *clusterPSRef;

    // PS cluster position, extrapolated to ECAL
    double zPS = clusterPS.position().Z();
    double xPS = clusterPS.position().X();
    double yPS = clusterPS.position().Y();

    double etaPS = fabs(clusterPS.positionREP().eta());
    double deltaX = 0.;
    double deltaY = 0.;
    double xPSonEcal = xPS;
    double yPSonEcal = yPS;

    if (clusterPS.layer() == PFLayer::PS1) { // PS1

      // vertical strips, measure x with pitch precision
      deltaX = resPSpitch_;
      deltaY = resPSlength_;
      xPSonEcal *= ps1ToEcal_;
      yPSonEcal *= ps1ToEcal_;

    } else { // PS2

      // horizontal strips, measure y with pitch precision
      deltaY = resPSpitch_;
      deltaX = resPSlength_;
      xPSonEcal *= ps2ToEcal_;
      yPSonEcal *= ps2ToEcal_;

    }
 
    
    // Estimate the maximal envelope in phi/eta that will be used to find rechit candidates.
    // Same envelope for cap et barrel rechits.
    
    
    double maxEcalRadius = getCristalXYMaxSize() / 2.;

    // The inflation factor includes the approximate projection from Preshower to ECAL
    double inflation = 2.4 - (etaPS-1.6);
    double rangeX = maxEcalRadius * (1 + (0.05 + 1.0 / maxEcalRadius * deltaX / 2.)) * inflation; 
    double rangeY = maxEcalRadius * (1 + (0.05 + 1.0 / maxEcalRadius * deltaY / 2.)) * inflation; 
    
    // We search for all candidate detIds, ie all detids of cells contained in the maximal size envelope.
    std::vector<KDTreeNodeInfoAOD> detIds;
    KDTreeBox trackBox(xPSonEcal - rangeX, xPSonEcal + rangeX, 
		  yPSonEcal - rangeY, yPSonEcal + rangeY);

    if (zPS < 0)
      treeNeg_.search(trackBox, detIds);
    else
      treePos_.search(trackBox, detIds);


    for(std::vector<KDTreeNodeInfoAOD>::const_iterator dId = detIds.begin(); 
	dId != detIds.end(); ++dId) {
 
     // get the position and corners - check for depth corrections (AA) <---------------------------------- !!!!
		const CaloCellGeometry* cellGeometry = caloGeometry_->getGeometry(*(dId->ptr));
		const GlobalPoint gpPosition = cellGeometry->getPosition();
		CaloCellGeometry::CornersVec cellGeoCorners = cellGeometry->getCorners();
		// convert to types used in the PF framework
		// this unscaled position is scaled by zPS/clusterz below
		const math::XYZPoint posxyzU(gpPosition.x(), gpPosition.y(), gpPosition.z()); 
		std::vector< math::XYZPoint > corners;
	   for (uint i=0; i<4; ++i) {
			corners.push_back(math::XYZPoint(cellGeoCorners[3-i].x(),
														cellGeoCorners[3-i].y(),
														cellGeoCorners[3-i].z()
														));
		} // end of corner re-ordering according to the convention in PFRecHit

     if(corners.size() != 4) continue;

      // Find all clusters associated to a given detid
      DetId2BlockEltMap::iterator ret = detid2ClusterLinks_.find(dId->ptr);
      
      for(BlockEltSet::const_iterator clusterIt = ret->second.begin(); 
	  clusterIt != ret->second.end(); clusterIt++) {
	
	reco::PFClusterRef clusterref = (*clusterIt)->clusterRef();
	double clusterz = clusterref->position().Z();

	// scale the unscaled position obtained from the geometry
	const math::XYZPoint& posxyz = posxyzU * zPS / clusterz;

	double x[5];
	double y[5];
	for ( unsigned jc=0; jc<4; ++jc ) {
	  math::XYZPoint cornerpos = corners[jc] * zPS / clusterz;
	  x[jc] = cornerpos.X() + (cornerpos.X()-posxyz.X()) * (0.05 +1.0/fabs((cornerpos.X()-posxyz.X()))*deltaX/2.);
	  y[jc] = cornerpos.Y() + (cornerpos.Y()-posxyz.Y()) * (0.05 +1.0/fabs((cornerpos.Y()-posxyz.Y()))*deltaY/2.);
	}

	x[4] = x[0];
	y[4] = y[0];
	
	bool isinside = TMath::IsInside(xPS,
					yPS,
					5,x,y);
	
	// Check if the track and the cluster are linked
	if( isinside )
	  target2ClusterLinks_[*it].insert(*clusterIt);	
      }
    }
    
  }
}

void
KDTreeLinkerPSEcalAOD::updatePFBlockEltWithLinks()
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
KDTreeLinkerPSEcalAOD::clear()
{
  targetSet_.clear();
  fieldClusterSet_.clear();

  detidsNegSet_.clear();
  detidsPosSet_.clear();

  detid2ClusterLinks_.clear();
  target2ClusterLinks_.clear();

  treeNeg_.clear();
  treePos_.clear();
}
