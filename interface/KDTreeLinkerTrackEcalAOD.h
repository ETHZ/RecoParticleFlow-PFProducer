#ifndef KDTreeLinkerTrackEcalAOD_h
#define KDTreeLinkerTrackEcalAOD_h

#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerBase.h"
//#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerToolsAOD.h"
#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerTools.h"
#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerAlgoAOD.h"


// This class is used to find all links between Tracks and ECAL clusters
// using a KDTree algorithm.
// It is used in PFBlockAlgo.cc in the function links().
class KDTreeLinkerTrackEcalAOD : public KDTreeLinkerBase
{
 public:
  KDTreeLinkerTrackEcalAOD();
  ~KDTreeLinkerTrackEcalAOD();

  // set the geometry (AA)
  void setCaloGeometry(const CaloGeometry* caloGeo) { caloGeometry_ = caloGeo; }

  // With this method, we create the list of psCluster that we want to link.
  void insertTargetElt(reco::PFBlockElement		*track);

  // Here, we create the list of ecalCluster that we want to link. From ecalCluster
  // and fraction, we will create a second list of detids that will be used to
  // build the KDTree.
  void insertFieldClusterElt(reco::PFBlockElement	*ecalCluster);  

  // The KDTree building from detids list.
  void buildTree();
  
  // Here we will iterate over all tracks. For each track intersection point with ECAL, 
  // we will search the closest detids in the KDTree, from detids we will find the 
  // ecalClusters and after that we will check the links between the track and 
  // all closest ecalClusters.  
  void searchLinks();
    
  // Here, we will store all PS/ECAL founded links in the PFBlockElement class
  // of each psCluster in the PFmultilinks field.
  void updatePFBlockEltWithLinks();
  
  // Here we free all allocated structures.
  void clear();
 
 private:
  // Data used by the KDTree algorithm : sets of Tracks and ECAL clusters.
  BlockEltSet		targetSet_;
  BlockEltSet		fieldClusterSet_;

  // Sets of deids that compose the ECAL clusters. 
  DetIdSet		detidsSet_;
  
  // Map of linked Track/ECAL clusters.
  BlockElt2BlockEltMap	target2ClusterLinks_;

  // Map of the ECAL clusters associated to a detid.
  DetId2BlockEltMap	detid2ClusterLinks_;
    
  // KD trees
  KDTreeLinkerAlgoAOD	tree_;

  // calorimeter geometry to get cell coordinates (AA)
  const CaloGeometry* caloGeometry_;

};


#endif /* !KDTreeLinkerTrackEcalAOD_h */
