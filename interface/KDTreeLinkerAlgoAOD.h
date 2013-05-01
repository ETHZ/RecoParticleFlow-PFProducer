#ifndef KDTreeLinkerAlgoAOD_h
#define KDTreeLinkerAlgoAOD_h

//#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerToolsAOD.h"
#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerTools.h"


#include <vector>

// Class that implements the KDTree partition of 2D space and 
// a closest point search algorithme.
class KDTreeLinkerAlgoAOD
{
 public:
  KDTreeLinkerAlgoAOD();
  
  // Dtor calls clear()
  ~KDTreeLinkerAlgoAOD();
  
  // Here we build the KD tree from the "eltList" in the space define by "region".
  void build(std::vector<KDTreeNodeInfoAOD>	&eltList,
	     const KDTreeBox			&region);
  
  // Here we search in the KDTree for all points that would be 
  // contained in the given searchbox. The founded points are stored in resRecHitList.
  void search(const KDTreeBox			&searchBox,
	      std::vector<KDTreeNodeInfoAOD>	&resDetIdList);
  
  // This method clears all allocated structures.
  void clear();
  
 private:
  // The KDTree root
  KDTreeNodeAOD*	root_;
  
  // The node pool allow us to do just 1 call to new for each tree building.
  KDTreeNodeAOD*	nodePool_;
  int		nodePoolSize_;
  int		nodePoolPos_;

 private:
  // Basic swap function.
  void swap(KDTreeNodeInfoAOD &e1, KDTreeNodeInfoAOD &e2);

  // Get next node from the node pool.
  KDTreeNodeAOD* getNextNode();

  //Fast median search with Wirth algorithm in eltList between low and high indexes.
  int medianSearch(std::vector<KDTreeNodeInfoAOD>	&eltList,
		   int				low,
		   int				high,
		   int				treeDepth);

  // Recursif kdtree builder. Is called by build()
  KDTreeNodeAOD *recBuild(std::vector<KDTreeNodeInfoAOD>	&eltList, 
		       int				low,
		       int				hight,
		       int				depth,
		       const KDTreeBox&			region);

  // Recursif kdtree search. Is called by search()
  void recSearch(const KDTreeNodeAOD		*current,
		 const KDTreeBox		&trackBox,
		 std::vector<KDTreeNodeInfoAOD>	&detIds);    

  // Add all elements of an subtree to the closest elements. Used during the recSearch().
  void addSubtree(const KDTreeNodeAOD		*current, 
		  std::vector<KDTreeNodeInfoAOD>	&detIds);

  // This method frees the KDTree.     
  void clearTree();
};

#endif
