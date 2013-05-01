#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerAlgoAOD.h"

KDTreeLinkerAlgoAOD::KDTreeLinkerAlgoAOD()
  : root_ (0),
    nodePool_(0),
    nodePoolSize_(-1),
    nodePoolPos_(-1)
{
}

KDTreeLinkerAlgoAOD::~KDTreeLinkerAlgoAOD()
{
  clear();
}

void
KDTreeLinkerAlgoAOD::build(std::vector<KDTreeNodeInfoAOD>	&eltList, 
			const KDTreeBox			&region)
{
  if (eltList.size()) {
    nodePoolSize_ = eltList.size() * 2 - 1;
    nodePool_ = new KDTreeNodeAOD[nodePoolSize_];

    // Here we build the KDTree
    root_ = recBuild(eltList, 0, eltList.size(), 0, region);
  }
}


KDTreeNodeAOD*
KDTreeLinkerAlgoAOD::recBuild(std::vector<KDTreeNodeInfoAOD>	&eltList, 
			   int				low, 
			   int				high, 
			   int				depth,
			   const KDTreeBox&		region)
{
  int portionSize = high - low;

  // By construction, portionSize > 0 can't happend.
  assert(portionSize > 0);

  if (portionSize == 1) { // Leaf case
   
    KDTreeNodeAOD *leaf = getNextNode();
    leaf->setAttributs(region, eltList[low]);
    return leaf;

  } else { // Node case
    
    // The even depth is associated to dim1 dimension
    // The odd one to dim2 dimension
    int medianId = medianSearch(eltList, low, high, depth);

    // We create the node
    KDTreeNodeAOD *node = getNextNode();
    node->setAttributs(region);


    // Here we split into 2 halfplanes the current plane
    KDTreeBox leftRegion = region;
    KDTreeBox rightRegion = region;
    if (depth & 1) {

      double medianVal = eltList[medianId].dim2;
      leftRegion.dim2max = medianVal;
      rightRegion.dim2min = medianVal;

    } else {

      double medianVal = eltList[medianId].dim1;
      leftRegion.dim1max = medianVal;
      rightRegion.dim1min = medianVal;

    }

    ++depth;
    ++medianId;

    // We recursively build the son nodes
    node->left = recBuild(eltList, low, medianId, depth, leftRegion);
    node->right = recBuild(eltList, medianId, high, depth, rightRegion);

    return node;
  }
}


//Fast median search with Wirth algorithm in eltList between low and high indexes.
int
KDTreeLinkerAlgoAOD::medianSearch(std::vector<KDTreeNodeInfoAOD>	&eltList,
			       int				low,
			       int				high,
			       int				treeDepth)
{
  //We should have at least 1 element to calculate the median...
  assert(low < high);

  int nbrElts = high - low;
  int median = (nbrElts & 1)	? nbrElts / 2 
				: nbrElts / 2 - 1;
  median += low;

  int l = low;
  int m = high - 1;
  
  while (l < m) {
    KDTreeNodeInfoAOD elt = eltList[median];
    int i = l;
    int j = m;

    do {

      // The even depth is associated to dim1 dimension
      // The odd one to dim2 dimension
      if (treeDepth & 1) {
	while (eltList[i].dim2 < elt.dim2) i++;
	while (eltList[j].dim2 > elt.dim2) j--;
      } else {
	while (eltList[i].dim1 < elt.dim1) i++;
	while (eltList[j].dim1 > elt.dim1) j--;
      }

      if (i <= j){
	swap(eltList[i], eltList[j]);
	i++; 
	j--;
      }
    } while (i <= j);
    if (j < median) l = i;
    if (i > median) m = j;
  }

  return median;
}

void 
KDTreeLinkerAlgoAOD::swap(KDTreeNodeInfoAOD	&e1, 
		       KDTreeNodeInfoAOD	&e2)
{
  KDTreeNodeInfoAOD tmp = e1;
  e1 = e2;
  e2 = tmp;
}

void
KDTreeLinkerAlgoAOD::search(const KDTreeBox		&trackBox,
			 std::vector<KDTreeNodeInfoAOD>	&detIds)
{
  if (root_)
    recSearch(root_, trackBox, detIds);
}

void 
KDTreeLinkerAlgoAOD::recSearch(const KDTreeNodeAOD		*current,
			    const KDTreeBox		&trackBox,
			    std::vector<KDTreeNodeInfoAOD>	&detIds)
{
  // By construction, current can't be null
  assert(current != 0);

  // By Construction, a node can't have just 1 son.
  assert (!(((current->left == 0) && (current->right != 0)) ||
	    ((current->left != 0) && (current->right == 0))));
    
  if ((current->left == 0) && (current->right == 0)) {//leaf case
  
    // If point inside the rectangle/area
    if ((current->di.dim1 >= trackBox.dim1min) && (current->di.dim1 <= trackBox.dim1max) &&
	(current->di.dim2 >= trackBox.dim2min) && (current->di.dim2 <= trackBox.dim2max))
      detIds.push_back(current->di);

  } else {

    //if region( v->left ) is fully contained in the rectangle
    if ((current->left->region.dim1min >= trackBox.dim1min) && 
	(current->left->region.dim1max <= trackBox.dim1max) &&
	(current->left->region.dim2min >= trackBox.dim2min) && 
	(current->left->region.dim2max <= trackBox.dim2max))
      addSubtree(current->left, detIds);
    
    else { //if region( v->left ) intersects the rectangle
      
      if (!((current->left->region.dim1min >= trackBox.dim1max) || 
	    (current->left->region.dim1max <= trackBox.dim1min) ||
	    (current->left->region.dim2min >= trackBox.dim2max) || 
	    (current->left->region.dim2max <= trackBox.dim2min)))
	recSearch(current->left, trackBox, detIds);
    }
    
    //if region( v->right ) is fully contained in the rectangle
    if ((current->right->region.dim1min >= trackBox.dim1min) && 
	(current->right->region.dim1max <= trackBox.dim1max) &&
	(current->right->region.dim2min >= trackBox.dim2min) && 
	(current->right->region.dim2max <= trackBox.dim2max))
      addSubtree(current->right, detIds);

    else { //if region( v->right ) intersects the rectangle
     
      if (!((current->right->region.dim1min >= trackBox.dim1max) || 
	    (current->right->region.dim1max <= trackBox.dim1min) ||
	    (current->right->region.dim2min >= trackBox.dim2max) || 
	    (current->right->region.dim2max <= trackBox.dim2min)))
	recSearch(current->right, trackBox, detIds);
    } 
  }
}

void
KDTreeLinkerAlgoAOD::addSubtree(const KDTreeNodeAOD		*current, 
		   std::vector<KDTreeNodeInfoAOD>	&detIds)
{
  // By construction, current can't be null
  assert(current != 0);

  if ((current->left == 0) && (current->right == 0)) // leaf
    detIds.push_back(current->di);
  else { // node
    addSubtree(current->left, detIds);
    addSubtree(current->right, detIds);
  }
}


void 
KDTreeLinkerAlgoAOD::clearTree()
{
  delete[] nodePool_;
  nodePool_ = 0;
  root_ = 0;
  nodePoolSize_ = -1;
  nodePoolPos_ = -1;
}

void 
KDTreeLinkerAlgoAOD::clear()
{
  if (root_)
    clearTree();
}


KDTreeNodeAOD* 
KDTreeLinkerAlgoAOD::getNextNode()
{
  ++nodePoolPos_;

  // The tree size is exactly 2 * nbrElts - 1 and this is the total allocated memory.
  // If we have used more than that....there is a big problem.
  assert(nodePoolPos_ < nodePoolSize_);

  return &(nodePool_[nodePoolPos_]);
}
