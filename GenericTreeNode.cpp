#include "config.h"
#include "GenericTreeNode.h"

namespace Tree {

GenericTreeNode *BinaryTreeNode::clone() const {
  return new BinaryTreeNode(*this);
};

void BinaryTreeNode::pup(PUP::er &p, int depth) {
  GenericTreeNode::pup(p);
  int isNull;
  for (int i=0; i<2; ++i) {
    isNull = (children[i]==NULL || depth==0) ? 0 : 1;
    p | isNull;
    CkAssert(isNull==0 || isNull==1);
    if (isNull != 0 && depth != 0) {
      if (p.isUnpacking()) children[i] = new BinaryTreeNode();
      children[i]->pup(p, depth-1);
      if (p.isUnpacking()) children[i]->parent = this;
    }
  }
};

}
