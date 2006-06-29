#include "GenericTreeNode.h"

namespace Tree {

GenericTreeNode *BinaryTreeNode::createNew() const {
  return new BinaryTreeNode();
};

GenericTreeNode *BinaryTreeNode::clone() const {
  //BinaryTreeNode *tmp = new BinaryTreeNode();
  //*tmp = *this;
  //return tmp;
  return new BinaryTreeNode(*this);
};

void BinaryTreeNode::pup(PUP::er &p, int depth) {
  //CkPrintf("Pupper of BinaryTreeNode(%d) called for %s (%d)\n",depth,p.isPacking()?"Packing":p.isUnpacking()?"Unpacking":"Sizing",p.isSizing()?((PUP::sizer*)&p)->size():((PUP::mem*)&p)->size());
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

int NodeKeyClass::staticCompare(const void *k1,const void *k2,size_t ){
  return ((const NodeKeyClass*)k1)->compare(*(const NodeKeyClass*)k2);
}

CkHashCode NodeKeyClass::staticHash(const void *v,size_t){
  return ((const NodeKeyClass*)v)->hash();
}

}
