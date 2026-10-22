#include "ParallelGravity.h"
#include "TreePieceReplica.h"

TreePieceReplica::TreePieceReplica() {
}

TreePieceReplica::TreePieceReplica(CkMigrateMessage *m) : CBase_TreePieceReplica(m) {
}

void TreePieceReplica::fillRequestNodeFromReplica(
		CkCacheRequestMsg<KeyType> *msg) {
  const Tree::GenericTreeNode* node = lookupNode(msg->key);
	//CkPrintf("Recvd fillrequest node from replica\n");
  if(node != NULL) {
    if(_cache) {
#if 1 || defined CACHE_BUFFER_MSGS
      int count = ((Tree::BinaryTreeNode*)node)->countDepth(_cacheLineDepth);
      // PAD_reply extra bytes are allocated in front of each node to
      // store the msg pointer, and packNodes() puts the nodes
      // ALIGN_DEFAULT apart; same layout as TreePiece::fillRequestNode().
      // See the free() and the unpackSingle() method in CacheInterface.cpp.
      CkCacheFillMsg<KeyType> *reply = new (count * ALIGN_DEFAULT(sizeof(Tree::BinaryTreeNode)+PAD_reply)) CkCacheFillMsg<KeyType>(msg->key);
      ((Tree::BinaryTreeNode*)node)->packNodes((Tree::BinaryTreeNode*)(reply->data+PAD_reply), _cacheLineDepth, PAD_reply);

#else
      PUP::sizer p1;
      node->pup(p1, msg->depth);
      FillNodeMsg *reply = new (p1.size(), 0) FillNodeMsg(thisIndex);

      /// @TODO: check that at destination of "remoteIndex" are correct
      PUP::toMem p2((void*)reply->nodes);
      node->pup(p2, msg->depth);
#endif
      cacheNode[msg->replyTo].recvData(reply);
    } else {
      CkAbort("Non cached version not anymore supported, feel free to fix it!");
    }
    delete msg;
  }
  else {	// Handle NULL nodes
      // Somehow, we don't don't own this node.  Fall back to using
      // the TreePieces to find it.
      // Go to the datamanager to find a local treepiece to hand this
      // request off to.
      //
      // DataManager *dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
      // TreePiece *tpLocal = dm->registeredTreePieces[0].treePiece;
      //
      // Unfortunately, registeredTreepieces gets cleared immediately
      // after the tree is built.  Until that is fixed, lets pick a
      // random piece (I know 0 exists) to handle this node.
      treeProxy[0].fillRequestNode(msg);
  }
}

void TreePieceReplica::recvTreePiece(TreeReplicaMsg *msg) {
	//CkPrintf("Recv TP %d at PE %d\n", msg->idx, CkMyPe());
	Tree::BinaryTreeNode *node = (Tree::BinaryTreeNode *) ((char *)msg->data);
	node->unpackNodes();
	treeProxy[msg->iFrom].recvAck();
	fillNodeLookupTable(node);
        msgRecvd.push_back(msg);
}

void TreePieceReplica::fillNodeLookupTable(Tree::BinaryTreeNode *node) {
	NodeLookupType::iterator iter = nodeLookupTable.find(node->getKey());
	if (iter == nodeLookupTable.end()) {
		nodeLookupTable[node->getKey()] = node;
	}
	if (node->children[0] != NULL) {
		fillNodeLookupTable(node->children[0]);
	}
	if (node->children[1] != NULL) {
		fillNodeLookupTable(node->children[1]);
	}
}

void TreePieceReplica::clearTable(bool active, const CkCallback &cb) {
        for(int i = 0; i < msgRecvd.size(); i++)
            delete msgRecvd[i];
        msgRecvd.clear();
	nodeLookupTable.clear();
	bActive = active;
	contribute(cb);
}
