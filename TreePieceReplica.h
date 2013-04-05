#ifndef TREEPIECE_REPLICA_H
#define TREEPIECE_REPLICA_H

#include "ParallelGravity.decl.h"

class TreePieceReplica : public CBase_TreePieceReplica {
	private:
  	typedef std::map<NodeKey, GenericTreeNode *> NodeLookupType;
		NodeLookupType nodeLookupTable;

	public:
		TreePieceReplica();
		TreePieceReplica(CkMigrateMessage *);

		void fillRequestNodeFromReplica(CkCacheRequestMsg<KeyType> *msg);

		void recvTreePiece(CkCacheFillMsg<KeyType> *msg);

		const GenericTreeNode *lookupNode(Tree::NodeKey key){
			return keyToNode(key);
		};

		inline GenericTreeNode *keyToNode(const Tree::NodeKey k){
			NodeLookupType::iterator iter = nodeLookupTable.find(k);
			if (iter != nodeLookupTable.end()) return iter->second;
			else return NULL;
		}

		void fillNodeLookupTable(Tree::BinaryTreeNode *node);
		void clearTable(CkCallback &cb);
};


#endif

