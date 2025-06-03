#ifndef TREEPIECE_REPLICA_H
#define TREEPIECE_REPLICA_H

#include "ParallelGravity.decl.h"

class TreeReplicaMsg : public CMessage_TreeReplicaMsg {
public:
    NodeKey key;                ///< Key of the root of this replica
    char *data;                 ///< Packed nodes of the replica
    int iFrom;                  ///< The TreePiece sending this replica
    TreeReplicaMsg (NodeKey k, int i) : key(k), iFrom(i) {}
};

class TreePieceReplica : public CBase_TreePieceReplica {
	private:
  	typedef std::map<NodeKey, GenericTreeNode *> NodeLookupType;
		NodeLookupType nodeLookupTable;
        std::vector<TreeReplicaMsg *> msgRecvd; ///< Hold received message
                                                /// pointers for memory management

	public:
		TreePieceReplica();
		TreePieceReplica(CkMigrateMessage *);

		void fillRequestNodeFromReplica(CkCacheRequestMsg<KeyType> *msg);

            void recvTreePiece(TreeReplicaMsg *msg);

		const GenericTreeNode *lookupNode(Tree::NodeKey key){
			return keyToNode(key);
		};

		inline GenericTreeNode *keyToNode(const Tree::NodeKey k){
			NodeLookupType::iterator iter = nodeLookupTable.find(k);
			if (iter != nodeLookupTable.end()) return iter->second;
			else return NULL;
		}

		void fillNodeLookupTable(Tree::BinaryTreeNode *node);
		void clearTable(const CkCallback &cb);
};


#endif

