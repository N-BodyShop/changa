#ifndef __CACHEINTERFACE_H__
#define __CACHEINTERFACE_H__

/** @file CacheInterface.h
 *
 *  Declares the interfaces used by the CacheManager: the software
 *  cache for requesting off processor particle and node data.
 */

#include <CkCache.h>
#if CHANGA_SMPCACHE
#include <CkTreeCache.h>
#include <atomic>
#endif
#include "config.h"
#include "gravity.h"
#include "GenericTreeNode.h"
#include "keytype.h"

/*********************************************************
 * Gravity interface: Particles
 *********************************************************/

/// @brief The data in a GravityParticle cache entry.
class CacheParticle {
public:
  /// Message containing the data for this entry.
  CkCacheFillMsg<KeyType> *msg;
  /// Index of the first particle in the home processor's myParticles array.
  int begin;
  /// Index of the last particle in the home processor's myParticles array.
  int end;
  /// The rest of the structure is an array of particles.  Declared as
  /// length 1, but can be arbitrary length.  It is assumed that these
  /// particles are contiguous in the myParticles array.
  ExternalGravityParticle part[1];
};

/// @brief Cache interface to particles for the gravity calculation.
/// This is a read-only cache of particles.
class EntryTypeGravityParticle : public CkCacheEntryType<KeyType> {
public:
  EntryTypeGravityParticle();
  /// @brief Request a bucket of particles from a TreePiece.
  void * request(CkArrayIndexMax&, KeyType);
  /// @brief Return data from fufilled cache request.
  void * unpack(CkCacheFillMsg<KeyType> *, int, CkArrayIndexMax &);
  /// @brief Do nothing: this is a read-only cache.
  void writeback(CkArrayIndexMax&, KeyType, void *);
  /// @brief free cached data.
  void free(void *);
  /// @brief return size of cached data.
  int size(void *);
  
  /// @brief callback to TreePiece after data is received.
  static void callback(CkArrayID, CkArrayIndexMax&, KeyType, CkCacheUserData &, void*, int);
};

/*********************************************************
 * Smooth interface: Particles
 *********************************************************/

/// @brief particle data in the smooth particle cache messages
class CacheSmoothParticle {
public:
    int begin;  ///< Beginning particle number
    int end;    ///< ending Particle number
    int nActual; ///< actual number of particles sent
    KeyType key; ///< Key of this bucket (for writeback)
    GravityParticle *partCached;        ///< particle data
    extraSPHData *extraSPHCached;       ///< particle extraData
    ExternalSmoothParticle partExt[1];  ///< particle data in the message
};

/// @brief Cache interface to the particles for smooth calculations.
/// This cache is a writeback cache.
class EntryTypeSmoothParticle : public CkCacheEntryType<KeyType> {
    // N.B. can't have helpful attributes because of the static function.
public:
  EntryTypeSmoothParticle();
  /// @brief Request a bucket of particles from a TreePiece.
  void * request(CkArrayIndexMax&, KeyType);
  /// @brief Return data from fufilled cache request.
  void * unpack(CkCacheFillMsg<KeyType> *, int, CkArrayIndexMax &);
  void writeback(CkArrayIndexMax&, KeyType, void *);
  /// @brief free cached data.
  void free(void *);
  /// @brief return size of cached data.
  int size(void *);
  
  /// @brief callback to TreePiece after data is received.
  static void callback(CkArrayID, CkArrayIndexMax&, KeyType, CkCacheUserData &, void*, int);
};

/*********************************************************
 * Gravity interface: Nodes
 *********************************************************/

/// @brief Cache interface to the Tree Nodes.
class EntryTypeGravityNode : public CkCacheEntryType<KeyType> {
  void *vptr; // For saving a copy of the virtual function table.
	      // It's use will be compiler dependent.
  void unpackSingle(CkCacheFillMsg<KeyType> *, Tree::BinaryTreeNode *, int, CkArrayIndexMax &, bool);
public:
  EntryTypeGravityNode();
  void * request(CkArrayIndexMax&, KeyType);
  void * unpack(CkCacheFillMsg<KeyType> *, int, CkArrayIndexMax &);
  void writeback(CkArrayIndexMax&, KeyType, void *);
  void free(void *);
  int size(void *);
  
  static void callback(CkArrayID, CkArrayIndexMax&, KeyType, CkCacheUserData &, void*, int);
};

#if CHANGA_SMPCACHE
/*********************************************************
 * Process-shared node cache (CkTreeCacheManager, charm ck-libs/cache)
 *********************************************************/

/// @brief Placeholder the node cache puts in an empty child slot of the
/// DataManager's merged tree while the node is being fetched. Type stays
/// Invalid, so BinaryTreeNode::getChildren() reports it as NULL to the
/// tree walks; only the cache manager ever dereferences it.
class NodeCachePlaceholder : public Tree::BinaryTreeNode {
public:
  std::atomic<void*> parked_head{nullptr};  ///< parked requestors (CkTreeCache)
  std::atomic<int> request_latch{0};        ///< first PE to set it sends the request
  std::atomic<Tree::BinaryTreeNode*> replacement{nullptr}; ///< the node that took its place
  int chunk;
  NodeCachePlaceholder(Tree::NodeKey k, Tree::BinaryTreeNode *p, int c) : chunk(c) {
    key = k;
    parent = p;
  }
};

/// @brief Traits binding CkTreeCacheManager to ChaNGa's BinaryTreeNode
/// (contract documented in CkTreeCache.h and TreeCacheCore.h).
struct GravityNodeTraits {
  typedef Tree::BinaryTreeNode Node;
  typedef KeyType Key;
  static Key key(const Node *n) { return n->getKey(); }
  static Node *parent(const Node *n) { return (Node *)n->parent; }
  static void setParent(Node *n, Node *p) { n->parent = p; }
  static int branchFactor() { return 2; }
  static Node *rawChild(const Node *n, int i) { return __atomic_load_n(&n->children[i], __ATOMIC_ACQUIRE); }
  static Node *exchangeChild(Node *n, int i, Node *c) { return __atomic_exchange_n(&n->children[i], c, __ATOMIC_ACQ_REL); }
  static bool casChild(Node *n, int i, Node *&expected, Node *desired) {
    return __atomic_compare_exchange_n(&n->children[i], &expected, desired, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE);
  }
  static void wireChild(Node *n, int i, Node *c) { n->children[i] = c; }
  static bool canHaveChildren(const Node *n) {
    switch (n->getType()) {
    case Tree::Bucket: case Tree::NonLocalBucket: case Tree::CachedBucket:
    case Tree::Empty: case Tree::CachedEmpty:
      return false;
    default:
      return true;
    }
  }
  static bool isPlaceholder(const Node *n) { return n->getType() == Tree::Invalid; }
  static Node *makePlaceholder(Key k, Node *p, int chunk) { return new NodeCachePlaceholder(k, p, chunk); }
  static void freePlaceholder(Node *n) { delete static_cast<NodeCachePlaceholder *>(n); }
  static int placeholderChunk(const Node *n) { return static_cast<const NodeCachePlaceholder *>(n)->chunk; }
  static std::atomic<void*> &parkedHead(Node *n) { return static_cast<NodeCachePlaceholder *>(n)->parked_head; }
  static std::atomic<int> &requestLatch(Node *n) { return static_cast<NodeCachePlaceholder *>(n)->request_latch; }
  static std::atomic<Node*> &replacement(Node *n) { return static_cast<NodeCachePlaceholder *>(n)->replacement; }
  static Node *root();   ///< the DataManager's merged tree (CacheInterface.cpp)
};

typedef CProxy_CkTreeCacheManager<KeyType, GravityNodeTraits> CProxy_NodeCache;
#else
/// The node cache: charm's per-PE CkCacheManager (configure --enable-smpcache=no).
typedef CProxy_CkCacheManager<KeyType> CProxy_NodeCache;
#endif

#endif


