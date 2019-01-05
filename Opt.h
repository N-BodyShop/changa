#ifndef __OPT_H__
#define __OPT_H__

/// @file Opt.h
///
/// Declare Opt classes for walk actions
///
#include "codes.h"

/// Base class for optimizing walk actions.

class Opt{
  protected:
  int action_array[2][NUM_NODE_TYPES+1];
  OptType type;
  Opt(OptType _type) : type(_type) {}
  public:

  /// Key method of Opt: given a opening decision, and a node, return
  /// an action to perform.  The action will depend on the type of
  /// node considered and which walk is being done.

  int action(int openDecision, GenericTreeNode *node){
    return action_array[openDecision][node->getType()];
  }
  OptType getSelfType() {return type;}
};

/// Class for optimizing remote gravity walk actions.
class RemoteOpt : public Opt{
  public:
  RemoteOpt() : Opt(Remote){
    // don't need to open
    //
    // local data
    action_array[0][Bucket] = DUMP;
    action_array[0][Internal] = DUMP;
    /*
    action_array[0][Boundary] = COMPUTE;
    action_array[0][NonLocal] = COMPUTE;
    action_array[0][NonLocalBucket] = COMPUTE;
    action_array[0][Cached] = COMPUTE;
    action_array[0][CachedBucket] = COMPUTE;
    */
    action_array[0][Boundary] = DUMP;
    // XXX - changed to COMPUTE from DUMP. see below also
    // CAMBRIDGE: action_array[0][NonLocal] change from COMPUTE to DUMP
    action_array[0][NonLocal] = DUMP;
    action_array[0][NonLocalBucket] = DUMP;
    
    // changed these two to COMPUTE from DEFER
    action_array[0][Cached] = COMPUTE;
    action_array[0][CachedBucket] = COMPUTE;

    action_array[0][CachedEmpty] = DUMP;

    action_array[0][Empty] = DUMP;
    action_array[0][Top] = ERROR;	
    action_array[0][Invalid] = ERROR;
     
    // need to open
    //
    // local data
    action_array[1][Internal] = DUMP;
    action_array[1][Bucket] = DUMP;
    
    // non-local data
    action_array[1][NonLocal] = KEEP;
    action_array[1][NonLocalBucket] = KEEP_REMOTE_BUCKET;
    action_array[1][Boundary] = KEEP;

    action_array[1][Cached] = KEEP;
    action_array[1][CachedBucket] = KEEP_REMOTE_BUCKET;

    // in this case, ancestors of CachedEmpty nodes must be checked as well
    action_array[1][CachedEmpty] = DUMP;

    action_array[1][Top] = ERROR; 	
    action_array[1][Invalid] = ERROR;
    action_array[1][Empty] = DUMP;

  }

};

/// Class for optimizing local gravity walk actions.
class LocalOpt : public Opt{
  public:
  LocalOpt() : Opt(Local){
    // don't need to open
    // these nodes are your concern
    action_array[0][Internal] = COMPUTE;
    action_array[0][Bucket] = COMPUTE;

    action_array[0][Boundary] = COMPUTE;
    // XXX - changed to DUMP from COMPUTE
    // CAMBRIDGE: action_array[0][NonLocal] change from DUMP to COMPUTE
    action_array[0][NonLocal] = COMPUTE;
    action_array[0][NonLocalBucket] = COMPUTE;	
    // changed these to DUMP from COMPUTE - remote does computation now
    action_array[0][Cached] = DUMP;	
    action_array[0][CachedBucket] = DUMP;

    // these nodes are no one's concern
    action_array[0][Empty] = DUMP;
    action_array[0][CachedEmpty] = DUMP;
    action_array[0][Top] = ERROR;
    action_array[0][Invalid] = ERROR;
    //--------------
    // need to open node
    // local data
    action_array[1][Internal] = KEEP;
    action_array[1][Bucket] = KEEP_LOCAL_BUCKET;
    action_array[1][Boundary] = KEEP;

    // remote data
    action_array[1][NonLocal] = DUMP;
    action_array[1][NonLocalBucket] = DUMP;
    // remote opt KEEPs Cached and KEEP_REMOTE_BUCKETs CachedBucket
    action_array[1][CachedBucket] = DUMP;
    action_array[1][Cached] = DUMP;

    // discard
    action_array[1][Empty] = DUMP;
    action_array[1][CachedEmpty] = DUMP;
    action_array[1][Top] = ERROR;
    action_array[1][Invalid] = ERROR;
  }
  
};

/// Class for optimizing experimental "push" gravity walk actions.

class PushGravityOpt : public Opt{
  public:
  PushGravityOpt() : Opt(PushGravity){
    // don't need to open node
    // these nodes are your concern
    action_array[0][Internal] = COMPUTE;
    action_array[0][Bucket] = COMPUTE;

    action_array[0][Boundary] = KEEP; 
    // XXX - changed to DUMP from COMPUTE
    action_array[0][NonLocal] = DUMP; 
    action_array[0][NonLocalBucket] = DUMP;	
    // changed these to DUMP from COMPUTE - remote does computation now
    action_array[0][Cached] = DUMP;	
    action_array[0][CachedBucket] = DUMP;

    // these nodes are no one's concern
    action_array[0][Empty] = DUMP;
    action_array[0][CachedEmpty] = DUMP;
    action_array[0][Top] = ERROR;
    action_array[0][Invalid] = ERROR;
    //--------------
    // need to open node
    // local data
    action_array[1][Internal] = KEEP;
    action_array[1][Bucket] = KEEP_LOCAL_BUCKET;
    action_array[1][Boundary] = KEEP;

    // remote data
    action_array[1][NonLocal] = DUMP;
    action_array[1][NonLocalBucket] = DUMP;
    // remote opt KEEPs Cached and KEEP_REMOTE_BUCKETs CachedBucket
    action_array[1][CachedBucket] = DUMP;
    action_array[1][Cached] = DUMP;

    // discard
    action_array[1][Empty] = DUMP;
    action_array[1][CachedEmpty] = DUMP;
    action_array[1][Top] = ERROR;
    action_array[1][Invalid] = ERROR;
  }
  
};

/// Optimization for Prefetch walk.

class PrefetchOpt : public Opt{
  public: 
  PrefetchOpt() : Opt(Pref){
    for(int i = 1; i <= NUM_NODE_TYPES; i++) 
      action_array[0][i] = DUMP;
    action_array[0][Invalid] = ERROR;
    action_array[0][Top] = ERROR;

    action_array[1][Internal] = DUMP;
    action_array[1][Bucket] = DUMP;

    action_array[1][Boundary] = KEEP;
    // we KEEP NonLocal nodes, even though we don't have their
    // children - the tw object will request the children if need be
    action_array[1][NonLocal] = KEEP;  
    action_array[1][NonLocalBucket] = KEEP_REMOTE_BUCKET;
    action_array[1][Cached] = KEEP;
    action_array[1][CachedBucket] = KEEP_REMOTE_BUCKET;

    action_array[1][CachedEmpty] = DUMP;
    action_array[1][Empty] = DUMP;
    action_array[1][Invalid] = ERROR; // node not properly initialized
    action_array[1][Top] = ERROR;    // not handled anywhere in the original code, so ERROR
  }
};
#endif
