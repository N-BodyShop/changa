#define NUM_NODE_TYPES 11

class Opt{
  protected:
  int action_array[2][NUM_NODE_TYPES+1];
  OptType type;
  Opt(OptType _type) : type(_type) {}
  public:

  int action(bool openDecision, GenericTreeNode *node){
    return action_array[openDecision][node->getType()];
  }
  OptType getSelfType() {return type;}
};

class RemoteOpt : public Opt{
  public:
  RemoteOpt() : Opt(Remote){
    // don't need to open
    //
    // local data
    action_array[false][Bucket] = DUMP;
    action_array[false][Internal] = DUMP;
    /*
    action_array[false][Boundary] = COMPUTE;
    action_array[false][NonLocal] = COMPUTE;
    action_array[false][NonLocalBucket] = COMPUTE;
    action_array[false][Cached] = COMPUTE;
    action_array[false][CachedBucket] = COMPUTE;
    */
    action_array[false][Boundary] = DUMP;
    action_array[false][NonLocal] = DUMP;
    action_array[false][NonLocalBucket] = DUMP;
    
    // changed these two to COMPUTE from DEFER
    action_array[false][Cached] = COMPUTE;
    action_array[false][CachedBucket] = COMPUTE;

    action_array[false][CachedEmpty] = DUMP;

    action_array[false][Empty] = DUMP;
    action_array[false][Top] = ERROR;	
    action_array[false][Invalid] = ERROR;
     
    // need to open
    //
    // local data
    action_array[true][Internal] = DUMP;
    action_array[true][Bucket] = DUMP;
    
    // non-local data
    action_array[true][NonLocal] = KEEP;
    action_array[true][NonLocalBucket] = KEEP_REMOTE_BUCKET;
    action_array[true][Boundary] = KEEP;

    action_array[true][Cached] = KEEP;
    action_array[true][CachedBucket] = KEEP_REMOTE_BUCKET;

    // in this case, ancestors of CachedEmpty nodes must be checked as well
    action_array[true][CachedEmpty] = DUMP;

    action_array[true][Top] = ERROR; 	
    action_array[true][Invalid] = ERROR;
    action_array[true][Empty] = DUMP;

  }

};

class LocalOpt : public Opt{
  public:
  LocalOpt() : Opt(Local){
    // don't need to open
    // these nodes are your concern
    action_array[false][Internal] = COMPUTE;
    action_array[false][Bucket] = COMPUTE;

    action_array[false][Boundary] = COMPUTE; 
    action_array[false][NonLocal] = COMPUTE; 
    action_array[false][NonLocalBucket] = COMPUTE;	
    // changed these to DUMP from COMPUTE - remote does computation now
    action_array[false][Cached] = DUMP;	
    action_array[false][CachedBucket] = DUMP;

    // these nodes are no one's concern
    action_array[false][Empty] = DUMP;
    action_array[false][CachedEmpty] = DUMP;
    action_array[false][Top] = ERROR;
    action_array[false][Invalid] = ERROR;
    //--------------
    // need to open node
    // local data
    action_array[true][Internal] = KEEP;
    action_array[true][Bucket] = KEEP_LOCAL_BUCKET;
    action_array[true][Boundary] = KEEP;

    // remote data
    action_array[true][NonLocal] = DUMP;
    action_array[true][NonLocalBucket] = DUMP;
    // remote opt KEEPs Cached and KEEP_REMOTE_BUCKETs CachedBucket
    action_array[true][CachedBucket] = DUMP;
    action_array[true][Cached] = DUMP;

    // discard
    action_array[true][Empty] = DUMP;
    action_array[true][CachedEmpty] = DUMP;
    action_array[true][Top] = ERROR;
    action_array[true][Invalid] = ERROR;
  }

};

class PrefetchOpt : public Opt{
  public: 
  PrefetchOpt() : Opt(Pref){
    for(int i = 1; i <= NUM_NODE_TYPES; i++) 
      action_array[false][i] = DUMP;
    action_array[false][Invalid] = ERROR;
    action_array[false][Top] = ERROR;

    action_array[true][Internal] = DUMP;
    action_array[true][Bucket] = DUMP;

    action_array[true][Boundary] = KEEP;
    // we KEEP NonLocal nodes, even though we don't have their
    // children - the tw object will request the children if need be
    action_array[true][NonLocal] = KEEP;  
    action_array[true][NonLocalBucket] = KEEP_REMOTE_BUCKET;
    action_array[true][Cached] = KEEP;
    action_array[true][CachedBucket] = KEEP_REMOTE_BUCKET;

    action_array[true][CachedEmpty] = DUMP;
    action_array[true][Empty] = DUMP;
    action_array[true][Invalid] = ERROR; // node not properly initialized
    action_array[true][Top] = ERROR;    // not handled anywhere in the original code, so ERROR
  }
};
