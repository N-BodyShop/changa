#ifndef _DECL_NChilada_H_
#define _DECL_NChilada_H_
#include "charm++.h"
/* DECLS: readonly CkArrayID treePieceID;
 */

/* DECLS: readonly CkGroupID dataManagerID;
 */

/* DECLS: readonly int verbosity;
 */


/* DECLS: chare Sorter: Chare{
Sorter(void);
void startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb);
void collectEvaluations(CkReductionMsg* impl_msg);
};
 */
 class Sorter;
/* --------------- index object ------------------ */
class CkIndex_Sorter{
  public:
    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: Sorter(void);
 */
    static int __idx_Sorter_void;
    static int ckNew(void) { return __idx_Sorter_void; }
    static void _call_Sorter_void(void* impl_msg,Sorter* impl_obj);

/* DECLS: void startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb);
 */
    static int __idx_startSorting_marshall2;
    static int startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb) { return __idx_startSorting_marshall2; }
    static void _call_startSorting_marshall2(void* impl_msg,Sorter* impl_obj);

/* DECLS: void collectEvaluations(CkReductionMsg* impl_msg);
 */
    static int __idx_collectEvaluations_CkReductionMsg;
    static int collectEvaluations(CkReductionMsg* impl_msg) { return __idx_collectEvaluations_CkReductionMsg; }
    static void _call_collectEvaluations_CkReductionMsg(void* impl_msg,Sorter* impl_obj);
};
/* --------------- element proxy ------------------ */
class CProxy_Sorter:public CProxy_Chare{
  public:
    CProxy_Sorter(void) {};
    CProxy_Sorter(CkChareID __cid) : CProxy_Chare(__cid){  }
    CProxy_Sorter(const Chare *c) : CProxy_Chare(c){  }
    CK_DISAMBIG_CHARE(CProxy_Chare)
    void ckDelegate(CkGroupID to) {
      CProxy_Chare::ckDelegate(to);
    }
    void ckUndelegate(void) {
      CProxy_Chare::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxy_Chare::pup(p);
    }
    void ckSetChareID(const CkChareID &c) {
      CProxy_Chare::ckSetChareID(c);
    }
    Sorter *ckLocal(void) const
     { return (Sorter *)CkLocalChare(&ckGetChareID()); }
/* DECLS: Sorter(void);
 */
    static CkChareID ckNew(int onPE=CK_PE_ANY);
    static void ckNew(CkChareID* pcid, int onPE=CK_PE_ANY);

/* DECLS: void startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb);
 */
    void startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb);

/* DECLS: void collectEvaluations(CkReductionMsg* impl_msg);
 */
    void collectEvaluations(CkReductionMsg* impl_msg);
};
PUPmarshall(CProxy_Sorter);
typedef CBaseT<Chare,CProxy_Sorter>  CBase_Sorter;

/* DECLS: nodegroup DataManager: NodeGroup{
DataManager(const CkArrayID &treePieceID);
void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb);
void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb);
void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb);
};
 */
 class DataManager;
/* --------------- index object ------------------ */
class CkIndex_DataManager{
  public:
    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: DataManager(const CkArrayID &treePieceID);
 */
    static int __idx_DataManager_marshall1;
    static int ckNew(const CkArrayID &treePieceID) { return __idx_DataManager_marshall1; }
    static void _call_DataManager_marshall1(void* impl_msg,DataManager* impl_obj);

/* DECLS: void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb);
 */
    static int __idx_loadParticles_marshall2;
    static int loadParticles(const std::string &filename, int nbodies, const CkCallback &cb) { return __idx_loadParticles_marshall2; }
    static void _call_loadParticles_marshall2(void* impl_msg,DataManager* impl_obj);

/* DECLS: void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb);
 */
    static int __idx_acceptCandidateKeys_marshall3;
    static int acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb) { return __idx_acceptCandidateKeys_marshall3; }
    static void _call_acceptCandidateKeys_marshall3(void* impl_msg,DataManager* impl_obj);

/* DECLS: void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb);
 */
    static int __idx_acceptFinalKeys_marshall4;
    static int acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb) { return __idx_acceptFinalKeys_marshall4; }
    static void _call_acceptFinalKeys_marshall4(void* impl_msg,DataManager* impl_obj);
};
/* --------------- element proxy ------------------ */
class CProxyElement_DataManager: public CProxyElement_NodeGroup{
  public:
    CProxyElement_DataManager(void) {}
    CProxyElement_DataManager(const IrrGroup *g) : CProxyElement_NodeGroup(g){  }
    CProxyElement_DataManager(CkGroupID _gid,int _onPE,CkGroupID dTo) : CProxyElement_NodeGroup(_gid,_onPE,dTo){  }
    CProxyElement_DataManager(CkGroupID _gid,int _onPE) : CProxyElement_NodeGroup(_gid,_onPE){  }
   CK_DISAMBIG_GROUP_ELEMENT(CProxyElement_NodeGroup)
    void ckDelegate(CkGroupID to) {
      CProxyElement_NodeGroup::ckDelegate(to);
    }
    void ckUndelegate(void) {
      CProxyElement_NodeGroup::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxyElement_NodeGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxyElement_NodeGroup::ckSetGroupID(g);
    }
    DataManager* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static DataManager* ckLocalBranch(CkGroupID gID) {
      return (DataManager*)CkLocalNodeBranch(gID);
    }
/* DECLS: DataManager(const CkArrayID &treePieceID);
 */

/* DECLS: void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb);
 */
    void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb)
    {
    ckCheck();
  //Marshall: const std::string &filename, int nbodies, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)filename;
    implP|nbodies;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)filename;
    implP|nbodies;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
      if (ckIsDelegated()) {
         CkNodeGroupMsgPrep(CkIndex_DataManager::__idx_loadParticles_marshall2, impl_msg, ckGetGroupID());
         ckDelegatedTo()->NodeGroupSend(CkIndex_DataManager::__idx_loadParticles_marshall2, impl_msg, ckGetGroupPe(), ckGetGroupID());
      } else CkSendMsgNodeBranch(CkIndex_DataManager::__idx_loadParticles_marshall2, impl_msg, ckGetGroupPe(), ckGetGroupID());
    }

/* DECLS: void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb);
 */
    void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb)
    {
    ckCheck();
  //Marshall: const Key *keys, int n, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  int impl_off_keys, impl_cnt_keys;
  impl_off_keys=impl_off=CK_ALIGN(impl_off,sizeof(Key));
  impl_off+=(impl_cnt_keys=sizeof(Key)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_keys;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_keys;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_keys,keys,impl_cnt_keys);
      if (ckIsDelegated()) {
         CkNodeGroupMsgPrep(CkIndex_DataManager::__idx_acceptCandidateKeys_marshall3, impl_msg, ckGetGroupID());
         ckDelegatedTo()->NodeGroupSend(CkIndex_DataManager::__idx_acceptCandidateKeys_marshall3, impl_msg, ckGetGroupPe(), ckGetGroupID());
      } else CkSendMsgNodeBranch(CkIndex_DataManager::__idx_acceptCandidateKeys_marshall3, impl_msg, ckGetGroupPe(), ckGetGroupID());
    }

/* DECLS: void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb);
 */
    void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb)
    {
    ckCheck();
  //Marshall: const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  int impl_off_keys, impl_cnt_keys;
  impl_off_keys=impl_off=CK_ALIGN(impl_off,sizeof(Key));
  impl_off+=(impl_cnt_keys=sizeof(Key)*(n));
  int impl_off_responsible, impl_cnt_responsible;
  impl_off_responsible=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_responsible=sizeof(int)*(n - 1));
  int impl_off_bins, impl_cnt_bins;
  impl_off_bins=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_bins=sizeof(int)*(n - 1));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_keys;
    implP|impl_off_responsible;
    implP|impl_off_bins;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_keys;
    implP|impl_off_responsible;
    implP|impl_off_bins;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_keys,keys,impl_cnt_keys);
  memcpy(impl_buf+impl_off_responsible,responsible,impl_cnt_responsible);
  memcpy(impl_buf+impl_off_bins,bins,impl_cnt_bins);
      if (ckIsDelegated()) {
         CkNodeGroupMsgPrep(CkIndex_DataManager::__idx_acceptFinalKeys_marshall4, impl_msg, ckGetGroupID());
         ckDelegatedTo()->NodeGroupSend(CkIndex_DataManager::__idx_acceptFinalKeys_marshall4, impl_msg, ckGetGroupPe(), ckGetGroupID());
      } else CkSendMsgNodeBranch(CkIndex_DataManager::__idx_acceptFinalKeys_marshall4, impl_msg, ckGetGroupPe(), ckGetGroupID());
    }
};
PUPmarshall(CProxyElement_DataManager);
/* ---------------- collective proxy -------------- */
class CProxy_DataManager: public CProxy_NodeGroup{
  public:
    CProxy_DataManager(void) {}
    CProxy_DataManager(const IrrGroup *g) : CProxy_NodeGroup(g){  }
    CProxy_DataManager(CkGroupID _gid,CkGroupID dTo) : CProxy_NodeGroup(_gid,dTo){  }
    CProxy_DataManager(CkGroupID _gid) : CProxy_NodeGroup(_gid){  }
    CProxyElement_DataManager operator[](int onPE) const
      {return CProxyElement_DataManager(ckGetGroupID(),onPE,ckDelegatedIdx());}
   CK_DISAMBIG_GROUP(CProxy_NodeGroup)
    void ckDelegate(CkGroupID to) {
      CProxy_NodeGroup::ckDelegate(to);
    }
    void ckUndelegate(void) {
      CProxy_NodeGroup::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxy_NodeGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxy_NodeGroup::ckSetGroupID(g);
    }
    DataManager* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static DataManager* ckLocalBranch(CkGroupID gID) {
      return (DataManager*)CkLocalNodeBranch(gID);
    }
/* DECLS: DataManager(const CkArrayID &treePieceID);
 */
    static CkGroupID ckNew(const CkArrayID &treePieceID);
    CProxy_DataManager(const CkArrayID &treePieceID);

/* DECLS: void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb);
 */
    void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb)
    {
    ckCheck();
  //Marshall: const std::string &filename, int nbodies, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)filename;
    implP|nbodies;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)filename;
    implP|nbodies;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
      if (ckIsDelegated()) {
         CkNodeGroupMsgPrep(CkIndex_DataManager::__idx_loadParticles_marshall2, impl_msg, ckGetGroupID());
         ckDelegatedTo()->NodeGroupBroadcast(CkIndex_DataManager::__idx_loadParticles_marshall2, impl_msg, ckGetGroupID());
      } else CkBroadcastMsgNodeBranch(CkIndex_DataManager::__idx_loadParticles_marshall2, impl_msg, ckGetGroupID());
    }

/* DECLS: void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb);
 */
    void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb)
    {
    ckCheck();
  //Marshall: const Key *keys, int n, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  int impl_off_keys, impl_cnt_keys;
  impl_off_keys=impl_off=CK_ALIGN(impl_off,sizeof(Key));
  impl_off+=(impl_cnt_keys=sizeof(Key)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_keys;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_keys;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_keys,keys,impl_cnt_keys);
      if (ckIsDelegated()) {
         CkNodeGroupMsgPrep(CkIndex_DataManager::__idx_acceptCandidateKeys_marshall3, impl_msg, ckGetGroupID());
         ckDelegatedTo()->NodeGroupBroadcast(CkIndex_DataManager::__idx_acceptCandidateKeys_marshall3, impl_msg, ckGetGroupID());
      } else CkBroadcastMsgNodeBranch(CkIndex_DataManager::__idx_acceptCandidateKeys_marshall3, impl_msg, ckGetGroupID());
    }

/* DECLS: void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb);
 */
    void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb)
    {
    ckCheck();
  //Marshall: const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  int impl_off_keys, impl_cnt_keys;
  impl_off_keys=impl_off=CK_ALIGN(impl_off,sizeof(Key));
  impl_off+=(impl_cnt_keys=sizeof(Key)*(n));
  int impl_off_responsible, impl_cnt_responsible;
  impl_off_responsible=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_responsible=sizeof(int)*(n - 1));
  int impl_off_bins, impl_cnt_bins;
  impl_off_bins=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_bins=sizeof(int)*(n - 1));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_keys;
    implP|impl_off_responsible;
    implP|impl_off_bins;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_keys;
    implP|impl_off_responsible;
    implP|impl_off_bins;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_keys,keys,impl_cnt_keys);
  memcpy(impl_buf+impl_off_responsible,responsible,impl_cnt_responsible);
  memcpy(impl_buf+impl_off_bins,bins,impl_cnt_bins);
      if (ckIsDelegated()) {
         CkNodeGroupMsgPrep(CkIndex_DataManager::__idx_acceptFinalKeys_marshall4, impl_msg, ckGetGroupID());
         ckDelegatedTo()->NodeGroupBroadcast(CkIndex_DataManager::__idx_acceptFinalKeys_marshall4, impl_msg, ckGetGroupID());
      } else CkBroadcastMsgNodeBranch(CkIndex_DataManager::__idx_acceptFinalKeys_marshall4, impl_msg, ckGetGroupID());
    }
};
PUPmarshall(CProxy_DataManager);
typedef CBaseT<NodeGroup,CProxy_DataManager>  CBase_DataManager;

/* DECLS: array TreePiece: ArrayElement{
TreePiece(CkMigrateMessage* impl_msg);
TreePiece(void);
void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb);
void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb);
void assignKeys(CkReductionMsg* impl_msg);
void evaluateBoundaries(const CkCallback &cb);
void unshuffleParticles(const CkCallback &cb);
void acceptSortedParticles(const FullParticle *particles, int n);
void shareBoundaries(CkReductionMsg* impl_msg);
void acceptBoundaryKey(const Key &k);
void startTreeBuild(const CkCallback &cb);
void report(void);
};
 */
 class TreePiece;
/* --------------- index object ------------------ */
class CkIndex_TreePiece{
  public:
    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: TreePiece(CkMigrateMessage* impl_msg);
 */
    static int __idx_TreePiece_CkMigrateMessage;
    static int ckNew(CkMigrateMessage* impl_msg) { return __idx_TreePiece_CkMigrateMessage; }
    static void _call_TreePiece_CkMigrateMessage(void* impl_msg,TreePiece* impl_obj);

/* DECLS: TreePiece(void);
 */
    static int __idx_TreePiece_void;
    static int ckNew(void) { return __idx_TreePiece_void; }
    static void _call_TreePiece_void(void* impl_msg,TreePiece* impl_obj);

/* DECLS: void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb);
 */
    static int __idx_registerWithDataManager_marshall2;
    static int registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb) { return __idx_registerWithDataManager_marshall2; }
    static void _call_registerWithDataManager_marshall2(void* impl_msg,TreePiece* impl_obj);

/* DECLS: void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb);
 */
    static int __idx_receiveParticles_marshall3;
    static int receiveParticles(const FullParticle *particles, int n, const CkCallback &cb) { return __idx_receiveParticles_marshall3; }
    static void _call_receiveParticles_marshall3(void* impl_msg,TreePiece* impl_obj);

/* DECLS: void assignKeys(CkReductionMsg* impl_msg);
 */
    static int __idx_assignKeys_CkReductionMsg;
    static int assignKeys(CkReductionMsg* impl_msg) { return __idx_assignKeys_CkReductionMsg; }
    static void _call_assignKeys_CkReductionMsg(void* impl_msg,TreePiece* impl_obj);

/* DECLS: void evaluateBoundaries(const CkCallback &cb);
 */
    static int __idx_evaluateBoundaries_marshall5;
    static int evaluateBoundaries(const CkCallback &cb) { return __idx_evaluateBoundaries_marshall5; }
    static void _call_evaluateBoundaries_marshall5(void* impl_msg,TreePiece* impl_obj);

/* DECLS: void unshuffleParticles(const CkCallback &cb);
 */
    static int __idx_unshuffleParticles_marshall6;
    static int unshuffleParticles(const CkCallback &cb) { return __idx_unshuffleParticles_marshall6; }
    static void _call_unshuffleParticles_marshall6(void* impl_msg,TreePiece* impl_obj);

/* DECLS: void acceptSortedParticles(const FullParticle *particles, int n);
 */
    static int __idx_acceptSortedParticles_marshall7;
    static int acceptSortedParticles(const FullParticle *particles, int n) { return __idx_acceptSortedParticles_marshall7; }
    static void _call_acceptSortedParticles_marshall7(void* impl_msg,TreePiece* impl_obj);

/* DECLS: void shareBoundaries(CkReductionMsg* impl_msg);
 */
    static int __idx_shareBoundaries_CkReductionMsg;
    static int shareBoundaries(CkReductionMsg* impl_msg) { return __idx_shareBoundaries_CkReductionMsg; }
    static void _call_shareBoundaries_CkReductionMsg(void* impl_msg,TreePiece* impl_obj);

/* DECLS: void acceptBoundaryKey(const Key &k);
 */
    static int __idx_acceptBoundaryKey_marshall9;
    static int acceptBoundaryKey(const Key &k) { return __idx_acceptBoundaryKey_marshall9; }
    static void _call_acceptBoundaryKey_marshall9(void* impl_msg,TreePiece* impl_obj);

/* DECLS: void startTreeBuild(const CkCallback &cb);
 */
    static int __idx_startTreeBuild_marshall10;
    static int startTreeBuild(const CkCallback &cb) { return __idx_startTreeBuild_marshall10; }
    static void _call_startTreeBuild_marshall10(void* impl_msg,TreePiece* impl_obj);

/* DECLS: void report(void);
 */
    static int __idx_report_void;
    static int report(void) { return __idx_report_void; }
    static void _call_report_void(void* impl_msg,TreePiece* impl_obj);
};
/* --------------- element proxy ------------------ */
 class CProxyElement_TreePiece : public CProxyElement_ArrayElement{
  public:
    CProxyElement_TreePiece(void) {}
    CProxyElement_TreePiece(const ArrayElement *e) : CProxyElement_ArrayElement(e){  }
    void ckDelegate(CkGroupID to) {
      CProxyElement_ArrayElement::ckDelegate(to);
    }
    void ckUndelegate(void) {
      CProxyElement_ArrayElement::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxyElement_ArrayElement::pup(p);
    }
    CK_DISAMBIG_ARRAY_ELEMENT(CProxyElement_ArrayElement)
    TreePiece *ckLocal(void) const
      { return (TreePiece *)CProxyElement_ArrayElement::ckLocal(); }
    CProxyElement_TreePiece(const CkArrayID &aid,const CkArrayIndex1D &idx,CkGroupID dTo)
        :CProxyElement_ArrayElement(aid,idx,dTo) {}
    CProxyElement_TreePiece(const CkArrayID &aid,const CkArrayIndex1D &idx)
        :CProxyElement_ArrayElement(aid,idx) {}
/* DECLS: TreePiece(CkMigrateMessage* impl_msg);
 */

/* DECLS: TreePiece(void);
 */
    void insert(int onPE=-1);
/* DECLS: void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb);
 */
    void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb) ;

/* DECLS: void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb);
 */
    void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb) ;

/* DECLS: void assignKeys(CkReductionMsg* impl_msg);
 */
    void assignKeys(CkReductionMsg* impl_msg) ;

/* DECLS: void evaluateBoundaries(const CkCallback &cb);
 */
    void evaluateBoundaries(const CkCallback &cb) ;

/* DECLS: void unshuffleParticles(const CkCallback &cb);
 */
    void unshuffleParticles(const CkCallback &cb) ;

/* DECLS: void acceptSortedParticles(const FullParticle *particles, int n);
 */
    void acceptSortedParticles(const FullParticle *particles, int n) ;

/* DECLS: void shareBoundaries(CkReductionMsg* impl_msg);
 */
    void shareBoundaries(CkReductionMsg* impl_msg) ;

/* DECLS: void acceptBoundaryKey(const Key &k);
 */
    void acceptBoundaryKey(const Key &k) ;

/* DECLS: void startTreeBuild(const CkCallback &cb);
 */
    void startTreeBuild(const CkCallback &cb) ;

/* DECLS: void report(void);
 */
    void report(void) ;
};
PUPmarshall(CProxyElement_TreePiece);
/* ---------------- collective proxy -------------- */
 class CProxy_TreePiece : public CProxy_ArrayElement{
  public:
    CProxy_TreePiece(void) {}
    CProxy_TreePiece(const ArrayElement *e) : CProxy_ArrayElement(e){  }
    void ckDelegate(CkGroupID to) {
      CProxy_ArrayElement::ckDelegate(to);
    }
    void ckUndelegate(void) {
      CProxy_ArrayElement::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxy_ArrayElement::pup(p);
    }
    CK_DISAMBIG_ARRAY(CProxy_ArrayElement)
    static CkArrayID ckNew(void) {return ckCreateEmptyArray();}
//Generalized array indexing:
    CProxyElement_TreePiece operator [] (const CkArrayIndex1D &idx) const
        {return CProxyElement_TreePiece(ckGetArrayID(), idx, ckDelegatedIdx());}
    CProxyElement_TreePiece operator() (const CkArrayIndex1D &idx) const
        {return CProxyElement_TreePiece(ckGetArrayID(), idx, ckDelegatedIdx());}
    CProxyElement_TreePiece operator [] (int idx) const 
        {return CProxyElement_TreePiece(ckGetArrayID(), CkArrayIndex1D(idx), ckDelegatedIdx());}
    CProxyElement_TreePiece operator () (int idx) const 
        {return CProxyElement_TreePiece(ckGetArrayID(), CkArrayIndex1D(idx), ckDelegatedIdx());}
    CProxy_TreePiece(const CkArrayID &aid,CkGroupID dTo) 
        :CProxy_ArrayElement(aid,dTo) {}
    CProxy_TreePiece(const CkArrayID &aid) 
        :CProxy_ArrayElement(aid) {}
/* DECLS: TreePiece(CkMigrateMessage* impl_msg);
 */

/* DECLS: TreePiece(void);
 */
    static CkArrayID ckNew(const CkArrayOptions &opts);

/* DECLS: void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb);
 */
    void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb) ;

/* DECLS: void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb);
 */
    void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb) ;

/* DECLS: void assignKeys(CkReductionMsg* impl_msg);
 */
    void assignKeys(CkReductionMsg* impl_msg) ;

/* DECLS: void evaluateBoundaries(const CkCallback &cb);
 */
    void evaluateBoundaries(const CkCallback &cb) ;

/* DECLS: void unshuffleParticles(const CkCallback &cb);
 */
    void unshuffleParticles(const CkCallback &cb) ;

/* DECLS: void acceptSortedParticles(const FullParticle *particles, int n);
 */
    void acceptSortedParticles(const FullParticle *particles, int n) ;

/* DECLS: void shareBoundaries(CkReductionMsg* impl_msg);
 */
    void shareBoundaries(CkReductionMsg* impl_msg) ;

/* DECLS: void acceptBoundaryKey(const Key &k);
 */
    void acceptBoundaryKey(const Key &k) ;

/* DECLS: void startTreeBuild(const CkCallback &cb);
 */
    void startTreeBuild(const CkCallback &cb) ;

/* DECLS: void report(void);
 */
    void report(void) ;
};
PUPmarshall(CProxy_TreePiece);
/* ---------------- section proxy -------------- */
 class CProxySection_TreePiece : public CProxySection_ArrayElement{
  public:
    CProxySection_TreePiece(void) {}
    void ckDelegate(CkGroupID to) {
      CProxySection_ArrayElement::ckDelegate(to);
    }
    void ckUndelegate(void) {
      CProxySection_ArrayElement::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxySection_ArrayElement::pup(p);
    }
    CK_DISAMBIG_ARRAY_SECTION(CProxySection_ArrayElement)
//Generalized array indexing:
    CProxyElement_TreePiece operator [] (const CkArrayIndex1D &idx) const
        {return CProxyElement_TreePiece(ckGetArrayID(), idx, ckDelegatedIdx());}
    CProxyElement_TreePiece operator() (const CkArrayIndex1D &idx) const
        {return CProxyElement_TreePiece(ckGetArrayID(), idx, ckDelegatedIdx());}
    CProxyElement_TreePiece operator [] (int idx) const 
        {return CProxyElement_TreePiece(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], ckDelegatedIdx());}
    CProxyElement_TreePiece operator () (int idx) const 
        {return CProxyElement_TreePiece(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], ckDelegatedIdx());}
    CProxySection_TreePiece(const CkArrayID &aid, CkArrayIndexMax *elems, int nElems, CkGroupID dTo) 
        :CProxySection_ArrayElement(aid,elems,nElems,dTo) {}
    CProxySection_TreePiece(const CkArrayID &aid, CkArrayIndexMax *elems, int nElems) 
        :CProxySection_ArrayElement(aid,elems,nElems) {}
    CProxySection_TreePiece(const CkSectionID &sid)       :CProxySection_ArrayElement(sid) {}
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndexMax *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
/* DECLS: TreePiece(CkMigrateMessage* impl_msg);
 */

/* DECLS: TreePiece(void);
 */

/* DECLS: void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb);
 */
    void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb) ;

/* DECLS: void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb);
 */
    void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb) ;

/* DECLS: void assignKeys(CkReductionMsg* impl_msg);
 */
    void assignKeys(CkReductionMsg* impl_msg) ;

/* DECLS: void evaluateBoundaries(const CkCallback &cb);
 */
    void evaluateBoundaries(const CkCallback &cb) ;

/* DECLS: void unshuffleParticles(const CkCallback &cb);
 */
    void unshuffleParticles(const CkCallback &cb) ;

/* DECLS: void acceptSortedParticles(const FullParticle *particles, int n);
 */
    void acceptSortedParticles(const FullParticle *particles, int n) ;

/* DECLS: void shareBoundaries(CkReductionMsg* impl_msg);
 */
    void shareBoundaries(CkReductionMsg* impl_msg) ;

/* DECLS: void acceptBoundaryKey(const Key &k);
 */
    void acceptBoundaryKey(const Key &k) ;

/* DECLS: void startTreeBuild(const CkCallback &cb);
 */
    void startTreeBuild(const CkCallback &cb) ;

/* DECLS: void report(void);
 */
    void report(void) ;
};
PUPmarshall(CProxySection_TreePiece);
typedef CBaseT<ArrayElementT<CkIndex1D>,CProxy_TreePiece>  CBase_TreePiece;

extern void _registerNChilada(void);
#endif
