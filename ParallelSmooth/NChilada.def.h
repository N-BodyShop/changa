/* DEFS: readonly CkArrayID treePieceID;
 */
extern CkArrayID treePieceID;
extern "C" void __xlater_roPup_treePieceID(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|treePieceID;
}

/* DEFS: readonly CkGroupID dataManagerID;
 */
extern CkGroupID dataManagerID;
extern "C" void __xlater_roPup_dataManagerID(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|dataManagerID;
}

/* DEFS: readonly int verbosity;
 */
extern int verbosity;
extern "C" void __xlater_roPup_verbosity(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|verbosity;
}


/* DEFS: chare Sorter: Chare{
Sorter(void);
void startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb);
void collectEvaluations(CkReductionMsg* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_Sorter::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: Sorter(void);
 */
CkChareID CProxy_Sorter::ckNew(int impl_onPE)
{
  void *impl_msg = CkAllocSysMsg();
  CkChareID impl_ret;
  CkCreateChare(CkIndex_Sorter::__idx, CkIndex_Sorter::__idx_Sorter_void, impl_msg, &impl_ret, impl_onPE);
  return impl_ret;
}
void CProxy_Sorter::ckNew(CkChareID* pcid, int impl_onPE)
{
  void *impl_msg = CkAllocSysMsg();
  CkCreateChare(CkIndex_Sorter::__idx, CkIndex_Sorter::__idx_Sorter_void, impl_msg, pcid, impl_onPE);
}
 int CkIndex_Sorter::__idx_Sorter_void=0;
void CkIndex_Sorter::_call_Sorter_void(void* impl_msg,Sorter * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  new (impl_obj) Sorter();
}

/* DEFS: void startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb);
 */
void CProxy_Sorter::startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb)
{
  ckCheck();
  //Marshall: const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkGroupID &)dataManagerID;
    implP|nChares;
    implP|toler;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkGroupID &)dataManagerID;
    implP|nChares;
    implP|toler;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Sorter::__idx_startSorting_marshall2, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(CkIndex_Sorter::__idx_startSorting_marshall2, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Sorter::__idx_startSorting_marshall2, impl_msg, &ckGetChareID());
}
 int CkIndex_Sorter::__idx_startSorting_marshall2=0;
void CkIndex_Sorter::_call_startSorting_marshall2(void* impl_msg,Sorter * impl_obj)
{
  //Unmarshall pup'd fields: const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  CkGroupID dataManagerID; implP|dataManagerID;
  int nChares; implP|nChares;
  double toler; implP|toler;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->startSorting(dataManagerID, nChares, toler, cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void collectEvaluations(CkReductionMsg* impl_msg);
 */
void CProxy_Sorter::collectEvaluations(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Sorter::__idx_collectEvaluations_CkReductionMsg, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(CkIndex_Sorter::__idx_collectEvaluations_CkReductionMsg, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Sorter::__idx_collectEvaluations_CkReductionMsg, impl_msg, &ckGetChareID());
}
 int CkIndex_Sorter::__idx_collectEvaluations_CkReductionMsg=0;
void CkIndex_Sorter::_call_collectEvaluations_CkReductionMsg(void* impl_msg,Sorter * impl_obj)
{
  impl_obj->collectEvaluations((CkReductionMsg*)impl_msg);
}
#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_Sorter::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size);
  _REGISTER_BASE(__idx, CkIndex_Chare::__idx);
// REG: Sorter(void);
  __idx_Sorter_void = CkRegisterEp("Sorter(void)",
     (CkCallFnPtr)_call_Sorter_void, 0, __idx);
  CkRegisterDefaultCtor(__idx, __idx_Sorter_void);

// REG: void startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb);
  __idx_startSorting_marshall2 = CkRegisterEp("startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb)",
     (CkCallFnPtr)_call_startSorting_marshall2, CkMarshallMsg::__idx, __idx);

// REG: void collectEvaluations(CkReductionMsg* impl_msg);
  __idx_collectEvaluations_CkReductionMsg = CkRegisterEp("collectEvaluations(CkReductionMsg* impl_msg)",
     (CkCallFnPtr)_call_collectEvaluations_CkReductionMsg, CMessage_CkReductionMsg::__idx, __idx);
}
#endif

/* DEFS: nodegroup DataManager: NodeGroup{
DataManager(const CkArrayID &treePieceID);
void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb);
void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb);
void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_DataManager::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: DataManager(const CkArrayID &treePieceID);
 */

/* DEFS: void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb);
 */

/* DEFS: void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb);
 */

/* DEFS: void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb);
 */
/* DEFS: DataManager(const CkArrayID &treePieceID);
 */
CkGroupID CProxy_DataManager::ckNew(const CkArrayID &treePieceID)
{
  //Marshall: const CkArrayID &treePieceID
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayID &)treePieceID;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayID &)treePieceID;
  }
  return CkCreateNodeGroup(CkIndex_DataManager::__idx, CkIndex_DataManager::__idx_DataManager_marshall1, impl_msg);
}
  CProxy_DataManager::CProxy_DataManager(const CkArrayID &treePieceID)
{
  //Marshall: const CkArrayID &treePieceID
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayID &)treePieceID;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayID &)treePieceID;
  }
  ckSetGroupID(CkCreateNodeGroup(CkIndex_DataManager::__idx, CkIndex_DataManager::__idx_DataManager_marshall1, impl_msg));
}
 int CkIndex_DataManager::__idx_DataManager_marshall1=0;
void CkIndex_DataManager::_call_DataManager_marshall1(void* impl_msg,DataManager * impl_obj)
{
  //Unmarshall pup'd fields: const CkArrayID &treePieceID
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  CkArrayID treePieceID; implP|treePieceID;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  new (impl_obj) DataManager(treePieceID);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb);
 */
 int CkIndex_DataManager::__idx_loadParticles_marshall2=0;
void CkIndex_DataManager::_call_loadParticles_marshall2(void* impl_msg,DataManager * impl_obj)
{
  //Unmarshall pup'd fields: const std::string &filename, int nbodies, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  std::string filename; implP|filename;
  int nbodies; implP|nbodies;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->loadParticles(filename, nbodies, cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb);
 */
 int CkIndex_DataManager::__idx_acceptCandidateKeys_marshall3=0;
void CkIndex_DataManager::_call_acceptCandidateKeys_marshall3(void* impl_msg,DataManager * impl_obj)
{
  //Unmarshall pup'd fields: const Key *keys, int n, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  int impl_off_keys; implP|impl_off_keys;
  int n; implP|n;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  Key *keys=(Key *)(impl_buf+impl_off_keys);
  impl_obj->acceptCandidateKeys(keys, n, cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb);
 */
 int CkIndex_DataManager::__idx_acceptFinalKeys_marshall4=0;
void CkIndex_DataManager::_call_acceptFinalKeys_marshall4(void* impl_msg,DataManager * impl_obj)
{
  //Unmarshall pup'd fields: const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  int impl_off_keys; implP|impl_off_keys;
  int impl_off_responsible; implP|impl_off_responsible;
  int impl_off_bins; implP|impl_off_bins;
  int n; implP|n;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  Key *keys=(Key *)(impl_buf+impl_off_keys);
  int *responsible=(int *)(impl_buf+impl_off_responsible);
  int *bins=(int *)(impl_buf+impl_off_bins);
  impl_obj->acceptFinalKeys(keys, responsible, bins, n, cb);
  delete (CkMarshallMsg *)impl_msg;
}
#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_DataManager::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size);
  _REGISTER_BASE(__idx, CkIndex_NodeGroup::__idx);
// REG: DataManager(const CkArrayID &treePieceID);
  __idx_DataManager_marshall1 = CkRegisterEp("DataManager(const CkArrayID &treePieceID)",
     (CkCallFnPtr)_call_DataManager_marshall1, CkMarshallMsg::__idx, __idx);

// REG: void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb);
  __idx_loadParticles_marshall2 = CkRegisterEp("loadParticles(const std::string &filename, int nbodies, const CkCallback &cb)",
     (CkCallFnPtr)_call_loadParticles_marshall2, CkMarshallMsg::__idx, __idx);

// REG: void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb);
  __idx_acceptCandidateKeys_marshall3 = CkRegisterEp("acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb)",
     (CkCallFnPtr)_call_acceptCandidateKeys_marshall3, CkMarshallMsg::__idx, __idx);

// REG: void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb);
  __idx_acceptFinalKeys_marshall4 = CkRegisterEp("acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb)",
     (CkCallFnPtr)_call_acceptFinalKeys_marshall4, CkMarshallMsg::__idx, __idx);
}
#endif

/* DEFS: array TreePiece: ArrayElement{
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
void report(const CkCallback &cb);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_TreePiece::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: TreePiece(CkMigrateMessage* impl_msg);
 */

/* DEFS: TreePiece(void);
 */
void CProxyElement_TreePiece::insert(int onPE)
{ 
  void *impl_msg = CkAllocSysMsg();
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_TreePiece::__idx_TreePiece_void,onPE);
}

/* DEFS: void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb);
 */
void CProxyElement_TreePiece::registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkGroupID &dataManagerID, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkGroupID &)dataManagerID;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkGroupID &)dataManagerID;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_registerWithDataManager_marshall2);
}

/* DEFS: void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb);
 */
void CProxyElement_TreePiece::receiveParticles(const FullParticle *particles, int n, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const FullParticle *particles, int n, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  int impl_off_particles, impl_cnt_particles;
  impl_off_particles=impl_off=CK_ALIGN(impl_off,sizeof(FullParticle));
  impl_off+=(impl_cnt_particles=sizeof(FullParticle)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_particles;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_particles;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_particles,particles,impl_cnt_particles);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_receiveParticles_marshall3);
}

/* DEFS: void assignKeys(CkReductionMsg* impl_msg);
 */
void CProxyElement_TreePiece::assignKeys(CkReductionMsg* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_assignKeys_CkReductionMsg);
}

/* DEFS: void evaluateBoundaries(const CkCallback &cb);
 */
void CProxyElement_TreePiece::evaluateBoundaries(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_evaluateBoundaries_marshall5);
}

/* DEFS: void unshuffleParticles(const CkCallback &cb);
 */
void CProxyElement_TreePiece::unshuffleParticles(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_unshuffleParticles_marshall6);
}

/* DEFS: void acceptSortedParticles(const FullParticle *particles, int n);
 */
void CProxyElement_TreePiece::acceptSortedParticles(const FullParticle *particles, int n) 
{
  ckCheck();
  //Marshall: const FullParticle *particles, int n
  int impl_off=0,impl_arrstart=0;
  int impl_off_particles, impl_cnt_particles;
  impl_off_particles=impl_off=CK_ALIGN(impl_off,sizeof(FullParticle));
  impl_off+=(impl_cnt_particles=sizeof(FullParticle)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_particles;
    implP|n;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_particles;
    implP|n;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_particles,particles,impl_cnt_particles);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_acceptSortedParticles_marshall7);
}

/* DEFS: void shareBoundaries(CkReductionMsg* impl_msg);
 */
void CProxyElement_TreePiece::shareBoundaries(CkReductionMsg* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_shareBoundaries_CkReductionMsg);
}

/* DEFS: void acceptBoundaryKey(const Key &k);
 */
void CProxyElement_TreePiece::acceptBoundaryKey(const Key &k) 
{
  ckCheck();
  //Marshall: const Key &k
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(Key &)k;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(Key &)k;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_acceptBoundaryKey_marshall9);
}

/* DEFS: void startTreeBuild(const CkCallback &cb);
 */
void CProxyElement_TreePiece::startTreeBuild(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_startTreeBuild_marshall10);
}

/* DEFS: void report(const CkCallback &cb);
 */
void CProxyElement_TreePiece::report(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_report_marshall11);
}
/* DEFS: TreePiece(CkMigrateMessage* impl_msg);
 */
 int CkIndex_TreePiece::__idx_TreePiece_CkMigrateMessage=0;
void CkIndex_TreePiece::_call_TreePiece_CkMigrateMessage(void* impl_msg,TreePiece * impl_obj)
{
  new (impl_obj) TreePiece((CkMigrateMessage*)impl_msg);
}

/* DEFS: TreePiece(void);
 */
CkArrayID CProxy_TreePiece::ckNew(const CkArrayOptions &opts)
{ 
  void *impl_msg = CkAllocSysMsg();
   return ckCreateArray((CkArrayMessage *)impl_msg,CkIndex_TreePiece::__idx_TreePiece_void,opts);
}
 int CkIndex_TreePiece::__idx_TreePiece_void=0;
void CkIndex_TreePiece::_call_TreePiece_void(void* impl_msg,TreePiece * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  new (impl_obj) TreePiece();
}

/* DEFS: void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb);
 */
void CProxy_TreePiece::registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkGroupID &dataManagerID, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkGroupID &)dataManagerID;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkGroupID &)dataManagerID;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_TreePiece::__idx_registerWithDataManager_marshall2);
}
 int CkIndex_TreePiece::__idx_registerWithDataManager_marshall2=0;
void CkIndex_TreePiece::_call_registerWithDataManager_marshall2(void* impl_msg,TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const CkGroupID &dataManagerID, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  CkGroupID dataManagerID; implP|dataManagerID;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->registerWithDataManager(dataManagerID, cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb);
 */
void CProxy_TreePiece::receiveParticles(const FullParticle *particles, int n, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const FullParticle *particles, int n, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  int impl_off_particles, impl_cnt_particles;
  impl_off_particles=impl_off=CK_ALIGN(impl_off,sizeof(FullParticle));
  impl_off+=(impl_cnt_particles=sizeof(FullParticle)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_particles;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_particles;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_particles,particles,impl_cnt_particles);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_TreePiece::__idx_receiveParticles_marshall3);
}
 int CkIndex_TreePiece::__idx_receiveParticles_marshall3=0;
void CkIndex_TreePiece::_call_receiveParticles_marshall3(void* impl_msg,TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const FullParticle *particles, int n, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  int impl_off_particles; implP|impl_off_particles;
  int n; implP|n;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  FullParticle *particles=(FullParticle *)(impl_buf+impl_off_particles);
  impl_obj->receiveParticles(particles, n, cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void assignKeys(CkReductionMsg* impl_msg);
 */
void CProxy_TreePiece::assignKeys(CkReductionMsg* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_TreePiece::__idx_assignKeys_CkReductionMsg);
}
 int CkIndex_TreePiece::__idx_assignKeys_CkReductionMsg=0;
void CkIndex_TreePiece::_call_assignKeys_CkReductionMsg(void* impl_msg,TreePiece * impl_obj)
{
  impl_obj->assignKeys((CkReductionMsg*)impl_msg);
}

/* DEFS: void evaluateBoundaries(const CkCallback &cb);
 */
void CProxy_TreePiece::evaluateBoundaries(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_TreePiece::__idx_evaluateBoundaries_marshall5);
}
 int CkIndex_TreePiece::__idx_evaluateBoundaries_marshall5=0;
void CkIndex_TreePiece::_call_evaluateBoundaries_marshall5(void* impl_msg,TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->evaluateBoundaries(cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void unshuffleParticles(const CkCallback &cb);
 */
void CProxy_TreePiece::unshuffleParticles(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_TreePiece::__idx_unshuffleParticles_marshall6);
}
 int CkIndex_TreePiece::__idx_unshuffleParticles_marshall6=0;
void CkIndex_TreePiece::_call_unshuffleParticles_marshall6(void* impl_msg,TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->unshuffleParticles(cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void acceptSortedParticles(const FullParticle *particles, int n);
 */
void CProxy_TreePiece::acceptSortedParticles(const FullParticle *particles, int n) 
{
  ckCheck();
  //Marshall: const FullParticle *particles, int n
  int impl_off=0,impl_arrstart=0;
  int impl_off_particles, impl_cnt_particles;
  impl_off_particles=impl_off=CK_ALIGN(impl_off,sizeof(FullParticle));
  impl_off+=(impl_cnt_particles=sizeof(FullParticle)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_particles;
    implP|n;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_particles;
    implP|n;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_particles,particles,impl_cnt_particles);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_TreePiece::__idx_acceptSortedParticles_marshall7);
}
 int CkIndex_TreePiece::__idx_acceptSortedParticles_marshall7=0;
void CkIndex_TreePiece::_call_acceptSortedParticles_marshall7(void* impl_msg,TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const FullParticle *particles, int n
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  int impl_off_particles; implP|impl_off_particles;
  int n; implP|n;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  FullParticle *particles=(FullParticle *)(impl_buf+impl_off_particles);
  impl_obj->acceptSortedParticles(particles, n);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void shareBoundaries(CkReductionMsg* impl_msg);
 */
void CProxy_TreePiece::shareBoundaries(CkReductionMsg* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_TreePiece::__idx_shareBoundaries_CkReductionMsg);
}
 int CkIndex_TreePiece::__idx_shareBoundaries_CkReductionMsg=0;
void CkIndex_TreePiece::_call_shareBoundaries_CkReductionMsg(void* impl_msg,TreePiece * impl_obj)
{
  impl_obj->shareBoundaries((CkReductionMsg*)impl_msg);
}

/* DEFS: void acceptBoundaryKey(const Key &k);
 */
void CProxy_TreePiece::acceptBoundaryKey(const Key &k) 
{
  ckCheck();
  //Marshall: const Key &k
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(Key &)k;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(Key &)k;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_TreePiece::__idx_acceptBoundaryKey_marshall9);
}
 int CkIndex_TreePiece::__idx_acceptBoundaryKey_marshall9=0;
void CkIndex_TreePiece::_call_acceptBoundaryKey_marshall9(void* impl_msg,TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const Key &k
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  Key k; implP|k;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->acceptBoundaryKey(k);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void startTreeBuild(const CkCallback &cb);
 */
void CProxy_TreePiece::startTreeBuild(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_TreePiece::__idx_startTreeBuild_marshall10);
}
 int CkIndex_TreePiece::__idx_startTreeBuild_marshall10=0;
void CkIndex_TreePiece::_call_startTreeBuild_marshall10(void* impl_msg,TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->startTreeBuild(cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void report(const CkCallback &cb);
 */
void CProxy_TreePiece::report(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_TreePiece::__idx_report_marshall11);
}
 int CkIndex_TreePiece::__idx_report_marshall11=0;
void CkIndex_TreePiece::_call_report_marshall11(void* impl_msg,TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->report(cb);
  delete (CkMarshallMsg *)impl_msg;
}
/* DEFS: TreePiece(CkMigrateMessage* impl_msg);
 */

/* DEFS: TreePiece(void);
 */

/* DEFS: void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb);
 */
void CProxySection_TreePiece::registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkGroupID &dataManagerID, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkGroupID &)dataManagerID;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkGroupID &)dataManagerID;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_registerWithDataManager_marshall2);
}

/* DEFS: void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb);
 */
void CProxySection_TreePiece::receiveParticles(const FullParticle *particles, int n, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const FullParticle *particles, int n, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  int impl_off_particles, impl_cnt_particles;
  impl_off_particles=impl_off=CK_ALIGN(impl_off,sizeof(FullParticle));
  impl_off+=(impl_cnt_particles=sizeof(FullParticle)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_particles;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_particles;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_particles,particles,impl_cnt_particles);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_receiveParticles_marshall3);
}

/* DEFS: void assignKeys(CkReductionMsg* impl_msg);
 */
void CProxySection_TreePiece::assignKeys(CkReductionMsg* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_assignKeys_CkReductionMsg);
}

/* DEFS: void evaluateBoundaries(const CkCallback &cb);
 */
void CProxySection_TreePiece::evaluateBoundaries(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_evaluateBoundaries_marshall5);
}

/* DEFS: void unshuffleParticles(const CkCallback &cb);
 */
void CProxySection_TreePiece::unshuffleParticles(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_unshuffleParticles_marshall6);
}

/* DEFS: void acceptSortedParticles(const FullParticle *particles, int n);
 */
void CProxySection_TreePiece::acceptSortedParticles(const FullParticle *particles, int n) 
{
  ckCheck();
  //Marshall: const FullParticle *particles, int n
  int impl_off=0,impl_arrstart=0;
  int impl_off_particles, impl_cnt_particles;
  impl_off_particles=impl_off=CK_ALIGN(impl_off,sizeof(FullParticle));
  impl_off+=(impl_cnt_particles=sizeof(FullParticle)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_particles;
    implP|n;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_particles;
    implP|n;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_particles,particles,impl_cnt_particles);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_acceptSortedParticles_marshall7);
}

/* DEFS: void shareBoundaries(CkReductionMsg* impl_msg);
 */
void CProxySection_TreePiece::shareBoundaries(CkReductionMsg* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_shareBoundaries_CkReductionMsg);
}

/* DEFS: void acceptBoundaryKey(const Key &k);
 */
void CProxySection_TreePiece::acceptBoundaryKey(const Key &k) 
{
  ckCheck();
  //Marshall: const Key &k
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(Key &)k;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(Key &)k;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_acceptBoundaryKey_marshall9);
}

/* DEFS: void startTreeBuild(const CkCallback &cb);
 */
void CProxySection_TreePiece::startTreeBuild(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_startTreeBuild_marshall10);
}

/* DEFS: void report(const CkCallback &cb);
 */
void CProxySection_TreePiece::report(const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_TreePiece::__idx_report_marshall11);
}
#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_TreePiece::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size);
  _REGISTER_BASE(__idx, CkIndex_ArrayElement::__idx);
// REG: TreePiece(CkMigrateMessage* impl_msg);
  __idx_TreePiece_CkMigrateMessage = CkRegisterEp("TreePiece(CkMigrateMessage* impl_msg)",
     (CkCallFnPtr)_call_TreePiece_CkMigrateMessage, 0, __idx);
  CkRegisterMigCtor(__idx, __idx_TreePiece_CkMigrateMessage);

// REG: TreePiece(void);
  __idx_TreePiece_void = CkRegisterEp("TreePiece(void)",
     (CkCallFnPtr)_call_TreePiece_void, 0, __idx);
  CkRegisterDefaultCtor(__idx, __idx_TreePiece_void);

// REG: void registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb);
  __idx_registerWithDataManager_marshall2 = CkRegisterEp("registerWithDataManager(const CkGroupID &dataManagerID, const CkCallback &cb)",
     (CkCallFnPtr)_call_registerWithDataManager_marshall2, CkMarshallMsg::__idx, __idx);

// REG: void receiveParticles(const FullParticle *particles, int n, const CkCallback &cb);
  __idx_receiveParticles_marshall3 = CkRegisterEp("receiveParticles(const FullParticle *particles, int n, const CkCallback &cb)",
     (CkCallFnPtr)_call_receiveParticles_marshall3, CkMarshallMsg::__idx, __idx);

// REG: void assignKeys(CkReductionMsg* impl_msg);
  __idx_assignKeys_CkReductionMsg = CkRegisterEp("assignKeys(CkReductionMsg* impl_msg)",
     (CkCallFnPtr)_call_assignKeys_CkReductionMsg, CMessage_CkReductionMsg::__idx, __idx);

// REG: void evaluateBoundaries(const CkCallback &cb);
  __idx_evaluateBoundaries_marshall5 = CkRegisterEp("evaluateBoundaries(const CkCallback &cb)",
     (CkCallFnPtr)_call_evaluateBoundaries_marshall5, CkMarshallMsg::__idx, __idx);

// REG: void unshuffleParticles(const CkCallback &cb);
  __idx_unshuffleParticles_marshall6 = CkRegisterEp("unshuffleParticles(const CkCallback &cb)",
     (CkCallFnPtr)_call_unshuffleParticles_marshall6, CkMarshallMsg::__idx, __idx);

// REG: void acceptSortedParticles(const FullParticle *particles, int n);
  __idx_acceptSortedParticles_marshall7 = CkRegisterEp("acceptSortedParticles(const FullParticle *particles, int n)",
     (CkCallFnPtr)_call_acceptSortedParticles_marshall7, CkMarshallMsg::__idx, __idx);

// REG: void shareBoundaries(CkReductionMsg* impl_msg);
  __idx_shareBoundaries_CkReductionMsg = CkRegisterEp("shareBoundaries(CkReductionMsg* impl_msg)",
     (CkCallFnPtr)_call_shareBoundaries_CkReductionMsg, CMessage_CkReductionMsg::__idx, __idx);

// REG: void acceptBoundaryKey(const Key &k);
  __idx_acceptBoundaryKey_marshall9 = CkRegisterEp("acceptBoundaryKey(const Key &k)",
     (CkCallFnPtr)_call_acceptBoundaryKey_marshall9, CkMarshallMsg::__idx, __idx);

// REG: void startTreeBuild(const CkCallback &cb);
  __idx_startTreeBuild_marshall10 = CkRegisterEp("startTreeBuild(const CkCallback &cb)",
     (CkCallFnPtr)_call_startTreeBuild_marshall10, CkMarshallMsg::__idx, __idx);

// REG: void report(const CkCallback &cb);
  __idx_report_marshall11 = CkRegisterEp("report(const CkCallback &cb)",
     (CkCallFnPtr)_call_report_marshall11, CkMarshallMsg::__idx, __idx);
}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registerNChilada(void)
{
  static int _done = 0; if(_done) return; _done = 1;
  CkRegisterReadonly("treePieceID","CkArrayID",sizeof(treePieceID),(void *) &treePieceID,__xlater_roPup_treePieceID);

  CkRegisterReadonly("dataManagerID","CkGroupID",sizeof(dataManagerID),(void *) &dataManagerID,__xlater_roPup_dataManagerID);

  CkRegisterReadonly("verbosity","int",sizeof(verbosity),(void *) &verbosity,__xlater_roPup_verbosity);

      registerNChiladaReductions();

/* REG: chare Sorter: Chare{
Sorter(void);
void startSorting(const CkGroupID &dataManagerID, int nChares, double toler, const CkCallback &cb);
void collectEvaluations(CkReductionMsg* impl_msg);
};
*/
  CkIndex_Sorter::__register("Sorter", sizeof(Sorter));

/* REG: nodegroup DataManager: NodeGroup{
DataManager(const CkArrayID &treePieceID);
void loadParticles(const std::string &filename, int nbodies, const CkCallback &cb);
void acceptCandidateKeys(const Key *keys, int n, const CkCallback &cb);
void acceptFinalKeys(const Key *keys, const int *responsible, const int *bins, int n, const CkCallback &cb);
};
*/
  CkIndex_DataManager::__register("DataManager", sizeof(DataManager));

/* REG: array TreePiece: ArrayElement{
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
void report(const CkCallback &cb);
};
*/
  CkIndex_TreePiece::__register("TreePiece", sizeof(TreePiece));

}
#endif
