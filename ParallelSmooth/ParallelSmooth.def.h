
/* DEFS: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
threaded void loadSortBuild(void);
threaded void doCalculations(void);
threaded void makeHistogram(CkReductionMsg* impl_msg);
threaded void doStressTest(CkReductionMsg* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_Main::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: Main(CkArgMsg* impl_msg);
 */
CkChareID CProxy_Main::ckNew(CkArgMsg* impl_msg, int impl_onPE)
{
  CkChareID impl_ret;
  CkCreateChare(CkIndex_Main::__idx, CkIndex_Main::__idx_Main_CkArgMsg, impl_msg, &impl_ret, impl_onPE);
  return impl_ret;
}
void CProxy_Main::ckNew(CkArgMsg* impl_msg, CkChareID* pcid, int impl_onPE)
{
  CkCreateChare(CkIndex_Main::__idx, CkIndex_Main::__idx_Main_CkArgMsg, impl_msg, pcid, impl_onPE);
}
  CProxy_Main::CProxy_Main(CkArgMsg* impl_msg, int impl_onPE)
{
  CkChareID impl_ret;
  CkCreateChare(CkIndex_Main::__idx, CkIndex_Main::__idx_Main_CkArgMsg, impl_msg, &impl_ret, impl_onPE);
  ckSetChareID(impl_ret);
}
 int CkIndex_Main::__idx_Main_CkArgMsg=0;
void CkIndex_Main::_call_Main_CkArgMsg(void* impl_msg,Main * impl_obj)
{
  new (impl_obj) Main((CkArgMsg*)impl_msg);
}

/* DEFS: threaded void loadSortBuild(void);
 */
void CProxy_Main::loadSortBuild(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::__idx_loadSortBuild_void, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(CkIndex_Main::__idx_loadSortBuild_void, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Main::__idx_loadSortBuild_void, impl_msg, &ckGetChareID());
}
 int CkIndex_Main::__idx_loadSortBuild_void=0;
void CkIndex_Main::_call_loadSortBuild_void(void* impl_msg,Main * impl_obj)
{
  CthAwaken(CthCreate((CthVoidFn)_callthr_loadSortBuild_void, new CkThrCallArg(impl_msg,impl_obj), 0));
}
void CkIndex_Main::_callthr_loadSortBuild_void(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  Main *impl_obj = (Main *) impl_arg->obj;
  delete impl_arg;
  CkFreeSysMsg(impl_msg);
  impl_obj->loadSortBuild();
}

/* DEFS: threaded void doCalculations(void);
 */
void CProxy_Main::doCalculations(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::__idx_doCalculations_void, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(CkIndex_Main::__idx_doCalculations_void, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Main::__idx_doCalculations_void, impl_msg, &ckGetChareID());
}
 int CkIndex_Main::__idx_doCalculations_void=0;
void CkIndex_Main::_call_doCalculations_void(void* impl_msg,Main * impl_obj)
{
  CthAwaken(CthCreate((CthVoidFn)_callthr_doCalculations_void, new CkThrCallArg(impl_msg,impl_obj), 0));
}
void CkIndex_Main::_callthr_doCalculations_void(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  Main *impl_obj = (Main *) impl_arg->obj;
  delete impl_arg;
  CkFreeSysMsg(impl_msg);
  impl_obj->doCalculations();
}

/* DEFS: threaded void makeHistogram(CkReductionMsg* impl_msg);
 */
void CProxy_Main::makeHistogram(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::__idx_makeHistogram_CkReductionMsg, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(CkIndex_Main::__idx_makeHistogram_CkReductionMsg, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Main::__idx_makeHistogram_CkReductionMsg, impl_msg, &ckGetChareID());
}
 int CkIndex_Main::__idx_makeHistogram_CkReductionMsg=0;
void CkIndex_Main::_call_makeHistogram_CkReductionMsg(void* impl_msg,Main * impl_obj)
{
  CthAwaken(CthCreate((CthVoidFn)_callthr_makeHistogram_CkReductionMsg, new CkThrCallArg(impl_msg,impl_obj), 0));
}
void CkIndex_Main::_callthr_makeHistogram_CkReductionMsg(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  Main *impl_obj = (Main *) impl_arg->obj;
  delete impl_arg;
  impl_obj->makeHistogram((CkReductionMsg*)impl_msg);
}

/* DEFS: threaded void doStressTest(CkReductionMsg* impl_msg);
 */
void CProxy_Main::doStressTest(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::__idx_doStressTest_CkReductionMsg, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(CkIndex_Main::__idx_doStressTest_CkReductionMsg, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Main::__idx_doStressTest_CkReductionMsg, impl_msg, &ckGetChareID());
}
 int CkIndex_Main::__idx_doStressTest_CkReductionMsg=0;
void CkIndex_Main::_call_doStressTest_CkReductionMsg(void* impl_msg,Main * impl_obj)
{
  CthAwaken(CthCreate((CthVoidFn)_callthr_doStressTest_CkReductionMsg, new CkThrCallArg(impl_msg,impl_obj), 0));
}
void CkIndex_Main::_callthr_doStressTest_CkReductionMsg(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  Main *impl_obj = (Main *) impl_arg->obj;
  delete impl_arg;
  impl_obj->doStressTest((CkReductionMsg*)impl_msg);
}
#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_Main::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size);
  _REGISTER_BASE(__idx, CkIndex_Chare::__idx);
// REG: Main(CkArgMsg* impl_msg);
  __idx_Main_CkArgMsg = CkRegisterEp("Main(CkArgMsg* impl_msg)",
     (CkCallFnPtr)_call_Main_CkArgMsg, CMessage_CkArgMsg::__idx, __idx);
  CkRegisterMainChare(__idx, __idx_Main_CkArgMsg);

// REG: threaded void loadSortBuild(void);
  __idx_loadSortBuild_void = CkRegisterEp("loadSortBuild(void)",
     (CkCallFnPtr)_call_loadSortBuild_void, 0, __idx);

// REG: threaded void doCalculations(void);
  __idx_doCalculations_void = CkRegisterEp("doCalculations(void)",
     (CkCallFnPtr)_call_doCalculations_void, 0, __idx);

// REG: threaded void makeHistogram(CkReductionMsg* impl_msg);
  __idx_makeHistogram_CkReductionMsg = CkRegisterEp("makeHistogram(CkReductionMsg* impl_msg)",
     (CkCallFnPtr)_call_makeHistogram_CkReductionMsg, CMessage_CkReductionMsg::__idx, __idx);

// REG: threaded void doStressTest(CkReductionMsg* impl_msg);
  __idx_doStressTest_CkReductionMsg = CkRegisterEp("doStressTest(CkReductionMsg* impl_msg)",
     (CkCallFnPtr)_call_doStressTest_CkReductionMsg, CMessage_CkReductionMsg::__idx, __idx);
}
#endif

/* DEFS: array Smooth_TreePiece: TreePiece{
Smooth_TreePiece(CkMigrateMessage* impl_msg);
Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic);
void findSmoothingRadius(int n, const CkCallback &cb);
void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb);
void handleTreeRequest(const TreeRequest_Pointer &req);
void receiveResponse(const Response_Pointer &resp);
void saveInformation(const std::string &prefix, const CkCallback &cb);
void minmaxDensity(const CkCallback &cb);
void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb);
void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_Smooth_TreePiece::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: Smooth_TreePiece(CkMigrateMessage* impl_msg);
 */

/* DEFS: Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic);
 */
void CProxyElement_Smooth_TreePiece::insert(const PeriodicBoundaryConditions<double> &periodic, int onPE)
{ 
  //Marshall: const PeriodicBoundaryConditions<double> &periodic
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(PeriodicBoundaryConditions<double> &)periodic;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(PeriodicBoundaryConditions<double> &)periodic;
  }
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_Smooth_TreePiece::__idx_Smooth_TreePiece_marshall1,onPE);
}

/* DEFS: void findSmoothingRadius(int n, const CkCallback &cb);
 */
void CProxyElement_Smooth_TreePiece::findSmoothingRadius(int n, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: int n, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_findSmoothingRadius_marshall2);
}

/* DEFS: void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb);
 */
void CProxyElement_Smooth_TreePiece::performSmoothOperation(const SmoothOperation &op, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const SmoothOperation &op, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_performSmoothOperation_marshall3);
}

/* DEFS: void handleTreeRequest(const TreeRequest_Pointer &req);
 */
void CProxyElement_Smooth_TreePiece::handleTreeRequest(const TreeRequest_Pointer &req) 
{
  ckCheck();
  //Marshall: const TreeRequest_Pointer &req
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(TreeRequest_Pointer &)req;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(TreeRequest_Pointer &)req;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_handleTreeRequest_marshall4);
}

/* DEFS: void receiveResponse(const Response_Pointer &resp);
 */
void CProxyElement_Smooth_TreePiece::receiveResponse(const Response_Pointer &resp) 
{
  ckCheck();
  //Marshall: const Response_Pointer &resp
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(Response_Pointer &)resp;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(Response_Pointer &)resp;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_receiveResponse_marshall5);
}

/* DEFS: void saveInformation(const std::string &prefix, const CkCallback &cb);
 */
void CProxyElement_Smooth_TreePiece::saveInformation(const std::string &prefix, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const std::string &prefix, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)prefix;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)prefix;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_saveInformation_marshall6);
}

/* DEFS: void minmaxDensity(const CkCallback &cb);
 */
void CProxyElement_Smooth_TreePiece::minmaxDensity(const CkCallback &cb) 
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
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_minmaxDensity_marshall7);
}

/* DEFS: void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb);
 */
void CProxyElement_Smooth_TreePiece::makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|numDensityBins;
    implP|minDensity;
    implP|maxDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|numDensityBins;
    implP|minDensity;
    implP|maxDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_makeDensityHistogram_marshall8);
}

/* DEFS: void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb);
 */
void CProxyElement_Smooth_TreePiece::densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const SmoothOperation &op, double minDensity, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    implP|minDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    implP|minDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_densityCutOperation_marshall9);
}
/* DEFS: Smooth_TreePiece(CkMigrateMessage* impl_msg);
 */
 int CkIndex_Smooth_TreePiece::__idx_Smooth_TreePiece_CkMigrateMessage=0;
void CkIndex_Smooth_TreePiece::_call_Smooth_TreePiece_CkMigrateMessage(void* impl_msg,Smooth_TreePiece * impl_obj)
{
  new (impl_obj) Smooth_TreePiece((CkMigrateMessage*)impl_msg);
}

/* DEFS: Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic);
 */
CkArrayID CProxy_Smooth_TreePiece::ckNew(const PeriodicBoundaryConditions<double> &periodic, const CkArrayOptions &opts)
{ 
  //Marshall: const PeriodicBoundaryConditions<double> &periodic
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(PeriodicBoundaryConditions<double> &)periodic;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(PeriodicBoundaryConditions<double> &)periodic;
  }
   return ckCreateArray((CkArrayMessage *)impl_msg,CkIndex_Smooth_TreePiece::__idx_Smooth_TreePiece_marshall1,opts);
}
 int CkIndex_Smooth_TreePiece::__idx_Smooth_TreePiece_marshall1=0;
void CkIndex_Smooth_TreePiece::_call_Smooth_TreePiece_marshall1(void* impl_msg,Smooth_TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const PeriodicBoundaryConditions<double> &periodic
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  PeriodicBoundaryConditions<double> periodic; implP|periodic;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  new (impl_obj) Smooth_TreePiece(periodic);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void findSmoothingRadius(int n, const CkCallback &cb);
 */
void CProxy_Smooth_TreePiece::findSmoothingRadius(int n, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: int n, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Smooth_TreePiece::__idx_findSmoothingRadius_marshall2);
}
 int CkIndex_Smooth_TreePiece::__idx_findSmoothingRadius_marshall2=0;
void CkIndex_Smooth_TreePiece::_call_findSmoothingRadius_marshall2(void* impl_msg,Smooth_TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: int n, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  int n; implP|n;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->findSmoothingRadius(n, cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb);
 */
void CProxy_Smooth_TreePiece::performSmoothOperation(const SmoothOperation &op, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const SmoothOperation &op, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Smooth_TreePiece::__idx_performSmoothOperation_marshall3);
}
 int CkIndex_Smooth_TreePiece::__idx_performSmoothOperation_marshall3=0;
void CkIndex_Smooth_TreePiece::_call_performSmoothOperation_marshall3(void* impl_msg,Smooth_TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const SmoothOperation &op, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  SmoothOperation op; implP|op;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->performSmoothOperation(op, cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void handleTreeRequest(const TreeRequest_Pointer &req);
 */
void CProxy_Smooth_TreePiece::handleTreeRequest(const TreeRequest_Pointer &req) 
{
  ckCheck();
  //Marshall: const TreeRequest_Pointer &req
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(TreeRequest_Pointer &)req;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(TreeRequest_Pointer &)req;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Smooth_TreePiece::__idx_handleTreeRequest_marshall4);
}
 int CkIndex_Smooth_TreePiece::__idx_handleTreeRequest_marshall4=0;
void CkIndex_Smooth_TreePiece::_call_handleTreeRequest_marshall4(void* impl_msg,Smooth_TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const TreeRequest_Pointer &req
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  TreeRequest_Pointer req; implP|req;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->handleTreeRequest(req);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void receiveResponse(const Response_Pointer &resp);
 */
void CProxy_Smooth_TreePiece::receiveResponse(const Response_Pointer &resp) 
{
  ckCheck();
  //Marshall: const Response_Pointer &resp
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(Response_Pointer &)resp;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(Response_Pointer &)resp;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Smooth_TreePiece::__idx_receiveResponse_marshall5);
}
 int CkIndex_Smooth_TreePiece::__idx_receiveResponse_marshall5=0;
void CkIndex_Smooth_TreePiece::_call_receiveResponse_marshall5(void* impl_msg,Smooth_TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const Response_Pointer &resp
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  Response_Pointer resp; implP|resp;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->receiveResponse(resp);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void saveInformation(const std::string &prefix, const CkCallback &cb);
 */
void CProxy_Smooth_TreePiece::saveInformation(const std::string &prefix, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const std::string &prefix, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)prefix;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)prefix;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Smooth_TreePiece::__idx_saveInformation_marshall6);
}
 int CkIndex_Smooth_TreePiece::__idx_saveInformation_marshall6=0;
void CkIndex_Smooth_TreePiece::_call_saveInformation_marshall6(void* impl_msg,Smooth_TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const std::string &prefix, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  std::string prefix; implP|prefix;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->saveInformation(prefix, cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void minmaxDensity(const CkCallback &cb);
 */
void CProxy_Smooth_TreePiece::minmaxDensity(const CkCallback &cb) 
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
  ckBroadcast(impl_amsg, CkIndex_Smooth_TreePiece::__idx_minmaxDensity_marshall7);
}
 int CkIndex_Smooth_TreePiece::__idx_minmaxDensity_marshall7=0;
void CkIndex_Smooth_TreePiece::_call_minmaxDensity_marshall7(void* impl_msg,Smooth_TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->minmaxDensity(cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb);
 */
void CProxy_Smooth_TreePiece::makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|numDensityBins;
    implP|minDensity;
    implP|maxDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|numDensityBins;
    implP|minDensity;
    implP|maxDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Smooth_TreePiece::__idx_makeDensityHistogram_marshall8);
}
 int CkIndex_Smooth_TreePiece::__idx_makeDensityHistogram_marshall8=0;
void CkIndex_Smooth_TreePiece::_call_makeDensityHistogram_marshall8(void* impl_msg,Smooth_TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  int numDensityBins; implP|numDensityBins;
  double minDensity; implP|minDensity;
  double maxDensity; implP|maxDensity;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->makeDensityHistogram(numDensityBins, minDensity, maxDensity, cb);
  delete (CkMarshallMsg *)impl_msg;
}

/* DEFS: void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb);
 */
void CProxy_Smooth_TreePiece::densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const SmoothOperation &op, double minDensity, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    implP|minDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    implP|minDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Smooth_TreePiece::__idx_densityCutOperation_marshall9);
}
 int CkIndex_Smooth_TreePiece::__idx_densityCutOperation_marshall9=0;
void CkIndex_Smooth_TreePiece::_call_densityCutOperation_marshall9(void* impl_msg,Smooth_TreePiece * impl_obj)
{
  //Unmarshall pup'd fields: const SmoothOperation &op, double minDensity, const CkCallback &cb
  char *impl_buf=((CkMarshallMsg *)impl_msg)->msgBuf;
  PUP::fromMem implP(impl_buf);
  SmoothOperation op; implP|op;
  double minDensity; implP|minDensity;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  //Unmarshall arrays:
  impl_obj->densityCutOperation(op, minDensity, cb);
  delete (CkMarshallMsg *)impl_msg;
}
/* DEFS: Smooth_TreePiece(CkMigrateMessage* impl_msg);
 */

/* DEFS: Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic);
 */

/* DEFS: void findSmoothingRadius(int n, const CkCallback &cb);
 */
void CProxySection_Smooth_TreePiece::findSmoothingRadius(int n, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: int n, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_findSmoothingRadius_marshall2);
}

/* DEFS: void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb);
 */
void CProxySection_Smooth_TreePiece::performSmoothOperation(const SmoothOperation &op, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const SmoothOperation &op, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_performSmoothOperation_marshall3);
}

/* DEFS: void handleTreeRequest(const TreeRequest_Pointer &req);
 */
void CProxySection_Smooth_TreePiece::handleTreeRequest(const TreeRequest_Pointer &req) 
{
  ckCheck();
  //Marshall: const TreeRequest_Pointer &req
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(TreeRequest_Pointer &)req;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(TreeRequest_Pointer &)req;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_handleTreeRequest_marshall4);
}

/* DEFS: void receiveResponse(const Response_Pointer &resp);
 */
void CProxySection_Smooth_TreePiece::receiveResponse(const Response_Pointer &resp) 
{
  ckCheck();
  //Marshall: const Response_Pointer &resp
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(Response_Pointer &)resp;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(Response_Pointer &)resp;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_receiveResponse_marshall5);
}

/* DEFS: void saveInformation(const std::string &prefix, const CkCallback &cb);
 */
void CProxySection_Smooth_TreePiece::saveInformation(const std::string &prefix, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const std::string &prefix, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)prefix;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)prefix;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_saveInformation_marshall6);
}

/* DEFS: void minmaxDensity(const CkCallback &cb);
 */
void CProxySection_Smooth_TreePiece::minmaxDensity(const CkCallback &cb) 
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
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_minmaxDensity_marshall7);
}

/* DEFS: void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb);
 */
void CProxySection_Smooth_TreePiece::makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|numDensityBins;
    implP|minDensity;
    implP|maxDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|numDensityBins;
    implP|minDensity;
    implP|maxDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_makeDensityHistogram_marshall8);
}

/* DEFS: void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb);
 */
void CProxySection_Smooth_TreePiece::densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb) 
{
  ckCheck();
  //Marshall: const SmoothOperation &op, double minDensity, const CkCallback &cb
  int impl_off=0,impl_arrstart=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    implP|minDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=new (impl_off,0)CkMarshallMsg;
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(SmoothOperation &)op;
    implP|minDensity;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Smooth_TreePiece::__idx_densityCutOperation_marshall9);
}
#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_Smooth_TreePiece::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size);
  _REGISTER_BASE(__idx, CkIndex_TreePiece::__idx);
// REG: Smooth_TreePiece(CkMigrateMessage* impl_msg);
  __idx_Smooth_TreePiece_CkMigrateMessage = CkRegisterEp("Smooth_TreePiece(CkMigrateMessage* impl_msg)",
     (CkCallFnPtr)_call_Smooth_TreePiece_CkMigrateMessage, 0, __idx);
  CkRegisterMigCtor(__idx, __idx_Smooth_TreePiece_CkMigrateMessage);

// REG: Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic);
  __idx_Smooth_TreePiece_marshall1 = CkRegisterEp("Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic)",
     (CkCallFnPtr)_call_Smooth_TreePiece_marshall1, CkMarshallMsg::__idx, __idx);

// REG: void findSmoothingRadius(int n, const CkCallback &cb);
  __idx_findSmoothingRadius_marshall2 = CkRegisterEp("findSmoothingRadius(int n, const CkCallback &cb)",
     (CkCallFnPtr)_call_findSmoothingRadius_marshall2, CkMarshallMsg::__idx, __idx);

// REG: void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb);
  __idx_performSmoothOperation_marshall3 = CkRegisterEp("performSmoothOperation(const SmoothOperation &op, const CkCallback &cb)",
     (CkCallFnPtr)_call_performSmoothOperation_marshall3, CkMarshallMsg::__idx, __idx);

// REG: void handleTreeRequest(const TreeRequest_Pointer &req);
  __idx_handleTreeRequest_marshall4 = CkRegisterEp("handleTreeRequest(const TreeRequest_Pointer &req)",
     (CkCallFnPtr)_call_handleTreeRequest_marshall4, CkMarshallMsg::__idx, __idx);

// REG: void receiveResponse(const Response_Pointer &resp);
  __idx_receiveResponse_marshall5 = CkRegisterEp("receiveResponse(const Response_Pointer &resp)",
     (CkCallFnPtr)_call_receiveResponse_marshall5, CkMarshallMsg::__idx, __idx);

// REG: void saveInformation(const std::string &prefix, const CkCallback &cb);
  __idx_saveInformation_marshall6 = CkRegisterEp("saveInformation(const std::string &prefix, const CkCallback &cb)",
     (CkCallFnPtr)_call_saveInformation_marshall6, CkMarshallMsg::__idx, __idx);

// REG: void minmaxDensity(const CkCallback &cb);
  __idx_minmaxDensity_marshall7 = CkRegisterEp("minmaxDensity(const CkCallback &cb)",
     (CkCallFnPtr)_call_minmaxDensity_marshall7, CkMarshallMsg::__idx, __idx);

// REG: void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb);
  __idx_makeDensityHistogram_marshall8 = CkRegisterEp("makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb)",
     (CkCallFnPtr)_call_makeDensityHistogram_marshall8, CkMarshallMsg::__idx, __idx);

// REG: void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb);
  __idx_densityCutOperation_marshall9 = CkRegisterEp("densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb)",
     (CkCallFnPtr)_call_densityCutOperation_marshall9, CkMarshallMsg::__idx, __idx);
}
#endif



#ifndef CK_TEMPLATES_ONLY
void _registerParallelSmooth(void)
{
  static int _done = 0; if(_done) return; _done = 1;
      _registerNChilada();

/* REG: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
threaded void loadSortBuild(void);
threaded void doCalculations(void);
threaded void makeHistogram(CkReductionMsg* impl_msg);
threaded void doStressTest(CkReductionMsg* impl_msg);
};
*/
  CkIndex_Main::__register("Main", sizeof(Main));

/* REG: array Smooth_TreePiece: TreePiece{
Smooth_TreePiece(CkMigrateMessage* impl_msg);
Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic);
void findSmoothingRadius(int n, const CkCallback &cb);
void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb);
void handleTreeRequest(const TreeRequest_Pointer &req);
void receiveResponse(const Response_Pointer &resp);
void saveInformation(const std::string &prefix, const CkCallback &cb);
void minmaxDensity(const CkCallback &cb);
void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb);
void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb);
};
*/
  CkIndex_Smooth_TreePiece::__register("Smooth_TreePiece", sizeof(Smooth_TreePiece));

      initRequestResponsePUP();

      initSPHPUP();

}
extern "C" void CkRegisterMainModule(void) {
  _registerParallelSmooth();
  _REGISTER_DONE();
}
#endif
