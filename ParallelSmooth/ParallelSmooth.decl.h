#ifndef _DECL_ParallelSmooth_H_
#define _DECL_ParallelSmooth_H_
#include "charm++.h"
#include "NChilada.decl.h"

/* DECLS: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
threaded void loadSortBuild(void);
threaded void doCalculations(void);
threaded void makeHistogram(CkReductionMsg* impl_msg);
threaded void doStressTest(CkReductionMsg* impl_msg);
};
 */
 class Main;
/* --------------- index object ------------------ */
class CkIndex_Main{
  public:
    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: Main(CkArgMsg* impl_msg);
 */
    static int __idx_Main_CkArgMsg;
    static int ckNew(CkArgMsg* impl_msg) { return __idx_Main_CkArgMsg; }
    static void _call_Main_CkArgMsg(void* impl_msg,Main* impl_obj);

/* DECLS: threaded void loadSortBuild(void);
 */
    static int __idx_loadSortBuild_void;
    static int loadSortBuild(void) { return __idx_loadSortBuild_void; }
    static void _call_loadSortBuild_void(void* impl_msg,Main* impl_obj);
    static void _callthr_loadSortBuild_void(CkThrCallArg *);

/* DECLS: threaded void doCalculations(void);
 */
    static int __idx_doCalculations_void;
    static int doCalculations(void) { return __idx_doCalculations_void; }
    static void _call_doCalculations_void(void* impl_msg,Main* impl_obj);
    static void _callthr_doCalculations_void(CkThrCallArg *);

/* DECLS: threaded void makeHistogram(CkReductionMsg* impl_msg);
 */
    static int __idx_makeHistogram_CkReductionMsg;
    static int makeHistogram(CkReductionMsg* impl_msg) { return __idx_makeHistogram_CkReductionMsg; }
    static void _call_makeHistogram_CkReductionMsg(void* impl_msg,Main* impl_obj);
    static void _callthr_makeHistogram_CkReductionMsg(CkThrCallArg *);

/* DECLS: threaded void doStressTest(CkReductionMsg* impl_msg);
 */
    static int __idx_doStressTest_CkReductionMsg;
    static int doStressTest(CkReductionMsg* impl_msg) { return __idx_doStressTest_CkReductionMsg; }
    static void _call_doStressTest_CkReductionMsg(void* impl_msg,Main* impl_obj);
    static void _callthr_doStressTest_CkReductionMsg(CkThrCallArg *);
};
/* --------------- element proxy ------------------ */
class CProxy_Main:public CProxy_Chare{
  public:
    CProxy_Main(void) {};
    CProxy_Main(CkChareID __cid) : CProxy_Chare(__cid){  }
    CProxy_Main(const Chare *c) : CProxy_Chare(c){  }
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
    Main *ckLocal(void) const
     { return (Main *)CkLocalChare(&ckGetChareID()); }
/* DECLS: Main(CkArgMsg* impl_msg);
 */
    static CkChareID ckNew(CkArgMsg* impl_msg, int onPE=CK_PE_ANY);
    static void ckNew(CkArgMsg* impl_msg, CkChareID* pcid, int onPE=CK_PE_ANY);
    CProxy_Main(CkArgMsg* impl_msg, int onPE=CK_PE_ANY);

/* DECLS: threaded void loadSortBuild(void);
 */
    void loadSortBuild(void);

/* DECLS: threaded void doCalculations(void);
 */
    void doCalculations(void);

/* DECLS: threaded void makeHistogram(CkReductionMsg* impl_msg);
 */
    void makeHistogram(CkReductionMsg* impl_msg);

/* DECLS: threaded void doStressTest(CkReductionMsg* impl_msg);
 */
    void doStressTest(CkReductionMsg* impl_msg);
};
PUPmarshall(CProxy_Main);
typedef CBaseT<Chare,CProxy_Main>  CBase_Main;

/* DECLS: array Smooth_TreePiece: TreePiece{
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
 class Smooth_TreePiece;
/* --------------- index object ------------------ */
class CkIndex_Smooth_TreePiece{
  public:
    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: Smooth_TreePiece(CkMigrateMessage* impl_msg);
 */
    static int __idx_Smooth_TreePiece_CkMigrateMessage;
    static int ckNew(CkMigrateMessage* impl_msg) { return __idx_Smooth_TreePiece_CkMigrateMessage; }
    static void _call_Smooth_TreePiece_CkMigrateMessage(void* impl_msg,Smooth_TreePiece* impl_obj);

/* DECLS: Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic);
 */
    static int __idx_Smooth_TreePiece_marshall1;
    static int ckNew(const PeriodicBoundaryConditions<double> &periodic) { return __idx_Smooth_TreePiece_marshall1; }
    static void _call_Smooth_TreePiece_marshall1(void* impl_msg,Smooth_TreePiece* impl_obj);

/* DECLS: void findSmoothingRadius(int n, const CkCallback &cb);
 */
    static int __idx_findSmoothingRadius_marshall2;
    static int findSmoothingRadius(int n, const CkCallback &cb) { return __idx_findSmoothingRadius_marshall2; }
    static void _call_findSmoothingRadius_marshall2(void* impl_msg,Smooth_TreePiece* impl_obj);

/* DECLS: void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb);
 */
    static int __idx_performSmoothOperation_marshall3;
    static int performSmoothOperation(const SmoothOperation &op, const CkCallback &cb) { return __idx_performSmoothOperation_marshall3; }
    static void _call_performSmoothOperation_marshall3(void* impl_msg,Smooth_TreePiece* impl_obj);

/* DECLS: void handleTreeRequest(const TreeRequest_Pointer &req);
 */
    static int __idx_handleTreeRequest_marshall4;
    static int handleTreeRequest(const TreeRequest_Pointer &req) { return __idx_handleTreeRequest_marshall4; }
    static void _call_handleTreeRequest_marshall4(void* impl_msg,Smooth_TreePiece* impl_obj);

/* DECLS: void receiveResponse(const Response_Pointer &resp);
 */
    static int __idx_receiveResponse_marshall5;
    static int receiveResponse(const Response_Pointer &resp) { return __idx_receiveResponse_marshall5; }
    static void _call_receiveResponse_marshall5(void* impl_msg,Smooth_TreePiece* impl_obj);

/* DECLS: void saveInformation(const std::string &prefix, const CkCallback &cb);
 */
    static int __idx_saveInformation_marshall6;
    static int saveInformation(const std::string &prefix, const CkCallback &cb) { return __idx_saveInformation_marshall6; }
    static void _call_saveInformation_marshall6(void* impl_msg,Smooth_TreePiece* impl_obj);

/* DECLS: void minmaxDensity(const CkCallback &cb);
 */
    static int __idx_minmaxDensity_marshall7;
    static int minmaxDensity(const CkCallback &cb) { return __idx_minmaxDensity_marshall7; }
    static void _call_minmaxDensity_marshall7(void* impl_msg,Smooth_TreePiece* impl_obj);

/* DECLS: void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb);
 */
    static int __idx_makeDensityHistogram_marshall8;
    static int makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb) { return __idx_makeDensityHistogram_marshall8; }
    static void _call_makeDensityHistogram_marshall8(void* impl_msg,Smooth_TreePiece* impl_obj);

/* DECLS: void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb);
 */
    static int __idx_densityCutOperation_marshall9;
    static int densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb) { return __idx_densityCutOperation_marshall9; }
    static void _call_densityCutOperation_marshall9(void* impl_msg,Smooth_TreePiece* impl_obj);
};
/* --------------- element proxy ------------------ */
 class CProxyElement_Smooth_TreePiece : public CProxyElement_TreePiece{
  public:
    CProxyElement_Smooth_TreePiece(void) {}
    CProxyElement_Smooth_TreePiece(const ArrayElement *e) : CProxyElement_TreePiece(e){  }
    void ckDelegate(CkGroupID to) {
      CProxyElement_TreePiece::ckDelegate(to);
    }
    void ckUndelegate(void) {
      CProxyElement_TreePiece::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxyElement_TreePiece::pup(p);
    }
    CK_DISAMBIG_ARRAY_ELEMENT(CProxyElement_TreePiece)
    Smooth_TreePiece *ckLocal(void) const
      { return (Smooth_TreePiece *)CProxyElement_TreePiece::ckLocal(); }
    CProxyElement_Smooth_TreePiece(const CkArrayID &aid,const CkArrayIndex1D &idx,CkGroupID dTo)
        :CProxyElement_TreePiece(aid,idx,dTo) {}
    CProxyElement_Smooth_TreePiece(const CkArrayID &aid,const CkArrayIndex1D &idx)
        :CProxyElement_TreePiece(aid,idx) {}
/* DECLS: Smooth_TreePiece(CkMigrateMessage* impl_msg);
 */

/* DECLS: Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic);
 */
    void insert(const PeriodicBoundaryConditions<double> &periodic, int onPE=-1);
/* DECLS: void findSmoothingRadius(int n, const CkCallback &cb);
 */
    void findSmoothingRadius(int n, const CkCallback &cb) ;

/* DECLS: void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb);
 */
    void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb) ;

/* DECLS: void handleTreeRequest(const TreeRequest_Pointer &req);
 */
    void handleTreeRequest(const TreeRequest_Pointer &req) ;

/* DECLS: void receiveResponse(const Response_Pointer &resp);
 */
    void receiveResponse(const Response_Pointer &resp) ;

/* DECLS: void saveInformation(const std::string &prefix, const CkCallback &cb);
 */
    void saveInformation(const std::string &prefix, const CkCallback &cb) ;

/* DECLS: void minmaxDensity(const CkCallback &cb);
 */
    void minmaxDensity(const CkCallback &cb) ;

/* DECLS: void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb);
 */
    void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb) ;

/* DECLS: void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb);
 */
    void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb) ;
};
PUPmarshall(CProxyElement_Smooth_TreePiece);
/* ---------------- collective proxy -------------- */
 class CProxy_Smooth_TreePiece : public CProxy_TreePiece{
  public:
    CProxy_Smooth_TreePiece(void) {}
    CProxy_Smooth_TreePiece(const ArrayElement *e) : CProxy_TreePiece(e){  }
    void ckDelegate(CkGroupID to) {
      CProxy_TreePiece::ckDelegate(to);
    }
    void ckUndelegate(void) {
      CProxy_TreePiece::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxy_TreePiece::pup(p);
    }
    CK_DISAMBIG_ARRAY(CProxy_TreePiece)
    static CkArrayID ckNew(void) {return ckCreateEmptyArray();}
//Generalized array indexing:
    CProxyElement_Smooth_TreePiece operator [] (const CkArrayIndex1D &idx) const
        {return CProxyElement_Smooth_TreePiece(ckGetArrayID(), idx, ckDelegatedIdx());}
    CProxyElement_Smooth_TreePiece operator() (const CkArrayIndex1D &idx) const
        {return CProxyElement_Smooth_TreePiece(ckGetArrayID(), idx, ckDelegatedIdx());}
    CProxyElement_Smooth_TreePiece operator [] (int idx) const 
        {return CProxyElement_Smooth_TreePiece(ckGetArrayID(), CkArrayIndex1D(idx), ckDelegatedIdx());}
    CProxyElement_Smooth_TreePiece operator () (int idx) const 
        {return CProxyElement_Smooth_TreePiece(ckGetArrayID(), CkArrayIndex1D(idx), ckDelegatedIdx());}
    CProxy_Smooth_TreePiece(const CkArrayID &aid,CkGroupID dTo) 
        :CProxy_TreePiece(aid,dTo) {}
    CProxy_Smooth_TreePiece(const CkArrayID &aid) 
        :CProxy_TreePiece(aid) {}
/* DECLS: Smooth_TreePiece(CkMigrateMessage* impl_msg);
 */

/* DECLS: Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic);
 */
    static CkArrayID ckNew(const PeriodicBoundaryConditions<double> &periodic, const CkArrayOptions &opts);

/* DECLS: void findSmoothingRadius(int n, const CkCallback &cb);
 */
    void findSmoothingRadius(int n, const CkCallback &cb) ;

/* DECLS: void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb);
 */
    void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb) ;

/* DECLS: void handleTreeRequest(const TreeRequest_Pointer &req);
 */
    void handleTreeRequest(const TreeRequest_Pointer &req) ;

/* DECLS: void receiveResponse(const Response_Pointer &resp);
 */
    void receiveResponse(const Response_Pointer &resp) ;

/* DECLS: void saveInformation(const std::string &prefix, const CkCallback &cb);
 */
    void saveInformation(const std::string &prefix, const CkCallback &cb) ;

/* DECLS: void minmaxDensity(const CkCallback &cb);
 */
    void minmaxDensity(const CkCallback &cb) ;

/* DECLS: void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb);
 */
    void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb) ;

/* DECLS: void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb);
 */
    void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb) ;
};
PUPmarshall(CProxy_Smooth_TreePiece);
/* ---------------- section proxy -------------- */
 class CProxySection_Smooth_TreePiece : public CProxySection_TreePiece{
  public:
    CProxySection_Smooth_TreePiece(void) {}
    void ckDelegate(CkGroupID to) {
      CProxySection_TreePiece::ckDelegate(to);
    }
    void ckUndelegate(void) {
      CProxySection_TreePiece::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxySection_TreePiece::pup(p);
    }
    CK_DISAMBIG_ARRAY_SECTION(CProxySection_TreePiece)
//Generalized array indexing:
    CProxyElement_Smooth_TreePiece operator [] (const CkArrayIndex1D &idx) const
        {return CProxyElement_Smooth_TreePiece(ckGetArrayID(), idx, ckDelegatedIdx());}
    CProxyElement_Smooth_TreePiece operator() (const CkArrayIndex1D &idx) const
        {return CProxyElement_Smooth_TreePiece(ckGetArrayID(), idx, ckDelegatedIdx());}
    CProxyElement_Smooth_TreePiece operator [] (int idx) const 
        {return CProxyElement_Smooth_TreePiece(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], ckDelegatedIdx());}
    CProxyElement_Smooth_TreePiece operator () (int idx) const 
        {return CProxyElement_Smooth_TreePiece(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], ckDelegatedIdx());}
    CProxySection_Smooth_TreePiece(const CkArrayID &aid, CkArrayIndexMax *elems, int nElems, CkGroupID dTo) 
        :CProxySection_TreePiece(aid,elems,nElems,dTo) {}
    CProxySection_Smooth_TreePiece(const CkArrayID &aid, CkArrayIndexMax *elems, int nElems) 
        :CProxySection_TreePiece(aid,elems,nElems) {}
    CProxySection_Smooth_TreePiece(const CkSectionID &sid)       :CProxySection_TreePiece(sid) {}
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndexMax *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
/* DECLS: Smooth_TreePiece(CkMigrateMessage* impl_msg);
 */

/* DECLS: Smooth_TreePiece(const PeriodicBoundaryConditions<double> &periodic);
 */

/* DECLS: void findSmoothingRadius(int n, const CkCallback &cb);
 */
    void findSmoothingRadius(int n, const CkCallback &cb) ;

/* DECLS: void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb);
 */
    void performSmoothOperation(const SmoothOperation &op, const CkCallback &cb) ;

/* DECLS: void handleTreeRequest(const TreeRequest_Pointer &req);
 */
    void handleTreeRequest(const TreeRequest_Pointer &req) ;

/* DECLS: void receiveResponse(const Response_Pointer &resp);
 */
    void receiveResponse(const Response_Pointer &resp) ;

/* DECLS: void saveInformation(const std::string &prefix, const CkCallback &cb);
 */
    void saveInformation(const std::string &prefix, const CkCallback &cb) ;

/* DECLS: void minmaxDensity(const CkCallback &cb);
 */
    void minmaxDensity(const CkCallback &cb) ;

/* DECLS: void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb);
 */
    void makeDensityHistogram(int numDensityBins, double minDensity, double maxDensity, const CkCallback &cb) ;

/* DECLS: void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb);
 */
    void densityCutOperation(const SmoothOperation &op, double minDensity, const CkCallback &cb) ;
};
PUPmarshall(CProxySection_Smooth_TreePiece);
typedef CBaseT<TreePiece,CProxy_Smooth_TreePiece>  CBase_Smooth_TreePiece;



extern void _registerParallelSmooth(void);
extern "C" void CkRegisterMainModule(void);
#endif
