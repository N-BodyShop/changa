/** @file NChilada.cpp
 Implementation of some basic functions of the NChilada library.
 @author Graeme Lufkin (gwl@u.washington.edu)
*/

#include <iostream>

#include "conv-ccs.h"
#include "pup_stl.h"
#include "pup_toNetwork4.h"

#include "NChilada.h"

#include "Sorter.h"
#include "DataManager.h"
#include "TreePiece.h"
#include "Tree.h"

using std::string;

int verbosity;

CkArrayID treePieceID;
CkGroupID dataManagerID;

Space3D<double> space;

CkReduction::reducerType boxGrowthReduction;
CkReduction::reducerType minmaxReduction;
CkReduction::reducerType callbackReduction;

/// Combine reduction messages to grow a box
CkReductionMsg* boxGrowth(int nMsg, CkReductionMsg** msgs) {
	if(nMsg < 1)
		cerr << "Small reduction message!" << endl;
	
	OrientedBox<double>* pbox = static_cast<OrientedBox<double> *>(msgs[0]->getData());
	OrientedBox<double>* msgpbox;
	for(int i = 1; i < nMsg; i++) {
		msgpbox = static_cast<OrientedBox<double> *>(msgs[i]->getData());
		pbox->grow(msgpbox->lesser_corner);
		pbox->grow(msgpbox->greater_corner);
	}
	
	return CkReductionMsg::buildNew(sizeof(OrientedBox<double>), pbox);
}

/// Combine reduction messages to get min/max pair
template <typename T>
CkReductionMsg* minmax(int nMsg, CkReductionMsg** msgs) {
	if(nMsg < 1)
		cerr << "Small reduction message!" << endl;
	
	T* pminmax = static_cast<T *>(msgs[0]->getData());
	T* msgpminmax;
	for(int i = 1; i < nMsg; i++) {
		msgpminmax = static_cast<T *>(msgs[i]->getData());
		if(msgpminmax[0] < pminmax[0])
			pminmax[0] = msgpminmax[0];
		if(msgpminmax[1] > pminmax[1])
			pminmax[1] = msgpminmax[1];
	}
	
	return CkReductionMsg::buildNew(2 * sizeof(T), pminmax);
}

/// Return a single object, given many copies of it
template <typename T>
CkReductionMsg* same(int nMsg, CkReductionMsg** msgs) {
	return CkReductionMsg::buildNew(sizeof(T), static_cast<T *>(msgs[0]->getData()));
}

/// Register my reductions
void registerNChiladaReductions() {
	boxGrowthReduction = CkReduction::addReducer(boxGrowth);
	minmaxReduction = CkReduction::addReducer(minmax<double>);
	callbackReduction = CkReduction::addReducer(same<CkCallback>);
}

/*
void perParticleCompiled(void* param, void* msg) {
	CcsDelayedReply* token = (CcsDelayedReply *) param;
	CkReductionMsg* m = (CkReductionMsg *) msg;
	
	string s;
	if(m->getSize() == 0)
		s = "Success";
	else
		s = "Failure";
	CcsSendDelayedReply(*token, s.length(), s.c_str());
	delete m;
	cout << "Function applied" << endl;
}

extern "C" void applyPerParticleHandler(byte* msg) {
	cout << "CCS: Got request to apply per-particle function" << endl;
	//get string out of msg
	PUP_toNetwork4_unpack p(msg + CmiMsgHeaderSizeBytes);
	string s;
	p | s;
	cout << "CCS: Code is:\n" << s << endl;
	CProxy_DataManager dataManager(dataManagerID);
	CcsDelayedReply token = CcsDelayReply();
	dataManager.compilePerParticleFunction(s, CkCallback(perParticleCompiled, &token));
}

void initializeNChiladaCcs() {
	cout << "Registering handlers" << endl;
	CcsRegisterHandler("applyPerParticleHandler", (CmiHandler) applyPerParticleHandler);
}
*/
#include "NChilada.def.h"
