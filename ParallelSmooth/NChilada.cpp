/** \file NChilada.cpp
 Implementation of some basic functions of the NChilada library.
 \author Graeme Lufkin (gwl@u.washington.edu)
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

CkReduction::reducerType boxGrowthReduction;

/// Combine reduction messages to grow a box
CkReductionMsg* boxGrowth(int nMsg, CkReductionMsg** msgs) {
	if(nMsg < 1)
		cerr << "Small reduction message!" << endl;
	
	OrientedBox<float>* pbox = static_cast<OrientedBox<float> *>(msgs[0]->getData());
	OrientedBox<float>* msgpbox;
	for(int i = 1; i < nMsg; i++) {
		msgpbox = static_cast<OrientedBox<float> *>(msgs[i]->getData());
		pbox->grow(msgpbox->lesser_corner);
		pbox->grow(msgpbox->greater_corner);
	}
	
	return CkReductionMsg::buildNew(sizeof(OrientedBox<float>), pbox);
}

/// Register my reductions
void registerNChiladaReductions() {
	boxGrowthReduction = CkReduction::addReducer(boxGrowth);
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
