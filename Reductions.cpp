/** @file Reductions.cpp
 */
 
#include "Reductions.h"

#include "OrientedBox.h"

#include "dumpframe.h"

CkReduction::reducerType growOrientedBox_float;
CkReduction::reducerType growOrientedBox_double;

CkReduction::reducerType minmax_int;
CkReduction::reducerType minmax_long;
CkReduction::reducerType minmax_float;
CkReduction::reducerType minmax_double;

CkReduction::reducerType max_count;

CkReduction::reducerType callbackReduction;
CkReduction::reducerType boxReduction;

CkReduction::reducerType dfImageReduction;

/// Combine reduction messages to grow a box
template <typename T>
CkReductionMsg* boxGrowth(int nMsg, CkReductionMsg** msgs) {

  CkAssert(msgs[0]->getSize() == sizeof(OrientedBox<T>)); 
	OrientedBox<T>* pbox = static_cast<OrientedBox<T> *>(msgs[0]->getData());
	OrientedBox<T>* msgpbox;
	for(int i = 1; i < nMsg; i++) {
          CkAssert(msgs[i]->getSize() == sizeof(OrientedBox<T>)); 
		msgpbox = static_cast<OrientedBox<T> *>(msgs[i]->getData());
		pbox->grow(*msgpbox);
	}
	
	return CkReductionMsg::buildNew(sizeof(OrientedBox<T>), pbox);
}

/// Combine reduction messages to get min/max pair
template <typename T>
CkReductionMsg* minmax(int nMsg, CkReductionMsg** msgs) {
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

/// Reduction for determining both the maximum timestep rung
/// (iMaxRung) and the number of particles in that rung (nMaxRung).
///
CkReductionMsg* max_count_reduce(int nMsg, CkReductionMsg** msgs) {
    int64_t* pmaxcount = static_cast<int64_t *>(msgs[0]->getData());
    int iMaxRung = pmaxcount[0];  // maxmimum rung
    int64_t nMaxRung = pmaxcount[1];  // count in maximum rung
    for(int i = 1; i < nMsg; i++) {
	pmaxcount = static_cast<int64_t *>(msgs[i]->getData());
	if(pmaxcount[0] > iMaxRung) {
	    iMaxRung = pmaxcount[0];
	    nMaxRung = pmaxcount[1];
	    }
	else if(pmaxcount[0] == iMaxRung) {
	    nMaxRung += pmaxcount[1];
	    }
	}
    int64_t newcount[2];
    newcount[0] = iMaxRung;
    newcount[1] = nMaxRung;
    return CkReductionMsg::buildNew(2 * sizeof(int64_t), newcount);
}

/// Return a single object, given many copies of it
template <typename T>
CkReductionMsg* same(int nMsg, CkReductionMsg** msgs) {
	return CkReductionMsg::buildNew(sizeof(T), static_cast<T *>(msgs[0]->getData()));
}

// Merge images for DumpFrame
// Messages consist of a struct inDumpFrame header followed by the
// image data.
CkReductionMsg* dfImageReducer(int nMsg, CkReductionMsg** msgs) {
    struct inDumpFrame *in = (struct inDumpFrame *)msgs[0]->getData();
    int nImage1 = in->nxPix*in->nyPix*sizeof(DFIMAGE);
    
    CkAssert(msgs[0]->getSize() == (sizeof(struct inDumpFrame) + nImage1));
    void *ImageOut = ((char *)msgs[0]->getData()) + sizeof(struct inDumpFrame);
    
    for(int i = 1; i < nMsg; i++) {
	void *Image2 = ((char *)msgs[i]->getData())
	    + sizeof(struct inDumpFrame);
	
	int nImage2 =  in->nxPix*in->nyPix*sizeof(DFIMAGE);
	CkAssert(msgs[i]->getSize() == (sizeof(struct inDumpFrame) + nImage2));
	dfMergeImage( in, ImageOut, &nImage1, Image2, &nImage2);
	}
    
    return CkReductionMsg::buildNew(msgs[0]->getSize(), msgs[0]->getData());
}

void registerReductions() {
	growOrientedBox_float = CkReduction::addReducer(boxGrowth<float>);
	growOrientedBox_double = CkReduction::addReducer(boxGrowth<double>);
	
	minmax_int = CkReduction::addReducer(minmax<int>);
	minmax_long = CkReduction::addReducer(minmax<int64_t>);
	minmax_float = CkReduction::addReducer(minmax<float>);
	minmax_double = CkReduction::addReducer(minmax<double>);
	max_count = CkReduction::addReducer(max_count_reduce);
	callbackReduction = CkReduction::addReducer(same<CkCallback>);
	boxReduction = CkReduction::addReducer(same<OrientedBox<float> >);
	dfImageReduction = CkReduction::addReducer(dfImageReducer);
	
}

#include "Reductions.def.h"
