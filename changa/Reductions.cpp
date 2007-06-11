/** @file Reductions.cpp
 */
 
#include "Reductions.h"

#include "OrientedBox.h"

#include "dumpframe.h"

CkReduction::reducerType growOrientedBox_float;
CkReduction::reducerType growOrientedBox_double;

CkReduction::reducerType minmax_int;
CkReduction::reducerType minmax_float;
CkReduction::reducerType minmax_double;

CkReduction::reducerType callbackReduction;
CkReduction::reducerType boxReduction;

CkReduction::reducerType dfImageReduction;

/// Combine reduction messages to grow a box
template <typename T>
CkReductionMsg* boxGrowth(int nMsg, CkReductionMsg** msgs) {
	OrientedBox<T>* pbox = static_cast<OrientedBox<T> *>(msgs[0]->getData());
	OrientedBox<T>* msgpbox;
	for(int i = 1; i < nMsg; i++) {
		msgpbox = static_cast<OrientedBox<T> *>(msgs[i]->getData());
		pbox->grow(msgpbox->lesser_corner);
		pbox->grow(msgpbox->greater_corner);
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

/// Return a single object, given many copies of it
template <typename T>
CkReductionMsg* same(int nMsg, CkReductionMsg** msgs) {
	return CkReductionMsg::buildNew(sizeof(T), static_cast<T *>(msgs[0]->getData()));
}

// Merge images for DumpFrame
// Messages consist of a struct inDumpFrame header followed by the
// image data.
CkReductionMsg* dfImageReducer(int nMsg, CkReductionMsg** msgs) {
    CkAssert(msgs[0]->getSize() == (sizeof(struct inDumpFrame)
				    + DF_NBYTEDUMPFRAME));
    struct inDumpFrame *in = (struct inDumpFrame *)msgs[0]->getData();
    void *ImageOut = ((char *)msgs[0]->getData()) + sizeof(struct inDumpFrame);
    
    for(int i = 1; i < nMsg; i++) {
	void *Image2 = ((char *)msgs[i]->getData())
	    + sizeof(struct inDumpFrame);
	
	int nImage1 = DF_NBYTEDUMPFRAME;
	int nImage2 = DF_NBYTEDUMPFRAME;
	dfMergeImage( in, ImageOut, &nImage1, Image2, &nImage2);
	}
    
    return CkReductionMsg::buildNew(msgs[0]->getSize(), msgs[0]->getData());
}

void registerReductions() {
	growOrientedBox_float = CkReduction::addReducer(boxGrowth<float>);
	growOrientedBox_double = CkReduction::addReducer(boxGrowth<double>);
	
	minmax_int = CkReduction::addReducer(minmax<int>);
	minmax_float = CkReduction::addReducer(minmax<float>);
	minmax_double = CkReduction::addReducer(minmax<double>);
	callbackReduction = CkReduction::addReducer(same<CkCallback>);
	boxReduction = CkReduction::addReducer(same<OrientedBox<float> >);
	dfImageReduction = CkReduction::addReducer(dfImageReducer);
	
}

#include "Reductions.def.h"
