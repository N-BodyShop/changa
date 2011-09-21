/*
 *  Bucket.C
 *  
 *
 *  Created by Aditya Devarakonda on 7/25/11.
 *
 *
 */

#include "Bucket.h"

	LBBucket::LBBucket(){
		load = 0;
		hilbertID = 0;
		numobjs = 0;
	}

	bool LBBucket::operator<=(const LBBucket &b) const{
		return hilbertID <= b.hilbertID;
	}
	bool LBBucket::operator>=(const LBBucket &b) const{
		return hilbertID >= b.hilbertID;
	}
	bool LBBucket::operator>(const LBBucket &b) const{
		return !(hilbertID <= b.hilbertID);
	}
	bool LBBucket::operator<(const LBBucket &b) const{
		return !(hilbertID >= b.hilbertID);
	}
	
	void LBBucket::setLoad(float l){
		//CkPrintf("Load = %f\n",l);
		CkAssert(l > 0);
		load += l;
	}

	void LBBucket::setTP(TPObject* tpobj, int len){
		CkAssert(tpobj != NULL);
		tp = tpobj;
		numobjs = len;
	}
		
	void LBBucket::setIndex(int i){
		tpStartIndex = i;
	}

	int LBBucket::getIndex(){
		return tpStartIndex;
	}

	float LBBucket::getLoad(){
		return load;
	}

	long LBBucket::getHilbertID(){
		return hilbertID;
	}

	int LBBucket::getNumTPs(){
		return numobjs;
	}

	TPObject* LBBucket::getTPs(){
		return tp;
	}
