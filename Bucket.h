/*
 *  Bucket.h
 *  
 *
 *  Created by Aditya Devarakonda on 7/25/11.
 *  
 *
 */

#ifndef _BUCKET_H_
#define _BUCKET_H_
#include "MapStructures.h"
class LBBucket{

public:
	float load;
	CmiUInt8 hilbertID;
	TPObject *tp;
	int numobjs;
	int tpStartIndex;
	float x,y,z;
	
	LBBucket();
	bool operator<=(const LBBucket &b) const;
	bool operator>=(const LBBucket &b) const;
	bool operator>(const LBBucket &b) const;
	bool operator<(const LBBucket &b) const;
	void setLoad(float l);
	void setTP(TPObject* tpobj, int len);
	void setCentroid(float xc, float yc, float zc);
	void setIndex(int i);
	int getIndex();
	float getLoad();
	long getHilbertID();
	int getNumTPs();
	TPObject* getTPs();
	float getx();
	float gety();
	float getz();
};

#endif
