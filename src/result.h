/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef  RESULT_H
#define RESULT_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <map>
#include <unordered_map>
#include <iomanip>
#include "common.h"
#include "options.h"
#include "writer.h"
#include "barcodeProcessor.h"
//#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>

using namespace std;

class Result {
public:
	Result(Options* opt, int threadId, bool paired = true);
	~Result();
	Writer* getWriter() { return mWriter; }
	static Result* merge(vector<Result*>& list);
	void print();
	void dumpDNBs(string& mappedDNBOutFile);
	void setBarcodeProcessor(unordered_map<uint64_t, Position1>* bpmap);
	void setBarcodeProcessor(folly::F14ValueMap<uint64_t, Position1>* bpmap);
	void setBarcodeProcessor(folly::F14ValueMap<uint64_t, Position1>* bpmap,BloomFilter* pbf);
#if defined(OPT_MEMORY_PLAN_A) || defined(OPT_MEMORY_PLAN_B)
	void setBarcodeProcessor(folly::F14ValueMap<uint64_t, Position1>* bpmap,BloomFilter* pbf,unordered_map<uint64_t,int>* dnbCount);
#endif
	void setBarcodeProcessor(unordered_map<uint64_t, Position2>* bpmap);
	void setBarcodeProcessor(folly::F14ValueMap<uint64_t, Position2>* bpmap);
private:
	void setBarcodeProcessor();
public:
	Options* mOptions;
	bool mPaired;
	long mTotalRead;
	long mFxiedFilterRead;
	long mDupRead;
	long mLowQuaRead;
	long mWithoutPositionReads;
	long overlapReadsWithMis = 0;
	long overlapReadsWithN = 0;
	//added by gc
	uint64_t seqQ10=0;
	uint64_t seqQ20=0;
	uint64_t seqQ30=0;
	uint64_t readLen=0;
	Writer* mWriter;
	BarcodeProcessor* mBarcodeProcessor;
	int mThreadId;
  	uint64_t tooManyNinBarcode=0;
	uint64_t totalBases=0;
	uint64_t containPolyA=0;
	uint64_t filterByPolyA=0;
};

#endif // ! RESULT_H
