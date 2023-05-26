/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef  BARCODE_POSITION_CONFIG_H
#define BARCODE_POSITION_CONFIG_H

#include <stdio.h>
#include <unordered_map>
#include <unordered_set>
#include "common.h"
#include "options.h"
#include "barcodePositionMap.h"

using namespace std;

class BarcodePositionConfig {
public:
	BarcodePositionConfig(Options* opt, int threadId);
	~BarcodePositionConfig();
	static BarcodePositionConfig* merge(BarcodePositionConfig** configs);
	void addBarcode(uint64_t barcodeInt, Position1& pos);
	void print();
public:
	Options* mOptions;
	long totalReads;
	long lowQReads;
	long dupReads;
	long withoutPositionRead;
	int mThreadId;
	unordered_set<uint64_t> barcodeSet;
	unordered_map<uint64_t, Position1> bpmap;
};

#endif 
