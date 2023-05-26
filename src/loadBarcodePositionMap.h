/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef LOAD_BARCODE_POSITION_MAP_H
#define LOAD_BARCODE_POSITION_MAP_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <atomic>
#include <unordered_set>
#include <functional>
#include <cstring>
#include <unistd.h>
#include "options.h"
#include "read.h"
#include "fastqreader.h"
#include "barcodePositionMap.h"
#include "slideMask.h"
#include "barcodePositionConfig.h"

using namespace std;

typedef struct ReadPack {
	Read** data;
	int count;
} ReadPack;

typedef struct ReadRepository {
	ReadPack** packBuffer;
	atomic_long readPos;
	atomic_long writePos;
}ReadRepository;

class LoadBarcodePositionMap {
public:
	// try mutiple thread for load barcodePositionMap, but faild
	LoadBarcodePositionMap(Options* opt);
	~LoadBarcodePositionMap();
	BarcodePositionMap* load();
private:
	bool processSingleEnd(ReadPack* pack, BarcodePositionConfig* config);
	void initPackRepository();
	void destoryPackRepository();
	void producePack(ReadPack* pack);
	void consumePack(BarcodePositionConfig* config);
	void producerTask();
	void consumerTask(BarcodePositionConfig* config);
public:
	BarcodePositionMap* barcodePositionMap;
	SlideMask* sm;
private:
	Options* mOptions;
	ReadRepository mRepo;
	atomic_bool mProduceFinished;
	atomic_int mFinishedThreads;
	std::mutex mInputMtx;
	gzFile mZipFile;
	int barcodeStart;
	int barcodeLen;
	bool rc;
	const string SEQ500 = "SEQ500";
	bool isSeq500 = false;
};

#endif