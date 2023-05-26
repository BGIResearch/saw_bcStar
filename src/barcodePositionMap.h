/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef  BARCODEPOSITIONMAP_H
#define BARCODEPOSITIONMAP_H

#include <string>
#include <cstring>
#include "fastqreader.h"
#include "read.h"
#include "util.h"
#include "slideMask.h"
#include "options.h"
#include <folly/container/F14Map.h>
//#include "heatMap.h"
#include "chipMaskHDF5.h"
#include <unordered_map>
#include <iomanip>
#include <set>

//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>
//#include <boost/serialization/string.hpp>
//#include <boost/serialization/unordered_map.hpp>
//#include <opencv2/opencv.hpp>

using namespace std;

class BarcodePositionMap {
public:
	//BarcodePositionMap(vector<string>& InFile, string MaskFile, int BarcodeStart, int BarcodeLen, int Segment, int TurnFovDegree, bool IsSeq500, int RC, string readidSep = "/");
	//BarcodePositionMap(string InFile, string MaskFile, int BarcodeStart, int BarcodeLen, int Segment, int TurnFovDegree, bool IsSeq500, int RC);
	//BarcodePositionMap(string InFile, int BarcodeStart, int BarcodeLen);
	BarcodePositionMap(Options* opt);
	~BarcodePositionMap();
private:
	void initiate();
	void getSuffixLen();
	int readQualityStat(string& readQ, int index);
	bool barcodeFilter(string& readSeq, int index);
	bool isEst(Position1& pos1, Position1& pos2);
public:
	long getBarcodeTypes();
	void dumpbpmap(string& mapOutFile);
	void loadbpmap();
	unordered_map<uint64_t, Position1>* getBpmap() { return &bpmap; };
	uint64_t touint64_t(Position1 pos);
	Position1 toPosition1(uint64_t);
	bool isNeighbor(uint64_t pos1,uint64_t pos2);
	bool isNeighbor(set<uint64_t> manyPos);
public:
	unordered_map<uint64_t,  Position1> bpmap;
    unordered_map<uint32_t,  Position1>* splitBpmap;
#if defined(OPT_MEMORY_PLAN_A) || defined(OPT_MEMORY_PLAN_B)
//  	folly::F14ValueMap<uint64_t,uint64_t> f14bpmap;
  	unordered_map<uint64_t,int> dnbCount;
#endif
  	folly::F14ValueMap<uint64_t,Position1> f14bpmap;

	unordered_map<uint64_t,  Position2> bpmap2;
//	unordered_map<uint32_t,  Position2>* splitBpmap;
//	folly::F14ValueMap<uint64_t,Position2> f14bpmap2;
	BloomFilter *bloomFilter;
    int splitNum;
	long* totalReads;
	long* readsWithN;
	long* dupReads;
	long* ESTdupReads;
	long* readsWithoutPos;
	long** polyReads;
	long* readsQ10;
	long* readsQ20;
	long* readsQ30;
	long* totalBase;
	long* totalBarcodes;
	int inFastqNumber;
	vector<string> inFile;
	Options* mOptions;
	set<uint64_t> dupBarcode;
	string maskFile;
	int barcodeStart;
	int barcodeLen;
	int segment;
	int turnFovDegree;
	bool isSeq500;	
	string firstFastq;
	string readidSep;
	int suffixLen = 0;
	int readLen;
	int indexLen = 8;
	//rc==0 only get forward barcode, rc==1 only get reverse complement barcode, rc==2 get both forward and reverse complement barcode
	int rc;
	slideRange sliderange;
	SlideMask* sm;
	//uint64_t polyTint;
    void loadbpmap(string barcodePositionMapFile);
};

#endif // ! BARCODEPOSITIONMAP_H