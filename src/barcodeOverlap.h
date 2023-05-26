/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef BARCODEOVERLAP_H
#define BARCODEOVERLAP_H

#include "util.h"
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <iomanip>
#include "fastqreader.h"
#include "options.h"
#include "common.h"
#include "barcodeMap.h"

class BarcodeOverlap {
public:
	BarcodeOverlap(Options* opt);
	~BarcodeOverlap();
	void overlapStat();
private:
	void barcodeMapDump(unordered_map<uint64_t, uint16>* barcodeMap, string outputFile);
	void misMaskGenerate(int mismatch);
	uint64_t getMisOverlap(uint64_t barodeInt);
	void logWrite();
	void setToArray(set<uint64_t>& mset, uint64_t* marray);
public:
	Options* mOptions;
	string map1File;
	string map2File;
	string overlapMapOut;
	string logOut;
	BarcodeMap* barcodeMap;
	unordered_map<uint64_t, uint16> overlapMap;
	unordered_map<uint64_t, uint16>misOverlapMap;
	uint64_t* baseMisCount;
	uint64_t baseTransCount[4][4] = {{0}};
	uint64_t* misMask;
	int misMaskLen = 0;
	int mismatch = 0;
	long totalReads = 0;
	long overlapReads = 0;
	long misOverlapReads = 0;
	long overlapReadsUniq = 0;
	long misOverlapReadsUniq = 0;
	long overlapBarcode = 0;
	long misOverlapBarcode = 0;
	long overlapBarcodeUniq = 0;
	long misOverlapBarcodeUniq = 0;
	long totalBarcode = 0;
	long dupMisOverlapBarcode = 0;

};

#endif // !BARCODEOVERLAP_H
