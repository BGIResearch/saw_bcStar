/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef  BARCODEMAP_H
#define BARCODEMAP_H

#include "util.h"
#include <string>
#include <unordered_map>
#include <time.h>
#include <set>
#include "fastqreader.h"
#include "options.h"
#include "common.h"
//#include "heatMap.h"

class BarcodeMap {
public:
	BarcodeMap(int BarcodeStart, int BarcodeLength, long MapSize);
	~BarcodeMap();
	void fastqToBarcodeMap(string fastqFile, int segment = 1, bool isRC=false);
	void loadBarcodeMap(string mapFile);
	void loadBarcodeMapFromTxt(string mapFile);
	long dumpBarcodeMap(string mapFile);
	long dumpBarcodeMapToTxt(string mapFile);
public:
	int barcodeStart;
	int barcodeLength;
	long mapSize = 100000000;
	unordered_map<uint64_t, uint16> bmap;
	long totalReads = 0;
	long readsWithN = 0;
	long estDup = 0;
	uint64_t barcodeInt;
	uint16 bcount;
	string barcodeSeq;
};

#endif 

