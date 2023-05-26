/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef BARCODESTAT_H
#define BARCODESTAT_H

#include "util.h"
#include <unordered_map>
#include <string>
#include <fstream>
#include <time.h>
#include <iomanip>
#include "common.h"
#include "barcodeMap.h"

class BarcodeStat {
public:
	BarcodeStat(Options* opt);
	~BarcodeStat();
public:
	void barcodeStat();
private:
	Options* mOptions;
	string fastqFile;
public:
	BarcodeMap* barcodeMap;
	string mapOut;
	string logOut;
	uint64_t barcodeInt;
	uint16 barcodeCount;
	
};

#endif // !BARCODESTAT_H
