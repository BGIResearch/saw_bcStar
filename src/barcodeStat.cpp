/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "barcodeStat.h"

BarcodeStat::BarcodeStat(Options* opt) 
{
	mOptions = opt;
	mapOut = opt->out;
	logOut = opt->out + ".stat";
}

BarcodeStat::~BarcodeStat()
{
}

void BarcodeStat::barcodeStat()
{
	barcodeMap = new BarcodeMap(mOptions->barcodeStart, mOptions->barcodeLen, mOptions->mapSize);
	barcodeMap->fastqToBarcodeMap(mOptions->in, mOptions->barcodeStat.segment, mOptions->barcodeStat.rc);
	long uniqReads = 0;
	if (ends_with(mapOut, ".bin"))
		uniqReads = barcodeMap->dumpBarcodeMap(mapOut);
	else
		uniqReads = barcodeMap->dumpBarcodeMapToTxt(mapOut);
	double totalReads = double(barcodeMap->totalReads);
	float lowQRate = double(barcodeMap->readsWithN) / (totalReads * mOptions->barcodeStat.segment) * 100;
	
	ofstream logWriter(logOut);
	logWriter << setiosflags(ios::fixed) << setprecision(2);
	logWriter << mOptions->in << endl;
	logWriter << "Reads:\t" << barcodeMap->totalReads << endl;
	logWriter << "barcodeSegement:\t" << mOptions->barcodeStat.segment << endl;
	logWriter << "LowQ:\t" << barcodeMap->readsWithN <<"\t" << lowQRate << "%" << endl;
	float uniqRate = double(uniqReads) / (totalReads * mOptions->barcodeStat.segment) * 100;
	long dupReads = totalReads * mOptions->barcodeStat.segment- barcodeMap->readsWithN - uniqReads;
	float dupRate = 100 - lowQRate - uniqRate;
	float estDupRate = barcodeMap->estDup / (totalReads * mOptions->barcodeStat.segment) * 100;
	logWriter << "Unique:\t" << uniqReads << "\t" << uniqRate << "%" << endl;
	logWriter << "Dup:\t" << dupReads << "\t" << dupRate << "%" << endl;
	logWriter << "BarcodeTypes:\t" << barcodeMap->bmap.size() << endl;
	logWriter << "EstNeiDup:\t" << barcodeMap->estDup << "\t" << estDupRate << "%" << endl;
	logWriter.close();
	delete barcodeMap;
}
