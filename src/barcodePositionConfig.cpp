/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "barcodePositionConfig.h"

BarcodePositionConfig::BarcodePositionConfig(Options* opt, int threadId)
{
	mOptions = opt;
	totalReads = 0;
	lowQReads = 0;
	dupReads = 0;
	withoutPositionRead = 0;
	mThreadId = threadId;
}

BarcodePositionConfig::~BarcodePositionConfig()
{
	if (!barcodeSet.empty()) {
		barcodeSet.clear();
		unordered_set<uint64_t>().swap(barcodeSet);
	}
	if (!bpmap.empty()) {
		bpmap.clear();
		unordered_map<uint64_t, Position1>().swap(bpmap);
	}
}

BarcodePositionConfig* BarcodePositionConfig::merge(BarcodePositionConfig** configs)
{
	if (sizeof(configs) == 0) {
		return nullptr;
	}
	BarcodePositionConfig* config = new BarcodePositionConfig(configs[0]->mOptions, 0);
	for (uint64_t i = 0; i < sizeof(configs) / sizeof(BarcodePositionConfig*); i++) {
		config->totalReads += configs[i]->totalReads;
		config->lowQReads += configs[i]->lowQReads;
		config->dupReads += configs[i]->dupReads;
		config->withoutPositionRead += configs[i]->withoutPositionRead;
		for (auto iter = configs[i]->bpmap.begin(); iter != configs[i]->bpmap.end(); iter++) {
			if (config->barcodeSet.count(iter->first) == 0) {
				config->barcodeSet.insert(iter->first);
				config->bpmap[iter->first] = iter->second;
			}
			else if(config->bpmap.count(iter->first)>0){
				config->bpmap.erase(iter->first);
			}
		}
		configs[i]->barcodeSet.clear();
		unordered_set<uint64_t>().swap(configs[i]->barcodeSet);
		configs[i]->bpmap.clear();
		unordered_map<uint64_t, Position1>().swap(configs[i]->bpmap);
	}
	return config;
}

void BarcodePositionConfig::addBarcode(uint64_t barcodeInt, Position1& pos)
{
	if (barcodeSet.count(barcodeInt) == 0) {
		barcodeSet.insert(barcodeInt);
		bpmap[barcodeInt] = pos;
	}
	else if (bpmap.count(barcodeInt) > 0) {
		dupReads++;
		bpmap.erase(barcodeInt);
	}
}

void BarcodePositionConfig::print()
{
	cout << "Total Reads:\t" << totalReads << endl;
	cout << "LowQ reads:\t" << lowQReads << endl;
	cout << "Duplicated reads:\t" << dupReads << endl;
	cout << "Barcode types:\t" << barcodeSet.size() << endl;
	cout << "Unique barcode types:\t" << bpmap.size() << endl;
}
