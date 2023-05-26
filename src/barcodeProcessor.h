/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef  BARCODE_PROCESSOR_H
#define BARCODE_PROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <chrono>
#include "read.h"
#include "util.h"
#include "barcodePositionMap.h"
#include "options.h"
#include "loadBarcodePositionMap.h"
#include "util.h"
//#include <boost/archive/binary_iarchive.hpp>

using namespace std;

class BarcodeProcessor {
public:
	BarcodeProcessor(Options* opt, unordered_map<uint64_t, Position1>* mbpmap);
	BarcodeProcessor(Options* opt, folly::F14ValueMap<uint64_t, Position1>* mbpmap);
  	BarcodeProcessor(Options* opt, folly::F14ValueMap<uint64_t, Position1>* mbpmap,BloomFilter* pBloomFilter);
//  	BarcodeProcessor(Options* opt, folly::F14ValueMap<uint64_t, uint64_t>* mbpmap,BloomFilter* pBloomFilter);
#if defined(OPT_MEMORY_PLAN_A) || defined(OPT_MEMORY_PLAN_B)
  	BarcodeProcessor(Options* opt, folly::F14ValueMap<uint64_t, Position1>* mbpmap,BloomFilter* pBloomFilter,unordered_map<uint64_t,int>* dnbCount);
#endif
	BarcodeProcessor(Options* opt, unordered_map<uint64_t, Position2>* mbpmap);
	BarcodeProcessor(Options* opt, folly::F14ValueMap<uint64_t, Position2>* mbpmap);
	BarcodeProcessor();
	~BarcodeProcessor();
	bool process(Read* read1, Read* read2);
	bool process(Read* read1,string tagMode,int Nindex);
    uint64_t getBarcodeSeq(string name,string& barcode,string& umi,string tagMode);
	void dumpDNBmap(string& dnbMapFile);

    bool barcodeStatAndFilter(string& barcodeQ);

private:
	void addPositionToName(Read* r, Position1* position, pair<string, string>* umi = NULL);
	void addPositionToName(Read* r, Position1* position, string umi);
  void addPositionToName(Read* r, Position2* position, pair<string, string>* umi = NULL);
  void addPositionToName(Read* r, Position2* position, string umi);
	void  getUMI(Read* r, pair<string, string>& umi, bool isRead2=false);
	void decodePosition(const uint32_t codePos, pair<uint16, uint16>& decodePos);
	void decodePosition(const uint64_t codePos, pair<uint32_t, uint32_t>& decodePos);
	uint32_t encodePosition(int fovCol, int fovRow);
	uint64_t encodePosition(uint32_t x, uint32_t y);
	long getBarcodeTypes();
	Position1* getPosition(uint64_t barcodeInt);
#ifdef OPT_MEMORY_PLAN_B
	Position1* getPosition(uint64_t barcodeInt,uint64_t& foundBarcode);
#endif
#ifdef OPT_MEMORY_PLAN_B
  Position1* getF14Position(uint64_t barcodeInt,uint64_t& foundBarcode);
  Position1* getF14Position(uint64_t barcodeInt,uint64_t &foundBarcode,int Nindex);
#else
  	Position1* getF14Position(uint64_t barcodeInt);
#endif
//  uint64_t getF14Position(uint64_t barcodeInt);
	Position1* getPosition(string& barcodeString);
  Position2* getPosition2(uint64_t barcodeInt);
  Position2* getF14Position2(uint64_t barcodeInt);
  Position2* getPosition2(string& barcodeString);
	void misMaskGenerate();
	string positionToString(Position1* position);
	string positionToString(Position2* position);
	string positionToString(Position* position);
	unordered_map<uint64_t, Position1>::iterator getMisOverlap(uint64_t barcodeInt);
#ifdef OPT_MEMORY_PLAN_B
	unordered_map<uint64_t, Position1>::iterator getMisOverlap(uint64_t barcodeInt,uint64_t& foundBarcode);
#endif
	unordered_map<uint64_t, Position2>::iterator getMisOverlap2(uint64_t barcodeInt);
//  uint64_t getF14MisOverlap(uint64_t barcodeInt);
#ifdef OPT_MEMORY_PLAN_B
  Position1* getF14MisOverlap(uint64_t barcodeInt,uint64_t& foundBarcode);
#else
  Position1* getF14MisOverlap(uint64_t barcodeInt);
#endif
//  uint64_t getF14MisOverlapWithBloomFiler(uint64_t barcodeInt);
#ifdef OPT_MEMORY_PLAN_B
  Position1 *getF14MisOverlapWithBloomFiler(uint64_t barcodeInt,uint64_t& foundBarcode);
#else
  Position1 *getF14MisOverlapWithBloomFiler(uint64_t barcodeInt);
#endif
	Position1* getNOverlap(string& barcodeString, uint8 Nindex);
	Position2* getF14MisOverlap2(uint64_t barcodeInt);
	Position2* getNOverlap2(string& barcodeString, uint8 Nindex);
	int getNindex(string& barcodeString);
	void addDNB(uint64_t barcodeInt);
	bool barcodeStatAndFilter(pair<string, string>& barcode);

    bool umiStatAndFilter (pair<string, string>& umi);
private:
	uint64_t* misMask;
	int misMaskLen;
	int* misMaskLens;
	const char q10 = '+';
	const char q20 = '5';
	const char q30 = '?';
	string polyT = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
	uint64_t polyTInt;
public:
	Options* mOptions;
	unordered_map<uint64_t, Position1>* bpmap;
#if defined(OPT_MEMORY_PLAN_A) || defined(OPT_MEMORY_PLAN_B)
//  folly::F14ValueMap<uint64_t,uint64_t>* f14bpmap;
  unordered_map<uint64_t,int>* dnbCount;
#endif
  	folly::F14ValueMap<uint64_t,Position1>* f14bpmap;

  unordered_map<uint64_t, Position2>* bpmap2;
  folly::F14ValueMap<uint64_t,Position2>* f14bpmap2;
  BloomFilter *bloomFilter;
//  set<uint64_t> checknames;
//  set<string> checknames2;
	long totalReads = 0;
	long mMapToSlideRead = 0;
	long overlapReads = 0;
	long overlapReadsWithMis = 0;
	long overlapReadsWithN = 0;
	long barcodeQ10 = 0;
	long barcodeQ20 = 0;
	long barcodeQ30 = 0;
	long umiQ10 = 0;
	long umiQ20 = 0;
	long umiQ30 = 0;
	long umiQ10FilterReads = 0;
	long umiNFilterReads = 0;
	long umiPloyAFilterReads = 0;
	unordered_map<uint64_t, int> mDNB;
	int mismatch;
	int barcodeLen;
  uint64_t queryTime=0;
    uint64_t getBarcodeInt(string name,string& umi);
  void getUmi(string name, string &umi);
//  Position1 *getPosition(uint64_t barcodeInt, string readName);
  
//  Position1 *getF14MisOverlapWithBloomFiler(unsigned long long int barcodeInt);
};


#endif // ! BARCODE_PROCESSOR_H
