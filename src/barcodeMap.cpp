/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "barcodeMap.h"

BarcodeMap::BarcodeMap(int BarcodeStart, int BarcodeLength, long MapSize)
{
	barcodeStart = BarcodeStart;
	barcodeLength = BarcodeLength;
	if (MapSize > 0) {
		mapSize = MapSize;
		bmap.rehash(mapSize);
	}
}

BarcodeMap::~BarcodeMap()
{
	bmap.clear();
}

void BarcodeMap::fastqToBarcodeMap(string fastqFile, int segment, bool isRC)
{
	string* readSeq;
	string barcodeSeq;
	uint64_t oldBarcodeInt = MAX_BARCODE;
	set<uint64_t> barcodeSet;
	FastqReader fastqReader(fastqFile);
	cout << "fastq read begin..." << endl;
	clock_t start;
	start = time(NULL);
	while (true) {
		Read* read = fastqReader.read();
		if (!read) {
			delete read;
			break;
		}
		totalReads++;
		readSeq = &read->mSeq.mStr;
		barcodeSet.clear();
		int barcodeWithN = 0;
		for (int i = 0; i < segment; i++) {
			barcodeSeq = readSeq->substr(barcodeStart + barcodeLength * i, barcodeLength);
		  string::size_type Npos = barcodeSeq.find("N");
			if (Npos != readSeq->npos) {
				readsWithN++;
				barcodeWithN++;
				continue;
			}
			barcodeInt = seqEncode(barcodeSeq.c_str(), 0 , barcodeLength, isRC);
			barcodeSet.insert(barcodeInt);
			if (segment == 1 && barcodeInt == oldBarcodeInt)
				estDup++;
			oldBarcodeInt = barcodeInt;
		}
		if (barcodeSet.size() == 0) {
			delete read;
			continue;
		}
		for (auto iter = barcodeSet.begin(); iter != barcodeSet.end(); iter++) {
			if (bmap.count(*iter) == 0) {
				bmap[*iter] = 1 ;
			}
			else {
				bmap[*iter]++;
			}
		}
		if (segment > 1 && (int)barcodeSet.size() < segment-barcodeWithN) {
			estDup++;
		}
		delete read;
	}
	cout << "barcode load time: " << (time(NULL) - start)<< " Sec" << endl;
}


void BarcodeMap::loadBarcodeMap(string mapFile)
{
	cout << "map file read begin..." << endl;
	time_t start = time(NULL);
	ifstream mapFileStream;
	mapFileStream.open(mapFile, ios::binary);
	while (true) {		
		mapFileStream.read((char*)&barcodeInt, sizeof(uint64_t));
		mapFileStream.read((char*)&bcount, sizeof(uint16));
		if (mapFileStream.eof())
			break;
		bmap[barcodeInt] = bcount;
		totalReads += bcount;
	}
	mapFileStream.close();
	cout << "barcode load time: " << (time(NULL)-start) << " Sec" << endl;
}

void BarcodeMap::loadBarcodeMapFromTxt(string mapFile)
{
	cout << "map file readbegin..." << endl;
	time_t start = time(NULL);
	ifstream mapFileStream;
	mapFileStream.open(mapFile);
	string line = "";
	int index = 0;
	while (true) {
		getline(mapFileStream, line);
		if (mapFileStream.eof())
			break;
		index = line.find("\t");
		if (starts_with(line, "#") || index == -1)
			continue;
		barcodeSeq = line.substr(0, barcodeLength);
		if (barcodeSeq.find("N") != barcodeSeq.npos && barcodeSeq.find("N") != barcodeSeq.rfind("N"))
			continue;
		barcodeInt = seqEncode(barcodeSeq.c_str(), 0 , barcodeLength);
		bcount = (uint16)stoi(line.substr(index + 1, line.size()));
		bmap[barcodeInt] = bcount;
	}
	mapFileStream.close();
	cout << "barcode load time: " << (time(NULL)-start) << " Sec" <<endl;

}

long BarcodeMap::dumpBarcodeMap(string mapFile)
{
	long uniqBarcode = 0;
	cout << "write map to bin file begin..." << endl;
	time_t start = time(NULL);
	ofstream binFileStream(mapFile, ios::out | ios::binary);
	unordered_map<uint64_t, uint16>::iterator mapIter;
	mapIter = bmap.begin();
	while (mapIter != bmap.end()) {
		binFileStream.write((char*)&mapIter->first, sizeof(uint64_t));
		binFileStream.write((char*)&mapIter->second, sizeof(uint16));
		if (mapIter->second == 1)
			uniqBarcode++;
		mapIter++;
	}
	binFileStream.close();
	cout << "barcode dump time: " << (time(NULL)-start) << " Sec" << endl;
	return uniqBarcode;
}

long BarcodeMap::dumpBarcodeMapToTxt(string mapFile)
{
	long uniqBarcode = 0;
	cout <<"write map to txtfile begin..." << endl;
	time_t start = time(NULL);
	ofstream mapOutWriter(mapFile);
	unordered_map<uint64_t, uint16>::iterator mapIter;
	mapIter = bmap.begin();
	while (mapIter != bmap.end()) {
		string barcodeString = seqDecode(mapIter->first, barcodeLength);
		mapOutWriter << barcodeString << "\t" << mapIter->second << endl;
		if (mapIter->second == 1)
			uniqBarcode++;
		mapIter++;
	}
	mapOutWriter.close();
	cout << "barcode dump time: " << (time(NULL) - start)  << "Sec" << endl;
	return uniqBarcode;
}
