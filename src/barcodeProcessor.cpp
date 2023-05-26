/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#include "barcodeProcessor.h"
BarcodeProcessor::BarcodeProcessor(Options *opt, unordered_map<uint64_t, Position1> *mbpmap)
{
	mOptions = opt;
	bpmap = mbpmap;
	mismatch = opt->transBarcodeToPos.mismatch;
	barcodeLen = opt->barcodeLen;
	polyTInt = seqEncode(polyT.c_str(), 0, barcodeLen, mOptions->rc);
	misMaskGenerate();
}
BarcodeProcessor::BarcodeProcessor(Options *opt, folly::F14ValueMap<uint64_t, Position1> *mbpmap)
{
	mOptions = opt;
	f14bpmap = mbpmap;
	mismatch = opt->transBarcodeToPos.mismatch;
	barcodeLen = opt->barcodeLen;
	polyTInt = seqEncode(polyT.c_str(), 0, barcodeLen, mOptions->rc);
	misMaskGenerate();
}
BarcodeProcessor::BarcodeProcessor(Options *opt, folly::F14ValueMap<uint64_t, Position1> *mbpmap, BloomFilter *pBloomFilter)
{
	mOptions = opt;
	f14bpmap = mbpmap;
	bloomFilter = pBloomFilter;
	mismatch = opt->transBarcodeToPos.mismatch;
	barcodeLen = opt->barcodeLen;
	polyTInt = seqEncode(polyT.c_str(), 0, barcodeLen, mOptions->rc);
	misMaskGenerate();
}
// BarcodeProcessor::BarcodeProcessor(Options* opt, folly::F14ValueMap<uint64_t, uint64_t>* mbpmap,BloomFilter* pBloomFilter)
//{
//   mOptions = opt;
//   f14bpmap = mbpmap;
//   bloomFilter=pBloomFilter;
//   mismatch = opt->transBarcodeToPos.mismatch;
//   barcodeLen = opt->barcodeLen;
//   polyTInt = seqEncode(polyT.c_str(), 0, barcodeLen, mOptions->rc);
//   misMaskGenerate();
// }
#if defined(OPT_MEMORY_PLAN_A) || defined(OPT_MEMORY_PLAN_B)
BarcodeProcessor::BarcodeProcessor(Options *opt, folly::F14ValueMap<uint64_t, Position1> *mbpmap, BloomFilter *pBloomFilter, unordered_map<uint64_t, int> *dnbCount)
{
	mOptions = opt;
	f14bpmap = mbpmap;
	this->dnbCount = dnbCount;
	bloomFilter = pBloomFilter;
	mismatch = opt->transBarcodeToPos.mismatch;
	barcodeLen = opt->barcodeLen;
	polyTInt = seqEncode(polyT.c_str(), 0, barcodeLen, mOptions->rc);
	misMaskGenerate();
}
#endif
BarcodeProcessor::BarcodeProcessor(Options *opt, unordered_map<uint64_t, Position2> *mbpmap)
{
	mOptions = opt;
	bpmap2 = mbpmap;
	mismatch = opt->transBarcodeToPos.mismatch;
	barcodeLen = opt->barcodeLen;
	polyTInt = seqEncode(polyT.c_str(), 0, barcodeLen, mOptions->rc);
	misMaskGenerate();
}
BarcodeProcessor::BarcodeProcessor(Options *opt, folly::F14ValueMap<uint64_t, Position2> *mbpmap)
{
	mOptions = opt;
	f14bpmap2 = mbpmap;
	mismatch = opt->transBarcodeToPos.mismatch;
	barcodeLen = opt->barcodeLen;
	polyTInt = seqEncode(polyT.c_str(), 0, barcodeLen, mOptions->rc);
	misMaskGenerate();
}

BarcodeProcessor::BarcodeProcessor()
{
}

BarcodeProcessor::~BarcodeProcessor()
{
}

bool BarcodeProcessor::process(Read *read1, Read *read2)
{
	totalReads++;
	string barcode = read1->mSeq.mStr.substr(mOptions->barcodeStart, mOptions->barcodeLen);
	string barcodeQ = read1->mQuality.substr(mOptions->barcodeStart, mOptions->barcodeLen);
	barcodeStatAndFilter(barcodeQ);
	Position1 *position = getPosition(barcode);
	pair<string, string> umi;
	bool umiPassFilter = true;
	if (mOptions->transBarcodeToPos.umiStart >= 0)
	{
		if (mOptions->transBarcodeToPos.umiRead == 1)
		{
			getUMI(read1, umi);
		}
		else
		{
			getUMI(read2, umi, true);
		}
		umiPassFilter = umiStatAndFilter(umi);
	}

	if (position != nullptr)
	{
		if (mOptions->transBarcodeToPos.umiStart >= 0)
		{
			addPositionToName(read2, position, &umi);
		}
		else
		{
			addPositionToName(read2, position);
		}
		if (!mOptions->transBarcodeToPos.mappedDNBOutFile.empty())
			addDNB(encodePosition((uint32_t)position->x, (uint32_t)position->y));
		mMapToSlideRead++;
		if (umiPassFilter)
			return true;
	}
	return false;
}
// added by gc
// processed fq from writeFq
bool BarcodeProcessor::process(Read *read1, string tagMode, int Nindex)
{
	//  leaktracer_startMonitoringThisThread();
	if (read1 == nullptr)
	{
		return false;
	}
	totalReads++;
	//    string barcode = read1->mSeq.mStr.substr(mOptions->barcodeStart, mOptions->barcodeLen);
	string barcode, umi;

	uint64_t barcodeInt = getBarcodeSeq(read1->mName, barcode, umi, tagMode);
	uint64_t foundBarcode = barcodeInt;
	//    cout<<barcodeInt<<endl;
	//    string barcodeQ = read1->mQuality.substr(mOptions->barcodeStart, mOptions->barcodeLen);
	//    barcodeStatAndFilter(barcodeQ);
	//    Position1* position = getPosition(barcodeInt,read1->mName);

	Position1 *position = nullptr;
	auto tStart = chrono::high_resolution_clock::now();
	if (Nindex == -1)
	{
		position = getF14Position(barcodeInt, foundBarcode);
	}
	else
	{
		position = getF14Position(barcodeInt, foundBarcode, Nindex);
	}
	//	  position = getF14Position(barcodeInt,foundBarcode);
	auto tEnd = chrono::high_resolution_clock::now();
	queryTime += chrono::duration<double, nano>(tEnd - tStart).count();
	//    pair<string, string> umi;
	//    bool umiPassFilter = true;
	//    if (mOptions->transBarcodeToPos.umiStart >= 0) {
	//        if (mOptions->transBarcodeToPos.umiRead == 1) {
	//            getUMI(read1, umi);
	//        }
	//        else {
	//            getUMI(read2, umi, true);
	//        }
	//        umiPassFilter = umiStatAndFilter(umi);
	//    }

	if (position != nullptr)
	{
		//	  	cout<<barcodeInt<<endl;
		if (mOptions->transBarcodeToPos.umiStart >= 0)
		{
			addPositionToName(read1, position, umi);
		}
		else
		{
			addPositionToName(read1, position);
		}
		if (!mOptions->transBarcodeToPos.mappedDNBOutFile.empty())
		{
#ifdef OPT_MEMORY_PLAN_A
			//		  coor=((position->x)>>coorBits) | position->y;
			// set 64th bit as a overflow flag,63-39 bits as count, 38-1 bits as position
			//		posLock[position->y&0b111].lock();
			uint64_t coor = encodePosition((uint32_t)position->x, (uint32_t)position->y);
			posLock.lock();
			if (dnbCount->find(coor) != dnbCount->end())
			{
				(*dnbCount)[coor]++;
			}
			else
			{
				(*dnbCount)[coor] = 1;
			}
			//		if(coor>>63){
			//		  totalDnb[encodePosition((uint32_t)position->x, (uint32_t)position->y)]+=1;
			//		}else{
			//		  uint64_t curCount=(coor & countMask)>>countBits;
			//		  curCount++;
			//		  if((curCount>>25)==1) {
			//			coor |= flagMask;
			//			totalDnb[encodePosition((uint32_t)position->x, (uint32_t)position->y)]=(1<<25);
			//		  }else {
			//			coor = (coor & noCountMask) | (curCount << countBits);
			//		  }
			//		  (&f14bpmap)[barcodeInt]=coor;
			//		}
			posLock.unlock();
#endif
#ifdef OPT_MEMORY_PLAN_B

			uint64_t coor = encodePosition((uint32_t)position->x, (uint32_t)position->y);
			posLock.lock();
			if (mOptions->useF14)
			{
				if (f14bpmap->find(foundBarcode) == f14bpmap->end())
				{
					string errMsg = "code error";
					sawErrCode(err_sw_exception, errMsg);
					cerr << "Error, code error" << endl;
					exit(1);
				}
				if ((*f14bpmap)[foundBarcode].fullTag == 0)
				{
					(*f14bpmap)[foundBarcode].count++;
					//		  cout<<foundBarcode<<"\t"<<(*f14bpmap)[foundBarcode].x<<"\t"<<(*f14bpmap)[foundBarcode].y<<"\t"<<(*f14bpmap)[foundBarcode].count<<endl;
					if ((*f14bpmap)[foundBarcode].count >= (1 << 25))
					{
						(*f14bpmap)[foundBarcode].fullTag = 1;
						(*dnbCount)[coor] = 1 >> 25;
					}
				}
				else
				{
					(*dnbCount)[coor]++;
				}
			}
			else
			{
				if (bpmap->find(foundBarcode) == bpmap->end())
				{
					string errMsg = "code error";
					sawErrCode(err_sw_exception, errMsg);
					cerr << "Error, code error" << endl;
					exit(1);
				}
				if ((*bpmap)[foundBarcode].fullTag == 0)
				{
					(*bpmap)[foundBarcode].count++;
					//		  cout<<foundBarcode<<"\t"<<(*bpmap)[foundBarcode].x<<"\t"<<(*bpmap)[foundBarcode].y<<"\t"<<(*bpmap)[foundBarcode].count<<endl;
					if ((*bpmap)[foundBarcode].count >= (1 << 25))
					{
						(*bpmap)[foundBarcode].fullTag = 1;
						(*dnbCount)[coor] = 1 >> 25;
					}
				}
				else
				{
					(*dnbCount)[coor]++;
				}
			}
			posLock.unlock();
#else
			addDNB(encodePosition(position->x, position->y));
#endif
		}
		mMapToSlideRead++;
		//        delete position;
		return true;
		//        if (umiPassFilter)
		//            return true;
	}
	else
	{
		//        cout<<read1->mName<<endl;
	}
	//  leaktracer_stopAllMonitoring();
	//  leaktracer_writeLeaksToFile("leak.out");
	//  exit(1);
	return false;
}
uint64_t BarcodeProcessor::getBarcodeSeq(string name, string &barcode, string &umi, string tagMode)
{
	vector<string> eles;
	split(name, eles, " ");
	if (eles.size() != 3)
	{
		string errMsg = "barcode format error";
		sawErrCode(err_sw_exception, errMsg);
		cerr << "Error,barcode format error" << endl;
		exit(1);
	}
	char *value;
	long barcodeInt = strtol(eles[1].c_str(), &value, 16);
	long umiInt = strtol(eles[2].c_str(), &value, 16);
	for (int i = 0; i < mOptions->barcodeLen; i++)
	{
		int n = (barcodeInt >> (2 * (mOptions->barcodeLen - i - 1))) & 3;
		barcode += ATCG_BASES2[n];
	}
	for (int i = 0; i < mOptions->transBarcodeToPos.umiLen; i++)
	{
		int n = (umiInt >> (2 * (mOptions->transBarcodeToPos.umiLen - i - 1))) & 3;
		umi += ATCG_BASES2[n];
	}
	// change umi to 0x or not
	if (tagMode == "spatial")
	{
		umi = eles[2];
	}
	return barcodeInt;
}
void BarcodeProcessor::getUmi(string name, string &umi)
{
	vector<string> eles;
	split(name, eles, " ");
	if (eles.size() != 3)
	{
		string errMsg = "barcode format error";
		sawErrCode(err_sw_exception, errMsg);
		cerr << "Error,barcode format error" << endl;
		exit(1);
	}
	char *value;
	//    int barcodeInt = 0;
	//  long barcodeInt=strtol(eles[1].c_str(), &value, 16);
	//    sscanf(eles[1].c_str(), "%x", &barcodeInt);
	//    int umiInt=0;
	//    sscanf(eles[2].c_str(), "%x", &umiInt);
	long umiInt = strtol(eles[2].c_str(), &value, 16);
	//    long barcodeInt=strtol(eles[1].c_str(),&value,16);
	//    long umiInt=strtol(eles[2].c_str(),&value,16);
	for (int i = 0; i < mOptions->transBarcodeToPos.umiLen; i++)
	{
		int n = (umiInt >> (2 * (mOptions->transBarcodeToPos.umiLen - i - 1))) & 3;
		umi += ATCG_BASES2[n];
	}
	//  vector<string>().swap(eles);
	//    string returnV(value);
}
uint64_t BarcodeProcessor::getBarcodeInt(string name, string &umi)
{
	vector<string> eles;
	split(name, eles, " ");
	if (eles.size() != 3)
	{
		string errMsg = "barcode format error";
		sawErrCode(err_sw_exception, errMsg);
		cerr << "Error,barcode format error" << endl;
		exit(1);
	}
	char *value;
	//    int barcodeInt = 0;
	long barcodeInt = strtol(eles[1].c_str(), &value, 16);
	long umiInt = strtol(eles[2].c_str(), &value, 16);
	for (int i = 0; i < mOptions->transBarcodeToPos.umiLen; i++)
	{
		int n = (umiInt >> (2 * (mOptions->transBarcodeToPos.umiLen - i - 1))) & 3;
		umi += ATCG_BASES2[n];
	}
	return barcodeInt;
	//    string returnV(value);
}
void BarcodeProcessor::addPositionToName(Read *r, Position1 *position, pair<string, string> *umi)
{
	string position_tag = positionToString(position);
	string::size_type readTagPos = r->mName.find("/");
	string readName;
	if (readTagPos != string::npos)
	{
		readName = r->mName.substr(0, readTagPos);
	}
	else
	{
		readName = r->mName;
	}
	if (umi == NULL)
		r->mName = readName + "|||CB:Z:" + position_tag;
	else
	{
		r->mName = readName + "|||CB:Z:" + position_tag + "|||UR:Z:" + umi->first + "|||UY:Z:" + umi->second;
	}
}
void BarcodeProcessor::addPositionToName(Read *r, Position1 *position, string umi)
{
	string position_tag = positionToString(position);
	//    int readTagPos = r->mName.find("/");
	uint64_t readTagPos = r->mName.find_first_of(" ");
	string readName;
	if (readTagPos != string::npos)
	{
		readName = r->mName.substr(0, readTagPos);
	}
	else
	{
		readName = r->mName;
	}
	r->mName = readName + "|||CB:Z:" + position_tag + "|||UR:Z:" + umi;
}
void BarcodeProcessor::addPositionToName(Read *r, Position2 *position, pair<string, string> *umi)
{
	string position_tag = positionToString(position);
	string::size_type readTagPos = r->mName.find("/");
	string readName;
	if (readTagPos != string::npos)
	{
		readName = r->mName.substr(0, readTagPos);
	}
	else
	{
		readName = r->mName;
	}
	if (umi == NULL)
		r->mName = readName + "|||CB:Z:" + position_tag;
	else
	{
		r->mName = readName + "|||CB:Z:" + position_tag + "|||UR:Z:" + umi->first + "|||UY:Z:" + umi->second;
	}
}
void BarcodeProcessor::addPositionToName(Read *r, Position2 *position, string umi)
{
	string position_tag = positionToString(position);
	//    int readTagPos = r->mName.find("/");
	uint64_t readTagPos = r->mName.find_first_of(" ");
	string readName;
	if (readTagPos != string::npos)
	{
		readName = r->mName.substr(0, readTagPos);
	}
	else
	{
		readName = r->mName;
	}
	r->mName = readName + "|||CB:Z:" + position_tag + "|||UR:Z:" + umi;
}
void BarcodeProcessor::getUMI(Read *r, pair<string, string> &umi, bool isRead2)
{
	string umiSeq = r->mSeq.mStr.substr(mOptions->transBarcodeToPos.umiStart, mOptions->transBarcodeToPos.umiLen);
	string umiQ = r->mQuality.substr(mOptions->transBarcodeToPos.umiStart, mOptions->transBarcodeToPos.umiLen);
	umi.first = umiSeq;
	umi.second = umiQ;
	if (isRead2)
	{
		r->mSeq.mStr = r->mSeq.mStr.substr(0, mOptions->transBarcodeToPos.umiStart);
		r->mQuality = r->mQuality.substr(0, mOptions->transBarcodeToPos.umiStart);
	}
}

void BarcodeProcessor::decodePosition(const uint32_t codePos, pair<uint16, uint16> &decodePos)
{
	decodePos.first = codePos >> 16;
	decodePos.second = codePos & 0x0000FFFF;
}

void BarcodeProcessor::decodePosition(const uint64_t codePos, pair<uint32_t, uint32_t> &decodePos)
{
	decodePos.first = codePos >> 32;
	decodePos.second = codePos & 0x00000000FFFFFFFF;
}

uint32_t BarcodeProcessor::encodePosition(int fovCol, int fovRow)
{
	uint32_t encodePos = (fovCol << 16) | fovRow;
	return encodePos;
}

uint64_t BarcodeProcessor::encodePosition(uint32_t x, uint32_t y)
{
	uint64_t encodePos = ((uint64_t)x << 32) | (uint64_t)y;
	return encodePos;
}

long BarcodeProcessor::getBarcodeTypes()
{
	return bpmap->size();
}

Position1 *BarcodeProcessor::getPosition(uint64_t barcodeInt)
{
	unordered_map<uint64_t, Position1>::iterator iter = bpmap->find(barcodeInt);
	if (iter != bpmap->end())
	{
		overlapReads++;
		return &iter->second;
	}
	else if (mismatch > 0)
	{
		iter = getMisOverlap(barcodeInt);
		if (iter != bpmap->end())
		{
			overlapReadsWithMis++;
			return &iter->second;
		}
		else
		{
			return nullptr;
		}
	}
	return nullptr;
}
#ifdef OPT_MEMORY_PLAN_B
Position1 *BarcodeProcessor::getPosition(uint64_t barcodeInt, uint64_t &foundBarcode)
{
	unordered_map<uint64_t, Position1>::iterator iter = bpmap->find(barcodeInt);
	if (iter != bpmap->end())
	{
		overlapReads++;
		foundBarcode = barcodeInt;
		return &iter->second;
	}
	else if (mismatch > 0)
	{
		iter = getMisOverlap(barcodeInt, foundBarcode);
		if (iter != bpmap->end())
		{
			overlapReadsWithMis++;
			return &iter->second;
		}
		else
		{
			return nullptr;
		}
	}
	return nullptr;
}
#endif
#ifdef OPT_MEMORY_PLAN_B
Position1 *BarcodeProcessor::getF14Position(uint64_t barcodeInt, uint64_t &foundBarcode, int Nindex)
{
	if (Nindex == -1)
	{
		cerr << "Error, code error." << __FILE__ << ":" << __LINE__ << endl;
		exit(1);
	}
	// N has the same encode (11) with G
	int misCount = 0;
	auto iter = f14bpmap->find(barcodeInt);
	auto overlapIter = f14bpmap->end();
	if (iter != f14bpmap->end())
	{
		misCount++;
		overlapIter = iter;
		foundBarcode = barcodeInt;
	}
	for (uint64_t j = 1; j < 4; j++)
	{
		uint64_t misBarcodeInt = barcodeInt ^ (j << (Nindex * 2));
		iter = f14bpmap->find(misBarcodeInt);
		if (iter != f14bpmap->end())
		{
			misCount++;
			if (misCount > 1)
			{
				return nullptr;
			}
			foundBarcode = misBarcodeInt;
			overlapIter = iter;
		}
	}
	if (misCount == 1)
	{
		overlapReadsWithN++;
		return &overlapIter->second;
	}
	return nullptr;
}
Position1 *BarcodeProcessor::getF14Position(uint64_t barcodeInt, uint64_t &foundBarcode)
{
	auto iter = f14bpmap->find(barcodeInt);
	if (iter != f14bpmap->end())
	{
		foundBarcode = barcodeInt;
		overlapReads++;
		return &iter->second;
		//	Position1* rPos=nullptr;
		//	rPos->x=((iter->second) & coorXMask)>>coorBits;
		//	rPos->y=(iter->second) & coorMask;
		//	return rPos;
	}
	else if (mismatch > 0)
	{
		Position1 *rv;
		if (mismatch == 1)
		{
#ifdef LOADH5_OPENMP
#ifdef OPT_MEMORY_IN_H5
			rv = getF14MisOverlap(barcodeInt, foundBarcode);
#else
			if (mOptions->useBf)
			{
				rv = getF14MisOverlapWithBloomFiler(barcodeInt, foundBarcode);
			}
			else
			{
				rv = getF14MisOverlap(barcodeInt, foundBarcode);
			}
#endif
#else
			rv = getF14MisOverlap(barcodeInt);
#endif
		}
		else
		{
			rv = getF14MisOverlap(barcodeInt, foundBarcode);
		}
		if (rv != nullptr)
		{
			overlapReadsWithMis++;
			return rv;
		}
		else
		{
#ifdef TESTLRMOVE1
			//	  move 1 bp on left. Suppose bc length is 25
			for (uint64_t i = 0; i < 4; i++)
			{
				uint64_t newBc = (barcodeInt << 2) & 0x3ffffffffffff;
				newBc = newBc | i;
				//		cout << "move 1 bp on right\t"<<seqDecode(barcodeInt,25)<<"\t"<<seqDecode(newBc,25)<<endl;
				if (f14bpmap->find(newBc) != f14bpmap->end())
				{
					foundBarcode = newBc;
					posLock.lock();
					cout << "move 1 bp on right\t" << seqDecode(barcodeInt, 25) << "\t" << seqDecode(newBc, 25) << endl;
					posLock.unlock();
					return &(f14bpmap->find(newBc)->second);
				}
			}
			//	  move 1 bp on right. Suppose bc length is 25
			for (uint64_t i = 0; i < 4; i++)
			{
				uint64_t newBc = barcodeInt >> 2;
				//		uint64_t firstBp=i<<48;
				newBc = newBc | (i << 48);
				//		cout << "move 1 bp on left\t"<<seqDecode(barcodeInt,25)<<"\t"<<seqDecode(newBc,25)<<endl;

				if (f14bpmap->find(newBc) != f14bpmap->end())
				{
					foundBarcode = newBc;
					posLock.lock();
					cout << "move 1 bp on left\t" << seqDecode(barcodeInt, 25) << "\t" << seqDecode(newBc, 25) << endl;
					posLock.unlock();
					return &(f14bpmap->find(newBc)->second);
				}
			}
			//	  exit(1);
			return nullptr;
#else

			return nullptr;
#endif
		}
	}
	return nullptr;
}
#else
Position1 *BarcodeProcessor::getF14Position(uint64_t barcodeInt)
{
	auto iter = f14bpmap->find(barcodeInt);
	if (iter != f14bpmap->end())
	{
		overlapReads++;
		return &iter->second;
		//	Position1* rPos=nullptr;
		//	rPos->x=((iter->second) & coorXMask)>>coorBits;
		//	rPos->y=(iter->second) & coorMask;
		//	return rPos;
	}
	else if (mismatch > 0)
	{
		Position1 *rv;
		if (mismatch == 1)
		{
#ifdef LOADH5_OPENMP
			if (mOptions->useBf)
			{
				rv = getF14MisOverlapWithBloomFiler(barcodeInt);
			}
			else
			{
				rv = getF14MisOverlap(barcodeInt);
			}
#else
			rv = getF14MisOverlap(barcodeInt);
#endif
		}
		else
		{
			rv = getF14MisOverlap(barcodeInt);
		}
		if (rv != nullptr)
		{
			overlapReadsWithMis++;
			return rv;
		}
		else
		{
			return nullptr;
		}
	}
	return nullptr;
}
#endif
// uint64_t BarcodeProcessor::getF14Position(uint64_t barcodeInt)
//{
//   uint64_t rv;
//   auto iter = f14bpmap->find(barcodeInt);
//   if (iter!=f14bpmap->end()) {
//	overlapReads++;
//	return iter->second;
//   }
//   else if (mismatch > 0) {
////	Position1* rv= nullptr;
//	if(mismatch==1) {
//#ifdef LOADH5_OPENMP
//	  if(mOptions->useBf) {
//		rv=getF14MisOverlapWithBloomFiler(barcodeInt);
//	  }else{
//		rv=getF14MisOverlap(barcodeInt);
//	  }
//	  if(rv!=-1){
//		overlapReadsWithMis++;
//		return rv;
//	  }
//#else
//	  rv = getF14MisOverlap(barcodeInt);
//#endif
//	}else{
//	  rv = getF14MisOverlap(barcodeInt);
//	  return rv;
//	}
//  }
//  return -1;
//}
Position2 *BarcodeProcessor::getPosition2(uint64_t barcodeInt)
{
	unordered_map<uint64_t, Position2>::iterator iter = bpmap2->find(barcodeInt);
	if (iter != bpmap2->end())
	{
		overlapReads++;
		return &iter->second;
	}
	else if (mismatch > 0)
	{
		iter = getMisOverlap2(barcodeInt);
		if (iter != bpmap2->end())
		{
			overlapReadsWithMis++;
			return &iter->second;
		}
		else
		{
			return nullptr;
		}
	}
	return nullptr;
}
Position2 *BarcodeProcessor::getF14Position2(uint64_t barcodeInt)
{
	auto iter = f14bpmap2->find(barcodeInt);
	if (iter != f14bpmap2->end())
	{
		overlapReads++;
		return &iter->second;
	}
	else if (mismatch > 0)
	{
		Position2 *rv = getF14MisOverlap2(barcodeInt);
		if (rv != nullptr)
		{
			overlapReadsWithMis++;
			return rv;
		}
		else
		{
			return nullptr;
		}
	}
	return nullptr;
}
// Position1* BarcodeProcessor::getPosition(uint64_t barcodeInt,string readName)
//{
//   unordered_map<uint64_t, Position1>::iterator iter = bpmap->find(barcodeInt);
//   if(checknames.find(barcodeInt)!=checknames.end()) {
//	vector<string> eles;
//	split(readName, eles, " ");
//	if(checknames2.find(eles[0])!=checknames2.end()) {
//	  cout << readName << "\t" << barcodeInt << "\t";
//	  if(iter!=bpmap->end()) {
//		cout << "perfert found\t" << &iter->second << endl;
//	  } else {
//		iter = getMisOverlap(barcodeInt);
//		if(iter!=bpmap->end()) {
//		  cout << "mismatch found\t" << (&iter->second)->x << "\t" << (&iter->second)->y << "\t";
//		  uint64_t misBarcodeInt;
//		  int misCount = 0;
//		  int misMaskIndex = 0;
//		  unordered_map<uint64_t, Position1>::iterator iter;
//		  unordered_map<uint64_t, Position1>::iterator overlapIter;
//
//		  for (int mis = 0; mis < mismatch; mis++) {
//			misCount = 0;
//			while (misMaskIndex < misMaskLens[mis]) {
//			  misBarcodeInt = barcodeInt^misMask[misMaskIndex];
//			  misMaskIndex++;
//			  iter = bpmap->find(misBarcodeInt);
//			  if(iter!=bpmap->end()) {
//				overlapIter = iter;
//				misCount++;
//				if(misCount > 1) {
//				  break;
//				}
//			  }
//			}
//			if(misCount==1) {
//			  cout << misBarcodeInt << endl;
//			  break;
//			}
//		  }
//		} else {
//		  cout << "mismatch not found" << endl;
//		}
//	  }
//	}
//   }
//   if (iter!=bpmap->end()) {
//	overlapReads++;
//	return &iter->second;
//   }
//   else if (mismatch > 0) {
//	iter = getMisOverlap(barcodeInt);
//	if (iter != bpmap->end()) {
//	  overlapReadsWithMis++;
//	  return &iter->second;
//	}
//	else {
//	  return nullptr;
//	}
//   }
//   return nullptr;
// }

Position1 *BarcodeProcessor::getPosition(string &barcodeString)
{
	int Nindex = getNindex(barcodeString);
	if (Nindex == -1)
	{
		uint64_t barcodeInt = seqEncode2(barcodeString.c_str(), 0, barcodeLen);
		if (barcodeInt == polyTInt)
		{
			return nullptr;
		}
		return getPosition(barcodeInt);
	}
	else if (Nindex == -2)
	{
		return nullptr;
	}
	else if (mismatch > 0)
	{
		return getNOverlap(barcodeString, Nindex);
	}
	return nullptr;
}

void BarcodeProcessor::misMaskGenerate()
{
	misMaskLen = possibleMis(barcodeLen, mismatch);
	misMaskLens = new int[mismatch];
	for (int i = 0; i < mismatch; i++)
	{
		misMaskLens[i] = possibleMis(barcodeLen, i + 1);
	}

	misMask = (uint64_t *)malloc(misMaskLen * sizeof(uint64_t));
	set<uint64_t> misMaskSet;
	int index = 0;
	if (mismatch > 0)
	{
		for (int i = 0; i < barcodeLen; i++)
		{
			for (uint64_t j = 1; j < 4; j++)
			{
				uint64_t misMaskInt = j << i * 2;
				misMask[index] = misMaskInt;
				index++;
			}
		}
		if (mOptions->verbose)
		{
			string msg = "1 mismatch mask barcode number: " + to_string(index);
			//			loginfo(msg);
		}
	}
	if (mismatch == 2)
	{
		misMaskSet.clear();
		for (int i = 0; i < barcodeLen; i++)
		{
			for (uint64_t j = 1; j < 4; j++)
			{
				uint64_t misMaskInt1 = j << i * 2;
				for (int k = 0; k < barcodeLen; k++)
				{
					if (k == i)
					{
						continue;
					}
					for (uint64_t j2 = 1; j2 < 4; j2++)
					{
						uint64_t misMaskInt2 = j2 << k * 2;
						uint64_t misMaskInt = misMaskInt1 | misMaskInt2;
						misMaskSet.insert(misMaskInt);
					}
				}
			}
		}
		for (auto iter = misMaskSet.begin(); iter != misMaskSet.end(); iter++)
		{
			misMask[index] = *iter;
			index++;
		}
		if (mOptions->verbose)
		{
			string msg = "2 mismatch mask barcode number: " + to_string(misMaskSet.size());
			//			loginfo(msg);
		}
	}
	if (mismatch == 3)
	{
		misMaskSet.clear();
		for (int i = 0; i < barcodeLen; i++)
		{
			for (uint64_t j = 1; j < 4; j++)
			{
				uint64_t misMaskInt1 = j << i * 2;
				for (int k = 0; k < barcodeLen; k++)
				{
					if (k == i)
					{
						continue;
					}
					for (uint64_t j2 = 1; j2 < 4; j2++)
					{
						uint64_t misMaskInt2 = j2 << k * 2;
						for (int h = 0; h < barcodeLen; h++)
						{
							if (h == k || h == i)
							{
								continue;
							}
							for (uint64_t j3 = 1; j3 < 4; j3++)
							{
								uint64_t misMaskInt3 = j3 << h * 2;
								uint64_t misMaskInt = misMaskInt1 | misMaskInt2 | misMaskInt3;
								misMaskSet.insert(misMaskInt);
							}
						}
					}
				}
			}
		}
		for (auto iter = misMaskSet.begin(); iter != misMaskSet.end(); iter++)
		{
			misMask[index] = *iter;
			index++;
		}
		if (mOptions->verbose)
		{
			string msg = "3 mismatch mask barcode number: " + to_string(misMaskSet.size());
			//			loginfo(msg);
		}
		misMaskSet.clear();
	}
	if (mOptions->verbose)
	{
		string msg = "total mismatch mask length: " + to_string(misMaskLen);
		//		loginfo(msg);
	}
}

string BarcodeProcessor::positionToString(Position *position)
{
	stringstream positionString;
	positionString << position->x << "_" << position->y;
	return positionString.str();
}

string BarcodeProcessor::positionToString(Position1 *position)
{
	stringstream positionString;
	positionString << position->x << "_" << position->y;
	return positionString.str();
}

string BarcodeProcessor::positionToString(Position2 *position)
{
	stringstream positionString;
	positionString << position->getX() << "_" << position->getY();
	return positionString.str();
}

unordered_map<uint64_t, Position1>::iterator BarcodeProcessor::getMisOverlap(uint64_t barcodeInt)
{
	uint64_t misBarcodeInt;
	int misCount = 0;
	int misMaskIndex = 0;
	unordered_map<uint64_t, Position1>::iterator iter;
	unordered_map<uint64_t, Position1>::iterator overlapIter;

	for (int mis = 0; mis < mismatch; mis++)
	{
		misCount = 0;
		while (misMaskIndex < misMaskLens[mis])
		{
			misBarcodeInt = barcodeInt ^ misMask[misMaskIndex];
			misMaskIndex++;
			iter = bpmap->find(misBarcodeInt);
			if (iter != bpmap->end())
			{
				overlapIter = iter;
				misCount++;
				if (misCount > 1)
				{
					return bpmap->end();
				}
			}
		}
		if (misCount == 1)
		{
			return overlapIter;
		}
	}
	return bpmap->end();
}
#ifdef OPT_MEMORY_PLAN_B
unordered_map<uint64_t, Position1>::iterator BarcodeProcessor::getMisOverlap(uint64_t barcodeInt, uint64_t &foundBarcode)
{
	uint64_t misBarcodeInt;
	int misCount = 0;
	int misMaskIndex = 0;
	unordered_map<uint64_t, Position1>::iterator iter;
	unordered_map<uint64_t, Position1>::iterator overlapIter;

	for (int mis = 0; mis < mismatch; mis++)
	{
		misCount = 0;
		while (misMaskIndex < misMaskLens[mis])
		{
			misBarcodeInt = barcodeInt ^ misMask[misMaskIndex];
			misMaskIndex++;
			iter = bpmap->find(misBarcodeInt);
			if (iter != bpmap->end())
			{
				foundBarcode = misBarcodeInt;
				overlapIter = iter;
				misCount++;
				if (misCount > 1)
				{
					return bpmap->end();
				}
			}
		}
		if (misCount == 1)
		{
			return overlapIter;
		}
	}
	return bpmap->end();
}
#endif
// uint64_t BarcodeProcessor::getF14MisOverlap(uint64_t barcodeInt)
//{
//   uint64_t misBarcodeInt;
//   int misCount = 0;
//   int misMaskIndex = 0;
//   auto iter=f14bpmap->end();
//   auto overlapIter=f14bpmap->end();
//
//   for (int mis = 0; mis < mismatch; mis++){
//	misCount = 0;
//	while (misMaskIndex < misMaskLens[mis]) {
//	  misBarcodeInt = barcodeInt ^ misMask[misMaskIndex];
//	  misMaskIndex++;
//	  iter = f14bpmap->find(misBarcodeInt);
//	  if (iter != f14bpmap->end()) {
//		overlapIter = iter;
//		misCount++;
//		if (misCount > 1) {
//		  return -1;
//		}
//	  }
//	}
//	if (misCount == 1) {
//	  return overlapIter->second;
//	}
//   }
//   return -1;
// }
#ifdef OPT_MEMORY_PLAN_B
Position1 *BarcodeProcessor::getF14MisOverlap(uint64_t barcodeInt, uint64_t &foundBarcode)
{
	uint64_t misBarcodeInt;
	int misCount = 0;
	int misMaskIndex = 0;
	auto iter = f14bpmap->end();
	auto overlapIter = f14bpmap->end();

	for (int mis = 0; mis < mismatch; mis++)
	{
		misCount = 0;
		while (misMaskIndex < misMaskLens[mis])
		{
			misBarcodeInt = barcodeInt ^ misMask[misMaskIndex];
			misMaskIndex++;
			iter = f14bpmap->find(misBarcodeInt);
			if (iter != f14bpmap->end())
			{
				foundBarcode = misBarcodeInt;
				overlapIter = iter;
				misCount++;
				if (misCount > 1)
				{
					return nullptr;
				}
			}
		}
		if (misCount == 1)
		{
			return &overlapIter->second;
		}
	}
	return nullptr;
}
#else
Position1 *BarcodeProcessor::getF14MisOverlap(uint64_t barcodeInt)
{
	uint64_t misBarcodeInt;
	int misCount = 0;
	int misMaskIndex = 0;
	auto iter = f14bpmap->end();
	auto overlapIter = f14bpmap->end();

	for (int mis = 0; mis < mismatch; mis++)
	{
		misCount = 0;
		while (misMaskIndex < misMaskLens[mis])
		{
			misBarcodeInt = barcodeInt ^ misMask[misMaskIndex];
			misMaskIndex++;
			iter = f14bpmap->find(misBarcodeInt);
			if (iter != f14bpmap->end())
			{
				overlapIter = iter;
				misCount++;
				if (misCount > 1)
				{
					return nullptr;
				}
			}
		}
		if (misCount == 1)
		{
			return &overlapIter->second;
		}
	}
	return nullptr;
}
#endif
#ifdef OPT_MEMORY_PLAN_B
Position1 *BarcodeProcessor::getF14MisOverlapWithBloomFiler(uint64_t barcodeInt, uint64_t &foundBarcode)
{
	int misCount = 0;
	/*
	 *  处理mismatch == 1 的情况
	 */
	auto iter = f14bpmap->end();
	auto overlapIter = f14bpmap->end();
	// 处理低32位
	for (int i = 0; i < 16 * 3; i++)
	{
		uint64 misBarcodeInt = barcodeInt ^ misMask[i];
		if (bloomFilter->get_Classification(misBarcodeInt))
		{
			iter = f14bpmap->find(misBarcodeInt);
			if (iter != f14bpmap->end())
			{
				foundBarcode = misBarcodeInt;
				overlapIter = iter;
				misCount++;
				if (misCount > 1)
				{
					return nullptr;
				}
			}
		}
	}
	//  if (misCount == 1) {
	//	return &overlapIter->second;
	//  }
	// 处理高位
	//  iter=f14bpmap->end();
	//  overlapIter=f14bpmap->end();
	//  misCount=0;
	if (bloomFilter->get_Classification(barcodeInt))
	{
		for (int i = 16 * 3; i < misMaskLen; i++)
		{
			uint64 misBarcodeInt = barcodeInt ^ misMask[i];
			iter = f14bpmap->find(misBarcodeInt);
			if (iter != f14bpmap->end())
			{
				foundBarcode = misBarcodeInt;
				overlapIter = iter;
				misCount++;
				if (misCount > 1)
				{
					return nullptr;
				}
			}
		}
	}
	if (misCount == 1)
	{
		return &overlapIter->second;
	}
	return nullptr;
}
#else
Position1 *BarcodeProcessor::getF14MisOverlapWithBloomFiler(uint64_t barcodeInt)
{
	int misCount = 0;
	/*
	 *  处理mismatch == 1 的情况
	 */
	auto iter = f14bpmap->end();
	auto overlapIter = f14bpmap->end();
	// 处理低32位
	for (int i = 0; i < 16 * 3; i++)
	{
		uint64 misBarcodeInt = barcodeInt ^ misMask[i];
		if (bloomFilter->get_Classification(misBarcodeInt))
		{
			iter = f14bpmap->find(misBarcodeInt);
			if (iter != f14bpmap->end())
			{
				overlapIter = iter;
				misCount++;
				if (misCount > 1)
				{
					return nullptr;
				}
			}
		}
	}
	//  if (misCount == 1) {
	//	return &overlapIter->second;
	//  }
	// 处理高位
	//  iter=f14bpmap->end();
	//  overlapIter=f14bpmap->end();
	//  misCount=0;
	if (bloomFilter->get_Classification(barcodeInt))
	{
		for (int i = 16 * 3; i < misMaskLen; i++)
		{
			uint64 misBarcodeInt = barcodeInt ^ misMask[i];
			iter = f14bpmap->find(misBarcodeInt);
			if (iter != f14bpmap->end())
			{
				overlapIter = iter;
				misCount++;
				if (misCount > 1)
				{
					return nullptr;
				}
			}
		}
	}
	if (misCount == 1)
	{
		return &overlapIter->second;
	}
	return nullptr;
}
#endif
// uint64_t BarcodeProcessor::getF14MisOverlapWithBloomFiler(uint64_t barcodeInt){
//   int misCount = 0;
///*
// *  处理mismatch == 1 的情况
// */
//  auto iter=f14bpmap->end();
//  auto overlapIter=f14bpmap->end();
//  // 处理低32位
//  for (int i = 0; i < 16 * 3; i++) {
//	uint64 misBarcodeInt = barcodeInt ^ misMask[i];
//	if (bloomFilter->get_Classification(misBarcodeInt)){
//	  iter = f14bpmap->find(misBarcodeInt);
//	  if (iter != f14bpmap->end()) {
//		overlapIter = iter;
//		misCount++;
//		if (misCount > 1) {
//		  return -1;
//		}
//	  }
//	}
//  }
//  if (bloomFilter->get_Classification(barcodeInt)) {
//	for (int i = 16 * 3; i < misMaskLen; i++) {
//	  uint64 misBarcodeInt = barcodeInt ^ misMask[i];
//	  iter = f14bpmap->find(misBarcodeInt);
//	  if (iter != f14bpmap->end()) {
//		overlapIter = iter;
//		misCount++;
//		if (misCount > 1) {
//		  return -1;
//		}
//	  }
//	}
//  }
//  if (misCount == 1) {
//	return overlapIter->second;
//  }
//  return -1;
//}
unordered_map<uint64_t, Position2>::iterator BarcodeProcessor::getMisOverlap2(uint64_t barcodeInt)
{
	uint64_t misBarcodeInt;
	int misCount = 0;
	int misMaskIndex = 0;
	unordered_map<uint64_t, Position2>::iterator iter;
	unordered_map<uint64_t, Position2>::iterator overlapIter;

	for (int mis = 0; mis < mismatch; mis++)
	{
		misCount = 0;
		while (misMaskIndex < misMaskLens[mis])
		{
			misBarcodeInt = barcodeInt ^ misMask[misMaskIndex];
			misMaskIndex++;
			iter = bpmap2->find(misBarcodeInt);
			if (iter != bpmap2->end())
			{
				overlapIter = iter;
				misCount++;
				if (misCount > 1)
				{
					return bpmap2->end();
				}
			}
		}
		if (misCount == 1)
		{
			return overlapIter;
		}
	}
	return bpmap2->end();
}
Position2 *BarcodeProcessor::getF14MisOverlap2(uint64_t barcodeInt)
{
	uint64_t misBarcodeInt;
	int misCount = 0;
	int misMaskIndex = 0;
	auto iter = f14bpmap2->end();
	auto overlapIter = f14bpmap2->end();

	for (int mis = 0; mis < mismatch; mis++)
	{
		misCount = 0;
		while (misMaskIndex < misMaskLens[mis])
		{
			misBarcodeInt = barcodeInt ^ misMask[misMaskIndex];
			misMaskIndex++;
			iter = f14bpmap2->find(misBarcodeInt);
			if (iter != f14bpmap2->end())
			{
				overlapIter = iter;
				misCount++;
				if (misCount > 1)
				{
					return nullptr;
				}
			}
		}
		if (misCount == 1)
		{
			return &overlapIter->second;
		}
	}
	return nullptr;
}

Position1 *BarcodeProcessor::getNOverlap(string &barcodeString, uint8 Nindex)
{
	// N has the same encode (11) with G
	int misCount = 0;
	uint64_t barcodeInt = seqEncode(barcodeString.c_str(), 0, barcodeString.length());
	unordered_map<uint64_t, Position1>::iterator iter;
	unordered_map<uint64_t, Position1>::iterator overlapIter;
	iter = bpmap->find(barcodeInt);
	if (iter != bpmap->end())
	{
		misCount++;
		overlapIter = iter;
	}
	for (uint64_t j = 1; j < 4; j++)
	{
		uint64_t misBarcodeInt = barcodeInt ^ (j << Nindex * 2);
		iter = bpmap->find(misBarcodeInt);
		if (iter != bpmap->end())
		{
			misCount++;
			if (misCount > 1)
			{
				return nullptr;
			}
			overlapIter = iter;
		}
	}
	if (misCount == 1)
	{
		overlapReadsWithN++;
		return &overlapIter->second;
	}
	return nullptr;
}

int BarcodeProcessor::getNindex(string &barcodeString)
{
	string::size_type Nindex = barcodeString.find("N");
	if (Nindex == barcodeString.npos)
	{
		return -1;
	}
	else if (Nindex != barcodeString.rfind("N"))
	{
		return -2;
	}
	return Nindex;
}

void BarcodeProcessor::addDNB(uint64_t barcodeInt)
{
	if (mDNB.count(barcodeInt) > 0)
	{
		mDNB[barcodeInt]++;
	}
	else
	{
		mDNB[barcodeInt] = 1;
	}
}

bool BarcodeProcessor::barcodeStatAndFilter(pair<string, string> &barcode)
{
	for (int i = 0; i < barcodeLen; i++)
	{
		if (barcode.second[i] >= q30)
		{
			barcodeQ30++;
			barcodeQ20++;
			barcodeQ10++;
		}
		else if (barcode.second[i] >= q20)
		{
			barcodeQ20++;
			barcodeQ10++;
		}
		else if (barcode.second[i] >= q10)
		{
			barcodeQ10++;
		}
	}
	return true;
}

bool BarcodeProcessor::barcodeStatAndFilter(string &barcodeQ)
{
	for (int i = 0; i < barcodeLen; i++)
	{
		if (barcodeQ[i] >= q30)
		{
			barcodeQ30++;
			barcodeQ20++;
			barcodeQ10++;
		}
		else if (barcodeQ[i] >= q20)
		{
			barcodeQ20++;
			barcodeQ10++;
		}
		else if (barcodeQ[i] >= q10)
		{
			barcodeQ10++;
		}
	}
	return true;
}

bool BarcodeProcessor::umiStatAndFilter(pair<string, string> &umi)
{
	int q10BaseCount = 0;
	for (int i = 0; i < mOptions->transBarcodeToPos.umiLen; i++)
	{
		if (umi.second[i] >= q30)
		{
			umiQ30++;
			umiQ20++;
			umiQ10++;
		}
		else if (umi.second[i] >= q20)
		{
			umiQ20++;
			umiQ10++;
		}
		else if (umi.second[i] >= q10)
		{
			umiQ10++;
		}
		else
		{
			q10BaseCount++;
		}
	}
	if (umi.first.find("N") != string::npos)
	{
		umiNFilterReads++;
		return false;
	}
	else if (seqEncode(umi.first.c_str(), 0, mOptions->transBarcodeToPos.umiLen) == 0)
	{
		umiPloyAFilterReads++;
		return false;
	}
	else if (q10BaseCount > 1)
	{
		umiQ10FilterReads++;
		return false;
	}
	else
	{
		return true;
	}
}

void BarcodeProcessor::dumpDNBmap(string &dnbMapFile)
{
	ofstream writer;
	//	if (ends_with(dnbMapFile, ".bin")) {
	//		mDNB.reserve(mDNB.size());
	//		writer.open(dnbMapFile, ios::out | ios::binary);
	////		boost::archive::binary_oarchive oa(writer);
	////		oa << mDNB;
	//	}
	//	else {
	writer.open(dnbMapFile);
	unordered_map<uint64_t, int>::iterator mapIter = mDNB.begin();
	while (mapIter != mDNB.end())
	{
		uint32_t x = mapIter->first >> 32;
		uint32_t y = mapIter->first & 0x00000000FFFFFFFF;
		writer << x << "\t" << y << "\t" << mapIter->second << "\n";
		mapIter++;
	}
	//	}
	writer.close();
}
