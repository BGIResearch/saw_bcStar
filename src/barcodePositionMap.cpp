/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#include "barcodePositionMap.h"
#include <chrono>
#include <unistd.h>
//#include "/Users/berry/Documents/software/hps-1.0.0/src/hps.h"
//#include <folly/container/F14Map.h>
//#include <folly/container/detail/F14Table.h>
// using namespace::folly;
BarcodePositionMap::BarcodePositionMap(Options *opt)
{
	mOptions = opt;
	split(opt->in, inFile, ",");
	maskFile = opt->maskFile;
	barcodeStart = opt->barcodeStart;
	barcodeLen = opt->barcodeLen;
	segment = opt->barcodeSegment;
	turnFovDegree = opt->turnFovDegree;
	isSeq500 = opt->isSeq500;
	rc = opt->rc;
#ifndef OPT_MEMORY_IN_H5
	bloomFilter = new BloomFilter();
#endif
	//#ifdef D_PRINT
	//  	cout<<"before construct...\t"<<bloomFilter->hashtableClassification[0]<<"\t"<<bloomFilter->hashtableClassification[(1ll<<26)-1]<<endl;
	//#endif
	// polyTint = getPolyTint(barcodeLen);
	readidSep = opt->barcodeStat.readidSep;
	inFastqNumber = inFile.size();
	firstFastq = *inFile.begin();
	initiate();
	splitNum = 64;
	int splitUseBaseNum = 8; // defined in writeFq
	if (mOptions->splitBarcode)
	{
		splitBpmap = new unordered_map<uint32_t, Position1>[(1 << splitUseBaseNum) - 1];
	}
	//    if(opt->mapSize!=100000000){
	//        bpmap=unordered_map<uint64_t,Position1>(opt->mapSize);
	//    }
	if (ends_with(firstFastq, "fq") || ends_with(firstFastq, "fastq") || ends_with(firstFastq, "fq.gz") || ends_with(firstFastq, "fastq.gz"))
	{
		// if (inFastqNumber == 1){
		//	loadMapFromFastq();
		// }else {
		//	loadMapFromMultiFastq();
		// }
	}
	else
	{
		for (vector<string>::iterator ix = inFile.begin(); ix != inFile.end(); ix++)
		{
			bool firstSeqInputAsList = false;
			ifstream tmpIn(*ix);
			if (!tmpIn)
			{
				string errMsg = "barcodePositionMapFile does not exists: " + inFile[0];
				sawErrCode(err_fileOpen_failed, errMsg);
				cerr << "barcodePositionMapFile does not exists: " << inFile[0] << endl;
				exit(1);
			}
			string tmpLine;
			getline(tmpIn, tmpLine);
			if (!tmpLine.empty())
			{
				ifstream tmpIn2(tmpLine);
				if (tmpIn2)
				{
					firstSeqInputAsList = true;
				}
				tmpIn2.close();
			}
			tmpIn.close();
			if (firstSeqInputAsList)
			{
				ifstream inList(*ix);
				string filePath;
				while (getline(inList, filePath))
				{
					loadbpmap(filePath);
				}
				inList.close();
			}
			else
			{
				loadbpmap(*ix);
			}
		}
		//		loadbpmap();
	}
	// debug loadbpmap
	//	exit(1);
	//#ifdef D_PRINT
	//  cout<<"after construct...\t"<<bloomFilter->hashtableClassification[0]<<"\t"<<bloomFilter->hashtableClassification[(1ll<<26)-1]<<endl;
	//#endif
}

BarcodePositionMap::~BarcodePositionMap()
{
	bpmap.clear();
	unordered_map<uint64_t, Position1>().swap(bpmap);
	dupBarcode.clear();
	set<uint64_t>().swap(dupBarcode);
	inFile.clear();
	vector<string>().swap(inFile);
	delete[] totalReads;
	delete[] readsWithN;
	delete[] dupReads;
	delete[] ESTdupReads;
	for (int i = 0; i < inFastqNumber; i++)
	{
		delete[] polyReads[i];
	}
	delete[] polyReads;
	delete[] readsQ10;
	delete[] readsQ20;
	delete[] readsQ30;
	delete[] totalBase;
	delete[] totalBarcodes;
}

void BarcodePositionMap::initiate()
{
	totalReads = new long[inFastqNumber]();
	readsWithN = new long[inFastqNumber]();
	readsWithoutPos = new long[inFastqNumber]();
	polyReads = new long *[inFastqNumber]();
	for (int i = 0; i < inFastqNumber; i++)
	{
		polyReads[i] = new long[4]();
	}
	dupReads = new long[inFastqNumber]();
	ESTdupReads = new long[inFastqNumber]();
	readsQ10 = new long[inFastqNumber]();
	readsQ20 = new long[inFastqNumber]();
	readsQ30 = new long[inFastqNumber]();
	totalBase = new long[inFastqNumber]();
	totalBarcodes = new long[inFastqNumber]();
}


long BarcodePositionMap::getBarcodeTypes()
{
	return bpmap.size();
}

void BarcodePositionMap::dumpbpmap(string &mapOutFile)
{
	time_t start = time(NULL);
	cout << "##########dump barcodeToPosition map begin..." << endl;
	if (ends_with(mapOutFile, ".bin"))
	{
		bpmap.reserve(bpmap.size());
		ofstream writer(mapOutFile, ios::out | ios::binary);
		//		boost::archive::binary_oarchive oa(writer);
		//		oa << bpmap;
		// while (mapIter != bpmap.end()) {
		//	writer.write((char*)&mapIter->first, sizeof(uint64_t));
		//	writer.write((char*)&mapIter->second, sizeof(Position));
		//	mapIter++;
		//}
		writer.close();
	}
	else if (ends_with(mapOutFile, "h5") || ends_with(mapOutFile, "hdf5"))
	{
		ChipMaskHDF5 chipMaskH5(mapOutFile);
		chipMaskH5.creatFile();
		int slidePitch = sm->getSlidePitch(mOptions->platform);
		chipMaskH5.writeDataSet(mOptions->chipID, sliderange, bpmap, barcodeLen, mOptions->barcodeSegment, slidePitch, mOptions->compression);
	}
	else if (ends_with(mapOutFile, ".hps"))
	{
	}
	else
	{
		ofstream writer(mapOutFile);
		unordered_map<uint64_t, Position1>::iterator mapIter = bpmap.begin();
		while (mapIter != bpmap.end())
		{
			writer << seqDecode(mapIter->first, barcodeLen) << "\t" << mapIter->second.x << "\t" << mapIter->second.y << endl;
			mapIter++;
		}
		writer.close();
	}
	cout << "##########dump barcodeToPosition map finished, time used: " << time(NULL) - start << " seconds" << endl;
}

void BarcodePositionMap::loadbpmap(string barcodePositionMapFile)
{
	time_t start = time(NULL);
	//    string barcodePositionMapFile = inFile.at(0);
	if (!file_exists(barcodePositionMapFile))
	{
		string errMsg = "barcodePositionMapFile does not exists: " + barcodePositionMapFile;
		sawErrCode(err_fileOpen_failed, errMsg);
		cerr << "barcodePositionMapFile does not exists: " << barcodePositionMapFile << endl;
		exit(1);
	}

	cout << "###############load barcodeToPosition map begin..." << endl;
	// cout << "###############barcode map file: " << barcodePositionMapFile << endl;
#ifdef D_PRINT
	pid_t pd = getpid();
	if (barcodePositionMapFile.size() > 1)
	{
		cout << "current memory use " << (float)get_proc_VmRSS(pd) / (1024 * 1024) << "\t" << (float)get_proc_virtualMem(pd) / (1024 * 1024) << "\t" << (float)get_proc_virtualMemPeak(pd) / (1024 * 1024) << endl;
	}
	//		exit(1);
#endif
	if (ends_with(barcodePositionMapFile, ".bin"))
	{
		cout << "bin file:\t" << barcodePositionMapFile << endl;
		ifstream mapReader(barcodePositionMapFile, ios::in | ios::binary);
		if (!mapReader.is_open())
		{
			throw invalid_argument("Could not open the file: " + barcodePositionMapFile);
		}
		uint64_t barcodeInt = 0;
		uint32_t posX, posY;
		Position1 v;
		if (mOptions->useF14)
		{
			if (mOptions->bcNum > 0)
			{
				f14bpmap.reserve(mOptions->bcNum);
			}
		}
		while (!mapReader.eof())
		{
			mapReader.read((char *)&barcodeInt, sizeof(barcodeInt));
			mapReader.read((char *)&posX, sizeof(posX));
			mapReader.read((char *)&posY, sizeof(posY));
			v.x = posX;
			v.y = posY;
#ifdef OPT_MEMORY_PLAN_B
			v.count = 0;
			v.fullTag = 0;
#endif
			if (mOptions->useF14)
			{
				f14bpmap[barcodeInt] = v;
			}
			else
			{
				bpmap[barcodeInt] = v;
			}
		}
		mapReader.close();
#ifdef D_PRINT
		cout << "current memory use " << (float)get_proc_VmRSS(pd) / (1024 * 1024) << "\t" << (float)get_proc_virtualMem(pd) / (1024 * 1024) << "\t" << (float)get_proc_virtualMemPeak(pd) / (1024 * 1024) << endl;
#endif
	}
	else if (ends_with(barcodePositionMapFile, "h5") || ends_with(barcodePositionMapFile, "hdf5"))
	{
		ChipMaskHDF5 chipMaskH5(barcodePositionMapFile);
		chipMaskH5.openFile();
#ifdef D_PRINT
		pid_t pd = getpid();

		cout << "current memory use " << (float)get_proc_VmRSS(pd) / (1024 * 1024) << "\t" << (float)get_proc_virtualMem(pd) / (1024 * 1024) << "\t" << (float)get_proc_virtualMemPeak(pd) / (1024 * 1024) << endl;
#endif
#ifdef LOADH5_OPENMP
		// use openmp to parallelize
		mOptions->useBf = true;
		if (mOptions->bcNum > 0)
		{
			f14bpmap.reserve(mOptions->bcNum);
		}
		//	    cout<<"load h5 and use f14"<<endl;
		if (chipMaskH5.readDataSetParallelize(f14bpmap, bloomFilter) < 0)
		{
			cerr << "Warning, h5 file is not a chunked dataset, it will be read by another method" << endl;
			mOptions->useBf = false;
			chipMaskH5.readDataSet(f14bpmap);
		}
#else
		if (mOptions->useF14)
		{
#ifdef D_PRINT
			cout << "using F14..." << endl;
#endif
			chipMaskH5.readDataSet(f14bpmap);
		}
		else
		{
#ifdef D_PRINT
			cout << "using unordered_map..." << endl;
#endif
			chipMaskH5.readDataSet(bpmap);
		}
#endif
//	  exit(1);
#ifdef D_PRINT
		cout << "current memory use " << (float)get_proc_VmRSS(pd) / (1024 * 1024) << "\t" << (float)get_proc_virtualMem(pd) / (1024 * 1024) << "\t" << (float)get_proc_virtualMemPeak(pd) / (1024 * 1024) << endl;
#endif
		//            exit(1);
	}
	else if (ends_with(barcodePositionMapFile, ".hps"))
	{
		// added by gc
		// parse .hps file

		//        ifstream in_file(barcodePositionMapFile,ifstream::binary);
		//        in_file.seekg(0,in_file.end);
		//        size_t serialized_size = in_file.tellg();
		//        bpmap=hps::from_stream<unordered_map<uint64_t,  Position1>>(in_file);
		//        cout<<"size(B): "<<serialized_size<<endl;
		//        in_file.close();
	}
	else if (barcodePositionMapFile.find("prefix") == 0 || barcodePositionMapFile.find("/prefix") != string::npos)
	{
		ifstream mapReader(barcodePositionMapFile);
		string line;
		uint64_t barcodeInt;
		Position1 position;
		//        double nonMap=0;
		//        double maptime=0;
		//        F14FastMap<uint64_t,  Position1> f14map;
		//                F14BasicMap<uint64_t>* f14map=new F14BasicMap<uint64_t>(10,3);
		while (getline(mapReader, line))
		{
			if (line.empty())
			{
				cerr << "barcodePositionMap file read finished." << endl;
				break;
			}
			vector<string> splitLine;
			// std::getline(mapReader, line);
			split(line, splitLine, "\t");
			if (splitLine.size() != 3)
			{
				//			    cerr<<line<<endl;
				cerr << "Error,parse barcode input error" << endl;
				//			    exit(1);
			}
			barcodeInt = seqEncode(splitLine[0].c_str(), barcodeStart, barcodeLen);
			//            auto& position=f14map[barcodeInt];
			position.x = atoi(splitLine[1].c_str());
			position.y = atoi(splitLine[2].c_str());
			//            f14map[barcodeInt]=position;
			//            f14map[barcodeInt]=touint64_t(position);
			//            f14map.insert(barcodeInt,position);
			bpmap[barcodeInt] = position;
		}
		cout << "bpmap load suceessfully. Map size:" << bpmap.size() << endl;
		//        cout<<"nonMap time:"<<nonMap/1000000000<<" s"<<endl;
		//        cout<<"Map time:"<<maptime/1000000000<<" s"<<endl;
		mapReader.close();
	}
	else
	{
		// added by gc
		// load ldxidb
		sm = new SlideMask(maskFile, isSeq500, turnFovDegree);

		uint64_t barcodeInt;
		Position1 position;
		string line;
		ifstream mapReader(barcodePositionMapFile);
		uint64_t dupBarcodeNum = 0;
		unordered_map<uint64_t, set<uint64_t>> potentialDup;
		while (std::getline(mapReader, line))
		{
			if (line.empty())
			{
				cerr << "barcodePositionMap file read finished." << endl;
				break;
			}
			vector<string> splitLine;
			// std::getline(mapReader, line);
			split(line, splitLine, "\t");
			if (splitLine.size() != 4)
			{
				//			    cerr<<line<<endl;
				cerr << "Error,parse barcode input error" << endl;
				//			    exit(1);
			}
			// todo row/10 for debug
			//			int lastNum=splitLine[3][splitLine[3].size()-1]-'0';
			//			pair<int32_t,int32_t> pos=sm->getIdx(atoi(splitLine[1].c_str())+lastNum*10000000,atoi(splitLine[2].c_str()),atoi(splitLine[3].c_str())/10);
			pair<int32_t, int32_t> pos = sm->getIdx(atoi(splitLine[1].c_str()), atoi(splitLine[2].c_str()), atoi(splitLine[3].c_str()));
			position.x = pos.first;
			position.y = pos.second;
			//			if (splitLine.size() < 3){
			//				break;
			//			}
			//			else if (splitLine.size() == 3){
			//				position.x = std::stoi(splitLine[1]);
			//				position.y = std::stoi(splitLine[2]);
			//			}else {
			//				position.x = std::stoi(splitLine[3]);
			//				position.y = std::stoi(splitLine[4]);
			//			}
			string realBarcode = splitLine[0].substr(barcodeStart, barcodeLen);
			barcodeInt = seqEncode2(realBarcode.c_str(), barcodeStart, barcodeLen);
			if (bpmap.find(barcodeInt) != bpmap.end())
			{
				potentialDup[barcodeInt].insert(touint64_t(position));
				potentialDup[barcodeInt].insert(touint64_t(bpmap[barcodeInt]));
			}
			else
			{
				bpmap[barcodeInt] = position;
			}
			//			cout << "barcode: " << barcodeInt << " position: " << position.x << " " << position.y <<endl;
		}
		delete sm;
		// check which barcodes are real duplicates
		for (unordered_map<uint64_t, set<uint64_t>>::iterator ix = potentialDup.begin(); ix != potentialDup.end(); ix++)
		{
			if (isNeighbor(ix->second))
			{
				bpmap[ix->first] = toPosition1(*((ix->second).begin()));
			}
			else
			{
				dupBarcodeNum++;
			}
		}
		potentialDup.clear();
		unordered_map<uint64_t, set<uint64_t>>().swap(potentialDup);
		cout << "bpmap load suceessfully. Map size:" << bpmap.size() << endl;
		mapReader.close();
	}
	cout << "###############load barcodeToPosition map finished, time used: " << time(NULL) - start << " seconds" << endl;
	cout << resetiosflags(ios::fixed) << setprecision(2);
	cout << "getBarcodePositionMap_uniqBarcodeTypes: ";
	if (mOptions->useF14)
	{
		cout << f14bpmap.size();
		if (f14bpmap.size() == 0)
		{
			string errMsg = "no barcode found";
			sawErrCode(err_fileParse_failed, errMsg);
			cerr << "no barcode found" << endl;
			exit(1);
		}
	}
	else
	{
		cout << bpmap.size();
		if (bpmap.size() == 0)
		{
			string errMsg = "no barcode found";
			sawErrCode(err_fileParse_failed, errMsg);
			cerr << "no barcode found" << endl;
			exit(1);
		}
	}
	cout << endl;
}
void BarcodePositionMap::getSuffixLen()
{
	FastqReader reader(firstFastq);
	Read *read = reader.read();
	string readName = read->mName;
	string::size_type sepPos = readName.find(readidSep);
	if (sepPos != readName.npos)
	{
		suffixLen = readName.size() - sepPos;
	}
	else
	{
		suffixLen = 0;
	}
}

int BarcodePositionMap::readQualityStat(string &readQ, int index)
{
	int lowQ10 = 0;
	for (uint32_t i = 0; i < readQ.size(); i++)
	{
		totalBase[index]++;
		if (readQ[i] >= 30 + 33)
		{
			readsQ30[index]++;
			readsQ20[index]++;
			readsQ10[index]++;
		}
		else if (readQ[i] >= 20 + 33)
		{
			readsQ20[index]++;
			readsQ10[index]++;
		}
		else if (readQ[i] >= 10 + 33)
		{
			readsQ10[index]++;
		}
		else
		{
			lowQ10++;
		}
	}
	return lowQ10;
}

/**
 * baseCount: ['A', 'C', 'T', 'G']
 */
bool BarcodePositionMap::barcodeFilter(string &readSeq, int index)
{
	int baseCount[5] = {0, 0, 0, 0, 0};
	int readLen = readSeq.size();
	for (int i = 0; i < readLen; i++)
	{
		switch (readSeq[i])
		{
		case 'A':
			baseCount[0]++;
			break;
		case 'C':
			baseCount[1]++;
			break;
		case 'T':
			baseCount[2]++;
			break;
		case 'G':
			baseCount[3]++;
			break;
		default:
			baseCount[4]++;
		}
	}
	if (baseCount[4] > 0)
	{
		readsWithN[index]++;
		return true;
	}
	for (int i = 0; i < 4; i++)
	{
		float polyRate = baseCount[i] / readLen;
		if (polyRate >= 0.8)
		{
			polyReads[index][i]++;
			return true;
		}
	}
	return false;
}
bool BarcodePositionMap::isEst(Position1 &pos1, Position1 &pos2)
{
	if (fabs(pos1.x - pos2.x) > EST_DNB_DISTANCE || fabs(pos1.y - pos2.y) > EST_DNB_DISTANCE)
	{
		return false;
	}
	else
	{
		return true;
	}
}
uint64_t BarcodePositionMap::touint64_t(Position1 pos)
{
	uint64_t v = uint64_t(pos.x) << 32;
	v += pos.y;
	return v;
}

bool BarcodePositionMap::isNeighbor(uint64_t a, uint64_t b)
{
	Position1 pos1 = toPosition1(a);
	Position1 pos2 = toPosition1(b);
	if (fabs(pos1.x - pos2.x) <= 1 && fabs(pos1.y - pos2.y) <= 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool BarcodePositionMap::isNeighbor(set<uint64_t> manyPos)
{
	if (manyPos.size() > 4)
	{
		return false;
	}
	else
	{
		for (set<uint64_t>::iterator ix = manyPos.begin(); ix != manyPos.end(); ix++)
		{
			set<uint64_t>::iterator ix2 = ix;
			ix2++;
			for (; ix2 != manyPos.end(); ix2++)
			{
				if (!isNeighbor(*ix, *ix2))
				{
					return false;
				}
			}
		}
		return true;
	}
}

Position1 BarcodePositionMap::toPosition1(uint64_t a)
{
	Position1 v;
	v.y = a & 0xFFFFFFFF;
	v.x = a >> 32 & 0xFFFFFFFF;
	return v;
}
