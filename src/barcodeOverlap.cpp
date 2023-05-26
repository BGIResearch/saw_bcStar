/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#include "barcodeOverlap.h"

BarcodeOverlap::BarcodeOverlap(Options *opt)
{
	mOptions = opt;
	map1File = mOptions->in;
	map2File = mOptions->barcodeOverlap.in2;
	overlapMapOut = mOptions->out;
	logOut = mOptions->out + ".stat";
	mismatch = mOptions->barcodeOverlap.mismatch;
	baseMisCount = new uint64_t[mOptions->barcodeLen]{0};
	if (mismatch > 0)
	{
		misMaskGenerate(mismatch);
	}
	barcodeMap = new BarcodeMap(mOptions->barcodeStart, mOptions->barcodeLen, mOptions->mapSize);
	if (ends_with(map1File, ".bin"))
	{
		barcodeMap->loadBarcodeMap(map1File);
	}
	else
	{
		barcodeMap->loadBarcodeMapFromTxt(map1File);
	}
}

BarcodeOverlap::~BarcodeOverlap()
{
	delete barcodeMap;
	overlapMap.clear();
}

void BarcodeOverlap::overlapStat()
{
	uint64_t barcodeInt;
	//	uint64_t misBarcodeInt;
	string barcodeString;
	uint16 bcount;
	uint16 b1count = 0;
	ifstream map2Reader;
	// ofstream mapOutWriter(overlapMapOut);
	if (ends_with(map2File, ".bin"))
	{
		map2Reader.open(map2File, ios::binary);
	}
	else
	{
		map2Reader.open(map2File);
	}
	while (true)
	{
		if (ends_with(map2File, ".bin"))
		{
			map2Reader.read((char *)&barcodeInt, sizeof(uint64_t));
			map2Reader.read((char *)&bcount, sizeof(uint16));
		}
		else
		{
			map2Reader >> barcodeString >> bcount;
			if (barcodeString.find("N") != barcodeString.npos)
				continue;
			barcodeInt = seqEncode(barcodeString.c_str(), mOptions->barcodeStart, mOptions->barcodeLen);
		}
		if (map2Reader.eof())
		{
			break;
		}
		totalReads += bcount;
		if (barcodeMap->bmap.count(barcodeInt) > 0)
		{
			b1count = barcodeMap->bmap[barcodeInt];
			overlapReads += bcount;
			overlapBarcode++;
			if (overlapMap.count(barcodeInt) > 0)
			{
				overlapMap[barcodeInt] += bcount;
			}
			else
			{
				overlapMap[barcodeInt] = bcount;
				// barcodeString = seqDecode(barcodeInt, myOptions->barcodeLen);
				// mapOutWriter << barcodeString << "\t" << bcount << endl;
			}
			if (b1count == 1)
			{
				overlapBarcodeUniq++;
				overlapReadsUniq += bcount;
			}
		}
		else if (mismatch > 0)
		{
			uint64_t misBarcodeInt = getMisOverlap(barcodeInt);
			if (misBarcodeInt == MAX_BARCODE)
			{
				continue;
			}
			misOverlapReads += bcount;
			misOverlapBarcode++;
			misOverlapMap[barcodeInt] = bcount;
			if (overlapMap.count(misBarcodeInt) > 0)
			{
				overlapMap[misBarcodeInt] += bcount;
			}
			else
			{
				overlapMap[misBarcodeInt] = bcount;
			}
			if (b1count == 1)
			{
				misOverlapBarcodeUniq++;
				misOverlapReadsUniq += bcount;
			}
		}
	}
	map2Reader.close();
	// mapOutWriter.close();
	totalBarcode = overlapMap.size();
	barcodeMapDump(&overlapMap, overlapMapOut);
	if (mismatch > 0)
	{
		string misOverlapMapOut = overlapMapOut + ".mis";
		barcodeMapDump(&misOverlapMap, misOverlapMapOut);
	}
	logWrite();
}

void BarcodeOverlap::barcodeMapDump(unordered_map<uint64_t, uint16> *barcodeMap, string outputFile)
{
	cout << "write map to txtfile begin..." << endl;
	time_t start = time(NULL);
	ofstream mapOutWriter(outputFile);
	unordered_map<uint64_t, uint16>::iterator mapIter;
	mapIter = barcodeMap->begin();
	while (mapIter != barcodeMap->end())
	{
		string barcodeString = seqDecode(mapIter->first, mOptions->barcodeLen);
		mapOutWriter << barcodeString << "\t" << mapIter->second << endl;
		mapIter++;
	}
	mapOutWriter.close();
	cout << "barcode dump time: " << (time(NULL) - start) << "Sec" << endl;
}

void BarcodeOverlap::logWrite()
{
	ofstream logWriter(logOut);
	logWriter << setiosflags(ios::fixed) << setprecision(2);
	logWriter << map2File << endl;
	logWriter << "RawReads:\t" << totalReads << endl;
	logWriter << "TotalOverlapBarcode:\t" << totalBarcode << endl;
	logWriter << "OverlapBarcode:\t" << overlapBarcode << endl;
	logWriter << "OverlapReads:\t" << overlapReads << endl;
	logWriter << "OverlapBarcodesUniq:\t" << overlapBarcodeUniq << endl;
	logWriter << "OverlapReadsUniq:\t" << overlapReadsUniq << endl;
	logWriter << "OverlapReadsPercent:\t" << double(overlapReads) / double(totalReads) * 100 << "%" << endl;
	logWriter << "OverlapReadsUniqPercent:\t" << double(overlapReadsUniq) / double(totalReads) * 100 << "%" << endl;
	if (mismatch > 0)
	{
		logWriter << "##################Overlap with  mismatch########################" << endl;
		logWriter << "OverlapReadsWithMis:\t" << misOverlapReads << endl;
		logWriter << "OverlapReadsUniqWithMis:\t" << misOverlapReadsUniq << endl;
		logWriter << "OverlapBarcodeWithMis:\t" << misOverlapBarcode << endl;
		logWriter << "OverlapBarcodeUniqWithMis:\t" << misOverlapBarcodeUniq << endl;
		logWriter << "OverlapReadsWithMisPercent:\t" << double(misOverlapReads) / double(totalReads) * 100 << "%" << endl;
		logWriter << "OverlpaReadsUniqWithMisPercent:\t" << double(misOverlapReadsUniq) / double(totalReads) * 100 << "%" << endl;
		logWriter << "DupOverlapBarcodeWithMis:\t" << dupMisOverlapBarcode << endl;
		for (int i = 0; i < mOptions->barcodeLen; i++)
		{
			logWriter << baseMisCount[i] << "\t";
		}
		logWriter << endl;
		logWriter << "base transfer distribution:" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				logWriter << baseTransCount[i][j] << "\t";
			}
			logWriter << endl;
		}
	}
	logWriter.close();
}

void BarcodeOverlap::setToArray(set<uint64_t> &mset, uint64_t *marray)
{
	marray = new uint64_t[mset.size()];
	int i = 0;
	for (auto iter = mset.begin(); iter != mset.end(); iter++)
	{
		marray[i] = *iter;
		i++;
	}
}

void BarcodeOverlap::misMaskGenerate(int mismatch)
{
	int PossibleMisMaskLen = possibleMis(mOptions->barcodeLen, mismatch);
	if (PossibleMisMaskLen == 0)
	{
		string errMsg = "allowed mismatch is 1, 2, 3, please check the mismak options you give.";
		sawErrCode(err_param_invalid_bcmap, errMsg);
		error_exit("allowed mismatch is 1, 2, 3, please check the mismak options you give.");
	}
	misMask = (uint64_t *)malloc(PossibleMisMaskLen * sizeof(uint64_t));
	set<uint64_t> misMaskSet;
	int index = 0;
	if (mismatch > 0)
	{
		for (int i = 0; i < mOptions->barcodeLen; i++)
		{
			for (uint64_t j = 1; j < 4; j++)
			{
				uint64_t misMaskInt = j << i * 2;
				misMask[index] = misMaskInt;
				index++;
			}
		}
		cout << "1 mismatch mask barcode number: " << index << endl;
	}
	if (mismatch >= 2)
	{
		for (int i = 0; i < mOptions->barcodeLen; i++)
		{
			for (uint64_t j = 1; j < 4; j++)
			{
				uint64_t misMaskInt1 = j << i * 2;
				for (int k = 0; k < mOptions->barcodeLen; k++)
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
		cout << "2 mismatch mask barcode number: " << misMaskSet.size() << endl;
		misMaskSet.clear();
	}
	if (mismatch == 3)
	{
		for (int i = 0; i < mOptions->barcodeLen; i++)
		{
			for (uint64_t j = 1; j < 4; j++)
			{
				uint64_t misMaskInt1 = j << i * 2;
				for (int k = 0; k < mOptions->barcodeLen; k++)
				{
					if (k == i)
					{
						continue;
					}
					for (uint64_t j2 = 1; j2 < 4; j2++)
					{
						uint64_t misMaskInt2 = j2 << k * 2;
						for (int h = 0; h < mOptions->barcodeLen; h++)
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
		cout << "3 mismatch mask barcode number: " << misMaskSet.size() << endl;
		misMaskSet.clear();
	}
	misMaskLen = index;
	cout << "total mismatch mask length: " << index << endl;
}

uint64_t BarcodeOverlap::getMisOverlap(uint64_t barcodeInt)
{
	uint64_t misBarcodeInt;
	uint64_t returnBarcodeInt;
	int misPos;
	int misCount = 0;
	for (int i = 0; i < misMaskLen; i++)
	{
		misBarcodeInt = barcodeInt ^ misMask[i];
		if (barcodeMap->bmap.count(misBarcodeInt) > 0)
		{
			returnBarcodeInt = misBarcodeInt;
			misPos = i;
			misCount++;
			if (misCount > 1)
			{
				dupMisOverlapBarcode++;
				return MAX_BARCODE;
			}
		}
	}
	if (misCount == 1)
	{
		misPos = int(misPos / 3);
		if (misPos < mOptions->barcodeLen)
		{
			baseMisCount[misPos]++;
			baseTransCount[(barcodeInt >> (misPos * 2)) & 3][(returnBarcodeInt >> (misPos * 2)) & 3]++;
		}
		return returnBarcodeInt;
	}
	return MAX_BARCODE;
}
