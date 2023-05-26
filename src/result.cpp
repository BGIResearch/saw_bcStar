/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "result.h"

Result::Result(Options* opt, int threadId, bool paired)
{
	mOptions = opt;
	mPaired = paired;
	mTotalRead = 0;
	mFxiedFilterRead = 0;
	mWithoutPositionReads = 0;
	mLowQuaRead = 0;
	overlapReadsWithMis = 0;
	overlapReadsWithN = 0;

	mThreadId = threadId;
	if (opt->verbose) {
		string msg = "new result class in thread: " + to_string(mThreadId);
//		loginfo(msg);
	}
}

Result::~Result()
{
	delete mBarcodeProcessor;
}

Result* Result::merge(vector<Result*>& list)
{
	if (list.size() == 0)
		return nullptr;
	Result* result = new Result(list[0]->mOptions, 0, list[0]->mPaired);
	result->setBarcodeProcessor();

	for (uint64_t i = 0; i < list.size(); i++) {
		result->mTotalRead += list[i]->mTotalRead;
		result->totalBases+=list[i]->totalBases;
		result->mFxiedFilterRead += list[i]->mFxiedFilterRead;
		result->mWithoutPositionReads += list[i]->mWithoutPositionReads;
		result->mLowQuaRead += list[i]->mLowQuaRead;
		result->seqQ10+=list[i]->seqQ10;
		result->seqQ20+=list[i]->seqQ20;
		result->seqQ30+=list[i]->seqQ30;
	  	result->tooManyNinBarcode+=list[i]->tooManyNinBarcode;
		  result->containPolyA+=list[i]->containPolyA;
		  result->filterByPolyA+=list[i]->filterByPolyA;
		if(result->readLen==0) {
		  result->readLen = list[i]->readLen;
		}
		result->mBarcodeProcessor->totalReads += list[i]->mBarcodeProcessor->totalReads;
		result->mBarcodeProcessor->mMapToSlideRead += list[i]->mBarcodeProcessor->mMapToSlideRead;
		result->mBarcodeProcessor->overlapReads += list[i]->mBarcodeProcessor->overlapReads;
		result->mBarcodeProcessor->overlapReadsWithMis += list[i]->mBarcodeProcessor->overlapReadsWithMis;
		result->mBarcodeProcessor->overlapReadsWithN += list[i]->mBarcodeProcessor->overlapReadsWithN;
		result->mBarcodeProcessor->barcodeQ10 += list[i]->mBarcodeProcessor->barcodeQ10;
		result->mBarcodeProcessor->barcodeQ20 += list[i]->mBarcodeProcessor->barcodeQ20;
		result->mBarcodeProcessor->barcodeQ30 += list[i]->mBarcodeProcessor->barcodeQ30;
		result->mBarcodeProcessor->umiQ10 += list[i]->mBarcodeProcessor->umiQ10;
		result->mBarcodeProcessor->umiQ20 += list[i]->mBarcodeProcessor->umiQ20;
		result->mBarcodeProcessor->umiQ30 += list[i]->mBarcodeProcessor->umiQ30;
		result->mBarcodeProcessor->umiNFilterReads += list[i]->mBarcodeProcessor->umiNFilterReads;
		result->mBarcodeProcessor->umiPloyAFilterReads += list[i]->mBarcodeProcessor->umiPloyAFilterReads;
		result->mBarcodeProcessor->umiQ10FilterReads += list[i]->mBarcodeProcessor->umiQ10FilterReads;
#ifndef OPT_MEMORY_PLAN_A
		if (!list[i]->mOptions->transBarcodeToPos.mappedDNBOutFile.empty()) {
			unordered_map<uint64_t, int>::iterator mergeIter;
			for (auto iter = list[i]->mBarcodeProcessor->mDNB.begin(); iter != list[i]->mBarcodeProcessor->mDNB.end(); iter++) {
				mergeIter = result->mBarcodeProcessor->mDNB.find(iter->first);
				if (mergeIter != result->mBarcodeProcessor->mDNB.end()) {
					mergeIter->second += iter->second;
				}
				else {
					result->mBarcodeProcessor->mDNB[iter->first] = iter->second;
				}
			}
			list[i]->mBarcodeProcessor->mDNB.clear();
			unordered_map<uint64_t, int>().swap(list[i]->mBarcodeProcessor->mDNB);
		}
#endif
	}
	return result;
}

void Result::print()
{
	
	cout << fixed << setprecision(2);
	cout << "total_reads:\t" << mTotalRead << endl;
//	cout << "fixed_sequence_contianing_reads:\t" << mFxiedFilterRead << "\t" << (double)mFxiedFilterRead / (double)mTotalRead * 100 << "%" << endl;
//	cout << "pass_filter_reads:\t" << mTotalRead << endl;
	cout << "reads_with_polyA:\t"<<containPolyA<<"\t"<<(double)containPolyA/(double)mTotalRead*100<<"%"<<endl;
	cout << "reads_filteredByPolyA:\t"<<filterByPolyA<<"\t"<<(double)filterByPolyA/(double)mTotalRead*100<<"%"<<endl;
	cout << "mapped_reads:\t" << mBarcodeProcessor->mMapToSlideRead << "\t" << (double)mBarcodeProcessor->mMapToSlideRead / (double)(mTotalRead - mFxiedFilterRead)* 100 << "%" << endl;
	cout << "reads_with_adapter:\t0\t0.00%" << endl;
	cout << "reads_with_dnb:\t0\t0.00%" << endl;
	double overlapReadsWithMisRate = (double)mBarcodeProcessor->overlapReadsWithMis / (double)mTotalRead * 100;
	double overlapReadsWithNRate = (double)mBarcodeProcessor->overlapReadsWithN / (double)mTotalRead * 100;
	cout << "barcode_exactlyOverlap_reads:\t" << mBarcodeProcessor->overlapReads << "\t" << (double)mBarcodeProcessor->overlapReads / (double)mTotalRead * 100 << "%" << endl;
	cout << "barcode_misOverlap_reads:\t" << mBarcodeProcessor->overlapReadsWithMis << "\t" << overlapReadsWithMisRate << "%" << endl;
	cout << "barcode_withN_reads:\t" << mBarcodeProcessor->overlapReadsWithN << "\t" << overlapReadsWithNRate << "%" << endl;
  //	cout<<"filtered barcode with too many N bases:\t"<<tooManyNinBarcode<<endl;
	double barcodeQ10 = (double)mBarcodeProcessor->barcodeQ10 / (double)(mTotalRead * mOptions->barcodeLen) * 100;
	double barcodeQ20 = (double)mBarcodeProcessor->barcodeQ20 / (double)(mTotalRead * mOptions->barcodeLen) * 100;
	double barcodeQ30 = (double)mBarcodeProcessor->barcodeQ30 / (double)(mTotalRead * mOptions->barcodeLen) * 100;
	double umiQ10 = 0;
	double umiQ20 = 0;
	double umiQ30 = 0;
	if (mOptions->transBarcodeToPos.umiLen>0){
		umiQ10 = (double)mBarcodeProcessor->umiQ10 / (double)(mTotalRead * mOptions->transBarcodeToPos.umiLen) * 100;
		umiQ20 = (double)mBarcodeProcessor->umiQ20 / (double)(mTotalRead * mOptions->transBarcodeToPos.umiLen) * 100;
		umiQ30 = (double)mBarcodeProcessor->umiQ30 / (double)(mTotalRead * mOptions->transBarcodeToPos.umiLen) * 100;
	}
	
	cout << "Q10_bases_in_barcode:\t" << barcodeQ10 << "%" << endl;
	cout << "Q20_bases_in_barcode:\t" << barcodeQ20 << "%" << endl;
	cout << "Q30_bases_in_barcode:\t" << barcodeQ30 << "%" << endl;
	cout << "Q10_bases_in_umi:\t" << umiQ10 << "%" << endl;
	cout << "Q20_bases_in_umi:\t" << umiQ20 << "%" << endl;
	cout << "Q30_bases_in_umi:\t" << umiQ30 << "%" << endl;
//	double seqQ10Ratio=(double)seqQ10/(double)(mTotalRead * readLen) * 100;
//	double seqQ20Ratio=(double)seqQ20/(double)(mTotalRead * readLen) * 100;
//	double seqQ30Ratio=(double)seqQ30/(double)(mTotalRead * readLen) * 100;
	double seqQ10Ratio=(double)seqQ10/(double)totalBases * 100;
	double seqQ20Ratio=(double)seqQ20/(double)totalBases * 100;
	double seqQ30Ratio=(double)seqQ30/(double)totalBases * 100;
	cout << "Q10_bases_in_seq:\t" << seqQ10Ratio << "%" << endl;
	cout << "Q20_bases_in_seq:\t" << seqQ20Ratio << "%" << endl;
	cout << "Q30_bases_in_seq:\t" << seqQ30Ratio << "%" << endl;
	long umiTotalFilterReads = mBarcodeProcessor->umiNFilterReads + mBarcodeProcessor->umiPloyAFilterReads + mBarcodeProcessor->umiQ10FilterReads;
	double umiFilterRate = (double)umiTotalFilterReads/(double)mTotalRead * 100;
	double umiNFilterRate = (double)mBarcodeProcessor->umiNFilterReads/(double)mTotalRead * 100;
	double umiPolyAFilterRate = (double)mBarcodeProcessor->umiPloyAFilterReads/(double)mTotalRead * 100;
	double umiQ10FilterRate = (double)mBarcodeProcessor->umiQ10FilterReads/(double)mTotalRead * 100;
	cout << "umi_filter_reads:\t" << umiTotalFilterReads << "\t" << umiFilterRate << "%" << endl;
	cout << "umi_with_N_reads:\t" << mBarcodeProcessor->umiNFilterReads << "\t" << umiNFilterRate << "%" << endl;
	cout << "umi_with_polyA_reads:\t" << mBarcodeProcessor->umiPloyAFilterReads << "\t" << umiPolyAFilterRate << "%" << endl;
	cout << "umi_with_low_quality_base_reads:\t" << mBarcodeProcessor->umiQ10FilterReads << "\t" << umiQ10FilterRate << "%" << endl;
}

void Result::dumpDNBs(string& mapOutFile)
{
	mBarcodeProcessor->dumpDNBmap(mapOutFile);
}

void Result::setBarcodeProcessor(unordered_map<uint64_t, Position1>* bpmap)
{
	mBarcodeProcessor = new BarcodeProcessor(mOptions, bpmap);
}
void Result::setBarcodeProcessor(folly::F14ValueMap<uint64_t, Position1>* bpmap)
{
  mBarcodeProcessor = new BarcodeProcessor(mOptions, bpmap);
}
void Result::setBarcodeProcessor(unordered_map<uint64_t, Position2>* bpmap)
{
  mBarcodeProcessor = new BarcodeProcessor(mOptions, bpmap);
}
void Result::setBarcodeProcessor(folly::F14ValueMap<uint64_t, Position2>* bpmap)
{
  mBarcodeProcessor = new BarcodeProcessor(mOptions, bpmap);
}
void Result::setBarcodeProcessor()
{
	mBarcodeProcessor = new BarcodeProcessor();
}
void Result::setBarcodeProcessor(folly::F14ValueMap<uint64_t, Position1> *bpmap, BloomFilter *pbf) {
  mBarcodeProcessor = new BarcodeProcessor(mOptions, bpmap,pbf);
}
#if defined(OPT_MEMORY_PLAN_A) || defined(OPT_MEMORY_PLAN_B)
void Result::setBarcodeProcessor(folly::F14ValueMap<uint64_t, Position1> *bpmap, BloomFilter *pbf, unordered_map<uint64_t, int> *dnbCount) {
  mBarcodeProcessor = new BarcodeProcessor(mOptions, bpmap,pbf,dnbCount);
}
#endif


