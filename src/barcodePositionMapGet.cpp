/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "barcodePositionMapGet.h"

BarcodePositionMapGet::BarcodePositionMapGet(Options* opt)
{
	mOptions = opt;
}

BarcodePositionMapGet::~BarcodePositionMapGet()
{
}

void BarcodePositionMapGet::process()
{
	BarcodePositionMap barcodePositionMap(mOptions);
	barcodePositionMap.dumpbpmap(mOptions->out);

	if (!mOptions->report.empty()){
		/*
		ofstream reportWriter(mOptions->report);
		reportWriter << fixed << setprecision(2);
		for (int i = 0; i< barcodePositionMap.inFastqNumber; i++){	
			reportWriter << "getBarcodePositionMap_fastqFile:\t" << barcodePositionMap.inFile[i] << endl; 	
			reportWriter << "getBarcodePositionMap_TotalReads:\t" << barcodePositionMap.totalReads[i] << endl;
			reportWriter << "getBarcodePositionMap_TotalBases:\t" << barcodePositionMap.totalBase[i] << endl;
			reportWriter << "getBarcodePositionMap_ReadsWithoutPosition:\t" << barcodePositionMap.readsWithoutPos[i] << "\t" << (double)barcodePositionMap.readsWithoutPos[i] / (double)barcodePositionMap.totalBarcodes[i] * 100 << "%" << endl;
			reportWriter << "getBarcodePositionMap_ReadsWithN:\t" << barcodePositionMap.readsWithN[i] << "\t" << (double)barcodePositionMap.readsWithN[i]/(double)barcodePositionMap.totalBarcodes[i]*100 << "%" << endl;
			for (int j = 0; j < 4; j++){
				reportWriter << "getBarcodePositionMap_poly" << ATCG_BASES[j] << ":\t" << barcodePositionMap.polyReads[i][j] << "\t" << (double)barcodePositionMap.polyReads[i][j]/(double)(barcodePositionMap.totalReads[i]*barcodePositionMap.segment)*100 << "%" << endl;
			}
			reportWriter << "getBarcodePositionMap_dup:\t" << barcodePositionMap.dupReads[i] << "\t" << (double)barcodePositionMap.dupReads[i]/(double)barcodePositionMap.totalBarcodes[i]*100 << "%" << endl;
			reportWriter << "getBarcodePositionMap_ESTdup:\t" << barcodePositionMap.ESTdupReads[i] << "\t" << (double)barcodePositionMap.ESTdupReads[i]/(double)barcodePositionMap.totalBarcodes[i]*100 << "%" << endl;
			
			reportWriter << "getBarcodePositionMap_readQ10:\t" << (double)barcodePositionMap.readsQ10[i]/(double)barcodePositionMap.totalBase[i]*100 << "%" << endl;
			reportWriter << "getBarcodePositionMap_readQ20:\t" << (double)barcodePositionMap.readsQ20[i]/(double)barcodePositionMap.totalBase[i]*100 << "%" << endl;
			reportWriter << "getBarcodePositionMap_readQ30:\t" << (double)barcodePositionMap.readsQ30[i]/(double)barcodePositionMap.totalBase[i]*100 << "%" << endl;	
		}
		reportWriter << "getBarcodePositionMap_uniqBarcodeTypes:\t" << barcodePositionMap.bpmap.size() << endl;
		reportWriter << "getBarcidePositionMap_rc:\t" << mOptions->rcString << endl;
		reportWriter << "getBarcodePositionMap_barcodeStart:\t" << mOptions->barcodeStart << endl;
		reportWriter << "getBarcodePositionMap_barcodeLen:\t" << mOptions->barcodeLen << endl;
		reportWriter << "getBarcodePositionMap_barcodeSegment:\t" << mOptions->barcodeSegment << endl;
		reportWriter.close();
		*/
		//string htmlReportFile = mOptions->report+".html";
		
		HtmlReporter htmlReporter(mOptions);
		htmlReporter.printReport(mOptions->report, barcodePositionMap);

	}
}
