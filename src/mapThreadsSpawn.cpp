/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#include "mapThreadsSpawn.h"
#include "ThreadControl.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include "cmdline.h"
#include "barcodeToPositionMulti.h"
#include "JoinedMultiThreadPara.h"

void mapThreadsSpawn(Parameters &P, ReadAlignChunk **RAchunk)
{
	Options opt = P.bcOpt;
	// preprocess in ST_barcodeMapping
#ifdef TIME_STAT
	auto tStart = chrono::high_resolution_clock::now();
#endif
	BarcodeToPositionMulti barcodeToPosMulti(&opt);
#ifdef TIME_STAT
	auto tEnd = chrono::high_resolution_clock::now();
	P.loadMaskTime = chrono::duration<double, nano>(tEnd - tStart).count();
	cout << "time cost of loading mask file(s):\t" << P.loadMaskTime / 1000000000 << endl;
#endif
	barcodeToPosMulti.initOutput();
	barcodeToPosMulti.initPackRepositoey();

//  	MPI_Finalize();
//    std::thread producer(std::bind(&BarcodeToPositionMulti::producerTask, &barcodeToPosMulti));
#ifdef D_PRINT
	if (opt.isMgzInput)
	{
		cout << "mgz input" << endl;
	}
	else
	{
		cout << "gz input" << endl;
	}
#endif
	Result **results = new Result *[P.runThreadN];
	//    BarcodeProcessor** barcodeProcessors = new BarcodeProcessor*[P.runThreadN];
	for (int t = 0; t < P.runThreadN; t++)
	{
		results[t] = new Result(barcodeToPosMulti.mOptions, true);
		if (opt.useF14)
		{
#ifdef LOADH5_OPENMP
#if defined(OPT_MEMORY_PLAN_A) || defined(OPT_MEMORY_PLAN_B)
			results[t]->setBarcodeProcessor(&(barcodeToPosMulti.mbpmap->f14bpmap), barcodeToPosMulti.mbpmap->bloomFilter, &(barcodeToPosMulti.mbpmap->dnbCount));
#else
			results[t]->setBarcodeProcessor(&(barcodeToPosMulti.mbpmap->f14bpmap), barcodeToPosMulti.mbpmap->bloomFilter);
#endif
#else
			results[t]->setBarcodeProcessor(&(barcodeToPosMulti.mbpmap->f14bpmap));
#endif
		}
		else
		{
			results[t]->setBarcodeProcessor(&(barcodeToPosMulti.mbpmap->bpmap));
		}
	}

	//    std::thread** threads = new thread * [P.runThreadN];
	//    barcodeToPosMulti.process();
	string file1, file2;
	if (opt.isMgzInput)
	{
		file1 = P.reader->mLeft->mFilename;
		if (P.reader->mRight)
		{
			file2 = P.reader->mRight->mFilename;
		}
	}
	for (int ithread = 1; ithread < P.runThreadN; ithread++)
	{ // spawn threads
		//        JoinedMultiThreadPara* tmpPara;
		//        tmpPara->bmPara->stBcPara=results[ithread];
		//        tmpPara->bmPara->opt=opt;
		//        tmpPara->starPara=RAchunk[ithread];
		JoinedMultiThreadPara *tmpPara;
		if (opt.isMgzInput)
		{
			mgzipParser *threadMgzipParser = new mgzipParser(file1, file2, ithread, P.runThreadN);
			tmpPara = new JoinedMultiThreadPara(RAchunk[ithread], new BarcodeMappingThreadPara(results[ithread], opt, &barcodeToPosMulti, threadMgzipParser));
		}
		else
		{
			tmpPara = new JoinedMultiThreadPara(RAchunk[ithread], new BarcodeMappingThreadPara(results[ithread], opt, &barcodeToPosMulti));
		}

		int threadStatus = pthread_create(&g_threadChunks.threadArray[ithread], NULL, &g_threadChunks.threadRAprocessChunks, (void *)tmpPara);
		//        int threadStatus=pthread_create(&g_threadChunks.threadArray[ithread], NULL, &g_threadChunks.threadRAprocessChunks, (void *) RAchunk[ithread]);
		if (threadStatus > 0)
		{ // something went wrong with one of threads
			ostringstream errOut;
			string errMsg = "EXITING because of FATAL ERROR: phtread error while creating thread # " + to_string(ithread) + ", error code: " + to_string(threadStatus);
			sawErrCode(err_otherAPI_failed, errMsg);
			errOut << "EXITING because of FATAL ERROR: phtread error while creating thread # " << ithread << ", error code: " << threadStatus;
			exitWithError(errOut.str(), std::cerr, P.inOut->logMain, 1, P);
		};
		pthread_mutex_lock(&g_threadChunks.mutexLogMain);
		P.inOut->logMain << "Created thread # " << ithread << "\n"
						 << flush;
		pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
	};
	JoinedMultiThreadPara *t0para;
	if (opt.isMgzInput)
	{
		mgzipParser *threadMgzipParser = new mgzipParser(file1, file2, 0, P.runThreadN);
		t0para = new JoinedMultiThreadPara(RAchunk[0], new BarcodeMappingThreadPara(results[0], opt, &barcodeToPosMulti, threadMgzipParser));
	}
	else
	{
		t0para = new JoinedMultiThreadPara(RAchunk[0], new BarcodeMappingThreadPara(results[0], opt, &barcodeToPosMulti));
	}
	RAchunk[0]->processChunks(t0para->bmPara); // start main thread
	//    t0para->bmPara->stBcPara=results[0];
	//    t0para->bmPara->opt=opt;
	//    t0para->starPara=RAchunk[0];

	for (int ithread = 1; ithread < P.runThreadN; ithread++)
	{ // wait for all threads to complete
		int threadStatus = pthread_join(g_threadChunks.threadArray[ithread], NULL);
		if (threadStatus > 0)
		{ // something went wrong with one of threads
			ostringstream errOut;
			string errMsg = "EXITING because of FATAL ERROR: phtread error while joining thread # " + to_string(ithread) + ", error code: " + to_string(threadStatus);
			sawErrCode(err_otherAPI_failed, errMsg);
			errOut << "EXITING because of FATAL ERROR: phtread error while joining thread # " << ithread << ", error code: " << threadStatus;
			exitWithError(errOut.str(), std::cerr, P.inOut->logMain, 1, P);
		};
		pthread_mutex_lock(&g_threadChunks.mutexLogMain);
		P.inOut->logMain << "Joined thread # " << ithread << "\n"
						 << flush;
		pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
	};

	uint64_t totalQueryTime = 0;
#ifdef TIME_STAT
	uint64_t totalBcTime = 0;
	uint64_t totalMapTime = 0;
#endif
	vector<Result *> resultList;
	for (int t = 0; t < P.runThreadN; t++)
	{
		resultList.push_back(results[t]);
		totalQueryTime += results[t]->mBarcodeProcessor->queryTime;
#ifdef TIME_STAT
		totalBcTime += RAchunk[t]->bcTime;
		totalMapTime += RAchunk[t]->mapTime;
		cout << "accumulate time cost from all threads of BarcodeMapping(s):\t" << totalBcTime / 1000000000 << endl;
		cout << "accumulate time cost from all threads of StarMapping(s):\t" << totalMapTime / 1000000000 << endl;
#endif
	}

	//    cout<<"total query time:\t"<<totalQueryTime<<endl;
	//  pid_t pd=getpid();
	//  cout<<"current memory use "<<get_proc_virtualmem(pd)<<endl;
	Result *finalResult = Result::merge(resultList);
	finalResult->print();
	double bcMapRate = (double)finalResult->mBarcodeProcessor->mMapToSlideRead / (double)(finalResult->mTotalRead - finalResult->mFxiedFilterRead);
	if (bcMapRate < opt.bcMappingRateCutoff)
	{
		cerr << fixed << setprecision(2);
		string errMsg = "barcode mapping rate is too low";
		sawErrCode(err_cidRate_low, errMsg);
		cerr << "Error, barcode mapping rate is too low," << bcMapRate << endl;
		exit(1);
	}
	cout << resetiosflags(ios::fixed) << setprecision(2);
#if defined(OPT_MEMORY_PLAN_A) || defined(OPT_MEMORY_PLAN_B)
	uint64_t dnbCount = 0;
#endif
	if (!opt.transBarcodeToPos.mappedDNBOutFile.empty())
	{
#if !defined(OPT_MEMORY_PLAN_A) && !defined(OPT_MEMORY_PLAN_B)
		finalResult->dumpDNBs(opt.transBarcodeToPos.mappedDNBOutFile);
#else
#ifdef OPT_MEMORY_PLAN_A
		ofstream writer;
		writer.open(opt.transBarcodeToPos.mappedDNBOutFile);
		//	  for(auto f14Iter=barcodeToPosMulti.mbpmap->f14bpmap.begin();f14Iter!=barcodeToPosMulti.mbpmap->f14bpmap.end();f14Iter++){
		//		if(!((f14Iter->second)>>63)) {
		//		  if((((f14Iter->second)&countMask)>>countBits)>0) {
		//			// output count in barcodePos
		//			writer << (uint32_t)((f14Iter->second&coorXMask) >> coorBits) << "\t" << (uint32_t)(f14Iter->second&coorMask) << (uint32_t)((f14Iter->second&countMask) >> countBits) << "\n";
		//			dnbCount++;
		//		  }
		//		}else{
		//		  // output count overflow
		//		  uint64_t k=(uint64_t)(((f14Iter->second & coorXMask)>>coorBits)>>32) | (uint64_t)(f14Iter->second & coorMask);
		//		  writer << (uint32_t)((f14Iter->second & coorXMask)>>coorBits) << "\t" << (uint32_t)(f14Iter->second & coorMask)<<"\t"<<totalDnb[k]<<"\n";
		//		  dnbCount++;
		//		}
		//	  }
		for (unordered_map<uint64_t, int>::iterator ix = barcodeToPosMulti.mbpmap->dnbCount.begin(); ix != barcodeToPosMulti.mbpmap->dnbCount.end(); ix++)
		{
			writer << (uint32_t)((ix->first) >> 32) << "\t" << ((ix->first) & 0xffffffff) << "\t" << ix->second << endl;
			//		  writer << (uint32_t)((ix->first&coorXMask) >> coorBits)<<"\t"<<uint32_t(ix->first&coorMask)<<"\t"<<ix->second<<endl;
			dnbCount++;
		}
		writer.close();
		//	  barcodeToPosMulti.mbpmap->dnbCount.clear();
		unordered_map<uint64_t, int>().swap(barcodeToPosMulti.mbpmap->dnbCount);
		malloc_trim(0);
#endif
#ifdef OPT_MEMORY_PLAN_B
		ofstream writer;
		writer.open(opt.transBarcodeToPos.mappedDNBOutFile);
		if (opt.useF14)
		{
			for (auto f14Iter = barcodeToPosMulti.mbpmap->f14bpmap.begin(); f14Iter != barcodeToPosMulti.mbpmap->f14bpmap.end(); f14Iter++)
			{
				if (!(f14Iter->second.fullTag))
				{
					if (f14Iter->second.count > 0)
					{
						// output count in barcodePos
						writer << (uint32_t)(f14Iter->second.x) << "\t" << (uint32_t)(f14Iter->second.y) << "\t" << (uint32_t)(f14Iter->second.count) << "\n";
						dnbCount++;
					}
				}
				else
				{
					// output count overflow
					for (unordered_map<uint64_t, int>::iterator ix = barcodeToPosMulti.mbpmap->dnbCount.begin(); ix != barcodeToPosMulti.mbpmap->dnbCount.end(); ix++)
					{
						writer << (uint32_t)((ix->first) >> 32) << "\t" << ((ix->first) & 0xffffffff) << "\t" << ix->second << endl;
						//		  writer << (uint32_t)((ix->first&coorXMask) >> coorBits)<<"\t"<<uint32_t(ix->first&coorMask)<<"\t"<<ix->second<<endl;
						dnbCount++;
					}
				}
			}
		}
		else
		{
			for (auto bpIter = barcodeToPosMulti.mbpmap->bpmap.begin(); bpIter != barcodeToPosMulti.mbpmap->bpmap.end(); bpIter++)
			{
				if (!(bpIter->second.fullTag))
				{
					if (bpIter->second.count > 0)
					{
						// output count in barcodePos
						writer << (uint32_t)(bpIter->second.x) << "\t" << (uint32_t)(bpIter->second.y) << "\t" << (uint32_t)(bpIter->second.count) << "\n";
						dnbCount++;
					}
				}
				else
				{
					// output count overflow
					for (unordered_map<uint64_t, int>::iterator ix = barcodeToPosMulti.mbpmap->dnbCount.begin(); ix != barcodeToPosMulti.mbpmap->dnbCount.end(); ix++)
					{
						writer << (uint32_t)((ix->first) >> 32) << "\t" << ((ix->first) & 0xffffffff) << "\t" << ix->second << endl;
						//		  writer << (uint32_t)((ix->first&coorXMask) >> coorBits)<<"\t"<<uint32_t(ix->first&coorMask)<<"\t"<<ix->second<<endl;
						dnbCount++;
					}
				}
			}
		}
		writer.close();
		//	  barcodeToPosMulti.mbpmap->dnbCount.clear();
		unordered_map<uint64_t, int>().swap(barcodeToPosMulti.mbpmap->dnbCount);
		malloc_trim(0);
#endif
#endif
#if !defined(OPT_MEMORY_PLAN_A) && !defined(OPT_MEMORY_PLAN_B)
		cout << "mapped_dnbs: " << finalResult->mBarcodeProcessor->mDNB.size() << endl;
#else
		cout << "mapped_dnbs: " << dnbCount << endl;
#endif
	}
	if (opt.useF14)
	{
		barcodeToPosMulti.mbpmap->f14bpmap.clear();
		folly::F14ValueMap<uint64_t, Position1>().swap(barcodeToPosMulti.mbpmap->f14bpmap);
	}
	else
	{
		//	  t0para->bmPara->stBcPara->mBarcodeProcessor->bpmap->clear();
		//	  unordered_map<uint64_t, Position1>().swap(*(t0para->bmPara->stBcPara->mBarcodeProcessor->bpmap));
		barcodeToPosMulti.mbpmap->bpmap.clear();
		unordered_map<uint64_t, Position1> tmpMap;
		tmpMap.swap(barcodeToPosMulti.mbpmap->bpmap);
	}
	malloc_trim(0);
#ifdef D_PRINT
	cout << "before freemem:\t" << (float)get_proc_virtualMem(getpid()) / (1024 * 1024) << "\t" << (float)get_proc_VmRSS(getpid()) / (1024 * 1024) << endl;
#endif

//  	delete[] barcodeToPosMulti.mbpmap->bloomFilter->hashtableClassification;
#if !defined(OPT_MEMORY_PLAN_A) && !defined(OPT_MEMORY_PLAN_B)
	finalResult->mBarcodeProcessor->mDNB.clear();
	unordered_map<uint64_t, int>().swap(finalResult->mBarcodeProcessor->mDNB);
#endif

//    sleep(5);
#ifdef D_PRINT
	cout << "before freemem:\t" << (float)get_proc_virtualMem(getpid()) / (1024 * 1024) << "\t" << (float)get_proc_VmRSS(getpid()) / (1024 * 1024) << endl;
#endif
	//    cout<<get_local_time()<<"...mapThreadsSpawn done"<<endl;
}

Options getStBCpara(Parameters &P)
{

	string paraLine;
	Options opt;
	cmdline::parser cmd;
	// input/output
	cmd.add<string>("in", 'i', "the first sequencing fastq file path of read or bin file of  first sequencing", "");
	cmd.add<string>("in1", 'I', "the second sequencing fastq file path of read1 or bin file of second sequencing", false, "");
	cmd.add<string>("in2", 0, "the second sequencing fastq file or bin file path of read1", false, "");
	cmd.add<string>("encodeRule", 0, "encode Rule of 4 bases, default:ACGT which means A=>0,C=>1,G=>2,T=>3, and ACTG means T=>2,G=>3", false, "ACGT");
	cmd.add<string>("barcodeReadsCount", 0, "the mapped barcode list file with reads count per barcode.", false, "");
	cmd.add<string>("maskFile", 'M', "mask file path", false, "");
	cmd.add<string>("platform", 0, "sequence platform [SEQ500, T1, T5, T10]", false, "T1");
	cmd.add<string>("out", 'O', "output file prefix or fastq output file of read1", true, "");
	cmd.add<string>("report", 0, "logging file path.", false, "");
	cmd.add<int>("compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4.", false, 4);
	cmd.add<string>("unmappedOut", 0, "output file path for barcode unmapped reads, if this path isn't given, discard the reads.", false, "");
	cmd.add<string>("q10Tiff", 0, "generate q10 heatmap tiff when get barcode to position map.", false, "");
	cmd.add<string>("fov", 0, "fov range in format minCol-maxCol_minRow-maxRow.", false, "0-17_0-96");
	cmd.add<string>("chipID", 0, "chip id", false, "");
	cmd.add<int>("barcodeLen", 'l', "barcode length, default  is 25", false, 25);
	cmd.add<int>("barcodeStart", 0, "barcode start position", false, 0);
	cmd.add<int>("umiRead", 0, "read1 or read2 contains the umi sequence.", false, 1);
	cmd.add<int>("umiStart", 0, "umi start position. if the start postion is negative number, no umi sequence will be found", false, 40);
	cmd.add<int>("umiLen", 0, "umi length.", false, 10);
	cmd.add<string>("fixedSequence", 0, "fixed sequence in read1 that will be filtered.", false, "");
	cmd.add<int>("fixedStart", 0, "fixed sequence start position can by specied.", false, -1);
	cmd.add<string>("fixedSequenceFile", 0, "file contianing the fixed sequences and the start position, one sequence per line in the format: TGCCTCTCAG\t-1. when position less than 0, means wouldn't specified", false, "");
	cmd.add<int>("turnFovDegree", 0, "degree of every fov will be turned.", false, 0);
	cmd.add<long>("mapSize", 0, "bucket size of the new unordered_map.", false, 0);
	cmd.add<int>("mismatch", 0, "max mismatch is allowed for barcode overlap find.", false, 0);
	cmd.add<int>("segment", 0, "barcode segment number of first sequencing read.", false, 1);
	cmd.add<string>("rc", 0, "true means get the reverse complement barcode to barcode map. all means get both forward and reverse complement barcode to barcode map", false, "false");
	cmd.add<string>("readidSep", 0, "number of characters will be trimed from readid to get read position in slide. If not given this will be automatically get from fastq file.", false, "/");
	cmd.add<int>("action", 0, "chose one action you want to run [barcode_stat = 1, barcode_overlap = 2, get_barcode_position_map = 3, map_barcode_to_slide = 4, merge_barcode_list = 5].", false, 1);
	cmd.add<int>("thread", 'w', "number of thread that will be used to run.", false, 2);
	cmd.add("verbose", 'V', "output verbose log information (i.e. when every 1M reads are processed).");
	cmd.add("useF14", 0, "use F14 to replace unordered_map");
	cmd.add<uint64_t>("bcNum", 0, "barcode number in maskFile.", false, 0);
	cmd.add<float>("bcMappingRateCutoff", 0, "cutoff of barcode mapping rate.", false, 0);
	cmd.add<int>("polyAnum", 0, "cutoff of continuous A bases number which would be considered as a polyA.", false, 0);
	cmd.add<int>("mismatchInPolyA", 0, "mismatch in polyA.", false, 0);
	//    cmd.add("isMgzInput",0,"the format of input fastq files is mGzip");
	string paraFile = P.stBCparaFile;
	ifstream paraIn(paraFile);
	if (!paraIn)
	{
		string errMsg = "fail to open such file," + paraFile;
		sawErrCode(err_fileOpen_failed, errMsg);
		cerr << "Error,fail to open such file," << paraFile << endl;
		cerr << cmd.usage() << endl;
		exit(1);
	}
	opt.in = cmd.get<string>("in");
	opt.out = cmd.get<string>("out");

	opt.encodeRule = cmd.get<string>("encodeRule");

	opt.compression = cmd.get<int>("compression");
	opt.barcodeLen = cmd.get<int>("barcodeLen");
	opt.barcodeStart = cmd.get<int>("barcodeStart");
	opt.platform = cmd.get<string>("platform");
	opt.turnFovDegree = cmd.get<int>("turnFovDegree");
	opt.mapSize = cmd.get<long>("mapSize");
	opt.actionInt = cmd.get<int>("action");
	opt.verbose = cmd.exist("verbose");
	opt.useF14 = true;
	opt.bcNum = cmd.get<uint64_t>("bcNum");
	opt.bcMappingRateCutoff = cmd.get<float>("bcMappingRateCutoff");
	opt.polyAnum = cmd.get<int>("polyAnum");
	opt.mismatchInPolyA = cmd.get<int>("mismatchInPolyA");
	//    opt.isMgzInput=cmd.exist("isMgzInput");
	opt.isMgzInput = false;
	opt.thread = cmd.get<int>("thread");
	opt.rcString = cmd.get<string>("rc");
	opt.rc = opt.transRC(opt.rcString);
	opt.barcodeSegment = cmd.get<int>("segment");
	opt.maskFile = cmd.get<string>("maskFile");
	opt.report = cmd.get<string>("report");
	opt.drawHeatMap.fovRange = cmd.get<string>("fov");
	opt.chipID = cmd.get<string>("chipID");
	opt.drawHeatMap.maskFile = cmd.get<string>("maskFile");
	opt.drawHeatMap.q10TiffFile = cmd.get<string>("q10Tiff");
	opt.barcodeOverlap.in2 = cmd.get<string>("in1");
	opt.barcodeOverlap.mismatch = cmd.get<int>("mismatch");
	opt.barcodeStat.rcString = cmd.get<string>("rc");
	opt.barcodeStat.rc = opt.transRC(opt.barcodeStat.rcString);
	opt.barcodeStat.segment = cmd.get<int>("segment");
	opt.barcodeStat.readidSep = cmd.get<string>("readidSep");
	opt.transBarcodeToPos.in = cmd.get<string>("in");
	opt.transBarcodeToPos.in1 = cmd.get<string>("in1");
	opt.transBarcodeToPos.in2 = cmd.get<string>("in2");
	opt.transBarcodeToPos.mismatch = cmd.get<int>("mismatch");
	opt.transBarcodeToPos.unmappedOutFile = cmd.get<string>("unmappedOut");
	opt.transBarcodeToPos.umiRead = cmd.get<int>("umiRead");
	opt.transBarcodeToPos.umiStart = cmd.get<int>("umiStart");
	opt.transBarcodeToPos.umiLen = cmd.get<int>("umiLen");
	opt.transBarcodeToPos.mappedDNBOutFile = cmd.get<string>("barcodeReadsCount");
	opt.transBarcodeToPos.fixedSequence = cmd.get<string>("fixedSequence");
	opt.transBarcodeToPos.fixedStart = cmd.get<int>("fixedStart");
	opt.transBarcodeToPos.fixedSequenceFile = cmd.get<string>("fixedSequenceFile");
	while (getline(paraIn, paraLine))
	{
		// if(paraLine.empty() || paraLine.find("#")==0 || paraLine.find("//")==0){
		// 	continue;
		// }
		vector<string> eles;
		line_split(paraLine, '=',eles);
		if (eles[0] == "in")
		{
			opt.in = eles[1];
			opt.transBarcodeToPos.in = eles[1];
		}
		if (eles[0] == "out")
		{
			opt.out = eles[1];
		}

		if (eles[0] == "encodeRule")
		{
			opt.encodeRule = eles[1];
		}

		if (eles[0] == "compression")
		{
			opt.compression = atoi(eles[1].c_str());
		}
		if (eles[0] == "barcodeLen")
		{
			opt.barcodeLen = atoi(eles[1].c_str());
		}
		if (eles[0] == "barcodeStart")
		{
			opt.barcodeStart = atoi(eles[1].c_str());
		}
		if (eles[0] == "platform")
		{
			opt.platform = eles[1];
		}
		if (eles[0] == "turnFovDegree")
		{
			opt.turnFovDegree = atoi(eles[1].c_str());
		}
		if (eles[0] == "mapSize")
		{
			opt.mapSize = atoi(eles[1].c_str());
		}
		if (eles[0] == "action")
		{
			opt.actionInt = atoi(eles[1].c_str());
		}
		if (eles[0] == "verbose")
		{
			opt.verbose = true;
		}
		if (eles[0] == "useF14")
		{
			opt.useF14 = true;
		}
		if (eles[0] == "bcNum")
		{
			opt.bcNum = atol(eles[1].c_str());
		}
		if (eles[0] == "bcMappingRateCutoff")
		{
			// cromwell
			// check the value need to reset to 1 or not
			// if exist call-ref_read dir && call-ref_read is a cacheCopy && find low bcMappingRate in bcStar log, then change the value to 1
			opt.bcMappingRateCutoff = atof(eles[1].c_str());
			if (opt.bcMappingRateCutoff < 0 || opt.bcMappingRateCutoff > 1)
			{
				string errMsg = "bcMappingRateCutoff should be a float which in [0,1]";
				sawErrCode(err_param_invalid_bcmap, errMsg);
				cerr << "bcMappingRateCutoff should be a float which in [0,1]" << endl;
				exit(1);
			}
			string currentShardIdx = "";
			string expectedCromwellPath = "../../../";
			string refReadRcPath;
			string refReadRcPath1 = expectedCromwellPath + "call-ref_read/cacheCopy/execution/rc";
			string refReadRcPath2 = expectedCromwellPath + "call-RefRead/cacheCopy/execution/rc";
			if (access(refReadRcPath1.c_str(), F_OK) == 0)
			{
				refReadRcPath = refReadRcPath1;
			}
			else if (access(refReadRcPath2.c_str(), F_OK) == 0)
			{
				refReadRcPath = refReadRcPath2;
			}
			char *buf1 = new char[500];
			char *buf2 = new char[500];
			char *buf3 = new char[500];
			string absRefReadRcPath = "";
			if (realpath(refReadRcPath.c_str(), buf3) != NULL)
			{
				absRefReadRcPath = buf3;
			}
			if (access(absRefReadRcPath.c_str(), F_OK) == 0)
			{
				if (realpath("./script", buf1) != NULL)
				{
					string realPath = buf1;
					vector<string> eles;
					line_split(realPath, '/', eles);
					currentShardIdx = eles[eles.size() - 3];
					if (currentShardIdx.find("shard") != 0)
					{
						cerr << "Error, code error" << endl;
						exit(1);
					}
				}
				string lastCromwellPath = "";
				if (readlink(refReadRcPath.c_str(), buf2, 500) != -1)
				{
					string rawPath = buf2;
					vector<string> eles;
					line_split(rawPath, '/', eles);
					for (uint64_t i = 0; i < eles.size() - 3; i++)
					{
						lastCromwellPath += "/" + eles[i];
					}
				}
				string lastBcStarPath = "";
				string lastBcStarPath1 = lastCromwellPath + "/call-barcodeMappingAndStar";
				string lastBcStarPath2 = lastCromwellPath + "/call-BarcodeMappingAndStar";
				if (access(lastBcStarPath1.c_str(), F_OK) == 0)
				{
					lastBcStarPath = lastBcStarPath1;
				}
				else if (access(lastBcStarPath2.c_str(), F_OK) == 0)
				{
					lastBcStarPath = lastBcStarPath2;
				}
				if (!lastBcStarPath.empty())
				{
					string lastStdErr = lastBcStarPath + "/" + currentShardIdx + "/execution/stderr";
					if (access(lastStdErr.c_str(), F_OK) == 0)
					{
						ifstream stderrIn(lastStdErr);
						if (stderrIn)
						{
							string tmpLine;
							while (getline(stderrIn, tmpLine))
							{
								if (tmpLine.find("Error, barcode mapping rate is too low") != string::npos)
								{
									cerr << "change bc mapRate to 0" << endl;
									opt.bcMappingRateCutoff = 0.0;
									break;
								}
							}
						}
						stderrIn.close();
					}
				}
			}
			delete[] buf1;
			delete[] buf2;
		}
		//		if(eles[0]=="isMgzInput"){
		//		  opt.isMgzInput=true;
		//		}
		if (eles[0] == "mismatchInPolyA")
		{
			opt.mismatchInPolyA = atoi(eles[1].c_str());
		}
		if (eles[0] == "polyAnum")
		{
			opt.polyAnum = atoi(eles[1].c_str());
		}
		if (eles[0] == "thread")
		{
			opt.thread = atoi(eles[1].c_str());
		}
		if (eles[0] == "rc")
		{
			opt.rcString = atoi(eles[1].c_str());
		}
		if (eles[0] == opt.rcString)
		{
			opt.rc = atoi(eles[1].c_str());
		}
		if (eles[0] == "segment")
		{
			opt.barcodeSegment = atoi(eles[1].c_str());
		}
		if (eles[0] == "maskFile")
		{
			opt.maskFile = eles[1];
		}
		if (eles[0] == "report")
		{
			opt.report = eles[1];
		}
		// if (eles[0] == "fov")
		// {
		// 	opt.drawHeatMap.fovRange = eles[1];
		// }
		// if (eles[0] == "chipID")
		// {
		// 	opt.chipID = eles[1];
		// }
		if (eles[0] == "maskFile")
		{
			opt.drawHeatMap.maskFile = eles[1];
		}
		// if (eles[0] == "q10TiffFile")
		// {
		// 	opt.drawHeatMap.q10TiffFile = eles[1];
		// }
		if (eles[0] == "in1")
		{
			opt.barcodeOverlap.in2 = eles[1];
		}
		if (eles[0] == "mismatch")
		{
			opt.barcodeOverlap.mismatch = atoi(eles[1].c_str());
		}
		if (eles[0] == "rc")
		{
			opt.rcString = eles[1];
			opt.rc = opt.transRC(opt.rcString);
			opt.barcodeStat.rcString = opt.rcString;
			opt.barcodeStat.rc = opt.transRC(opt.barcodeStat.rcString);
		}
		if (eles[0] == "segment")
		{
			opt.barcodeStat.segment = atoi(eles[1].c_str());
		}
		if (eles[0] == "readidSep")
		{
			opt.barcodeStat.readidSep = eles[1];
		}
		if (eles[0] == "in1")
		{
			opt.transBarcodeToPos.in1 = eles[1];
		}
		if (eles[0] == "in2")
		{
			opt.transBarcodeToPos.in2 = eles[1];
		}
		if (eles[0] == "mismatch")
		{
			opt.transBarcodeToPos.mismatch = atoi(eles[1].c_str());
		}
		if (eles[0] == "unmappedOut")
		{
			opt.transBarcodeToPos.unmappedOutFile = eles[1];
		}
		if (eles[0] == "umiRead")
		{
			opt.transBarcodeToPos.umiRead = atoi(eles[1].c_str());
		}
		if (eles[0] == "umiStart")
		{
			opt.transBarcodeToPos.umiStart = atoi(eles[1].c_str());
		}
		if (eles[0] == "umiLen")
		{
			opt.transBarcodeToPos.umiLen = atoi(eles[1].c_str());
		}
		if (eles[0] == "barcodeReadsCount")
		{
			opt.transBarcodeToPos.mappedDNBOutFile = eles[1];
		}
		if (eles[0] == "fixedSequence")
		{
			opt.transBarcodeToPos.fixedSequence = eles[1];
		}
		if (eles[0] == "fixedStart")
		{
			opt.transBarcodeToPos.fixedStart = atoi(eles[1].c_str());
		}
		if (eles[0] == "fixedSequenceFile")
		{
			opt.transBarcodeToPos.fixedSequenceFile = eles[1];
		}
	}
	opt.init();
	opt.validate();
	return opt;
}